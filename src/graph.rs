//! `DbgGraph`: bidirected de Bruijn graph built on top of petgraph's `StableGraph`.

use std::collections::HashMap;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::hash::BuildHasherDefault;
use std::io::Write;

use nohash_hasher::NoHashHasher;
use petgraph::algo::connected_components as petgraph_connected_components;
use petgraph::algo::tarjan_scc;
use petgraph::dot::{Config, Dot};
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use petgraph::Direction::{Incoming, Outgoing};

use crate::node::{EmptyEdge, NodeStruct};
use crate::types::{
    CarryType, EdgeId, EdgeType, EdgeWeight, HashInfoSimple, Idx, IndexedEdge, NodeId,
};

/// Inner petgraph type alias.
type Inner = petgraph::stable_graph::StableGraph<NodeStruct, EmptyEdge, petgraph::Directed, Idx>;
type BackendNodeIndex = petgraph::stable_graph::NodeIndex<Idx>;
type BackendEdgeIndex = petgraph::stable_graph::EdgeIndex<Idx>;
type ExportEdgeKey = (NodeId, EdgeType, NodeId);

#[inline]
fn to_backend_node(node: NodeId) -> BackendNodeIndex {
    BackendNodeIndex::new(node.0)
}

#[inline]
fn from_backend_node(node: BackendNodeIndex) -> NodeId {
    NodeId(node.index())
}

#[inline]
fn to_backend_edge(edge: EdgeId) -> BackendEdgeIndex {
    BackendEdgeIndex::new(edge.0)
}

#[inline]
fn from_backend_edge(edge: BackendEdgeIndex) -> EdgeId {
    EdgeId(edge.index())
}

/// Candidate outgoing edge from the start of a potential bubble.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BubbleStartEdge {
    pub target: NodeId,
    pub edge_type: EdgeType,
}

/// A connection invariant violated by a graph.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GraphValidationIssue {
    /// A directed edge has fewer reverse-type partners than its multiplicity requires.
    UnpairedEdge {
        edge_id: EdgeId,
        from: NodeId,
        to: NodeId,
        edge_type: EdgeType,
    },
    /// More than one identical directed edge connects the same pair of nodes.
    DuplicateEdge {
        from: NodeId,
        to: NodeId,
        edge_type: EdgeType,
        count: usize,
    },
}

/// All connection invariant violations found by [`DbgGraph::validate`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct GraphValidationReport {
    /// Issues are reported in deterministic connection/edge order.
    pub issues: Vec<GraphValidationIssue>,
}

impl GraphValidationReport {
    /// Whether validation found no issues.
    pub fn is_empty(&self) -> bool {
        self.issues.is_empty()
    }
}

impl fmt::Display for GraphValidationReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "graph validation failed with {} issue(s)",
            self.issues.len()
        )
    }
}

impl std::error::Error for GraphValidationReport {}

/// Bidirected de Bruijn graph.
pub struct DbgGraph {
    inner: Inner,
    k: usize,
}

#[inline]
fn assert_valid_k(k: usize) {
    assert!(k > 0, "k-mer length must be at least 1");
}

impl Default for DbgGraph {
    fn default() -> Self {
        Self::new(1)
    }
}

// ─── Construction ────────────────────────────────────────────────────────────

impl DbgGraph {
    /// Create an empty graph for the given k-mer length.
    pub fn new(k: usize) -> Self {
        assert_valid_k(k);

        DbgGraph {
            inner: Inner::default(),
            k,
        }
    }

    /// Build a de Bruijn graph from the k-mer map produced by preprocessing.
    ///
    /// `map` is a `HashMap<canonical_hash, HashInfoSimple, ...>` as returned by
    /// `preprocessing_standalone` / `preprocessing_wasm`.
    ///
    /// The map must be closed under neighbour references: every hash appearing in
    /// a `pre` or `post` entry must also be present as a key in `map`. Missing
    /// neighbour entries cause this constructor to panic.
    pub fn from_kmer_map(
        k: usize,
        map: &HashMap<u64, HashInfoSimple, BuildHasherDefault<NoHashHasher<u64>>>,
    ) -> Self {
        assert_valid_k(k);

        let mut g = DbgGraph {
            inner: Inner::with_capacity(map.len(), map.len() * 2),
            k,
        };

        let mut tmpdict: HashMap<u64, BackendNodeIndex, BuildHasherDefault<NoHashHasher<u64>>> =
            HashMap::with_capacity_and_hasher(map.len(), BuildHasherDefault::default());

        for (h, hi) in map {
            // First, check if the node exists.
            if !tmpdict.contains_key(h) {
                let nid = g.inner.add_node(NodeStruct {
                    counts: hi.counts,
                    abs_ind: vec![*h],
                    innerdir: None,
                });
                tmpdict.insert(*h, nid);
            }

            // ...and also all of its preceding edges...
            for hpre in hi.pre.iter() {
                if let std::collections::hash_map::Entry::Vacant(e) = tmpdict.entry(hpre.0) {
                    // First, we add the node.
                    let nid2 = g.inner.add_node(NodeStruct {
                        counts: map.get(&hpre.0).unwrap().counts,
                        abs_ind: vec![hpre.0],
                        innerdir: None,
                    });
                    e.insert(nid2);
                }
                // Then, we add the edge.
                if !g
                    .inner
                    .edges_connecting(*tmpdict.get(&hpre.0).unwrap(), *tmpdict.get(h).unwrap())
                    .any(|e| e.weight().t == hpre.1)
                {
                    g.inner.add_edge(
                        *tmpdict.get(&hpre.0).unwrap(),
                        *tmpdict.get(h).unwrap(),
                        EmptyEdge { t: hpre.1 },
                    );
                }
            }

            // ...and the forward ones.
            for hpost in hi.post.iter() {
                if let std::collections::hash_map::Entry::Vacant(e) = tmpdict.entry(hpost.0) {
                    // First, we add the node.
                    let nid2 = g.inner.add_node(NodeStruct {
                        counts: map.get(&hpost.0).unwrap().counts,
                        abs_ind: vec![hpost.0],
                        innerdir: None,
                    });
                    e.insert(nid2);
                }
                // Then, we add the edge.
                if !g
                    .inner
                    .edges_connecting(*tmpdict.get(h).unwrap(), *tmpdict.get(&hpost.0).unwrap())
                    .any(|e| e.weight().t == hpost.1)
                {
                    g.inner.add_edge(
                        *tmpdict.get(h).unwrap(),
                        *tmpdict.get(&hpost.0).unwrap(),
                        EmptyEdge { t: hpost.1 },
                    );
                }
            }
        }

        log::debug!(
            "DbgGraph::from_kmer_map: {} nodes, {} edges",
            g.inner.node_count(),
            g.inner.edge_count()
        );

        g
    }

    /// Build a de Bruijn graph from aligned k-mer vectors.
    ///
    /// Neighbours for each source are stored in one vector. The first
    /// `predecessor_counts[source]` entries are incoming; the remainder are outgoing. All input
    /// vectors are consumed so their allocations can be released as construction progresses.
    pub fn from_indexed_kmers(
        k: usize,
        canonical_hashes: Vec<u64>,
        counts: Vec<EdgeWeight>,
        neighbours: Vec<Vec<IndexedEdge>>,
        predecessor_counts: Vec<u8>,
    ) -> Self {
        assert_valid_k(k);

        let node_count = canonical_hashes.len();
        assert_eq!(
            counts.len(),
            node_count,
            "indexed k-mer hashes and counts have different lengths"
        );
        assert_eq!(
            neighbours.len(),
            node_count,
            "indexed k-mer hashes and neighbour lists have different lengths"
        );
        assert_eq!(
            predecessor_counts.len(),
            node_count,
            "indexed k-mer hashes and predecessor counts have different lengths"
        );

        for (source, (edges, predecessor_count)) in
            neighbours.iter().zip(predecessor_counts.iter()).enumerate()
        {
            assert!(
                usize::from(*predecessor_count) <= edges.len(),
                "predecessor count exceeds neighbour-list length for k-mer index {source}"
            );
            for edge in edges {
                assert!(
                    edge.target < node_count,
                    "neighbour target {} is outside the indexed k-mer table of length {node_count}",
                    edge.target
                );
            }
        }

        let mut g = DbgGraph {
            inner: Inner::with_capacity(node_count, node_count.saturating_mul(2)),
            k,
        };

        for (expected_index, (hash, count)) in canonical_hashes.into_iter().zip(counts).enumerate()
        {
            let actual = g.inner.add_node(NodeStruct {
                counts: count,
                abs_ind: vec![hash],
                innerdir: None,
            });
            assert_eq!(
                actual.index(),
                expected_index,
                "sequential graph-node insertion did not preserve indexed k-mer order"
            );
        }

        for (current, (edges, predecessor_count)) in
            neighbours.into_iter().zip(predecessor_counts).enumerate()
        {
            let current = BackendNodeIndex::new(current);
            let split = usize::from(predecessor_count);

            for edge in &edges[..split] {
                let predecessor = BackendNodeIndex::new(edge.target);
                if !g
                    .inner
                    .edges_connecting(predecessor, current)
                    .any(|existing| existing.weight().t == edge.edge_type)
                {
                    g.inner
                        .add_edge(predecessor, current, EmptyEdge { t: edge.edge_type });
                }
            }

            for edge in &edges[split..] {
                let successor = BackendNodeIndex::new(edge.target);
                if !g
                    .inner
                    .edges_connecting(current, successor)
                    .any(|existing| existing.weight().t == edge.edge_type)
                {
                    g.inner
                        .add_edge(current, successor, EmptyEdge { t: edge.edge_type });
                }
            }
        }

        log::debug!(
            "DbgGraph::from_indexed_kmers: {} nodes, {} edges",
            g.inner.node_count(),
            g.inner.edge_count()
        );

        g
    }

    /// Insert a bare node (used in tests and programmatic construction).
    pub fn add_node(&mut self, node: NodeStruct) -> NodeId {
        from_backend_node(self.inner.add_node(node))
    }
}

// ─── Graph info ──────────────────────────────────────────────────────────────

impl DbgGraph {
    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.inner.node_count()
    }

    /// Number of edges.
    pub fn edge_count(&self) -> usize {
        self.inner.edge_count()
    }

    /// k-mer length.
    pub fn k(&self) -> usize {
        self.k
    }

    /// Whether node index `n` is still present.
    pub fn contains_node(&self, n: NodeId) -> bool {
        self.inner.contains_node(to_backend_node(n))
    }

    /// Iterator over all live node indices.
    pub fn node_indices(&self) -> impl Iterator<Item = NodeId> + '_ {
        self.inner.node_indices().map(from_backend_node)
    }

    /// Count nodes without incoming or outgoing neighbours.
    pub fn isolated_node_count(&self) -> usize {
        self.inner
            .node_indices()
            .filter(|node| {
                self.inner.neighbors_directed(*node, Incoming).count() == 0
                    && self.inner.neighbors_directed(*node, Outgoing).count() == 0
            })
            .count()
    }

    /// Strongly connected components, using petgraph's current backend implementation.
    pub fn strongly_connected_components(&self) -> Vec<Vec<NodeId>> {
        tarjan_scc(&self.inner)
            .into_iter()
            .map(|component| component.into_iter().map(from_backend_node).collect())
            .collect()
    }
}

// ─── Node / edge weight access ───────────────────────────────────────────────

impl DbgGraph {
    /// Immutable reference to node weight.
    pub fn node_weight(&self, n: NodeId) -> Option<&NodeStruct> {
        self.inner.node_weight(to_backend_node(n))
    }

    /// Mutable reference to node weight.
    pub fn node_weight_mut(&mut self, n: NodeId) -> Option<&mut NodeStruct> {
        self.inner.node_weight_mut(to_backend_node(n))
    }

    /// Immutable reference to edge weight.
    pub fn edge_weight(&self, e: EdgeId) -> Option<&EmptyEdge> {
        self.inner.edge_weight(to_backend_edge(e))
    }

    /// Mutable reference to edge weight.
    pub fn edge_weight_mut(&mut self, e: EdgeId) -> Option<&mut EmptyEdge> {
        self.inner.edge_weight_mut(to_backend_edge(e))
    }

    /// Endpoints of an edge.
    pub fn edge_endpoints(&self, e: EdgeId) -> Option<(NodeId, NodeId)> {
        self.inner
            .edge_endpoints(to_backend_edge(e))
            .map(|(from, to)| (from_backend_node(from), from_backend_node(to)))
    }

    /// Return the total number of k-mers represented by an ordered path.
    ///
    /// Returns `None` if the path contains a removed node or if the total
    /// overflows `usize`.
    pub fn path_kmer_length(&self, path: &[NodeId]) -> Option<usize> {
        path.iter().try_fold(0usize, |total, &node| {
            let node_len = self.node_weight(node)?.abs_ind.len();
            total.checked_add(node_len)
        })
    }

    /// Validate all directed connections in the graph.
    ///
    /// Every edge must have a reverse edge with the type returned by [`EdgeType::rev`], with
    /// multiplicity taken into account. Identical parallel edges are reported as duplicates.
    /// Self-loops are valid when their reverse-type multiplicity is present.
    ///
    /// The validator checks graph topology only; it cannot verify k-mer overlap because the graph
    /// stores hashes rather than the underlying sequences.
    pub fn validate(&self) -> Result<(), GraphValidationReport> {
        let mut connections: BTreeMap<(NodeId, NodeId, EdgeType), Vec<EdgeId>> = BTreeMap::new();

        for edge in self.inner.edge_references() {
            let from = from_backend_node(edge.source());
            let to = from_backend_node(edge.target());
            let edge_id = from_backend_edge(edge.id());
            connections
                .entry((from, to, edge.weight().t))
                .or_default()
                .push(edge_id);
        }

        let mut issues = Vec::new();
        for (&(from, to, edge_type), edge_ids) in &connections {
            if edge_ids.len() > 1 {
                issues.push(GraphValidationIssue::DuplicateEdge {
                    from,
                    to,
                    edge_type,
                    count: edge_ids.len(),
                });
            }

            let reverse_count = connections
                .get(&(to, from, edge_type.rev()))
                .map_or(0, Vec::len);
            for &edge_id in edge_ids.iter().skip(reverse_count) {
                issues.push(GraphValidationIssue::UnpairedEdge {
                    edge_id,
                    from,
                    to,
                    edge_type,
                });
            }
        }

        if issues.is_empty() {
            Ok(())
        } else {
            Err(GraphValidationReport { issues })
        }
    }
}

// ─── Traversal – bidirected aware ────────────────────────────────────────────

impl DbgGraph {
    /// Forward neighbours from a given canonicality origin.
    ///
    /// Replaces `out_neighbours_bi` / `out_neighbours_min` / `out_neighbours_max`.
    pub fn forward_neighbors(&self, n: NodeId, carry: CarryType) -> Vec<(NodeId, EdgeType)> {
        self.inner
            .edges_directed(to_backend_node(n), Outgoing)
            .filter(|e| e.weight().t.get_from_and_to().0 == carry)
            .map(|e| (from_backend_node(e.target()), e.weight().t))
            .collect()
    }

    /// Backward neighbours arriving at a given canonicality.
    ///
    /// Replaces `in_neighbours_bi` / `in_neighbours_min` / `in_neighbours_max`.
    pub fn backward_neighbors(&self, n: NodeId, carry: CarryType) -> Vec<(NodeId, EdgeType)> {
        self.inner
            .edges_directed(to_backend_node(n), Incoming)
            .filter(|e| e.weight().t.get_from_and_to().1 == carry)
            .map(|e| (from_backend_node(e.source()), e.weight().t))
            .collect()
    }

    /// First outgoing edge type for a node, in backend iteration order.
    ///
    /// The returned edge is whichever edge petgraph yields first; no semantic ordering is
    /// guaranteed.
    #[inline]
    pub fn first_outgoing_edge_type(&self, node: NodeId) -> Option<EdgeType> {
        self.inner
            .edges_directed(to_backend_node(node), Outgoing)
            .next()
            .map(|edge| edge.weight().t)
    }

    /// Outgoing edges whose source carry matches `carry`, including edge id.
    pub fn outgoing_edges_by_carry(
        &self,
        node: NodeId,
        carry: CarryType,
    ) -> Vec<(EdgeId, NodeId, EdgeType)> {
        self.inner
            .edges_directed(to_backend_node(node), Outgoing)
            .filter(|edge| edge.weight().t.get_from_and_to().0 == carry)
            .map(|edge| {
                (
                    from_backend_edge(edge.id()),
                    from_backend_node(edge.target()),
                    edge.weight().t,
                )
            })
            .collect()
    }

    /// Outgoing bubble candidates whose source carry matches `carry`.
    pub fn bubble_start_edges_by_carry(
        &self,
        node: NodeId,
        carry: CarryType,
    ) -> Vec<BubbleStartEdge> {
        self.inner
            .edges_directed(to_backend_node(node), Outgoing)
            .filter(|edge| edge.weight().t.get_from_and_to().0 == carry)
            .map(|edge| BubbleStartEdge {
                target: from_backend_node(edge.target()),
                edge_type: edge.weight().t,
            })
            .collect()
    }

    /// All outgoing edges as target and edge type.
    pub fn outgoing_edges(&self, node: NodeId) -> Vec<(NodeId, EdgeType)> {
        self.inner
            .edges_directed(to_backend_node(node), Outgoing)
            .map(|edge| (from_backend_node(edge.target()), edge.weight().t))
            .collect()
    }

    /// All incoming edges as source and edge type.
    pub fn incoming_edges(&self, node: NodeId) -> Vec<(NodeId, EdgeType)> {
        self.inner
            .edges_directed(to_backend_node(node), Incoming)
            .map(|edge| (from_backend_node(edge.source()), edge.weight().t))
            .collect()
    }

    /// All outgoing non-self-loop neighbours (carry-agnostic).
    ///
    /// Replaces `get_good_neighbours_bi`.
    pub fn all_neighbors(&self, n: NodeId) -> Vec<(NodeId, EdgeType)> {
        self.inner
            .edges_directed(to_backend_node(n), Outgoing)
            .filter(|e| e.source() != e.target())
            .map(|e| (from_backend_node(e.target()), e.weight().t))
            .collect()
    }

    /// Count of forward edges from a given canonicality.
    ///
    /// Replaces `out_degree_bi` / `out_degree_min` / `out_degree_max`.
    pub fn forward_degree(&self, n: NodeId, carry: CarryType) -> usize {
        self.inner
            .edges_directed(to_backend_node(n), Outgoing)
            .filter(|e| e.weight().t.get_from_and_to().0 == carry)
            .count()
    }

    /// Count of backward edges to a given canonicality.
    ///
    /// Replaces `in_degree_bi` / `in_degree_min` / `in_degree_max`.
    pub fn backward_degree(&self, n: NodeId, carry: CarryType) -> usize {
        self.inner
            .edges_directed(to_backend_node(n), Incoming)
            .filter(|e| e.weight().t.get_from_and_to().1 == carry)
            .count()
    }

    /// Total outgoing edge count.
    pub fn out_degree(&self, n: NodeId) -> usize {
        self.inner
            .neighbors_directed(to_backend_node(n), petgraph::EdgeDirection::Outgoing)
            .count()
    }

    /// Total incoming edge count.
    pub fn in_degree(&self, n: NodeId) -> usize {
        self.inner
            .neighbors_directed(to_backend_node(n), petgraph::EdgeDirection::Incoming)
            .count()
    }

    /// Count outgoing non-self-loop edges.
    ///
    /// This is the ordinary connection count used by path and collapse helpers. It is not the
    /// total edge count used by `ambiguous_nodes()`.
    ///
    /// Replaces `get_good_connections_degree`.
    pub fn nonself_degree(&self, n: NodeId) -> usize {
        self.inner
            .edges_directed(to_backend_node(n), Outgoing)
            .filter(|e| e.source() != e.target())
            .count()
    }

    // ── Old names kept as thin wrappers for minimal call-site churn ──────────

    /// Alias for `forward_neighbors`.
    #[inline]
    pub fn out_neighbours_bi(&self, n: NodeId, ty: CarryType) -> Vec<(NodeId, EdgeType)> {
        self.forward_neighbors(n, ty)
    }

    /// Alias for `forward_neighbors(n, CarryType::Min)`.
    #[inline]
    pub fn out_neighbours_min(&self, n: NodeId) -> Vec<(NodeId, EdgeType)> {
        self.forward_neighbors(n, CarryType::Min)
    }

    /// Alias for `forward_neighbors(n, CarryType::Max)`.
    #[inline]
    pub fn out_neighbours_max(&self, n: NodeId) -> Vec<(NodeId, EdgeType)> {
        self.forward_neighbors(n, CarryType::Max)
    }

    /// Alias for `backward_neighbors`.
    #[inline]
    pub fn in_neighbours_bi(&self, n: NodeId, ty: CarryType) -> Vec<(NodeId, EdgeType)> {
        self.backward_neighbors(n, ty)
    }

    /// Alias for `backward_neighbors(n, CarryType::Min)`.
    #[inline]
    pub fn in_neighbours_min(&self, n: NodeId) -> Vec<(NodeId, EdgeType)> {
        self.backward_neighbors(n, CarryType::Min)
    }

    /// Alias for `backward_neighbors(n, CarryType::Max)`.
    #[inline]
    pub fn in_neighbours_max(&self, n: NodeId) -> Vec<(NodeId, EdgeType)> {
        self.backward_neighbors(n, CarryType::Max)
    }

    /// Alias for `forward_degree`.
    #[inline]
    pub fn out_degree_bi(&self, n: NodeId, ty: CarryType) -> usize {
        self.forward_degree(n, ty)
    }

    /// Alias for `forward_degree(n, CarryType::Min)`.
    #[inline]
    pub fn out_degree_min(&self, n: NodeId) -> usize {
        self.forward_degree(n, CarryType::Min)
    }

    /// Alias for `forward_degree(n, CarryType::Max)`.
    #[inline]
    pub fn out_degree_max(&self, n: NodeId) -> usize {
        self.forward_degree(n, CarryType::Max)
    }

    /// Alias for `backward_degree`.
    #[inline]
    pub fn in_degree_bi(&self, n: NodeId, ty: CarryType) -> usize {
        self.backward_degree(n, ty)
    }

    /// Alias for `backward_degree(n, CarryType::Min)`.
    #[inline]
    pub fn in_degree_min(&self, n: NodeId) -> usize {
        self.backward_degree(n, CarryType::Min)
    }

    /// Alias for `backward_degree(n, CarryType::Max)`.
    #[inline]
    pub fn in_degree_max(&self, n: NodeId) -> usize {
        self.backward_degree(n, CarryType::Max)
    }

    /// Alias for `all_neighbors`.
    #[inline]
    pub fn get_good_neighbours_bi(&self, n: NodeId) -> Vec<(NodeId, EdgeType)> {
        self.all_neighbors(n)
    }

    /// Alias for `nonself_degree`.
    #[inline]
    pub fn get_good_connections_degree(&self, n: NodeId) -> usize {
        self.nonself_degree(n)
    }
}

// ─── Topology ────────────────────────────────────────────────────────────────

impl DbgGraph {
    /// Nodes that form junctions or path starts/ends in a bidirected graph.
    ///
    /// Nodes without self-loops retain the ordinary rule: exactly two outgoing non-self edges
    /// with one Min-origin edge are treated as a straight intermediate path; all other non-empty
    /// patterns are ambiguous.
    ///
    /// A node with a self-loop is ambiguous only when it also has at least one outgoing Min-origin
    /// edge to a different node. A self-loop by itself, or a self-loop with only non-Min outgoing
    /// edges, is not ambiguous. Self-loop edges remain excluded from ordinary neighbour lists.
    ///
    /// Replaces `get_ambiguous_nodes_bi`.
    pub fn ambiguous_nodes(&self) -> BTreeSet<NodeId> {
        self.inner
            .node_indices()
            .filter(|n| {
                let mut nonself_connections = 0usize;
                let mut nonself_min_connections = 0usize;
                let mut has_self_loop = false;

                for edge in self.inner.edges_directed(*n, Outgoing) {
                    if edge.source() == edge.target() {
                        has_self_loop = true;

                        // A self-loop is ambiguous only if a previously observed
                        // non-self Min edge already exists.
                        if nonself_min_connections > 0 {
                            return true;
                        }
                        continue;
                    }

                    nonself_connections += 1;

                    if edge.weight().t.get_from_and_to().0 == CarryType::Min {
                        nonself_min_connections += 1;

                        // A non-self Min edge plus a self-loop makes this node
                        // an ambiguity boundary, regardless of other edges.
                        if has_self_loop {
                            return true;
                        }
                    }
                }

                // A self-loop by itself, or with only non-Min non-self edges,
                // is not ambiguous under the current criterion.
                if has_self_loop {
                    return false;
                }

                if nonself_connections == 0 {
                    return false;
                }

                if nonself_connections == 2 {
                    return nonself_min_connections != 1;
                }

                true
            })
            .map(from_backend_node)
            .collect()
    }

    /// Entry-point nodes (all incoming edges share the same destination canonicality, or no incoming edges).
    ///
    /// Replaces `externals_bi`.
    pub fn externals(&self) -> Vec<NodeId> {
        self.inner
            .node_indices()
            .filter(|n| {
                let mut it = self.inner.edges_directed(*n, petgraph::Direction::Incoming);
                let ct: CarryType;
                if let Some(e) = it.next() {
                    ct = e.weight().t.get_from_and_to().1;
                } else {
                    return true;
                }
                for e in it {
                    if e.weight().t.get_from_and_to().1 != ct {
                        return false;
                    }
                }
                true
            })
            .map(from_backend_node)
            .collect()
    }

    /// Whether node `n` has any self-loop edge.
    ///
    /// Replaces `node_has_self_loops`.
    pub fn has_self_loop(&self, n: NodeId) -> bool {
        self.inner
            .edges_directed(to_backend_node(n), Outgoing)
            .any(|e| e.target() == to_backend_node(n))
    }

    /// Number of weakly connected components in the graph.
    ///
    /// `StableGraph` cannot be passed directly to petgraph's component algorithm because it does
    /// not implement `NodeCompactIndexable`: removed nodes leave holes in its index space. The
    /// conversion below compacts node indices before running the union-find algorithm.
    pub fn connected_components(&self) -> usize {
        petgraph_connected_components(&petgraph::graph::Graph::from(self.inner.clone()))
    }

    // ── Old names ────────────────────────────────────────────────────────────

    /// Alias for `ambiguous_nodes`.
    #[inline]
    pub fn get_ambiguous_nodes_bi(&self) -> BTreeSet<NodeId> {
        self.ambiguous_nodes()
    }

    /// Alias for `externals`.
    #[inline]
    pub fn externals_bi(&self) -> Vec<NodeId> {
        self.externals()
    }

    /// Alias for `has_self_loop`.
    #[inline]
    pub fn node_has_self_loops(&self, n: NodeId) -> bool {
        self.has_self_loop(n)
    }
}

// ─── Mutation ────────────────────────────────────────────────────────────────

impl DbgGraph {
    /// Remove a node and all its incident edges, returning the node data.
    pub fn remove_node(&mut self, n: NodeId) -> Option<NodeStruct> {
        self.inner.remove_node(to_backend_node(n))
    }

    /// Add a single directed edge.
    pub fn add_edge(&mut self, from: NodeId, to: NodeId, edge: EdgeType) {
        self.inner.add_edge(
            to_backend_node(from),
            to_backend_node(to),
            EmptyEdge { t: edge },
        );
    }

    /// Add both directions of a bidirected edge pair.
    ///
    /// This method does not check whether an edge pair already connects these nodes; callers must
    /// avoid duplicate insertion.
    pub fn add_bi_edge(&mut self, from: NodeId, to: NodeId, edge: EdgeType) {
        self.inner.add_edge(
            to_backend_node(from),
            to_backend_node(to),
            EmptyEdge { t: edge },
        );
        self.inner.add_edge(
            to_backend_node(to),
            to_backend_node(from),
            EmptyEdge { t: edge.rev() },
        );
    }

    /// Remove a specific edge by index.
    pub fn remove_edge(&mut self, e: EdgeId) {
        self.inner.remove_edge(to_backend_edge(e));
    }

    /// Remove all self-loop edges.
    pub fn remove_self_loops(&mut self) {
        self.inner.retain_edges(|g, e| {
            let (n1, n2) = g.edge_endpoints(e).unwrap();
            n1 != n2
        });
    }

    /// Remove all nodes whose count is below `min`.
    pub fn retain_nodes_by_count(&mut self, min: u32) {
        self.inner
            .retain_nodes(|g, n| g.node_weight(n).unwrap().counts >= min);
    }

    /// All edge indices connecting `from` → `to` (regardless of type).
    pub fn edges_between(&self, from: NodeId, to: NodeId) -> Vec<EdgeId> {
        self.inner
            .edges_connecting(to_backend_node(from), to_backend_node(to))
            .map(|e| from_backend_edge(e.id()))
            .collect()
    }

    /// Remove all incoming and outgoing edges of node `n`.
    pub fn remove_all_edges_of(&mut self, n: NodeId) {
        let edges: Vec<EdgeId> = self
            .inner
            .edges_directed(to_backend_node(n), Outgoing)
            .chain(self.inner.edges_directed(to_backend_node(n), Incoming))
            .map(|e| from_backend_edge(e.id()))
            .collect();
        for e in edges {
            self.inner.remove_edge(to_backend_edge(e));
        }
    }

    /// Set an edge type by edge index.
    #[inline]
    pub fn set_edge_type(&mut self, edge: EdgeId, edge_type: EdgeType) {
        self.inner.edge_weight_mut(to_backend_edge(edge)).unwrap().t = edge_type;
    }

    /// Set the first edge type between two nodes, in backend iteration order.
    ///
    /// The edge is selected by its endpoints only; `edge_type` is the replacement
    /// type, not a filter. Callers must ensure that a connecting edge exists and
    /// that no competing parallel edge makes the selection ambiguous.
    ///
    /// # Panics
    ///
    /// Panics if no edge connects `from` to `to`.
    #[inline]
    pub fn set_first_edge_type_between(&mut self, from: NodeId, to: NodeId, edge_type: EdgeType) {
        let edge = self
            .inner
            .edges_connecting(to_backend_node(from), to_backend_node(to))
            .next()
            .expect("set_first_edge_type_between requires an existing connecting edge")
            .id();
        self.inner.edge_weight_mut(edge).unwrap().t = edge_type;
    }

    /// Find the incoming edge and modify edge orientations after shrinking through a non-direct internal edge.
    ///
    /// The caller must ensure that `internal_edge_ty` is non-direct (`MinToMax` or `MaxToMin`),
    /// that `in_edge_ty` ends at the source carry of `internal_edge_ty`, and that exactly one
    /// directed edge of type `in_edge_ty` exists from `prev_node` to `base_node`.
    /// `base_node` and `prev_node` must be distinct. This method does not use `in_edge_ty` to
    /// disambiguate multiple edges between the same endpoints.
    ///
    /// The shrinker establishes these conditions by filtering the incoming edge carry and checking
    /// single-edge path degrees before calling this method.
    ///
    /// # Panics
    ///
    /// Panics if the preconditions are violated, including when `base_node` and `prev_node` are
    /// the same node, when the connecting edge is absent, or when multiple connecting edges are
    /// present.
    #[inline]
    pub fn modify_edges_when_shrinking_between(
        &mut self,
        base_node: NodeId,
        prev_node: NodeId,
        internal_edge_ty: EdgeType,
        in_edge_ty: EdgeType,
    ) {
        assert!(
            base_node != prev_node,
            "modify_edges_when_shrinking_between requires distinct base_node and prev_node; got {base_node:?}"
        );

        let edges = self.edges_between(prev_node, base_node);
        if edges.len() > 1 {
            panic!("More than one linking outgoing edge, this should not happen unless there are multiple connections to the same node.");
        }
        let in_edge_ind = edges.first().copied().expect(
            "modify_edges_when_shrinking_between requires an existing edge from prev_node to base_node",
        );
        self.modify_edges_when_shrinking(
            base_node,
            prev_node,
            internal_edge_ty,
            in_edge_ind,
            in_edge_ty,
        );
    }

    /// Modify edge orientations after shrinking through a non-direct internal edge.
    #[inline]
    fn modify_edges_when_shrinking(
        &mut self,
        base_node: NodeId,
        prev_node: NodeId,
        internal_edge_ty: EdgeType,
        in_edge_ind: EdgeId,
        in_edge_ty: EdgeType,
    ) {
        log::trace!("base_node: {:?} prev_node: {:?} internal_edge_ty: {:?} in_edge_ind: {:?} in_edge_ty: {:?}",
            base_node, prev_node, internal_edge_ty, in_edge_ind, in_edge_ty,
        );
        log::trace!(
            "edges coming to the prev_node from the base_node: {:?}",
            self.edges_between(base_node, prev_node)
        );
        log::trace!(
            "edges coming to the base_node from the prev_node: {:?}",
            self.edges_between(prev_node, base_node)
        );

        match internal_edge_ty {
            EdgeType::MinToMax => {
                match in_edge_ty {
                    EdgeType::MinToMin => {
                        self.set_edge_type(in_edge_ind, EdgeType::MinToMax);
                        self.set_first_edge_type_between(base_node, prev_node, EdgeType::MinToMax);
                    }
                    EdgeType::MaxToMin => {
                        self.set_edge_type(in_edge_ind, EdgeType::MaxToMax);
                        self.set_first_edge_type_between(base_node, prev_node, EdgeType::MinToMin);
                    }
                    _ => panic!("Not expected edge type"),
                }
                self.node_weight_mut(base_node)
                    .unwrap()
                    .set_internal_edge(EdgeType::MaxToMax);
            }
            EdgeType::MaxToMin => {
                match in_edge_ty {
                    EdgeType::MaxToMax => {
                        self.set_edge_type(in_edge_ind, EdgeType::MaxToMin);
                        self.set_first_edge_type_between(base_node, prev_node, EdgeType::MaxToMin);
                    }
                    EdgeType::MinToMax => {
                        self.set_edge_type(in_edge_ind, EdgeType::MinToMin);
                        self.set_first_edge_type_between(base_node, prev_node, EdgeType::MaxToMax);
                    }
                    _ => panic!("Not expected edge type"),
                }
                self.node_weight_mut(base_node)
                    .unwrap()
                    .set_internal_edge(EdgeType::MinToMin);
            }
            _ => panic!("Value not expected"),
        }
    }
}

// ─── Merge primitive ─────────────────────────────────────────────────────────

impl DbgGraph {
    /// Merge `child` into `parent`, remove `child` from the graph, and return the child's data.
    ///
    /// `parent` and `child` must be distinct, existing nodes. Any required edge rewiring must be
    /// completed before calling this method, because removing `child` also removes its incident
    /// edges.
    ///
    /// The caller is responsible for edge rewiring and finalising `set_mean_counts`,
    /// `set_internal_edge`, and `invert_if_needed` on `parent`.
    ///
    /// # Panics
    ///
    /// Panics if either node ID is invalid or if `parent` and `child` refer to the same node.
    pub fn merge_nodes(&mut self, parent: NodeId, child: NodeId, edge: EdgeType) -> NodeStruct {
        let child_data = self.inner.remove_node(to_backend_node(child)).unwrap();
        self.inner
            .node_weight_mut(to_backend_node(parent))
            .unwrap()
            .merge(&child_data, edge);
        child_data
    }
}

// ─── Output ──────────────────────────────────────────────────────────────────

impl DbgGraph {
    /// Write graph to a DOT format writer.
    pub fn write_to_dot<W: Write>(&self, f: &mut W) {
        let output = self.get_dot_string();
        f.write_all(output.as_bytes())
            .unwrap_or_else(|error| panic!("failed to write DOT graph: {error}"));
        f.flush()
            .unwrap_or_else(|error| panic!("failed to flush DOT graph: {error}"));
    }

    /// Return graph as a DOT format string.
    pub fn get_dot_string(&self) -> String {
        let mut representatives = self.representative_edge_indices();
        representatives.sort_unstable();

        let mut graphfordot = self.inner.clone();
        graphfordot.retain_edges(|_, edge| representatives.binary_search(&edge).is_ok());
        format!(
            "{:?}",
            Dot::with_attr_getters(
                &graphfordot,
                &[Config::NodeNoLabel, Config::EdgeNoLabel],
                &|_, e| {
                    let (source_carry, target_carry) = e.weight().t.get_from_and_to();

                    let source_sign = match source_carry {
                        CarryType::Min => "+",
                        CarryType::Max => "-",
                    };
                    let target_sign = match target_carry {
                        CarryType::Min => "+",
                        CarryType::Max => "-",
                    };

                    format!(
                        "label = \"{:?} ({source_sign},{target_sign})\"",
                        e.weight().t
                    )
                },
                &|_, n| format!(
                    "label = \"{:?} | {:?}\" counts = {:?} kmers = {:?}",
                    n.1.counts,
                    n.1.abs_ind.len(),
                    n.1.counts,
                    n.1.abs_ind.len()
                ),
            )
        )
    }

    /// Write graph to a GFAv1 format writer.
    pub fn write_to_gfa<W: Write>(&self, f: &mut W) {
        let output = self.get_gfa_string();
        f.write_all(output.as_bytes())
            .unwrap_or_else(|error| panic!("failed to write GFA1 graph: {error}"));
        f.flush()
            .unwrap_or_else(|error| panic!("failed to flush GFA1 graph: {error}"));
    }

    /// Return graph as a GFAv1 string.
    pub fn get_gfa_string(&self) -> String {
        let mut output = "H\tVN:Z:1.0\n".to_owned();

        // Given the duality of the graph, to do the exportation of it as GFA, it is enough to assume that e.g. the canonical hashes represent the direct strand.
        // With that, the resulting nodes and edges will represent, by construction, correctly both strands.

        // Add nodes.
        //
        // `LN` is the segment LENGTH and `KC` the k-mer COUNT: these two arguments used to be the other
        // way round, so every tool reading these files (Bandage, …) saw the coverage as the length and
        // vice versa. `get_gfa2_string` below has always been right.
        self.inner.node_indices().for_each(|ni| {
            let tmpw = self.inner.node_weight(ni).unwrap();
            output.push_str(&format!(
                "S\t{}\t*\tLN:i:{}\tKC:i:{}\n",
                ni.index(),
                self.k + tmpw.abs_ind.len() - 1,
                tmpw.counts
            ));
        });

        // Add one record for each reciprocal edge class. The orientation signs
        // preserve which end of each oriented segment the link attaches to.
        self.representative_edge_indices()
            .into_iter()
            .for_each(|ei| {
                let (sid, tid) = self.inner.edge_endpoints(ei).unwrap();
                let (st, tt) = self.inner.edge_weight(ei).unwrap().t.get_from_and_to();

                let ssign: &str = match st {
                    CarryType::Min => "+",
                    CarryType::Max => "-",
                };
                let tsign: &str = match tt {
                    CarryType::Min => "+",
                    CarryType::Max => "-",
                };

                output.push_str(&format!(
                    "L\t{}\t{}\t{}\t{}\t{}M\tID:Z:{}\n",
                    sid.index(),
                    ssign,
                    tid.index(),
                    tsign,
                    self.k - 1,
                    ei.index(),
                ));
            });

        output
    }

    /// Write graph to a GFAv2 format writer.
    pub fn write_to_gfa2<W: Write>(&self, f: &mut W) {
        let output = self.get_gfa2_string();
        f.write_all(output.as_bytes())
            .unwrap_or_else(|error| panic!("failed to write GFA2 graph: {error}"));
        f.flush()
            .unwrap_or_else(|error| panic!("failed to flush GFA2 graph: {error}"));
    }

    /// Return graph as a GFAv2 string.
    pub fn get_gfa2_string(&self) -> String {
        let mut output = "H\tVN:Z:2.0\n".to_owned();

        // Given the duality of the graph, to do the exportation of it as GFA, it is enough to assume that e.g. the canonical hashes represent the direct strand.
        // With that, the resulting nodes and edges will represent, by construction, correctly both strands.

        // Add nodes
        self.inner.node_indices().for_each(|ni| {
            output.push_str(&format!(
                "S\t{}\t{}\t*\n",
                ni.index(),
                self.k + self.inner.node_weight(ni).unwrap().abs_ind.len() - 1
            ));
        });

        // Add one record for each reciprocal edge class. The orientation signs
        // and endpoint intervals preserve the attachment ends in GFA2.
        self.representative_edge_indices()
            .into_iter()
            .for_each(|ei| {
                let (sid, tid) = self.inner.edge_endpoints(ei).unwrap();
                let (st, tt) = self.inner.edge_weight(ei).unwrap().t.get_from_and_to();

                let ssign: &str;
                let tsign: &str;
                let sbeg: String;
                let send: String;
                let tbeg: String;
                let tend: String;

                match st {
                    CarryType::Min => {
                        ssign = "+";
                        let tmplen = self.inner.node_weight(sid).unwrap().abs_ind.len();
                        sbeg = format!("{}", tmplen);
                        send = format!("{}$", self.k + tmplen - 1);
                    }
                    CarryType::Max => {
                        ssign = "-";
                        sbeg = "0".to_owned();
                        send = format!("{}", self.k - 1);
                    }
                }

                match tt {
                    CarryType::Min => {
                        tsign = "+";
                        tbeg = "0".to_owned();
                        tend = format!("{}", self.k - 1);
                    }
                    CarryType::Max => {
                        tsign = "-";
                        let tmplen = self.inner.node_weight(tid).unwrap().abs_ind.len();
                        tbeg = format!("{}", tmplen);
                        tend = format!("{}$", self.k + tmplen - 1);
                    }
                }

                output.push_str(&format!(
                    "E\t{}\t{}{}\t{}{}\t{}\t{}\t{}\t{}\t{}M\n",
                    ei.index(),
                    sid.index(),
                    ssign,
                    tid.index(),
                    tsign,
                    sbeg,
                    send,
                    tbeg,
                    tend,
                    self.k - 1,
                ));
            });

        output
    }

    fn representative_edge_indices(&self) -> Vec<BackendEdgeIndex> {
        let mut seen: BTreeSet<ExportEdgeKey> = BTreeSet::new();
        let mut representatives = Vec::new();

        for edge in self.inner.edge_references() {
            let key = (
                from_backend_node(edge.source()),
                edge.weight().t,
                from_backend_node(edge.target()),
            );
            let reverse_key = (key.2, key.1.rev(), key.0);
            let canonical_key = if key <= reverse_key { key } else { reverse_key };

            if seen.insert(canonical_key) {
                representatives.push(edge.id());
            }
        }

        representatives
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use nohash_hasher::NoHashHasher;
    use std::collections::HashMap;
    use std::hash::BuildHasherDefault;

    #[test]
    fn test_new_graph() {
        let graph = DbgGraph::new(31);
        assert_eq!(graph.k(), 31);
        assert_eq!(graph.node_count(), 0);
        assert_eq!(graph.edge_count(), 0);
    }

    #[test]
    fn test_add_node() {
        let mut graph = DbgGraph::new(31);
        let node_data = NodeStruct {
            counts: 1,
            abs_ind: vec![0],
            innerdir: None,
        };
        let node_idx = graph.add_node(node_data.clone());

        assert_eq!(graph.node_count(), 1);
        assert!(graph.contains_node(node_idx));
        assert_eq!(graph.node_weight(node_idx), Some(&node_data));
    }

    #[test]
    fn test_path_kmer_length_sums_live_nodes() {
        let mut graph = DbgGraph::new(31);
        let first = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![0; 3],
            innerdir: None,
        });
        let second = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![0; 7],
            innerdir: None,
        });

        assert_eq!(graph.path_kmer_length(&[first, second]), Some(10));
    }

    #[test]
    fn test_path_kmer_length_rejects_removed_nodes() {
        let mut graph = DbgGraph::new(31);
        let node = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![0],
            innerdir: None,
        });
        graph.remove_node(node);

        assert_eq!(graph.path_kmer_length(&[node]), None);
    }

    #[test]
    fn test_default_graph_uses_k_one_and_exports_empty_nodes() {
        let mut graph = DbgGraph::default();
        let node1 = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![],
            innerdir: None,
        });
        let node2 = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![],
            innerdir: None,
        });
        graph.add_bi_edge(node1, node2, EdgeType::MinToMin);

        assert_eq!(graph.k(), 1);
        assert!(graph.get_gfa_string().contains("LN:i:0"));
        assert!(graph.get_gfa_string().contains("\t0M"));
        assert!(graph.get_gfa2_string().contains("\t0$"));
    }

    #[test]
    #[should_panic(expected = "k-mer length must be at least 1")]
    fn test_new_rejects_zero_k() {
        let _ = DbgGraph::new(0);
    }

    #[test]
    #[should_panic(expected = "k-mer length must be at least 1")]
    fn test_from_kmer_map_rejects_zero_k() {
        let map: HashMap<u64, HashInfoSimple, BuildHasherDefault<NoHashHasher<u64>>> =
            HashMap::default();
        let _ = DbgGraph::from_kmer_map(0, &map);
    }

    #[test]
    fn test_add_edge() {
        let mut graph = DbgGraph::new(31);
        let node1 = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![0],
            innerdir: None,
        });
        let node2 = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![1],
            innerdir: None,
        });

        graph.add_edge(node1, node2, EdgeType::MinToMin);

        assert_eq!(graph.edge_count(), 1);
        assert_eq!(graph.out_degree(node1), 1);
        assert_eq!(graph.out_degree(node2), 0);
    }

    #[test]
    fn test_graph_from_empty_kmer_map() {
        let k = 31;
        let empty_map: HashMap<u64, HashInfoSimple, BuildHasherDefault<NoHashHasher<u64>>> =
            HashMap::default();

        let graph = DbgGraph::from_kmer_map(k, &empty_map);

        assert_eq!(graph.k(), k);
        assert_eq!(graph.node_count(), 0);
        assert_eq!(graph.edge_count(), 0);
    }

    #[test]
    fn test_graph_from_indexed_kmers_consumes_aligned_edges() {
        let graph = DbgGraph::from_indexed_kmers(
            31,
            vec![10, 20],
            vec![7, 11],
            vec![
                vec![
                    IndexedEdge {
                        target: 1,
                        edge_type: EdgeType::MaxToMax,
                    },
                    IndexedEdge {
                        target: 1,
                        edge_type: EdgeType::MinToMin,
                    },
                ],
                vec![
                    IndexedEdge {
                        target: 0,
                        edge_type: EdgeType::MinToMin,
                    },
                    IndexedEdge {
                        target: 0,
                        edge_type: EdgeType::MaxToMax,
                    },
                ],
            ],
            vec![1, 1],
        );

        assert_eq!(graph.k(), 31);
        assert_eq!(graph.node_count(), 2);
        assert_eq!(graph.edge_count(), 2);
        assert!(graph.validate().is_ok());
    }

    #[test]
    #[should_panic(expected = "outside the indexed k-mer table")]
    fn test_graph_from_indexed_kmers_rejects_invalid_target() {
        let _ = DbgGraph::from_indexed_kmers(
            31,
            vec![10],
            vec![7],
            vec![vec![IndexedEdge {
                target: 1,
                edge_type: EdgeType::MinToMin,
            }]],
            vec![0],
        );
    }

    #[test]
    #[should_panic(expected = "predecessor count exceeds neighbour-list length")]
    fn test_graph_from_indexed_kmers_rejects_invalid_split() {
        let _ = DbgGraph::from_indexed_kmers(31, vec![10], vec![7], vec![Vec::new()], vec![1]);
    }

    #[test]
    fn test_node_degrees() {
        let mut graph = DbgGraph::new(31);
        let node1 = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![0],
            innerdir: None,
        });
        let node2 = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![1],
            innerdir: None,
        });
        let node3 = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![2],
            innerdir: None,
        });

        graph.add_edge(node1, node2, EdgeType::MinToMin);
        graph.add_edge(node2, node1, EdgeType::MaxToMax);
        graph.add_edge(node2, node3, EdgeType::MinToMin);

        assert_eq!(graph.out_degree(node1), 1);
        assert_eq!(graph.in_degree(node1), 1);
        assert_eq!(graph.out_degree(node2), 2); // node2 has edges to node1 and node3
        assert_eq!(graph.in_degree(node2), 1);
        assert_eq!(graph.out_degree(node3), 0);
        assert_eq!(graph.in_degree(node3), 1);
    }

    #[test]
    fn test_forward_backward_neighbors() {
        let mut graph = DbgGraph::new(31);
        let node1 = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![0],
            innerdir: None,
        });
        let node2 = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![1],
            innerdir: None,
        });

        graph.add_edge(node1, node2, EdgeType::MinToMin);
        graph.add_edge(node2, node1, EdgeType::MaxToMax);

        let forward_neighbors = graph.forward_neighbors(node1, CarryType::Min);
        let backward_neighbors = graph.backward_neighbors(node1, CarryType::Max);

        assert_eq!(forward_neighbors.len(), 1);
        assert_eq!(backward_neighbors.len(), 1);
    }

    #[test]
    fn test_gfa_serialization_empty_graph() {
        let graph = DbgGraph::new(31);
        let gfa_string = graph.get_gfa_string();

        // Check that it's valid GFA format
        assert!(gfa_string.contains("H\tVN:Z:1.0"));
        // Empty graphs may not have segment lines, so just check it's not empty
        assert!(!gfa_string.is_empty());
    }

    #[test]
    fn test_graph_contains_node() {
        let mut graph = DbgGraph::new(31);
        let node_data = NodeStruct {
            counts: 1,
            abs_ind: vec![0],
            innerdir: None,
        };
        let node_idx = graph.add_node(node_data);

        assert!(graph.contains_node(node_idx));
        graph.remove_node(node_idx);
        assert!(!graph.contains_node(node_idx));
    }

    #[test]
    fn test_edge_types() {
        let mut graph = DbgGraph::new(31);
        let node1 = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![0],
            innerdir: None,
        });
        let node2 = graph.add_node(NodeStruct {
            counts: 1,
            abs_ind: vec![1],
            innerdir: None,
        });

        graph.add_edge(node1, node2, EdgeType::MinToMin);
        graph.add_edge(node2, node1, EdgeType::MaxToMax);

        let edges = graph.all_neighbors(node1);
        assert_eq!(edges.len(), 1);
        assert_eq!(edges[0].1, EdgeType::MinToMin);

        let edges = graph.all_neighbors(node2);
        assert_eq!(edges.len(), 1);
        assert_eq!(edges[0].1, EdgeType::MaxToMax);
    }
}
