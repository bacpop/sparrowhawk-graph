//! `DbgGraph`: bidirected de Bruijn graph built on top of petgraph's `StableGraph`.

use std::collections::BTreeSet;
use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use std::io::Write;

use nohash_hasher::NoHashHasher;
use petgraph::algo::connected_components as petgraph_connected_components;
use petgraph::dot::{Config, Dot};
use petgraph::visit::EdgeRef;
use petgraph::Direction::{Incoming, Outgoing};

use crate::node::{EmptyEdge, NodeStruct};
use crate::types::{CarryType, EdgeIndex, EdgeType, HashInfoSimple, Idx, NodeIndex};

/// Inner petgraph type alias.
type Inner = petgraph::stable_graph::StableGraph<NodeStruct, EmptyEdge, petgraph::Directed, Idx>;

/// Bidirected de Bruijn graph.
#[derive(Default)]
pub struct DbgGraph {
    inner: Inner,
    k: usize,
}

// ─── Construction ────────────────────────────────────────────────────────────

impl DbgGraph {
    /// Create an empty graph for the given k-mer length.
    pub fn new(k: usize) -> Self {
        DbgGraph {
            inner: Inner::default(),
            k,
        }
    }

    /// Build a de Bruijn graph from the k-mer map produced by preprocessing.
    ///
    /// `map` is a `HashMap<canonical_hash, HashInfoSimple, ...>` as returned by
    /// `preprocessing_standalone` / `preprocessing_wasm`.
    pub fn from_kmer_map(
        k: usize,
        map: &HashMap<u64, HashInfoSimple, BuildHasherDefault<NoHashHasher<u64>>>,
    ) -> Self {
        let mut g = DbgGraph {
            inner: Inner::with_capacity(map.len(), map.len() * 2),
            k,
        };

        let mut tmpdict: HashMap<u64, NodeIndex, BuildHasherDefault<NoHashHasher<u64>>> =
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

    /// Insert a bare node (used in tests and programmatic construction).
    pub fn add_node(&mut self, node: NodeStruct) -> NodeIndex {
        self.inner.add_node(node)
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
    pub fn contains_node(&self, n: NodeIndex) -> bool {
        self.inner.contains_node(n)
    }

    /// Iterator over all live node indices.
    pub fn node_indices(&self) -> impl Iterator<Item = NodeIndex> + '_ {
        self.inner.node_indices()
    }
}

// ─── Node / edge weight access ───────────────────────────────────────────────

impl DbgGraph {
    /// Immutable reference to node weight.
    pub fn node_weight(&self, n: NodeIndex) -> Option<&NodeStruct> {
        self.inner.node_weight(n)
    }

    /// Mutable reference to node weight.
    pub fn node_weight_mut(&mut self, n: NodeIndex) -> Option<&mut NodeStruct> {
        self.inner.node_weight_mut(n)
    }

    /// Immutable reference to edge weight.
    pub fn edge_weight(&self, e: EdgeIndex) -> Option<&EmptyEdge> {
        self.inner.edge_weight(e)
    }

    /// Mutable reference to edge weight.
    pub fn edge_weight_mut(&mut self, e: EdgeIndex) -> Option<&mut EmptyEdge> {
        self.inner.edge_weight_mut(e)
    }

    /// Endpoints of an edge.
    pub fn edge_endpoints(&self, e: EdgeIndex) -> Option<(NodeIndex, NodeIndex)> {
        self.inner.edge_endpoints(e)
    }
}

// ─── Traversal – bidirected aware ────────────────────────────────────────────

impl DbgGraph {
    /// Forward neighbours from a given canonicality origin.
    ///
    /// Replaces `out_neighbours_bi` / `out_neighbours_min` / `out_neighbours_max`.
    pub fn forward_neighbors(&self, n: NodeIndex, carry: CarryType) -> Vec<(NodeIndex, EdgeType)> {
        self.inner
            .edges_directed(n, Outgoing)
            .filter(|e| e.weight().t.get_from_and_to().0 == carry)
            .map(|e| (e.target(), e.weight().t))
            .collect()
    }

    /// Backward neighbours arriving at a given canonicality.
    ///
    /// Replaces `in_neighbours_bi` / `in_neighbours_min` / `in_neighbours_max`.
    pub fn backward_neighbors(&self, n: NodeIndex, carry: CarryType) -> Vec<(NodeIndex, EdgeType)> {
        self.inner
            .edges_directed(n, Incoming)
            .filter(|e| e.weight().t.get_from_and_to().1 == carry)
            .map(|e| (e.source(), e.weight().t))
            .collect()
    }

    /// All outgoing non-self-loop neighbours (carry-agnostic).
    ///
    /// Replaces `get_good_neighbours_bi`.
    pub fn all_neighbors(&self, n: NodeIndex) -> Vec<(NodeIndex, EdgeType)> {
        self.inner
            .edges_directed(n, Outgoing)
            .filter(|e| e.source() != e.target())
            .map(|e| (e.target(), e.weight().t))
            .collect()
    }

    /// Count of forward edges from a given canonicality.
    ///
    /// Replaces `out_degree_bi` / `out_degree_min` / `out_degree_max`.
    pub fn forward_degree(&self, n: NodeIndex, carry: CarryType) -> usize {
        self.inner
            .edges_directed(n, Outgoing)
            .filter(|e| e.weight().t.get_from_and_to().0 == carry)
            .count()
    }

    /// Count of backward edges to a given canonicality.
    ///
    /// Replaces `in_degree_bi` / `in_degree_min` / `in_degree_max`.
    pub fn backward_degree(&self, n: NodeIndex, carry: CarryType) -> usize {
        self.inner
            .edges_directed(n, Incoming)
            .filter(|e| e.weight().t.get_from_and_to().1 == carry)
            .count()
    }

    /// Total outgoing edge count.
    pub fn out_degree(&self, n: NodeIndex) -> usize {
        self.inner
            .neighbors_directed(n, petgraph::EdgeDirection::Outgoing)
            .count()
    }

    /// Total incoming edge count.
    pub fn in_degree(&self, n: NodeIndex) -> usize {
        self.inner
            .neighbors_directed(n, petgraph::EdgeDirection::Incoming)
            .count()
    }

    /// Count of non-self-loop outgoing edges.
    ///
    /// Replaces `get_good_connections_degree`.
    pub fn nonself_degree(&self, n: NodeIndex) -> usize {
        self.inner
            .edges_directed(n, Outgoing)
            .filter(|e| e.source() != e.target())
            .count()
    }

    // ── Old names kept as thin wrappers for minimal call-site churn ──────────

    /// Alias for `forward_neighbors`.
    #[inline]
    pub fn out_neighbours_bi(&self, n: NodeIndex, ty: CarryType) -> Vec<(NodeIndex, EdgeType)> {
        self.forward_neighbors(n, ty)
    }

    /// Alias for `forward_neighbors(n, CarryType::Min)`.
    #[inline]
    pub fn out_neighbours_min(&self, n: NodeIndex) -> Vec<(NodeIndex, EdgeType)> {
        self.forward_neighbors(n, CarryType::Min)
    }

    /// Alias for `forward_neighbors(n, CarryType::Max)`.
    #[inline]
    pub fn out_neighbours_max(&self, n: NodeIndex) -> Vec<(NodeIndex, EdgeType)> {
        self.forward_neighbors(n, CarryType::Max)
    }

    /// Alias for `backward_neighbors`.
    #[inline]
    pub fn in_neighbours_bi(&self, n: NodeIndex, ty: CarryType) -> Vec<(NodeIndex, EdgeType)> {
        self.backward_neighbors(n, ty)
    }

    /// Alias for `backward_neighbors(n, CarryType::Min)`.
    #[inline]
    pub fn in_neighbours_min(&self, n: NodeIndex) -> Vec<(NodeIndex, EdgeType)> {
        self.backward_neighbors(n, CarryType::Min)
    }

    /// Alias for `backward_neighbors(n, CarryType::Max)`.
    #[inline]
    pub fn in_neighbours_max(&self, n: NodeIndex) -> Vec<(NodeIndex, EdgeType)> {
        self.backward_neighbors(n, CarryType::Max)
    }

    /// Alias for `forward_degree`.
    #[inline]
    pub fn out_degree_bi(&self, n: NodeIndex, ty: CarryType) -> usize {
        self.forward_degree(n, ty)
    }

    /// Alias for `forward_degree(n, CarryType::Min)`.
    #[inline]
    pub fn out_degree_min(&self, n: NodeIndex) -> usize {
        self.forward_degree(n, CarryType::Min)
    }

    /// Alias for `forward_degree(n, CarryType::Max)`.
    #[inline]
    pub fn out_degree_max(&self, n: NodeIndex) -> usize {
        self.forward_degree(n, CarryType::Max)
    }

    /// Alias for `backward_degree`.
    #[inline]
    pub fn in_degree_bi(&self, n: NodeIndex, ty: CarryType) -> usize {
        self.backward_degree(n, ty)
    }

    /// Alias for `backward_degree(n, CarryType::Min)`.
    #[inline]
    pub fn in_degree_min(&self, n: NodeIndex) -> usize {
        self.backward_degree(n, CarryType::Min)
    }

    /// Alias for `backward_degree(n, CarryType::Max)`.
    #[inline]
    pub fn in_degree_max(&self, n: NodeIndex) -> usize {
        self.backward_degree(n, CarryType::Max)
    }

    /// Alias for `all_neighbors`.
    #[inline]
    pub fn get_good_neighbours_bi(&self, n: NodeIndex) -> Vec<(NodeIndex, EdgeType)> {
        self.all_neighbors(n)
    }

    /// Alias for `nonself_degree`.
    #[inline]
    pub fn get_good_connections_degree(&self, n: NodeIndex) -> usize {
        self.nonself_degree(n)
    }
}

// ─── Topology ────────────────────────────────────────────────────────────────

impl DbgGraph {
    /// Nodes that form junctions or path starts/ends in a bidirected graph.
    ///
    /// A node is *ambiguous* unless it has exactly two non-self-loop outgoing edges and both
    /// come from `Min`.
    ///
    /// Replaces `get_ambiguous_nodes_bi`.
    pub fn ambiguous_nodes(&self) -> BTreeSet<NodeIndex> {
        self.inner
            .node_indices()
            .filter(|n| {
                let conns = self.nonself_degree(*n);
                if conns == 0 {
                    return false;
                } else if conns == 2 {
                    // If conns == 2, we just need to avoid a perfect intermediate kmer in a sequence of them.
                    if self.out_degree_min(*n) == 1 {
                        return false;
                    }
                }
                true
            })
            .collect()
    }

    /// Entry-point nodes (all incoming edges share the same destination canonicality, or no incoming edges).
    ///
    /// Replaces `externals_bi`.
    pub fn externals(&self) -> Vec<NodeIndex> {
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
            .collect()
    }

    /// Whether node `n` has any self-loop edge.
    ///
    /// Replaces `node_has_self_loops`.
    pub fn has_self_loop(&self, n: NodeIndex) -> bool {
        self.inner
            .edges_directed(n, Outgoing)
            .any(|e| e.target() == n)
    }

    /// Number of weakly connected components in the graph.
    pub fn connected_components(&self) -> usize {
        petgraph_connected_components(&petgraph::graph::Graph::from(self.inner.clone()))
    }

    // ── Old names ────────────────────────────────────────────────────────────

    /// Alias for `ambiguous_nodes`.
    #[inline]
    pub fn get_ambiguous_nodes_bi(&self) -> BTreeSet<NodeIndex> {
        self.ambiguous_nodes()
    }

    /// Alias for `externals`.
    #[inline]
    pub fn externals_bi(&self) -> Vec<NodeIndex> {
        self.externals()
    }

    /// Alias for `has_self_loop`.
    #[inline]
    pub fn node_has_self_loops(&self, n: NodeIndex) -> bool {
        self.has_self_loop(n)
    }
}

// ─── Mutation ────────────────────────────────────────────────────────────────

impl DbgGraph {
    /// Remove a node and all its incident edges, returning the node data.
    pub fn remove_node(&mut self, n: NodeIndex) -> Option<NodeStruct> {
        self.inner.remove_node(n)
    }

    /// Add a single directed edge.
    pub fn add_edge(&mut self, from: NodeIndex, to: NodeIndex, edge: EdgeType) {
        self.inner.add_edge(from, to, EmptyEdge { t: edge });
    }

    /// Add both directions of a bidirected edge pair.
    pub fn add_bi_edge(&mut self, from: NodeIndex, to: NodeIndex, edge: EdgeType) {
        self.inner.add_edge(from, to, EmptyEdge { t: edge });
        self.inner.add_edge(to, from, EmptyEdge { t: edge.rev() });
    }

    /// Remove a specific edge by index.
    pub fn remove_edge(&mut self, e: EdgeIndex) {
        self.inner.remove_edge(e);
    }

    /// Remove all self-loop edges.
    pub fn remove_self_loops(&mut self) {
        self.inner.retain_edges(|g, e| {
            let (n1, n2) = g.edge_endpoints(e).unwrap();
            n1 != n2
        });
    }

    /// Remove all nodes whose count is below `min`.
    pub fn retain_nodes_by_count(&mut self, min: u16) {
        self.inner
            .retain_nodes(|g, n| g.node_weight(n).unwrap().counts >= min);
    }

    /// All edge indices connecting `from` → `to` (regardless of type).
    pub fn edges_between(&self, from: NodeIndex, to: NodeIndex) -> Vec<EdgeIndex> {
        self.inner
            .edges_connecting(from, to)
            .map(|e| e.id())
            .collect()
    }

    /// Remove all incoming and outgoing edges of node `n`.
    pub fn remove_all_edges_of(&mut self, n: NodeIndex) {
        let edges: Vec<EdgeIndex> = self
            .inner
            .edges_directed(n, Outgoing)
            .chain(self.inner.edges_directed(n, Incoming))
            .map(|e| e.id())
            .collect();
        for e in edges {
            self.inner.remove_edge(e);
        }
    }
}

// ─── Merge primitive ─────────────────────────────────────────────────────────

impl DbgGraph {
    /// Merge `child` into `parent`, remove `child` from the graph, and return the child's data.
    ///
    /// The caller is responsible for edge rewiring and finalising `set_mean_counts`,
    /// `set_internal_edge`, and `invert_if_needed` on `parent`.
    pub fn merge_nodes(
        &mut self,
        parent: NodeIndex,
        child: NodeIndex,
        edge: EdgeType,
    ) -> NodeStruct {
        let child_data = self.inner.remove_node(child).unwrap();
        self.inner
            .node_weight_mut(parent)
            .unwrap()
            .merge(&child_data, edge);
        child_data
    }
}

// ─── Inner graph access ──────────────────────────────────────────────────────

impl DbgGraph {
    /// Read-only access to the underlying petgraph `StableGraph`.
    ///
    /// Used by algorithms (e.g. `tarjan_scc`, `connected_components`) that operate directly on
    /// petgraph types and have no wrapper alternative.
    pub fn inner_graph(
        &self,
    ) -> &petgraph::stable_graph::StableGraph<NodeStruct, EmptyEdge, petgraph::Directed, Idx> {
        &self.inner
    }

    /// Mutable access to the underlying petgraph `StableGraph`.
    ///
    /// Use sparingly — prefer the higher-level methods above.
    pub fn inner_graph_mut(
        &mut self,
    ) -> &mut petgraph::stable_graph::StableGraph<NodeStruct, EmptyEdge, petgraph::Directed, Idx>
    {
        &mut self.inner
    }
}

// ─── Output ──────────────────────────────────────────────────────────────────

impl DbgGraph {
    /// Write graph to a DOT format writer.
    pub fn write_to_dot<W: Write>(&self, f: &mut W) {
        let output = self.get_dot_string();
        let _ = f.write(output.as_bytes());
    }

    /// Return graph as a DOT format string.
    pub fn get_dot_string(&self) -> String {
        let mut graphfordot = self.inner.clone();
        graphfordot.retain_nodes(|g, n| {
            g.neighbors_directed(n, petgraph::EdgeDirection::Outgoing)
                .count()
                != 0
                || g.neighbors_directed(n, petgraph::EdgeDirection::Incoming)
                    .count()
                    != 0
        });
        format!(
            "{:?}",
            Dot::with_attr_getters(
                &graphfordot,
                &[Config::NodeNoLabel, Config::EdgeNoLabel],
                &|_, e| format!("label = \"{:?}\"", e.weight().t),
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
        let _ = f.write(output.as_bytes());
    }

    /// Return graph as a GFAv1 string.
    pub fn get_gfa_string(&self) -> String {
        let mut output = "H\tVN:Z:1.0\n".to_owned();

        // Given the duality of the graph, to do the exportation of it as GFA, it is enough to assume that e.g. the canonical hashes represent the direct strand.
        // With that, the resulting nodes and edges will represent, by construction, correctly both strands.

        // Add nodes
        self.inner.node_indices().for_each(|ni| {
            let tmpw = self.inner.node_weight(ni).unwrap();
            output.push_str(&format!(
                "S\t{}\t*\tLN:i:{}\tKC:i:{}\n",
                ni.index(),
                tmpw.counts,
                self.k + tmpw.abs_ind.len() - 1
            ));
        });

        // Add edges
        self.inner.edge_indices().for_each(|ei| {
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
        let _ = f.write(output.as_bytes());
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

        // Add edges
        self.inner.edge_indices().for_each(|ei| {
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
        assert!(!graph.contains_node(NodeIndex::new(999)));
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
