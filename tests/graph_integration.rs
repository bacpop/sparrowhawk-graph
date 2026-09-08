// Integration tests for sparrowhawk-graph
// Tests the public API and module interactions

use nohash_hasher::NoHashHasher;
use sparrowhawk_graph::*;
use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use std::io::{self, Write};

struct TestWriter {
    bytes: Vec<u8>,
    max_write: usize,
    write_error: bool,
    flush_error: bool,
}

impl TestWriter {
    fn chunked(max_write: usize) -> Self {
        Self {
            bytes: Vec::new(),
            max_write,
            write_error: false,
            flush_error: false,
        }
    }

    fn failing_write() -> Self {
        Self {
            bytes: Vec::new(),
            max_write: 1,
            write_error: true,
            flush_error: false,
        }
    }

    fn failing_flush() -> Self {
        Self {
            bytes: Vec::new(),
            max_write: usize::MAX,
            write_error: false,
            flush_error: true,
        }
    }
}

impl Write for TestWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if self.write_error {
            return Err(io::Error::other("test write failure"));
        }
        let n = buf.len().min(self.max_write);
        self.bytes.extend_from_slice(&buf[..n]);
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.flush_error {
            Err(io::Error::other("test flush failure"))
        } else {
            Ok(())
        }
    }
}

#[test]
fn test_complete_graph_workflow() {
    // Test a complete workflow: create graph, add nodes/edges, serialize
    let mut graph = DbgGraph::new(31);

    // Add some nodes
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

    // Add edges to form a simple path
    graph.add_edge(node1, node2, EdgeType::MinToMin);
    graph.add_edge(node2, node3, EdgeType::MinToMin);

    // Verify graph structure
    assert_eq!(graph.node_count(), 3);
    assert_eq!(graph.edge_count(), 2);
    assert_eq!(graph.out_degree(node1), 1);
    assert_eq!(graph.in_degree(node3), 1);

    // Test GFA serialization
    let gfa_string = graph.get_gfa_string();
    assert!(gfa_string.contains("H\tVN:Z:1.0"));
    assert!(gfa_string.contains("S\t1\t*"));
}

#[test]
fn test_graph_from_kmer_map_integration() {
    // Test creating a graph from a k-mer map (simulated)
    let k = 21; // Smaller k for testing

    // Create a minimal k-mer map
    let mut kmer_map = HashMap::with_hasher(BuildHasherDefault::<NoHashHasher<u64>>::default());

    // Add some test k-mers (simplified - in reality these would be proper k-mer hashes)
    kmer_map.insert(
        12345,
        HashInfoSimple {
            hnc: 12345,
            b: 0,
            pre: vec![],
            post: vec![],
            counts: 1,
        },
    );
    kmer_map.insert(
        23456,
        HashInfoSimple {
            hnc: 23456,
            b: 0,
            pre: vec![],
            post: vec![],
            counts: 1,
        },
    );

    let graph = DbgGraph::from_kmer_map(k, &kmer_map);

    // Basic verification
    assert_eq!(graph.k(), k);
    // Note: The actual node count depends on how k-mers connect
    // This is just a basic integration test
}

#[test]
fn test_public_api_access() {
    // Test that all public API functions are accessible
    let graph = DbgGraph::new(31);

    // These should all compile and run without panicking
    assert_eq!(graph.k(), 31);
    assert_eq!(graph.node_count(), 0);
    assert_eq!(graph.edge_count(), 0);
    let gfa_string = graph.get_gfa_string();
    assert!(gfa_string.contains("H\tVN:Z:1.0"));

    // Test that public types are accessible
    let _node_struct = NodeStruct {
        counts: 1,
        abs_ind: vec![0],
        innerdir: None,
    };
    let _edge_type = EdgeType::MinToMin;
}

#[test]
fn test_graph_operations_sequence() {
    // Test a sequence of operations that might be used in real scenarios
    let mut graph = DbgGraph::new(25);

    // Build a small graph
    let nodes: Vec<NodeId> = (0..5)
        .map(|i| {
            graph.add_node(NodeStruct {
                counts: 1,
                abs_ind: vec![i as u64],
                innerdir: None,
            })
        })
        .collect();

    // Connect nodes in a chain
    for i in 0..4 {
        graph.add_edge(nodes[i], nodes[i + 1], EdgeType::MinToMin);
    }

    // Verify the chain structure
    for i in 0..4 {
        assert_eq!(graph.out_degree(nodes[i]), 1);
        assert_eq!(graph.in_degree(nodes[i + 1]), 1);
    }

    // Test neighbor queries
    let neighbors = graph.all_neighbors(nodes[2]);
    assert_eq!(neighbors.len(), 1); // Should have one outgoing neighbor
    assert_eq!(neighbors[0].0, nodes[3]); // Should be node 3
}

#[test]
fn test_empty_graph_operations() {
    // Test that operations on empty graphs don't panic
    let graph = DbgGraph::new(31);

    let gfa_string = graph.get_gfa_string();
    assert!(gfa_string.contains("H\tVN:Z:1.0"));
    assert_eq!(graph.node_count(), 0);
    assert_eq!(graph.edge_count(), 0);

    // Test GFA serialization on empty graph
    let gfa_string = graph.get_gfa_string();
    assert!(gfa_string.contains("H\tVN:Z:1.0"));
}

fn graph_test_node(counts: u32, abs_ind: u64) -> NodeStruct {
    NodeStruct {
        counts,
        abs_ind: vec![abs_ind],
        innerdir: None,
    }
}

#[test]
fn test_exporters_complete_partial_writes() {
    let graph = DbgGraph::new(31);

    let mut dot = TestWriter::chunked(2);
    graph.write_to_dot(&mut dot);
    assert_eq!(dot.bytes, graph.get_dot_string().as_bytes());

    let mut gfa = TestWriter::chunked(2);
    graph.write_to_gfa(&mut gfa);
    assert_eq!(gfa.bytes, graph.get_gfa_string().as_bytes());

    let mut gfa2 = TestWriter::chunked(2);
    graph.write_to_gfa2(&mut gfa2);
    assert_eq!(gfa2.bytes, graph.get_gfa2_string().as_bytes());
}

#[test]
#[should_panic(expected = "failed to write DOT graph")]
fn test_dot_export_reports_write_errors() {
    let graph = DbgGraph::new(31);
    graph.write_to_dot(&mut TestWriter::failing_write());
}

#[test]
#[should_panic(expected = "failed to write GFA1 graph")]
fn test_gfa_export_reports_write_errors() {
    let graph = DbgGraph::new(31);
    graph.write_to_gfa(&mut TestWriter::failing_write());
}

#[test]
#[should_panic(expected = "failed to write GFA2 graph")]
fn test_gfa2_export_reports_write_errors() {
    let graph = DbgGraph::new(31);
    graph.write_to_gfa2(&mut TestWriter::failing_write());
}

#[test]
#[should_panic(expected = "failed to flush GFA1 graph")]
fn test_gfa_export_reports_flush_errors() {
    let graph = DbgGraph::new(31);
    graph.write_to_gfa(&mut TestWriter::failing_flush());
}

#[test]
fn test_isolated_node_count() {
    let mut graph = DbgGraph::new(31);
    let node1 = graph.add_node(graph_test_node(1, 0));
    let node2 = graph.add_node(graph_test_node(1, 1));
    let node3 = graph.add_node(graph_test_node(1, 2));

    assert_eq!(graph.isolated_node_count(), 3);
    graph.add_edge(node1, node2, EdgeType::MinToMin);
    assert_eq!(graph.isolated_node_count(), 1);
    graph.add_edge(node3, node3, EdgeType::MinToMin);
    assert_eq!(graph.isolated_node_count(), 0);
}

#[test]
fn test_first_outgoing_edge_type() {
    let mut graph = DbgGraph::new(31);
    let node1 = graph.add_node(graph_test_node(1, 0));
    let node2 = graph.add_node(graph_test_node(1, 1));

    assert_eq!(graph.first_outgoing_edge_type(node1), None);
    graph.add_edge(node1, node2, EdgeType::MaxToMin);
    assert_eq!(
        graph.first_outgoing_edge_type(node1),
        Some(EdgeType::MaxToMin)
    );
}

#[test]
fn test_outgoing_edges_by_carry() {
    let mut graph = DbgGraph::new(31);
    let source = graph.add_node(graph_test_node(1, 0));
    let min_target = graph.add_node(graph_test_node(1, 1));
    let max_target = graph.add_node(graph_test_node(1, 2));

    graph.add_edge(source, min_target, EdgeType::MinToMax);
    graph.add_edge(source, max_target, EdgeType::MaxToMin);

    let min_edges = graph.outgoing_edges_by_carry(source, CarryType::Min);
    assert_eq!(min_edges.len(), 1);
    assert_eq!(min_edges[0].1, min_target);
    assert_eq!(min_edges[0].2, EdgeType::MinToMax);

    let max_edges = graph.outgoing_edges_by_carry(source, CarryType::Max);
    assert_eq!(max_edges.len(), 1);
    assert_eq!(max_edges[0].1, max_target);
    assert_eq!(max_edges[0].2, EdgeType::MaxToMin);
}

#[test]
fn test_set_edge_type() {
    let mut graph = DbgGraph::new(31);
    let source = graph.add_node(graph_test_node(1, 0));
    let target = graph.add_node(graph_test_node(1, 1));
    graph.add_edge(source, target, EdgeType::MinToMin);
    let edge = graph.edges_between(source, target)[0];

    graph.set_edge_type(edge, EdgeType::MaxToMin);
    assert_eq!(graph.edge_weight(edge).unwrap().t, EdgeType::MaxToMin);
}

#[test]
fn test_set_first_edge_type_between() {
    let mut graph = DbgGraph::new(31);
    let source = graph.add_node(graph_test_node(1, 0));
    let target = graph.add_node(graph_test_node(1, 1));
    graph.add_edge(source, target, EdgeType::MinToMin);
    graph.add_edge(source, target, EdgeType::MaxToMax);
    let first_edge = graph.edges_between(source, target)[0];

    graph.set_first_edge_type_between(source, target, EdgeType::MinToMax);
    assert_eq!(graph.edge_weight(first_edge).unwrap().t, EdgeType::MinToMax);
}

#[test]
fn test_modify_edges_when_shrinking_min_to_max_with_min_to_min_input() {
    let mut graph = DbgGraph::new(31);
    let prev = graph.add_node(graph_test_node(1, 0));
    let base = graph.add_node(graph_test_node(1, 1));

    graph.add_edge(prev, base, EdgeType::MinToMin);
    graph.add_edge(base, prev, EdgeType::MaxToMax);
    let outgoing = graph.edges_between(base, prev)[0];

    graph.modify_edges_when_shrinking_between(base, prev, EdgeType::MinToMax, EdgeType::MinToMin);

    assert_eq!(
        graph
            .edge_weight(graph.edges_between(prev, base)[0])
            .unwrap()
            .t,
        EdgeType::MinToMax
    );
    assert_eq!(graph.edge_weight(outgoing).unwrap().t, EdgeType::MinToMax);
    assert_eq!(
        graph.node_weight(base).unwrap().innerdir,
        Some(EdgeType::MaxToMax)
    );
}
