// Integration tests for sphk-graph
// Tests the public API and module interactions

use sphk_graph::*;
use std::collections::HashMap;
use nohash_hasher::NoHashHasher;
use std::hash::BuildHasherDefault;

#[test]
fn test_complete_graph_workflow() {
    // Test a complete workflow: create graph, add nodes/edges, serialize
    let mut graph = DbgGraph::new(31);
    
    // Add some nodes
    let node1 = graph.add_node(NodeStruct { counts: 1, abs_ind: vec![0], innerdir: None });
    let node2 = graph.add_node(NodeStruct { counts: 1, abs_ind: vec![1], innerdir: None });
    let node3 = graph.add_node(NodeStruct { counts: 1, abs_ind: vec![2], innerdir: None });
    
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
    kmer_map.insert(12345, HashInfoSimple { 
        hnc: 12345, 
        b: 0, 
        pre: vec![], 
        post: vec![], 
        counts: 1 
    });
    kmer_map.insert(23456, HashInfoSimple { 
        hnc: 23456, 
        b: 0, 
        pre: vec![], 
        post: vec![], 
        counts: 1 
    });
    
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
    let _node_struct = NodeStruct { counts: 1, abs_ind: vec![0], innerdir: None };
    let _edge_type = EdgeType::MinToMin;
}

#[test]
fn test_graph_operations_sequence() {
    // Test a sequence of operations that might be used in real scenarios
    let mut graph = DbgGraph::new(25);
    
    // Build a small graph
    let nodes: Vec<NodeIndex> = (0..5).map(|i| {
        graph.add_node(NodeStruct { counts: 1, abs_ind: vec![i as u64], innerdir: None })
    }).collect();
    
    // Connect nodes in a chain
    for i in 0..4 {
        graph.add_edge(nodes[i], nodes[i+1], EdgeType::MinToMin);
    }
    
    // Verify the chain structure
    for i in 0..4 {
        assert_eq!(graph.out_degree(nodes[i]), 1);
        assert_eq!(graph.in_degree(nodes[i+1]), 1);
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
    
    // These should return empty results rather than panic
    let neighbors = graph.all_neighbors(NodeIndex::new(0));
    assert!(neighbors.is_empty());
    
    // Test GFA serialization on empty graph
    let gfa_string = graph.get_gfa_string();
    assert!(gfa_string.contains("H\tVN:Z:1.0"));
}