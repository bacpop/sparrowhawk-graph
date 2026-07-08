//! `sparrowhawk-graph`: bidirected de Bruijn graph library for the Sparrowhawk assembler.

pub mod graph;
pub mod node;
pub mod types;

pub use graph::{BubbleStartEdge, DbgGraph};
pub use node::{EmptyEdge, NodeStruct};
pub use types::{CarryType, EdgeId, EdgeType, EdgeWeight, HashInfoSimple, Idx, KmerMap, NodeId};

/// Serialized representation of a single contig (ordered list of node data).
pub type SerializedContig = Vec<NodeStruct>;

/// Collection of serialized contigs.
pub type SerializedContigs = Vec<Vec<NodeStruct>>;

/// Return the number of k-mers represented by an ordered list of graph nodes.
pub fn get_nodelist_kmer_length(nodes: &[NodeStruct]) -> usize {
    nodes.iter().map(|node| node.abs_ind.len()).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_node(len: usize) -> NodeStruct {
        NodeStruct {
            counts: 1,
            abs_ind: vec![0u64; len],
            innerdir: None,
        }
    }

    #[test]
    fn get_nodelist_kmer_length_empty() {
        assert_eq!(get_nodelist_kmer_length(&[]), 0);
    }

    #[test]
    fn get_nodelist_kmer_length_single_node_no_kmers() {
        assert_eq!(get_nodelist_kmer_length(&[make_node(0)]), 0);
    }

    #[test]
    fn get_nodelist_kmer_length_single_node_five() {
        assert_eq!(get_nodelist_kmer_length(&[make_node(5)]), 5);
    }

    #[test]
    fn get_nodelist_kmer_length_multiple_nodes() {
        let nodes = vec![make_node(3), make_node(0), make_node(7)];
        assert_eq!(get_nodelist_kmer_length(&nodes), 10);
    }

    #[test]
    fn get_nodelist_kmer_length_large() {
        let nodes: Vec<_> = (0..100).map(|_| make_node(50)).collect();
        assert_eq!(get_nodelist_kmer_length(&nodes), 5000);
    }
}
