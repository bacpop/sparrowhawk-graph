//! `sparrowhawk-graph`: bidirected de Bruijn graph library for the Sparrowhawk assembler.

pub mod graph;
pub mod node;
pub mod types;

pub use graph::DbgGraph;
pub use node::{EmptyEdge, NodeStruct};
pub use types::{
    CarryType, EdgeIndex, EdgeType, EdgeWeight, HashInfoSimple, Idx, KmerMap, NodeIndex,
};

/// Serialized representation of a single contig (ordered list of node data).
pub type SerializedContig = Vec<NodeStruct>;

/// Collection of serialized contigs.
pub type SerializedContigs = Vec<Vec<NodeStruct>>;
