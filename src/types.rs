//! Core types for the sparrowhawk de Bruijn graph.

use std::collections::HashMap;
use std::hash::BuildHasherDefault;

use nohash_hasher::NoHashHasher;

/// Index type for both nodes and edges in the graph.
pub type Idx = usize;

/// Type for representing the weight (count) of a k-mer.
pub type EdgeWeight = u16;

/// Type denoting the index of a node in the graph.
pub type NodeIndex = petgraph::stable_graph::NodeIndex<Idx>;

/// Type denoting the index of an edge in the graph.
pub type EdgeIndex = petgraph::stable_graph::EdgeIndex<Idx>;

/// Type alias for the canonical k-mer hash map produced by preprocessing.
pub type KmerMap = HashMap<u64, HashInfoSimple, BuildHasherDefault<NoHashHasher<u64>>>;

/// Struct that contains the basic information for one k-mer.
pub struct HashInfoSimple {
    /// Non-canonical (maximum) hash of this k-mer pair.
    pub hnc: u64,
    /// First and last bases (packed: top 2 bits = last base of fwd, bottom 2 = first base of fwd).
    pub b: u8,
    /// Neighbours found, if any, previous to this k-mer.
    pub pre: Vec<(u64, EdgeType)>,
    /// Neighbours found, if any, posterior to this k-mer.
    pub post: Vec<(u64, EdgeType)>,
    /// Count of this k-mer.
    pub counts: u16,
}

/// Describes the type of an edge: which canonicality it originates from and which it points to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeType {
    /// Links canonical hash to canonical hash.
    MinToMin,
    /// Links non-canonical hash to non-canonical hash.
    MaxToMax,
    /// Links canonical hash to non-canonical hash.
    MinToMax,
    /// Links non-canonical hash to canonical hash.
    MaxToMin,
}

impl EdgeType {
    /// Returns the reverse edge type.
    pub fn rev(&self) -> EdgeType {
        match self {
            EdgeType::MinToMin => EdgeType::MaxToMax,
            EdgeType::MaxToMax => EdgeType::MinToMin,
            EdgeType::MinToMax | EdgeType::MaxToMin => *self,
        }
    }

    /// Returns the `CarryType` of the origin and destination of this edge.
    pub fn get_from_and_to(&self) -> (CarryType, CarryType) {
        match self {
            EdgeType::MinToMin => (CarryType::Min, CarryType::Min),
            EdgeType::MaxToMax => (CarryType::Max, CarryType::Max),
            EdgeType::MaxToMin => (CarryType::Max, CarryType::Min),
            EdgeType::MinToMax => (CarryType::Min, CarryType::Max),
        }
    }

    /// Constructs an `EdgeType` from two `CarryType` values.
    pub fn from_carrytypes(first: CarryType, second: CarryType) -> EdgeType {
        match (first, second) {
            (CarryType::Min, CarryType::Min) => EdgeType::MinToMin,
            (CarryType::Max, CarryType::Max) => EdgeType::MaxToMax,
            (CarryType::Min, CarryType::Max) => EdgeType::MinToMax,
            (CarryType::Max, CarryType::Min) => EdgeType::MaxToMin,
        }
    }

    /// Whether this edge connects hashes of the same canonicality (direct orientation).
    pub fn is_direct(&self) -> bool {
        matches!(self, EdgeType::MinToMin | EdgeType::MaxToMax)
    }
}

/// Canonicality of a hash.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CarryType {
    /// Canonical hash (the minimum hash of the canonical/non-canonical pair).
    Min,
    /// Non-canonical hash (the maximum hash of the pair).
    Max,
}
