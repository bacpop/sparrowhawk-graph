//! Core types for the sparrowhawk de Bruijn graph.

use std::collections::HashMap;
use std::hash::BuildHasherDefault;

use nohash_hasher::NoHashHasher;

/// Index type for both nodes and edges in the graph.
pub type Idx = usize;

/// Type for representing the weight (count) of a k-mer.
pub type EdgeWeight = u32;

/// A graph edge whose target is an index into an aligned k-mer table.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IndexedEdge {
    pub target: Idx,
    pub edge_type: EdgeType,
}

/// Opaque, copyable handle identifying a node in the graph.
///
/// The handle is stable for the lifetime of its node, can be stored and passed by value, and is
/// meaningful only for the graph in which it was created while that node exists.
#[repr(transparent)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeId(pub(crate) Idx);

/// Opaque, copyable handle identifying an edge in the graph.
///
/// The handle is stable for the lifetime of its edge, can be stored and passed by value, and is
/// meaningful only for the graph in which it was created while that edge exists.
#[repr(transparent)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdgeId(pub(crate) Idx);

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
    pub counts: u32,
}

/// Enum that describes the type of one edge of the graph (essentially,
/// from which hash it comes (either canonical/minimum or non-canonical/maximum)
/// and with what it is linked (again, either min/max).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum EdgeType {
    /// Links canonical hash to canonical hash
    MinToMin,
    /// Links non-canonical hash to non-canonical hash
    MaxToMax,
    /// Links canonical hash to non-canonical hash
    MinToMax,
    /// Links non-canonical hash to canonical hash
    MaxToMin,
}

impl EdgeType {
    /// Reverses the EdgeType, changing it to the DNA de Bruijn-graph edge that would exist
    /// and begin at the end of the original EdgeType and end at the beginning of the same one.
    pub fn rev(&self) -> EdgeType {
        match self {
            EdgeType::MinToMin => EdgeType::MaxToMax,
            EdgeType::MaxToMax => EdgeType::MinToMin,
            EdgeType::MinToMax | EdgeType::MaxToMin => *self,
        }
    }

    /// Returns the CarryType of the origin and end of the edges, i.e. the canonicality of the
    /// hashes that this edge connects.
    pub fn get_from_and_to(&self) -> (CarryType, CarryType) {
        match self {
            EdgeType::MinToMin => (CarryType::Min, CarryType::Min),
            EdgeType::MaxToMax => (CarryType::Max, CarryType::Max),
            EdgeType::MaxToMin => (CarryType::Max, CarryType::Min),
            EdgeType::MinToMax => (CarryType::Min, CarryType::Max),
        }
    }

    /// Constructs an EdgeType from the canonicality (CarryType) of the source and end of the edges.
    pub fn from_carrytypes(first: CarryType, second: CarryType) -> EdgeType {
        match (first, second) {
            (CarryType::Min, CarryType::Min) => EdgeType::MinToMin,
            (CarryType::Max, CarryType::Max) => EdgeType::MaxToMax,
            (CarryType::Min, CarryType::Max) => EdgeType::MinToMax,
            (CarryType::Max, CarryType::Min) => EdgeType::MaxToMin,
        }
    }

    /// Answers whether the EdgeType is of the direct types (i.e. linkes a canonical hash with another canonical one,
    /// or the equivalent with non-canonicals).
    pub fn is_direct(&self) -> bool {
        matches!(self, EdgeType::MinToMin | EdgeType::MaxToMax)
    }
}

/// This enum describes the canonicality of a hash
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CarryType {
    /// Canonical hash (the minimum hash of the pair)
    Min,
    /// Non-canonical hash (the maximum hash of the pair)
    Max,
}
