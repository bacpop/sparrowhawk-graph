//! Node and edge weight structs for the de Bruijn graph.

use core::fmt;

use crate::types::EdgeType;

/// Contains the information stored at a graph node.
#[derive(Clone, Debug)]
pub struct NodeStruct {
    /// k-mer count (or mean count after shrinkage).
    pub counts: u16,

    /// Ordered list of canonical k-mer hashes in this node (one per k-mer for unshrunk nodes).
    pub abs_ind: Vec<u64>,

    /// Internal edge direction for shrunk nodes (the orientation used when the node was created).
    pub innerdir: Option<EdgeType>,
}

impl NodeStruct {
    /// Merges `other` node's k-mer list into `self`, respecting orientation.
    ///
    /// `tytoother` is the edge type connecting `self` → `other`.
    pub fn merge(&mut self, other: &NodeStruct, tytoother: EdgeType) {
        if let Some(thisid) = self.innerdir {
            let thisct = thisid.get_from_and_to().0;

            if thisct == tytoother.get_from_and_to().0 {
                if let Some(otherid) = other.innerdir {
                    let otherct = otherid.get_from_and_to().0;
                    if tytoother.get_from_and_to().1 == otherct {
                        self.abs_ind.extend(&other.abs_ind);
                    } else {
                        let mut newvec = other.abs_ind.clone();
                        newvec.reverse();
                        self.abs_ind.extend(&newvec);
                    }
                } else {
                    self.abs_ind.extend(&other.abs_ind);
                }
            } else {
                if let Some(otherid) = other.innerdir {
                    let otherct = otherid.get_from_and_to().0;
                    if tytoother.get_from_and_to().1 == otherct {
                        let mut tmpvec = other.abs_ind.clone();
                        tmpvec.reverse();
                        self.abs_ind.splice(0..0, tmpvec.iter().cloned());
                    } else {
                        self.abs_ind.splice(0..0, other.abs_ind.iter().cloned());
                    }
                } else {
                    self.abs_ind.splice(0..0, other.abs_ind.iter().cloned());
                }
            }
        } else {
            if let Some(otherid) = other.innerdir {
                let otherct = otherid.get_from_and_to().0;
                if tytoother.get_from_and_to().1 == otherct {
                    self.abs_ind.extend(&other.abs_ind);
                } else {
                    let mut newvec = other.abs_ind.clone();
                    newvec.reverse();
                    self.abs_ind.extend(&newvec);
                }
            } else {
                self.abs_ind.extend(&other.abs_ind);
            }
        }
    }

    /// Reverses `abs_ind` and the inner edge if the current orientation does not match `outedge`.
    pub fn invert_if_needed(&mut self, outedge: EdgeType) {
        if let Some(id) = self.innerdir {
            if id.get_from_and_to().1 != outedge.get_from_and_to().0 {
                self.abs_ind.reverse();
                self.innerdir.unwrap().rev();
            }
        }
    }

    /// Sets `counts` to the rounded mean of `countsvec`.
    pub fn set_mean_counts(&mut self, countsvec: &[u16]) {
        self.counts = (countsvec.iter().map(|&e| e as u32).sum::<u32>() as f32
            / countsvec.len() as f32)
            .round() as u16;
    }

    /// Sets the internal edge type (must be `MinToMin` or `MaxToMax`).
    pub fn set_internal_edge(&mut self, ed: EdgeType) {
        match ed {
            EdgeType::MinToMin | EdgeType::MaxToMax => self.innerdir = Some(ed),
            _ => panic!("Non-valid internal edge!"),
        }
    }
}

impl fmt::Display for NodeStruct {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.counts)
    }
}

/// Edge weight in the graph — holds only the `EdgeType`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EmptyEdge {
    /// Type of the edge.
    pub t: EdgeType,
}

impl fmt::Display for EmptyEdge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "")
    }
}
