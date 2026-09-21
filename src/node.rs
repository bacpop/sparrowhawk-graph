//! Node and edge weight structs for the de Bruijn graph.

use core::fmt;

use crate::types::EdgeType;

/// Structure that contains the information of a node.
#[derive(Clone, Debug, PartialEq)]
pub struct NodeStruct {
    /// Value that reflects either the counts of one k-mer, or a
    /// proxy value for shrunk nodes.
    pub counts: u32,

    /// List of hashes os k-mers
    pub abs_ind: Vec<u64>,

    /// Inner edge for those edges result of a shrinkage
    pub innerdir: Option<EdgeType>,
}

impl NodeStruct {
    /// Allows merging two nodes and their inner information
    pub fn merge(&mut self, other: &NodeStruct, tytoother: EdgeType) {
        // First, we'll see if we have an internal edge already.
        if let Some(thisid) = self.innerdir {
            // This is a bit more complicated...
            let thisct = thisid.get_from_and_to().0;

            if thisct == tytoother.get_from_and_to().0 {
                // Great! We don't need to change nothing from this vector.
                // Let's see now the OTHER vector...
                if let Some(otherid) = other.innerdir {
                    // This is a bit more complicated...
                    let otherct = otherid.get_from_and_to().0;
                    if tytoother.get_from_and_to().1 == otherct {
                        // We can directly merge the vectors
                        self.abs_ind.extend(&other.abs_ind);
                    } else {
                        // This is a bit more complicated: we need to INVERT
                        // the other vector BEFORE joining them
                        let mut newvec = other.abs_ind.clone();
                        newvec.reverse();
                        self.abs_ind.extend(&newvec);
                    }
                } else {
                    // This is easy!
                    self.abs_ind.extend(&other.abs_ind);
                }
            } else {
                // This is not great, this node is currently in the opposite orientation to that
                // of the current merge (as defined by tytoother).
                // Thus, we need to insert whatever comes at the beginning of the vector

                if let Some(otherid) = other.innerdir {
                    // This is a bit more complicated...
                    let otherct = otherid.get_from_and_to().0;
                    if tytoother.get_from_and_to().1 == otherct {
                        // The other shrunk node is aligned with the edge that connects it with this one.
                        // As ours is not, we must prepend the REVERSED information of the other shrunk node.
                        let mut tmpvec = other.abs_ind.clone();
                        tmpvec.reverse();
                        self.abs_ind.splice(0..0, tmpvec.iter().cloned());
                    } else {
                        // This is a bit more complicated: we DON'T need to reverse
                        // the other vector BEFORE joining them, because they are already in the same order.
                        self.abs_ind.splice(0..0, other.abs_ind.iter().cloned());
                    }
                } else {
                    self.abs_ind.splice(0..0, other.abs_ind.iter().cloned());
                }
            }
        } else {
            // This is easier!
            // Let's see if the other node has an internal edge itself

            if let Some(otherid) = other.innerdir {
                // This is a bit more complicated...
                let otherct = otherid.get_from_and_to().0;
                if tytoother.get_from_and_to().1 == otherct {
                    // We can directly merge the vectors
                    self.abs_ind.extend(&other.abs_ind);
                } else {
                    // This is a bit more complicated: we need to INVERT
                    // the other vector BEFORE joining them
                    let mut newvec = other.abs_ind.clone();
                    newvec.reverse();
                    self.abs_ind.extend(&newvec);
                }
            } else {
                // This is very easy! We just have to join the vectors
                self.abs_ind.extend(&other.abs_ind);
            }
        }
    }

    /// Checks whether it is needed to reverse the list of hashes and the inner edge.
    pub fn invert_if_needed(&mut self, outedge: EdgeType) {
        if let Some(id) = self.innerdir {
            if id.get_from_and_to().1 != outedge.get_from_and_to().0 {
                self.abs_ind.reverse();
                self.innerdir = Some(id.rev());
            }
        }
    }

    /// Sets the counts of the node as the mean coverage of everything merged into it, each
    /// `(counts, k-mers represented)` contribution weighted by its second element. A one-k-mer node
    /// and a thousand-k-mer node are not equal evidence.
    ///
    /// # Panics
    ///
    /// Panics if `contributions` is empty, or if they represent no k-mers at all.
    pub fn set_mean_counts(&mut self, contributions: &[(u32, usize)]) {
        assert!(
            !contributions.is_empty(),
            "set_mean_counts requires at least one count"
        );

        // u64 is ample: the sum is bounded by u32::MAX times the k-mers in the whole graph.
        let (weighted, kmers) = contributions.iter().fold((0u64, 0u64), |(w, n), &(c, k)| {
            (w + u64::from(c) * k as u64, n + k as u64)
        });
        assert!(kmers > 0, "set_mean_counts requires at least one k-mer");

        // Rounds half up. A weighted mean of u32s cannot exceed u32::MAX, so the cast is lossless.
        self.counts = ((weighted + kmers / 2) / kmers) as u32;
    }

    /// Sets the type of the internal edge as the one you provide the function.
    ///
    /// Internal edges must be direct: only `MinToMin` and `MaxToMax` are valid.
    pub fn set_internal_edge(&mut self, ed: EdgeType) {
        match ed {
            EdgeType::MinToMin | EdgeType::MaxToMax => self.innerdir = Some(ed),
            _ => panic!(
                "invalid internal edge {ed:?}: internal edges must be direct (MinToMin or MaxToMax)"
            ),
        }
    }
}

impl fmt::Display for NodeStruct {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.counts)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn node() -> NodeStruct {
        NodeStruct {
            counts: 7,
            abs_ind: vec![1],
            innerdir: None,
        }
    }

    #[test]
    #[should_panic(expected = "set_mean_counts requires at least one count")]
    fn set_mean_counts_panics_on_empty_input() {
        let mut node = node();
        node.set_mean_counts(&[]);
    }

    #[test]
    fn set_mean_counts_preserves_nonempty_behaviour() {
        let mut node = node();
        // Equal weights, so this is still the plain mean it always was.
        node.set_mean_counts(&[(2, 1), (4, 1)]);
        assert_eq!(node.counts, 3);
    }

    #[test]
    #[should_panic(expected = "set_mean_counts requires at least one k-mer")]
    fn set_mean_counts_panics_when_nothing_is_represented() {
        let mut node = node();
        node.set_mean_counts(&[(10, 0)]);
    }

    /// The whole point of the weighting: a single leftover k-mer must not drag a long unitig's
    /// coverage halfway down to its own.
    #[test]
    fn set_mean_counts_weights_by_represented_kmers() {
        let mut node = node();
        node.set_mean_counts(&[(10, 1), (100, 9)]);
        assert_eq!(node.counts, 91, "an unweighted mean would give 55");
    }

    /// Shrinking reaches the same unitig by different merge orders, so the answer must not depend
    /// on how the contributions were grouped. An unweighted mean fails this.
    #[test]
    fn set_mean_counts_is_associative_across_regroupings() {
        let (a, b, c) = ((10u32, 1usize), (40, 3), (100, 6));

        let mut all_at_once = node();
        all_at_once.set_mean_counts(&[a, b, c]);

        // First {a, b} into one node of 4 k-mers, then that node with c.
        let mut ab = node();
        ab.set_mean_counts(&[a, b]);
        let mut in_stages = node();
        in_stages.set_mean_counts(&[(ab.counts, a.1 + b.1), c]);

        assert_eq!(all_at_once.counts, in_stages.counts);
    }

    #[test]
    #[should_panic(expected = "invalid internal edge MinToMax: internal edges must be direct")]
    fn set_internal_edge_rejects_non_direct_edges() {
        let mut node = node();
        node.set_internal_edge(EdgeType::MinToMax);
    }

    #[test]
    fn invert_if_needed_reverses_sequence_and_direct_internal_edge() {
        let mut node = NodeStruct {
            counts: 7,
            abs_ind: vec![1, 2, 3],
            innerdir: Some(EdgeType::MinToMin),
        };

        node.invert_if_needed(EdgeType::MaxToMin);

        assert_eq!(node.abs_ind, vec![3, 2, 1]);
        assert_eq!(node.innerdir, Some(EdgeType::MaxToMax));
    }

    #[test]
    fn invert_if_needed_reverses_sequence_and_self_reversing_internal_edge() {
        let mut node = NodeStruct {
            counts: 7,
            abs_ind: vec![1, 2, 3],
            innerdir: Some(EdgeType::MinToMax),
        };

        node.invert_if_needed(EdgeType::MinToMin);

        assert_eq!(node.abs_ind, vec![3, 2, 1]);
        assert_eq!(node.innerdir, Some(EdgeType::MinToMax));
    }

    #[test]
    fn invert_if_needed_keeps_aligned_orientation_unchanged() {
        let mut node = NodeStruct {
            counts: 7,
            abs_ind: vec![1, 2, 3],
            innerdir: Some(EdgeType::MinToMin),
        };

        node.invert_if_needed(EdgeType::MinToMax);

        assert_eq!(node.abs_ind, vec![1, 2, 3]);
        assert_eq!(node.innerdir, Some(EdgeType::MinToMin));
    }
}

/// This struct contains the information of an edge in the graph, which is "empty" because it only contains the type of the edge
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
