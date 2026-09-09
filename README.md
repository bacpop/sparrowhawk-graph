# sparrowhawk-graph
Bidirected de Bruijn graph library for the [sparrowhawk](https://github.com/bacpop/sparrowhawk-asm) bacterial genome assembler, written in Rust.

---
## Disclaimer :warning: :construction:
This is a **work in progress** project. This in particular implies:

- The API is not yet stable, and might change without notice.
- Code might be messy, and not even documented.
- General documentation might be short or even missing.
- Finding unexpected errors/behaviour or bugs should not be a surprise.

These (and potentially other) items will be progressively fixed before version 1.0.

---

## sparrowhawk?
Sparrowhawk was at one time the Archmage of [Earthsea](https://en.wikipedia.org/wiki/Earthsea).
Also, the [sparrowhawk](https://en.wikipedia.org/wiki/Eurasian_sparrowhawk) (*Accipiter nisus*) is a bird of prey native to Europe (and the island of Gont).

# Description

**Note:** this repository contains the de Bruijn graph library used by [sparrowhawk-asm](https://github.com/bacpop/sparrowhawk-asm), the Rust-based genomic assembler. If you are looking for the assembler itself, see that repository; for its web implementation, see [sparrowhawk](https://github.com/bacpop/sparrowhawk).

sparrowhawk-graph provides the node-based, bidirected de Bruijn graph on which sparrowhawk assembles genomes. It abstracts the graph backend (currently [petgraph](https://docs.rs/petgraph/latest/petgraph)'s `StableGraph`) away from the assembler, so that it can be changed in the future without touching the assembler itself.

Current **main features**:
- A bidirected de Bruijn graph, `DbgGraph`: each node represents a canonical k-mer pair (a k-mer and its reverse complement), and each adjacency is stored as a pair of directed edges typed by the canonicity of their endpoints (`EdgeType`: `MinToMin`, `MaxToMax`, `MinToMax`, `MaxToMin`), so both strands are represented at once.
- Construction from the canonical-hash k-mer map produced by sparrowhawk's preprocessing (`DbgGraph::from_kmer_map`).
- Strand-aware traversal and topology queries: forward/backward neighbours and degrees by canonicity, ambiguous (junction) nodes, external nodes, self-loops, and connected components.
- The mutation primitives used by the assembler's graph-simplification stages: node merging for path shrinking, bidirected edge insertion, edge retyping, and node/edge removal.
- Graph exportation in [DOT](https://en.wikipedia.org/wiki/DOT_%28graph_description_language%29) and [GFA](https://gfa-spec.github.io/GFA-spec/) versions 1.1 and 2.
- Only three dependencies ([petgraph](https://docs.rs/petgraph/latest/petgraph), [nohash-hasher](https://docs.rs/nohash-hasher/latest/nohash_hasher/), and [log](https://docs.rs/log/latest/log/)), and compilation both to native and WebAssembly targets.

# Installation
This crate is not yet published on crates.io: it is meant to be used as a git dependency. To get an in principle working version (up to some degree), always pin a versioned tag, e.g.:

```toml
[dependencies]
sparrowhawk-graph = { git = "https://github.com/bacpop/sparrowhawk-graph.git", tag = "vX.Y.Z" }
```

To compile it from source you will need the [rust toolchain](https://www.rust-lang.org/tools/install) installed in your system. Development has been done only on x86_64 GNU/Linux-based systems, and most surely will probably stay that way (i.e. no other systems have been tested). Again, always clone a versioned tag, e.g.:

```
git clone --branch vX.Y.Z https://github.com/bacpop/sparrowhawk-graph.git
cd sparrowhawk-graph
cargo build --release
```

# Usage
sparrowhawk-graph is a library, so there is no binary to run: its reference consumer is [sparrowhawk-asm](https://github.com/bacpop/sparrowhawk-asm). A minimal example of programmatic use:

```rust
use sparrowhawk_graph::{DbgGraph, EdgeType, NodeStruct};

let mut graph = DbgGraph::new(31);
let n1 = graph.add_node(NodeStruct { counts: 5, abs_ind: vec![0xBEEF], innerdir: None });
let n2 = graph.add_node(NodeStruct { counts: 7, abs_ind: vec![0xCAFE], innerdir: None });

// One biological adjacency = two paired directed edges, one per strand.
graph.add_bi_edge(n1, n2, EdgeType::MinToMin);

assert_eq!(graph.out_degree(n1), 1);
println!("{}", graph.get_gfa_string());
```

You can run the unit and integration tests with `cargo test`, and build the API documentation with `cargo doc --open`.
