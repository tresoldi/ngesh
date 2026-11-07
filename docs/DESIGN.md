# ngesh 2 design log

## Audit of 1.2.1 (baseline 866b900)

Verified in `random_tree.py`: a single gamma-CDF threshold per edge replaces a
continuous-time process; `mut_exp` depends on column index; `k_hgt_vec` uses `k`
instead of `k_hgt`; donor selection uses older endpoints (including ancestors)
and the arbitrary min+max distance inversion; fast generation rejects exact tip
counts. `simulate_bad_sampling` only deletes tips, with no observation layer.
`common.set_seeds` mutates two global generators. The brief's claim that this
code hashes arbitrary seeds with Python `hash()` was **not** reproduced: strings
and floats use SHA-256, while other values are passed to NumPy. ETE3, NumPy and
six are mandatory dependencies; packaging targets Python 3.6. We replace this
core rather than preserve misleading 1.x parameter semantics.

## G1: reproducibility and representation

* SHA-256 uses type tags and length framing, avoiding ambiguous concatenations,
  randomized Python hashes and collisions between string, bytes and integer seeds.
* Every simulation component owns a private Random instance. Only its `random()`
  primitive is used, rather than relying on version-dependent convenience methods.
  Component streams and replicate keys are derived, never obtained from globals.
* Exponential inversion, Box–Muller normal, Marsaglia–Tsang gamma, gamma-ratio
  beta and additive Poisson arrivals are transparent, documented algorithms.
  Discrete gamma categories use conditional means, rather than quantile midpoints,
  so category rates average exactly one.
* Primary trees have explicit edges, absolute node times and arbitrary attributes.
  Edges are alive on [parent.time, node.time); zero-length observed tips are not
  biological lineages. Networks store secondary edges separately. Tree distances
  deliberately refer to the primary tree, not a shortest path through a network.
* Unary pruning sums edge lengths, retaining absolute sample times and named
  sampled ancestors. Root stem length is meaningful and is serialized.
