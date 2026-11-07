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

## G2: stopping, loss, ancestors and copying

* Yule starts with one lineage at time zero, including its first waiting time
  (root stem). At n tips its expected height is H_(n-1)/birth. Crown-conditioned
  height has a different expectation and is not silently substituted.
* Exact-n generation ends at the first hitting time. n plus a deadline means
  reaching n **by** that time, with rejection of extinct/late attempts. Time-only
  generation conditions on >=2 surviving lineages. Limits bound failed attempts
  and explosive populations. Standard and fast names share one correct kernel.
* Polytomies add Poisson(lambda) offspring to a binary birth. At the final birth
  offspring are capped to the remaining tip budget: the exact-n model therefore
  has a truncated final polytomy. Rejection of overshoots was an alternative, but
  biases the entire trajectory rather than just its stated stopping transition.
* Sampling is Bernoulli(rho) per extant tip, after evolution. Fossil/manuscript
  observations are a Poisson process per lineage-time, split into a continuation
  and a zero-length sampled tip. Loss later in time does not erase earlier samples.
* Copying irregularity uses independent mean-one lognormal birth multipliers per
  new lineage. A birth can have an additional Poisson burst with an explicit
  probability. This models simultaneous extra copies, not undocumented calendar
  periods of scriptorium activity. High loss is a preset parameter, not an
  assumption that changes the process.
