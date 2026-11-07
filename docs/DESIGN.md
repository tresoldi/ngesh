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

## G3 lexical / G6 pilot

* Replacements are exponential-waiting-time events per character/branch, so
  arbitrarily many events occur and every event time is recorded in history.tsv.
  Synonym death is another exponential clock, inherited across speciation.
  Each concept has root ID 0, recurrent pool IDs 1..K, and monotonically allocated
  novel IDs above K. Recurrent pool draws may be silent; the event is still logged.
* Gamma rates have mean one before invariant masking. Invariants are a separate
  mixture, so the unconditional mean is rate*(1-p_invariant). Relaxed clocks are
  independent mean-one lognormal edge multipliers. Column indices affect labels,
  not the distribution of rates.
* Long TSV deliberately cannot distinguish binary absence from missing: both
  omit rows under the consumer contract. Full observed.json and truth.json are
  authoritative lossless matrices (null=missing, []=known absence). NEXUS and
  PHYLIP retain this distinction. Every evolution event is in history.tsv;
  horizontal events will additionally be in events.tsv.
* Integer-token multistate NEXUS supports lexical classes beyond 62 symbols.
  Relaxed PHYLIP uses deterministic whitespace-free aliases with a name map.
* A manifest hashes all other files, never itself. Grids stage each replicate in
  a temporary directory and rename on completion; resume verifies hashes, version,
  full config, preset and seed. Conflicting or corrupt outputs fail visibly.
* Python 3.10 lacks tomllib. Rather than add a runtime dependency or drop 3.10,
  it uses a strict documented scalar/array/table TOML subset. Python >=3.11 uses
  tomllib. Example configurations intentionally stay within the common subset.
* Summary quartet scores use at most 1000 evenly spaced lexicographic quartets,
  record their sample size, and consume no random draws. Full quartets can be
  requested explicitly. This keeps 1000-replicate benchmarking practical.

## G4/G5: chronological processes and observation

* One chronological scheduler handles births/edge starts, state changes, synonym
  deaths and transfers. Independent and borrowing runs use the identical kernel;
  the horizontal stream is separate, so a zero rate is exactly identical.
  Borrowing proposals have intensity b per character per lineage-time. Proposals
  with no eligible donor are skipped, so recorded event intensity is b only while
  another lineage is alive. Silent transfers are retained in the log.
* Donor weights are exp(-d/scale), normalized after subtracting the minimum
  distance to prevent underflow. Tree distance uses the two lineages at event time.
  Geography is two independent Brownian coordinates with SD diffusion*sqrt(dt);
  chronological queries simulate increments, with children inheriting the parental
  position. Geography has its own stream and does not perturb character events.
* Synonym deadlines are inherited at splits. An incoming transferred reading is
  established (its donor's future synonym deadlines are not inherited). Replace
  resets the cell; add unions it. This is a chosen transmission model, rather than
  a claim that all borrowed synonyms retain identical lifetimes.
* Contamination occurs on newly copied branches, draws a contemporaneous second
  exemplar excluding siblings, and replaces either Bernoulli(f) scattered loci
  or one contiguous block of round(f*L) loci. A transfer network splits edges at
  event times; secondary edges have length zero. It is an event display, not a
  probabilistic network with inferred inheritance probabilities.
* Missing cells use null, true selection/omission uses the empty set. Selection
  masks are inherited and can only shrink, explicitly generating nested item
  collections. Coverage is one Beta(a,b) draw per taxon. Block missingness is one
  contiguous block per taxon. Lumping maps each class globally to an earlier class;
  splitting creates new per-cell labels. Ascertainment records removed column
  indices; original full truth remains available while observed columns are mapped.
* The sampled-ancestor internal encoding restores the observation's original label
  after removing its zero-length tip; lineage IDs in history remain unique.

## G3 remaining domains

* Stochastic Dollo conditions each non-root-present trait on exactly one origin:
  its gain is uniform over total lineage-time (equivalent to a homogeneous
  Poisson origin process conditioned on one event). Descendants lose at mu times
  their character/clock multiplier. Trait count is fixed; we do not silently
  condition on retention in the observed tips. Invariants have no gain/loss.
* **Acceptance correction:** ordinary binary Fitch CI need not equal one under
  Dollo with repeated losses. A single origin followed by two independent losses
  is valid Dollo evolution but can require >1 unordered-parsimony step. Tests
  therefore check single-origin histories, valid losses, and CI=1 for the
  no-loss special case. Requiring CI=1 with nonzero repeated loss would change the
  requested stochastic Dollo model into an incompatible infinite-sites model.
* Binary Mk uses 0->1 gain_rate and 1->0 loss_rate, times the shared rate/clock.
  Uniformization uses max(gain,loss) and silent proposals, preserving exact
  transition probabilities even when transfers change the state mid-branch.
* Unordered symmetric k-state Mk has off-diagonal rates rate/(k-1); ordered Mk
  keeps only adjacent edges with that same rate, giving lower exit rates at the
  endpoints and uniform equilibrium. An explicit asymmetric matrix supplies
  nonnegative off-diagonal rates and zero diagonal; row sums determine exits.
  Root state is zero; equilibrium is an outcome, not a hidden root assumption.
* Mk correlated blocks use shared uniformization proposals at rho times the
  smallest proposal rate in the block; independent residual clocks keep each
  character's marginal generator exactly unchanged. Each shared proposal draws
  individual transitions, including silent proposals. This exposes dependence
  without inflating marginal rates.
* Textual replacement uses the lexical innovation mechanism plus recurrent easy
  readings. Omission is a separate competing hazard, omission_rate times the
  per-character rate. The empty state is genuine absence; later copying error
  can restore a reading. Binary coding lump/split is rejected: its semantics
  concern class labels, not the biological presence indicator.
* Presets remain provisional until fitted to real target summaries; no empirical
  support is implied by these starting parameter values.

* Evolution is run on the full lineage history before incomplete extant sampling
  and extinction pruning. Otherwise removed tips would incorrectly disappear as
  potential donors and change the process. A separate sampling stream determines
  retained extant observations. `history_tree.nwk` preserves latent lineages when
  different from the observed truth tree; transfer networks refer to that full
  history, while `truth.nwk` matches the retained observation taxa. Stem samples
  are included by introducing an explicit time-zero origin when necessary.

## G7: calibration evidence and release gate

* Sixteen lexical CLDF revisions were retrieved. Only derived summaries, input
  SHA-256s, revision IDs, citations and licenses are tracked. Unjudged cells are
  missing, and partial/alternate cognates form state sets. A 36-point grid with
  three replicates per point, equal family weights and four fractional features
  selects rate=.3, alpha=.7, missing=.15. Residual mismatches are published;
  this preliminary fit is not evidence of universal language-family realism.
* The private stemmatology source returned 404. The 640-MB gene-content and
  6.8-GB morphological inputs are assigned to owner-side extraction per the brief.
  Scripts operate on explicit local source formats and fail on unrecognized
  encodings; citations/licenses must come from source metadata. No synthetic
  numbers stand in for missing empirical targets. Ancestor/contamination rates
  cannot be identified from a matrix alone; reference histories are needed.
* Missing calibration blocks rc1/final, not implementation and documentation of
  the available models. Version remains 2.0.0a3. No PyPI publishing is authorized.
* Single-state CI uses unordered Sankoff; polymorphic data use binary membership
  recoding with an explicit encoding flag. The ambiguous-state shortcut was
  rejected because its minimum-step numerator can exceed its tree-step count.
* Copying donors must predate the copy; newly born simultaneous copies and
  sampled-ancestor continuation points are not extra copying exemplars. This
  also prevents artificial same-time cycles in displayed networks.

* Binary pairwise distances use bitsets and population counts, exactly matching
  observed-cell Hamming distances. This avoids billions of Python cell comparisons
  for 332-species orthogroup matrices. Extraction reuses immutable presence/absence
  state objects and hashes large input files in streaming chunks.
* Final local verification includes stdlib-only tests under `python -S`, actual
  NEXUS/PHYLIP round trips, complete simulation golden outputs under distinct
  PYTHONHASHSEED values, README doctests, all four CLI examples, Ruff and wheel/sdist
  builds. The hosted Python 3.10–3.13 CI matrix remains the cross-version check.

## Post-milestone verification correction

Mean-height rescaling in orchestration must reference retained observations, not
latent extinct tips. Sampling is therefore drawn (with its independent stream)
before evolution without pruning the history; the entire history is rescaled by
the retained tips' height/mean, then characters evolve on that time scale. A
regression test covers extinction, sampled ancestors and incomplete sampling.
Ancestor sampling intensity refers to the original tree-generation time units;
character rates refer to the final rescaled units. Published milestone tags
remain fixed; this correction is a subsequent branch commit.

## Review and empirical follow-up

The count stopping rule, capped last polytomy, sampled-ancestor representation,
and single-origin Dollo interpretation remain as documented: each exposes its
conditioning and rates. Python 3.10 retains the documented TOML subset because
tomllib is unavailable; validation now rejects malformed subset syntax. NEXUS
reading now handles interleaved matrices, comma-separated state sets and symbol
ranges encountered in the public morphology data.

The gene-content and morphology presets adopt the v2 candidate parameters after
their held-out summary errors fell substantially. The lexical candidate worsened
on its holdout, so the previous default stays. The 10% promotion threshold was
chosen after viewing results, and results remain exploratory. The yeast holdout
shares taxa with training; the lexical prior was already fitted to all families.
The stemma reference repository remains unavailable; parameters for surviving
ancestors and contamination remain uncalibrated.

## Benchmark and release review

The 70-dataset comparison uses the consumer's own validation/scoring code and
its NJ and UPGMA implementations. It is a smoke-scale comparison across seven
scenarios, not a statistically powered contest or a calibration objective.
Stemma parameters remain provisional and the morphology source has two damaged
matrices. The rc1 gate therefore remains unmet; no release-candidate tag is made.

## Review follow-up: Yule count stopping (partial)

The earlier first-hitting-time assumption produces zero-length terminal
branches. For constant-rate Yule trees without a deadline or polytomies, we now extend the first n-lineage
interval by an independent exponential wait of rate n*birth. This gives mean
origin height H_n/birth and positive terminal edges. This correction does not
implement the general sampling approach for birth-death or relaxed copying;
those models still require a full GSA implementation before release.
