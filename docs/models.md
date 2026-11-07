# Models and parameters

Times increase forward from an origin at zero. Lengths are elapsed time, not
expected substitutions. Rates have units of events per character per time unit
unless stated otherwise. See `TreeConfig`, `CharacterConfig`, `HorizontalConfig`
and `ObservationConfig` for the complete typed parameter list and defaults.

## Tree history

Each extant lineage gives birth at λ and dies at μ. The next waiting time is
Exponential(sum of all active birth and death rates). A lineage is selected in
proportion to its total event rate. Binary births add one extant lineage; optional
polytomies have 2+Poisson(polytomy) children. `copy_sigma` gives mean-one
lognormal birth-rate multipliers to each new lineage. With `burst_probability`,
a birth has Poisson(`burst_mean`) extra children.

`n_taxa` stops at the first exact extant count. At the final birth, offspring are
capped to that count. `max_time` additionally rejects histories that fail to hit
the count by the deadline. Use `stop="time"` and `max_time` for time-only
simulation (CLI: `ngesh tree --time-only --max-time 3`). Failed/extinct histories
are retried up to `max_attempts`; populations are bounded by `max_lineages`.
With no death, the expected exact-n origin height is
E[H_n] = sum(k=1..n-1) 1/(k λ). This includes the first stem waiting time.

`full` retains extinct tips. Extant tips are independently sampled with
probability `sampling` (ρ). `ancestor_rate` samples points along lineage-time as
a Poisson process, retaining a continuing lineage and a zero-length observation.
The manuscript preset combines high loss, irregular copying and observed ancestors.
The simulator evolves characters on the full history before pruning observations.

`shape` may be `balanced`, `caterpillar`, `birth_death`, or `newick` with `newick`
text. Controls have unit height; balanced splits differ by at most one tip.
`height` rescales origin height; `mean_height` rescales mean observed tip time.
They are mutually exclusive. Pruning preserves absolute observation times.

## Character processes

Each character has a rate r_c = rate × G_c, with G_c=1 by default. For `alpha`,
G_c ~ Gamma(alpha, scale=1/alpha). `categories>0` uses equal-probability gamma
bins and their conditional means. A separate `invariant` fraction has zero rate;
the mixture's unconditional mean is rate×(1-invariant). A branch multiplier is
L_b = exp(N(-clock_sigma²/2, clock_sigma²)), with E[L_b]=1. Event intensity is
r_c L_b. Event counts on a branch of length t are Poisson(r_c L_b t) for lexical
replacements. No cap on changes per branch is imposed. `max_events` bounds memory.

**Lexical.** Each concept begins at class 0. Replacement creates a unique class,
except with `p_parallel`, when a class from the recurrent pool 1..`pool_size` is
drawn uniformly. A recurrent draw can be silent. New classes have larger IDs.
`synonym_duration>0` retains old classes for independent exponential lifetimes
with that mean. Expiry times are inherited across splits. Without parallelism,
borrowing or synonymy, multistate characters have CI=1 on their truth tree.

**Gene content.** `model="dollo"` gives each trait initial presence probability
`root_frequency`. Non-root-present, noninvariant traits have exactly one origin,
uniform over total lineage-time: this is a homogeneous origin process conditioned
on one gain. Losses have intensity loss_rate×r_c×L_b. Traits may disappear entirely;
there is no hidden retention ascertainment. `model="binary"` instead has repeated
0→1 gains at gain_rate×r_c×L_b and 1→0 losses at loss_rate×r_c×L_b. Its equilibrium
presence is gain/(gain+loss). Dollo with repeated losses can have ordinary CI<1;
single origin, rather than CI=1, is its defining constraint.

**Morphology.** `model="mk"`, k=2..6. Unordered symmetric Q_ij=1/(k-1), i≠j,
multiplied by r_c L_b. Ordered Mk retains only adjacent transitions at that same
rate; endpoint exit rates are lower. Both have uniform stationary frequencies.
An asymmetric `rate_matrix` supplies off-diagonal rates and zero diagonal. All
rows must have k entries. Root state is zero. Uniformization uses the maximum
exit rate, including silent proposals. For symmetric unordered Mk,
P_ii(t)=1/k+(1-1/k)exp[-k r_c L_b t/(k-1)].

Correlated Mk blocks have `block_size` and `block_correlation` (ρ_block). Shared
uniformization proposals occur at ρ_block times the smallest proposal intensity
in the block; independent residual proposals preserve each marginal Q exactly.
Each character draws its transition independently at a shared proposal, permitting
simultaneous changes without inflating its marginal rate.

**Textual.** Unique new readings and recurrent easy errors use the lexical model.
Omission occurs at omission_rate×r_c×L_b and leaves an empty cell. A later copying
error can restore a reading. Contamination acts at copying time (below).

## Horizontal processes

`borrowing_rate` is a proposal intensity per character per lineage-time. A donor
must be a different lineage alive at that time; proposals without a donor are
skipped. Silent events remain recorded. Donor weight is exp(-d/scale), where d
is the patristic distance **at the event time**, or a distance in simulated
geography. Geography has two Brownian coordinates, each with increment variance
diffusion²×elapsed time; children inherit parental positions.

`mode="replace"` establishes the donor's current state set. `mode="add"` unions
presence states and is restricted to binary/Dollo models. Synonym deadlines are
reset by a transfer. No future donor state or ancestral endpoint is consulted.

At each copying split, `contamination_fraction` selects Bernoulli(f) scattered
loci or, with `contamination_blocks`, one block of round(f×L) loci. The second
exemplar exists at copying time and is not the recipient or its new sibling.
When no second exemplar exists, there is no contamination. Multiple loci share
one reticulation edge when transferred together.

## Observation

Observation has an independent stream and leaves truth unchanged. `missing` is
MCAR per cell. `coverage_a,b` draw per-taxon Beta coverage; a cell is retained
with that probability. `block_missing` masks one contiguous fraction of loci.
`selection` independently retains loci at each copying branch using inherited masks;
sampling points and continuation branches preserve their masks. Collections can
only shrink, producing nested, known absences rather than missing
cells. Lumping globally merges a class into an earlier class with its specified
probability; splitting assigns a fresh per-cell state. These class-label errors
are undefined for binary presence/absence and are rejected there.

`remove_invariant` removes constant observed state sets; `remove_singleton`
removes groups containing any state observed in only one taxon. All decisions
are logged. `taxon_sampling` is Bernoulli subsampling with consistent truth-tree
pruning. A run that loses every taxon fails explicitly.

## Summary statistics

Classes and singleton fractions count observed class memberships. Coverage
counts nonmissing cells, including known absences. Pairwise distance averages
per-group Jaccard distance over jointly observed cells; empty/empty is zero.
For binary data this is observed-cell Hamming distance. This is one documented
benchmark distance, not a recommendation for scoring nested manuscript collections.

For a quartet, sort its three pair sums S1≤S2≤S3. δ=(S3-S2)/(S3-S1), or zero for
an unresolved equal-sum quartet. Q-residual averages (S3-S2)² and divides by the
squared mean observed pair distance. Missing-distance quartets are omitted.
By default at most 1000 evenly spaced lexicographic quartets are used; counts and
sampling rule are recorded. Holland et al. (2002), *Delta plots: a tool for
analyzing phylogenetic distance data*, Molecular Biology and Evolution 19:2051–2059.

Unordered CI uses Sankoff dynamic programming (including multifurcations), with
missing cells unrestricted. Constant data have undefined CI (JSON null). For
polymorphic matrices, CI uses binary membership recoding and records that choice;
treating all states as ambiguous would give an invalid minimum-step numerator.
Ordered/asymmetric evolution is summarized with this same explicitly unordered CI.
