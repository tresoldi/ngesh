# Dataset bundle contract

Every replicate is a new directory. Deterministic ordering, UTF-8, LF endings and
JSON with sorted keys make checksums reproducible. `params.json` is written last.

| File | Meaning |
|---|---|
| `truth.nwk` | Primary tree for retained observations, branch lengths in time |
| `truth_sa_leaves.nwk` | Sampled ancestors as zero-length observed tips |
| `truth_sa_internal.nwk` | Same observations as labelled internal nodes |
| `history_tree.nwk` | Full latent lineage history, if sampling/pruning changed it |
| `truth.json` | Unfiltered, latent states for the initially retained taxa |
| `observed.json` | Final observed matrix, with filters applied |
| `data.tsv` | Sparse long format for class-membership consumers |
| `data.bin.nex`, `data.bin.phy` | Binary membership matrices, `?` for missing |
| `data.multi.nex` | Integer-token multistate NEXUS; Mk partition and ordering |
| `lineages.json` | Full lineage lifetimes, parent IDs, extinction and sampling flags |
| `columns.json` | Group ordering, observed taxon names, state-column mapping, PHYLIP aliases |
| `history.tsv` | Every evolutionary event, including silent events |
| `events.tsv`, `network.enwk` | Horizontal events and displayed event network, when any occurred |
| `observations.json` | Separate observation and ascertainment log |
| `summary.json` | Descriptive statistics and quartet sampling information |
| `params.json` | Effective parameters, seeds, version, rates and file hashes |

Matrices use `null` for missing, `[]` for known absence, and `[s,...]` for observed
state sets. `truth.json` has original groups; `observed.json` follows the retained
indices in `columns.json`. Taxon subsampling can remove rows only from observed
outputs. The network includes the full donor history, even if a donor or its
descendants are absent from observed data. Use `history_tree.nwk` with event
lineage IDs when truth pruning has collapsed branches.

The TSV header is exactly `Language<TAB>Concept<TAB>Cognateset`. Each later line
has exactly three columns by position: taxon, group, and `<group>.<integer>`.
Integers are nonnegative and refer to classes within a group. Prefixes are ignored
by GRAPE; classes from separate groups are never implicitly equated. Taxa and
states are sorted; groups follow their deterministic declared order. Each line,
including newline, is strictly under 4096 UTF-8 bytes. Labels with tabs or newlines
are rejected. Synonyms have extra rows; missing cells have no rows.

**Binary presence/absence:** each original trait is one group, present is state
`.1`; known absence and missing both have **no TSV row**. Sparse TSV alone cannot
recover that distinction, nor even recover taxa that have no present traits.
Use `observed.json`, NEXUS or PHYLIP, with `columns.json`, when absence matters.
Binary NEXUS/PHYLIP retain all original traits, even entirely absent ones. For
lexical/morphological data, binarization creates one column per observed state;
missing source cells become `?` across all those columns.

Multistate NEXUS uses integer `TOKENS` to allow more than 62 lexical classes;
polymorphisms use braces, genuine absence uses `-`. Consumer support for integer
TOKENS varies; binary formats are available for single-symbol consumers. PHYLIP
uses stable `T000001` aliases because whitespace in names is not portable.

`params.json` contains every effective model/observation/tree parameter and all
component seeds and their actual SHA-derived Random integer seeds, plus a SHA-256 for **every other file**, including summaries.
It excludes itself to avoid recursive hashing. `verify_bundle` checks this
manifest. The versioned seed framing distinguishes integer, string and byte
seeds. Grid replicate seeds derive from master seed, factorial combination index
and replicate index; combination axes are sorted. Resuming a different config,
seed, preset or version in the same destination fails instead of silently mixing
runs. Interrupted grid staging directories are disposable and are regenerated.
