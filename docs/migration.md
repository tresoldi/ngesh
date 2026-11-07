# Migration from 1.x

2.x is a breaking model replacement. It does not reproduce 1.x output digests.
There is no statistically faithful mapping from the old gamma-CDF threshold to
a continuous-time replacement rate, so `k`, `th`, `mut_exp`, `k_hgt`, `th_hgt`
are removed. Choose rates in events per character per time unit, and use gamma
`alpha` to express rate heterogeneity without column-position effects.

| 1.x | 2.x replacement |
|---|---|
| `gen_tree(...)` | `generate_tree(TreeConfig(...), seed=...)` from `ngesh.models` |
| ETE3 `Tree` | `ngesh.Tree`, `Node`, `Tree.from_newick`, `tree.newick()` |
| `add_characters` | `evolve(tree, CharacterConfig(...), seed=...)` |
| `simulate_bad_sampling` | `ObservationConfig(taxon_sampling=...)` |
| `tree2wordlist`, `tree2nexus` | `simulate(..., out=...)`, bundle writers |
| `sorted_newick` | Deterministic `tree.newick()` for a given ordered tree |
| `set_seeds` | Private streams; no global seeding API |
| INI files | TOML tables: tree, characters, horizontal, observation, output |
| GUI / random species names | Removed; assign arbitrary node names explicitly |

No old function name is retained with altered statistical meaning. `method`
accepts `standard` and `fast`, both now using the same exact-count event kernel.

Python 3.10 has no `tomllib`. Its dependency-free configuration reader supports
bare-key tables, scalar strings/numbers/booleans, and single-line arrays. It rejects
unsupported constructs such as inline tables, dates, quoted/dotted assignment
keys and multiline strings. Python 3.11+ uses `tomllib`. The shipped TOML examples
work on both. Use Python >=3.11 for the full TOML grammar.

CLI seed `--seed 42` is the string "42", distinct from Python API seed `42`.
Use `simulate(seed="42")` to match it. Byte seeds are recorded as hex in manifests.
Set `stop="time"` plus `max_time` in TOML to disable the default exact-50 stop.

Historical screenshots and the published JOSS manuscript remain historical
artifacts, not examples of 2.x behavior. See the design log for the defect audit
and the calibration document for evidence that is still needed before 2.0 final.
