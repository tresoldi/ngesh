# Realism calibration: evidence and remaining work

Status: lexical has a **preliminary coarse grid fit**. Stemma, gene content and
morphology have implemented simulators and extraction tools, but **no empirical
fit is claimed**. Their inputs are unavailable in this session. This prevents
claiming completion of G7 or advancing to 2.0.0rc1/2.0.0 final. No release was
published to PyPI.

## Reproducible inputs

Only small derived statistics, input hashes, exact Git revisions, bibliographic
citations and source licenses are committed in `calibration/targets/lexical/`.
No third-party forms, cognate tables, trees or sequence data are committed.
Extraction restricts attention to taxa and concepts with at least one cognate
judgment; cells without a judgment are missing. Partial cognates and alternate
forms become class memberships, not fabricated unique classes. This choice is
important for Burmish and other partially cognate datasets.

All 16 requested lexical sources were retrieved. Each JSON records the original
citation and license, which remain applicable to those sources. Most are CC BY
4.0; **Kitchen Semitic is CC BY-NC 4.0**. Derived summaries here do not relicense
source data under the software's MIT license. Source-specific citations and
licenses are machine-readable; see their repository metadata before reuse.

```sh
python scripts/calibrate.py fetch-lexical --work work/calibration --out calibration/targets/lexical
python scripts/calibrate.py fit --targets calibration/targets/lexical --preset lexical --out calibration/lexical-fit.json
```

Fetching preserves existing checkouts; delete a checkout explicitly if an updated
source revision is wanted. The fit records all 36 grid points, with three
replicates each. Rates are {0.03,0.1,0.3,0.6}, alpha {0.3,0.7,2}, MCAR missing
{0,0.15,0.4}. Simulations use 50 taxa and 200 concepts. We minimize the sum of
squared differences in four fractional features: mean pairwise distance,
singleton fraction, missing fraction, and mean classes per group divided by taxon
count. Families have equal weight, so large datasets do not dominate. Seeds and
the algorithm are fixed in the script. No continuous optimization is hidden.

The selected defaults are **rate=0.3, alpha=0.7, missing=0.15**. The fitted and
real features below show a meaningful remaining mismatch; these are benchmark
starting points, not a validated model of all language families. Tree-likeness
and rate/borrowing identifiability are not established by this fit. A larger
validation grid, family-matched taxon counts and held-out targets are recommended
before freezing final-release defaults.

| Feature | Real (equal family mean) | Simulated (three fitting replicates) |
|---|---:|---:|
| Pairwise distance | 0.52896 | 0.51130 |
| Singleton fraction | 0.47196 | 0.55713 |
| Missing fraction | 0.13016 | 0.15210 |
| Classes per group / taxon count | 0.28548 | 0.19830 |

## Real and simulated distributions

`calibration/simulated.json` contains five independent simulations of each preset,
including distance/coverage/class quantiles and quartet scores. These comparison
replicates are separate from the fitting replicates. Quantiles summarize per-cell,
per-taxon or per-pair distributions within each dataset, not uncertainty intervals
for fitted parameters. Unknown real distributions are marked pending, not
replaced by synthetic targets.

| Preset | Real pairwise-distance quantiles (min, p10, p25, p50, p75, p90, max) | Simulated quantiles (five-replicate mean) |
|---|---|---|
| genecontent | pending source extraction | 0.000, 0.112, 0.181, 0.218, 0.239, 0.268, 0.302 |
| lexical | 0.110, 0.353, 0.456, 0.553, 0.626, 0.674, 0.732 | 0.000, 0.381, 0.489, 0.552, 0.582, 0.605, 0.666 |
| morphology | pending source extraction | 0.000, 0.296, 0.396, 0.452, 0.488, 0.516, 0.600 |
| stemma | pending source extraction | 0.000, 0.733, 0.890, 0.958, 0.978, 0.986, 1.000 |

Lexical sources, with equal-family aggregation above:

| Dataset | Taxa | Concepts | Mean distance | Singletons | Missing |
|---|---:|---:|---:|---:|---:|
| [chaconarawakan](https://github.com/lexibank/chaconarawakan/tree/d7874df4a7d3de37c436017bcca9fda443641d4d) | 8 | 102 | 0.497 | 0.524 | 0.252 |
| [constenlachibchan](https://github.com/lexibank/constenlachibchan/tree/def48749fb269a0dccb7ae164d6ccdf85f38a559) | 25 | 110 | 0.836 | 0.695 | 0.044 |
| [dravlex](https://github.com/lexibank/dravlex/tree/37578075e5ccb09c43022f7f1282125e748de84d) | 20 | 100 | 0.631 | 0.640 | 0.073 |
| [dunnaslian](https://github.com/lexibank/dunnaslian/tree/2cbc18138908440e01a6e046cc5101e28855b344) | 32 | 146 | 0.634 | 0.398 | 0.160 |
| [iecor](https://github.com/lexibank/iecor/tree/dd1ebc997bc38c824fb18a1ef234bbf2c229420b) | 160 | 170 | 0.774 | 0.470 | 0.065 |
| [kitchensemitic](https://github.com/lexibank/kitchensemitic/tree/de5d0fd71ab09cb6b4d3d77c9b160bd520184883) | 25 | 97 | 0.633 | 0.495 | 0.118 |
| [leekoreanic](https://github.com/lexibank/leekoreanic/tree/b2ca6413274d592ee15deb59d35ef71e4957d608) | 15 | 246 | 0.170 | 0.219 | 0.372 |
| [liusinitic](https://github.com/lexibank/liusinitic/tree/54f6742d9fa60315ae41b91d0d1e02f04036efb5) | 19 | 203 | 0.315 | 0.563 | 0.000 |
| [mannburmish](https://github.com/lexibank/mannburmish/tree/513edae24c395440ceb23ead61b933b1f52c6d0c) | 7 | 391 | 0.170 | 0.002 | 0.242 |
| [savelyevturkic](https://github.com/lexibank/savelyevturkic/tree/4a540590580f0981d1c5679668cd08578f856de2) | 32 | 254 | 0.315 | 0.271 | 0.066 |
| [sidwellbahnaric](https://github.com/lexibank/sidwellbahnaric/tree/be17c6e8a96d4c1e457bfa90917dd75ccb6e51af) | 24 | 200 | 0.455 | 0.498 | 0.104 |
| [syrjaenenuralic](https://github.com/lexibank/syrjaenenuralic/tree/f4b9010d5d5d9fc5208d32a75bca0e488697357d) | 7 | 173 | 0.748 | 0.710 | 0.016 |
| [tuled](https://github.com/tupian-language-resources/tuled/tree/91c20e199ad3963d41c8987ac59431c4706f7645) | 89 | 447 | 0.594 | 0.606 | 0.423 |
| [utoaztecan](https://github.com/lexibank/utoaztecan/tree/13aaffc2f50124f3e77e6086d6a1ff5a9e5c1653) | 41 | 121 | 0.756 | 0.504 | 0.066 |
| [walworthpolynesian](https://github.com/lexibank/walworthpolynesian/tree/515aa7fd5210949e79b57d5dd5dfb8a8ae9644d7) | 31 | 210 | 0.451 | 0.544 | 0.050 |
| [wichmannmixezoquean](https://github.com/lexibank/wichmannmixezoquean/tree/6793ec2fdfd2e47bc0bbb63c26bc60d93995f35f) | 10 | 110 | 0.483 | 0.413 | 0.030 |

## Inputs needed from the owner

**Stemma:** `evotext/stemmatology` is private and returned 404 through the connected
GitHub account. Supply an archive or an accessible checkout containing
`binary/*.bin.nex`, `multistate/*.nex`, and reference stemmata under `testing/`.
The traditions are Heinrichi, Parzival, Notre Besoin and Legend. Example:

```sh
python scripts/calibrate.py extract binary-nexus work/stemmatology/binary/TRADITION.bin.nex --citation 'evotext/stemmatology, tradition and revision' --license 'COPY FROM SOURCE' --out calibration/targets/stemma/TRADITION.json
```

Use `nexus` for multistate files. Do not combine binary recodings and their
multistate source as independent calibration datasets. Reference stemmata and
known artificial-tradition provenance are also required to estimate observed
ancestor and contamination rates; matrices alone cannot identify these rates.
Those parameter fits remain pending rather than reporting arbitrary estimates.

**Gene content:** Shen et al. (2018), *Tempo and mode of genome evolution in the
budding yeast subphylum*, Cell 175:1533,
[doi:10.1016/j.cell.2018.10.023](https://doi.org/10.1016/j.cell.2018.10.023).
[Figshare article 5854692](https://figshare.com/articles/dataset/5854692) includes
`4_orthomcl_output.zip` (~640 MB) and `1_molecular_phylogenetics.zip`. The owner
can download and extract these outside version control. Find the OrthoMCL groups
file and confirm its species/gene identifier convention. The extraction script
requires an explicit species separator (default `|`) and fails if it is absent:

```sh
python scripts/calibrate.py extract orthogroups work/yeast/groups.txt --separator '|' --citation 'Shen et al. 2018; figshare 5854692; exact file' --license 'COPY FROM FIGSHARE METADATA' --out calibration/targets/genecontent/yeast.json
```

Each line must be `group: species|gene species|gene ...`. If source IDs use a
separate mapping, apply that documented mapping before extraction; do not guess
species from arbitrary prefixes. Check the resulting taxa against the 332-species
tree and document exclusions. Absent orthogroups are known zero, not missing.
Use the tree for validation of low-transfer defaults after extraction.

**Morphology:** Keating, Garwood & Sansom (2023), *Phylogenetic congruence, conflict
and consilience between molecular and morphological data*, BMC Ecology and
Evolution 23:30,
[doi:10.1186/s12862-023-02131-z](https://doi.org/10.1186/s12862-023-02131-z).
[Zenodo 10.5281/zenodo.6579584](https://doi.org/10.5281/zenodo.6579584),
`Keating_et_al_sup.zip` (~6.8 GB), contains 32 paired metazoan datasets. The owner
runs extraction locally. Use only the morphological STANDARD matrices; the reader
rejects DNA/RNA/protein/mixed matrices. If a STANDARD matrix includes multiple
partitions, specify the 1-based inclusive morphology range:

```sh
python scripts/calibrate.py extract nexus work/keating/morphology.nex --characters 1:152 --citation 'Keating et al. 2023; Zenodo 6579584; exact dataset' --license 'COPY FROM ZENODO METADATA' --out calibration/targets/morphology/DATASET.json
```

Omit `--characters` if all columns are morphological. Check each source's ordered
characters, polymorphisms and inapplicable-state encoding before extraction.
`?` is missing; `-` is treated as known empty/inapplicable under ngesh's generic
contract, so recode it to `?` first if the source defines it as unknown. Molecular
partitions must never enter these targets. Record any recoding alongside input
hashes and citations. The script accepts extracted matrices rather than silently
unpacking or loading multi-gigabyte archives into RAM.

Source license values above are deliberately placeholders requiring inspection
of the attached archive metadata. No license or real-data statistics are invented.
Once these targets and reference trees exist, extend/fit the domain-specific
parameter grids, compare held-out distributions and advance through rc1 review.

## Performance evidence

The pilot `2.0.0a1` run wrote all 1000 bundles with 50 taxa, 200 concepts, full
checksums and default summaries in 285.71 s wall / 285.62 s CPU on Python 3.12.14,
Linux x86-64, one Python process. This used the then-provisional lexical rate 0.1,
with no missingness; see `docs/benchmarks/lexical-a1.json`. A final calibrated
preset run is recorded separately after the remaining changes; neither timing
is a promise for all hardware or for exhaustive O(n^4) quartet statistics.

The calibrated `2.0.0a3` preset (rate .3, alpha .7, missing .15) wrote 1000 bundles
in **320.65 s wall / 320.51 s CPU (5.34 minutes)** on the same runtime. See
`docs/benchmarks/lexical-calibrated-a3.json`. This measurement preceded final
metadata additions; model parameters, dimensions, quartet limit and core process
were the same. Each run is a single sequential process; there are no worker pools
in the benchmark simulator.
