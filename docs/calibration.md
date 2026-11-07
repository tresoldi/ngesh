# Realism calibration: evidence and remaining work

Status: lexical, gene content and morphology have preliminary empirical fits.
Stemma calibration remains blocked by the inaccessible private reference repository.
Two morphology matrices in the published archive contain only NUL bytes. These
limitations prevent completion of G7 and advancing to 2.0.0rc1/2.0.0 final.

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

## Initial alpha3 comparison (historical)

`calibration/simulated.json` contains five independent simulations of each preset,
including distance/coverage/class quantiles and quartet scores. These comparison
replicates are separate from the fitting replicates. Quantiles summarize per-cell,
per-taxon or per-pair distributions within each dataset, not uncertainty intervals
for fitted parameters. At the time of this initial comparison, the biological targets had not yet been extracted; see the subsequent v2 results below.

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

## Subsequent held-out calibration and GRAPE pilot

The v2 fitting procedure uses dataset-level holdouts for lexical and morphology,
and a deterministic 20% character holdout on the same 332 yeast taxa. The latter
is not an independent biological-family validation. Equal family weights and a
150-character simulation budget limit domination by large datasets. Training
coverage mean and within-family variance determine Beta coverage parameters.
The previous lexical fit used all 16 families, including the later v2 holdouts,
so its comparison is not an independent test of that earlier fit.

| Preset | Previous held-out error | Candidate held-out error | Decision |
|---|---:|---:|---|
| Lexical | 0.02356 | 0.03608 | Retain prior default |
| Gene content | 0.21281 | 0.01655 | Adopt candidate |
| Morphology | 0.01951 | 0.01199 | Adopt candidate |

These are mean squared differences of documented descriptive features, not
accuracy estimates for a phylogenetic method. The promotion rule required at
least 10% reduction; it was recorded after seeing the results. No alternative
candidate was selected using the holdout. One fitting replicate and three
validation replicates give a coarse, preliminary result. The prior lexical
default was fitted to all families; its baseline is not independent. Borrowing,
contamination and absolute historical rates cannot be identified from these
statistics alone. Complete trials, splits, seeds, simulated distributions and
source hashes are in `calibration/fits/`.

The yeast target contains 193,052 reported orthogroups with at least one of
332 ingroup species; 11 published outgroups are excluded using the supplied
species roster. Unclustered genes are not represented. Simulations are projected
to groups present in at least one taxon for calibration; exported simulations
retain unobserved traits. Morphological targets cover 30 of 32 datasets: the
Caviidae matrix was recovered from the archive's TNT-derived NEXUS file; the
lemur and osteogloss matrices and their alternate files contain only NUL bytes.
For Mk calibration, inapplicable gaps are treated as unobserved.

The public archives were accessed with validated HTTP byte ranges. Only selected
members were downloaded; zip CRCs and member SHA-256 hashes were checked. The
published whole-archive checksums were recorded but not verified. Reproduce with:

```sh
python scripts/source_archives.py yeast --work work/sources
python scripts/source_archives.py morphology --work work/sources
python scripts/extract_public.py yeast --root work/sources/yeast --out calibration/targets/genecontent
python scripts/extract_public.py morphology --root work/sources/morphology --out calibration/targets/morphology
python scripts/fit_presets.py --preset genecontent --targets calibration/targets/genecontent --out calibration/fits/genecontent-v2.json
python scripts/fit_presets.py --preset morphology --targets calibration/targets/morphology --out calibration/fits/morphology-v2.json
python scripts/fit_presets.py --preset lexical --targets calibration/targets/lexical --out calibration/fits/lexical-v2.json
```

The GRAPE native CLI at revision `74fec5d78aeb4317253025d8cf05903bd1c6105f`
passed its two native tests and all 17 integration cases, including the expected
rejection of an entirely empty TSV. The pilot verifies roster preservation for
represented taxa, deterministic fixed-seed output, transfer donor lifetimes and
network parsing; it does not measure inference accuracy. Reports are in
`docs/pilot/`, and `scripts/grape_pilot.py` reproduces the checks with a separately
built GRAPE executable. A taxon with no TSV rows is not visible to GRAPE; use the
bundle's observed matrix and taxon roster to detect such omissions.
