# ngesh 2

ngesh simulates phylogenetic histories and known-truth character data for
benchmarking. It supports lexical cognate classes, manuscript readings, gene
presence/absence and morphological characters. Runtime and core tests use only
the Python standard library. Python 3.10–3.13 are supported.

**Current development release: 2.0.0a3.** This is a breaking rewrite of 1.2.1.
The lexical, gene-content and morphology defaults have preliminary empirical fits.
The stemma preset remains provisional pending access to the reference traditions. It is not yet a release candidate.

```sh
python -m pip install -e .
ngesh tree --n-taxa 50 --seed demo
ngesh simulate --preset lexical --config examples/lexical.toml --seed demo --out work/pilot
ngesh grid examples/grid.toml --out work/grid
```

Each replicate is a directory containing a time-scaled truth tree, lossless true
and observed matrices, sparse long-format `data.tsv`, binary NEXUS and PHYLIP,
multistate NEXUS, event history, parameters, checksums and descriptive statistics.
Transfers also produce an event table and displayed eNewick network. Sampled
ancestors have both labelled-internal and zero-length-leaf truth encodings.

```python
>>> from ngesh import simulate
>>> result = simulate(seed="example", config={"tree": {"n_taxa": 6}})
>>> len(result.truth)
6
>>> result.groups[0]
'C0001'

```

The returned `Evolution.truth` contains latent states before observation noise.
Bundles contain `observed.json` separately. CLI seeds are strings; API integer,
string and byte seeds are deliberately distinct. Reproducibility uses versioned
SHA-256 derivation and private RNG streams; it never reseeds global generators.

- [Models, equations and parameter meanings](docs/models.md)
- [Bundle format and missing/absent semantics](docs/formats.md)
- [Migration from 1.x](docs/migration.md)
- [Examples for all four domains](examples/README.md)
- [Calibration evidence and remaining inputs](docs/calibration.md)
- [Running design log](docs/DESIGN.md)

Run the core suite with `PYTHONPATH=src python -m unittest discover -s tests -v`.
Development checks: `pip install -e '.[dev]'`, `pytest`, `ruff check src tests scripts`,
`python -m build`. No PyPI publication is part of this work.

Please cite Tresoldi (2021), *Ngesh: a Python library for synthetic phylogenetic
data*, JOSS 6(66), 3173, doi:[10.21105/joss.03173](https://doi.org/10.21105/joss.03173).
See [CITATION.cff](CITATION.cff). Code is MIT licensed; calibration target metadata
records the separate licenses and citations of the source datasets.
