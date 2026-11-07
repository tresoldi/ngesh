# Runnable 2.x examples

Run from the repository after `python -m pip install -e .`:

```sh
ngesh simulate --preset lexical --config examples/lexical.toml --seed demo --out work/lexical
ngesh simulate --preset stemma --config examples/stemma.toml --seed demo --out work/stemma
ngesh simulate --preset genecontent --config examples/genecontent.toml --seed demo --out work/genecontent
ngesh simulate --preset morphology --config examples/morphology.toml --seed demo --out work/morphology
ngesh grid examples/grid.toml --out work/grid
```

Each `simulate` destination must be new. `grid` verifies completed bundles and
resumes. The test suite executes all four preset examples in temporary directories.
The 1.x INI and XML artifacts were removed; they describe the replaced model.
