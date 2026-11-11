# Changelog

## Unreleased

- Record where contact episodes' copied readings sit in the observed data:
  `episode_cells.tsv` (carriers, holders, donor-side holders per episode and
  character), and `effective` and `carriers` columns in `episodes.tsv`. The
  founders' living descendants miss most episodes when donor lines are lost, as
  in manuscript traditions. Runs without episodes are unchanged.

- Add contact episodes (`episode_rate`, `episode_duration`, `episode_intensity`,
  `episode_symmetric`): sustained transfer between two populations that
  continues through splits, logged in `episodes.tsv` and in an `episode` column of
  `events.tsv`. A zero rate reproduces earlier outputs exactly.

- Generate deep histories without recursion: the GSA snapshot in
  `_birth_death` copies the tree iteratively (`tree.copy_subtree`) instead of
  `copy.deepcopy`, which raised `RecursionError` for some 100-taxon stemma
  histories; `Tree.from_newick` parses with an explicit stack, so nesting depth
  is no longer bounded by the recursion limit. Outputs are unchanged.

## 2.0.0a4 — review fixes

- Use target-count interval sampling and distinguish crown and origin heights.
- Preserve sampled witnesses, apply standard ascertainment after subsampling,
  and write single-symbol Mk matrices.
- Share lineage clocks across sampling points, weight Dollo gains by exposure,
  label silent omissions and retry empty tip samples.
- Refit available calibration targets and rerun the GRAPE comparison.


## 2.0.0a3 — 2026-11-07

- Add single-origin Dollo, binary and ordered/asymmetric Mk, correlated blocks,
  textual errors and omissions, manuscript and other domain presets.
- Evolve full lineage histories before pruning observations.
- Add empirical lexical summaries and calibration extraction/fitting tools.
- Calibration for stemma, gene content and morphology awaits external inputs;
  2.0.0rc1 and 2.0.0 final are not claimed or published.

## 2.0.0a2 — 2026-11-07

- Add chronological borrowing, Brownian geography and copying contamination.
- Add missingness, nested selection, coding error, ascertainment and subsampling.

## 2.0.0a1 — 2026-11-07

- Replace ETE3/NumPy core with standard-library timed trees and stable private RNGs.
- Add exact-count birth–death, sampled ancestors and control trees.
- Add timed lexical evolution, gamma rates, bundles, CLI and resumable grids.
- Remove incompatible 1.x APIs and INI configuration; require Python >=3.10.


Version 1.2:
  - Added exact number of leaves as stop criterion

Version 1.1.1:
  - Small bugfixes, related to JOSS publication

Version 1.0:
  - Release for JOSS publication

Version 0.5:
  - Code improvements preparing for release and submission
  - Included new faster generation method contributed by @NicolaDM

Version 0.4.1:
  - Fixed bug with character generation during simulation of horizontal
    gene transfer
  - Fixed bug in upper boundary number of syllables for label generation

Version 0.4:
  - General improvements for submission
  - Extended function documentation and in-code comments
  - Documentation at ReadTheDocs.io
  - Complete output (whenever possible) for all output formats
  - Complete command-line access to all parameters, such as for hard
    politomies
  - Automatic code review, more tests
  - Removed dependency on the `abzu` library

Version 0.3.1:
  - Code improvements for next release and submission

Version 0.3:
  - General improvements to code quality
  - Full reproducibility from seeds for the pseudo-random generators,
    allowing string, ints, and floats
  - Changes for further integration with `abzu` and `alteruphono` for
    simulating linguistic data
