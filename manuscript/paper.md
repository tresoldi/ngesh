---
title: 'ngesh: a Python library for phylogenetic tree simulation'
tags:
  - Python
  - phylogenetics
  - random phylogenetic tree
  - simulation
authors:
  - name: Tiago Tresoldi
    orcid: 0000-0002-2863-1467
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: Department of Linguistic and Cultural Evolution, Max Planck Institute for the Science of Human History
   index: 1
date: 02 April 2020
bibliography: paper.bib
---

# Summary

This work describes [`ngesh`](https://pypi.org/project/ngesh/), a Python library
and related command-line tools for simulating phylogenetic
trees and related data (characters, states, branch length, etc.). It is intended
for benchmarking phylogenetic methods, especially in historical linguistics, and
for providing dummy trees for their development and debugging. The generation of
random phylogenetic trees also goes by the name of "simulation methods for
phylogenetic trees" or just "phylogenetic tree simulation".

In detail, with ngesh:

  - trees can be generated according to user-specified birth and death ratios
(and the death ratio can be set to zero, resulting in a birth-only tree)
  - speciation events default to two descendants, but the number of descendants
can be randomly drawn from a user-defined Poisson process (allowing to model
hard politomies)
  - trees will have random topologies and, if necessary, random branch-lengths
  - trees can be limited in terms of number of extant leaves, evolution time (as
related to the birth and death parameters), or both
  - non-extant leaves can be pruned from birth-death trees
  - character evolution can be simulated in relation to branch lengths, with
user-specified ratios for mutation and horizontal gene transfer, with different
rates of change for each character
  - trees can be generated from user-provided seeds, so that the random
generation can be maintained across executions (and, in most cases, the
execution should be reproducible also on different machines and different
vestions of Python)
  - nodes can optionally receive unique labels, either sequential ones (like
"L01", "L02", and "L03"), random human-readable names (like "Sume", "Fekobir",
and "Tukok"), or random biological names approximating the binomial nomenclature
standard (like "Sburas wioris", "Zurbata ceglaces", and "Spellis spusso")
  - trees can be returned as ETE tree objects or exported in a variety of
formats, such as Newick trees, ASCII representation, tabular textual listings,
etc.

# Background

Just to cite @Bailey:1990, @Foote:1999, @Harmon:2019, @Stadler:2011,
@Noutahi:2017 

# Installation, Usage, & Examples

The library can be installed with the standard `pip` tool for package
management:

```bash
$ pip install ngesh
```

The [documentation](https://ngesh.readthedocs.io/en/latest/) offers detailed
instructions on how to use the library. For most purposes,

# Alternatives

There are many tools for simulating phylogenetic processes in order to
obtain random phylogenetic trees. The most complete is probably the R
package `TreeSim` by Tanja Stadler, which includes many flexible tree
simulation functions. In R, one can also use the `rtree()` function from
package `ape` and the `birthdeath.tree()` one from package `geiger`,
as well as manually randomizing taxon placement in cladograms.

In Python, some code similar to `ngesh` and which served as initial
inspiration is provided by @Noutahi:2017.

For simpler simulations, the `.populate()` method of the `Tree` class in ETE
might be enough as well. Documentation on the method is available here.

A number of on-line tools are also available at the time of writing:

  - T-Rex (Tree and reticulogram REConstruction at the Université du Québec à Montréal (UQAM)
  - Anvi'o Server can be used on-line as a wrapper to T-Rex above
  - phyloT, which by randomly sampling taxonomic names, identifiers or protein accessions can be used for the same purpose


# Code and Documentation Availability

The `ngesh` source code is available on GitHub at
[https://github.com/tresoldi/ngesh](https://github.com/tresoldi/ngesh).

The documentation is available at
[https://ngesh.readthedocs.io/](https://ngesh.readthedocs.io/).

# Acknowledgements

The author has received funding from the European Research Council (ERC)
under the European Union’s Horizon 2020 research and innovation
programme (grant agreement
No. [ERC Grant #715618](https://cordis.europa.eu/project/rcn/206320/factsheet/en),
[Computer-Assisted Language Comparison](https://digling.org/calc/)).

# References
