---
title: 'Ngesh: a Python library for phylogenetic simulation'
tags:
  - Python
  - phylogenetics
  - random phylogenetic tree
  - phylogenetic tree simulation
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

This work presents [`ngesh`](https://pypi.org/project/ngesh/), a Python library
for simulating phylogenetic trees and data, primarily designed for
reserach in historical linguistics and stemmatics. 
It generates reproducible stochastic simulations,
according to various criteria and in a range of output formats,
for the development, debugging, and benchmarking of
phylogenetic methods and tools.

# Background

Following the popularization of techniques for Bayesian inference,
computational phylogenetics is being adopted in fields beyond biology,
including historical linguistics [@Bouckaert:2012] and stemmatics
[@Robinson:2016]. Despite some remarkable research outcomes,
such adoptions still prompt some concerns
on the adaptability of the underlying evolutionary metaphors into other disciplines,
calling for demonstrations of the efficiency of such methods.
Stochastic simulations, long advocated for
natural sciences in general [@Bailey:1990] and phylogenetics in specific
[@Foote:1999; @Harmon:2019], are a practical answer to this issue,
as they offer means to assess performance and applicability
by allowing to design and study extensive amounts of simulated evolutions, without
time limits for data collection and with complete knowledge of the
expected results.
As side-effects, simulations provide
basic fuzzy testing of software pipelines and can help in analyzing
which models and which sets of criteria better match
observed phenomena.

The `ngesh` library provides resources for such reproducible simulations,
generating trees and correlated data following user-defined parameters,
such as birth and death ratios, and constrains, such as
maximum branch length or number of extant nodes, optionally
manipulating the results in diverse ways, for example by pruning extinct leaves
or simulating uneven sampling.
Character evolution analogous to the topology can likewise be simulated,
with different rates for *ex novo* mutation and horizontal gene transfer.
Taxa can be labeled either in progression or randomly, allowing both human-readable names
(like "Sume", "Fekobir", and "Tukok") and binominal biological nomenclature
(like "Sburas wioris", "Zurbata ceglaces", and "Spellis spusso").
The resulting trees are standard ETE objects [@ETE:2016] and can be exported in a
variety of formats, including Newick trees, ASCII-art representation, and tabular
textual listings.

# Installation, Usage, & Examples

The library can be installed with the standard `pip` tool for Python package
management.
Trees can be generated from the command-line, defaulting to small
phylogenies in Newick format:

```bash
$ ngesh
(Ukis:1.11985,(Koge:0.880823,(Rozkob:0.789548,(Meu:0.706601,
(((Felbuh:0.189693,Kefa:0.189693)1:0.117347,((Epib:0.153782,
Vugog:0.153782)1:0.0884745,Puluk:0.242256)1:0.0647836)1:0.0469885,
Efam:0.354028)1:0.352573)1:0.0829465)1:0.0912757)1:0.23903);
```

The tool allows both configuration files and command-line flags overriding
settings. Here we define a model to generate Nexus data for a reproducible Yule
tree with a birth rate of 0.75, at least 5 leaves,
"human" labels, and 20 presence/absence features:

```bash
$ cat my_tree.conf
[Config]
labels=human
birth=0.75
death=0.0
output=nexus
min_leaves=5
num_chars=20
$ ngesh -c my_tree.conf --seed 12345
#NEXUS

begin data;
  dimensions ntax=6 nchar=33;
  format datatype=standard missing=? gap=-;
  matrix
Buza      111110110111011011010101000100110
Lenlar    111111010110111101100010010011001
Mukom     111110111011011011101001000100110
Pagil     111110110111011011100100100100110
Suglu     111110110111011011100011001001010
Wite      111110110111011011100101000100110
  ;
end;
```

Despite the convenience of its command-line tool, the package is designed
for usage as a library. The two main simulation functions are `gen_tree()`,
which returns a random tree, and `add_characters()`, which inserts character
evolution data in a tree object. As they are separate methods, it is
possible to generate a random tree without character information or
to simulate character evolution within existing trees, including
non-simulated ones.

```python
>>> import ngesh
>>> tree = ngesh.gen_tree(1.0, 0.5, max_time=0.5, labels="bio",
                          seed="abc")
>>> print(tree)

      /-Nedoros seveddi
   /-|
--|   \-Radabmo toras
  |
   \-Wirror ubislis
>>> print(tree.write())
((Nedoros seveddi:0.158084,Radabmo toras:0.158084)1:0.325586,
Wirror ubislis:0.48367);
>>> tree = ngesh.add_characters(tree, 15, 2.0, 0.5)
>>> print(ngesh.tree2nexus(tree))
#NEXUS

begin data;
  dimensions ntax=3 nchar=22;
  format datatype=standard missing=? gap=-;
  matrix
Nedoros_seveddi    1111011011101010101101
Radabmo_toras      1111011101011010101101
Wirror_ubislis     1111101011010101011110
  ;
end;
```

# Alternatives

A recommended alternative for simulating phylogenetic processes, despite
no specific support for historical linguistics or stemmatics, is the
R package `TreeSim` by @Stadler:2011;
methods of the `geiger` package [@Pennell:2014] and the `rtree()` function of
the `ape` package [@Paradis:2018]
might be satisfactory depending on the requirements.
In Python, @Noutahi:2017 gives a demonstration script comparable to `ngesh`,
and the `populate()` method of ETE’s `Tree` class [@ETE:2016]
can generate simpler simulations of character evolution.
In all languages, manual randomization of taxa
placement in existing cladograms is a well-known alternative to stochastic
simulation.

# Code and Documentation Availability

The `ngesh` source code is available on GitHub at
[https://github.com/tresoldi/ngesh](https://github.com/tresoldi/ngesh).

The full user documentation is available at
[https://ngesh.readthedocs.io/](https://ngesh.readthedocs.io/).

# Acknowledgements

The author has received funding from the European Research Council (ERC)
under the European Union’s Horizon 2020 research and innovation
programme (grant agreement
No. [ERC Grant #715618](https://cordis.europa.eu/project/rcn/206320/factsheet/en),
[Computer-Assisted Language Comparison](https://digling.org/calc/)).

# References
