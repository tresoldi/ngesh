---
title: 'Ngesh: a Python library for phylogenetic tree simulation'
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
simulating phylogenetic trees and data, particularly for historical linguistics. 
It provides simulated reproducible data for
the development, debugging, and benchmarking of
phylogenetic methods and tools, allowing simulations according to
different parameters and constrains and output in a variety of file formats.

# Background

Computational phylogenetics, especially after the popularization of
Bayesian inference, is being increasingly used in fields beyond biology,
such as historical linguistics and stemmatics (cite).
Such expansions has led to some uncertainties and resistence about its
adaptability, as the underlying evolutionary metaphors may not be
immediately transposable across fields, requiring its efficiency to
be demonstrated. Stochastic simulations, long recommended for
natural sciences (@Bailey:1990) and phylogenetics in general
(e.g. @Foote:1999, @Harmon:2019), becomes even more important, allowing
to empirically evaluate efficiency and accuracy, particularly in
humanity fields in which the number of documented evolutions to serve
as references may be limited. As a side-effect, simulations can serve as
an initial fuzzy testing of software pipelines, evaluating input, output,
parsing, and other functions.

The `ngesh` library is intended to allow such reproducible generation
according to user-provided seeds. Trees and related
data can be generated according to different parameters,
such as birth and death ratios, contrained in different ways, such as maximum
to a maximum branch length from root or number of extant nodes,
and manipulating output in different ways, such as pruning non-extant leaves
or simulating bad taxa sampling.
Character evolution corresponding to the topology can be likewise
generated, including user-specified ratios for mutation and horizontal
gene transfer. Taxa can be labeled either sequentially 
(such as "L01", "L02", "L03") or randomly with either human-readable
(like "Sume", "Fekobir", and "Tukok") or binominal biological nomenclature
(like "Sburas wioris", "Zurbata ceglaces", and "Spellis spusso").
The simulated trees can be used as ETE (@ete2016) tree objects or exported in a
variety of formats, such as Newick trees, ASCII representation, tabular
textual listings, etc.

# Installation, Usage, & Examples

The library can be installed with the standard `pip` tool for package
management, by issuing `pip install ngesh`.
Trees can be generated directly from the command line, defaulting to small
structures in Newick format:

```bash
$ ngesh
(Mamut:1.11985,(Koge:0.880823,(Rozkob:0.789548,(Meu:0.706601,(((Felbuh:0.189693,Kefa:0.189693)1:0.117347,((Epib:0.153782,Vugog:0.153782)1:0.0884745,Puluk:0.242256)1:0.0647836)1:0.0469885,Efam:0.354028)1:0.352573)1:0.0829465)1:0.0912757)1:0.23903);
```

The tool allows for both configuration files and command-line flags overriding
them. Here we generate the Nexus data for a reproducible Yule
tree (note the "12345" seed) with a birth ratio of 0.75, at least 5 leaves
with "human" labels, and 10 presence/absence characters.

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

The package is, however, designed to be used as a library.
The main functions for generation are `gen_tree()`, which returns a random tree
topology, and `add_characters()`, which simulates character evolution in a
provided tree. As they are separate tasks, it is possible to just generate
a random tree or to simulate character evolution in an user provided tree.

```python

In [1]: import ngesh

In [2]: tree = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels="human")

In [3]: print(tree)

      /-Butobfa
   /-|
  |  |   /-Defomze
  |   \-|
  |      \-Gegme
--|
  |      /-Bo
  |   /-|
  |  |   \-Peoni
   \-|
     |   /-Riuzo
      \-|
         \-Hoale

In [4]: tree = ngesh.add_characters(tree, 10, 3.0, 1.0)

In [5]: print(ngesh.tree2nexus(tree))
#NEXUS

begin data;
  dimensions ntax=7 nchar=15;
  format datatype=standard missing=? gap=-;
  matrix
Hoale      100111101101110
Butobfa    101011101110101
Defomze    101011110110101
Riuzo      100111101101110
Peoni      110011101110110
Bo         110011101110110
Gegme      101011101110101
  ;
end;
```

# Alternatives

The most complete alternative for simulating phylogenetic processes, even
though with no particular support for historical linguistics, is the
R package `TreeSim` by @Stadler:2011. Always in R, the `rtree()` function
of the `ape` package and the `birthdeath.tree()` one of the `geiger`
package might also be sufficient. In Python, code similar to `ngesh` and
which served as an initial inspiration is provided by @Noutahi:2017,
and for simpler simulations the `.populate()` method of the `Tree` class in
ETE can be used as well. In all languages, manual randomization of taxon
placement in existing cladograms is a well known alternative.

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
