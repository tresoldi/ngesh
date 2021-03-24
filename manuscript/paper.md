---
title: "Ngesh: a Python library for generating synthetic phylogenetic data"
tags:
  - Python
  - phylogenetics
  - random phylogenetic tree
  - phylogenetic tree simulation
  - synthetic data
authors:
  - name: Tiago Tresoldi
    orcid: 0000-0002-2863-1467
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Department of Linguistics and Philology, Uppsala University
   index: 1
 - name: Department of Linguistic and Cultural Evolution, Max Planck Institute for the Science of Human History
   index: 2
date: 02 April 2021
bibliography: paper.bib
---

# Summary

This work presents [`ngesh`](https://pypi.org/project/ngesh/), a Python library for simulating phylogenetic trees and
data,
designed for usage in development, debugging, and benchmarking of analysis pipelines and methods for phylogenetic
inference, particularly in historical linguistics and stemmatics.
The package generates reproducible stochastic simulations of evolution according to various criteria,
including character
mutation rates and probability of horizontal transfer, and its results can include the simulation of inadequate
data compilation and sampling. Different output formats are supported, both for visualization
(such as plain text and with integrated graphical viewers)
and for software interoperability (such as Newick and NEXUS).

# Background

Computational phylogenetics is being increasingly accepted in fields beyond biology, such as historical
linguistics [@Bouckaert:2012] and stemmatics [@Robinson:2016].
Stochastic simulations, long advocated for natural sciences in general [@Bailey:1990] and
genetics in specific [@Foote:1999; @Harmon:2019], are not yet used often enough in these
fields. However, they are highly desirable, allowing to
evaluate evolutionary analogies and models, as well as methods' performance,
through vast amounts of simulated histories, without limits imposed by
data availability and collection time and with quantifiable precision of results.
Simulations can also be used to perform basic fuzzy testing of software and to
support studies on which evolutionary models, processes, and evolutionary parameters
better match the observed phenomena.

The [`ngesh`](https://pypi.org/project/ngesh/) library is a tool that
allows to perform such simulations, designed for easy
integration into phylogenetic pipelines. It can generate reproducible trees and
correlated data from random seeds, following both user-established parameters,
such as ratios of birth and death, and constrains, such as branch
lengths and minimum number of taxa. The
library can label taxa progressive enumeration or with random names
that are easy to pronounce (e.g. "Sume" and "Fekobir") or which
imitate the binominal nomenclature (e.g. "Sburas wioris" and "Zurbata pusso").
Character evolution related to the
tree topology can likewise be simulated, including *ex novo* mutations and
horizontal gene transfers. Results can be manipulated in diverse
manners, for example by pruning extinct leaves or simulating uneven sampling. 
The simulated trees are standard ETE3 objects [@ETE:2016] and may
be exported into different formats such as Newick trees, ASCII-art representation,
and tabular lists.


As a building block for evaluating pipelines of analysis, the library
is a curated alternative to the basic technique of randomizing taxa
placement in existing cladograms, and to simpler tools such as the one by
@Noutahi:2017 or the `populate()` method of ETE3’s `Tree`
class [@ETE:2016]. Within its intended scope, it compares to
popular alternatives, including the R `TreeSim` [@Stadler:2011] and `geiger`
packages [@Pennell:2014] and the `rtree()` function of `ape` [@Paradis:2018],
however with specific support for historical linguistics and stemmatics.

# Installation, Usage, & Examples

Users can install the library with the standard `pip` tool for managing Python packages. Trees can be
generated from the command-line, defaulting to small phylogenies in Newick format:

```bash
$ ngesh
(Ukis:1.11985,(Koge:0.880823,(Rozkob:0.789548,(Meu:0.706601,
(((Felbuh:0.189693,Kefa:0.189693)1:0.117347,((Epib:0.153782,
Vugog:0.153782)1:0.0884745,Puluk:0.242256)1:0.0647836)1:0.0469885,
Efam:0.354028)1:0.352573)1:0.0829465)1:0.0912757)1:0.23903);
```

The tool supports both configuration files and command-line flags
that take precedence over the former. Here we specify a model to generate Nexus data
for a reproducible Yule tree, with a birth rate of 0.75, at least 5 leaves,
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

Despite the benefit of a stand-alone tool, the package is designed to be run as a library.
The two primary functions are `gen_tree()`, which returns a random tree, and
`add_characters()`, which adds character evolution data to a tree. Users can
generate random trees without character information or simulate character evolution within existing trees,
including non-simulated ones.

```python
>>> import ngesh
>>> tree = ngesh.gen_tree(1.0, 0.5, max_time=0.3, labels="bio",
                          seed="135")
>>> print(tree)

   /-Lubedsas larpes
--|
  |   /-Rasso wimapudda
   \-|
      \-Sbaes rapis
>>> print(tree.write())
(Lubedsas larpes:0.201311,(Rasso wimapudda:0.0894405,Sbaes rapis:0.0894405)
1:0.11187);
>>> tree = ngesh.add_characters(tree, 15, 2.0, 0.5)
```

Besides the `write()` method of the example above, which outputs Newick trees, results can be exported in NEXUS
format with the `tree2nexus()` function and in a tabular output expected by
tools such as BEASTling [@Maurits:2017], with the `tree2wordlist()` function.

# Code and Documentation Availability

The `ngesh` source code is available on GitHub at [https://github.com/tresoldi/ngesh](https://github.com/tresoldi/ngesh).

User documentation is available at [https://ngesh.readthedocs.io/](https://ngesh.readthedocs.io/).

# Acknowledgements

The author has received funding from the Riksbankens Jubileumsfond
(ID: MXM19-1087:1, ["Cultural Evolution of Texts"](https://www.rj.se/en/anslag/2019/cultural-evolution-of-texts/))
and from the European Research Council (ERC) under the European
Union’s Horizon 2020 research and innovation programme
(No. ERC #715618, ["Computer-Assisted Language Comparison"](https://digling.org/calc/)).

# References
