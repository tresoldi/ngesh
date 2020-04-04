---
title: 'ngesh: a Python library for phylogenetic tree simulation'
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

This work describes [`ngesh`](https://pypi.org/project/ngesh/), a Python library
for the task of phylogenetic tree and data simulation.
It is intended to provide dummy
trees and data for the development, debugging, and benchmarking of
phylogenetic methods and tools, particularly in the field of historical
linguistics. The library allows to reproducibly generate random trees
according to user-specified parameters, such as birth and death ratio, and
withing certain constraints, such as relative time depth or number or
extant nodes, allowing to output either as programming variables or
in a variety of file formats.

# Background

Métodos estocásticos têm sido utilizado em uma variedade campos além daquele
biológico para o qual foram inicialmente desenvolvidos, como a linguística
histórica (citar). A adaptaćão de tais métodos tem levado a algumas
incertezas sobre sua adaptabilidade, pois as metáforas evolutivas básicas
podem não ser imediatamente transponíveis; deste modo, faz-se ainda mais
importante a prática de simulaćões estocásticas que permitam desenvolver e
avaliar a eficiências e eficácia de tais métodos, desde a verificaćão
de métodos de input e output, como o parsing de formatos ou a geraćão de
visualizaćões gráficas, avaliando o desempenho de pipelines de modo
independete da análise real, até a verifićão da robustez dos métodos,
como na equivalente prática de desenvolvimento de software de fuzzy
testing. Para estas finalidades é útil a random phylogenetic tree generation,
which also goes by the names of "phylogenetic tree simulation".
Just to cite @Bailey:1990, @Foote:1999, @Harmon:2019, @Stadler:2011,
@Noutahi:2017 

O ngesh é uma biblioteca e ferramentas de linha de comando que permitem
tal geraćão. Trees can be generated according to different parameters
such as birth and death ratio (allowing birth-only trees if desired),
contrained according to different parameters (such as maximum branch length
from root or number of extant nodes),
output distribution for speciation events (allowing to model hard politomies),
pruning of non-extant leaves, etc.
Character evolution corresponding to the random topology can be likewise
generated, including user-specified ratios for mutation and horizontal
gene transfer (including different rates of change for each character).
For usability, trees can be generated according to user-provided seeds and
can optionally receive unique labels, being either sequences (such as
"L01", "L02", "L03", random human-readable names (like "Sume", "Fekobir",
and "Tukok"), or random biological names approximating the binomial nomenclature
standard (like "Sburas wioris", "Zurbata ceglaces", and "Spellis spusso").
The simulated trees can be used as ETE tree objects or exported in a
variety of formats, such as Newick trees, ASCII representation, tabular
textual listings, etc.

# Installation, Usage, & Examples

The library can be installed with the standard `pip` tool for package
management:

```bash
$ pip install ngesh
```

The [documentation](https://ngesh.readthedocs.io/en/latest/) offers detailed
instructions on how to use the library. Simpler trees can be generated
directly from the command line, which by default will return a different
random small birth-death tree in Newick format at each call:

```bash
$ ngesh
(output)
```

The same command line tool can refer to values provided in a textual
configuration file. Here, we generate the Nexus data for a reproducible Yule
tree (note the 12345 seed) with a birth ratio of 0.75, at least 8 leaves
with "human" labels, and 10 presence/absence characters:
Parameters set in a configuration file can be overridden at the command line.
As ASCII representation of the same tree could be obtained with the
`--ascii` flag.

```bash
$ cat
```

The package is, however, designed to be used as a library. If you have PyQt5
installed (which is not listed as a dependency), the following code will pop
up the ETE Tree Viewer on a random tree:

```bash
$ python -c "import ngesh ; ngesh.display_random_tree()"
```

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

Particularly for specific biological simulations, a number of on-line tools can
also be used, such as 
`T-Rex` (Tree and reticulogram REConstruction at the Université du Québec à Montréal (UQAM),
`Anvi'o Server` can be used on-line as a wrapper to T-Rex above,
and `phyloT`, which by randomly sampling taxonomic names, identifiers or
protein accessions can be used for the same purpose.


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
