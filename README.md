# ngesh

`ngesh` is a Python library for generating random phylogenetic trees and related data.
It is intended for benchmarking phylogenetic methods and for providing dummy trees
during the development or debugging of phylogenetic methods. The generation of
random phylogenetic trees also goes by the name of "simulation methods for phylogenetic
trees" or just "simulating phylogenetic trees".

*Please note that this is both a work in progress and a library intended
to be simple, for complex methods see the bibliography and the alternatives
here listed*

In detail with ngesh:

* trees are returned either as strings in Newick representation or as
Python's [ETE](http://etetoolkit.org/) tree objects
* trees can be generated from user-specified seeds, so that the
generation can be reproduced
* trees can be generated according to user-provided birth and death ratios
* the death ratio can be set to zero
* non-extant leaves can be pruned
* speciation events default to two descendants, but the number of descendants
can be randomly drawn from a user-defined Poisson process (so that it
is possible to model hard politomies)
* trees can be limited in terms of number of extant leaves, evolution time
(as related to the birth and death parameters), or both
* nodes can optionally receive unique labels, either sequential ones,
random human-readable names, or random biological names following the
binomial nomenclature standard

## Installation

(lorem ipsum)

## How to use

(lorem ipsum)

## How to cite

(lorem ipsum)

## How does ngesh work?

For each tree, an `event_rate` is computed from the sum of the `birth` and
`death` rates. At each iteration, which takes places after an
random expovariant time from the `event_rate`, one of the extant nodes is
selected for an "event": either a birth or a death from the
proportion of each rate. All other extant leaves have their distances
updated with the event time.

The random labels follow the expected methods for random text generation
from a set of patterns, taking care to generate names as universally
readable (if not pronounceable) as possible.

## What does "ngesh" mean?

Technically it is just a name, but it was derived from one of the Sumerian words
for "tree", [ĝeš](http://psd.museum.upenn.edu/epsd/epsd/e2052.html), albeit
with an uncommon transcription. The name comes from the library once being
a module of a larger system for simulating language evolution and benchmarking
related tools, called [Enki](https://en.wikipedia.org/wiki/Enki) after the
Sumerian god of (among many other things) language and mischief.

The intended pronounciation, as in the most accepted reconstructions, is /ŋeʃ/. 
But don't strees over it: it is just a unique name.

## Alternatives

There are many tools for simulating phylogenetic processes in order to obtain
random phylogenetic trees. The most complete is probably the R package
[`TreeSim`](https://cran.r-project.org/web/packages/TreeSim/index.html)
by Tanja Stadler, which includes many flexible tree simulation functions. In
R, one can also use the `rtree` function from `ape` and the
`birthdeath.tree` from `geiger`, as well as manually randomizing taxon
placement in cladograms.

In Python, some code similar to ngesh and which served as initial inspiration
is provided by Marc-Rolland Noutahi on the blog post
[How to simulate a phylogenetic tree ? (part 1)](https://mrnoutahi.com/2017/12/05/How-to-simulate-a-tree/).

A number of on-line tools are also available at the time of writing:

https://lukejharmon.github.io/pcm/chapter10_birthdeath/#ref-Stadler2011-xu

* [T-Rex (Tree and reticulogram REConstruction](http://www.trex.uqam.ca/index.php?action=randomtreegenerator&project=trex)
at the Université du Québec à Montréal (UQAM)
* [Anvi'o Server](https://anvi-server.org/meren/random_phylogenetic_tree_w500_nodes) can
be used on-line as a wrapper to T-Rex above
* [phyloT](https://phylot.biobyte.de/), which by randomly sampling taxonomic names,
identifiers or protein accessions can be used for the same purpose

(add reference bibliography)

## References

https://lukejharmon.github.io/pcm/chapter10_birthdeath/

## Author

Tiago Tresoldi (tresoldi@shh.mpg.de)

The author was supported during development by the 
[ERC Grant #715618](https://cordis.europa.eu/project/rcn/206320/factsheet/en)
for the project [CALC](http://calc.digling.org)
(Computer-Assisted Language Comparison: Reconciling Computational and Classical
Approaches in Historical Linguistics).
