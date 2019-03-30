# ngesh

`ngesh` is a Python library for generating random phylogenetic trees and related data.
It is intended for benchmarking phylogenetic methods and for providing dummy trees
during the development or debugging of phylogenetic methods. The generation of
random phylogenetic trees also goes by the name of "simulation methods for phylogenetic
trees" or just "simulating phylogenetic trees".

In detail:

* trees are returned either as strings in Newick representation or as
Python's [ETE](http://etetoolkit.org/) tree objects
* trees can be generated from user-specified seeds, so that the
generation can be reproduced
* trees can be generated according to user-provided birth and death ratios
* the death ratio can be set to zero
* non-extant leaves can be pruned
* speciation events default to two descendants, but the number of descendants
can be randomly drawn from a user-defined Poisson process
* trees can be limited in terms of number of extant leaves, evolution time
(in relation to the birth death parameters), or both
* nodes can optionally receive unique labels, either sequential ones,
random human-readable names, or random biological names following the
binomial nomenclature standard.

In its simplest usage...

## How does ngesh work?

For each tree, an `event_rate` is computed from the sum of the `birth` and
`death` rates. At each iteration, which takes places after an
random expovariant time from the `event_rate`, one of the extant nodes is
selected for an "event", randomly either a birth or a death from the
proportion of their rates. All other extant leaves have their distances
ipdated with the event time.

The random labels follow the expected methods for random text generation
from a set of patterns, taking care to make names as universally
readable (if not pronounceable) as possible.

## What does "ngesh" mean?

Technically it is just a name, but it was derived from one of the Sumerian words
for "tree", [ĝeš](http://psd.museum.upenn.edu/epsd/epsd/e2052.html), albeit
with an uncommon transcription. The name comes from the library once being
a module of a larger system for simulating language evolution and benchmarking
related tools, called [Enki](https://en.wikipedia.org/wiki/Enki) after the
Sumerian god of (among many other things) language and mischief.

The intended pronountiation, as in the most accepted reconstructions, is /ŋeʃ/. 
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

* [T-Rex (Tree and reticulogram REConstruction](http://www.trex.uqam.ca/index.php?action=randomtreegenerator&project=trex)
at the Université du Québec à Montréal (UQAM)
* [phylot](https://phylot.biobyte.de/), which by randomly sampling taxonomic names,
identifiers or protein accessions can be used for the same purpose

## Author

Tiago Tresoldi <tresoldi@shh.mpg.de>

The author was supported by the 
[ERC Grant #715618](https://cordis.europa.eu/project/rcn/206320/factsheet/en)
for the project [CALC](http://calc.digling.org)
(Computer-Assisted Language Comparison: Reconciling Computational and Classical
Approaches in Historical Linguistics).
