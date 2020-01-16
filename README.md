# Ngesh, simulation of random phylogenetic trees with characters

[![Build Status](https://travis-ci.org/tresoldi/ngesh.svg?branch=master)](https://travis-ci.org/tresoldi/ngesh)
[![codecov](https://codecov.io/gh/tresoldi/ngesh/branch/master/graph/badge.svg)](https://codecov.io/gh/tresoldi/ngesh)
[![Codacy
Badge](https://api.codacy.com/project/badge/Grade/16ece2c98e3e4f319cb134bef2ade19c)](https://www.codacy.com/manual/tresoldi/ngesh?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=tresoldi/ngesh&amp;utm_campaign=Badge_Grade)
[![PyPI](https://img.shields.io/pypi/v/ngesh.svg)](https://pypi.org/project/ngesh)
[![DOI](https://zenodo.org/badge/178537103.svg)](https://zenodo.org/badge/latestdoi/178537103)

`ngesh` is a Python library for simulating random phylogenetic trees and
related data (characters, states, branch length, etc.).
It is intended for benchmarking phylogenetic methods, especially in
historical linguistics, and for providing dummy trees
for their development and debugging. The generation of
random phylogenetic trees also goes by the name of "simulation methods
for phylogenetic trees" or just "phylogenetic tree simulation".

*Please remember that `ngesh` is a work-in-progress and a library intended to
be a simple drop-in for cases where random trees are needed; for complex
methods, see the alternatives listed below or consult the
bibliographic references.*

In detail, with `ngesh`:

* trees can be generated according to user-specified birth and death ratios (and
the death ratio can be set to zero, resulting in a birth-only tree)
* speciation events default to two descendants, but the number of descendants
can be randomly drawn from a user-defined Poisson process (allowing
to model hard politomies)
* trees will have random topologies and, if necessary, random branch-lengths
* trees can be limited in terms of number of extant leaves, evolution time
(as related to the birth and death parameters), or both
* non-extant leaves can be pruned from birth-death trees
* character evolution can be simulated in relation to branch lengths,
with user-specified ratios for mutation and horizontal gene transfer,
with different rates of change for each character
* trees can be generated from user-provided seeds, so that the random
generation can be maintained across executions (and, in most cases, the
execution should be reproducible *also* on different machines and different
vestions of Python)
* nodes can optionally receive unique labels, either sequential ones
(like "L01", "L02", and "L03"),
random human-readable names (like "Sume", "Fekobir", and "Tukok"),
or random biological names approximating the
binomial nomenclature standard (like "Sburas wioris", "Zurbata ceglaces",
and "Spellis spusso")
* trees can be returned as [ETE](http://etetoolkit.org/) tree objects or
exported in a variety of formats, such as Newick trees, ASCII representation,
tabular textual listings, etc.

### Changelog

Version 0.3.1:
 - Code improvements for next release and submission

Version 0.3:

- General improvements to code quality
- Full reproducibility from seeds for the pseudo-random generators,
  allowing string, ints, and floats
- Changes for further integration with `abzu` and `alteruphono` for
  simulating linguistic data

## How does ngesh work?

For each tree, an `event_rate` is computed from the sum of the `birth` and
`death` rates. At each iteration, which takes place after an
random expovariant time from the `event_rate`, one of the extant nodes is
selected for an "event": either a birth or a death from the
proportion of each rate. All other extant leaves have their distances
updated with the event time.

The random labels follow the expected methods for random text generation
from a set of patterns, taking care to generate names as universally
readable (if not pronounceable) as possible.

*missing on character generation*

## Installation

In any standard Python environment, `ngesh` can be installed with:

```
pip install ngesh
```

The `pip` installation will also fetch the dependencies `ete3` and
`numpy`, if necessary.

## How to use

You can test your installation from the command line with the `ngesh` command, which
will return a different random small birth-death tree in Newick format each time it
is called:

```
$ ngesh
(Saorus getes:1.31562,((Voces earas:1.07567,(Dallao spettus:0.703609,Sburas wioris:0.703609)1:0.372063)1:0.464667,(Zurbaza ceglaces:0.527431,(Amduo vizoris:0.345862,Uras wiurus:0.345862)1:0.18551)1:1.00897)1:2.1707);

$ ngesh
((Ollio zavis:0.698453,(Spectuo sicui:0.596731,((Ronis mivulis:0.0431014,Vaporus conomattas:0.0431014)1:0.413634,Rizarus urrus:0.456735)1:0.139996)1:0.101722)1:3.17827,(Deses mepus:2.22061,(Ovegpuves wiumoras:1.88469,(Easas ecdebus:0.201891,Muggas lupas:0.201891)1:1.6828)1:0.335918)1:1.65611);
```

The same command line tool can refer to values provided in a textual
configuration file. Here, we generate the Nexus data for a
reproducible Yule tree (note the `12345` seed)
with a birth ratio of 0.75, at least 8 leaves with `"human"` labels,
and 10 presence/absence characters:

```
$ cat my_tree.conf 
[Config]
labels=human
birth=0.666
death=0.0
output=nexus
min_leaves=8
num_chars=10

$ ngesh -c mine.conf --seed 12345
#NEXUS

begin data;
  dimensions ntax=15 nchar=25;
  format datatype=standard missing=? gap=-;
  matrix
Foro        1011000101010010100100100
Meno        1100100101010010010010001
Vuea        1100110001010100001001001
Vegevo      1100100101010010010010100
Bufuri      1100110001010100001001001
Novake      1100110001010100001001001
Fonulip     1100110001010100001001001
Omih        1101001001010011000100100
Onegro      1101000011010010001100100
Rolsoa      1100100100110010010100100
Wigu        1101001001010010001100100
Teozu       1101001001010011000100010
Kabu        1100100101001010010100100
Timebbed    1100100101010010010010100
Okuna       1100110001010100001001001
  ;
end;
```

Parameters set in a configuration file can be overridden at the command line.
The ASCII representation of the topology of the same tree can be obtained
with:

```
$ ngesh -c mine.conf --seed 12345 -o ascii

         /-Omih
        |
      /-|      /-Onegro
     |  |   /-|
     |   \-|   \-Wigu
     |     |
     |      \-Teozu
     |
   /-|   /-Kabu
  |  |  |
  |  |  |            /-Novake
  |  |  |           |
  |  |  |         /-|      /-Okuna
  |  |  |        |  |   /-|
  |   \-|        |   \-|   \-Fonulip
  |     |      /-|     |
  |     |     |  |      \-Bufuri
  |     |     |  |
--|     |   /-|   \-Vuea
  |     |  |  |
  |     |  |  |   /-Meno
  |      \-|   \-|
  |        |      \-Rolsoa
  |        |
  |        |   /-Vegevo
  |         \-|
  |            \-Timebbed
  |
   \-Foro
```

The package is, however, designed to be used as a library. If you have
PyQt5 installed (which is not listed as a dependency), the following code
will pop up the ETE Tree Viewer on a random tree:

```
python3 -c "import ngesh ; ngesh.display_random_tree()"
```

![random tree](https://raw.githubusercontent.com/tresoldi/ngesh/master/doc/tree001.png){width=100%}

The main functions for generation are `gen_tree()`, which returns a random
tree topology, and `add_characters()`, which simulates character evolution
in a provided tree. As they are separate tasks, it is possible to just
generate a random tree or to simulate character evolution in an user
provided tree.

The full documentation is available in the functions docstring (which
can be visualized with `print(ngesh.gen_tree.__doc__)` and
`print(ngesh.add_characters.__doc__)`) or
[directly in the source code](https://github.com/tresoldi/ngesh/blob/master/ngesh/random_tree.py).
The code snipped below shows a basic tree generation, character evolution,
and output flow; the parameters for generation are the same listed in
the docstrings and in the following below.

```
$ ipython3

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

### Integrating with other software

Integration is easy due to the various export functions. For example, it
is possible to generate random trees with characters for which we know
all details on evolution and parameters, and generate Nexus files that
can be fed to phylogenetic software such as
[MrBayes](http://nbisweden.github.io/MrBayes/) or
[BEAST2](https://www.beast2.org/)
to either check how they perform or
how good is our generation in terms of real data.

Let's simulate phylogenetic data for an analysis using BEAST2 through
[BEASTling](https://github.com/lmaurits/BEASTling). We start with
a birth-death tree (lambda=0.9, mu=0.3), with at least 15 leaves, and 100
characters whose evolution is modelled with the default parameters
and a string seed `"jena"` for reproducibility; the tree data is exported
in `"wordlist"` format:

```
$ cat examples/example_ngesh.conf 
[Config]
labels=human
birth=0.9
death=0.3
output=nexus
min_leaves=15
num_chars=100

$ ngesh -c examples/example.conf --seed jena > examples/example.csv

$ head examples/example.csv 
Language_ID,Feature_ID,Value
Dotare,feature_0,0
Dotare,feature_1,1
Dotare,feature_2,2
Dotare,feature_3,3
Dotare,feature_4,4
Dotare,feature_5,5
Dotare,feature_6,6
Dotare,feature_7,7
Dotare,feature_8,8
```

We can now use a minimal BEASTling configuration and generate an XML
input for BEAST2. Let's assume we want to test how well our pipeline
performs when assuming a Yule tree when the data actually includes
extinct taxa. The results here presented are not expected to perfect,
as we will use
a short chain length to make it faster and a model which is different
from the assumptions used for generation (besides the fact of the
default parameters for
horizontal gene transfer being a bit too aggressive).

```
$ cat examples/example_beastling.conf 
[admin]
basename=example
[MCMC]
chainlength=500000
[model example]
model=covarion
data=example.csv

$ beastling example_beastling.conf 

$ beast example.xml 
```

We can proceed normally here: use BEAST2's `treeannotator` (or similar
software) to generate a summary tree,
which we store in `examples/summary.nex`,
and plot the results with `figtree` (or, again, similar software).

Let's plot our summary tree and compare the results with the
actual topology (which we can regenerate with the earlier seed).

![summary tree](https://raw.githubusercontent.com/tresoldi/ngesh/master/doc/summary.nex.png){width=100%}

```
$ ngesh -c examples/example_ngesh.conf --seed jena --output newick > examples/example.nw
```

![original tree](https://raw.githubusercontent.com/tresoldi/ngesh/master/doc/example.nw.png){width=100%}

The results are not excellent given the limits we set for quick demonstration,
but it still capture major information and subgroupings (as clearer by
the radial layout below) -- manual data exploration show that at least some
of the errors, including the group in the first split, are due to horizontal
gene transfer. For an analysis of
the inference performance we would need to improve the parameters above
and repeat the analysis on a range of random trees, including studying the
log of character changes (including borrowings) involved in this particular
random tree.

![summary tree radial](https://raw.githubusercontent.com/tresoldi/ngesh/master/doc/summary.nex2.png){width=100%}

*TODO: Compare trees (Robinson-Foulds symmetric difference?)*

The files used and generated in this example
can be found in the `/examples` directory.

## Parameters for tree generation

The parameters for tree generation, as also given by `ngesh -h`, are:

* `birth`: The tree birth rate (l)
* `death`: The tree death rate (mu)
* `max_time`: The stopping criterion for maximum evolution time
* `min_leaves`: The stopping criterion for minimum number of leaves
* `labels`: The model for textual generation of random labels
(`None`, `"enum"` for a simple enumeration, `"human"` for randomly
generated names, and `"bio"` for randomly generated specie names)
* `num_chars`: The number of characters to be simulated
* `k_mut`: The character mutation gamma `k` parameter
* `th_mut`: The character mutation gamma `th` parameter
* `k_hgt`: The character HGT gamma `k` parameter
* `th_hgt`: The character HGT gamma `th` parameter
* `e`: The character general mutation `e` parameter

## What does "ngesh" mean?

Technically it is just an unique name, but it was derived from one of the Sumerian words
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
[`TreeSim`](https://CRAN.R-project.org/package=TreeSim)
by Tanja Stadler, which includes many flexible tree simulation functions. In
R, one can also use the `rtree()` function from package `ape` and the
`birthdeath.tree()` one from package `geiger`, as well as manually randomizing taxon
placement in cladograms.

In Python, some code similar to `ngesh` and which served as initial inspiration
is provided by Marc-Rolland Noutahi on the blog post
[How to simulate a phylogenetic tree ? (part 1)](https://mrnoutahi.com/2017/12/05/How-to-simulate-a-tree/).

For simpler simulations, the `.populate()` method of the `Tree` class
in ETE might be enough as well. Documentation on the method is
available
[here](http://etetoolkit.org/docs/latest/reference/reference_tree.html#ete3.TreeNode.populate).

A number of on-line tools are also available at the time of writing:

* [T-Rex (Tree and reticulogram REConstruction](http://www.trex.uqam.ca/index.php?action=randomtreegenerator&project=trex)
at the Université du Québec à Montréal (UQAM)
* [Anvi'o Server](https://anvi-server.org/meren/random_phylogenetic_tree_w500_nodes) can
be used on-line as a wrapper to T-Rex above
* [phyloT](https://phylot.biobyte.de/), which by randomly sampling taxonomic names,
identifiers or protein accessions can be used for the same purpose

## TODO

* Shorter-term
    * Write better documentation of function parameters
    * Add all still unavailable parameters to the command line tool (e.g.,
      setting hard politomies)
    * Automatically generate developer documentation (possibly with Sphinx)
    * Allow generation of unlabelled trees from the command-line (a text
      generation model is currently mandatory)
    * Look for default parameters that are more closely related to
      linguistic trees (or, at least, to the Indo-European one)
    * Check if all outputs are complete (e.g., characters are currently
      missing in the Newick format)
    * Add a command-line option (or a new tool) that allows to write the
      output to one file and the reference tree to a second one (possibly
      with the log of character evolution)
    * Check for alternatives for the exponential correction of a character
      resistance to mutation (e.g. Zipf law), including separating mutation
      and borrowing rates
    * Allow to exclude non extant taxa from horizontal gene transfer
      events
    * Add stopping criterion on the global number of nodes (in complement
      to the number of *extant* nodes, currently implemented), either absolute
      or as a range

* Longer-term
    * Simulation of data problems (incomplete sampling, errors in
      sequencing/cognate judgment, etc.)
    * Variable birth/death ratios
    * Rewrite the random text generation functions, possibly as actual
      Python generators
    * Consider replacing or complementing `expovariate()` in birth/death
      events with actual random Poisson sampling, allowing additional models
    * Build a simple website
    * Implement parallel character evolution as controlled by a parameter
    * Rewrite functions with too many arguments to accept dictionaries of
      parameters
    * Implement more models for random character generation, especially those
      frm genetics (first candidate, a General Time Reversable model with
      a proportion of invariable sites and a gamma-shaped distribution of
      rates across sites)
    * Simulate a "donation power" for taxa, making borrowing events globally
      more likely from a given donor (analogous to cultural influence in
      linguistics)
    * Allow to guarantee that borrowing events will always result in
      altered states (it is currently possible that an event will borrow
      an equal state for a given character, especially considering that we
      favor borrowing from closer taxa)
    * Implement character simulation for other datatypes, particularly from
      genetics (currently only standard binary presence/absence)

## Gallery

![random tree](https://raw.githubusercontent.com/tresoldi/ngesh/master/doc/tree001.png){width=100%}
![random tree](https://raw.githubusercontent.com/tresoldi/ngesh/master/doc/tree002.png){width=100%}
![random tree](https://raw.githubusercontent.com/tresoldi/ngesh/master/doc/tree003.png){width=100%}

## How to cite

If you use `ngesh`, please cite it as:

> Tresoldi, Tiago (2019). Ngesh, a tool for simulating random phylogenetic trees.
Version 0.3. Jena. Available at: https://github.com/tresoldi/ngesh

In BibTeX:

```
@misc{Tresoldi2019ngesh,
  author = {Tresoldi, Tiago},
  title = {Ngesh, a tool for simulating random phylogenetic trees. Version 0.3},
  howpublished = {\url{https://github.com/tresoldi/ngesh}},
  address = {Jena},
  year = {2019},
  doi = {10.5281/zenodo.2619311},
}
```

## References

* Bailey, N. T. J. (1964). *The elements of stochastic processes with applications to the natural sciences*. John Wiley & Sons.

* Foote, M., J. P. Hunter, C. M. Janis, and J. J. Sepkoski Jr. (1999). *Evolutionary and preservational constraints on origins of biologic groups: Divergence times of eutherian mammals*. Science 283:1310–1314.

* Harmon, Luke J (2019). *Phylogenetic Comparative Methods -- learning from trees*.
Available at: [https://lukejharmon.github.io/pcm/chapter10_birthdeath/](https://lukejharmon.github.io/pcm/chapter10_birthdeath/). Access date: 2019-03-31.

* Noutahi, Marc-Rolland (2017). *How to simulate a phylogenetic tree? (part 1)*. Available at:
[https://mrnoutahi.com/2017/12/05/How-to-simulate-a-tree/](https://mrnoutahi.com/2017/12/05/How-to-simulate-a-tree/). Access date: 2019-03-31

* Stadler, Tanja (2011). *Simulating Trees with a Fixed Number of Extant Species*. Systematic Biology 60.5:676-684. DOI: [https://doi.org/10.1093/sysbio/syr029](https://doi.org/10.1093/sysbio/syr029)

## Author

Tiago Tresoldi (tresoldi@shh.mpg.de)

The author was supported during development by the 
[ERC Grant #715618](https://cordis.europa.eu/project/rcn/206320/factsheet/en)
for the project [CALC](http://calc.digling.org)
(Computer-Assisted Language Comparison: Reconciling Computational and Classical
Approaches in Historical Linguistics), led by
[Johann-Mattis List](http://www.lingulist.de).
