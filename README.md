# Ngesh, a library for phylogenetic tree simulation

[![PyPI](https://img.shields.io/pypi/v/ngesh.svg)](https://pypi.org/project/ngesh)
[![CI](https://github.com/tresoldi/ngesh/actions/workflows/CI.yml/badge.svg)](https://github.com/tresoldi/ngesh/actions/workflows/CI.yml)
[![Codacy
Badge](https://api.codacy.com/project/badge/Grade/16ece2c98e3e4f319cb134bef2ade19c)](https://www.codacy.com/manual/tresoldi/ngesh?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=tresoldi/ngesh&amp;utm_campaign=Badge_Grade)
[![Documentation Status](https://readthedocs.org/projects/ngesh/badge/?version=latest)](https://ngesh.readthedocs.io/en/latest/?badge=latest)

`ngesh` is a Python library and command-line tool
for simulating phylogenetic trees and related data (characters, states,
branch length, etc.).
It is intended for benchmarking phylogenetic methods, especially in
historical linguistics and stemmatology. The generation of
stochastic phylogenetic trees also goes by the name of "simulation methods
for phylogenetic trees", "synthetic data", or just "phylogenetic tree simulation".

![ngesh](https://raw.githubusercontent.com/tresoldi/ngesh/master/docs/banner.png)

Among the highlights of the package, with `ngesh`:

* any hashable element can be provided as a seed for the pseudo-random number
  generators, guaranteeing that the generated trees are reproducible (including
  across different systems)
* trees can be generated according to user-specified parameters such as birth and death ratios (and
  the death ratio can be set to zero, resulting in a birth-only tree)
* trees will have random topologies and, if necessary, random branch-lengths
* trees can be limited in terms of number of extant leaves, evolution time
  (as related to the birth and death parameters), or both
* non-extant leaves can be pruned from birth-death trees
* speciation events default to two descendants, but the number of descendants
  can be randomly drawn from a user-defined Poisson process (allowing
  modelling of hard politomies)
* character evolution can be simulated in relation to branch lengths,
  with user-specified ratios for mutation and for horizontal gene transfer,
  with different rates of change for each character
* nodes can optionally receive unique labels, either sequential ones
  (like "L01", "L02", and "L03"), random pronounceable names (like "Sume", "Fekobir", and "Tukok"),
  or random biological names approximating the binomial nomenclature standard
  (like "Sburas wioris", "Zurbata ceglaces", and "Spellis spusso")
* trees can be returned as [ETE](http://etetoolkit.org/) tree objects or
  exported in a variety of formats, such as Newick trees, ASCII representation,
  tabular textual listings, etc.

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

```bash
$ ngesh
((Vovrera:0.149348,(Wigag:3.11592,(Pallo:2.68125,Zoei:1.85803)1:1.29704)1:0.204529)1:0.607805,(((Avi:0.347942,Uemi:0.0137646)1:1.41697,(((Kufo:0.817012,
(Gapurem:0.0203582,Hukub:0.0203582)1:0.796654)1:0.395727,Tablo:0.00846148)1:0.484705,(Kaza:0.140656,((Tozea:0.240634,Pebigmom:0.240634)1:1.13579,(Kata:0
.109977,((Fabom:0.04242,Upik:0.04242)1:0.549364,(Amue:0.182635,Lunida:0.182635)1:0.409149)1:0.366701)1:0.417941)1:0.162968)1:0.158051)1:1.47281)1:1.0326
,(Kunizob:0.650455,Madku:0.221172)1:1.22008)1:0.587783);


$ ngesh
((((Povi:0.325601,Udo:0.325601)1:0.0750448,Hiruta:0.400646)1:0.181454,(Voebi:0.0293506,Sodi:0.0293506)1:0.55275)1:0.258834,((Vandemif:0.0160558,(((Dubik
:0.0543122,Fuvu:0.0543122)1:0.36458,Hitfuv:0.418892)1:0.0388987,Pizuna:0.457791)1:0.0535386)1:0.179893,(Uo:0.67132,Zegna:0.163427)1:0.0199021)1:0.149711
);
```

The same command line tool can refer to values provided in a textual
configuration file. Here, we generate the Nexus data for a
reproducible Yule tree (note the `12345` seed)
with a birth ratio of 0.75, at least 8 leaves with `"human"` labels,
and 10 presence/absence characters:

```bash
$ cat ngesh_demo.conf
[Config]
labels=human
birth=0.666
death=0.0
output=nexus
min_leaves=8
num_chars=10

$ ngesh -c ngesh_demo.conf --seed 123
#NEXUS

begin data;
  dimensions ntax=16 nchar=38;
  format datatype=standard missing=? gap=-;
  matrix
Abel        10001001011000010000010010010000100000
Azogu       10001001011000010000010010010000100000
Bou         10001001100010100000010010010000000010
Dipu        10001001010001000010000110010000000001
Gezepsem    10001001100010100000010010010000000010
Gupote      10001001010010010000010010010000000100
Hefi        10100100010010010001000001010001000000
Lerzo       10001001010001000010000110010000000001
Magumel     10001001010010010000010010010000000010
Pao         01001010010100001000100010001000100000
Sanigo      10010100010010000100001000100010010000
Tuzizo      10001001100010100000010010010000000010
Wialum      10001001011000010000010010000100100000
Zudal       10001001010010010000010010010000100000
Zukar       10001001011000010000010010000100100000
Zusu        10010100010010000100001000100010001000
  ;
end;
```

Parameters set in a configuration file can be overridden at the command line.
The ASCII representation of the topology of the same tree can be obtained
with:

```bash
$ ngesh -c ngesh_demo.conf --seed 123 -o ascii

         /-Zudal
        |
        |               /-Azogu
        |              |
        |            /-|      /-Wialum
        |           |  |   /-|
        |           |   \-|   \-Zukar
        |         /-|     |
        |        |  |      \-Abel
        |        |  |
      /-|        |  |   /-Dipu
     |  |        |   \-|
     |  |      /-|      \-Lerzo
     |  |     |  |
     |  |     |  |         /-Bou
     |  |     |  |      /-|
     |  |     |  |   /-|   \-Gezepsem
     |  |   /-|  |  |  |
   /-|  |  |  |   \-|   \-Tuzizo
  |  |  |  |  |     |
  |  |   \-|  |      \-Magumel
  |  |     |  |
  |  |     |   \-Pao
  |  |     |
--|  |      \-Gupote
  |  |
  |  |   /-Zusu
  |   \-|
  |      \-Sanigo
  |
   \-Hefi
```

The package is, however, designed to be used as a library. If you have
PyQt5 installed (which is *not* listed as a dependency and must be installed
separately), the following code will pop up the ETE Tree Viewer on a random tree:

```bash
python3 -c "import ngesh ; ngesh.show_random_tree()"
```

![random tree](https://raw.githubusercontent.com/tresoldi/ngesh/master/docs/tree001.png)

The main functions for generation are `gen_tree()` ([doc](https://ngesh.readthedocs.io/en/latest/source/ngesh.html#ngesh.random_tree.gen_tree)),
which returns a random
tree topology, and
`add_characters()` ([doc](https://ngesh.readthedocs.io/en/latest/source/ngesh.html#ngesh.random_tree.add_characters)),
which simulates character evolution
in a provided tree. As they are separate tasks, it is possible to just
generate a random tree or to simulate character evolution in an user
provided tree.

The code snipped below shows a basic tree generation, character evolution,
and output flow; the parameters for generation are the same listed in
the docstrings and in the following below.

```python
>>> import ngesh
>>> tree = ngesh.gen_tree(1.0, 0.5, max_time=3.0, labels="human")
>>> print(tree)

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

>>> tree = ngesh.add_characters(tree, 10, 3.0, 1.0)
>>> print(ngesh.tree2nexus(tree))
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

## Parameters for tree generation

The parameters for tree generation, as also given by the command `ngesh -h`, are:

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
and a string seed `"uppsala"` for reproducibility; the tree data is exported
in `"wordlist"` format:

```bash
$ cat examples/example_ngesh.conf
[Config]
labels=human
birth=0.9
death=0.3
output=nexus
min_leaves=15
num_chars=100

$ ngesh -c examples/example_ngesh.conf --seed uppsala > examples/example.csv

$ head -n 20 examples/example.csv
Language_ID,Feature_ID,Value
Akup,feature_0,0
Buter,feature_0,0
Dufou,feature_0,0
Emot,feature_0,0
Kiu,feature_0,0
Kovala,feature_0,0
Lusei,feature_0,0
Oso,feature_0,0
Puota,feature_0,0
Relenin,feature_0,976
Sotok,feature_0,0
Tetosur,feature_0,0
Usimi,feature_0,976
Voe,feature_0,0
Vusodur,feature_0,0
Zeba,feature_0,0
Zufe,feature_0,0
Akup,feature_1,1
Buter,feature_1,1
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

```bash
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

![summary tree](https://raw.githubusercontent.com/tresoldi/ngesh/master/docs/summary.nex.png)

```bash
$ ngesh -c examples/example_ngesh.conf --seed uppsala --output newick > examples/example.nw
```

![original tree](https://raw.githubusercontent.com/tresoldi/ngesh/master/docs/example.nw.png)

The results are not excellent given the limits we set for quick demonstration,
but it still capture major information and subgroupings (as clearer by
the radial layout below) -- manual data exploration show that at least some
of the errors, including the group in the first split, are due to horizontal
gene transfer. For an analysis of
the inference performance we would need to improve the parameters above
and repeat the analysis on a range of random trees, including studying the
log of character changes (including borrowings) involved in this particular
random tree.

![summary tree radial](https://raw.githubusercontent.com/tresoldi/ngesh/master/docs/summary.nex2.png)

Trees can, as expected, be compared with common methods of tree
comparison, such as [Robinson–Foulds metric](https://en.wikipedia.org/wiki/Robinson%E2%80%93Foulds_metric).
All packages and programming languages for this purpose should be
able to read the trees exported in Newick or NEXUS format; however,
as `ngesh` trees are actually ETE3 trees, we can also do it directly
from Python:

```python
d = tree1.robinson_foulds(tree_2)
```

The files used and generated in this example can be found in the `/examples` directory.

## What does "ngesh" mean?

Technically it is just an unique name, but it was originally derived from one of the Sumerian words
for "tree", [ĝeš](http://psd.museum.upenn.edu/epsd/epsd/e2052.html). The name
was chosen because the library was first planned as part of
a larger system for simulating language evolution and benchmarking
related tools, named [Enki](https://en.wikipedia.org/wiki/Enki) after the
Sumerian god of (among many other things) language and "randomness".

The intended pronunciation, as in the most accepted reconstructions, is /ŋeʃ/.
But don't stress over it, and feel free to call it /n̩.gɛʃ/, as
most people seem to do: it is just a unique name.

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
The `toytree` and `dendropy` packages also offer comparable functionality.

A number of on-line tools are also available at the time of writing:

* [T-Rex (Tree and reticulogram REConstruction](http://www.trex.uqam.ca/index.php?action=randomtreegenerator&project=trex)
at the Université du Québec à Montréal (UQAM)
* [Anvi'o Server](https://anvi-server.org/meren/random_phylogenetic_tree_w500_nodes) can
be used on-line as a wrapper to T-Rex above
* [phyloT](https://phylot.biobyte.de/), which by randomly sampling taxonomic names,
identifiers or protein accessions can be used for the same purpose

## Gallery

![random tree](https://raw.githubusercontent.com/tresoldi/ngesh/master/docs/tree001.png)
![random tree](https://raw.githubusercontent.com/tresoldi/ngesh/master/docs/tree002.png)
![random tree](https://raw.githubusercontent.com/tresoldi/ngesh/master/docs/tree003.png)

## References

* Bailey, N. T. J. (1964). *The elements of stochastic processes with applications to the natural sciences*. John Wiley & Sons.

* Foote, M., J. P. Hunter, C. M. Janis, and J. J. Sepkoski Jr. (1999). *Evolutionary and preservational constraints on origins of biologic groups: Divergence times of eutherian mammals*. Science 283:1310–1314.

* Harmon, Luke J (2019). *Phylogenetic Comparative Methods -- learning from trees*.
Available at: [https://lukejharmon.github.io/pcm/chapter10_birthdeath/](https://lukejharmon.github.io/pcm/chapter10_birthdeath/). Access date: 2019-03-31.

* Noutahi, Marc-Rolland (2017). *How to simulate a phylogenetic tree? (part 1)*. Available at:
[https://mrnoutahi.com/2017/12/05/How-to-simulate-a-tree/](https://mrnoutahi.com/2017/12/05/How-to-simulate-a-tree/). Access date: 2019-03-31

* Stadler, Tanja (2011). *Simulating Trees with a Fixed Number of Extant Species*. Systematic Biology 60.5:676-684. DOI: [https://doi.org/10.1093/sysbio/syr029](https://doi.org/10.1093/sysbio/syr029)

The `ngesh` banner was designed by Tiago Tresoldi on basis of the
vignette "Sherwood Forest" by J. Needham
published in Needham, J. (1895) *Studies of trees in pencil and in water colors*. First
series. London, Glasgow, Edinburgh: Blackie & Son. (under public domain and
available on [archive.org](https://archive.org/details/studiesoftreesin00need/page/n3/mode/2up)).

## Community guidelines

While the author can be contacted directly for support, it is recommended that
third parties use GitHub standard features, such as issues and pull requests, to
contribute, report problems, or seek support.

Contributing guidelines, including a code of conduct, can be found in the
`CONTRIBUTING.md` file.

## Author and citation

The library is developed by Tiago Tresoldi (tiago.tresoldi@lingfil.uu.se). The library is developed in the context of
the [Cultural Evolution of Texts](https://github.com/evotext/), with funding from the
[Riksbankens Jubileumsfond](https://www.rj.se/) (grant agreement ID:
[MXM19-1087:1](https://www.rj.se/en/anslag/2019/cultural-evolution-of-texts/)).

During the first stages of development, the author received funding from the
[European Research Council](https://erc.europa.eu/) (ERC) under the European Union’s Horizon 2020
research and innovation programme (grant agreement
No. [ERC Grant #715618](https://cordis.europa.eu/project/rcn/206320/factsheet/en),
[Computer-Assisted Language Comparison](https://digling.org/calc/)).

If you use `ngesh`, please cite it as:

> Tresoldi, Tiago (2021). Ngesh, a tool for simulating random phylogenetic trees.
Version 0.5. Uppsala: Uppsala universitet. Available at: https://github.com/tresoldi/ngesh

In BibTeX:

```
@misc{Tresoldi2020ngesh,
  author = {Tresoldi, Tiago},
  title = {Ngesh, a tool for simulating random phylogenetic trees. Version 0.5},
  howpublished = {\url{https://github.com/tresoldi/ngesh}},
  address = {Uppsala},
  publisher = {Uppsala universitet},
  year = {2021},
}
```
