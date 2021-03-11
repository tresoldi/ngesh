User guide
==========

The ``ngesh`` library is designed to generate simulated data in
reproducible manner: the pseudo-random generator seeds are used not only
in a way that guarantees reproducibility (that is, the same seed
returning the same tree) but also in an incremental way that will return
similar trees when different parameters are used to the same seed.

A pure random tree can be generated with a call to the
``ngesh.gen_tree()`` function. Either a minimum number of leaves
(``min_leaves``) or maximum time (``max_time``) must be specified as a
stopping criterion.

.. code:: python

   >>> import ngesh
   >>> tree1 = ngesh.gen_tree(min_leaves=5)
   >>> print(tree1)

         /-L1
      /-|
     |   \-L2
   --|
     |   /-L3
      \-|
        |   /-L4
         \-|
            \-L5

   >>> tree2 = ngesh.gen_tree(max_time=1.0)
   >>> print(tree2)

         /-L1
      /-|
   --|   \-L2
     |
      \-L3

The main parameters are the ``birth`` and ``death`` rates, which default
to 1.0 and 0.5 respectively and are the first and second parameter. A
Yule tree can be simulated by setting the ``death`` rate to 0.0 (i.e.,
``ngesh.gen_tree(1.0, 0.0)``).

For reproducibility, the main parameter is the ``seed``, which takes any
Python object that can be represented as a string. The equivalence of
the generated tree can be verified by generating trees from the same
seed with different labels:

.. code:: python

   >>> tree3 = ngesh.gen_tree(1.0, 0.33, min_leaves=5, labels="human", seed="ngesh")
   >>> print(tree3)

         /-Fupsikmu
      /-|
     |   \-Enafu
   --|
     |   /-Zuhu
      \-|
        |   /-Poizes
         \-|
            \-Buge

   >>> tree4 = ngesh.gen_tree(1.0, 0.33, min_leaves=5, labels="bio", seed="ngesh")
   >>> print(tree4)

         /-Ennapus satvo
      /-|
     |   \-Sbuges asecis
   --|
     |   /-Spoizes rovues
      \-|
        |   /-Spupsicmus essis
         \-|
            \-Zusis spesus

The progression of the random generation can be verified by expanding
the number of leaves in the tree with human labels above:

.. code:: python

   >>> tree5 = ngesh.gen_tree(1.0, 0.33, min_leaves=6, labels="human", seed="ngesh")
   >>> print(tree5)

         /-Fupsikmu
      /-|
     |   \-Enafu
   --|
     |   /-Zuhu
     |  |
      \-|      /-Poizes
        |   /-|
         \-|   \-Buge
           |
            \-Hesi

Note that a single new taxon was created, “Buge”, while the rest of tree
is equivalent to one above.

All trees are normal ETE3 ``Tree`` objects, which means that all methods
from ETE3, including those for visualization and export, can be used
directly. For example, we can easily iterate over all the descendants to
grab the node distances, that is, the age in relation to the root. Note
that internal nodes are not named by default:

.. code:: python

   >>> for node in tree5.iter_descendants():
   >>>     print([node.name, node.dist])

   ["", 0.8181260108242158]
   ["", 1.3503435886693707]
   ["Fupsikmu", 1.1736487266924596]
   ["Enafu", 1.1736487266924596]
   ["Zuhu", 0.6414311488473046]
   ["", 0.5001786833989259]
   ["", 0.07984912684272014]
   ["Hesi", 0.1412524654483787]
   ["Poizes", 0.06140333860565855]
   ["Buge", 0.06140333860565855]

Random characters, matching the topology, can be added with the
``ngesh.add_characters()`` function. A new tree is returned (characters
are not added in-place) and it is necessary to specify, at least, the
number of characters to be simulated along with the ``k`` and ``theta``
arguments for the gamma distribution related to mutation events. Note
that the function also allows to simulate events equivalent to
horizontal gene transfer.

.. code:: python

   >>> tree_char = ngesh.add_characters(tree5, 10, 2.0, 1.0)
   >>> for node in tree_char.iter_descendants():
   >>>     if node.name:
   >>>         print("%10s - " % node.name, " ".join(["%02i" % v for v in node.chars]))

     Fupsikmu -  15 28 02 17 04 05 06 07 29 09 10 11 30 13 19
        Enafu -  15 16 02 17 04 24 06 25 26 09 27 11 12 13 19
         Zuhu -  00 01 02 31 20 05 06 21 32 09 22 23 33 13 14
         Hesi -  00 01 02 03 20 05 06 21 08 09 22 23 12 13 14
       Poizes -  00 01 02 03 20 05 06 21 08 09 22 23 12 13 14
         Buge -  00 01 02 03 20 05 06 21 08 09 22 23 12 13 14

The trees can be exported to different formats, as described in the
modules documentation and perfomed in the tests.
