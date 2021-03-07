import hashlib
from ete3 import Tree
import ngesh

# Some pre-generated newick trees for testing
_TREES = [
    "(((Nucroto zolos:0.339415,Coddopus zoggaus:0.339415)1:2.29301,Aporos "
    "oiasis:2.63243)1:0.511706,Spetitis mubvoppis:3.14413);",
    "((((Ataba eolus:0.274414,Dasoros audus:0.274414)1:1.59309,(Uvuros "
    "spalus:1.63679,(Zilavis sicagas:0.158265,Uzazopus aolo:0.158265)"
    "1:1.47853)1:0.230708)1:0.976915,Spempo gipus:2.84442)1:1.08647,"
    "(((Cobbas linis:0.242355,Ciggus sopebbas:0.242355)1:1.33741,"
    "(Vaoras ovamla:0.235349,Nirceo spemgazzo:0.235349)1:1.34442)"
    "1:0.796904,(Wiopepus spiparzas:1.86067,Eavoros airos:0.906654)"
    "1:0.515998)1:1.55422);",
    "((((Egasis reggasas:0.0747242,(Niponis ecales:0.381168,Vossignis "
    "spepupapes:0.381168)1:0.0385491)1:0.697337,Tuapebbos "
    "eppozas:0.457931)1:2.02569,(Vucotes sparas:3.09293,((Ruttinges "
    "daerus:0.300024,(Ovugas setlanes:0.297199,Rezas "
    "emdis:0.260688)1:1.20787)1:1.35227,(Gemdas naoris:1.77171,(Tirreres "
    "sbissas:0.161549,Speppizipus spebbes:0.33078)1:1.44093)1:1.08563)"
    "1:0.235588)1:0.0498126)1:0.593566,Seti spines:0.874473);",
    "(((Uorus sbepis:1.00768,Winozbus veses:1.00768)1:2.02024,(((((Ragis "
    "zubarivos:0.239467,(Zoprotollas naessadda:0.032829,Sauvo "
    "tibbunnis:0.032829)1:0.206638)1:0.298691,Mirba "
    "ummes:0.538158)1:0.751761,Gitopes guaras:1.28992)1:0.219143,(Uvoros "
    "sbego:0.429929,Rouslis novuddenus:0.429929)1:1.07913)1:0.0145254,"
    "(((Spaupis gieves:0.155583,Ruasus coeddas:0.155583)1:0.222901,Tipapo "
    "vati:0.378484)1:0.044333,Zuvcus luruzuco:0.422817)1:1.10077)"
    "1:1.50433)1:0.97142,((Nunpis ozabmis:0.491355,Sbigirrecas sbilpo:"
    "0.491355)1:0.33275,Agolus galgis:0.824105)1:3.17523);",
]


def test_add_characters():
    # gamma parameters
    NUM_CONCEPTS = 10
    k = 4.0  # shape
    th = 1.0  # scale
    z = 1.045  # "zipf" correction
    z = 1.0  # TODO: fix test later

    # Add characters to all trees, for coverage
    trees = [
        ngesh.add_characters(
            Tree(newick), NUM_CONCEPTS, k=k, th=th, mut_exp=z, seed="myseed"
        )
        for newick in _TREES
    ]

    # Assert the first one
    digest = hashlib.sha256(str(ngesh.tree2wordlist(trees[0])).encode("utf-8")).digest()
    assert (
        digest
        == b"\xdc<\x1f\x10N\xbf\xcc\xc4|l26\x10\xfc\xbaN\xb7\\c\x8bB\xca\x95.\xcbH\x82T\xa3\xbd\xff\x15"
    )


def test_add_characters_with_hgt():
    tree_hgt = Tree(_TREES[0])
    ngesh.add_characters(
        tree_hgt,
        10,
        k=4.0,
        th=1.0,
        mut_exp=1.0,  # z=1.045, TODO: fix test later
        k_hgt=2.0,
        th_hgt=1.1,
        seed="myseed",
    )
    digest = hashlib.sha256(str(ngesh.tree2wordlist(tree_hgt)).encode("utf-8")).digest()

    assert (
        digest
        == b"\xc0\x84\x8a\xa8\x0b\xd6\xf8\x1b\x9c\xc7\xe6\xaf]\xe7\xee\x105\x989\x9b\xcfke\xe5\xf3\x03\xdc\x17U\x14 \xac"
    )
