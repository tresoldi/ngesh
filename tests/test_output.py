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


def test_tree_output():
    # Add characters to all test trees
    trees = [
        ngesh.add_characters(
            Tree(newick), 100, k=4.0, th=1.0, mut_exp=1.05, seed="myseed"
        )
        for newick in _TREES
    ]

    # Assert the first one
    digest_nx = hashlib.sha256(str(ngesh.tree2nexus(trees[0])).encode("utf-8")).digest()
    digest_wl = hashlib.sha256(
        str(ngesh.tree2wordlist(trees[0])).encode("utf-8")
    ).digest()

    assert (
        digest_nx
        == b"E\xf8\x97\xb6*\x7f\xf4_j\x89\x02dn\x1d\xbe\xb0\xb6\xcd\xd9.\xca:\x9ft\xe2m\xc5y\xa5\xaa\x0fa"
    )
    assert (
        digest_wl
        == b"\xb6\xc8i\xf4!\xf4l\x91\xb9\x8d\xb5Kae\x1aF\x9c\xfd\n \x06\xf2D\x1e<\xdd(U]6(\xbf"
    )

    # Test output for a tree without characters
    tree_nochar = Tree(_TREES[0])
    digest_nx_nochar = hashlib.sha256(
        str(ngesh.tree2nexus(tree_nochar)).encode("utf-8")
    ).digest()
    assert (
        digest_nx_nochar
        == b"A\x06N\x97n1iD\x01n\x07T\x99al\xa8)m\n\x93\x9ajKp\x07\x9bc\x8a\x03\x1a\x8d\xb7"
    )
