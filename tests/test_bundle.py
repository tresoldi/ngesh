import json
from pathlib import Path
import tempfile
import unittest

from ngesh import simulate
from ngesh.__main__ import run_grid
from ngesh.bundle import read_long, read_matrix, verify_bundle
from ngesh.config import toml_subset


class BundleTests(unittest.TestCase):
    def test_roundtrip_and_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp)/'bundle'
            result = simulate(seed='roundtrip',config={'tree':{'n_taxa':8},'characters':{'n_characters':12,'synonym_duration':2}},out=dest)
            self.assertEqual(read_long(dest/'data.tsv',sorted(result.truth),result.groups),result.truth)
            self.assertEqual(read_matrix(dest/'truth.json'),result.truth)
            self.assertTrue(verify_bundle(dest))
            for line in (dest/'data.tsv').read_text().splitlines()[1:]:
                taxon, group, state = line.split('\t')
                self.assertTrue(state.startswith(group+'.'))
                self.assertTrue(state.rsplit('.',1)[1].isdigit())
                self.assertLess(len(line.encode())+1,4096)
            (dest/'data.tsv').write_text('corrupt')
            self.assertFalse(verify_bundle(dest))

    def test_grid_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = {'seed':'grid','replicates':2,'base':{'tree':{'n_taxa':4},'characters':{'n_characters':3}},'grid':{'characters__rate':[.1,.2]}}
            out = Path(tmp)
            self.assertEqual(run_grid(cfg,out),{'completed':4,'skipped':0})
            self.assertEqual(run_grid(cfg,out),{'completed':0,'skipped':4})
            seeds = [json.loads(p.read_text())['master_seed']['value'] for p in out.glob('*/params.json')]
            self.assertEqual(len(set(seeds)),4)
            cfg['seed']='other'
            with self.assertRaises(ValueError):
                run_grid(cfg,out)

    def test_subset_toml(self):
        text = 'seed="a#b" # comment\n[tree]\nn_taxa=4\n[grid]\ncharacters__rate=[0.1,0.2]\n'
        import sys
        parsed = toml_subset(text)
        if sys.version_info >= (3,11):
            import tomllib
            self.assertEqual(parsed,tomllib.loads(text))
        self.assertEqual(parsed['seed'],'a#b')
