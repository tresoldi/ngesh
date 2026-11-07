"""Public orchestration and stable component/replicate seeds."""
from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

from .bundle import write_bundle
from .characters import CharacterConfig, Evolution, evolve
from .models import TreeConfig, generate_tree
from .rng import RandomStream, Seed, derive_seed
from .horizontal import HorizontalConfig
from .observation import ObservationConfig, observe


PRESETS: dict[str, dict[str, Any]] = {
    'lexical': {'tree':{}, 'characters':{'model':'lexical','rate':.1,'alpha':.7}},
    'genecontent': {'tree':{},'characters':{'model':'dollo','rate':1.0,'loss_rate':.1,'n_characters':500,'root_frequency':.5},'horizontal':{'mode':'add'}},
    'morphology': {'tree':{},'characters':{'model':'mk','rate':.2,'k':4,'n_characters':100,'alpha':.7}},
    'stemma': {'tree':{'birth':1.0,'death':.8,'ancestor_rate':.2,'copy_sigma':.6},
               'characters':{'model':'textual','rate':.2,'p_parallel':.2,'omission_rate':.05}},
}


def simulate(*, seed: Seed = 0, preset: str = 'lexical', config: dict[str, Any] | None = None,
             out: str | Path | None = None) -> Evolution:
    if preset not in PRESETS:
        raise ValueError(f'unknown preset: {preset}')
    config = config or {}
    unknown = set(config) - {'tree','characters','output','horizontal','observation'}
    if unknown:
        raise ValueError(f'unknown configuration sections: {sorted(unknown)}')
    base = PRESETS[preset]
    tree_cfg = TreeConfig(**(base.get('tree',{}) | config.get('tree',{})))
    char_cfg = CharacterConfig(**(base.get('characters',{}) | config.get('characters',{})))
    seeds = {component:derive_seed(seed,component) for component in ['tree','characters','observation','horizontal','geography','tree_sampling']}
    tree = generate_tree(replace(tree_cfg,full=True,sampling=1),seed=seeds['tree'])
    hcfg = HorizontalConfig(**(base.get('horizontal',{}) | config.get('horizontal',{})))
    ocfg = ObservationConfig(**(base.get('observation',{}) | config.get('observation',{})))
    if char_cfg.model in {'binary','dollo'} and (ocfg.lumping or ocfg.splitting):
        raise ValueError('class lumping/splitting is not defined for binary presence/absence')
    result = evolve(tree,char_cfg,seed=seeds['characters'],horizontal=hcfg,
                    horizontal_rng=RandomStream(seeds['horizontal']),geography_rng=RandomStream(seeds['geography']))
    srng = RandomStream(seeds['tree_sampling'])
    keep = {n.name for n in tree.tips() if n.attrs.get('sampled_tip') or
            (n.attrs.get('extinct',False) and tree_cfg.full) or
            (not n.attrs.get('extinct',False) and (tree_cfg.sampling==1 or srng.random()<tree_cfg.sampling))}
    if not keep:
        raise ValueError('tree sampling retained no observations')
    if keep!=set(result.truth):
        result.history_tree = result.tree
        result.tree = result.tree.prune(keep)
        result.truth = {t:row for t,row in result.truth.items() if t in keep}
    observed = observe(result,ocfg,rng=RandomStream(seeds['observation']))
    if out is not None:
        master = {'type':type(seed).__name__,'value':seed.hex() if isinstance(seed,bytes) else seed}
        output = config.get('output',{})
        if set(output) - {'include_ci','max_quartets'}:
            raise ValueError('unknown output configuration')
        write_bundle(out,result,params={'preset':preset,'master_seed':master,
                     'derived_seeds':seeds,'tree_config':asdict(tree_cfg),
                     'horizontal_config':asdict(hcfg),'observation_config':asdict(ocfg),'config':config},
                     observed=observed.matrix, truth_tree=observed.tree, observation_log=observed.log,
                     retained_groups=observed.retained_groups,
                     include_ci=output.get('include_ci',True),max_quartets=output.get('max_quartets',1000))
    return result
