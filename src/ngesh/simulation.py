"""Public orchestration and stable component/replicate seeds."""
from __future__ import annotations

from dataclasses import asdict
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
    seeds = {component:derive_seed(seed,component) for component in ['tree','characters','observation','horizontal','geography']}
    tree = generate_tree(tree_cfg,seed=seeds['tree'])
    hcfg = HorizontalConfig(**(base.get('horizontal',{}) | config.get('horizontal',{})))
    ocfg = ObservationConfig(**(base.get('observation',{}) | config.get('observation',{})))
    result = evolve(tree,char_cfg,seed=seeds['characters'],horizontal=hcfg,
                    horizontal_rng=RandomStream(seeds['horizontal']),geography_rng=RandomStream(seeds['geography']))
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
