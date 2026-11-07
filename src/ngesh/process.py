"""Chronological CTMC scheduler shared by independent and horizontal evolution."""
from __future__ import annotations

import bisect
import heapq
import math
from typing import Any

from .characters import CharacterConfig, Event, Evolution, character_rates
from .horizontal import HorizontalConfig
from .rng import RandomStream
from .tree import Node, Tree


def run_process(tree: Tree, cfg: CharacterConfig, rng: RandomStream,
                horizontal: HorizontalConfig, hrng: RandomStream,
                grng: RandomStream) -> Evolution:
    cfg.validate()
    horizontal.validate()
    rates = character_rates(cfg,rng)
    count = cfg.n_characters
    # Names are stable lineage IDs; sample observations retain tip names.
    tip_names = [n.name for n in tree.tips()]
    if len(set(tip_names)) != len(tip_names):
        raise ValueError('tip names must be unique')
    reserved = set(tip_names)
    for i,node in enumerate(tree.traverse()):
        if not node.name or (node.children and node.name in reserved):
            node.name = f'I{i:06d}'
        if node.children:
            reserved.add(node.name)
    if len({n.name for n in tree.traverse()}) != len(list(tree.traverse())):
        raise ValueError('lineage names must be unique')
    nodes = list(tree.traverse())
    branch_rates = {n.name:rng.lognormal(cfg.clock_sigma) for n in nodes}
    counters = [cfg.pool_size] * count
    states: dict[Node,list[frozenset[int]]] = {}
    deadlines: dict[Node,list[dict[int,float]]] = {}
    active: set[Node] = set()
    positions: dict[Node,tuple[float,float,float]] = {}
    events: list[Event] = []
    heap: list[tuple[float,int,int,str,Node,int,Any]] = []
    serial = 0

    def queue(t: float, priority: int, kind: str, node: Node, char: int = -1, value: Any = None) -> None:
        nonlocal serial
        serial += 1
        heapq.heappush(heap,(t,priority,serial,kind,node,char,value))

    def schedule_change(t: float, node: Node, c: int) -> None:
        rate = rates[c] * branch_rates[node.name]
        if rate:
            when = t + rng.exponential(rate)
            if when < node.time:
                queue(when,4,'change',node,c)

    def schedule_borrow(t: float, node: Node) -> None:
        rate = horizontal.borrowing_rate * count
        if rate:
            when = t + hrng.exponential(rate)
            if when < node.time:
                queue(when,5,'borrow',node)

    def position(node: Node, t: float) -> tuple[float,float]:
        x,y,last = positions[node]
        if t > last:
            sd = horizontal.diffusion * math.sqrt(t-last)
            if sd:
                x += sd * grng.normal()
                y += sd * grng.normal()
            positions[node] = x,y,t
        return x,y

    def choose_donor(node: Node, t: float, copying: bool = False) -> Node | None:
        candidates = sorted((n for n in active if n is not node and
                        n.parent is not None and n.parent.time <= t < n.time and
                        (not copying or n.parent is not node.parent)),key=lambda n:n.name)
        if not candidates:
            return None
        distances = []
        if horizontal.distance == 'geography':
            x,y = position(node,t)
            distances = [math.dist((x,y),position(d,t)) for d in candidates]
        else:
            distances = [tree.distance_at(node,d,t) for d in candidates]
        # Subtract min distance before exp, avoiding all-zero underflow.
        minimum = min(distances)
        weights = [math.exp(-(d-minimum)/horizontal.scale) for d in distances]
        cumulative = []
        total = 0.0
        for weight in weights:
            total += weight
            cumulative.append(total)
        return candidates[min(bisect.bisect_right(cumulative,hrng.random()*total),len(candidates)-1)]

    for node in nodes:
        start = node.parent.time if node.parent else 0.0
        # Starts at a split must follow parent completion, including zero edges.
        # Preorder insertion resolves chains of zero-length boundaries.
        queue(start,1,'start',node)
        if node.time > start:
            queue(node.time,0,'end',node)
    while heap:
        t,_,_,kind,node,c,value = heapq.heappop(heap)
        if kind == 'start':
            if node.parent is None:
                states[node] = [frozenset({0}) for _ in rates]
                deadlines[node] = [{} for _ in rates]
                positions[node] = (0,0,t)
            else:
                states[node] = states[node.parent][:]
                deadlines[node] = [dict(d) for d in deadlines[node.parent]]
                x,y = position(node.parent,t)
                positions[node] = (x,y,t)
            if t < node.time and not node.attrs.get('sampled_tip'):
                active.add(node)
                for i in range(count):
                    schedule_change(t,node,i)
                    for state,when in deadlines[node][i].items():
                        if when < node.time:
                            queue(when,3,'expire',node,i,state)
                schedule_borrow(t,node)
                if horizontal.contamination_fraction and node.parent and len(node.parent.children)>1:
                    queue(t,2,'contaminate',node)
            continue
        if kind == 'end':
            active.discard(node)
            if horizontal.distance == 'geography':
                position(node,t)
            continue
        if kind in {'borrow','contaminate'}:
            donor = choose_donor(node,t,copying=kind=='contaminate')
            if donor is not None:
                if kind == 'borrow':
                    chars = [hrng.index(count)]
                elif horizontal.contamination_blocks:
                    width = round(horizontal.contamination_fraction*count)
                    start = hrng.index(count-width+1)
                    chars = list(range(start,start+width))
                else:
                    chars = [i for i in range(count) if hrng.random()<horizontal.contamination_fraction]
                for i in chars:
                    old = states[node][i]
                    new = states[donor][i]
                    if horizontal.mode == 'add' and kind == 'borrow':
                        new = old | new
                    states[node][i] = new
                    deadlines[node][i] = {}  # incoming reading becomes established
                    events.append(Event(t,node.name,i,tuple(sorted(old)),tuple(sorted(new)),kind,donor.name))
            if kind == 'borrow':
                schedule_borrow(t,node)
            continue
        old = states[node][c]
        if kind == 'expire':
            if deadlines[node][c].get(value) != t:
                continue
            new = old - {value}
            deadlines[node][c].pop(value,None)
            event_kind = 'synonym_loss'
        else:
            if cfg.p_parallel and rng.random()<cfg.p_parallel:
                new_id = 1+rng.index(cfg.pool_size)
            else:
                counters[c] += 1
                new_id = counters[c]
            if cfg.synonym_duration:
                for s in old-{new_id}:
                    if s not in deadlines[node][c]:
                        when = t+rng.exponential(1/cfg.synonym_duration)
                        deadlines[node][c][s] = when
                        if when < node.time:
                            queue(when,3,'expire',node,c,s)
                deadlines[node][c].pop(new_id,None)
                new = old | {new_id}
            else:
                new = frozenset({new_id})
                deadlines[node][c].clear()
            event_kind = 'replacement'
            schedule_change(t,node,c)
        states[node][c] = frozenset(new)
        events.append(Event(t,node.name,c,tuple(sorted(old)),tuple(sorted(new)),event_kind))
    all_states = {n.name:states[n] for n in nodes}
    truth = {n.name:states[n][:] for n in tree.tips()}
    return Evolution(tree,cfg,[f'C{i+1:04d}' for i in range(count)],truth,events,rates,branch_rates,all_states)
