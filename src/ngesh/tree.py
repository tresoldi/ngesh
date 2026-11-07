"""Rooted trees, timed edges, and eNewick DAG serialization.

Distances are time units. A node's time is its parent's time plus its length.
Network nodes have one primary parent and optional secondary incoming edges.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import re
from typing import Any, Iterator


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float = 0.0
    time: float = 0.0
    attrs: dict[str, Any] = field(default_factory=dict)
    children: list[Node] = field(default_factory=list, repr=False)
    parent: Node | None = field(default=None, repr=False)

    def add(self, child: Node) -> Node:
        if child.parent is not None:
            raise ValueError("child already has a primary parent")
        child.parent = self
        self.children.append(child)
        return child

    def ancestors(self) -> Iterator[Node]:
        node: Node | None = self
        while node is not None:
            yield node
            node = node.parent


@dataclass
class Tree:
    root: Node
    extra_edges: list[tuple[Node, Node, float]] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.update_times()
        if self.extra_edges:
            nodes = list(self.traverse())
            members = set(nodes)
            incoming = {n: 0 for n in nodes}
            outgoing = {n: list(n.children) for n in nodes}
            for n in nodes:
                for child in n.children:
                    incoming[child] += 1
            for a, b, length in self.extra_edges:
                if a not in members or b not in members or not math.isfinite(length) or length < 0:
                    raise ValueError("invalid secondary network edge")
                if not math.isclose(a.time + length, b.time, abs_tol=1e-10):
                    raise ValueError("network incoming edges disagree on node time")
                outgoing[a].append(b)
                incoming[b] += 1
            ready = [n for n in nodes if incoming[n] == 0]
            seen = 0
            while ready:
                n = ready.pop()
                seen += 1
                for child in outgoing[n]:
                    incoming[child] -= 1
                    if incoming[child] == 0:
                        ready.append(child)
            if seen != len(nodes):
                raise ValueError("cyclic network")

    def traverse(self, order: str = "preorder") -> Iterator[Node]:
        if order not in {"preorder", "postorder"}:
            raise ValueError("order must be preorder or postorder")
        stack = [(self.root, False)]
        while stack:
            node, visited = stack.pop()
            if order == "preorder" or visited:
                yield node
            if not visited:
                if order == "postorder":
                    stack.append((node, True))
                stack.extend((child, False) for child in reversed(node.children))

    def tips(self, extant_only: bool = False) -> list[Node]:
        return [
            n
            for n in self.traverse()
            if not n.children and (not extant_only or not n.attrs.get("extinct", False))
        ]

    def update_times(self) -> None:
        for node in self.traverse():
            if not math.isfinite(node.length) or node.length < 0:
                raise ValueError("branch lengths must be finite and nonnegative")
            node.time = (node.parent.time if node.parent else 0.0) + node.length

    def alive(self, time: float) -> list[Node]:
        """Half-open edges [start,end); excludes zero-length sampled tips."""
        return [
            n
            for n in self.traverse()
            if (n.parent.time if n.parent else 0.0) <= time < n.time
            and not n.attrs.get("sampled_tip", False)
        ]

    def distance(self, a: Node, b: Node) -> float:
        ancestors = set(a.ancestors())
        mrca = next((n for n in b.ancestors() if n in ancestors), None)
        if mrca is None:
            raise ValueError("nodes do not belong to the same tree")
        return a.time + b.time - 2 * mrca.time

    def distance_at(self, a: Node, b: Node, time: float) -> float:
        ancestors = set(a.ancestors())
        mrca = next(n for n in b.ancestors() if n in ancestors)
        return 2 * (time - mrca.time)

    def prune(self, names: set[str]) -> Tree:
        """Copy the induced primary tree; preserve times and named ancestors."""
        if self.extra_edges:
            raise ValueError("prune the truth tree before constructing its network")
        copies: dict[Node, Node | None] = {}
        for node in self.traverse("postorder"):
            children = [copies[c] for c in node.children if copies[c] is not None]
            keep = node.name in names
            if not children and not keep:
                copies[node] = None
            elif len(children) == 1 and not keep:
                child = children[0]
                child.length += node.length
                copies[node] = child
            else:
                new = Node(node.name, node.length, attrs=dict(node.attrs))
                for child in children:
                    new.add(child)
                copies[node] = new
        root = copies[self.root]
        if root is None:
            raise ValueError("pruning removed every observation")
        result = Tree(root)
        for tip in result.tips():
            if tip.attrs.get("sampled_tip") and (
                tip.length != 0
                or tip.parent is None
                or not tip.parent.attrs.get("sampled_ancestor")
            ):
                tip.attrs["sampled_tip"] = False
        return result

    def rescale(self, target: float, mode: str = "height") -> None:
        if not math.isfinite(target) or target <= 0:
            raise ValueError("target must be finite and positive")
        tips = self.tips()
        if mode not in {"height", "mean", "crown"}:
            raise ValueError("mode must be height, mean or crown")
        if mode == "crown":
            common = set(tips[0].ancestors())
            for tip in tips[1:]:
                common.intersection_update(tip.ancestors())
            crown = next(n for n in tips[0].ancestors() if n in common)
            base = max(n.time for n in tips) - crown.time
        else:
            base = (
                max(n.time for n in tips)
                if mode == "height"
                else sum(n.time for n in tips) / len(tips)
            )
        if base == 0:
            raise ValueError("cannot rescale a zero-height tree")
        ratio = target / base
        for node in self.traverse():
            node.length *= ratio
        self.extra_edges = [(a, b, d * ratio) for a, b, d in self.extra_edges]
        self.update_times()

    def newick(self) -> str:
        hybrid = {
            node: f"#H{i + 1}"
            for i, node in enumerate(dict.fromkeys(b for _, b, _ in self.extra_edges))
        }
        extra: dict[Node, list[tuple[Node, float]]] = {}
        for a, b, length in self.extra_edges:
            extra.setdefault(a, []).append((b, length))

        def label(name: str) -> str:
            if not name:
                return ""
            return name if re.fullmatch(r"[\w.\-]+", name) else "'" + name.replace("'", "''") + "'"

        rendered: dict[Node, str] = {}
        for node in self.traverse("postorder"):
            children = [rendered[c] for c in node.children]
            children.extend(hybrid[b] + ":" + format(d, ".17g") for b, d in extra.get(node, []))
            body = "(" + ",".join(children) + ")" if children else ""
            rendered[node] = (
                body + label(node.name) + hybrid.get(node, "") + ":" + format(node.length, ".17g")
            )
        return rendered[self.root] + ";"

    @classmethod
    def from_newick(cls, text: str) -> Tree:
        """Parse quoted labels, lengths, comments, and #H hybrid references.

        Inheritance probabilities and support annotations are not interpreted.
        """
        tokens = re.findall(r"'(?:[^']|'')*'|\[[^\]]*\]|[(),:;]|[^\s(),:;\[\]]+", text)
        tokens = [t for t in tokens if not t.startswith("[")]
        pos = 0
        hybrids: dict[str, Node] = {}
        occurrences: list[tuple[Node, str]] = []

        def parse() -> Node:
            nonlocal pos
            node = Node()
            if pos < len(tokens) and tokens[pos] == "(":
                pos += 1
                node.add(parse())
                while pos < len(tokens) and tokens[pos] == ",":
                    pos += 1
                    node.add(parse())
                if pos >= len(tokens) or tokens[pos] != ")":
                    raise ValueError("unclosed Newick clade")
                pos += 1
            if pos < len(tokens) and tokens[pos] not in "(),:;":
                name = tokens[pos]
                pos += 1
                if name.startswith("'"):
                    node.name = name[1:-1].replace("''", "'")
                    if pos < len(tokens) and tokens[pos].startswith("#H"):
                        name = tokens[pos]
                        pos += 1
                    else:
                        name = ""
                match = re.search(r"#H[\w]+$", name)
                if match:
                    tag = match.group()
                    if not node.name:
                        node.name = name[: match.start()]
                    occurrences.append((node, tag))
                    if node.children or node.name:
                        if tag in hybrids:
                            raise ValueError("multiple hybrid definitions")
                        hybrids[tag] = node
                elif name:
                    node.name = name
            if pos < len(tokens) and tokens[pos] == ":":
                pos += 1
                try:
                    node.length = float(tokens[pos])
                except (ValueError, IndexError) as error:
                    raise ValueError("invalid Newick length") from error
                pos += 1
            return node

        try:
            root = parse()
        except RecursionError as error:
            raise ValueError("Newick nesting exceeds parser limit") from error
        if pos != len(tokens) - 1 or tokens[pos] != ";":
            raise ValueError("Newick must end with one semicolon")
        for node, tag in occurrences:
            hybrids.setdefault(tag, node)
        edges = []
        for node, tag in occurrences:
            target = hybrids[tag]
            if node is target:
                continue
            if node.children or node.parent is None:
                raise ValueError("invalid hybrid reference")
            if target in set(node.parent.ancestors()):
                raise ValueError("cyclic network")
            edges.append((node.parent, target, node.length))
            node.parent.children.remove(node)
        return cls(root, edges)
