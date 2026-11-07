"""Standard-library readers for the NEXUS/PHYLIP matrices ngesh exports.

NEXUS import also accepts common interleaved single-symbol matrices, ambiguous
state sets and character subsets for calibration. It does not infer partitions.
"""

from __future__ import annotations

from pathlib import Path
import re

from .characters import Cell, Matrix


def _cells(
    text: str, tokens: bool, symbols: list[str], missing: str, gap: str, presence: bool
) -> list[Cell]:
    pieces = (
        re.findall(r"\{[^}]*\}|\([^)]*\)|[^\s]+", text)
        if tokens
        else re.findall(r"\{[^}]*\}|\([^)]*\)|[^\s]", text)
    )
    result = []
    for piece in pieces:
        if piece == missing:
            result.append(None)
        elif piece == gap:
            result.append(frozenset())
        else:
            raw = piece[1:-1] if piece.startswith(("{", "(")) else piece
            values = raw.split() if tokens else list(raw.replace(" ", ""))
            states = set()
            for value in values:
                if tokens:
                    state = int(value)
                else:
                    if value not in symbols:
                        raise ValueError(f"unknown NEXUS symbol: {value}")
                    state = symbols.index(value)
                if not presence or state != 0:
                    states.add(state)
            result.append(frozenset(states))
    return result


def read_nexus(
    path: str | Path, *, presence: bool = False, characters: list[int] | None = None
) -> Matrix:
    text = Path(path).read_text(encoding="utf-8-sig")
    text = re.sub(r"\[[\s\S]*?\]", "", text)
    block = re.search(r"\bmatrix\b([\s\S]*?);", text, re.I)
    if not block:
        raise ValueError("no NEXUS matrix")
    fmt = re.search(r"\bformat\b([\s\S]*?);", text, re.I)
    options = fmt.group(1) if fmt else ""
    if re.search(r"datatype\s*=\s*(dna|rna|protein|mixed)", options, re.I):
        raise ValueError("select a morphological STANDARD matrix, not molecular data")
    sym = re.search(r'symbols\s*=\s*"([^"]*)"', options, re.I)
    tokens = bool(re.search(r"\btokens\b", options, re.I))
    symbols = (
        sym.group(1).split()
        if sym and tokens
        else list(sym.group(1).replace(" ", ""))
        if sym
        else list("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
    )
    miss = re.search(r"missing\s*=\s*(\S)", options, re.I)
    gap = re.search(r"gap\s*=\s*(\S)", options, re.I)
    dimensions = re.search(r"\bnchar\s*=\s*(\d+)", text, re.I)
    result: Matrix = {}
    taxa_order = []
    continuation = 0
    for line in block.group(1).splitlines():
        line = line.strip()
        if not line:
            continuation = 0
            continue
        if dimensions and int(dimensions.group(1)) == 0:
            name = line[1:-1].replace("''", "'") if line.startswith("'") else line
            result[name] = []
            continue
        match = re.match(r"('(?:[^']|'')*'|\S+)\s+(.+)$", line)
        if match:
            name = match.group(1)
            name = name[1:-1].replace("''", "'") if name.startswith("'") else name
            data = match.group(2)
            if name not in result:
                taxa_order.append(name)
                result[name] = []
        elif taxa_order:
            name = taxa_order[continuation % len(taxa_order)]
            continuation += 1
            data = line
        else:
            raise ValueError("NEXUS row has no taxon label")
        result[name].extend(
            _cells(
                data,
                tokens,
                symbols,
                miss.group(1) if miss else "?",
                gap.group(1) if gap else "-",
                presence,
            )
        )
    dimensions = re.search(r"\bnchar\s*=\s*(\d+)", text, re.I)
    widths = {len(row) for row in result.values()}
    if len(widths) != 1 or (dimensions and widths != {int(dimensions.group(1))}):
        raise ValueError("NEXUS matrix dimensions do not match")
    if characters is not None:
        result = {t: [row[i] for i in characters] for t, row in result.items()}
    return result


def read_phylip(path: str | Path, *, names: dict[str, str] | None = None) -> Matrix:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    ntax, nchar = map(int, lines[0].split())
    result = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        name, _, data = line.partition(" ")
        sequence = data.replace(" ", "")
        if len(sequence) != nchar or any(c not in "01?" for c in sequence):
            raise ValueError("invalid binary PHYLIP row")
        result[(names or {}).get(name, name)] = [
            None if c == "?" else frozenset({1}) if c == "1" else frozenset() for c in sequence
        ]
    if len(result) != ntax:
        raise ValueError("PHYLIP taxon count mismatch")
    return result
