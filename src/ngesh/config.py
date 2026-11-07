"""TOML configuration. Python 3.11+ uses tomllib.

3.10 supports the documented scalar/array/table subset with a strict local
parser; unsupported TOML constructs fail rather than being misinterpreted.
"""

from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Any


def toml_subset(text: str) -> dict[str, Any]:
    """3.10 subset: bare-key tables, scalar literals and single-line arrays.

    No inline tables, dates, multiline strings, quoted/dotted keys, or arrays
    of tables. This avoids a runtime dependency solely for configuration.
    """
    root: dict[str, Any] = {}
    current = root
    tables: set[str] = set()

    def parse(value: str) -> Any:
        value = value.strip()
        if value in {"true", "false"}:
            return value == "true"
        if value.startswith("[") and value.endswith("]"):
            # Tokenize array commas outside quoted strings/nested arrays.
            parts, start, depth, quote, escape = [], 1, 0, "", False
            for i, ch in enumerate(value[1:-1], 1):
                if quote:
                    if ch == quote and not escape:
                        quote = ""
                    escape = quote == '"' and ch == "\\" and not escape
                elif ch in "\"'":
                    quote = ch
                elif ch == "[":
                    depth += 1
                elif ch == "]":
                    depth -= 1
                elif ch == "," and depth == 0:
                    parts.append(value[start:i])
                    start = i + 1
            parts.append(value[start:-1])
            if quote or depth != 0:
                raise ValueError("unclosed TOML array value")
            if len(parts) == 1 and not parts[0].strip():
                return []
            if not parts[-1].strip():
                parts.pop()
            if any(not p.strip() for p in parts):
                raise ValueError("empty TOML array element")
            return [parse(p) for p in parts]
        if value.startswith(('"', "'")):
            if value.startswith("'"):
                if not value.endswith("'") or "'" in value[1:-1]:
                    raise ValueError("invalid literal TOML string")
                return value[1:-1]
            return json.loads(value)
        if re.fullmatch(
            r"[+-]?(?:0|[1-9](?:_?\d)*)(?:\.\d(?:_?\d)*)?(?:[eE][+-]?\d(?:_?\d)*)?", value
        ):
            return float(value) if any(x in value for x in ".eE") else int(value)
        raise ValueError(f"unsupported TOML value on Python 3.10: {value}")

    for raw in text.splitlines():
        # Strip comments only outside quoted strings.
        quote, escape, end = "", False, len(raw)
        for i, ch in enumerate(raw):
            if quote:
                if ch == quote and not escape:
                    quote = ""
                escape = quote == '"' and ch == "\\" and not escape
            elif ch in "\"'":
                quote = ch
            elif ch == "#":
                end = i
                break
        line = raw[:end].strip()
        if not line:
            continue
        if re.fullmatch(r"\[[A-Za-z0-9_.-]+\]", line):
            if line in tables:
                raise ValueError("duplicate TOML table")
            tables.add(line)
            current = root
            for key in line[1:-1].split("."):
                if not re.fullmatch(r"[A-Za-z0-9_-]+", key) or not isinstance(current, dict):
                    raise ValueError("invalid TOML table")
                if key in current and not isinstance(current[key], dict):
                    raise ValueError("TOML table conflicts with a scalar")
                current = current.setdefault(key, {})
            continue
        key, sep, value = line.partition("=")
        key = key.strip()
        if not sep or not re.fullmatch(r"[A-Za-z0-9_-]+", key) or key in current:
            raise ValueError("invalid or duplicate TOML key")
        current[key] = parse(value)
    return root


def load_config(path: str | Path) -> dict[str, Any]:
    try:
        import tomllib
    except ImportError:
        return toml_subset(Path(path).read_text(encoding="utf-8"))
    with Path(path).open("rb") as handle:
        return tomllib.load(handle)
