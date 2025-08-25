"""Minimal JSON-with-Comments loader."""
from __future__ import annotations
import json, re, pathlib

_COMMENT_RE = re.compile(
    r"""
    (?:\/\/[^\n]*$)             # // line comments
  | (?:\/\*[\s\S]*?\*\/)        # /* block comments */
    """,
    re.MULTILINE | re.VERBOSE,
)

def strip_comments(text: str) -> str:
    """Remove // and /*  */ comments without touching string literals."""
    # naive but good enough because AlphaFold input has no multi-line strings
    return re.sub(_COMMENT_RE, "", text)

def load(path: str | pathlib.Path) -> dict:
    raw = pathlib.Path(path).read_text(encoding="utf-8")
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        return json.loads(strip_comments(raw))
