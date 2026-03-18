#!/usr/bin/env python3
"""CLI bridge exposing ExtrafastInfi enumeration for the JS test suite."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    from ExtrafastInfi import (  # type: ignore
        canonicalize_smiles,
        generate_functionalized_isomers,
        validate_structure_possibility,
    )
except Exception as exc:  # pragma: no cover - surfaced to Node caller
    print(json.dumps({"error": f"Failed to import ExtrafastInfi: {exc}"}))
    sys.exit(0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="ExtrFastInfi SMILES export helper")
    parser.add_argument("--carbon_count", type=int, required=True)
    parser.add_argument("--double_bonds", type=int, default=0)
    parser.add_argument("--triple_bonds", type=int, default=0)
    parser.add_argument("--rings", type=int, default=0)
    parser.add_argument("--functional_groups", type=str, default="[]")
    parser.add_argument("--carbon_types", type=str, default='["primary","secondary","tertiary"]')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        functional_groups = json.loads(args.functional_groups) or []
        carbon_types = json.loads(args.carbon_types) or ["primary", "secondary", "tertiary"]
    except json.JSONDecodeError as exc:
        print(json.dumps({"error": f"Invalid JSON for arguments: {exc}"}))
        return

    ok, reason = validate_structure_possibility(
        args.carbon_count,
        functional_groups,
        args.double_bonds,
        args.triple_bonds,
        carbon_types,
        args.rings,
    )
    if not ok:
        print(json.dumps({"error": reason}))
        return

    try:
        smiles = generate_functionalized_isomers(
            args.carbon_count,
            functional_groups,
            args.double_bonds,
            args.triple_bonds,
            args.rings,
            carbon_types,
        )
    except Exception as exc:  # pragma: no cover - surfaced to JS caller
        print(json.dumps({"error": f"Generation failed: {exc}"}))
        return

    canonical = sorted(
        {canonicalize_smiles(entry) for entry in smiles if entry}
    )
    print(json.dumps({"smiles": canonical}))


if __name__ == "__main__":
    main()
