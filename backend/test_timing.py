#!/usr/bin/env python3
"""Simple timing test for molecule generation."""

import time
import json
from ExtrafastInfi import generate_functionalized_isomers, canonicalize_smiles

def test_generation():
    print("Testing Python molecule generation...")

    # Test parameters for ~100k molecules
    carbon_count = 12
    functional_groups = ["OH", "COOH", "NH2"]
    double_bonds = 2
    triple_bonds = 0
    rings = 1
    carbon_types = ["primary", "secondary", "tertiary"]

    print(f"Parameters: {carbon_count}C, {double_bonds}DB, {triple_bonds}TB, {rings}R, {functional_groups}, {carbon_types}")

    start_time = time.time()
    smiles_list = generate_functionalized_isomers(
        carbon_count, functional_groups, double_bonds, triple_bonds, rings, carbon_types
    )
    end_time = time.time()

    generation_time = end_time - start_time
    molecule_count = len(smiles_list)

    print(".2f")
    print(f"Molecules generated: {molecule_count}")
    print(".1f")

    # Show first 5 SMILES
    print("\nFirst 5 SMILES:")
    for i, smiles in enumerate(smiles_list[:5]):
        print(f"  {i+1}. {smiles}")

    return generation_time, molecule_count, smiles_list

if __name__ == "__main__":
    test_generation()
