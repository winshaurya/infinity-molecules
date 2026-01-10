#!/usr/bin/env python3
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'backend'))

from core_logic import generate_functionalized_isomers

print("Testing molecule generation...")

# Test 1: Simple case - 2 carbons with OH group
print("\n=== Test 1: 2 carbons + OH ===")
result1 = generate_functionalized_isomers(
    n_carbons=2,
    functional_groups=['OH'],
    n_double_bonds=0,
    n_triple_bonds=0,
    n_rings=0,
    carbon_types=['primary', 'secondary']
)
print(f"Generated {len(result1)} molecules: {result1}")

# Test 2: Just carbons, no functional groups
print("\n=== Test 2: 2 carbons, no functional groups ===")
result2 = generate_functionalized_isomers(
    n_carbons=2,
    functional_groups=[],
    n_double_bonds=0,
    n_triple_bonds=0,
    n_rings=0,
    carbon_types=['primary', 'secondary']
)
print(f"Generated {len(result2)} molecules: {result2}")

# Test 3: 0 carbons (should generate functional groups only)
print("\n=== Test 3: 0 carbons + COOH ===")
result3 = generate_functionalized_isomers(
    n_carbons=0,
    functional_groups=['COOH'],
    n_double_bonds=0,
    n_triple_bonds=0,
    n_rings=0,
    carbon_types=[]
)
print(f"Generated {len(result3)} molecules: {result3}")

# Test 4: Check what happens with empty parameters
print("\n=== Test 4: Empty parameters ===")
try:
    result4 = generate_functionalized_isomers(
        n_carbons=0,
        functional_groups=[],
        n_double_bonds=0,
        n_triple_bonds=0,
        n_rings=0,
        carbon_types=[]
    )
    print(f"Generated {len(result4)} molecules: {result4}")
except Exception as e:
    print(f"Error: {e}")

# Test what happens with the exact parameters that might be sent from frontend
print("\n=== Test 5: Frontend-like parameters ===")
# Simulate what happens when user selects 2 carbons, no functional groups, no unsaturations
frontend_params = {
    "carbon_count": 2,
    "double_bonds": 0,
    "triple_bonds": 0,
    "rings": 0,
    "carbon_types": ["primary", "secondary", "tertiary"],
    "functional_groups": []  # Empty list when no functional groups selected
}

result5 = generate_functionalized_isomers(
    n_carbons=frontend_params["carbon_count"],
    functional_groups=frontend_params["functional_groups"],
    n_double_bonds=frontend_params["double_bonds"],
    n_triple_bonds=frontend_params["triple_bonds"],
    n_rings=frontend_params["rings"],
    carbon_types=frontend_params["carbon_types"]
)
print(f"Frontend params generated {len(result5)} molecules: {result5}")

print("\n=== Summary ===")
print(f"Test 1 (2C+OH): {len(result1)} molecules")
print(f"Test 2 (2C only): {len(result2)} molecules")
print(f"Test 3 (0C+COOH): {len(result3)} molecules")
print(f"Test 5 (Frontend-like): {len(result5)} molecules")
