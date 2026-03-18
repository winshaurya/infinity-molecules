import re
import csv
import threading
import queue
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, SDWriter, rdMolDescriptors, Descriptors
from rdkit.Chem.rdchem import BondType
import multiprocessing
from multiprocessing import Pool, cpu_count
from functools import lru_cache, partial
import time
import psutil
import random
from typing import List, Set, Tuple, Dict
import itertools
from concurrent.futures import ProcessPoolExecutor, as_completed
import sys
import os
import shutil
from collections import OrderedDict, defaultdict
import networkx as nx
import tkinter.font as tkfont
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

DUPLICATE_REGISTRY = {}
FILE_EXISTENCE_CACHE = {}
_canonical_cache = {}

class UniqueNomenclatureSystem:
    def __init__(self):
        self.used_names = set()
        self.structure_cache = {}

    def get_branching_pattern(self, mol):
        if mol.GetNumAtoms() == 0:
            return "NULL"
        G = nx.Graph()
        for atom in mol.GetAtoms():
            G.add_node(atom.GetIdx(), symbol=atom.GetSymbol())
        for bond in mol.GetBonds():
            G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), 
                      bond_type=bond.GetBondType())
        return nx.weisfeiler_lehman_graph_hash(G)

    def get_functional_group_positions(self, mol, functional_groups):
        positions = []
        smarts_patterns = {
            'OH': '[OX2H]',
            'COOH': 'C(=O)[OH]',
            'CHO': '[CX3H1](=O)[#6]',
            'Ketone': '[CX3](=O)[#6]',
            'CN': 'C#N',
            'NC': '[N+]#[C-]',
            'F': '[F]',
            'Cl': '[Cl]',
            'Br': '[Br]',
            'I': '[I]',
            'NO2': '[N+](=O)[O-]',
            'NH2': '[NX3;H2]',
            'Ether': '[OD2]([#6])[#6]',
            'Imine': '[CX3]([#6])[#6]=[NX2]',
            'Amide': 'C(=O)[NX3]',
            'Azide': '[N-]=[N+]=N',
            'S_Bivalent': '[SX2H]',
            'S_Tetravalent': '[SX3](=O)',
            'S_Hexavalent': '[SX4](=O)(=O)',
            'S_Chain_Bi': '[SX2]([#6])[#6]',
            'S_Chain_Tetra': '[SX3](=O)([#6])[#6]',
            'S_Chain_Hexa': '[SX4](=O)(=O)([#6])[#6]'
        }
        for fg_name, smarts in smarts_patterns.items():
            if fg_name in functional_groups:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern:
                    matches = mol.GetSubstructMatches(pattern)
                    positions.append(f"{fg_name}_{len(matches)}")
        return sorted(positions)

    def get_unsaturation_pattern(self, mol):
        double_bonds = 0
        triple_bonds = 0
        aromatic_bonds = 0
        for bond in mol.GetBonds():
            bond_type = bond.GetBondType()
            if bond_type == BondType.DOUBLE:
                double_bonds += 1
            elif bond_type == BondType.TRIPLE:
                triple_bonds += 1
            elif bond_type == BondType.AROMATIC:
                aromatic_bonds += 1
        pattern = []
        if double_bonds > 0:
            pattern.append(f"DB{double_bonds}")
        if triple_bonds > 0:
            pattern.append(f"TB{triple_bonds}")
        if aromatic_bonds > 0:
            pattern.append(f"AR{aromatic_bonds}")
        return "_".join(pattern) if pattern else "SAT"

    def get_ring_pattern(self, mol):
        if not mol.GetRingInfo().NumRings():
            return "ACYCLIC"
        ring_info = mol.GetRingInfo()
        ring_sizes = []
        for ring in ring_info.AtomRings():
            ring_sizes.append(len(ring))
        size_count = defaultdict(int)
        for size in ring_sizes:
            size_count[size] += 1
        pattern_parts = []
        for size in sorted(size_count.keys()):
            count = size_count[size]
            pattern_parts.append(f"R{size}_{count}")
        return "_".join(pattern_parts)

    def generate_unique_name(self, mol, formula, functional_groups, carbon_count, 
                           n_double_bonds, n_triple_bonds, n_rings, index):
        base_name = f"{formula}_{index:04d}"
        structural_parts = []
        unsat_pattern = self.get_unsaturation_pattern(mol)
        structural_parts.append(unsat_pattern)
        ring_pattern = self.get_ring_pattern(mol)
        structural_parts.append(ring_pattern)
        fg_positions = self.get_functional_group_positions(mol, functional_groups)
        if fg_positions:
            structural_parts.extend(fg_positions)
        structural_str = "_".join(structural_parts)
        unique_name = f"{base_name}_{structural_str}"
        final_name = unique_name
        suffix = 1
        while final_name in self.used_names:
            final_name = f"{unique_name}_V{suffix}"
            suffix += 1
        self.used_names.add(final_name)
        return final_name

nomenclature_system = UniqueNomenclatureSystem()

def validate_structure_possibility(carbon_count, functional_groups, n_double_bonds, n_triple_bonds, carbon_types, n_rings=0):
    if 'Ether' in functional_groups:
        ether_count = functional_groups.count('Ether')
        if carbon_count < 2:
            return False, f"IMPOSSIBLE: Ethers require at least 2 carbon atoms (you have {carbon_count})"
        if carbon_count == 2:
            if ether_count > 1:
                return False, f"IMPOSSIBLE: For 2 carbon atoms, maximum 1 ether group is possible (you requested {ether_count})"
            if carbon_count == 2 and n_double_bonds > 0 and ether_count > 0:
                if n_double_bonds >= 1 and ether_count >= 1:
                    return False, "IMPOSSIBLE: For 2 carbon atoms, cannot have both double bonds and ethers"
        if carbon_count >= 3:
            max_ethers = carbon_count - 1
            if ether_count > max_ethers:
                return False, f"IMPOSSIBLE: For {carbon_count} carbon atoms, maximum {max_ethers} ether groups are possible (you requested {ether_count})"
            total_unsaturations = n_double_bonds + n_triple_bonds
            if total_unsaturations > 0 and ether_count > 0:
                pass        
    sulfur_chain_types = ['S_Chain_Bi', 'S_Chain_Tetra', 'S_Chain_Hexa']
    sulfur_chains_present = [fg for fg in functional_groups if fg in sulfur_chain_types]
    if sulfur_chains_present:
        total_sulfur_chains = len(sulfur_chains_present)
        if carbon_count < 2:
            return False, f"IMPOSSIBLE: Sulfur chains require at least 2 carbon atoms (you have {carbon_count})"
        if total_sulfur_chains > carbon_count - 1:
            return False, f"IMPOSSIBLE: Maximum sulfur chains for {carbon_count} carbons is {carbon_count - 1} (you requested {total_sulfur_chains})"
    if 'Ketone' in functional_groups:
        if carbon_count < 2:
            return False, f"IMPOSSIBLE: Ketones require at least 2 carbon atoms (you have {carbon_count})"
        if carbon_count == 2:
            if n_double_bonds > 0 or n_triple_bonds > 0:
                return False, "IMPOSSIBLE: For 2 carbons with ketone, no double/triple bonds allowed (would destroy the C-C bond needed for ketone)"
    if 'Azide' in functional_groups:
        azide_count = functional_groups.count('Azide')
        if carbon_count == 0:
            return False, "IMPOSSIBLE: Azides require at least 1 carbon atom (you have 0)"
    if carbon_count > 0:
        total_unsaturations = n_double_bonds + n_triple_bonds
        if n_rings > 0:
            max_possible_unsaturations = carbon_count + n_rings - 1
        else:
            max_possible_unsaturations = carbon_count - 1
        if total_unsaturations > max_possible_unsaturations:
            return False, f"IMPOSSIBLE: For {carbon_count} carbon(s) and {n_rings} ring(s), maximum unsaturations is {max_possible_unsaturations} (you requested {total_unsaturations})"
    else:
        if n_double_bonds > 0 or n_triple_bonds > 0 or n_rings > 0:
            return False, "IMPOSSIBLE: Cannot have double/triple bonds or rings with zero carbons"
    if carbon_count == 0:
        if not functional_groups:
            return False, "IMPOSSIBLE: Zero carbons with no functional groups is nothing"
        valid_zero_carbon_groups = {
            'COOH', 'CHO',  
            'COOR_CH3', 'COOR_C2H5', 'COOR_C3H7', 'COOR_CH(CH3)2',
            'COX_Cl', 'COX_Br', 'COX_F', 'COX_I',  
            'OX_Cl', 'OX_Br', 'OX_F', 'OX_I',  
            'NH2',  
        }
        invalid_groups_for_zero_carbon = {
            'Amide', 'Cl', 'Br', 'I', 'F', 'CN', 'NC', 
            'OCN', 'NCO', 'Imine', 'NO2', 'Ketone', 'Ether', 'OH', 'Azide',
            'S_Bivalent', 'S_Tetravalent', 'S_Hexavalent',
            'S_Chain_Bi', 'S_Chain_Tetra', 'S_Chain_Hexa'
        }
        invalid_groups = []
        for fg in functional_groups:
            if fg in invalid_groups_for_zero_carbon:
                invalid_groups.append(fg)
            elif fg not in valid_zero_carbon_groups:
                invalid_groups.append(fg)
        if invalid_groups:
            return False, f"IMPOSSIBLE: For zero carbons, these functional groups are invalid: {', '.join(set(invalid_groups))}"
    if carbon_count > 0:
        total_fg_count = len(functional_groups)
        if total_fg_count == 0:
            pass
        else:
            if carbon_count == 1:
                if total_fg_count > 4:
                    return False, f"IMPOSSIBLE: Single carbon can have maximum 4 functional groups (you selected {total_fg_count})"
                if 'Ketone' in functional_groups:
                    return False, "IMPOSSIBLE: Ketones require at least 2 carbon atoms (you have 1)"
            else:
                max_fg_estimate = carbon_count * 2.5
                if carbon_count <= 3:
                    max_fg_estimate = carbon_count * 4  
                elif carbon_count <= 6:
                    max_fg_estimate = carbon_count * 3
                if total_fg_count > max_fg_estimate:
                    return False, f"IMPOSSIBLE: {carbon_count} carbons can accommodate approximately {int(max_fg_estimate)} functional groups maximum (you selected {total_fg_count})"
    if carbon_count > 0 and not carbon_types:
        return False, "IMPOSSIBLE: At least one carbon type must be selected when carbon count > 0"
    if n_rings > 0:
        if carbon_count < 3:
            return False, f"IMPOSSIBLE: Need at least 3 carbons to form a ring (you have {carbon_count})"
        if n_rings > carbon_count - 1:
            return False, f"IMPOSSIBLE: Maximum rings for {carbon_count} carbons is {carbon_count - 2}"
    if carbon_count == 1 and 'Ketone' in functional_groups:
        return False, "IMPOSSIBLE: Ketones cannot be formed with only 1 carbon atom"
    return True, ""

def is_duplicate_or_exists(mol_key, output_path):
    if mol_key in DUPLICATE_REGISTRY:
        return True
    if output_path in FILE_EXISTENCE_CACHE:
        return FILE_EXISTENCE_CACHE[output_path]
    if os.path.exists(output_path):
        FILE_EXISTENCE_CACHE[output_path] = True
        return True
    FILE_EXISTENCE_CACHE[output_path] = False
    return False

def sanitize_filename(name):
    translation_table = str.maketrans({c: '_' for c in r'[^a-zA-Z0-9_\-\.]'})
    return name.translate(translation_table)

@lru_cache(maxsize=50000)
def get_degree_sequence_hash(degree_seq):
    return hash(tuple(sorted(degree_seq)))

def canonicalize_smiles(smiles: str) -> str:
    if smiles in _canonical_cache:
        return _canonical_cache[smiles]
    
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        _canonical_cache[smiles] = ""
        return ""
    
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    _canonical_cache[smiles] = canonical
    return canonical

def get_element_counts(mol):
    element_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        element_counts[symbol] = element_counts.get(symbol, 0) + 1
    return element_counts

def get_functional_group_type(mol):
    patterns = {
        'Alcohol': '[OX2H]',
        'Carboxylic_Acid': 'C(=O)[OH]',
        'Aldehyde': '[CX3H1](=O)[#6]',
        'Ketone': '[CX3](=O)[#6]',
        'Ester': 'C(=O)[OX2]',
        'Nitrile': 'C#N',
        'Isocyanide': '[N+]#[C-]',
        'Cyanate': 'O-C#N',
        'Isocyanate': 'N=C=O',
        'Fluoride': '[F]',
        'Chloride': '[Cl]',
        'Bromide': '[Br]',
        'Iodide': '[I]',
        'Nitro': '[N+](=O)[O-]',
        'Ether': '[OD2]([#6])[#6]',
        'Imine': '[CX3]([#6])[#6]=[NX2]',
        'Amide': 'C(=O)[NX3]',
        'Acid_Halide': 'C(=O)[F,Cl,Br,I]',
        'Hypohalite': 'O[F,Cl,Br,I]',
        'Amine': '[NX3;H2,H1;!$(N=*)]',
        'Azide': '[N-]=[N+]=N',
        'S_Bivalent': '[SX2H]',
        'S_Tetravalent': '[SX3](=O)',
        'S_Hexavalent': '[SX4](=O)(=O)',
        'S_Chain_Bi': '[SX2]([#6])[#6]',
        'S_Chain_Tetra': '[SX3](=O)([#6])[#6]',
        'S_Chain_Hexa': '[SX4](=O)(=O)([#6])[#6]'
    }
    detected_groups = []
    for group, smarts in patterns.items():
        try:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                detected_groups.append(group)
        except:
            continue
    if detected_groups:
        return ', '.join(detected_groups)
    else:
        return 'Hydrocarbon'

def get_functional_group_name(func_code):
    func_name_map = {
        'OH': 'Alco',
        'COOH': 'CarboAc',
        'CHO': 'Alde',
        'Ketone': 'Keto',
        'COOR_CH3': 'MeEst',
        'COOR_C2H5': 'EthEst',
        'COOR_C3H7': 'PropEst',
        'COOR_CH(CH3)2': 'IsoPropEst',
        'CN': 'Nitrile',
        'NC': 'IsoCyan',
        'OCN': 'Cyanate',
        'NCO': 'IsoCyanate',
        'F': 'Fluoro',
        'Cl': 'Chloro',
        'Br': 'Bromo',
        'I': 'Iodo',
        'NO2': 'Nitro',
        'Ether': 'Ether',
        'Imine': 'Imine',
        'Amide': 'Amide',
        'COX_Cl': 'AcidChlor',
        'COX_Br': 'AcidBrom',
        'COX_F': 'AcidFluor',
        'COX_I': 'AcidIod',
        'OX_Cl': 'HypoChlor',
        'OX_Br': 'HypoBrom',
        'OX_F': 'HypoFluor',
        'OX_I': 'HypoIod',
        'NH2': 'Amine',
        'Azide': 'Azide', 
        'S_Bivalent': 'S_Bi',
        'S_Tetravalent': 'S_Tetra', 
        'S_Hexavalent': 'S_Hexa',
        'S_Chain_Bi': 'SChain_Bi',
        'S_Chain_Tetra': 'SChain_Tetra',
        'S_Chain_Hexa': 'SChain_Hexa'
    }
    return func_name_map.get(func_code, func_code)

def get_ring_size_distribution(mol):
    """Get the distribution of ring sizes in a molecule"""
    if not mol.GetRingInfo().NumRings():
        return "ACYCLIC"
    
    ring_info = mol.GetRingInfo()
    ring_sizes = []
    for ring in ring_info.AtomRings():
        ring_sizes.append(len(ring))
    
    # Count rings by size
    size_count = defaultdict(int)
    for size in ring_sizes:
        size_count[size] += 1
    
    # Create a simplified classification
    if not size_count:
        return "ACYCLIC"
    
    # If multiple ring sizes, use the smallest one for classification
    smallest_ring = min(size_count.keys())
    return f"R{smallest_ring}"

def create_directory_structure(base_dir, ring_size="ACYCLIC"):
    try:
        # Create ring-specific subdirectory if not acyclic
        if ring_size != "ACYCLIC":
            base_dir = os.path.join(base_dir, f"Ring_Size_{ring_size}")
        
        main_db_dir = os.path.join(os.path.dirname(base_dir), "Main Database")
        mol_dir = os.path.join(base_dir, 'MOL Files')
        sdf_dir = os.path.join(base_dir, 'SDF Files')
        main_mol_dir = os.path.join(main_db_dir, 'MOL Files')
        main_sdf_dir = os.path.join(main_db_dir, 'SDF Files')
        
        os.makedirs(mol_dir, exist_ok=True)
        os.makedirs(sdf_dir, exist_ok=True)
        os.makedirs(main_mol_dir, exist_ok=True)
        os.makedirs(main_sdf_dir, exist_ok=True)
        
        for dir_path in [mol_dir, sdf_dir, main_mol_dir, main_sdf_dir]:
            if not os.path.exists(dir_path):
                print(f"WARNING: Directory not created: {dir_path}")
        
        return {
            'element': {
                'mol_files': mol_dir,
                'sdf_files': sdf_dir
            },
            'main_db': {
                'mol_files': main_mol_dir,
                'sdf_files': main_sdf_dir
            }
        }
    except Exception as e:
        print(f"ERROR creating directory structure: {e}")
        raise

def check_file_permissions(directory):
    test_file = os.path.join(directory, "test_write.permission")
    try:
        with open(test_file, 'w') as f:
            f.write("test")
        os.remove(test_file)
        return True
    except Exception as e:
        print(f"NO WRITE PERMISSION in {directory}: {e}")
        return False

def debug_molecule_processing(mol, unique_compound_id, paths):
    print(f"\n--- Debugging {unique_compound_id} ---")
    print(f"Mol valid: {mol is not None}")
    if mol:
        print(f"Num atoms: {mol.GetNumAtoms()}")
        print(f"Num bonds: {mol.GetNumBonds()}")
        print(f"SMILES: {Chem.MolToSmiles(mol)}")
    print(f"Paths to create:")
    for file_type, path in paths.items():
        print(f"  {file_type}: {path}")
        print(f"  Directory exists: {os.path.exists(os.path.dirname(path))}")
    print("--- End Debug ---\n")

def generate_sulfur_chain_structures(backbone_smiles, sulfur_chain_type, count_needed, carbon_count, n_double_bonds, n_triple_bonds):
    if carbon_count < 2 or count_needed == 0:
        return [backbone_smiles]
    mol = Chem.MolFromSmiles(backbone_smiles)
    if not mol:
        return []
    cc_bonds = []
    for bond in mol.GetBonds():
        if (bond.GetBeginAtom().GetSymbol() == 'C' and 
            bond.GetEndAtom().GetSymbol() == 'C' and
            bond.GetBondType() == BondType.SINGLE):
            cc_bonds.append((bond.GetIdx(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    if len(cc_bonds) < count_needed:
        return []
    all_bond_combinations = list(itertools.combinations(cc_bonds, count_needed))
    sulfur_smiles_set = set()
    for bond_combination in all_bond_combinations:
        try:
            mol_copy = Chem.RWMol(mol)
            sorted_bonds = sorted(bond_combination, key=lambda x: x[0], reverse=True)
            for bond_idx, begin_idx, end_idx in sorted_bonds:
                mol_copy.RemoveBond(begin_idx, end_idx)
                if sulfur_chain_type == 'S_Chain_Bi':
                    s_idx = mol_copy.AddAtom(Chem.Atom('S'))
                    mol_copy.AddBond(begin_idx, s_idx, BondType.SINGLE)
                    mol_copy.AddBond(s_idx, end_idx, BondType.SINGLE)
                elif sulfur_chain_type == 'S_Chain_Tetra':
                    s_idx = mol_copy.AddAtom(Chem.Atom('S'))
                    o_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    mol_copy.AddBond(begin_idx, s_idx, BondType.SINGLE)
                    mol_copy.AddBond(s_idx, o_idx, BondType.DOUBLE)
                    mol_copy.AddBond(s_idx, end_idx, BondType.SINGLE)
                elif sulfur_chain_type == 'S_Chain_Hexa':
                    s_idx = mol_copy.AddAtom(Chem.Atom('S'))
                    o1_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    o2_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    mol_copy.AddBond(begin_idx, s_idx, BondType.SINGLE)
                    mol_copy.AddBond(s_idx, o1_idx, BondType.DOUBLE)
                    mol_copy.AddBond(s_idx, o2_idx, BondType.DOUBLE)
                    mol_copy.AddBond(s_idx, end_idx, BondType.SINGLE)
            modified_mol = mol_copy.GetMol()
            Chem.SanitizeMol(modified_mol)
            final_double_bonds = sum(1 for bond in modified_mol.GetBonds() if bond.GetBondType() == BondType.DOUBLE)
            final_triple_bonds = sum(1 for bond in modified_mol.GetBonds() if bond.GetBondType() == BondType.TRIPLE)
            expected_double_bonds = n_double_bonds
            if sulfur_chain_type == 'S_Chain_Tetra':
                expected_double_bonds += count_needed
            elif sulfur_chain_type == 'S_Chain_Hexa':
                expected_double_bonds += count_needed * 2
            if final_double_bonds == expected_double_bonds and final_triple_bonds == n_triple_bonds:
                sulfur_smiles = Chem.MolToSmiles(modified_mol, isomericSmiles=True, canonical=True)
                sulfur_smiles_set.add(sulfur_smiles)
        except Exception as e:
            print(f"Error creating {count_needed} {sulfur_chain_type} chains: {e}")
            continue
    return list(sulfur_smiles_set) if sulfur_smiles_set else []

def generate_mixed_sulfur_chain_structures(backbone_smiles, sulfur_chain_counts, carbon_count, n_double_bonds, n_triple_bonds):
    if carbon_count < 2 or sum(sulfur_chain_counts.values()) == 0:
        return [backbone_smiles]
    mol = Chem.MolFromSmiles(backbone_smiles)
    if not mol:
        return []
    cc_bonds = []
    for bond in mol.GetBonds():
        if (bond.GetBeginAtom().GetSymbol() == 'C' and 
            bond.GetEndAtom().GetSymbol() == 'C' and
            bond.GetBondType() == BondType.SINGLE):
            cc_bonds.append((bond.GetIdx(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    total_chains_needed = sum(sulfur_chain_counts.values())
    if len(cc_bonds) < total_chains_needed:
        return []
    sulfur_smiles_set = set()
    all_bond_combinations = list(itertools.combinations(cc_bonds, total_chains_needed))
    for bond_combination in all_bond_combinations:
        chain_types = []
        for chain_type, count in sulfur_chain_counts.items():
            chain_types.extend([chain_type] * count)
        unique_assignments = set(itertools.permutations(chain_types))
        for chain_assignment in unique_assignments:
            try:
                mol_copy = Chem.RWMol(mol)
                sorted_bonds = sorted(bond_combination, key=lambda x: x[0], reverse=True)
                for (bond_idx, begin_idx, end_idx), chain_type in zip(sorted_bonds, chain_assignment):
                    mol_copy.RemoveBond(begin_idx, end_idx)
                    if chain_type == 'S_Chain_Bi':
                        s_idx = mol_copy.AddAtom(Chem.Atom('S'))
                        mol_copy.AddBond(begin_idx, s_idx, BondType.SINGLE)
                        mol_copy.AddBond(s_idx, end_idx, BondType.SINGLE)
                    elif chain_type == 'S_Chain_Tetra':
                        s_idx = mol_copy.AddAtom(Chem.Atom('S'))
                        o_idx = mol_copy.AddAtom(Chem.Atom('O'))
                        mol_copy.AddBond(begin_idx, s_idx, BondType.SINGLE)
                        mol_copy.AddBond(s_idx, o_idx, BondType.DOUBLE)
                        mol_copy.AddBond(s_idx, end_idx, BondType.SINGLE)
                    elif chain_type == 'S_Chain_Hexa':
                        s_idx = mol_copy.AddAtom(Chem.Atom('S'))
                        o1_idx = mol_copy.AddAtom(Chem.Atom('O'))
                        o2_idx = mol_copy.AddAtom(Chem.Atom('O'))
                        mol_copy.AddBond(begin_idx, s_idx, BondType.SINGLE)
                        mol_copy.AddBond(s_idx, o1_idx, BondType.DOUBLE)
                        mol_copy.AddBond(s_idx, o2_idx, BondType.DOUBLE)
                        mol_copy.AddBond(s_idx, end_idx, BondType.SINGLE)
                modified_mol = mol_copy.GetMol()
                Chem.SanitizeMol(modified_mol)
                final_double_bonds = sum(1 for bond in modified_mol.GetBonds() if bond.GetBondType() == BondType.DOUBLE)
                final_triple_bonds = sum(1 for bond in modified_mol.GetBonds() if bond.GetBondType() == BondType.TRIPLE)
                expected_double_bonds = n_double_bonds
                for chain_type in chain_assignment:
                    if chain_type == 'S_Chain_Tetra':
                        expected_double_bonds += 1
                    elif chain_type == 'S_Chain_Hexa':
                        expected_double_bonds += 2
                if final_double_bonds == expected_double_bonds and final_triple_bonds == n_triple_bonds:
                    sulfur_smiles = Chem.MolToSmiles(modified_mol, isomericSmiles=True, canonical=True)
                    sulfur_smiles_set.add(sulfur_smiles)
            except Exception as e:
                print(f"Error creating mixed sulfur chains {chain_assignment}: {e}")
                continue
    return list(sulfur_smiles_set) if sulfur_smiles_set else []

def add_sulfur_functional_groups(mol_copy, carbon_idx, func_group):
    if func_group == 'S_Bivalent':
        s_idx = mol_copy.AddAtom(Chem.Atom('S'))
        h_idx = mol_copy.AddAtom(Chem.Atom('H'))
        mol_copy.AddBond(carbon_idx, s_idx, BondType.SINGLE)
        mol_copy.AddBond(s_idx, h_idx, BondType.SINGLE)
    elif func_group == 'S_Tetravalent':
        s_idx = mol_copy.AddAtom(Chem.Atom('S'))
        o_idx = mol_copy.AddAtom(Chem.Atom('O'))
        mol_copy.AddBond(carbon_idx, s_idx, BondType.SINGLE)
        mol_copy.AddBond(s_idx, o_idx, BondType.DOUBLE)
        h_idx = mol_copy.AddAtom(Chem.Atom('H'))
        mol_copy.AddBond(s_idx, h_idx, BondType.SINGLE)
    elif func_group == 'S_Hexavalent':
        s_idx = mol_copy.AddAtom(Chem.Atom('S'))
        o1_idx = mol_copy.AddAtom(Chem.Atom('O'))
        o2_idx = mol_copy.AddAtom(Chem.Atom('O'))
        mol_copy.AddBond(carbon_idx, s_idx, BondType.SINGLE)
        mol_copy.AddBond(s_idx, o1_idx, BondType.DOUBLE)
        mol_copy.AddBond(s_idx, o2_idx, BondType.DOUBLE)

def standalone_process_isomer(args, base_name, output_root_dir, functional_groups, carbon_count, n_double_bonds, n_triple_bonds, n_rings):
    idx, smiles, compound_id = args
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        
        canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        functional_group_type = get_functional_group_type(mol)
        
        # Determine ring size for directory organization
        ring_size = get_ring_size_distribution(mol)
        
        unique_compound_id = nomenclature_system.generate_unique_name(
            mol, formula, functional_groups, carbon_count, n_double_bonds, n_triple_bonds, n_rings, idx
        )
        
        # Pass ring_size to directory creation
        dirs = create_directory_structure(output_root_dir, ring_size)
        
        mol_filename = f"{unique_compound_id}.mol"
        sdf_filename = f"{unique_compound_id}.sdf"
        
        element_mol_path = os.path.join(dirs['element']['mol_files'], mol_filename)
        element_sdf_path = os.path.join(dirs['element']['sdf_files'], sdf_filename)
        
        mol_key = canonical_smiles
        
        if (is_duplicate_or_exists(mol_key, element_mol_path) or
            is_duplicate_or_exists(mol_key, element_sdf_path)):
            return None
        
        try:
            Chem.SanitizeMol(mol)
            AllChem.Compute2DCoords(mol)
            print(f"Processing {unique_compound_id}: {canonical_smiles} | Ring Size: {ring_size}")
            
            try:
                Chem.MolToMolFile(mol, element_mol_path)
                print(f"  ✓ MOL created: {element_mol_path}")
            except Exception as mol_error:
                print(f"  ✗ MOL error: {mol_error}")
                return None
            
            try:
                writer = SDWriter(element_sdf_path)
                if writer:
                    mol.SetProp("_Name", unique_compound_id)
                    mol.SetProp("Formula", formula)
                    mol.SetProp("SMILES", canonical_smiles)
                    mol.SetProp("Functional_Group", functional_group_type)
                    mol.SetProp("Ring_Size", ring_size)  # Add ring size property
                    writer.write(mol)
                    writer.close()
                    print(f"  ✓ SDF created: {element_sdf_path}")
                else:
                    print(f"  ✗ SDF writer failed for: {element_sdf_path}")
                    return None
            except Exception as sdf_error:
                print(f"  ✗ SDF error: {sdf_error}")
                return None
                
        except Exception as e:
            print(f"ERROR saving files for {unique_compound_id}: {e}")
            return None
        
        # Also copy to main database with ring-based organization
        main_mol_path = os.path.join(dirs['main_db']['mol_files'], mol_filename)
        main_sdf_path = os.path.join(dirs['main_db']['sdf_files'], sdf_filename)
        
        if mol_key in DUPLICATE_REGISTRY:
            return None
        
        try:
            shutil.copy2(element_mol_path, main_mol_path)
            shutil.copy2(element_sdf_path, main_sdf_path)
            print(f"  ✓ Files copied to main database")
        except Exception as e:
            print(f"ERROR copying to main database: {e}")
            return None
        
        DUPLICATE_REGISTRY[mol_key] = True
        
        # Clear memory by deleting the molecule object
        del mol
        
        return {
            'index': idx,
            'formula': formula,
            'smiles': canonical_smiles,
            'element_counts': get_element_counts(Chem.MolFromSmiles(canonical_smiles)),
            'functional_group_type': functional_group_type,
            'compound_id': unique_compound_id,
            'ring_size': ring_size,
            'paths': {
                'element': {
                    'mol': element_mol_path,
                    'sdf': element_sdf_path
                },
                'main_db': {
                    'mol': main_mol_path,
                    'sdf': main_sdf_path
                }
            }
        }
        
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {e}")
        return None

def add_unsaturations(backbone_smiles, n_double_bonds, n_triple_bonds):
    mol = Chem.MolFromSmiles(backbone_smiles)
    if not mol:
        return []
    cc_single_bonds = []
    for bond in mol.GetBonds():
        if (bond.GetBeginAtom().GetSymbol() == 'C' and 
            bond.GetEndAtom().GetSymbol() == 'C' and 
            bond.GetBondType() == BondType.SINGLE):
            cc_single_bonds.append(bond.GetIdx())
    total_unsaturations = n_double_bonds + n_triple_bonds
    if len(cc_single_bonds) < total_unsaturations:
        return []
    all_unsaturation_combinations = list(itertools.combinations(cc_single_bonds, total_unsaturations))
    unsaturated_smiles = set()
    for combination in all_unsaturation_combinations:
        for double_bond_indices in itertools.combinations(combination, n_double_bonds):
            triple_bond_indices = [idx for idx in combination if idx not in double_bond_indices]
            mol_copy = Chem.RWMol(mol)
            for bond_idx in double_bond_indices:
                bond = mol_copy.GetBondWithIdx(bond_idx)
                bond.SetBondType(BondType.DOUBLE)
            for bond_idx in triple_bond_indices:
                bond = mol_copy.GetBondWithIdx(bond_idx)
                bond.SetBondType(BondType.TRIPLE)
            try:
                modified_mol = mol_copy.GetMol()
                Chem.SanitizeMol(modified_mol)
                can = Chem.MolToSmiles(modified_mol, isomericSmiles=True, canonical=True)
                unsaturated_smiles.add(can)
            except:
                continue
    return list(unsaturated_smiles)

def generate_alkane_backbone(n_carbons):
    if n_carbons == 0:
        return [""]
    elif n_carbons == 1:
        return ["C"]
    elif n_carbons == 2:
        return ["CC"]
    elif n_carbons == 3:
        return ["CCC", "C(C)C"]
    elif n_carbons == 4:
        return ["CCCC", "C(C)CC", "CC(C)C", "C(C)(C)C"]
    backbones = set()
    trees = list(nx.nonisomorphic_trees(n_carbons))
    for tree in trees:
        mol = Chem.RWMol()
        node_map = {}
        for node in tree.nodes():
            node_map[node] = mol.AddAtom(Chem.Atom('C'))
        for edge in tree.edges():
            mol.AddBond(node_map[edge[0]], node_map[edge[1]], BondType.SINGLE)
        try:
            mol = mol.GetMol()
            Chem.SanitizeMol(mol)
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            backbones.add(smiles)
        except:
            continue
    return list(backbones)

def generate_hydrocarbon_isomers(n_carbons, n_double_bonds=0, n_triple_bonds=0, n_rings=0):
    if n_rings == 0:
        if n_carbons == 0:
            if n_double_bonds == 0 and n_triple_bonds == 0:
                return [""]
            else:
                return []
        if n_carbons == 1:
            if n_double_bonds == 1 and n_triple_bonds == 0:
                return ["C=C"]
            elif n_double_bonds == 0 and n_triple_bonds == 1:
                return ["C#C"]
            elif n_double_bonds == 0 and n_triple_bonds == 0:
                return ["C"]
            else:
                return []
        if n_double_bonds > 0 or n_triple_bonds > 0:
            saturated_backbones = generate_alkane_backbone(n_carbons)
            unsaturated_backbones = []
            for backbone in saturated_backbones:
                exact_unsaturated = add_exact_unsaturations(backbone, n_double_bonds, n_triple_bonds)
                unsaturated_backbones.extend(exact_unsaturated)
            backbones = list(set(unsaturated_backbones))
        else:
            backbones = generate_alkane_backbone(n_carbons)
        return backbones
    else:
        return generate_cyclic_isomers(n_carbons, n_rings, n_double_bonds, n_triple_bonds)

def add_exact_unsaturations(backbone_smiles, n_double_bonds, n_triple_bonds):
    if n_double_bonds == 0 and n_triple_bonds == 0:
        return [backbone_smiles]
    mol = Chem.MolFromSmiles(backbone_smiles)
    if not mol:
        return []
    cc_single_bonds = []
    for bond in mol.GetBonds():
        if (bond.GetBeginAtom().GetSymbol() == 'C' and 
            bond.GetEndAtom().GetSymbol() == 'C' and 
            bond.GetBondType() == BondType.SINGLE):
            cc_single_bonds.append(bond.GetIdx())
    total_unsaturations = n_double_bonds + n_triple_bonds
    if len(cc_single_bonds) < total_unsaturations:
        return []
    all_unsaturation_combinations = list(itertools.combinations(cc_single_bonds, total_unsaturations))
    exact_unsaturated_smiles = set()
    for combination in all_unsaturation_combinations:
        for double_bond_indices in itertools.combinations(combination, n_double_bonds):
            triple_bond_indices = [idx for idx in combination if idx not in double_bond_indices]
            mol_copy = Chem.RWMol(mol)
            for bond_idx in double_bond_indices:
                bond = mol_copy.GetBondWithIdx(bond_idx)
                bond.SetBondType(BondType.DOUBLE)
            for bond_idx in triple_bond_indices:
                bond = mol_copy.GetBondWithIdx(bond_idx)
                bond.SetBondType(BondType.TRIPLE)
            try:
                modified_mol = mol_copy.GetMol()
                Chem.SanitizeMol(modified_mol)
                final_double_bonds = sum(1 for bond in modified_mol.GetBonds() if bond.GetBondType() == BondType.DOUBLE)
                final_triple_bonds = sum(1 for bond in modified_mol.GetBonds() if bond.GetBondType() == BondType.TRIPLE)
                if final_double_bonds == n_double_bonds and final_triple_bonds == n_triple_bonds:
                    can = Chem.MolToSmiles(modified_mol, isomericSmiles=True, canonical=True)
                    exact_unsaturated_smiles.add(can)
            except:
                continue
    return list(exact_unsaturated_smiles)

def cyclomatic_number(G):
    return G.number_of_edges() - G.number_of_nodes() + nx.number_connected_components(G)

def graph_to_smiles_from_nx(G):
    mol = Chem.RWMol()
    node_map = {}
    for n in G.nodes():
        node_map[n] = mol.AddAtom(Chem.Atom('C'))
    for u, v in G.edges():
        mol.AddBond(node_map[u], node_map[v], BondType.SINGLE)
    try:
        m = mol.GetMol()
        Chem.SanitizeMol(m)
        return Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
    except:
        return ""

def generate_cyclic_isomers(n_carbons, n_rings, n_double_bonds=0, n_triple_bonds=0):
    smiles_set = set()
    try:
        if n_carbons <= 7:
            for g in nx.graph_atlas_g():
                if (g.number_of_nodes() == n_carbons and 
                    nx.is_connected(g) and 
                    cyclomatic_number(g) == n_rings):
                    smi = graph_to_smiles_from_nx(g)
                    if smi:
                        mol = Chem.MolFromSmiles(smi)
                        if mol and sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'C') == n_carbons:
                            if n_double_bonds > 0 or n_triple_bonds > 0:
                                unsaturated_variants = add_unsaturations(smi, n_double_bonds, n_triple_bonds)
                                for usmi in unsaturated_variants:
                                    smiles_set.add(usmi)
                            else:
                                smiles_set.add(smi)
        else:
            trees = list(nx.nonisomorphic_trees(n_carbons))
            for tree in trees:
                G = tree.copy()
                possible_edges = []
                for i in range(n_carbons):
                    for j in range(i+1, n_carbons):
                        if not G.has_edge(i, j):
                            possible_edges.append((i, j))
                for edges_to_add in itertools.combinations(possible_edges, n_rings):
                    G_copy = G.copy()
                    for edge in edges_to_add:
                        G_copy.add_edge(edge[0], edge[1])
                    if (nx.is_connected(G_copy) and 
                        cyclomatic_number(G_copy) == n_rings):
                        smi = graph_to_smiles_from_nx(G_copy)
                        if smi:
                            mol = Chem.MolFromSmiles(smi)
                            if mol and sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'C') == n_carbons:
                                if n_double_bonds > 0 or n_triple_bonds > 0:
                                    unsaturated_variants = add_unsaturations(smi, n_double_bonds, n_triple_bonds)
                                    for usmi in unsaturated_variants:
                                        smiles_set.add(usmi)
                                else:
                                    smiles_set.add(smi)
    except Exception as e:
        print(f"Error generating cyclic isomers: {e}")
    return list(smiles_set)

def generate_ether_structures(backbone_smiles, carbon_count, n_double_bonds, n_triple_bonds, ether_count=1):
    if carbon_count < 2:
        return []
    mol = Chem.MolFromSmiles(backbone_smiles)
    if not mol:
        return []
    ether_smiles_set = set()
    current_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == BondType.DOUBLE)
    current_triple_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == BondType.TRIPLE)
    if current_double_bonds != n_double_bonds or current_triple_bonds != n_triple_bonds:
        return []
    cc_single_bonds = []
    for bond in mol.GetBonds():
        if (bond.GetBeginAtom().GetSymbol() == 'C' and 
            bond.GetEndAtom().GetSymbol() == 'C' and
            bond.GetBondType() == BondType.SINGLE):
            cc_single_bonds.append((bond.GetIdx(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    if ether_count > len(cc_single_bonds):
        return []
    all_bond_combinations = list(itertools.combinations(cc_single_bonds, ether_count))
    for bond_combination in all_bond_combinations:
        try:
            mol_copy = Chem.RWMol(mol)
            sorted_bonds = sorted(bond_combination, key=lambda x: x[0], reverse=True)
            for bond_idx, begin_idx, end_idx in sorted_bonds:
                mol_copy.RemoveBond(begin_idx, end_idx)
                oxygen_idx = mol_copy.AddAtom(Chem.Atom('O'))
                mol_copy.AddBond(begin_idx, oxygen_idx, BondType.SINGLE)
                mol_copy.AddBond(oxygen_idx, end_idx, BondType.SINGLE)
            modified_mol = mol_copy.GetMol()
            Chem.SanitizeMol(modified_mol)
            final_double_bonds = sum(1 for bond in modified_mol.GetBonds() if bond.GetBondType() == BondType.DOUBLE)
            final_triple_bonds = sum(1 for bond in modified_mol.GetBonds() if bond.GetBondType() == BondType.TRIPLE)
            if final_double_bonds == n_double_bonds and final_triple_bonds == n_triple_bonds:
                ether_smiles = Chem.MolToSmiles(modified_mol, isomericSmiles=True, canonical=True)
                ether_smiles_set.add(ether_smiles)
        except Exception as e:
            print(f"Error creating {ether_count} ether(s) from bonds {[b[0] for b in bond_combination]}: {e}")
            continue
    return list(ether_smiles_set)

def add_functional_groups(backbone_smiles, functional_groups, carbon_types=None, n_double_bonds=0, n_triple_bonds=0):
    if carbon_types is None:
        carbon_types = ['primary', 'secondary', 'tertiary']
    if not functional_groups:
        return [backbone_smiles]
    if backbone_smiles == "":
        functional_smiles_list = []
        func_group_smiles = {
            'COOH': 'C(=O)O',
            'CHO': 'C=O',
            'COOR_CH3': 'C(=O)OC',
            'COOR_C2H5': 'C(=O)OCC',
            'COOR_C3H7': 'C(=O)OCCC',
            'COOR_CH(CH3)2': 'C(=O)OC(C)C',
            'COX_Cl': 'C(=O)Cl',
            'COX_Br': 'C(=O)Br',
            'COX_F': 'C(=O)F',
            'COX_I': 'C(=O)I',
            'OX_Cl': 'OCl',
            'OX_Br': 'OBr',
            'OX_F': 'OF',
            'OX_I': 'OI',
            'NH2': 'N',
        }
        for func_group in functional_groups:
            if func_group in func_group_smiles:
                functional_smiles_list.append(func_group_smiles[func_group])
        return functional_smiles_list
    mol = Chem.MolFromSmiles(backbone_smiles)
    if not mol:
        return []
    backbone_carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if 'Ether' in functional_groups and backbone_carbon_count >= 2:
        ether_count = functional_groups.count('Ether')
        ether_structures = generate_ether_structures(
            backbone_smiles, backbone_carbon_count, n_double_bonds, n_triple_bonds, ether_count
        )
        remaining_groups = [fg for fg in functional_groups if fg != 'Ether']
        if remaining_groups:
            all_functionalized = []
            for ether_smiles in ether_structures:
                functionalized = add_functional_groups(ether_smiles, remaining_groups, carbon_types, n_double_bonds, n_triple_bonds)
                all_functionalized.extend(functionalized)
            return all_functionalized
        else:
            return ether_structures
    sulfur_chain_types = ['S_Chain_Bi', 'S_Chain_Tetra', 'S_Chain_Hexa']
    sulfur_chains_present = [fg for fg in functional_groups if fg in sulfur_chain_types]
    if sulfur_chains_present and backbone_carbon_count >= 2:
        sulfur_chain_counts = {}
        for chain_type in sulfur_chain_types:
            count = functional_groups.count(chain_type)
            if count > 0:
                sulfur_chain_counts[chain_type] = count
        total_sulfur_chains = sum(sulfur_chain_counts.values())
        if len(sulfur_chain_counts) > 1:
            mixed_structures = generate_mixed_sulfur_chain_structures(
                backbone_smiles, sulfur_chain_counts, backbone_carbon_count, n_double_bonds, n_triple_bonds
            )
            current_structures = mixed_structures
        else:
            chain_type = list(sulfur_chain_counts.keys())[0]
            count_needed = sulfur_chain_counts[chain_type]
            current_structures = generate_sulfur_chain_structures(
                backbone_smiles, chain_type, count_needed, backbone_carbon_count, n_double_bonds, n_triple_bonds
            )
        if not current_structures:
            return []
        remaining_groups = [fg for fg in functional_groups if fg not in sulfur_chain_types]
        if remaining_groups:
            all_functionalized = []
            for sulfur_smiles in current_structures:
                functionalized = add_functional_groups(sulfur_smiles, remaining_groups, carbon_types, n_double_bonds, n_triple_bonds)
                all_functionalized.extend(functionalized)
            return all_functionalized
        else:
            return current_structures
    ketone_count = functional_groups.count('Ketone')
    if ketone_count > 0 and backbone_carbon_count >= 2:
        ketone_smiles_set = set()
        cc_bonds = []
        for bond in mol.GetBonds():
            if (bond.GetBeginAtom().GetSymbol() == 'C' and 
                bond.GetEndAtom().GetSymbol() == 'C'):
                cc_bonds.append(bond.GetIdx())
        for bond_idx in cc_bonds:
            try:
                mol_copy = Chem.RWMol(mol)
                bond = mol_copy.GetBondWithIdx(bond_idx)
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                mol_copy.RemoveBond(begin_idx, end_idx)
                carbonyl_idx = mol_copy.AddAtom(Chem.Atom('C'))
                oxygen_idx = mol_copy.AddAtom(Chem.Atom('O'))
                mol_copy.AddBond(begin_idx, carbonyl_idx, BondType.SINGLE)
                mol_copy.AddBond(carbonyl_idx, oxygen_idx, BondType.DOUBLE)
                mol_copy.AddBond(carbonyl_idx, end_idx, BondType.SINGLE)
                remaining_groups = functional_groups.copy()
                remaining_groups.remove('Ketone')
                if remaining_groups:
                    ketone_smiles = Chem.MolToSmiles(mol_copy.GetMol())
                    mixed_compounds = add_functional_groups(ketone_smiles, remaining_groups, carbon_types)
                    ketone_smiles_set.update(mixed_compounds)
                else:
                    try:
                        modified_mol = mol_copy.GetMol()
                        Chem.SanitizeMol(modified_mol)
                        can = Chem.MolToSmiles(modified_mol)
                        ketone_smiles_set.add(can)
                    except:
                        continue
            except Exception as e:
                print(f"Error creating ketone: {e}")
                continue
        return list(ketone_smiles_set)
    if not functional_groups:
        return []
    carbon_atoms = []
    carbon_valence = {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            if len(mol.GetAtoms()) == 1:
                carbon_type = 'primary'
            else:
                carbon_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'C')
                if carbon_neighbors == 1:
                    carbon_type = 'primary'
                elif carbon_neighbors == 2:
                    carbon_type = 'secondary'
                elif carbon_neighbors == 3:
                    carbon_type = 'tertiary'
                else:
                    carbon_type = 'primary'
            if carbon_type in carbon_types:
                carbon_atoms.append(atom.GetIdx())
                carbon_valence[atom.GetIdx()] = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
    if not carbon_atoms or not functional_groups:
        return []
    functional_smiles_list = []
    for carbon_assignment in itertools.product(carbon_atoms, repeat=len(functional_groups)):
        mol_copy = Chem.RWMol(mol)
        skip_this_assignment = False
        carbon_valence_needed = {idx: 0 for idx in carbon_atoms}
        for i, carbon_idx in enumerate(carbon_assignment):
            func_group = functional_groups[i]
            if func_group in ['OH', 'COOH', 'CHO', 'CN', 'NC', 'OCN', 'NCO', 'F', 'Cl', 'Br', 'I', 
                             'NO2', 'Amide', 'NH2', 'S_Bivalent', 'S_Tetravalent', 'S_Hexavalent', 'Azide']:
                valence_needed = 1
            elif func_group == 'Imine':
                valence_needed = 2
            else:
                valence_needed = 1
            carbon_valence_needed[carbon_idx] += valence_needed
            if carbon_valence[carbon_idx] + carbon_valence_needed[carbon_idx] > 4:
                skip_this_assignment = True
                break
        if skip_this_assignment:
            continue
        try:
            for i, carbon_idx in enumerate(carbon_assignment):
                func_group = functional_groups[i]
                carbon_atom = mol_copy.GetAtomWithIdx(carbon_idx)
                current_valence = sum(bond.GetBondTypeAsDouble() for bond in carbon_atom.GetBonds())
                if func_group in ['OH', 'COOH', 'CHO', 'CN', 'NC', 'OCN', 'NCO', 'F', 'Cl', 'Br', 'I', 
                                 'NO2', 'Amide', 'NH2', 'S_Bivalent', 'S_Tetravalent', 'S_Hexavalent', 'Azide']:
                    valence_needed = 1
                elif func_group == 'Imine':
                    valence_needed = 2
                else:
                    valence_needed = 1
                if current_valence + valence_needed > 4:
                    raise ValueError(f"Cannot add {func_group} to carbon {carbon_idx}: valence would exceed 4")
                if func_group == 'OH':
                    o_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    h_idx = mol_copy.AddAtom(Chem.Atom('H'))
                    mol_copy.AddBond(carbon_idx, o_idx, BondType.SINGLE)
                    mol_copy.AddBond(o_idx, h_idx, BondType.SINGLE)
                elif func_group == 'COOH':
                    c_idx = mol_copy.AddAtom(Chem.Atom('C'))
                    o1_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    o2_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    h_idx = mol_copy.AddAtom(Chem.Atom('H'))
                    mol_copy.AddBond(carbon_idx, c_idx, BondType.SINGLE)
                    mol_copy.AddBond(c_idx, o1_idx, BondType.DOUBLE)
                    mol_copy.AddBond(c_idx, o2_idx, BondType.SINGLE)
                    mol_copy.AddBond(o2_idx, h_idx, BondType.SINGLE)
                elif func_group == 'CHO':
                    c_idx = mol_copy.AddAtom(Chem.Atom('C'))
                    o_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    h_idx = mol_copy.AddAtom(Chem.Atom('H'))
                    mol_copy.AddBond(carbon_idx, c_idx, BondType.SINGLE)
                    mol_copy.AddBond(c_idx, o_idx, BondType.DOUBLE)
                    mol_copy.AddBond(c_idx, h_idx, BondType.SINGLE)
                elif func_group.startswith('COOR'):
                    r_group = func_group.split('_')[1] if '_' in func_group else 'CH3'
                    c_idx = mol_copy.AddAtom(Chem.Atom('C'))
                    o1_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    o2_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    mol_copy.AddBond(carbon_idx, c_idx, BondType.SINGLE)
                    mol_copy.AddBond(c_idx, o1_idx, BondType.DOUBLE)
                    mol_copy.AddBond(c_idx, o2_idx, BondType.SINGLE)
                    if r_group == 'CH3':
                        r_c_idx = mol_copy.AddAtom(Chem.Atom('C'))
                        mol_copy.AddBond(o2_idx, r_c_idx, BondType.SINGLE)
                    elif r_group == 'C2H5':
                        r_c1_idx = mol_copy.AddAtom(Chem.Atom('C'))
                        r_c2_idx = mol_copy.AddAtom(Chem.Atom('C'))
                        mol_copy.AddBond(o2_idx, r_c1_idx, BondType.SINGLE)
                        mol_copy.AddBond(r_c1_idx, r_c2_idx, BondType.SINGLE)
                    elif r_group == 'C3H7':
                        r_c1_idx = mol_copy.AddAtom(Chem.Atom('C'))
                        r_c2_idx = mol_copy.AddAtom(Chem.Atom('C'))
                        r_c3_idx = mol_copy.AddAtom(Chem.Atom('C'))
                        mol_copy.AddBond(o2_idx, r_c1_idx, BondType.SINGLE)
                        mol_copy.AddBond(r_c1_idx, r_c2_idx, BondType.SINGLE)
                        mol_copy.AddBond(r_c2_idx, r_c3_idx, BondType.SINGLE)
                    elif r_group in ['CH(CH3)2', '-C(CH3)2']:
                        r_c1_idx = mol_copy.AddAtom(Chem.Atom('C'))
                        r_c2_idx = mol_copy.AddAtom(Chem.Atom('C'))
                        r_c3_idx = mol_copy.AddAtom(Chem.Atom('C'))
                        mol_copy.AddBond(o2_idx, r_c1_idx, BondType.SINGLE)
                        mol_copy.AddBond(r_c1_idx, r_c2_idx, BondType.SINGLE)
                        mol_copy.AddBond(r_c1_idx, r_c3_idx, BondType.SINGLE)
                elif func_group == 'CN':
                    c_idx = mol_copy.AddAtom(Chem.Atom('C'))
                    n_idx = mol_copy.AddAtom(Chem.Atom('N'))
                    mol_copy.AddBond(carbon_idx, c_idx, BondType.SINGLE)
                    mol_copy.AddBond(c_idx, n_idx, BondType.TRIPLE)
                elif func_group == 'NC':
                    n_idx = mol_copy.AddAtom(Chem.Atom('N'))
                    c_idx = mol_copy.AddAtom(Chem.Atom('C'))
                    mol_copy.AddBond(carbon_idx, n_idx, BondType.SINGLE)
                    mol_copy.AddBond(n_idx, c_idx, BondType.TRIPLE)
                    mol_copy.GetAtomWithIdx(n_idx).SetFormalCharge(1)
                    mol_copy.GetAtomWithIdx(c_idx).SetFormalCharge(-1)
                elif func_group == 'OCN':
                    o_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    c_idx = mol_copy.AddAtom(Chem.Atom('C'))
                    n_idx = mol_copy.AddAtom(Chem.Atom('N'))
                    mol_copy.AddBond(carbon_idx, o_idx, BondType.SINGLE)
                    mol_copy.AddBond(o_idx, c_idx, BondType.SINGLE)
                    mol_copy.AddBond(c_idx, n_idx, BondType.TRIPLE)
                elif func_group == 'NCO':
                    n_idx = mol_copy.AddAtom(Chem.Atom('N'))
                    c_idx = mol_copy.AddAtom(Chem.Atom('C'))
                    o_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    mol_copy.AddBond(carbon_idx, n_idx, BondType.SINGLE)
                    mol_copy.AddBond(n_idx, c_idx, BondType.DOUBLE)
                    mol_copy.AddBond(c_idx, o_idx, BondType.DOUBLE)
                elif func_group in ['F', 'Cl', 'Br', 'I']:
                    halogen_idx = mol_copy.AddAtom(Chem.Atom(func_group))
                    mol_copy.AddBond(carbon_idx, halogen_idx, BondType.SINGLE)
                elif func_group == 'NO2':
                    n_idx = mol_copy.AddAtom(Chem.Atom('N'))
                    o1_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    o2_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    mol_copy.AddBond(carbon_idx, n_idx, BondType.SINGLE)
                    mol_copy.AddBond(n_idx, o1_idx, BondType.DOUBLE)
                    mol_copy.AddBond(n_idx, o2_idx, BondType.SINGLE)
                    mol_copy.GetAtomWithIdx(n_idx).SetFormalCharge(1)
                    mol_copy.GetAtomWithIdx(o1_idx).SetFormalCharge(0)
                    mol_copy.GetAtomWithIdx(o2_idx).SetFormalCharge(-1)
                elif func_group == 'Imine':
                    n_idx = mol_copy.AddAtom(Chem.Atom('N'))
                    mol_copy.AddBond(carbon_idx, n_idx, BondType.DOUBLE)
                elif func_group == 'Amide':
                    c_idx = mol_copy.AddAtom(Chem.Atom('C'))
                    o_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    n_idx = mol_copy.AddAtom(Chem.Atom('N'))
                    h_idx = mol_copy.AddAtom(Chem.Atom('H'))
                    mol_copy.AddBond(carbon_idx, c_idx, BondType.SINGLE)
                    mol_copy.AddBond(c_idx, o_idx, BondType.DOUBLE)
                    mol_copy.AddBond(c_idx, n_idx, BondType.SINGLE)
                    mol_copy.AddBond(n_idx, h_idx, BondType.SINGLE)
                elif func_group.startswith('COX'):
                    halogen = func_group.split('_')[1] if '_' in func_group else 'Cl'
                    c_idx = mol_copy.AddAtom(Chem.Atom('C'))
                    o_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    x_idx = mol_copy.AddAtom(Chem.Atom(halogen))
                    mol_copy.AddBond(carbon_idx, c_idx, BondType.SINGLE)
                    mol_copy.AddBond(c_idx, o_idx, BondType.DOUBLE)
                    mol_copy.AddBond(c_idx, x_idx, BondType.SINGLE)
                elif func_group.startswith('OX'):
                    halogen = func_group.split('_')[1] if '_' in func_group else 'Cl'
                    o_idx = mol_copy.AddAtom(Chem.Atom('O'))
                    x_idx = mol_copy.AddAtom(Chem.Atom(halogen))
                    mol_copy.AddBond(carbon_idx, o_idx, BondType.SINGLE)
                    mol_copy.AddBond(o_idx, x_idx, BondType.SINGLE)
                elif func_group == 'NH2':
                    n_idx = mol_copy.AddAtom(Chem.Atom('N'))
                    h1_idx = mol_copy.AddAtom(Chem.Atom('H'))
                    h2_idx = mol_copy.AddAtom(Chem.Atom('H'))
                    mol_copy.AddBond(carbon_idx, n_idx, BondType.SINGLE)
                    mol_copy.AddBond(n_idx, h1_idx, BondType.SINGLE)
                    mol_copy.AddBond(n_idx, h2_idx, BondType.SINGLE)
                elif func_group == 'Azide':
                    n1_idx = mol_copy.AddAtom(Chem.Atom('N'))
                    n2_idx = mol_copy.AddAtom(Chem.Atom('N'))
                    n3_idx = mol_copy.AddAtom(Chem.Atom('N'))
                    mol_copy.AddBond(carbon_idx, n1_idx, BondType.SINGLE)
                    mol_copy.AddBond(n1_idx, n2_idx, BondType.DOUBLE)
                    mol_copy.AddBond(n2_idx, n3_idx, BondType.SINGLE)
                elif func_group in ['S_Bivalent', 'S_Tetravalent', 'S_Hexavalent']:
                    add_sulfur_functional_groups(mol_copy, carbon_idx, func_group)
            modified_mol = mol_copy.GetMol()
            Chem.SanitizeMol(modified_mol)
            canonical_smiles = Chem.MolToSmiles(modified_mol)
            if canonical_smiles not in functional_smiles_list:
                functional_smiles_list.append(canonical_smiles)
        except Exception as e:
            print(f"Error adding functional group {func_group} to carbon {carbon_idx}: {e}")
            continue
    return functional_smiles_list

def generate_functionalized_isomers(n_carbons, functional_groups, n_double_bonds=0, n_triple_bonds=0, n_rings=0, carbon_types=None):
    if carbon_types is None:
        carbon_types = ['primary', 'secondary', 'tertiary']
    ether_count = functional_groups.count('Ether')
    if ether_count > 0:
        if n_carbons < 2:
            functional_groups = [fg for fg in functional_groups if fg != 'Ether']
        elif n_carbons == 2 and ether_count > 1:
            ethers_kept = 0
            new_groups = []
            for fg in functional_groups:
                if fg == 'Ether':
                    if ethers_kept < 1:
                        new_groups.append(fg)
                        ethers_kept += 1
                else:
                    new_groups.append(fg)
            functional_groups = new_groups
    if 'Ketone' in functional_groups and n_carbons < 2:
        functional_groups = [fg for fg in functional_groups if fg != 'Ketone']
    sulfur_chain_types = ['S_Chain_Bi', 'S_Chain_Tetra', 'S_Chain_Hexa']
    for sulfur_chain in sulfur_chain_types:
        if sulfur_chain in functional_groups and n_carbons < 2:
            functional_groups = [fg for fg in functional_groups if fg != sulfur_chain]
    backbones = generate_hydrocarbon_isomers(n_carbons, n_double_bonds, n_triple_bonds, n_rings)
    all_functional_smiles = []
    for backbone in backbones:
        functionalized = add_functional_groups(backbone, functional_groups, carbon_types, n_double_bonds, n_triple_bonds)
        all_functional_smiles.extend(functionalized)
    return list(set(all_functional_smiles))

def process_alkane_batch_parallel(args):
    batch_data, batch_id, base_name, output_root_dir, functional_groups, carbon_count, n_double_bonds, n_triple_bonds, n_rings = args
    batch_results = []
    for task in batch_data:
        result = standalone_process_isomer(task, base_name, output_root_dir, functional_groups, carbon_count, n_double_bonds, n_triple_bonds, n_rings)
        if result:
            batch_results.append(result)
    return batch_results

def create_alkane_batches(tasks, num_processes, base_name, output_root_dir, functional_groups, carbon_count, n_double_bonds, n_triple_bonds, n_rings):
    if not tasks:
        return []
    total_tasks = len(tasks)
    num_batches = num_processes * 4
    batch_size = max(1, total_tasks // num_batches)
    if batch_size < 5 and total_tasks > 50:
        batch_size = max(5, total_tasks // (num_processes * 2))
    batches = []
    for i in range(0, total_tasks, batch_size):
        batch = tasks[i:i + batch_size]
        if batch:
            batches.append((batch, len(batches), base_name, output_root_dir, functional_groups, carbon_count, n_double_bonds, n_triple_bonds, n_rings))
    return batches

def parallel_alkane_generation(tasks, base_name, output_root_dir, functional_groups, carbon_count, n_double_bonds, n_triple_bonds, n_rings, progress_callback=None, num_processes=None):
    if num_processes is None:
        num_processes = cpu_count()
        if num_processes >= 8:
            num_processes = num_processes
        else:
            num_processes = max(2, num_processes - 1)
    print(f"Using {num_processes} parallel processes for functionalized compound generation")
    batches = create_alkane_batches(tasks, num_processes, base_name, output_root_dir, functional_groups, carbon_count, n_double_bonds, n_triple_bonds, n_rings)
    if not batches:
        return []
    all_results = []
    processed_batches = 0
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        future_to_batch = {executor.submit(process_alkane_batch_parallel, batch): batch 
                           for batch in batches}
        for future in as_completed(future_to_batch):
            try:
                batch_results = future.result()
                all_results.extend(batch_results)
                processed_batches += 1
                if progress_callback:
                    progress = (processed_batches / len(batches)) * 100
                    progress_callback(int(progress), processed_batches, len(batches))
            except Exception as e:
                print(f"Error processing batch: {e}")
                continue
    return all_results

class UltraFastCompoundGeneratorGUI:
    def __init__(self, root):
        self.root = root
        root.title("CHEM-∞ : The Universal Molecular Generator")
        root.state('zoomed')
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()
        root.geometry(f"{screen_width}x{screen_height}+0+0")
        self.mainframe = ttk.Frame(root, padding=15)
        self.mainframe.pack(fill=tk.BOTH, expand=True)
        canvas = tk.Canvas(self.mainframe)
        scrollbar = ttk.Scrollbar(self.mainframe, orient="vertical", command=canvas.yview)
        self.scrollable_frame = ttk.Frame(canvas)
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        cpu_count_val = cpu_count()
        memory_gb = psutil.virtual_memory().total / (1024**3)
        system_info = f"System: {cpu_count_val} CPU cores, {memory_gb:.1f}GB RAM"
        ttk.Label(self.scrollable_frame, text=system_info, font=('Arial', 12, 'bold'), foreground="darkgreen").grid(row=0, column=0, columnspan=6, sticky=tk.W, pady=(0, 5))
        title_label = ttk.Label(self.scrollable_frame, text="CHEM-∞ : The Universal Molecular Generator", 
                               font=('Arial', 14, 'bold'), foreground="blue")
        title_label.grid(row=1, column=0, columnspan=6, pady=(0, 5))
        self.scrollable_frame.grid_columnconfigure(0, weight=1)
        self.scrollable_frame.grid_columnconfigure(1, weight=1)
        self.scrollable_frame.grid_columnconfigure(2, weight=1)
        self.scrollable_frame.grid_columnconfigure(3, weight=1)
        self.scrollable_frame.grid_columnconfigure(4, weight=1)
        self.scrollable_frame.grid_columnconfigure(5, weight=1)
        ttk.Label(self.scrollable_frame, text="Number of Carbons (0-20):", font=('Arial', 12, 'bold')).grid(row=2, column=0, sticky=tk.W, pady=2)
        self.carbon_spinbox = ttk.Spinbox(self.scrollable_frame, from_=0, to=20, width=8, font=('Arial', 12))
        self.carbon_spinbox.set(0)
        self.carbon_spinbox.grid(row=2, column=1, sticky=tk.W, pady=2)
        self.carbon_spinbox.bind('<FocusOut>', self.update_functional_group_states)
        self.carbon_spinbox.bind('<ButtonRelease>', self.update_functional_group_states)
        ttk.Label(self.scrollable_frame, text="Number of Double Bonds:", font=('Arial', 12, 'bold')).grid(row=3, column=0, sticky=tk.W, pady=2)
        self.double_bond_spinbox = ttk.Spinbox(self.scrollable_frame, from_=0, to=10, width=8, font=('Arial', 12))
        self.double_bond_spinbox.set(0)
        self.double_bond_spinbox.grid(row=3, column=1, sticky=tk.W, pady=2)
        self.double_bond_spinbox.bind('<FocusOut>', self.update_valency_status)
        self.double_bond_spinbox.bind('<ButtonRelease>', self.update_valency_status)
        ttk.Label(self.scrollable_frame, text="Number of Triple Bonds:", font=('Arial', 12, 'bold')).grid(row=4, column=0, sticky=tk.W, pady=2)
        self.triple_bond_spinbox = ttk.Spinbox(self.scrollable_frame, from_=0, to=10, width=8, font=('Arial', 12))
        self.triple_bond_spinbox.set(0)
        self.triple_bond_spinbox.grid(row=4, column=1, sticky=tk.W, pady=2)
        self.triple_bond_spinbox.bind('<FocusOut>', self.update_valency_status)
        self.triple_bond_spinbox.bind('<ButtonRelease>', self.update_valency_status)
        ttk.Label(self.scrollable_frame, text="Number of Rings:", font=('Arial', 12, 'bold')).grid(row=5, column=0, sticky=tk.W, pady=2)
        self.rings_spinbox = ttk.Spinbox(self.scrollable_frame, from_=0, to=10, width=8, font=('Arial', 12))
        self.rings_spinbox.set(0)
        self.rings_spinbox.grid(row=5, column=1, sticky=tk.W, pady=2)
        self.rings_spinbox.bind('<FocusOut>', self.update_valency_status)
        self.rings_spinbox.bind('<ButtonRelease>', self.update_valency_status)
        ttk.Label(self.scrollable_frame, text="Carbon Types:", font=('Arial', 12, 'bold')).grid(row=6, column=0, sticky=tk.W, pady=2)
        carbon_frame = ttk.Frame(self.scrollable_frame)
        carbon_frame.grid(row=6, column=1, columnspan=5, sticky=tk.W, pady=2)
        ttk.Label(self.scrollable_frame, text="Functional Groups (Counts):", font=('Arial', 12, 'bold')).grid(row=7, column=0, columnspan=6, sticky=tk.W, pady=(5, 2))
        func_group_frame = ttk.Frame(self.scrollable_frame)
        func_group_frame.grid(row=8, column=0, columnspan=6, sticky=tk.W+tk.E, pady=(2, 2))
        self.func_group_spinboxes = {}
        functional_groups = [
            ('Amide', 'Amide*'),
            ('Azide', 'Azide*'),
            ('Br', 'Bromide*'),
            ('CHO', 'Aldehyde'),
            ('Cl', 'Chloride*'),
            ('CN', 'Cyanide*'),
            ('COOR_C2H5', 'Ester (C2H5)'),
            ('COOR_C3H7', 'Ester (C3H7)'),
            ('COOR_CH(CH3)2', 'Ester (CH(CH3)2)'),
            ('COOR_CH3', 'Ester (CH3)'),
            ('COOH', 'Carboxylic Acid'),
            ('COX_Br', 'Acid Bromide'),
            ('COX_Cl', 'Acid Chloride'),
            ('COX_F', 'Acid Fluoride'),
            ('COX_I', 'Acid Iodide'),
            ('Ether', 'Ether*'),
            ('F', 'Fluoride*'),
            ('I', 'Iodide*'),
            ('Imine', 'Imine*'),
            ('Ketone', 'Ketone*'),
            ('NC', 'Isocyanide*'),
            ('NCO', 'Isocyanate*'),
            ('NH2', 'Amine'),
            ('NO2', 'Nitro*'),
            ('OCN', 'Cyanate*'),
            ('OH', 'Hydroxyl*'),
            ('OX_Br', 'Hypobromite'),
            ('OX_Cl', 'Hypochlorite'),
            ('OX_F', 'Hypofluorite'),
            ('OX_I', 'Hypoiodite'),
            ('S_Bivalent', 'Sulfur (Bivalent)*'),
            ('S_Tetravalent', 'Sulfur (Tetravalent)*'),
            ('S_Hexavalent', 'Sulfur (Hexavalent)*'),
            ('S_Chain_Bi', 'S-Chain (Bi)*'),
            ('S_Chain_Tetra', 'S-Chain (Tetra)*'),
            ('S_Chain_Hexa', 'S-Chain (Hexa)*')
        ]
        num_columns = 6
        for i, (group_code, group_label) in enumerate(functional_groups):
            row = i // num_columns
            col = (i % num_columns) * 2
            label = ttk.Label(func_group_frame, text=group_label, width=17, anchor=tk.W, font=('Arial', 12))
            label.grid(row=row, column=col, sticky=tk.W, padx=(5, 2), pady=1)
            self.func_group_spinboxes[f"{group_code}_label"] = label
            spinbox = ttk.Spinbox(func_group_frame, from_=0, to=10, width=5, command=self.update_valency_status, font=('Arial', 12))
            spinbox.set(0)
            spinbox.grid(row=row, column=col+1, sticky=tk.W, padx=(2, 5), pady=1)
            spinbox.configure(validate="key", validatecommand=(spinbox.register(self.validate_spinbox), '%P', group_code))
            spinbox.bind('<FocusOut>', self.update_valency_status)
            spinbox.bind('<ButtonRelease>', self.update_valency_status)
            self.func_group_spinboxes[group_code] = spinbox
        for i in range(num_columns * 2):
            func_group_frame.columnconfigure(i, weight=1 if i % 2 == 0 else 0)
        self.valency_status_var = tk.StringVar(value="Valency status: OK")
        self.valency_status_label = ttk.Label(self.scrollable_frame, textvariable=self.valency_status_var, font=('Arial', 12))
        self.valency_status_label.grid(row=23, column=0, columnspan=6, sticky=tk.W, pady=(5, 2))
        ttk.Label(self.scrollable_frame, text="CPU Cores to use:", font=('Arial', 12, 'bold')).grid(row=24, column=0, sticky=tk.W, pady=2)
        self.cpu_cores_var = tk.StringVar(value=str(cpu_count_val))
        cpu_frame = ttk.Frame(self.scrollable_frame)
        cpu_frame.grid(row=24, column=1, columnspan=5, sticky=tk.W, pady=2)
        ttk.Radiobutton(cpu_frame, text=f"All {cpu_count_val} cores", variable=self.cpu_cores_var, value=str(cpu_count_val)).pack(side=tk.LEFT)
        ttk.Radiobutton(cpu_frame, text=f"Half ({cpu_count_val//2})", variable=self.cpu_cores_var, value=str(cpu_count_val//2)).pack(side=tk.LEFT, padx=(10, 0))
        ttk.Radiobutton(cpu_frame, text="Custom:", variable=self.cpu_cores_var, value="custom").pack(side=tk.LEFT, padx=(10, 0))
        self.custom_cores_spinbox = ttk.Spinbox(cpu_frame, from_=1, to=cpu_count_val, width=5, font=('Arial', 12))
        self.custom_cores_spinbox.set(cpu_count_val)
        self.custom_cores_spinbox.pack(side=tk.LEFT, padx=(5, 0))
        ttk.Label(self.scrollable_frame, text="Select output directory:", font=('Arial', 12, 'bold')).grid(row=25, column=0, sticky=tk.W, pady=2)
        self.outdir_var = tk.StringVar()
        self.outdir_entry = ttk.Entry(self.scrollable_frame, textvariable=self.outdir_var, width=70, font=('Arial', 12))
        self.outdir_entry.grid(row=26, column=0, columnspan=5, sticky=tk.W + tk.E, pady=2)
        self.browse_btn = ttk.Button(self.scrollable_frame, text="Browse...", command=self.browse_outdir)
        self.browse_btn.grid(row=26, column=5, sticky=tk.W, padx=5, pady=2)
        button_frame = ttk.Frame(self.scrollable_frame)
        button_frame.grid(row=27, column=0, columnspan=6, pady=5)
        self.clear_cache_btn = ttk.Button(button_frame, text="Clear Cache", command=self.clear_cache)
        self.clear_cache_btn.pack(side=tk.LEFT, padx=5)
        self.generate_btn = ttk.Button(button_frame, text="GENERATE COMPOUND STRUCTURES", command=self.start_generation)
        self.generate_btn.pack(side=tk.LEFT, padx=5)
        self.status_var = tk.StringVar(value="Ready for compound structure generation - Please select output directory")
        self.status_label = ttk.Label(self.scrollable_frame, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W, font=('Arial', 12))
        self.status_label.grid(row=28, column=0, columnspan=6, sticky=tk.E + tk.W, pady=2)
        self.progress = ttk.Progressbar(self.scrollable_frame, orient=tk.HORIZONTAL, length=800, mode='determinate')
        self.progress.grid(row=29, column=0, columnspan=6, pady=3, sticky=tk.E + tk.W)
        self.perf_label = ttk.Label(self.scrollable_frame, text="", foreground="green", font=('Arial', 12, 'bold'))
        self.perf_label.grid(row=30, column=0, columnspan=6, pady=2)
        self.batch_label = ttk.Label(self.scrollable_frame, text="", foreground="blue", font=('Arial', 12))
        self.batch_label.grid(row=31, column=0, columnspan=6, pady=2)
        info_text = ("Burry the Worries of Drawing the Molecules\n"
                 "A code by Dr. Nitin Sapre, Computational Chemist\n"
                 "Create Chemical Database with MOL files and SDF files using parallel processing.")
        self.info_label = ttk.Label(self.scrollable_frame, text=info_text, font=('Arial', 12), foreground="blue")
        self.info_label.grid(row=32, column=0, columnspan=6, pady=2)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        for i in range(6):
            self.mainframe.columnconfigure(i, weight=1)
        self.result_queue = queue.Queue()
        self.generation_thread = None
        self.update_functional_group_states()
        self.update_valency_status()

    def update_functional_group_states(self, event=None):
        try:
                carbon_count = int(self.carbon_spinbox.get())
                disabled_groups_for_zero_carbon = {
                        'Amide', 'Cl', 'Br', 'I', 'F', 'CN', 'NC', 
                        'OCN', 'NCO', 'Imine', 'NO2', 'Ketone', 'Ether', 'OH', 'NH2', 'OX_Cl', 'OX_Br', 'OX_F', 'OX_I', 'Azide',
                        'S_Bivalent', 'S_Tetravalent', 'S_Hexavalent',
                        'S_Chain_Bi', 'S_Chain_Tetra', 'S_Chain_Hexa'
                }
                always_allowed_groups = {
                        'COOH', 'CHO', 'COOR_CH3', 'COOR_C2H5', 'COOR_C3H7', 'COOR_CH(CH3)2',
                        'COX_Cl', 'COX_Br', 'COX_F', 'COX_I'
                }
                for group_code, spinbox in self.func_group_spinboxes.items():
                        if group_code.endswith('_label'):
                                continue
                        if carbon_count == 0:
                                if group_code in disabled_groups_for_zero_carbon:
                                        spinbox.configure(state='disabled')
                                        try:
                                                current_val = spinbox.get()
                                                if current_val != "" and int(current_val) > 0:
                                                        spinbox.set(0)
                                        except ValueError:
                                                spinbox.set(0)
                                else:
                                        spinbox.configure(state='normal')
                        else:
                                if group_code == 'Ether':
                                        if carbon_count >= 2:
                                                spinbox.configure(state='normal')
                                        else:
                                                spinbox.configure(state='disabled')
                                                try:
                                                        current_val = spinbox.get()
                                                        if current_val != "" and int(current_val) > 0:
                                                                spinbox.set(0)
                                                except ValueError:
                                                        spinbox.set(0)
                                elif group_code in ['S_Chain_Bi', 'S_Chain_Tetra', 'S_Chain_Hexa', 'Ketone']:
                                        if carbon_count >= 2:
                                                spinbox.configure(state='normal')
                                        else:
                                                spinbox.configure(state='disabled')
                                                try:
                                                        current_val = spinbox.get()
                                                        if current_val != "" and int(current_val) > 0:
                                                                spinbox.set(0)
                                                except ValueError:
                                                        spinbox.set(0)
                                else:
                                        spinbox.configure(state='normal')
        except ValueError:
                for group_code, spinbox in self.func_group_spinboxes.items():
                        if not group_code.endswith('_label'):
                                spinbox.configure(state='disabled')

    def browse_outdir(self):
        directory = filedialog.askdirectory(title="Select Output directory")
        if directory:
            self.outdir_var.set(directory)

    def clear_cache(self):
        get_degree_sequence_hash.cache_clear()
        _canonical_cache.clear()
        DUPLICATE_REGISTRY.clear()
        FILE_EXISTENCE_CACHE.clear()
        nomenclature_system.used_names.clear()
        self.status_var.set("Cache cleared successfully")

    def get_cpu_cores(self):
        cores_setting = self.cpu_cores_var.get()
        if cores_setting == "custom":
            return int(self.custom_cores_spinbox.get())
        else:
            return int(cores_setting)

    def get_carbon_types(self):
        return ['primary', 'secondary', 'tertiary']

    def get_functional_groups(self):
        groups = []
        for group_code, spinbox in self.func_group_spinboxes.items():
            if group_code.endswith('_label'):
                continue
            try:
                value = spinbox.get()
                count = int(value) if value != "" else 0
            except ValueError:
                count = 0
            groups.extend([group_code] * count)
        return groups

    def get_max_valency(self):
        try:
            carbon_count = int(self.carbon_spinbox.get())
            n_double_bonds = int(self.double_bond_spinbox.get())
            n_triple_bonds = int(self.triple_bond_spinbox.get())
            n_rings = int(self.rings_spinbox.get())
            if carbon_count == 0:
                return 1
            max_valency = (2 * carbon_count + 2) - (2 * n_double_bonds) - (4 * n_triple_bonds)
            return max_valency
        except:
            return 0

    def get_current_valency(self):
        total_valency = 0
        for group_code, spinbox in self.func_group_spinboxes.items():
            if group_code.endswith('_label'):
                continue
            try:
                value = spinbox.get()
                count = int(value) if value != "" else 0
            except ValueError:
                count = 0
            total_valency += count
        return total_valency

    def update_valency_status(self, event=None):
        try:
            max_valency = self.get_max_valency()
            current_valency = self.get_current_valency()
            if max_valency is None or current_valency is None:
                self.valency_status_var.set("Valency status: Calculating...")
                self.generate_btn.configure(state='disabled')
                return
            if current_valency > max_valency:
                self.valency_status_var.set(f"Valency status: EXCEEDED (Current: {current_valency}, Max: {max_valency})")
                self.valency_status_label.configure(foreground='red')
                self.generate_btn.configure(state='disabled')
                for group_code, spinbox in self.func_group_spinboxes.items():
                    if group_code.endswith('_label'):
                        continue
                    try:
                        value = spinbox.get()
                        current_val = int(value) if value != "" else 0
                    except ValueError:
                        current_val = 0
                    if current_val == 0 and spinbox['state'] == 'normal':
                        spinbox.configure(state='disabled')
                    elif current_val > 0:
                        spinbox.configure(state='normal')
            else:
                if current_valency == max_valency:
                    self.valency_status_var.set(f"Valency status: MAXIMUM REACHED (Current: {current_valency}, Max: {max_valency})")
                    self.valency_status_label.configure(foreground='orange')
                    self.generate_btn.configure(state='normal')
                    for group_code, spinbox in self.func_group_spinboxes.items():
                        if group_code.endswith('_label'):
                            continue
                        try:
                            value = spinbox.get()
                            current_val = int(value) if value != "" else 0
                        except ValueError:
                            current_val = 0
                        if current_val == 0 and spinbox['state'] == 'normal':
                            spinbox.configure(state='disabled')
                        elif current_val > 0:
                            spinbox.configure(state='normal')
                else:
                    self.valency_status_var.set(f"Valency status: OK (Current: {current_valency}, Max: {max_valency})")
                    self.valency_status_label.configure(foreground='green')
                    self.generate_btn.configure(state='normal')
                    self.update_functional_group_states()
        except Exception as e:
            self.valency_status_var.set(f"Valency status: Error calculating ({str(e)})")
            self.valency_status_label.configure(foreground='red')
            self.generate_btn.configure(state='disabled')

    def validate_spinbox(self, new_value, group_code):
        if not new_value.isdigit() and new_value != "":
            return False
        if new_value == "":
            return True
        new_count = int(new_value)
        try:
            current_value = self.func_group_spinboxes[group_code].get()
            current_count = int(current_value) if current_value != "" else 0
        except ValueError:
            current_count = 0
        if new_count < current_count:
            return True
        max_valency = self.get_max_valency()
        current_valency = self.get_current_valency()
        if current_valency + (new_count - current_count) > max_valency:
            messagebox.showwarning("Maximum Valency Reached", 
                                                     f"Cannot add more functional groups. The maximum allowed for these parameters is {max_valency}.")
            return False
        return True

    def calculate_unsaturation(self, carbon_count, n_double_bonds, n_triple_bonds, n_rings):
        total_unsaturation = n_double_bonds + 2 * n_triple_bonds + n_rings
        unsat_type = []
        if n_double_bonds > 0:
            unsat_type.append(f"DB{n_double_bonds}")
        if n_triple_bonds > 0:
            unsat_type.append(f"TB{n_triple_bonds}")
        if n_rings > 0:
            unsat_type.append(f"R{n_rings}")
        return total_unsaturation, "_".join(unsat_type) if unsat_type else "SAT"

    def start_generation(self):
        if self.generation_thread and self.generation_thread.is_alive():
            messagebox.showwarning("Generation in Progress", "Please wait for the current generation to complete.")
            return
        max_valency = self.get_max_valency()
        current_valency = self.get_current_valency()
        if current_valency > max_valency:
            messagebox.showerror("Valency Exceeded", 
                               f"Cannot generate structures: Current valency ({current_valency}) exceeds maximum valency ({max_valency})")
            return
        try:
            carbon_count = int(self.carbon_spinbox.get())
            if carbon_count < 0 or carbon_count > 20:
                raise ValueError("Carbon count must be between 0 and 20")
            n_double_bonds = int(self.double_bond_spinbox.get())
            if n_double_bonds < 0:
                raise ValueError("Double bond count cannot be negative")
            n_triple_bonds = int(self.triple_bond_spinbox.get())
            if n_triple_bonds < 0:
                raise ValueError("Triple bond count cannot be negative")
            n_rings = int(self.rings_spinbox.get())
            if n_rings < 0:
                raise ValueError("Number of rings cannot be negative")
            carbon_types = self.get_carbon_types()
            functional_groups = self.get_functional_groups()
            is_valid, error_message = validate_structure_possibility(
                carbon_count, functional_groups, n_double_bonds, n_triple_bonds, carbon_types, n_rings
            )
            if not is_valid:
                messagebox.showerror("Structure Generation Impossible", 
                                   f"Cannot generate any structures with these parameters:\n\n{error_message}\n\n" +
                                   "Please modify your parameters and try again.")
                return
        except ValueError as e:
            messagebox.showerror("Invalid Input", f"Invalid parameters: {e}")
            return
        output_dir = self.outdir_var.get().strip()
        if not output_dir:
            messagebox.showerror("Invalid Input", "Please select an output directory")
            return
        num_cores = self.get_cpu_cores()
        self.generate_btn.configure(state='disabled')
        self.generation_thread = threading.Thread(
            target=self.generate_structures_thread,
            args=(carbon_count, functional_groups, n_double_bonds, n_triple_bonds, n_rings, carbon_types, output_dir, num_cores),
            daemon=True
        )
        self.generation_thread.start()
        self.monitor_generation()

    def parallel_progress_callback(self, progress, processed_batches, total_batches):
        self.result_queue.put(("batch_progress", f"Parallel processing: {processed_batches}/{total_batches} batches complete"))
        overall_progress = 30 + int(progress * 0.5)
        self.result_queue.put(("progress", overall_progress))

    def generate_structures_thread(self, carbon_count, functional_groups, n_double_bonds, n_triple_bonds, n_rings, carbon_types, output_dir, num_cores):
        try:
            start_time = time.time()
            self.result_queue.put(("status", "Validating structure generation parameters..."))
            is_valid, error_message = validate_structure_possibility(
                carbon_count, functional_groups, n_double_bonds, n_triple_bonds, carbon_types, n_rings
            )
            if not is_valid:
                self.result_queue.put(("no_structures", error_message))
                self.result_queue.put(("status", "VALIDATION ERROR"))
                self.result_queue.put(("progress", 0))
                return
            total_unsat, unsat_type_str = self.calculate_unsaturation(
                carbon_count, n_double_bonds, n_triple_bonds, n_rings
            )
            db_str = f"_DB{n_double_bonds}" if n_double_bonds > 0 else ""
            tb_str = f"_TB{n_triple_bonds}" if n_triple_bonds > 0 else ""
            ring_str = f"_R{n_rings}" if n_rings > 0 else ""
            type_str = ""
            if carbon_count > 0:
                if 'primary' in carbon_types:
                    type_str += "P"
                if 'secondary' in carbon_types:
                    type_str += "S"
                if 'tertiary' in carbon_types:
                    type_str += "T"
            else:
                type_str = "NoC"
            func_counts = {}
            for fg in functional_groups:
                func_counts[fg] = func_counts.get(fg, 0) + 1
            func_str_parts = []
            for fg, count in func_counts.items():
                readable_name = get_functional_group_name(fg)
                if count == 1:
                    func_str_parts.append(readable_name)
                else:
                    func_str_parts.append(f"{readable_name}{count}")
            func_str = "_".join(func_str_parts)
            if len(func_str) > 20:
                func_str = func_str[:20] + "..."
            base_name = f"C{carbon_count}{db_str}{tb_str}{ring_str}_{type_str}_{func_str}_U{total_unsat}_{unsat_type_str}"
            main_output_dir = os.path.join(output_dir, base_name)
            self.result_queue.put(("status", f"Starting compound structure generation using {num_cores} CPU cores..."))
            self.result_queue.put(("progress", 5))
            smiles_start = time.time()
            raw_smiles_list = generate_functionalized_isomers(carbon_count, functional_groups, n_double_bonds, n_triple_bonds, n_rings, carbon_types)
            if not raw_smiles_list:
                error_msg = f"NO STRUCTURES GENERATED!\n\n"
                error_msg += f"Parameters used:\n"
                error_msg += f"• Carbons: {carbon_count}\n"
                error_msg += f"• Double bonds: {n_double_bonds}\n" 
                error_msg += f"• Triple bonds: {n_triple_bonds}\n"
                error_msg += f"• Rings: {n_rings}\n"
                error_msg += f"• Functional groups: {', '.join(functional_groups)}\n"
                error_msg += f"• Carbon types: {', '.join(carbon_types)}\n\n"
                error_msg += "This combination is chemically impossible or creates invalid structures.\n"
                error_msg += "Please try different parameters."
                self.result_queue.put(("no_structures", error_msg))
                self.result_queue.put(("status", "NO STRUCTURES GENERATED"))
                self.result_queue.put(("progress", 0))
                return
            unique_map = OrderedDict()
            for s in raw_smiles_list:
                can = canonicalize_smiles(s)
                if can:
                    unique_map.setdefault(can, None)
            smiles_list = list(unique_map.keys())
            smiles_gen_time = time.time() - smiles_start
            os.makedirs(main_output_dir, exist_ok=True)
            self.result_queue.put(("status", f"Generated {len(smiles_list)} unique compound isomers in {smiles_gen_time:.2f}s"))
            self.result_queue.put(("progress", 50))
            tasks = []
            for idx, smiles in enumerate(smiles_list, start=1):
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    continue
                formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                compound_id = f"{formula}_{idx:04d}"
                tasks.append((idx, smiles, compound_id))
                del mol  # Clear memory
            proc_start = time.time()
            self.result_queue.put(("status", f"Starting parallel structure processing with {num_cores} processes..."))
            all_results = parallel_alkane_generation(
                tasks,
                base_name,
                main_output_dir,
                functional_groups,
                carbon_count,
                n_double_bonds,
                n_triple_bonds,
                n_rings,
                progress_callback=self.parallel_progress_callback,
                num_processes=num_cores
            )
            proc_gen_time = time.time() - proc_start
            self.result_queue.put(("status", f"Processed {len(all_results)} structures in {proc_gen_time:.2f}s using {num_cores} cores"))
            self.result_queue.put(("progress", 80))
            self.result_queue.put(("status", "Generating CSV file..."))
            csv_path = os.path.join(main_output_dir, f"{base_name}.csv")
            header = ["No.", "Molecular Formula", "SMILES", "Functional Group", "Ring Size", "Unique Compound ID"]
            with open(csv_path, mode='w', newline='', encoding='utf-8') as csv_file:
                csv_writer = csv.writer(csv_file)
                csv_writer.writerow(header)
                valid_results = sorted([r for r in all_results if r is not None], key=lambda x: x['index'])
                for result in valid_results:
                    row = [
                        result['index'],
                        result['formula'],
                        result.get('smiles', 'N/A'),
                        result['functional_group_type'],
                        result.get('ring_size', 'ACYCLIC'),
                        result['compound_id']
                    ]
                    csv_writer.writerow(row)
            self.result_queue.put(("progress", 85))
            total_time = time.time() - start_time
            summary_filename = f"{base_name}_Summary.txt"
            summary_path = os.path.join(main_output_dir, summary_filename)
            with open(summary_path, 'w', encoding='utf-8') as f:
                f.write(f"COMPLETE PARALLEL COMPOUND GENERATION SUMMARY - C{carbon_count}\n")
                f.write("=" * 80 + "\n\n")
                f.write(f"GENERATE BUTTON SAFETY: Button disabled when valency exceeded\n")
                f.write(f"- Prevents generation of impossible chemical structures\n")
                f.write(f"- Automatic valency calculation and validation\n")
                f.write(f"- User-friendly visual feedback\n\n")
                f.write(f"PARAMETERS:\n")
                f.write("-" * 15 + "\n")
                f.write(f"Carbon atoms: {carbon_count}\n")
                f.write(f"Double bonds: {n_double_bonds}\n")
                f.write(f"Triple bonds: {n_triple_bonds}\n")
                f.write(f"Rings: {n_rings}\n")
                if carbon_count > 0:
                    f.write(f"Carbon types: {', '.join(carbon_types)}\n")
                f.write(f"Functional groups: {', '.join(functional_groups)}\n\n")
                f.write(f"PERFORMANCE METRICS:\n")
                f.write("-" * 30 + "\n")
                f.write(f"CPU Cores Used: {num_cores} / {cpu_count()}\n")
                f.write(f"Total Structures Generated: {len(valid_results)}\n")
                f.write(f"SMILES Generation Time: {smiles_gen_time:.3f} seconds\n")
                f.write(f"Parallel Processing Time: {proc_gen_time:.3f} seconds\n")
                f.write(f"Total Generation Time: {total_time:.3f} seconds\n")
                if total_time > 0:
                    f.write(f"Generation Rate: {len(valid_results)/total_time:.1f} structures/sec\n")
                else:
                    f.write("Generation Rate: N/A\n")
                if proc_gen_time > 0:
                    f.write(f"Parallel Efficiency: {len(valid_results)/proc_gen_time:.1f} structures/sec\n")
                    f.write(f"Core Utilization: {(len(valid_results)/proc_gen_time)/num_cores:.1f} structures/sec/core\n\n")
                else:
                    f.write("Parallel Efficiency: N/A\n\n")
                f.write(f"OUTPUT FILES GENERATED:\n")
                f.write("-" * 25 + "\n")
                f.write(f"- Data table: {base_name}.csv (with SMILES)\n")
                f.write(f"- MOL files: {len(valid_results)} files in element-based directories\n")
                f.write(f"- SDF files: {len(valid_results)} files in element-based directories\n")
                f.write(f"- Performance summary: {summary_filename}\n")
                f.write(f"\nAll files are organized in: {base_name}/\n")
                f.write(f"Files are also copied to: Main_Database/\n")
                # Add ring size distribution information
                if n_rings > 0:
                    ring_size_distribution = defaultdict(int)
                    for result in valid_results:
                        ring_size = result.get('ring_size', 'ACYCLIC')
                        ring_size_distribution[ring_size] += 1
                    f.write(f"\nRING SIZE DISTRIBUTION:\n")
                    f.write("-" * 25 + "\n")
                    for ring_size, count in sorted(ring_size_distribution.items()):
                        f.write(f"- {ring_size}: {count} structures\n")
            self.result_queue.put(("progress", 100))
            self.result_queue.put(("status", f"COMPOUND GENERATION FINISHED: {len(valid_results)} structures in {total_time:.2f}s using {num_cores} cores"))
            throughput = len(valid_results) / total_time if total_time > 0 else 0.0
            parallel_efficiency = (len(valid_results) / proc_gen_time) / num_cores if proc_gen_time > 0 else 0.0
            self.result_queue.put(("performance", f"Throughput: {throughput:.2f} structures/sec | Efficiency: {parallel_efficiency:.2f} struct/sec/core"))
        except Exception as e:
            self.result_queue.put(("status", f"Error: {e}"))
            self.result_queue.put(("progress", 0))

    def monitor_generation(self):
        try:
            while not self.result_queue.empty():
                msg_type, msg_value = self.result_queue.get_nowait()
                if msg_type == "status":
                    self.status_var.set(msg_value)
                elif msg_type == "progress":
                    self.progress['value'] = msg_value
                elif msg_type == "performance":
                    self.perf_label.config(text=msg_value)
                elif msg_type == "batch_progress":
                    self.batch_label.config(text=msg_value)
                elif msg_type == "no_structures":
                    messagebox.showerror("NO STRUCTURES POSSIBLE", msg_value)
        except Exception:
            pass
        if self.generation_thread and self.generation_thread.is_alive():
            self.root.after(200, self.monitor_generation)
        else:
            self.status_var.set("Finished. Ready for next generation.")
            self.update_valency_status()

def main():
    root = tk.Tk()
    app = UltraFastCompoundGeneratorGUI(root)
    root.mainloop()

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
