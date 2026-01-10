import re
import csv
import threading
import queue
import multiprocessing
from multiprocessing import Pool, cpu_count
from functools import lru_cache, partial
import time
import random
from typing import List, Set, Tuple, Dict
import itertools
from concurrent.futures import ProcessPoolExecutor, as_completed
import sys
import os
import shutil
from collections import OrderedDict, defaultdict
import networkx as nx
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw, AllChem, SDWriter, rdMolDescriptors, Descriptors
from rdkit.Chem.rdchem import BondType
import psutil
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
                    return False, "IMPOSSIBLE: Ketones cannot be formed with only 1 carbon atom"
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

    # If multiple ring sizes, use the smallest one for classification
    smallest_ring = min(size_count.keys())
    return f"R{smallest_ring}"

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
