'use strict';

const OCL = require('openchemlib');

const canonicalCache = new Map();
const rootedTreeCache = new Map();
const partitionsCache = new Map();
const alkaneBackboneCache = new Map();
const functionalGroupQueryCache = new Map();

const DEFAULT_CARBON_TYPES = ['primary', 'secondary', 'tertiary'];
const MAX_ASSIGNMENT_ATTEMPTS = 50000;
const MAX_RING_COMBINATIONS = 7500;
const WEISFEILER_ITERATIONS = 4;

const FUNCTIONAL_GROUP_PATTERNS = {
  Alcohol: 'CO',
  Carboxylic_Acid: 'C(=O)O',
  Aldehyde: 'C=O',
  Ketone: 'C(=O)C',
  Ester: 'C(=O)OC',
  Nitrile: 'C#N',
  Isocyanide: 'N#C',
  Cyanate: 'OC#N',
  Isocyanate: 'N=C=O',
  Fluoride: 'CF',
  Chloride: 'CCl',
  Bromide: 'CBr',
  Iodide: 'CI',
  Nitro: 'N(=O)O',
  Ether: 'COC',
  Imine: 'C=N',
  Amide: 'C(=O)N',
  Acid_Halide: 'C(=O)Cl',
  Hypohalite: 'OO',
  Amine: 'CN',
  Azide: 'N=N=N',
  S_Bivalent: 'CS',
  S_Tetravalent: 'CS(=O)',
  S_Hexavalent: 'CS(=O)(=O)',
  S_Chain_Bi: 'CSC',
  S_Chain_Tetra: 'CS(=O)C',
  S_Chain_Hexa: 'CS(=O)(=O)C'
};

const SULFUR_CHAIN_GROUPS = new Set(['S_Chain_Bi', 'S_Chain_Tetra', 'S_Chain_Hexa']);

rootedTreeCache.set(1, [{ size: 1, children: [], canonical: '()' }]);

function canonicalizeSmiles(smiles) {
  if (!smiles) {
    return '';
  }
  if (canonicalCache.has(smiles)) {
    return canonicalCache.get(smiles);
  }
  let canonical = '';
  try {
    const mol = OCL.Molecule.fromSmiles(smiles);
    const bondCount = mol.getAllBonds?.() ?? 0;
    for (let bond = 0; bond < bondCount; bond++) {
      if (typeof mol.setBondParityUnknownOrNone === 'function') {
        mol.setBondParityUnknownOrNone(bond);
      }
    }
    canonical = typeof mol.toSmiles === 'function' ? mol.toSmiles() : mol.toIsomericSmiles();
    canonical = canonical.replace(/[\\/]/g, '');
  } catch (err) {
    canonical = '';
  }
  canonicalCache.set(smiles, canonical);
  return canonical;
}

function createRootedTree(children) {
  const orderedChildren = children.slice().sort((a, b) => a.canonical.localeCompare(b.canonical));
  const canonical = `(${orderedChildren.map((child) => child.canonical).join('')})`;
  const size = 1 + orderedChildren.reduce((sum, child) => sum + child.size, 0);
  return { size, children: orderedChildren, canonical };
}

function getIntegerPartitions(n, minPart = 1) {
  const key = `${n}|${minPart}`;
  if (partitionsCache.has(key)) {
    return partitionsCache.get(key);
  }
  const partitions = [];
  if (n === 0) {
    partitions.push([]);
  } else {
    for (let i = minPart; i <= n; i++) {
      const remainder = getIntegerPartitions(n - i, i);
      for (const tail of remainder) {
        partitions.push([i, ...tail]);
      }
    }
  }
  partitionsCache.set(key, partitions);
  return partitions;
}

function* combinations(items, k, start = 0, prefix = []) {
  if (k === 0) {
    yield prefix.slice();
    return;
  }
  for (let i = start; i <= items.length - k; i++) {
    prefix.push(items[i]);
    yield* combinations(items, k - 1, i + 1, prefix);
    prefix.pop();
  }
}

function* combinationsWithReplacement(items, k, start = 0, prefix = []) {
  if (k === 0) {
    yield prefix.slice();
    return;
  }
  for (let i = start; i < items.length; i++) {
    prefix.push(items[i]);
    yield* combinationsWithReplacement(items, k - 1, i, prefix);
    prefix.pop();
  }
}

function groupPartition(parts) {
  const counters = new Map();
  for (const size of parts) {
    counters.set(size, (counters.get(size) || 0) + 1);
  }
  return Array.from(counters.entries()).map(([size, count]) => ({ size, count })).sort((a, b) => a.size - b.size);
}

function* enumerateChildAssignments(totalNodes) {
  for (const partition of getIntegerPartitions(totalNodes)) {
    const grouped = groupPartition(partition);
    yield* enumerateGroupedAssignments(grouped, 0, []);
  }
}

function* enumerateGroupedAssignments(groups, index, buffer) {
  if (index >= groups.length) {
    yield buffer.slice();
    return;
  }
  const { size, count } = groups[index];
  const candidates = getRootedTrees(size);
  for (const combo of combinationsWithReplacement(candidates, count)) {
    buffer.push(...combo);
    yield* enumerateGroupedAssignments(groups, index + 1, buffer);
    buffer.length -= combo.length;
  }
}

function getRootedTrees(size) {
  if (rootedTreeCache.has(size)) {
    return rootedTreeCache.get(size);
  }
  const trees = new Map();
  for (const children of enumerateChildAssignments(size - 1)) {
    const tree = createRootedTree(children);
    trees.set(tree.canonical, tree);
  }
  const list = Array.from(trees.values()).sort((a, b) => a.canonical.localeCompare(b.canonical));
  rootedTreeCache.set(size, list);
  return list;
}

function normalizeEdge(a, b) {
  return a < b ? [a, b] : [b, a];
}

function cloneRootedSubtree(tree, startIndex) {
  const current = startIndex;
  let cursor = current + 1;
  const edges = [];
  for (const child of tree.children) {
    const cloned = cloneRootedSubtree(child, cursor);
    edges.push(normalizeEdge(current, cloned.root));
    edges.push(...cloned.edges);
    cursor = cloned.nextIndex;
  }
  return { root: current, edges, nextIndex: cursor };
}

function buildTreeEdges(children, totalNodes) {
  if (totalNodes === 1) {
    return [];
  }
  let cursor = 1;
  const edges = [];
  for (const child of children) {
    const cloned = cloneRootedSubtree(child, cursor);
    edges.push(normalizeEdge(0, cloned.root));
    edges.push(...cloned.edges);
    cursor = cloned.nextIndex;
  }
  return edges;
}

function edgesToAdjacency(edges, nodeCount) {
  const adjacency = Array.from({ length: nodeCount }, () => []);
  for (const [a, b] of edges) {
    adjacency[a].push(b);
    adjacency[b].push(a);
  }
  return adjacency;
}

function findTreeCenters(adjacency) {
  const remaining = adjacency.length;
  if (remaining <= 2) {
    return Array.from({ length: remaining }, (_, idx) => idx);
  }
  const degrees = adjacency.map((neighbors) => neighbors.length);
  let leaves = degrees.map((deg, idx) => (deg <= 1 ? idx : null)).filter((idx) => idx !== null);
  let nodesLeft = remaining;
  while (nodesLeft > 2) {
    nodesLeft -= leaves.length;
    const nextLeaves = [];
    for (const leaf of leaves) {
      degrees[leaf] = 0;
      for (const neighbor of adjacency[leaf]) {
        if (degrees[neighbor] === 0) continue;
        degrees[neighbor] -= 1;
        if (degrees[neighbor] === 1) {
          nextLeaves.push(neighbor);
        }
      }
    }
    leaves = nextLeaves;
  }
  return leaves;
}

function encodeRooted(adjacency, node, parent) {
  const labels = [];
  for (const neighbor of adjacency[node]) {
    if (neighbor === parent) continue;
    labels.push(encodeRooted(adjacency, neighbor, node));
  }
  labels.sort();
  return `(${labels.join('')})`;
}

function canonicalTreeLabel(adjacency) {
  if (adjacency.length === 1) {
    return '()';
  }
  const centers = findTreeCenters(adjacency);
  if (centers.length === 1) {
    return encodeRooted(adjacency, centers[0], -1);
  }
  const encodings = [
    encodeRooted(adjacency, centers[0], centers[1]),
    encodeRooted(adjacency, centers[1], centers[0])
  ].sort();
  return `[${encodings.join('')}]`;
}

function generateTreeEdgeSets(size) {
  if (size === 0) {
    return [[]];
  }
  if (size === 1) {
    return [[]];
  }
  const unique = new Map();
  for (const children of enumerateChildAssignments(size - 1)) {
    const edges = buildTreeEdges(children, size);
    const adjacency = edgesToAdjacency(edges, size);
    const label = canonicalTreeLabel(adjacency);
    if (!unique.has(label)) {
      unique.set(label, edges);
    }
  }
  return Array.from(unique.values());
}

function createCarbonSkeleton(nAtoms, nBonds) {
  const mol = new OCL.Molecule(Math.max(nAtoms + 8, 64), Math.max(nBonds + 8, 64));
  mol.setFragment(false);
  for (let i = 0; i < nAtoms; i++) {
    mol.addAtom(6);
  }
  return mol;
}

function buildSmilesFromEdges(nAtoms, edges, orders = null) {
  if (nAtoms === 0) {
    return '';
  }
  try {
    const mol = createCarbonSkeleton(nAtoms, edges.length);
    edges.forEach(([a, b], idx) => {
      const order = orders ? orders[idx] : 1;
      const bondIdx = mol.addBond(a, b);
      mol.setBondOrder(bondIdx, order);
    });
    mol.ensureHelperArrays(OCL.Molecule.cHelperNeighbours);
    return mol.toIsomericSmiles();
  } catch (err) {
    return '';
  }
}

function generateAlkaneBackbone(nCarbons) {
  if (nCarbons === 0) {
    return [''];
  }
  if (alkaneBackboneCache.has(nCarbons)) {
    return alkaneBackboneCache.get(nCarbons).slice();
  }
  const edgeSets = generateTreeEdgeSets(nCarbons);
  const smilesSet = new Set();
  for (const edges of edgeSets) {
    const smiles = buildSmilesFromEdges(nCarbons, edges);
    if (smiles) {
      smilesSet.add(canonicalizeSmiles(smiles));
    }
  }
  const result = Array.from(smilesSet).sort();
  alkaneBackboneCache.set(nCarbons, result);
  return result.slice();
}

function enumerateNonEdges(nAtoms, edges) {
  const existing = new Set(edges.map(([a, b]) => `${Math.min(a, b)}-${Math.max(a, b)}`));
  const candidates = [];
  for (let i = 0; i < nAtoms; i++) {
    for (let j = i + 1; j < nAtoms; j++) {
      const key = `${i}-${j}`;
      if (!existing.has(key)) {
        candidates.push([i, j]);
      }
    }
  }
  return candidates;
}

function generateCyclicGraphs(nCarbons, nRings) {
  if (nCarbons < 3 || nRings <= 0) {
    return [];
  }
  const trees = generateTreeEdgeSets(nCarbons);
  const graphs = [];
  for (const treeEdges of trees) {
    const candidates = enumerateNonEdges(nCarbons, treeEdges);
    if (candidates.length < nRings) {
      continue;
    }
    let attempts = 0;
    for (const combo of combinations(candidates, nRings)) {
      attempts += 1;
      if (attempts > MAX_RING_COMBINATIONS) {
        break;
      }
      const augmented = treeEdges.slice();
      combo.forEach(([a, b]) => augmented.push(normalizeEdge(a, b)));
      graphs.push(augmented);
    }
  }
  return graphs;
}

function addHydrocarbonFromEdges(edges, nAtoms, nDouble, nTriple, smilesSet) {
  const edgeCount = edges.length;
  const totalUnsaturations = nDouble + nTriple;
  if (totalUnsaturations > edgeCount) {
    return;
  }
  if (totalUnsaturations === 0) {
    const smiles = buildSmilesFromEdges(nAtoms, edges);
    if (smiles) {
      smilesSet.add(canonicalizeSmiles(smiles));
    }
    return;
  }
  const isValenceValid = (orders) => {
    const valence = new Array(nAtoms).fill(0);
    edges.forEach(([a, b], idx) => {
      const order = orders ? orders[idx] : 1;
      valence[a] += order;
      valence[b] += order;
    });
    return valence.every((value) => value <= 4);
  };
  const indices = edges.map((_, idx) => idx);
  for (const doubleSelection of combinations(indices, nDouble)) {
    const doubleSet = new Set(doubleSelection);
    const remaining = indices.filter((idx) => !doubleSet.has(idx));
    for (const tripleSelection of combinations(remaining, nTriple)) {
      const orders = new Array(edgeCount).fill(1);
      for (const idx of doubleSelection) {
        orders[idx] = 2;
      }
      for (const idx of tripleSelection) {
        orders[idx] = 3;
      }
      if (!isValenceValid(orders)) {
        continue;
      }
      const smiles = buildSmilesFromEdges(nAtoms, edges, orders);
      if (smiles) {
        smilesSet.add(canonicalizeSmiles(smiles));
      }
    }
  }
}

function addExactUnsaturation(backboneSmiles, nDoubleBonds, nTripleBonds) {
  if (nDoubleBonds === 0 && nTripleBonds === 0) {
    return [backboneSmiles];
  }
  let mol;
  try {
    mol = OCL.Molecule.fromSmiles(backboneSmiles);
  } catch (err) {
    return [];
  }
  const totalAtoms = mol.getAllAtoms();
  const edges = [];
  const bondCount = mol.getAllBonds();
  for (let bond = 0; bond < bondCount; bond++) {
    const begin = mol.getBondAtom(0, bond);
    const end = mol.getBondAtom(1, bond);
    if (mol.getAtomicNo(begin) === 6 && mol.getAtomicNo(end) === 6) {
      edges.push([Math.min(begin, end), Math.max(begin, end)]);
    }
  }
  const smilesSet = new Set();
  addHydrocarbonFromEdges(edges, totalAtoms, nDoubleBonds, nTripleBonds, smilesSet);
  return Array.from(smilesSet).sort();
}

function generateHydrocarbonIsomers(nCarbons, nDoubleBonds = 0, nTripleBonds = 0, nRings = 0) {
  if (nCarbons < 0) {
    return [];
  }
  if (nCarbons === 0) {
    return nDoubleBonds === 0 && nTripleBonds === 0 && nRings === 0 ? [''] : [];
  }
  if (nCarbons === 1) {
    if (nDoubleBonds === 1 && nTripleBonds === 0 && nRings === 0) {
      return ['C=C'];
    }
    if (nTripleBonds === 1 && nDoubleBonds === 0 && nRings === 0) {
      return ['C#C'];
    }
    if (nDoubleBonds === 0 && nTripleBonds === 0 && nRings === 0) {
      return ['C'];
    }
    return [];
  }
  const smilesSet = new Set();
  if (nRings === 0) {
    const trees = generateTreeEdgeSets(nCarbons);
    for (const edges of trees) {
      addHydrocarbonFromEdges(edges, nCarbons, nDoubleBonds, nTripleBonds, smilesSet);
    }
  } else {
    const graphs = generateCyclicGraphs(nCarbons, nRings);
    for (const edges of graphs) {
      addHydrocarbonFromEdges(edges, nCarbons, nDoubleBonds, nTripleBonds, smilesSet);
    }
  }
  return Array.from(smilesSet).sort();
}

function ensureHelperArrays(mol, mask = OCL.Molecule.cHelperNeighbours) {
  const status = typeof mol.getHelperArrayStatus === 'function' ? mol.getHelperArrayStatus() : 0;
  if ((status & mask) === mask) {
    return;
  }
  mol.ensureHelperArrays(mask);
}

function parseMolGraph(mol) {
  ensureHelperArrays(mol, OCL.Molecule.cHelperNeighbours);
  const atomCount = mol.getAllAtoms();
  const atoms = [];
  const adjacency = Array.from({ length: atomCount }, () => []);
  for (let atom = 0; atom < atomCount; atom++) {
    atoms.push({ idx: atom, label: mol.getAtomLabel(atom) });
    const conn = mol.getConnAtoms(atom);
    for (let i = 0; i < conn; i++) {
      const neighbor = mol.getConnAtom(atom, i);
      const bond = mol.getConnBond(atom, i);
      const order = mol.getBondOrder(bond);
      adjacency[atom].push({ idx: neighbor, order });
    }
  }
  return { atoms, adjacency };
}

function computeAtomValence(graph, idx) {
  return graph.adjacency[idx].reduce((sum, neighbor) => sum + (neighbor.order || 1), 0);
}

function classifyCarbonAtoms(mol, allowedTypes) {
  const graph = parseMolGraph(mol);
  const carbons = [];
  graph.atoms.forEach((atom, idx) => {
    if (atom.label !== 'C') {
      return;
    }
    const carbonNeighbors = graph.adjacency[idx].filter((neighbor) => graph.atoms[neighbor.idx]?.label === 'C').length;
    let type = 'primary';
    if (carbonNeighbors === 2) {
      type = 'secondary';
    } else if (carbonNeighbors === 3) {
      type = 'tertiary';
    }
    if (allowedTypes.includes(type)) {
      carbons.push({ idx, type, valence: computeAtomValence(graph, idx) });
    }
  });
  return carbons;
}

function cloneEditableMol(mol) {
  const smiles = mol.toIsomericSmiles();
  return OCL.Molecule.fromSmiles(smiles);
}

function addSimpleAtom(mol, symbol) {
  const atomicNo = OCL.Molecule.getAtomicNoFromLabel(symbol);
  return mol.addAtom(atomicNo);
}

function addBondWithOrder(mol, begin, end, order = 1) {
  const bondIdx = mol.addBond(begin, end);
  mol.setBondOrder(bondIdx, order);
  return bondIdx;
}

function addExplicitHydrogen(mol, attachIdx) {
  const hIdx = addSimpleAtom(mol, 'H');
  addBondWithOrder(mol, attachIdx, hIdx, 1);
  return hIdx;
}

function setAtomCharge(mol, atomIdx, charge) {
  if (typeof mol.setAtomCharge === 'function') {
    mol.setAtomCharge(atomIdx, charge);
  }
}

function attachLinearChain(mol, startIdx, length) {
  let previous = startIdx;
  for (let i = 0; i < length; i++) {
    const carbon = addSimpleAtom(mol, 'C');
    addBondWithOrder(mol, previous, carbon, 1);
    previous = carbon;
  }
}

function attachAlkylGroup(mol, startIdx, descriptor) {
  switch (descriptor) {
    case 'CH3':
      attachLinearChain(mol, startIdx, 1);
      break;
    case 'C2H5':
      attachLinearChain(mol, startIdx, 2);
      break;
    case 'C3H7':
      attachLinearChain(mol, startIdx, 3);
      break;
    case 'CH(CH3)2': {
      const central = addSimpleAtom(mol, 'C');
      addBondWithOrder(mol, startIdx, central, 1);
      const branch1 = addSimpleAtom(mol, 'C');
      const branch2 = addSimpleAtom(mol, 'C');
      addBondWithOrder(mol, central, branch1, 1);
      addBondWithOrder(mol, central, branch2, 1);
      break;
    }
    default:
      attachLinearChain(mol, startIdx, 1);
      break;
  }
}

function attachEsterGroup(mol, carbonIdx, descriptor) {
  const carbonyl = addSimpleAtom(mol, 'C');
  const oDouble = addSimpleAtom(mol, 'O');
  const oSingle = addSimpleAtom(mol, 'O');
  addBondWithOrder(mol, carbonIdx, carbonyl, 1);
  addBondWithOrder(mol, carbonyl, oDouble, 2);
  addBondWithOrder(mol, carbonyl, oSingle, 1);
  attachAlkylGroup(mol, oSingle, descriptor);
}

function attachAcidHalideGroup(mol, carbonIdx, halogen) {
  const carbonyl = addSimpleAtom(mol, 'C');
  const oIdx = addSimpleAtom(mol, 'O');
  const halIdx = addSimpleAtom(mol, halogen);
  addBondWithOrder(mol, carbonIdx, carbonyl, 1);
  addBondWithOrder(mol, carbonyl, oIdx, 2);
  addBondWithOrder(mol, carbonyl, halIdx, 1);
}

function attachHypohaliteGroup(mol, carbonIdx, halogen) {
  const oIdx = addSimpleAtom(mol, 'O');
  const halIdx = addSimpleAtom(mol, halogen);
  addBondWithOrder(mol, carbonIdx, oIdx, 1);
  addBondWithOrder(mol, oIdx, halIdx, 1);
}

function insertSulfurChain(mol, beginAtom, endAtom, type) {
  const sIdx = addSimpleAtom(mol, 'S');
  addBondWithOrder(mol, beginAtom, sIdx, 1);
  addBondWithOrder(mol, sIdx, endAtom, 1);
  if (type === 'S_Chain_Tetra') {
    const oIdx = addSimpleAtom(mol, 'O');
    addBondWithOrder(mol, sIdx, oIdx, 2);
  } else if (type === 'S_Chain_Hexa') {
    const o1Idx = addSimpleAtom(mol, 'O');
    const o2Idx = addSimpleAtom(mol, 'O');
    addBondWithOrder(mol, sIdx, o1Idx, 2);
    addBondWithOrder(mol, sIdx, o2Idx, 2);
  }
}

function countUnsaturations(mol) {
  let doubles = 0;
  let triples = 0;
  const bondCount = mol.getAllBonds();
  for (let bond = 0; bond < bondCount; bond++) {
    const order = mol.getBondOrder(bond);
    if (order === 2) {
      doubles += 1;
    } else if (order === 3) {
      triples += 1;
    }
  }
  return { doubles, triples };
}

function* enumerateChainAssignments(chainTypes, counts, slots, buffer = []) {
  if (buffer.length === slots) {
    yield buffer.slice();
    return;
  }
  for (const type of chainTypes) {
    if ((counts[type] || 0) === 0) {
      continue;
    }
    counts[type] -= 1;
    buffer.push(type);
    yield* enumerateChainAssignments(chainTypes, counts, slots, buffer);
    buffer.pop();
    counts[type] += 1;
  }
}

function generateSulfurChainStructures(backboneSmiles, chainCounts, nDoubleBonds, nTripleBonds) {
  const totalChains = Object.values(chainCounts).reduce((sum, value) => sum + value, 0);
  if (totalChains === 0) {
    return [backboneSmiles];
  }
  let base;
  try {
    base = OCL.Molecule.fromSmiles(backboneSmiles);
  } catch (err) {
    return [];
  }
  const bondCount = base.getAllBonds();
  const ccSingleBonds = [];
  for (let bond = 0; bond < bondCount; bond++) {
    if (base.getBondOrder(bond) !== 1) {
      continue;
    }
    const begin = base.getBondAtom(0, bond);
    const end = base.getBondAtom(1, bond);
    if (base.getAtomicNo(begin) === 6 && base.getAtomicNo(end) === 6) {
      ccSingleBonds.push([Math.min(begin, end), Math.max(begin, end)]);
    }
  }
  if (ccSingleBonds.length < totalChains) {
    return [];
  }
  const chainTypes = Object.keys(chainCounts);
  const results = new Set();
  for (const combo of combinations(ccSingleBonds, totalChains)) {
    const counts = { ...chainCounts };
    for (const assignment of enumerateChainAssignments(chainTypes, counts, combo.length)) {
      const clone = cloneEditableMol(base);
      try {
        combo.forEach(([a, b], idx) => {
          const bondIdx = clone.getBond(a, b);
          if (bondIdx < 0) {
            throw new Error('missing bond');
          }
          clone.deleteBond(bondIdx);
          insertSulfurChain(clone, a, b, assignment[idx]);
        });
        ensureHelperArrays(clone, OCL.Molecule.cHelperNeighbours);
        const { doubles, triples } = countUnsaturations(clone);
        const expectedDouble = assignment.reduce((sum, type) => {
          if (type === 'S_Chain_Tetra') {
            return sum + 1;
          }
          if (type === 'S_Chain_Hexa') {
            return sum + 2;
          }
          return sum;
        }, nDoubleBonds);
        if (doubles !== expectedDouble || triples !== nTripleBonds) {
          continue;
        }
        results.add(canonicalizeSmiles(clone.toIsomericSmiles()));
      } catch (err) {
        // skip invalid sulfur chain placement
      }
    }
  }
  return Array.from(results).sort();
}

function addSulfurFunctionalGroup(mol, carbonIdx, group) {
  const sIdx = addSimpleAtom(mol, 'S');
  addBondWithOrder(mol, carbonIdx, sIdx, 1);
  if (group === 'S_Bivalent') {
    addExplicitHydrogen(mol, sIdx);
    return;
  }
  if (group === 'S_Tetravalent') {
    const oIdx = addSimpleAtom(mol, 'O');
    addBondWithOrder(mol, sIdx, oIdx, 2);
    addExplicitHydrogen(mol, sIdx);
    return;
  }
  if (group === 'S_Hexavalent') {
    const o1Idx = addSimpleAtom(mol, 'O');
    const o2Idx = addSimpleAtom(mol, 'O');
    addBondWithOrder(mol, sIdx, o1Idx, 2);
    addBondWithOrder(mol, sIdx, o2Idx, 2);
  }
}

function attachFunctionalGroup(mol, carbonIdx, group) {
  switch (group) {
    case 'OH': {
      const oIdx = addSimpleAtom(mol, 'O');
      const hIdx = addSimpleAtom(mol, 'H');
      addBondWithOrder(mol, carbonIdx, oIdx, 1);
      addBondWithOrder(mol, oIdx, hIdx, 1);
      break;
    }
    case 'COOH': {
      const cIdx = addSimpleAtom(mol, 'C');
      const o1Idx = addSimpleAtom(mol, 'O');
      const o2Idx = addSimpleAtom(mol, 'O');
      const hIdx = addSimpleAtom(mol, 'H');
      addBondWithOrder(mol, carbonIdx, cIdx, 1);
      addBondWithOrder(mol, cIdx, o1Idx, 2);
      addBondWithOrder(mol, cIdx, o2Idx, 1);
      addBondWithOrder(mol, o2Idx, hIdx, 1);
      break;
    }
    case 'CHO': {
      const cIdx = addSimpleAtom(mol, 'C');
      const oIdx = addSimpleAtom(mol, 'O');
      const hIdx = addSimpleAtom(mol, 'H');
      addBondWithOrder(mol, carbonIdx, cIdx, 1);
      addBondWithOrder(mol, cIdx, oIdx, 2);
      addBondWithOrder(mol, cIdx, hIdx, 1);
      break;
    }
    case 'CN': {
      const cIdx = addSimpleAtom(mol, 'C');
      const nIdx = addSimpleAtom(mol, 'N');
      addBondWithOrder(mol, carbonIdx, cIdx, 1);
      addBondWithOrder(mol, cIdx, nIdx, 3);
      break;
    }
    case 'NC': {
      const nIdx = addSimpleAtom(mol, 'N');
      const cIdx = addSimpleAtom(mol, 'C');
      addBondWithOrder(mol, carbonIdx, nIdx, 1);
      addBondWithOrder(mol, nIdx, cIdx, 3);
      setAtomCharge(mol, nIdx, 1);
      setAtomCharge(mol, cIdx, -1);
      break;
    }
    case 'OCN': {
      const oIdx = addSimpleAtom(mol, 'O');
      const cIdx = addSimpleAtom(mol, 'C');
      const nIdx = addSimpleAtom(mol, 'N');
      addBondWithOrder(mol, carbonIdx, oIdx, 1);
      addBondWithOrder(mol, oIdx, cIdx, 1);
      addBondWithOrder(mol, cIdx, nIdx, 3);
      break;
    }
    case 'NCO': {
      const nIdx = addSimpleAtom(mol, 'N');
      const cIdx = addSimpleAtom(mol, 'C');
      const oIdx = addSimpleAtom(mol, 'O');
      addBondWithOrder(mol, carbonIdx, nIdx, 1);
      addBondWithOrder(mol, nIdx, cIdx, 2);
      addBondWithOrder(mol, cIdx, oIdx, 2);
      break;
    }
    case 'NO2': {
      const nIdx = addSimpleAtom(mol, 'N');
      const o1Idx = addSimpleAtom(mol, 'O');
      const o2Idx = addSimpleAtom(mol, 'O');
      addBondWithOrder(mol, carbonIdx, nIdx, 1);
      addBondWithOrder(mol, nIdx, o1Idx, 2);
      addBondWithOrder(mol, nIdx, o2Idx, 1);
      setAtomCharge(mol, nIdx, 1);
      setAtomCharge(mol, o2Idx, -1);
      break;
    }
    case 'COOR_CH3':
    case 'COOR_C2H5':
    case 'COOR_C3H7':
    case 'COOR_CH(CH3)2': {
      const descriptor = group.slice('COOR_'.length);
      attachEsterGroup(mol, carbonIdx, descriptor);
      break;
    }
    case 'COX_Cl':
    case 'COX_Br':
    case 'COX_F':
    case 'COX_I': {
      const halogen = group.split('_')[1];
      attachAcidHalideGroup(mol, carbonIdx, halogen);
      break;
    }
    case 'OX_Cl':
    case 'OX_Br':
    case 'OX_F':
    case 'OX_I': {
      const halogen = group.split('_')[1];
      attachHypohaliteGroup(mol, carbonIdx, halogen);
      break;
    }
    case 'Imine': {
      const nIdx = addSimpleAtom(mol, 'N');
      addBondWithOrder(mol, carbonIdx, nIdx, 2);
      break;
    }
    case 'Amide': {
      const cIdx = addSimpleAtom(mol, 'C');
      const oIdx = addSimpleAtom(mol, 'O');
      const nIdx = addSimpleAtom(mol, 'N');
      addBondWithOrder(mol, carbonIdx, cIdx, 1);
      addBondWithOrder(mol, cIdx, oIdx, 2);
      addBondWithOrder(mol, cIdx, nIdx, 1);
      addExplicitHydrogen(mol, nIdx);
      break;
    }
    case 'NH2': {
      const nIdx = addSimpleAtom(mol, 'N');
      addBondWithOrder(mol, carbonIdx, nIdx, 1);
      addExplicitHydrogen(mol, nIdx);
      addExplicitHydrogen(mol, nIdx);
      break;
    }
    case 'Azide': {
      const n1Idx = addSimpleAtom(mol, 'N');
      const n2Idx = addSimpleAtom(mol, 'N');
      const n3Idx = addSimpleAtom(mol, 'N');
      addBondWithOrder(mol, carbonIdx, n1Idx, 1);
      addBondWithOrder(mol, n1Idx, n2Idx, 2);
      addBondWithOrder(mol, n2Idx, n3Idx, 1);
      break;
    }
    case 'Ketone': {
      const oIdx = addSimpleAtom(mol, 'O');
      addBondWithOrder(mol, carbonIdx, oIdx, 2);
      break;
    }
    case 'F':
    case 'Cl':
    case 'Br':
    case 'I': {
      const halIdx = addSimpleAtom(mol, group);
      addBondWithOrder(mol, carbonIdx, halIdx, 1);
      break;
    }
    case 'S_Bivalent':
    case 'S_Tetravalent':
    case 'S_Hexavalent':
      addSulfurFunctionalGroup(mol, carbonIdx, group);
      break;
    default:
      break;
  }
}

function iterateAssignments(pool, slots, callback) {
  if (slots <= 0) {
    callback([]);
    return 1;
  }
  const indices = new Array(slots).fill(0);
  let produced = 0;
  const limit = MAX_ASSIGNMENT_ATTEMPTS;
  while (produced < limit) {
    callback(indices.map((idx) => pool[idx]));
    produced += 1;
    let position = slots - 1;
    while (position >= 0) {
      indices[position] += 1;
      if (indices[position] < pool.length) {
        break;
      }
      indices[position] = 0;
      position -= 1;
    }
    if (position < 0) {
      break;
    }
  }
  return produced;
}

function buildZeroCarbonFragments(functionalGroups) {
  const fragments = [];
  const dictionary = {
    COOH: 'C(=O)O',
    CHO: 'C=O',
    COOR_CH3: 'C(=O)OC',
    COOR_C2H5: 'C(=O)OCC',
    COOR_C3H7: 'C(=O)OCCC',
    'COOR_CH(CH3)2': 'C(=O)OC(C)C',
    COX_Cl: 'C(=O)Cl',
    COX_Br: 'C(=O)Br',
    COX_F: 'C(=O)F',
    COX_I: 'C(=O)I',
    OX_Cl: 'OCl',
    OX_Br: 'OBr',
    OX_F: 'OF',
    OX_I: 'OI',
    NH2: 'N'
  };
  for (const group of functionalGroups) {
    if (dictionary[group]) {
      fragments.push(dictionary[group]);
    }
  }
  return Array.from(new Set(fragments));
}

function generateEtherStructures(backboneSmiles, nEthers, nDoubleBonds, nTripleBonds) {
  if (!nEthers) {
    return [backboneSmiles];
  }
  let base;
  try {
    base = OCL.Molecule.fromSmiles(backboneSmiles);
  } catch (err) {
    return [];
  }
  const ccSingleBonds = [];
  const bondCount = base.getAllBonds();
  for (let bond = 0; bond < bondCount; bond++) {
    const begin = base.getBondAtom(0, bond);
    const end = base.getBondAtom(1, bond);
    if (base.getBondOrder(bond) === 1 && base.getAtomicNo(begin) === 6 && base.getAtomicNo(end) === 6) {
      ccSingleBonds.push([Math.min(begin, end), Math.max(begin, end)]);
    }
  }
  if (ccSingleBonds.length < nEthers) {
    return [];
  }
  const results = new Set();
  for (const combo of combinations(ccSingleBonds, nEthers)) {
    const clone = cloneEditableMol(base);
    try {
      for (const [a, b] of combo) {
        const bondIdx = clone.getBond(a, b);
        if (bondIdx < 0) {
          throw new Error('bond missing');
        }
        clone.deleteBond(bondIdx);
        const oIdx = addSimpleAtom(clone, 'O');
        addBondWithOrder(clone, a, oIdx, 1);
        addBondWithOrder(clone, oIdx, b, 1);
      }
      ensureHelperArrays(clone, OCL.Molecule.cHelperNeighbours);
      results.add(canonicalizeSmiles(clone.toIsomericSmiles()));
    } catch (err) {
      // skip invalid structures
    }
  }
  return Array.from(results).sort();
}

function addFunctionalGroups(backboneSmiles, functionalGroups, carbonTypes = DEFAULT_CARBON_TYPES, nDoubleBonds = 0, nTripleBonds = 0) {
  if (!functionalGroups || functionalGroups.length === 0) {
    return [backboneSmiles];
  }
  if (!backboneSmiles) {
    return buildZeroCarbonFragments(functionalGroups);
  }
  const groups = functionalGroups.slice();
  const etherCount = groups.filter((fg) => fg === 'Ether').length;
  if (etherCount > 0) {
    const remaining = groups.filter((fg) => fg !== 'Ether');
    const etherStructures = generateEtherStructures(backboneSmiles, etherCount, nDoubleBonds, nTripleBonds);
    const targets = etherStructures.length ? etherStructures : [backboneSmiles];
    const downstream = new Set();
    for (const structure of targets) {
      const nested = addFunctionalGroups(structure, remaining, carbonTypes, nDoubleBonds, nTripleBonds);
      nested.forEach((smiles) => downstream.add(canonicalizeSmiles(smiles)));
    }
    return Array.from(downstream);
  }
  const sulfurChainGroups = groups.filter((fg) => SULFUR_CHAIN_GROUPS.has(fg));
  if (sulfurChainGroups.length) {
    const remaining = groups.filter((fg) => !SULFUR_CHAIN_GROUPS.has(fg));
    const counts = {};
    sulfurChainGroups.forEach((group) => {
      counts[group] = (counts[group] || 0) + 1;
    });
    const sulfurStructures = generateSulfurChainStructures(backboneSmiles, counts, nDoubleBonds, nTripleBonds);
    if (!sulfurStructures.length) {
      return [];
    }
    const downstream = new Set();
    for (const structure of sulfurStructures) {
      const nested = addFunctionalGroups(structure, remaining, carbonTypes, nDoubleBonds, nTripleBonds);
      nested.forEach((smiles) => downstream.add(canonicalizeSmiles(smiles)));
    }
    return Array.from(downstream);
  }
  let baseMol;
  try {
    baseMol = OCL.Molecule.fromSmiles(backboneSmiles);
  } catch (err) {
    return [];
  }
  const carbonAtoms = classifyCarbonAtoms(baseMol, carbonTypes);
  if (!carbonAtoms.length) {
    return [];
  }
  const unique = new Set();
  iterateAssignments(carbonAtoms, groups.length, (assignment) => {
    const candidate = cloneEditableMol(baseMol);
    try {
      groups.forEach((group, idx) => {
        attachFunctionalGroup(candidate, assignment[idx].idx, group);
      });
      ensureHelperArrays(candidate, OCL.Molecule.cHelperNeighbours);
      unique.add(canonicalizeSmiles(candidate.toIsomericSmiles()));
    } catch (err) {
      // invalid arrangement, skip
    }
  });
  return Array.from(unique).sort();
}

function generateFunctionalizedIsomers(nCarbons, functionalGroups, nDoubleBonds = 0, nTripleBonds = 0, nRings = 0, carbonTypes = null) {
  const carbonTypeList = carbonTypes && carbonTypes.length ? carbonTypes : DEFAULT_CARBON_TYPES;
  const groups = (functionalGroups || []).slice();
  if (nCarbons === 0) {
    return buildZeroCarbonFragments(groups);
  }
  if (groups.length === 0) {
    return generateHydrocarbonIsomers(nCarbons, nDoubleBonds, nTripleBonds, nRings);
  }
  if (nCarbons < 2) {
    for (let i = groups.length - 1; i >= 0; i--) {
      if (groups[i] === 'Ether' || SULFUR_CHAIN_GROUPS.has(groups[i])) {
        groups.splice(i, 1);
      }
    }
  } else if (nCarbons === 2) {
    let seen = 0;
    for (let i = groups.length - 1; i >= 0; i--) {
      if (groups[i] === 'Ether') {
        seen += 1;
        if (seen > 1) {
          groups.splice(i, 1);
        }
      }
    }
  }
  if (groups.length === 0) {
    return generateHydrocarbonIsomers(nCarbons, nDoubleBonds, nTripleBonds, nRings);
  }
  const backbones = generateHydrocarbonIsomers(nCarbons, nDoubleBonds, nTripleBonds, nRings);
  const decorated = new Set();
  for (const backbone of backbones) {
    const molecules = addFunctionalGroups(backbone, groups, carbonTypeList, nDoubleBonds, nTripleBonds);
    molecules.forEach((smiles) => decorated.add(canonicalizeSmiles(smiles)));
  }
  return Array.from(decorated).sort();
}

function validateStructurePossibility(carbonCount, functionalGroups, nDoubleBonds, nTripleBonds, carbonTypes, nRings = 0) {
  const groups = functionalGroups || [];
  const etherCount = groups.filter((fg) => fg === 'Ether').length;
  if (etherCount > 0 && carbonCount < 2) {
    return [false, `IMPOSSIBLE: Ethers require at least 2 carbon atoms (you have ${carbonCount})`];
  }
  const sulfurChainCount = groups.filter((fg) => SULFUR_CHAIN_GROUPS.has(fg)).length;
  if (sulfurChainCount > 0) {
    if (carbonCount < 2) {
      return [false, `IMPOSSIBLE: Sulfur chains require at least 2 carbon atoms (you have ${carbonCount})`];
    }
    const maxChains = Math.max(carbonCount - 1, 0);
    if (sulfurChainCount > maxChains) {
      return [false, `IMPOSSIBLE: Maximum sulfur chains for ${carbonCount} carbon(s) is ${maxChains}`];
    }
  }
  if (groups.includes('Ketone')) {
    if (carbonCount < 2) {
      return [false, `IMPOSSIBLE: Ketones require at least 2 carbon atoms (you have ${carbonCount})`];
    }
    if (carbonCount === 2 && (nDoubleBonds > 0 || nTripleBonds > 0)) {
      return [false, 'IMPOSSIBLE: For 2 carbons with ketone, no double/triple bonds allowed'];
    }
  }
  if (carbonCount === 0 && (nDoubleBonds > 0 || nTripleBonds > 0 || nRings > 0)) {
    return [false, 'IMPOSSIBLE: Cannot have double/triple bonds or rings with zero carbons'];
  }
  if (nRings > 0 && carbonCount < 3) {
    return [false, `IMPOSSIBLE: Need at least 3 carbons to form a ring (you have ${carbonCount})`];
  }
  if (carbonCount > 0) {
    const totalUnsaturation = nDoubleBonds + nTripleBonds;
    const maxUnsaturation = nRings > 0 ? carbonCount + nRings - 1 : carbonCount - 1;
    if (totalUnsaturation > maxUnsaturation) {
      return [false, `IMPOSSIBLE: For ${carbonCount} carbon(s) and ${nRings} ring(s), maximum unsaturations is ${maxUnsaturation}`];
    }
    if (!carbonTypes || carbonTypes.length === 0) {
      return [false, 'IMPOSSIBLE: At least one carbon type must be selected when carbon count > 0'];
    }
  }
  return [true, ''];
}

function getElementCounts(input) {
  let mol = input;
  if (typeof input === 'string') {
    try {
      mol = OCL.Molecule.fromSmiles(input);
    } catch (err) {
      return {};
    }
  }
  if (!mol) {
    return {};
  }
  const counts = {};
  const atomCount = mol.getAllAtoms();
  for (let atom = 0; atom < atomCount; atom++) {
    const label = mol.getAtomLabel(atom);
    counts[label] = (counts[label] || 0) + 1;
  }
  return counts;
}

function compileFunctionalGroupQueries() {
  if (functionalGroupQueryCache.size > 0) {
    return functionalGroupQueryCache;
  }
  Object.entries(FUNCTIONAL_GROUP_PATTERNS).forEach(([name, smiles]) => {
    try {
      const fragment = OCL.Molecule.fromSmiles(smiles);
      fragment.setFragment(true);
      functionalGroupQueryCache.set(name, fragment);
    } catch (err) {
      // ignore invalid fragment definitions
    }
  });
  return functionalGroupQueryCache;
}

function getFunctionalGroupType(input) {
  let mol = input;
  if (typeof input === 'string') {
    try {
      mol = OCL.Molecule.fromSmiles(input);
    } catch (err) {
      return 'Unknown';
    }
  }
  if (!mol) {
    return 'Unknown';
  }
  const matches = [];
  const queries = compileFunctionalGroupQueries();
  const searcher = new OCL.SSSearcher();
  for (const [label, fragment] of queries.entries()) {
    if (!fragment) {
      continue;
    }
    searcher.setMol(fragment, mol);
    if (searcher.isFragmentInMolecule()) {
      matches.push(label);
    }
  }
  return matches.length ? matches.join(', ') : 'Hydrocarbon';
}

const rdMolDescriptors = {
  CalcMolFormula: (input) => {
    const counts = getElementCounts(input);
    const ordering = ['C', 'H', 'O', 'N', 'S', 'F', 'Cl', 'Br', 'I'];
    let formula = '';
    for (const element of ordering) {
      if (counts[element]) {
        formula += element + (counts[element] > 1 ? counts[element] : '');
      }
    }
    Object.keys(counts)
      .filter((key) => !ordering.includes(key))
      .sort()
      .forEach((element) => {
        formula += element + (counts[element] > 1 ? counts[element] : '');
      });
    return formula;
  }
};

function computeWeisfeilerLehmanSignature(mol, iterations = WEISFEILER_ITERATIONS) {
  const graph = parseMolGraph(mol);
  let labels = graph.atoms.map((atom) => atom.label || '?');
  for (let iter = 0; iter < iterations; iter++) {
    labels = labels.map((label, idx) => {
      const neighborLabels = graph.adjacency[idx]
        .map((neighbor) => `${labels[neighbor.idx]}${neighbor.order}`)
        .sort()
        .join('');
      return `${label}[${neighborLabels}]`;
    });
  }
  return labels.sort().join('|');
}

module.exports = {
  canonicalizeSmiles,
  generateAlkaneBackbone,
  generateHydrocarbonIsomers,
  generateFunctionalizedIsomers,
  validateStructurePossibility,
  getFunctionalGroupType,
  getElementCounts,
  rdMolDescriptors,
  addExactUnsaturation,
  addFunctionalGroups,
  __testHooks: {
    computeWeisfeilerLehmanSignature,
    _generateTreeEdgeSets: generateTreeEdgeSets
  }
};
