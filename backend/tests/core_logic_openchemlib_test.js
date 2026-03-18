#!/usr/bin/env node
'use strict';

const assert = require('assert');
const {
  canonicalizeSmiles,
  generateHydrocarbonIsomers,
  generateFunctionalizedIsomers
} = require('../core_logic');

const EXPECTED_C6_ALKENES = [
  'CC(C)(C)C=C',
  'CC(C)=C(C)C',
  'CC(C)C(C)=C',
  'CC(C)C=CC',
  'CC(C)CC=C',
  'CCC(C)=CC',
  'CCC(C)C=C',
  'CCC(CC)=C',
  'CCC=C(C)C',
  'CCC=CCC',
  'CCCC(C)=C',
  'CCCC=CC',
  'CCCCC=C'
];

function testHydrocarbonEnumeration() {
  const actual = generateHydrocarbonIsomers(6, 1, 0, 0);
  assert.deepStrictEqual(actual, EXPECTED_C6_ALKENES, 'C6H12 alkene enumeration must match Extrafast reference set');
}

function testAlcoholFunctionalization() {
  const molecules = generateFunctionalizedIsomers(6, ['OH'], 1, 0, 0);
  assert.ok(molecules.includes('CCCCC=CO'), 'Primary alcohol should be present');
  assert.ok(molecules.includes('C=CCCCCO'), 'Terminal double-bond alcohol should be present');
  const canonicalSet = new Set(molecules.map((smiles) => canonicalizeSmiles(smiles)));
  assert.strictEqual(canonicalSet.size, molecules.length, 'Functionalization output must avoid duplicate canonical SMILES');
}

function testAcylAndHypohaliteGroups() {
  const ester = generateFunctionalizedIsomers(2, ['COOR_CH3']);
  assert.ok(ester.includes('CCC(OC)=O'), 'Methyl ester should be generated');

  const acidHalide = generateFunctionalizedIsomers(2, ['COX_Cl']);
  assert.ok(acidHalide.includes('CCC(Cl)=O'), 'Acid chloride should be generated');

  const hypohalite = generateFunctionalizedIsomers(2, ['OX_Br']);
  assert.ok(hypohalite.includes('CCOBr'), 'Hypobromite substitution should be generated');
}

function testSulfurChains() {
  const sulfurChains = generateFunctionalizedIsomers(4, ['S_Chain_Bi']);
  assert.ok(sulfurChains.includes('CCSCC'), 'Sulfur bridge should break a single C-C bond and insert S');
}

function testZeroCarbonFragments() {
  const fragments = generateFunctionalizedIsomers(0, ['COOH', 'OX_Cl']);
  assert.ok(fragments.includes('C(=O)O'), 'Carboxylic acid fragment should be emitted for zero-carbon inputs');
  assert.ok(fragments.includes('OCl'), 'Hypochlorite fragment should be emitted for zero-carbon inputs');
}

function testRingFunctionalization() {
  const molecules = generateFunctionalizedIsomers(7, ['Ketone', 'Cl'], 0, 0, 1);
  assert.ok(molecules.length > 0, 'Ringed ketone/halide combinations should produce candidates');
  assert.ok(
    molecules.some((smiles) => smiles.includes('Cl') && smiles.includes('=O')),
    'At least one molecule should retain both ketone and halide decoration'
  );
}

function testHighCarbonFunctionalization() {
  const molecules = generateFunctionalizedIsomers(6, ['Ether', 'NO2', 'COOR_CH3'], 1, 0, 0);
  assert.ok(molecules.length > 0, 'Decorated skeletons should be generated for mid-sized chains');
  const canonicalSet = new Set(molecules.map((smiles) => canonicalizeSmiles(smiles)));
  assert.strictEqual(canonicalSet.size, molecules.length, 'High-carbon enumeration must remain duplicate free');
}

function testMixedSulfurChainsLarge() {
  const molecules = generateFunctionalizedIsomers(5, ['S_Chain_Bi', 'S_Chain_Tetra'], 1, 0, 0);
  assert.ok(molecules.length > 0, 'Mixed sulfur chains with unsaturation should yield candidates');
  assert.ok(
    molecules.some((smiles) => smiles.includes('S') && smiles.includes('=O')),
    'Sulfur chain outputs should retain oxidized sulfur signatures'
  );
}

const tests = [
  { name: 'Hydrocarbon enumeration', fn: testHydrocarbonEnumeration },
  { name: 'Alcohol functionalization', fn: testAlcoholFunctionalization },
  { name: 'Ester / halide / hypohalite attachments', fn: testAcylAndHypohaliteGroups },
  { name: 'Sulfur chain insertion', fn: testSulfurChains },
  { name: 'Zero-carbon fragments', fn: testZeroCarbonFragments },
  { name: 'Ring ketone functionalization', fn: testRingFunctionalization },
  { name: 'High-carbon decorated skeleton', fn: testHighCarbonFunctionalization },
  { name: 'Mixed sulfur chain macro test', fn: testMixedSulfurChainsLarge }
];

(async () => {
  for (const test of tests) {
    try {
      test.fn();
      console.log(`✔ ${test.name}`);
    } catch (err) {
      console.error(`✖ ${test.name} failed:`, err.message || err);
      throw err;
    }
  }
  console.log('\n✅ Core logic OpenChemLib tests passed.');
})().catch((err) => {
  console.error('❌ Core logic tests failed:', err);
  process.exit(1);
});
