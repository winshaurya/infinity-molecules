#!/usr/bin/env node
'use strict';
const { spawnSync } = require('child_process');
const path = require('path');
const { generateFunctionalizedIsomers, validateStructurePossibility, canonicalizeSmiles } = require('../core_logic');

const SCENARIOS = [
  {
    name: 'C6 double-bond alcohol',
    carbon_count: 6,
    double_bonds: 1,
    triple_bonds: 0,
    rings: 0,
    functional_groups: ['OH']
  },
  {
    name: 'C8 ring ketone halide',
    carbon_count: 8,
    double_bonds: 1,
    triple_bonds: 0,
    rings: 1,
    functional_groups: ['Ketone', 'Cl']
  },
  {
    name: 'C10 ether + nitro + ester',
    carbon_count: 10,
    double_bonds: 2,
    triple_bonds: 0,
    rings: 0,
    functional_groups: ['Ether', 'NO2', 'COOR_CH3']
  },
  {
    name: 'C7 sulfur chain mix',
    carbon_count: 7,
    double_bonds: 0,
    triple_bonds: 0,
    rings: 0,
    functional_groups: ['S_Chain_Bi', 'S_Chain_Tetra', 'OH']
  }
];

function runPythonExport(params) {
  const script = path.join(__dirname, 'python_smiles_export.py');
  const args = [
    script,
    '--carbon_count', String(params.carbon_count),
    '--double_bonds', String(params.double_bonds || 0),
    '--triple_bonds', String(params.triple_bonds || 0),
    '--rings', String(params.rings || 0),
    '--functional_groups', JSON.stringify(params.functional_groups || []),
    '--carbon_types', JSON.stringify(params.carbon_types || ['primary', 'secondary', 'tertiary'])
  ];
  const res = spawnSync('python', args, { cwd: path.join(__dirname, '..') });
  if (res.error) {
    throw res.error;
  }
  const out = res.stdout.toString().trim();
  try {
    const obj = JSON.parse(out);
    if (obj.error) {
      throw new Error('python_error: ' + obj.error);
    }
    return obj.smiles || [];
  } catch (e) {
    throw new Error('python_output_parse_error: ' + e + '\nstdout=' + out + '\nstderr=' + res.stderr.toString());
  }
}

function canonicalizeSet(smilesList) {
  return Array.from(new Set(smilesList.map((s) => canonicalizeSmiles(s)).filter(Boolean))).sort();
}

function diffSets(jsCanon, pyCanon) {
  const jsSet = new Set(jsCanon);
  const pySet = new Set(pyCanon);
  return {
    onlyJS: jsCanon.filter((entry) => !pySet.has(entry)),
    onlyPY: pyCanon.filter((entry) => !jsSet.has(entry))
  };
}

async function main() {
  const failures = [];
  for (const scenario of SCENARIOS) {
    console.log(`\nScenario: ${scenario.name}`);
    const [ok, reason] = validateStructurePossibility(
      scenario.carbon_count,
      scenario.functional_groups,
      scenario.double_bonds || 0,
      scenario.triple_bonds || 0,
      scenario.carbon_types || ['primary', 'secondary', 'tertiary'],
      scenario.rings || 0
    );
    if (!ok) {
      throw new Error(`Validation failed for ${scenario.name}: ${reason}`);
    }

    const jsCanonical = canonicalizeSet(
      generateFunctionalizedIsomers(
        scenario.carbon_count,
        scenario.functional_groups,
        scenario.double_bonds || 0,
        scenario.triple_bonds || 0,
        scenario.rings || 0,
        scenario.carbon_types || ['primary', 'secondary', 'tertiary']
      )
    );

    let pyCanonical = [];
    try {
      pyCanonical = canonicalizeSet(runPythonExport(scenario));
    } catch (err) {
      console.error(`Python export failed for ${scenario.name}:`, err.message);
      process.exit(3);
    }

    console.log(' JS count:', jsCanonical.length);
    console.log(' PY count:', pyCanonical.length);

    const { onlyJS, onlyPY } = diffSets(jsCanonical, pyCanonical);
    if (onlyJS.length || onlyPY.length) {
      failures.push({ scenario: scenario.name, onlyJS, onlyPY });
      console.warn('  ⚠️  Divergence detected');
      if (onlyJS.length) {
        console.warn('    Only JS (first 10):', onlyJS.slice(0, 10));
      }
      if (onlyPY.length) {
        console.warn('    Only PY (first 10):', onlyPY.slice(0, 10));
      }
    } else {
      console.log('  ✅ Sets match');
    }
  }

  if (failures.length) {
    console.error(`\n❌ ${failures.length} scenario(s) diverged.`);
    process.exit(1);
  }

  console.log('\n✅ All comparison scenarios matched.');
  process.exit(0);
}

if (require.main === module) {
  main();
}
