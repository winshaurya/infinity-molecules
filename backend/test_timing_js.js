#!/usr/bin/env node
'use strict';

/**
 * JavaScript timing test for molecule generation
 */

const { hrtime } = require('process');
const { generateFunctionalizedIsomers } = require('./core_logic');

function testGeneration() {
  console.log("Testing JavaScript molecule generation...");

  // Test parameters for ~100k molecules - same as Python
  const carbonCount = 12;
  const functionalGroups = ["OH", "COOH", "NH2"];
  const doubleBonds = 2;
  const tripleBonds = 0;
  const rings = 1;
  const carbonTypes = ["primary", "secondary", "tertiary"];

  console.log(`Parameters: ${carbonCount}C, ${doubleBonds}DB, ${tripleBonds}TB, ${rings}R, ${functionalGroups}, ${carbonTypes}`);

  const startTime = hrtime.bigint();
  const smilesList = generateFunctionalizedIsomers(
    carbonCount, functionalGroups, doubleBonds, tripleBonds, rings, carbonTypes
  );
  const endTime = hrtime.bigint();

  const generationTime = Number(endTime - startTime) / 1e9;
  const moleculeCount = smilesList.length;

  console.log(`Generation time: ${generationTime.toFixed(2)} seconds`);
  console.log(`Molecules generated: ${moleculeCount}`);
  console.log(`Performance: ${(moleculeCount / generationTime).toFixed(1)} molecules/second`);

  // Show first 5 SMILES
  console.log("\nFirst 5 SMILES:");
  for (let i = 0; i < Math.min(5, smilesList.length); i++) {
    console.log(`  ${i + 1}. ${smilesList[i]}`);
  }

  return [generationTime, moleculeCount, smilesList];
}

if (require.main === module) {
  testGeneration();
}
