export function updateFunctionalGroupStates(carbonCount) {
  const disabled = {}

  const disabledForZeroCarbon = new Set([
    'Amide', 'Cl', 'Br', 'I', 'F', 'CN', 'NC', 'OCN', 'NCO', 'Imine', 'NO2', 'Ketone', 'Ether', 'OH', 'Azide',
    'S_Bivalent', 'S_Tetravalent', 'S_Hexavalent', 'S_Chain_Bi', 'S_Chain_Tetra', 'S_Chain_Hexa'
  ])

  const alwaysAllowed = new Set([
    'COOH', 'CHO', 'COOR_CH3', 'COOR_C2H5', 'COOR_C3H7', 'COOR_CH(CH3)2',
    'COX_Cl', 'COX_Br', 'COX_F', 'COX_I',
    'OX_Cl', 'OX_Br', 'OX_F', 'OX_I',
    'NH2'
  ])

  if (carbonCount === 0) {
    disabledForZeroCarbon.forEach(fg => {
      if (!alwaysAllowed.has(fg)) {
        disabled[fg] = true
      }
    })
  } else {
    if (carbonCount < 2) {
      disabled['Ether'] = true
      disabled['Ketone'] = true
      disabled['S_Chain_Bi'] = true
      disabled['S_Chain_Tetra'] = true
      disabled['S_Chain_Hexa'] = true
    }
    if (carbonCount === 1) {
      disabled['Ketone'] = true
    }
  }

  return disabled
}

export function updateValencyStatus(carbonCount, doubleBonds, tripleBonds, rings, functionalGroups) {
  if (carbonCount === 0) {
    return 'Valency status: OK'
  }

  let maxValency = 2 * carbonCount + 2 - 2 * doubleBonds - 4 * tripleBonds
  if (rings > 0) {
    maxValency += rings
  }

  let currentValency = 0
  Object.values(functionalGroups).forEach(count => {
    currentValency += count
  })

  if (currentValency > maxValency) {
    return `Valency status: EXCEEDED (Current: ${currentValency}, Max: ${maxValency})`
  } else if (currentValency === maxValency) {
    return `Valency status: MAXIMUM REACHED (Current: ${currentValency}, Max: ${maxValency})`
  } else {
    return `Valency status: OK (Current: ${currentValency}, Max: ${maxValency})`
  }
}

export function validateStructurePossibility({
  carbonCount = 0,
  functionalGroups = [],
  doubleBonds = 0,
  tripleBonds = 0,
  rings = 0,
  carbonTypes = []
}) {
  const countOccurrences = (name) => functionalGroups.filter(fg => fg === name).length
  const fail = (message) => ({ valid: false, message })

  if (functionalGroups.includes('Ether')) {
    const etherCount = countOccurrences('Ether')
    if (carbonCount < 2) {
      return fail(`IMPOSSIBLE: Ethers require at least 2 carbon atoms (you have ${carbonCount})`)
    }
    if (carbonCount === 2) {
      if (etherCount > 1) {
        return fail(`IMPOSSIBLE: For 2 carbon atoms, maximum 1 ether group is possible (you requested ${etherCount})`)
      }
      if (doubleBonds >= 1 && etherCount >= 1) {
        return fail('IMPOSSIBLE: For 2 carbon atoms, cannot have both double bonds and ethers')
      }
    }
    if (carbonCount >= 3) {
      const maxEthers = carbonCount - 1
      if (etherCount > maxEthers) {
        return fail(`IMPOSSIBLE: For ${carbonCount} carbon atoms, maximum ${maxEthers} ether groups are possible (you requested ${etherCount})`)
      }
    }
  }

  const sulfurChainTypes = ['S_Chain_Bi', 'S_Chain_Tetra', 'S_Chain_Hexa']
  const sulfurChainsPresent = functionalGroups.filter(fg => sulfurChainTypes.includes(fg))
  if (sulfurChainsPresent.length) {
    if (carbonCount < 2) {
      return fail(`IMPOSSIBLE: Sulfur chains require at least 2 carbon atoms (you have ${carbonCount})`)
    }
    if (sulfurChainsPresent.length > carbonCount - 1) {
      return fail(`IMPOSSIBLE: Maximum sulfur chains for ${carbonCount} carbons is ${carbonCount - 1} (you requested ${sulfurChainsPresent.length})`)
    }
  }

  if (functionalGroups.includes('Ketone')) {
    if (carbonCount < 2) {
      return fail(`IMPOSSIBLE: Ketones require at least 2 carbon atoms (you have ${carbonCount})`)
    }
    if (carbonCount === 2 && (doubleBonds > 0 || tripleBonds > 0)) {
      return fail('IMPOSSIBLE: For 2 carbons with ketone, no double/triple bonds allowed (would destroy the C-C bond needed for ketone)')
    }
    if (carbonCount === 1) {
      return fail('IMPOSSIBLE: Ketones cannot be formed with only 1 carbon atom')
    }
  }

  if (functionalGroups.includes('Azide') && carbonCount === 0) {
    return fail('IMPOSSIBLE: Azides require at least 1 carbon atom (you have 0)')
  }

  if (carbonCount > 0) {
    const totalUnsaturations = doubleBonds + tripleBonds
    const maxPossibleUnsaturations = rings > 0 ? carbonCount + rings - 1 : carbonCount - 1
    if (totalUnsaturations > maxPossibleUnsaturations) {
      return fail(`IMPOSSIBLE: For ${carbonCount} carbon(s) and ${rings} ring(s), maximum unsaturations is ${maxPossibleUnsaturations} (you requested ${totalUnsaturations})`)
    }
  } else {
    if (doubleBonds > 0 || tripleBonds > 0 || rings > 0) {
      return fail('IMPOSSIBLE: Cannot have double/triple bonds or rings with zero carbons')
    }
  }

  if (carbonCount === 0) {
    if (!functionalGroups.length) {
      return fail('IMPOSSIBLE: Zero carbons with no functional groups is nothing')
    }
    const validZeroCarbonGroups = new Set([
      'COOH', 'CHO',
      'COOR_CH3', 'COOR_C2H5', 'COOR_C3H7', 'COOR_CH(CH3)2',
      'COX_Cl', 'COX_Br', 'COX_F', 'COX_I',
      'OX_Cl', 'OX_Br', 'OX_F', 'OX_I',
      'NH2'
    ])
    const invalidGroupsForZeroCarbon = new Set([
      'Amide', 'Cl', 'Br', 'I', 'F', 'CN', 'NC',
      'OCN', 'NCO', 'Imine', 'NO2', 'Ketone', 'Ether', 'OH', 'Azide',
      'S_Bivalent', 'S_Tetravalent', 'S_Hexavalent',
      'S_Chain_Bi', 'S_Chain_Tetra', 'S_Chain_Hexa'
    ])
    const invalidGroups = []
    functionalGroups.forEach(fg => {
      if (invalidGroupsForZeroCarbon.has(fg)) {
        invalidGroups.push(fg)
      } else if (!validZeroCarbonGroups.has(fg)) {
        invalidGroups.push(fg)
      }
    })
    if (invalidGroups.length) {
      const unique = [...new Set(invalidGroups)]
      return fail(`IMPOSSIBLE: For zero carbons, these functional groups are invalid: ${unique.join(', ')}`)
    }
  }

  if (carbonCount > 0) {
    const totalFgCount = functionalGroups.length
    if (totalFgCount > 0) {
      if (carbonCount === 1) {
        if (totalFgCount > 4) {
          return fail(`IMPOSSIBLE: Single carbon can have maximum 4 functional groups (you selected ${totalFgCount})`)
        }
      } else {
        let maxFgEstimate = carbonCount * 2.5
        if (carbonCount <= 3) {
          maxFgEstimate = carbonCount * 4
        } else if (carbonCount <= 6) {
          maxFgEstimate = carbonCount * 3
        }
        if (totalFgCount > maxFgEstimate) {
          return fail(`IMPOSSIBLE: ${carbonCount} carbons can accommodate approximately ${Math.floor(maxFgEstimate)} functional groups maximum (you selected ${totalFgCount})`)
        }
      }
    }
  }

  if (carbonCount > 0 && (!carbonTypes || carbonTypes.length === 0)) {
    return fail('IMPOSSIBLE: At least one carbon type must be selected when carbon count > 0')
  }

  if (rings > 0) {
    if (carbonCount < 3) {
      return fail(`IMPOSSIBLE: Need at least 3 carbons to form a ring (you have ${carbonCount})`)
    }
    if (rings > carbonCount - 1) {
      return fail(`IMPOSSIBLE: Maximum rings for ${carbonCount} carbons is ${carbonCount - 2}`)
    }
  }

  return { valid: true }
}