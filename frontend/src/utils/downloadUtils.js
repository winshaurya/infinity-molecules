import { initRDKit, RDKitModule } from '@rdkit/rdkit';

let rdkitPromise = null;

export async function getRDKit() {
  if (!rdkitPromise) {
    rdkitPromise = initRDKit();
  }
  return await rdkitPromise;
}

export async function createDownloadBlob(smilesList, format, jobId, onProgress) {
  const rdkit = await getRDKit();

  if (format === 'csv') {
    // Create CSV with SMILES
    const csvContent = 'SMILES\n' + smilesList.join('\n');
    return new Blob([csvContent], { type: 'text/csv' });
  } else if (format === 'molsdf') {
    // Create ZIP with MOL, SDF, and CSV files
    const JSZip = (await import('jszip')).default;
    const zip = new JSZip();

    // Add CSV file
    const csvContent = 'SMILES\n' + smilesList.join('\n');
    zip.file(`molecules_${jobId}.csv`, csvContent);

    // Add SDF file
    let sdfContent = '';
    for (let i = 0; i < smilesList.length; i++) {
      const smiles = smilesList[i];
      try {
        const mol = rdkit.get_mol(smiles);
        if (mol) {
          const molBlock = mol.get_molblock();
          sdfContent += molBlock + '\n$$$$\n';
          mol.delete();
        }
      } catch (error) {
        console.error('Error converting SMILES to MOL:', smiles, error);
      }
      // Report progress for SDF creation (first 50% of progress)
      if (onProgress) {
        onProgress(Math.round((i + 1) / smilesList.length * 50));
      }
    }
    zip.file(`molecules_${jobId}.sdf`, sdfContent);

    // Add individual MOL files
    for (let i = 0; i < smilesList.length; i++) {
      const smiles = smilesList[i];
      try {
        const mol = rdkit.get_mol(smiles);
        if (mol) {
          const molBlock = mol.get_molblock();
          zip.file(`molecule_${(i + 1).toString().padStart(4, '0')}.mol`, molBlock);
          mol.delete();
        }
      } catch (error) {
        console.error('Error converting SMILES to MOL:', smiles, error);
      }
      // Report progress for MOL files creation (next 40% of progress)
      if (onProgress) {
        onProgress(Math.round(50 + (i + 1) / smilesList.length * 40));
      }
    }

    // Generate ZIP blob (last 10% of progress)
    if (onProgress) onProgress(90);
    const zipBlob = await zip.generateAsync({ type: 'blob' });
    if (onProgress) onProgress(100);
    return zipBlob;
  }

  throw new Error(`Unsupported format: ${format}`);
}
