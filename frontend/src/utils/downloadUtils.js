import { initRDKit, RDKitModule } from '@rdkit/rdkit';

let rdkitPromise = null;

export async function getRDKit() {
  if (!rdkitPromise) {
    rdkitPromise = initRDKit();
  }
  return await rdkitPromise;
}

export async function createDownloadBlob(smilesList, format, jobId) {
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
    for (const smiles of smilesList) {
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
    }
    zip.file(`molecules_${jobId}.sdf`, sdfContent);

    // Add individual MOL files
    smilesList.forEach((smiles, index) => {
      try {
        const mol = rdkit.get_mol(smiles);
        if (mol) {
          const molBlock = mol.get_molblock();
          zip.file(`molecule_${(index + 1).toString().padStart(4, '0')}.mol`, molBlock);
          mol.delete();
        }
      } catch (error) {
        console.error('Error converting SMILES to MOL:', smiles, error);
      }
    });

    // Generate ZIP blob
    const zipBlob = await zip.generateAsync({ type: 'blob' });
    return zipBlob;
  }

  throw new Error(`Unsupported format: ${format}`);
}
