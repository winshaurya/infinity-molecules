import { initRDKit, RDKitModule } from '@rdkit/rdkit';

let rdkitPromise = null;

export async function getRDKit() {
  if (!rdkitPromise) {
    rdkitPromise = (async () => {
      // If we are in production, we need to tell RDKit where the WASM is
      const options = {
        locateFile: (file) => {
          // Point to the local copy in public folder
          return `/${file}`;
        }
      };
      
      const module = await import('@rdkit/rdkit');
      const init = module.initRDKit || module.default?.initRDKit || module.default;
      if (typeof init !== 'function') {
        const msg = 'RDKit initialization function not found';
        console.error(msg, module);
        throw new Error(msg);
      }
      return await init(options);
    })();
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
    // Create ZIP with a structured directory: per-molecule MOL/SDF + metadata
    const JSZip = (await import('jszip')).default;
    const zip = new JSZip();
    const rootFolder = `molecules_${jobId}`;
    const molFolder = `${rootFolder}/MOL Files`;
    const sdfFolder = `${rootFolder}/SDF Files`;

    // Add CSV file at root
    const csvContent = 'SMILES\n' + smilesList.join('\n');
    zip.file(`${rootFolder}/molecules_${jobId}.csv`, csvContent);

    // Metadata file
    const metadataLines = [];
    metadataLines.push(`Job ID: ${jobId}`);
    metadataLines.push(`Generated: ${new Date().toISOString()}`);
    metadataLines.push(`Total molecules: ${smilesList.length}`);
    metadataLines.push('');
    metadataLines.push('Index,SMILES');

    // Create per-molecule files
    for (let i = 0; i < smilesList.length; i++) {
        const smiles = smilesList[i];
        const indexStr = (i + 1).toString().padStart(4, '0');
        metadataLines.push(`${indexStr},${smiles}`);

        try {
            // Check if smiles is valid before processing
            if (!smiles || typeof smiles !== 'string') {
                console.warn(`Invalid SMILES at index ${i}:`, smiles);
                continue;
            }
            
            const mol = rdkit.get_mol(smiles);
            if (mol) {
                const molBlock = mol.get_molblock();
                if (molBlock) {
                    // Save individual MOL file
                    zip.file(`${molFolder}/molecule_${indexStr}.mol`, molBlock);

                    // Save individual SDF file (single entry)
                    const sdfBlock = molBlock + '\n$$$$\n';
                    zip.file(`${sdfFolder}/molecule_${indexStr}.sdf`, sdfBlock);
                }
                mol.delete();
            }
        } catch (error) {
            console.error('Error converting SMILES to MOL:', smiles, error);
        }

      // Report progress (structured into 90% of work)
      if (onProgress) {
        onProgress(Math.round((i + 1) / smilesList.length * 90));
      }
    }

    // Write metadata file
    zip.file(`${rootFolder}/metadata.txt`, metadataLines.join('\n'));

    // Generate ZIP blob (last 10% of progress)
    if (onProgress) onProgress(90);
    const zipBlob = await zip.generateAsync({ type: 'blob' });
    if (onProgress) onProgress(100);
    return zipBlob;
  }

  throw new Error(`Unsupported format: ${format}`);
}
