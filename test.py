import os
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# Define the path to save the molecule image
mol_image_path = "D_Glucose.png"

# Input a SMILES string
smiles = 'OC[C@H]1NC(=O)[C@H](C(C)C)N(c2c3c(C1)c[nH]c3c(cc2)[C@](CCC=C(C)C)(C=C)C)C'

# Create a molecule object from the SMILES string
mol = Chem.MolFromSmiles(smiles)

# Check if the molecule was created successfully
if mol:
    # Generate a 2D depiction of the molecule and save it as a PNG image
    Draw.MolToFile(mol, mol_image_path)

    # Load and display the generated image
    if os.path.exists(mol_image_path):
        img = Image.open(mol_image_path)
        img.show()
    else:
        print("Image file was not generated.")
else:
    print("Failed to create molecule object from SMILES string.")
