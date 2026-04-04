from rdkit import Chem

# Tạo phân tử từ SMILES
mol = Chem.MolFromSmiles('CCO') # Ethanol

# Chuyển sang InChI
inchi = Chem.MolToInchi(mol)
print(f"Mã InChI: {inchi}")

# Chuyển sang InChIKey
inchikey = Chem.MolToInchiKey(mol)
print(f"Mã InChIKey: {inchikey}")
