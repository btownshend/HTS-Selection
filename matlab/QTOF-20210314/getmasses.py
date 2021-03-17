from rdkit import Chem
from rdkit.Chem import Descriptors

suppl = Chem.SDMolSupplier('../../data/HTBCFiles/ChemDIvFull.sdf')

for mol in suppl:
    plate=mol.GetProp('BATCH_PLATE')
    well=mol.GetProp('BATCH_WELL')
    smiles=Chem.MolToSmiles(mol)
    mass=Descriptors.ExactMolWt(mol)
    print(f"{mass}, {plate}-{well}")

