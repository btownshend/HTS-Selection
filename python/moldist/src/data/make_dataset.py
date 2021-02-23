import os

from rdkit.Chem import AllChem
import pubchempy as pc

def load():
    # Read in the CDIV molecules
    suppl = AllChem.SDMolSupplier(os.getenv("DATA") + "/HTBCFiles/ChemDivFull.sdf")

    # Extract only the 960 tested molecules
    # plates = ["CDIV%04d" % p for p in range(1, 121, 10)]
    # print(len(plates), plates)
    wells=["045B02","206E05","208H04","247B08","247C06","247E05","286E11","287E04","326C10","445F10","446A08","446G07","446G08","447E02","448H09","485D02","486E07","487E04","487F04","488D10","525D02","526A09","526B09","526D09","526H07","526H11","568B08","568C04","568C07","607C11"]
    # tested = [x for x in suppl if x.GetProp("BATCH_PLATE") in plates]
    tested = [x for x in suppl if x.GetProp("BATCH_PLATE")[-3:]+x.GetProp("BATCH_WELL") in wells]
    print(len(suppl), len(tested))

    # Set NAME property of molecules to match data from Matlab
    for mol in tested:
        if mol is None:
            continue
        plate = mol.GetProp("BATCH_PLATE")
        plate = int(plate[5:])
        well = mol.GetProp("BATCH_WELL")
        # if well[-2] == '0':
        # Remove leading 0 from column number
        # well = well[:-2] + well[-1]
        name = "%d%s" % (plate, well)
        # print(plate, well, name)
        mol.SetProp("NAME", name)

    return tested

def dumpInchi(molecules):
    for mol in molecules:
        AllChem.AddHs(mol)  # Not sure if this is needed, but seems idempotent
        inchi=AllChem.MolToInchi(mol)
        smiles=AllChem.MolToSmiles(mol)
        print(f"{mol.GetProp('NAME')}\n{inchi}\n{smiles}\n\n")
        #c=pc.get_compounds(inchi,'inchi')
        #print(c)
