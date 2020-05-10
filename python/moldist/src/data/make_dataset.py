import os

from rdkit.Chem import AllChem


def load():
    # Read in the CDIV molecules
    suppl = AllChem.SDMolSupplier(os.getenv("DATA") + "/HTBCFiles/ChemDivFull.sdf")

    # Extract only the 960 tested molecules
    plates = ["CDIV%04d" % p for p in range(1, 121, 10)]
    # print(len(plates), plates)
    tested = [x for x in suppl if x.GetProp("BATCH_PLATE") in plates]
    # print(len(suppl), len(tested))

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
