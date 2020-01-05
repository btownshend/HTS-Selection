import numpy as np
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures

def build(molecules, fpSize=2048):
    # Setup fingerprints for all the tested molecules
    fp = []
    for m in molecules:
        # fp.append(AllChem.RDKFingerprint(m, fpSize=fpSize))
        fp.append(AllChem.GetMorganFingerprintAsBitVect(m,radius=10,nBits=fpSize))

    # print(fp)
    X = np.array(list(fp))
    print("mean value of features =",np.mean(np.mean(X)))
    return X
