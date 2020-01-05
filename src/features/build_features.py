import numpy as np
from rdkit.Chem import AllChem


def build(molecules, fpSize=2048):
    # Setup fingerprints for all the tested molecules
    fp = []
    for m in molecules:
        fp.append(AllChem.RDKFingerprint(m, fpSize=fpSize))
    # print(fp)
    X = np.array(list(fp))
    print("mean value of features =",np.mean(np.mean(X)))
    return X
