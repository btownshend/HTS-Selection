import os
import numpy as np
from rdkit import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem import FragmentCatalog


def build(molecules, fpSize=2048):
    # Setup fingerprints for all the tested molecules
    fp = []
    for m in molecules:
        # fp.append(AllChem.RDKFingerprint(m, fpSize=fpSize))
        fp.append(AllChem.GetMorganFingerprintAsBitVect(m, radius=10, nBits=fpSize))

    # print(fp)
    X = np.array(list(fp))
    print("mean value of features =", np.mean(np.mean(X)))
    return X


def buildPharm(molecules):
    fdefName = '../data/external/MinimalFeatures.fdef'
    featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    # The fingerprints themselves are calculated using a signature (fingerprint) factory, which keeps track of all the parameters required to generate the pharmacophore:

    sigFactory = SigFactory(featFactory, minPointCount=2, maxPointCount=3)
    sigFactory.SetBins([(0, 2), (2, 5), (5, 8)])
    sigFactory.Init()
    sigFactory.GetSigSize()
    # The signature factory is now ready to be used to generate fingerprints, a task which is done using the rdkit.Chem.Pharm2D.Generate module:
    fp = []
    for m in molecules:
        # fp.append(AllChem.RDKFingerprint(m, fpSize=fpSize))
        fp.append(Generate.Gen2DFingerprint(m, sigFactory))
    return fp


def buildGobbi(molecules):
    fp = []
    for m in molecules:
        # fp.append(AllChem.RDKFingerprint(m, fpSize=fpSize))
        fp.append(Generate.Gen2DFingerprint(m, Gobbi_Pharm2D.factory))
    return fp


def buildFragment(molecules):
    fName = os.path.join(RDConfig.RDDataDir, 'FunctionalGroups.txt')
    fparams = FragmentCatalog.FragCatParams(1, 6, fName)
    print("Num functional groups:", fparams.GetNumFuncGroups())
    fcat = FragmentCatalog.FragCatalog(fparams)
    fcgen = FragmentCatalog.FragCatGenerator()
    for m in molecules:
        fcgen.AddFragsFromMol(m, fcat)
    print("Have", fcat.GetNumEntries(), "fragments in library")
    # Make fragment fingerprints
    fpgen = FragmentCatalog.FragFPGenerator()
    fp = []
    for m in molecules:
        # fp.append(AllChem.RDKFingerprint(m, fpSize=fpSize))
        fp.append(fpgen.GetFPForMol(m, fcat))
    #fp = np.array(list(fp))
    return fp, fcat
