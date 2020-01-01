#!/usr/bin/env python
# coding: utf-8

# In[2]:


import sklearn as sk
import rdkit as rd
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


# In[3]:


help(sk)


# In[4]:


help(rd.Chem)


# In[5]:


suppl = rd.Chem.SDMolSupplier('../data/raw/ChemDivFull.sdf')


# In[6]:


print(suppl)


# In[7]:


# Extract only the 960 tested molecules
plates = ["CDIV%04d"%p for p in range(1,121,10)]
print(len(plates),plates)


# In[8]:


tested = [x for x in suppl if x.GetProp("BATCH_PLATE") in plates]
print(len(suppl),len(tested))


# In[9]:


for mol in tested:
    if mol is None: continue
    if mol.GetProp("BATCH_WELL")=="A02":
        print(mol.GetNumAtoms(),mol.GetProp("_Name"),mol.GetProp("BATCH_PLATE"),mol.GetProp("BATCH_WELL"))


# In[10]:


for x in mol.GetPropNames(includePrivate=True,includeComputed=True):
    print(x,mol.GetProp(x))


# In[11]:


AllChem.Compute2DCoords(mol)


# In[12]:


print(mol)


# In[42]:


m2=rd.Chem.AddHs(mol)
AllChem.EmbedMolecule(m2)
AllChem.MMFFOptimizeMolecule(m2)
m2=rd.Chem.RemoveHs(m2)
img=Draw.MolsToGridImage([mol,m2],molsPerRow=2,subImgSize=(200,200),legends=["2d","3d"])
img.show()


# In[14]:


print(rd.Chem.MolToMolBlock(tested[1]))


# In[30]:


def getMol(mols,plate,well):
    mol=[x for x in mols if x.GetProp("BATCH_PLATE")=="CDIV%04d"%plate and x.GetProp("BATCH_WELL")==well]
    assert(len(mol)==1)
    return mol[0]


# In[38]:


# Compare some molecules which give the same response
m91d2=getMol(tested,91,'D02')
m91e2=getMol(tested,91,"E02")
img=Draw.MolsToGridImage([m91d2,m91e2],molsPerRow=2,subImgSize=(200,200),legends=["91E2","91D2"])
img.show()


# In[44]:


from rdkit.Chem import rdFMCS
res=rdFMCS.FindMCS([m91d2,m91e2])
res.smartsString


# In[45]:


res


# In[47]:


# Check chemical features
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import os
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
feats = [factory.GetFeaturesForMol(x) for x in [m91d2,m91e2]]


# In[51]:


len(feats[1])


# In[53]:


help(feats[0][0])


# In[61]:


for fs in feats:
    print("")
    for f in fs:
        print(f.GetFamily(),f.GetType(),f.GetAtomIds())


# In[62]:


# Generate pharmacophore fingerprints 
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
fdefName = 'data/MinimalFeatures.fdef'
featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)


# In[ ]:


from rdkit.Chem.Pharm2D.SigFactory import SigFactory
sigFactory = SigFactory(featFactory,minPointCount=2,maxPointCount=3)
sigFactory.SetBins([(0,2),(2,5),(5,8)])
sigFactory.Init()
sigFactory.GetSigSize()

