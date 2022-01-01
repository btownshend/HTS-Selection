#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sklearn as sk


# In[ ]:


import sklearn as sk
import rdkit as rd
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import matplotlib
# matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt


# In[1]:


suppl = rd.Chem.SDMolSupplier('../../../data/HTBCFiles/ChemDivFull.sdf')


# In[1]:


# Extract only the 960 tested molecules
plates = ["CDIV%04d"%p for p in range(1,121,10)]
print(len(plates),plates)


# In[4]:


tested = [x for x in suppl if x.GetProp("BATCH_PLATE") in plates]
print(len(suppl),len(tested))


# In[5]:


for mol in tested:
    if mol is None: continue
    plate=mol.GetProp("BATCH_PLATE")
    plate=int(plate[5:])
    well=mol.GetProp("BATCH_WELL")
    name="%d%s"%(plate,well.replace('0',''))
    print(plate,well,name)
    mol.SetProp("NAME",name)


# In[6]:


for x in tested[0].GetPropNames(includePrivate=True,includeComputed=True):
    print(x,mol.GetProp(x))


# In[7]:


AllChem.Compute2DCoords(mol)


# In[8]:


print(mol)


# In[9]:


m2=rd.Chem.AddHs(mol)
AllChem.EmbedMolecule(m2)
AllChem.MMFFOptimizeMolecule(m2)
m2=rd.Chem.RemoveHs(m2)
img=Draw.MolsToGridImage([mol,m2],molsPerRow=2,subImgSize=(200,200),legends=["2d","3d"])
img.show()


# In[10]:


print(rd.Chem.MolToMolBlock(tested[1]))


# In[11]:


def getMol(mols,plate,well):
    mol=[x for x in mols if x.GetProp("BATCH_PLATE")=="CDIV%04d"%plate and x.GetProp("BATCH_WELL")==well]
    assert(len(mol)==1)
    return mol[0]


# In[12]:


# Compare some molecules which give the same response
m91d2=getMol(tested,91,'D02')
m91e2=getMol(tested,91,"E02")
img=Draw.MolsToGridImage([m91d2,m91e2],molsPerRow=2,subImgSize=(200,200),legends=["91E2","91D2"])
img.show()


# In[13]:


from rdkit.Chem import rdFMCS
res=rdFMCS.FindMCS([m91d2,m91e2])
res.smartsString


# In[14]:


res


# In[15]:


# Check chemical features
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import os
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
feats = [factory.GetFeaturesForMol(x) for x in [m91d2,m91e2]]


# In[16]:


len(feats[1])


# In[17]:


help(feats[0][0])


# In[18]:


for fs in feats:
    print("")
    for f in fs:
        print(f.GetFamily(),f.GetType(),f.GetAtomIds())


# In[20]:


# Generate pharmacophore fingerprints 
# from rdkit import Chem
# from rdkit.Chem import ChemicalFeatures
# fdefName = 'data/MinimalFeatures.fdef'
# featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)


# In[21]:


# from rdkit.Chem.Pharm2D.SigFactory import SigFactory
# sigFactory = SigFactory(featFactory,minPointCount=2,maxPointCount=3)
# sigFactory.SetBins([(0,2),(2,5),(5,8)])
# sigFactory.Init()
# sigFactory.GetSigSize()


# In[22]:


# Generate RDK fingerprints
fp=[]
for m in [m91d2, m91e2]:
    fp.append(AllChem.RDKFingerprint(m, fpSize=2048))
print(fp)


# In[23]:


help(fp[0])


# In[24]:


d=(~fp[0]&~fp[1]).GetNumOnBits()
a=(fp[0]&~fp[1]).GetNumOnBits()
b=(~fp[0]&fp[1]).GetNumOnBits()
c=(fp[0]&fp[1]).GetNumOnBits()
print(a,b,c,d)
rd.DataStructs.FingerprintSimilarity(fp[0],fp[1])


# In[25]:


# Setup fingerprints for all the molecules 
fp=[]
for m in tested:
    fp.append(AllChem.RDKFingerprint(m, fpSize=2048))
print(fp)


# In[26]:


# Set categorization of molecules
# Property, HIT, will be set to true iff we found an aptamer for the molecule
import csv

allhits=[]
with open('../../../data/Hits/hasaptamer.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            print(f'\tPlate {row[0]}, well {row[1]} is a hit.')
            line_count += 1
            allhits.append("%s%s"%(row[0],row[1]))
    print(f'Processed {line_count} hits.')
# Set HIT property
for t in tested:
    t.SetBoolProp("HIT",t.GetProp("NAME") in allhits)


# In[27]:


# Set categorization of molecules
# Property, HIT, will be set to true iff we found an aptamer for the molecule
import csv
fold=[]
targets=[]
with open('../../../data/Hits/fold.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print(f'Column names are {", ".join(row)}')
            aptamers=row[1:]
            line_count += 1
        else:
            line_count += 1
            fold.append([float(x) for x in row[1:]])
            targets.append(row[0])
    print(f'Processed {line_count} targets.')
print(targets)
# Set fold property
for t in tested:
    t.SetBoolProp("HIT",t.GetProp("NAME") in allhits)


# In[28]:


# Setup ML input as X, output as Y
import numpy as np
X = np.array(list(fp))
from sklearn.preprocessing import StandardScaler
#st = StandardScaler()
#X = st.fit_transform(X)
y=[t.GetBoolProp("HIT") for t in tested]


# In[29]:


from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=0)
print('Fraction hits: train: %.3f, test: %.3f'%(np.mean(y_train),np.mean(y_test)))


# In[30]:


from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
clf = RandomForestClassifier(n_estimators=15)
scores = cross_val_score(clf, X, y, cv=5)
np.mean(scores)


# In[31]:


from sklearn.ensemble import RandomForestClassifier

rfc = RandomForestClassifier(random_state=42)
rfc.fit(X_train, y_train)
y_predict=rfc.predict(X_test)
y_train_predict=rfc.predict(X_train)


# In[32]:


import matplotlib.pyplot as plt
import sklearn.metrics as sm

ax = plt.gca()
sm.roc_auc_score(y_test, y_predict)


# In[33]:


print(sm.classification_report(y_test, y_predict))


# In[34]:


print(sm.confusion_matrix(y_test, y_predict))


# In[35]:


print(sm.classification_report(y_train, y_train_predict))


# In[36]:


print(sm.confusion_matrix(y_train, y_train_predict))


# In[37]:



from sklearn.linear_model import LogisticRegression

lr = LogisticRegression(random_state=42)
lr.fit(X_train, y_train)
y_predict=lr.predict(X_test)
y_train_predict=lr.predict(X_train)


# In[38]:


print(sm.classification_report(y_train, y_train_predict))


# In[40]:


lr.coef_.shape
plt.hist(lr.coef_[0],100)
plt.show()
plt.interactive()

