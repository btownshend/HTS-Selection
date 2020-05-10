#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from src.data import make_dataset
from src.features import build_features, build_targets
from src.models import split_data, train_model, predict_model, bitranking


# In[ ]:


# OPTIONAL: Load the "autoreload" extension so that code can change
get_ipython().run_line_magic('load_ext', 'autoreload')
# OPTIONAL: always reload modules so that as you change code in src, it gets loaded
get_ipython().run_line_magic('autoreload', '2')
seed=42


# In[5]:


# Load data
molecules=make_dataset.load()


# In[6]:


# Build the features
Xfull,fpcat=build_features.buildFragment(molecules)


# In[14]:


# Build the targets]
yfull,aptamers,tfull=build_targets.build_hitbyapt(molecules)
print("Loaded",len(aptamers),"aptamers and",len(tfull),"targets.")


# In[8]:


# Choose an aptamer to work on
# Setup training data to include only the ones for which we have definite hit/miss classification
apt=4
print("Testing",aptamers[apt])
tsel=[True if np.isfinite(yi) else False for yi in yfull[apt]]
y=yfull[apt][tsel]
X=[Xfull[i] for i in range(len(tsel)) if tsel[i]]
targets=[tfull[i] for i in range(len(tsel)) if tsel[i] ]


# In[9]:


# Create train/test sets
X_train, X_test, y_train, y_test, ind_train, ind_test = split_data.split(0.4,X,y,seed)
tgt_test=[targets[x] for x in ind_test]
tgt_train=[targets[x] for x in ind_train]
    


# In[10]:


# Train model
# Need to handle NaNs!
models=[train_model.train_rfc(X_train,y_train,seed,n_estimators=10)]


# In[11]:


# Test model
for model in models:
    print("------")
    print(model)
    print("Train:")
    yp_train=predict_model.predict(model,X_train,y_train,tgt_train)
    print("Test:")
    yp_test=predict_model.predict(model,X_test,y_test,tgt_test)
    print("LOO:")
    yp_loo=predict_model.predictLOO(model,X,y,targets)


# In[12]:


model.fit(X,y)
yfull = model.predict(Xfull)
for i in range(len(yfull)):
    if yfull[i]==1:
        if molecules[i].GetProp("NAME") not in targets:
            print(molecules[i].GetProp("NAME"),"may be an untested hit")    


# In[15]:


# Rank fragments based on observed activity
ntop=1
for apt in range(len(aptamers)):
    print("Testing",aptamers[apt])
    tsel=[True if np.isfinite(yi) else False for yi in yfull[apt]]
    y=yfull[apt][tsel]
    X=[Xfull[i] for i in range(len(tsel)) if tsel[i]]
    targets=[tfull[i] for i in range(len(tsel)) if tsel[i] ]

    frag=bitranking.getFragRanks(fpcat, X, y,ntop=ntop)
    print('Postives       : ',[targets[i] for i in range(len(targets)) if y[i]==1 and X[i].GetBit(int(frag))])
    print('False negatives: ',[targets[i] for i in range(len(targets)) if y[i]==1 and not X[i].GetBit(int(frag))])
    print('False positives: ',[targets[i] for i in range(len(targets)) if y[i]==0 and X[i].GetBit(int(frag))])


# In[21]:


m=molecules[648]
m2=Chem.AddHs(m)
print(CalcMolFormula(m2), Chem.MolToSmiles(m2))
print(Chem.Descriptors.ExactMolWt(m2))


# In[5]:


get_ipython().system('ls /tmp')

