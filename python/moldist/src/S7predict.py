#!/usr/bin/env python
# coding: utf-8

import numpy as np
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from src.features import load_singles
from src.features.build_features import buildFragment
from src.models import predict_model
from dotenv import load_dotenv, find_dotenv
from sacred import Experiment
from sacred.observers import FileStorageObserver

ex = Experiment('S7predict')
ex.observers.append(FileStorageObserver('my_runs'))

@ex.config
def my_config():
    foldthresh = 2.0  # Threshold for fold change to classify as a hit
    minhits = 2  # Minimum number of hits to run analysis of an aptamer
    seed = 2  # Random seed for reproducibility
    model = RandomForestClassifier(random_state=seed, n_estimators=100, class_weight="balanced")


@ex.automain
def my_main(foldthresh,minhits,model):
    # Load environment
    load_dotenv(find_dotenv())

    # Classifier to use for experiments
    print(f"Modelling using {model}")

    # Load fold change data exported from matlab
    # fold[napt][2][ncompound] - fold change confidence interval
    # aptamers[napt] - aptamer names (naseq)
    # targets[ncompound][2] - [name,mol] for each compound
    [fold, aptamers, targets] = load_singles.load_singles()
    print(f"Loaded {len(aptamers)} aptamers, {len(targets)} targets.")

    # List comppounds and their SMILES
    # for t in targets:
    #     print(t[0], AllChem.MolToSmiles(t[1]))

    # Build the features
    Xfull, fpcat = buildFragment([t[1] for t in targets])

    # Build the targets:
    # yfull[ncompound][napt] is 1 if fold change is definitely >foldthresh, 0 if it is definitely less,
    #                        and nan if indeterminate
    yfull = [[0 if i[1] < foldthresh else 1 if i[0] > foldthresh else np.nan for i in zip(ff[0], ff[1])] for ff in fold]
    # Transpose so we have yfull[napt][ncompound]
    yfull = np.array(yfull).T
    print(yfull.shape)

    # Loop over aptamers
    yest=[]
    cm=np.zeros((2,2))
    for apt in range(len(aptamers)):
        # Check for minimum number of hits
        npos = sum([yi == 1 for yi in yfull[apt]])
        if npos < minhits:
            print(f"Skipping {aptamers[apt]} since it has only {npos} hits")
            yest.append([])
            continue

        # Setup training data to include only the ones for which we have definite hit/miss classification
        tsel = [True if np.isfinite(yi) else False for yi in yfull[apt]]
        y = yfull[apt][tsel]
        X = [Xfull[i] for i in range(len(tsel)) if tsel[i]]
        tnames = [f"{targets[i][0]}[{fold[i][0][apt]:.2f},{fold[i][1][apt]:.2f}]" for i in range(len(tsel)) if tsel[i]]
        print(f"Testing {aptamers[apt]} against {len(tnames)} targets with {sum(y):.0f} positives")
        # Run leave-out-one prediction of the hits
        yp_loo = predict_model.predictLOO(model, X, y, tnames)
        yest.append(yp_loo)
        cm=cm+confusion_matrix(y, yp_loo)
        print(cm)

#%%
