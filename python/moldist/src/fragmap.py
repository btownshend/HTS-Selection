# %%
import os
import csv
import dotenv
import numpy as np

from src.data import make_dataset
from src.features import build_features, build_targets
from src.models import bitranking

project_dir = os.path.join(os.path.dirname(__file__), os.pardir)
dotenv_path = os.path.join(project_dir, '.env')
dotenv.load_dotenv(dotenv_path, verbose=True)

# %%
# Load data
molecules = make_dataset.load()

# %%
# Build the features using fragments from the structures
Xfull, fpcat = build_features.buildFragment(molecules)

# %%
# Build the targets]
yfull, aptamers, tfull = build_targets.build_hitbyapt(molecules)
print("Loaded", len(aptamers), "aptamers and", len(tfull), "targets.")

# %%
# Rank fragments based on observed activity
ntop = 20
ypred=[]  # Hit predictions
for apt in range(len(aptamers)):
    print("Testing", aptamers[apt])
    tsel = [True if np.isfinite(yi) else False for yi in yfull[apt]]
    y = yfull[apt][tsel]
    X = [Xfull[i] for i in range(len(tsel)) if tsel[i]]
    targets = [tfull[i] for i in range(len(tsel)) if tsel[i]]

    frag = bitranking.getFragRanks(fpcat, X, y, ntop=ntop)
    print("frag=",frag)
    pred=[1 if xx.GetBit(int(frag)) else 0 for xx in Xfull]
    ypred.append(pred)
    print('Postives       : ', [targets[i] for i in range(len(targets)) if y[i] == 1 and X[i].GetBit(int(frag))])
    print('False negatives: ', [targets[i] for i in range(len(targets)) if y[i] == 1 and not X[i].GetBit(int(frag))])
    print('False positives: ', [targets[i] for i in range(len(targets)) if y[i] == 0 and X[i].GetBit(int(frag))])
    print("ypred has ",len([i for i, val in enumerate(pred) if val==1] )," positives")
# %%
# Save ypred
with open('preds.csv','w') as csv_file:  # From matlab:   TRPSummary.dumphits()
    csvwriter = csv.writer(csv_file, delimiter=',')
    csvwriter.writerow(['Target']+aptamers)
    for i in range(len(tfull)):
        csvwriter.writerow(tfull[i:i+1]+[yp[i] for yp in ypred])
