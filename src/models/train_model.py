import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.svm import SVR, LinearSVR
from sklearn.tree import DecisionTreeRegressor


def crossval(mdl, X, y):
    scores = cross_val_score(mdl, X, y, cv=5, scoring='roc_auc')
    print(mdl)
    print("Cross-validation score = %.3f += %.3f" % (float(np.mean(scores)), float(np.std(scores))))


def train_rfc(X, y, random_state=42, n_estimators=5):
    mdl = RandomForestClassifier(random_state=random_state, n_estimators=n_estimators, class_weight="balanced")
    crossval(mdl, X, y)
    mdl.fit(X, y)
    return mdl


def train_lr(X, y, random_state=42):
    mdl = LogisticRegression(random_state=random_state, class_weight="balanced", solver='liblinear', multi_class='auto',
                             max_iter=1000)
    crossval(mdl, X, y)
    mdl.fit(X, y)
    # print("LR Coeff shape: ", mdl.coef_.shape)
    return mdl

def train_dtr(random_state=42):
    return DecisionTreeRegressor(random_state=random_state,max_leaf_nodes=10)

def train_svr(random_state=42):
    return SVR(gamma="auto")

def train_linearsvr(random_state=42):
    return LinearSVR(max_iter=10000,random_state=random_state)

def train_multi(tfun, X, y, random_state=42):
    mdl = [tfun(random_state=random_state) for _ in range(len(y[0]))]
    #crossval(mdl, X, y)
    for i in range(len(mdl)):
        mdl[i].fit(X, [y[j][i] for j in range(len(y))])
    return mdl
