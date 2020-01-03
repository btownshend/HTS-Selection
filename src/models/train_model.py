import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score


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
