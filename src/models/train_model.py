import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, cross_validate
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

def train_multi(tfun, X, y, random_state=42):
    mdl=[tfun(X,yi,random_state=random_state) for yi in y]
    for i in range(len(mdl)):
        mdl[i].fit(X, y[i])
    return mdl

def cv_model(mdl,X,y,scoring='neg_mean_squared_error'):
    print(mdl)
    cvres = cross_validate(mdl, X, y, cv=3,scoring=scoring,return_train_score=True)
    print(cvres)
    scores=cvres['test_score']
    tscores=cvres['train_score']
    print("Cross-validation scores: test=%.3f+-%.3f, train=%.3f+-%.3f"%(np.asscalar(np.mean(scores)),np.asscalar(np.std(scores)),np.asscalar(np.mean(tscores)),np.asscalar(np.std(tscores))))

def cv_model_sfm(mdl,X,y,max_features=3,scoring='neg_mean_squared_error'):
    smdl = SelectFromModel(mdl,max_features=max_features)
    smdl.fit(X,y)
    features=smdl.get_support(indices=True)
    if len(features)>0:
        print("features:", features)
        Xt=smdl.transform(X)
        cv_model(mdl,Xt,y,scoring=scoring)
    else:
        print("No features retainec")
