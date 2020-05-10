import numpy as np
import sklearn.metrics as sm
from sklearn.model_selection import LeaveOneOut


def predict(mdl, X, y_true, molecules):
    y = mdl.predict(X)
    print(sm.classification_report(y_true, y))
    print(sm.confusion_matrix(y_true, y))
    print("ROC_AUC Score = %.3f" % sm.roc_auc_score(y_true, y))
    for i in range(len(X)):
        if y[i] != y_true[i]:
            print("%6.6s: predict=%.1f, actual=%.1f" % (molecules[i], y[i], y_true[i]))
    return y


def predict_reg(mdl, X, y_true):
    y = np.transpose([mdl[i].predict(X) for i in range(len(mdl))])
    score = [mdl[i].score(X, y_true[i]) for i in range(len(mdl))]
    print("score=", score)
    print("mean(score)=%.3f" % np.mean(score))
    return y


def predictLOO(mdl, X, y_true, molecules):
    loo = LeaveOneOut()
    y = [0 for _ in y_true]
    print("Running Leave-out validation on", len(molecules), "molecules...")
    for train_index, test_index in loo.split(X):
        if test_index % 20 == 0:
            print("%d..." % test_index, end="")
        X_train, X_test = [X[i] for i in train_index], [X[i] for i in test_index]
        y_train, y_test = [y_true[i] for i in train_index], [y_true[i] for i in test_index]
        mdl.fit(X_train, y_train)
        p = mdl.predict(X_test)
        y[test_index[0]] = p[0]
        if p != y_test[0]:
            print("\n%6.6s: predict=%.1f, actual=%.1f" % (molecules[test_index[0]], p, y_test[0]))
    print('done')
    print(sm.classification_report(y_true, y))
    print(sm.confusion_matrix(y_true, y))
    print("ROC_AUC Score = %.3f" % sm.roc_auc_score(y_true, y))
    return y
