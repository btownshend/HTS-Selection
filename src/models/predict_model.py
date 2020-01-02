import sklearn.metrics as sm


def predict(mdl, X, y_true, molecules):
    y = mdl.predict(X)
    print(sm.classification_report(y_true, y))
    print(sm.confusion_matrix(y_true, y))
    print("ROC_AUC Score = %.3f" % sm.roc_auc_score(y_true, y))
    for i in range(len(X)):
        if y[i]!=y_true[i]:
            print("Molecule %3d %6.6s: predict=%.1f, actual=%.1f"%(i,molecules[i].GetProp("NAME"),y[i],y_true[i]))
    return y
