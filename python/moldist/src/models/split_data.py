from sklearn.model_selection import train_test_split
import numpy as np


def split(test_size, X, y, random_state=0):
    X_train, X_test, y_train, y_test, ind_train, ind_test = \
        train_test_split(X, y, range(len(X)), test_size=test_size, random_state=random_state, stratify=y)
    print('Fraction hits: train: %.3f, test: %.3f' % (float(np.mean(y_train)), float(np.mean(y_test))))
    return X_train, X_test, y_train, y_test, ind_train, ind_test
