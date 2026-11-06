from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis

MODEL_REGISTRY = {
    "mlp": MLPClassifier,
    "rf": RandomForestClassifier,
    "svm": SVC,
    "lr": LogisticRegression,
    "knn": KNeighborsClassifier,
    "lda": LinearDiscriminantAnalysis,
    "qda": QuadraticDiscriminantAnalysis
}

def build_model(type : str):
    """
    """
    if type == "svm":
        return MODEL_REGISTRY[type](probability = True)
    return MODEL_REGISTRY[type]()