from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression

MODEL_REGISTRY = {
    "mlp": MLPClassifier,
    "rf": RandomForestClassifier,
    "svm": SVC,
    "lr": LogisticRegression
}

def build_model(type : str):
    """
    """
    if type == "svm":
        return MODEL_REGISTRY[type](probabilites = True)
    return MODEL_REGISTRY[type]()