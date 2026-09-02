

def predict_pocket(model, X):
    """
    """
    Y = model.predict_proba(X)
    class_names = model.classes_
    output = [{cls: p for cls, p in zip(class_names, row)} for row in Y]
    return output

