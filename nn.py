import keras
import numpy as np
from collections import defaultdict

def myfunction(modelfile, X):
    '''
    Predicts Curvature (1/inch) from [Deflection (inch), Elevation (inch)]
    using a trained Keras neural network located in modelfile and data X.
    '''
    # NN was trained on preprocessed data, so input needs same processing.
    assert X.shape[1] == 2
    prep = defaultdict(dict,
            {'Deflection': {'mean': 0.0, 'std': 1.0},
             'Elevation': {'mean': 0.0, 'std': 100.0},
             'Curvature': {'mean': 0.0, 'std': 0.00032503688027413094}})
    X[:,0] = (X[:,0] - prep['Deflection']['mean']) / prep['Deflection']['std']
    X[:,1] = (X[:,1] - prep['Elevation']['mean']) / prep['Elevation']['std']

    # Load model and predict.
    model = keras.models.load_model(modelfile)
    yhat = model.predict(X)
    
    # Output needs to be processed as well.
    yhat = (yhat * prep['Curvature']['std']) + prep['Curvature']['mean']
    return yhat

# Example usage.
exampledata = np.array([[1., 60.], [-1., 60.]])
modelfile = './model_saved.h5'
yhat = myfunction(modelfile, exampledata)
print(yhat)
