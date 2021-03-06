from scipy import stats
import numpy as np 
import pandas as pd
import os
from tensorflow.keras.models import load_model
#import tensorflow as tf
import os

# This creates features out of the traces imported to it
# The traces need to be of a total length of 120 points.
def featureMaker(traces, steps):
    '''
    Input to this function is
    1: traces to calculate features
    2: Size of Steps between the times

    In This experiment we are going to attempt to extract some Features to see if this improves performance
    On the last experiment our dataset was of the for
    [A, B, C]
    A: Number of traces/samples ~10,000
    B: Number of TimeSteps 120
    C: Number of features 1

    For this experiment we want to reduce the number of timesteps while increasinG the number of Features
    A: Numver of traces/samples ~10,000
    B: Number of timeSteps 10
    C: Number of Features 3

    The features i would like to collect this time around are
    1: Mean over 12 points
    2: Standard Error over 12 points
    3: Derivative over 12 Points
    Function to create the above described features
    '''
    #some examples have na values get rid of
    traces = traces.dropna()
    
    #This need to be a 3 dimensional numpy array
    traces = np.asarray(traces)

    # Start my breaking up the traces into segments based on the
    # number of points in the trace.
    samples = traces.shape[0]
    timesteps = int(traces.shape[1] / steps)
    features = 4
    featureFrame = np.empty([samples, timesteps, features])
    rangeToCalc = np.arange(0, traces.shape[1]+1, steps).astype(int)

    for i in range(len(rangeToCalc)-1):
        meanFeat = traces[:,rangeToCalc[i]:rangeToCalc[i+1]].mean(axis=1)
        stdFeat = traces[:,rangeToCalc[i]:rangeToCalc[i+1]].std(axis=1)
        semFeat = stats.sem(traces[:,rangeToCalc[i]:rangeToCalc[i+1]], axis=1)
        derivFeat = np.mean(np.gradient(traces[:,rangeToCalc[i]:rangeToCalc[i+1]], axis=1), axis=1)

        featureFrame[:, i, 0] = meanFeat
        featureFrame[:, i, 1] = stdFeat
        featureFrame[:, i, 2] = semFeat
        featureFrame[:, i, 3] = derivFeat
    
    return featureFrame

# Function to create a specified number of features
def featureMaker2(traces, numWindows = 12):
    '''
    Input to this function is
    1: traces to calculate features
    2: Size of Steps between the times

    In This experiment we are going to attempt to extract some Features to see if this improves performance
    On the last experiment our dataset was of the for
    [A, B, C]
    A: Number of traces/samples ~10,000
    B: Number of TimeSteps 120
    C: Number of features 1

    For this experiment we want to reduce the number of timesteps while increasinG the number of Features
    A: Numver of traces/samples ~10,000
    B: Number of timeSteps 10
    C: Number of Features 3

    The features i would like to collect this time around are
    1: Mean over 12 points
    2: Standard Error over 12 points
    3: Derivative over 12 Points
    Function to create the above described features
    '''    
    #This need to be a 3 dimensional numpy array
    traces = np.asarray(traces)
    samples = traces.shape[0]
    dataPoints = traces.shape[1]

    timeSteps = np.linspace(start = int(0), stop = int(dataPoints), num=int(numWindows+1))
    timeSteps = np.floor(timeSteps).astype('int')
    # Start my breaking up the traces into segments based on the
    # number of points in the trace.
    features = 4
    featureFrame = np.empty([samples, len(timeSteps) - 1, features])

    for i in range(len(timeSteps) - 1):
        meanFeat = traces[:,timeSteps[i]:timeSteps[i+1]].mean(axis=1)
        stdFeat = traces[:,timeSteps[i]:timeSteps[i+1]].std(axis=1)
        semFeat = stats.sem(traces[:,timeSteps[i]:timeSteps[i+1]], axis=1)
        derivFeat = np.mean(np.gradient(traces[:,timeSteps[i]:timeSteps[i+1]], axis=1), axis=1)

        featureFrame[:, i, 0] = meanFeat
        featureFrame[:, i, 1] = stdFeat
        featureFrame[:, i, 2] = semFeat
        featureFrame[:, i, 3] = derivFeat

    return featureFrame

# Function to Run the models
def modelRunner(features, model):
    this_dir, this_filename = os.path.split(__file__)
    
    if model == 'aitc':
        DATA_PATH = os.path.join(this_dir, "models", "AITC_100uM_MAIN.h5")
    elif model == 'menthol':
        DATA_PATH = os.path.join(this_dir, "models", "Menthol_400uM_MAIN.h5")
    elif model == 'capsaicin':
        DATA_PATH = os.path.join(this_dir, "models", "Capsaicin_300nM_MAIN.h5")
    elif model == 'k40':
        DATA_PATH = os.path.join(this_dir, "models", "K_40mM_MAIN.h5")
    elif model == 'gfp':
        DATA_PATH = os.path.join(this_dir, "models", "gfp.h5")
    elif model == 'cy5':
        DATA_PATH = os.path.join(this_dir, "models", "cy5.h5")
    elif model == 'drop':
        DATA_PATH = os.path.join(this_dir, "models", "drop.h5")


    model = load_model(DATA_PATH)
    scores = model.predict(features)

    return scores

# Function to load the modesl    
def modelLoader(model):
    this_dir, this_filename = os.path.split(__file__)
    
    if model == 'aitc':
        DATA_PATH = os.path.join(this_dir, "models", "AITC_100uM_MAIN.h5")
    elif model == 'menthol':
        DATA_PATH = os.path.join(this_dir, "models", "Menthol_400uM_MAIN.h5")
    elif model == 'capsaicin':
        DATA_PATH = os.path.join(this_dir, "models", "Capsaicin_300nM_MAIN.h5")
    elif model == 'k40':
        DATA_PATH = os.path.join(this_dir, "models", "K_40mM_MAIN.h5")
    elif model == 'gfp':
        DATA_PATH = os.path.join(this_dir, "models", "gfp.h5")
    elif model == 'cy5':
        DATA_PATH = os.path.join(this_dir, "models", "cy5.h5")
    elif model == 'drop':
        DATA_PATH = os.path.join(this_dir, "models", "drop.h5")

    model = load_model(DATA_PATH)
    #scores = model.predict(features)

    return model
    
