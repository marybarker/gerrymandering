""" SAMPLE FILE FOR RUNNING MCMC """
""" change into the directory for that state (e.g. Pennsylvania/ folder) """
""" and run to get all functions loaded and a basic MCMC run on a random initial redistricting plan"""

import os
import time
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import math
from osgeo import ogr

SHAPEFILE_FILE='../setup_stuff.py'
MCMC_FILE='../cleanitallup.py'
THIS_STATE_SETUP_FILE='setup.py'

execfile(SHAPEFILE_FILE) #load all functions for dealing with shapefiles
execfile(MCMC_FILE) # load MCMC functions
execfile(THIS_STATE_SETUP_FILE) # setup datastructure for current state based on shapefiles


ndistricts=13 # how many congressional districts wanted 

# create random initial state 
initial_state=contiguousStart()

# save as a png file to see 
color_this_state(g,initial_state,'test_example.png'.1)


for i in range(1, 100):
    runningState = (currentNCstate.copy(), 0)
    updateGlobals(runningState[0])
    
    for j in range(300):
        if j > 200:
            goodnessWeights= [200,200,200]
        runningState = MH(runningState[0], 500, neighbor, goodness, switchDistrict)
        
    runningState[0].to_csv(foldername+'data/state%d.csv'%i, index=False)
    metrics.to_csv(foldername+'data/metrics%d.csv'%i, index=False)
    color_these_states(g, [runningState], foldername+'data/pictures/', i, 0.1)
    print 'finished with step %d. '%i, dfEquiv(runningState[0], currentNCstate)
