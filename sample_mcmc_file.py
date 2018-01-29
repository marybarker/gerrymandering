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
color_this_state(g,initial_state,'test_example.png',.1)

updateGlobals(initial_state)
initial_state_metrics=metrics.copy()

# run MCMC
new_state = MH(initial_state, 500, neighbor, goodness, switchDistrict)

# get new state of things and compare
updated_state_metrics=metrics.copy()

color_these_states(g, [(initial_state,0), new_state], '', 0, 0.1)



