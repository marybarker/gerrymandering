import os
os.getcwd()
os.chdir('Documents/tarleton_misc/gerrymandering/')
from mhstuff import *
from setup_stuff import * 

#State and precincts 
stateSHORT='PA'
os.chdir('Pennsylvania')
g = package_vtds("precinct/precinct.shp")

# VTD stats
blockstats = pd.read_csv("precinctstats.csv")
blockstats = blockstats.drop('Unnamed: 0', 1)
blockstats.set_index(blockstats.VTD)

# number of districts and VTDS
cdtable = pd.read_csv('../cdbystate1.txt', '\t')
ndistricts = int(cdtable[cdtable['STATE']==stateSHORT].CD)
nvtd = len(blockstats.VTD)

#Read adjacency frame
adjacencyFrame = pd.read_csv('PRECINCTconnections.csv')
adjacencyFrame = adjacencyFrame.drop('Unnamed: 0', 1)
adjacencyFrame.columns = ['low', 'high', 'length']
#adjacencyFrame.low  = [x[5:] for x in adjacencyFrame.low]
#adjacencyFrame.high = [x[5:] for x in adjacencyFrame.high]


numwalkers = 10
numsteps = 100
numplots = 10
for walker in range(numwalkers): 
    starting_state = contiguousStart(blockstats)
    starting_state.columns = ['key', 'value']
    foldername = 'rangledangle' + str(walker) + '/'
    runningState = (starting_state.copy(), 1)
    os.mkdir(foldername)
    plot_stuff.color_these_states(g, [runningState], foldername, 0)

    for step in range(numplots):
        runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
        plot_stuff.color_these_states(g, [runningState], foldername, step+1)
        pd.DataFrame(runningState[0]).to_csv(foldername+"state%d.csv"%(step+1))

