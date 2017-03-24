import os
os.getcwd()
os.chdir('/home/tsugrad/Documents/gerrymandering/')
from setup_stuff import * 
from mhstuff import *


#State and precincts 
stateSHORT='PA'
os.chdir('Pennsylvania')
g = package_vtds("precinct/precinct.shp")

# VTD stats
blockstats = pd.read_csv("vtdstats.csv")
blockstats = blockstats.drop('Unnamed: 0', 1)
blockstats = blockstats.set_index(blockstats.VTD)

# number of districts and VTDS
cdtable = pd.read_csv('../cdbystate1.txt', '\t')
ndistricts = int(cdtable[cdtable['STATE']==stateSHORT].CD)
nvtd = len(blockstats.VTD)
colors = colorDict(ndistricts)


#Read adjacency frame
adjacencyFrame = pd.read_csv('PRECINCTconnections.csv')
adjacencyFrame = adjacencyFrame.drop('Unnamed: 0', 1)
adjacencyFrame.columns = ['low', 'high', 'length']
#adjacencyFrame.low  = [x[5:] for x in adjacencyFrame.low]
#adjacencyFrame.high = [x[5:] for x in adjacencyFrame.high]


#mh_initialize_globals(nvtd, ndistricts, adjacencyFrame, blockstats, g)

numwalkers = 1
numsteps = 20
numplots = 20
for walker in range(numwalkers): 
    starting_state = contiguousStart()
    starting_state.columns = ['key', 'value']
    foldername = 'rangledangle' + str(walker) + '/'
    runningState = (starting_state.copy(), 1)
    os.mkdir(foldername)
    color_these_states(g, [runningState], foldername, 0)

    for step in range(numplots):
        runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
        color_these_states(g, [runningState], foldername, step+1)
        pd.DataFrame(runningState[0]).to_csv(foldername+"state%d.csv"%(step+1))


def semigoodness(state):
    #Haves
        #contiguousness
        #evenness of population
        #efficiency
        #bizarreness
    #Needs
        #Compactness
    
    stpops  = [population(state, i) for i in range(ndistricts)]
    stbiz   = [bizarreness(state, i) for i in range(ndistricts)]
    
    modTotalVar = sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in stpops])/(2*(1-float(1)/ndistricts))
    
    return -100*modTotalVar - 10*np.nansum(stbiz)


for step in range(numplots):
    runningState = MH(runningState[0], numsteps, neighbor, semigoodness, switchDistrict)
    color_these_states(g, [runningState], foldername, step+1)
    pd.DataFrame(runningState[0]).to_csv(foldername+"state%d.csv"%(step+1))





