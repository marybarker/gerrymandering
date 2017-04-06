import os
import time
os.getcwd()
#os.chdir('/Users/marybarker/Documents/tarleton_misc/gerrymandering/')
os.chdir('/home/odin/Documents/gerrymandering/gerrymandering/')
#from setup_stuff import * 
#from mhstuff import *


#State and precincts
stateSHORT='PA'
os.chdir('Pennsylvania')
g = package_vtds("precinct/precinct.shp")

# VTD stats
blockstats = pd.read_csv("vtdstats.csv")
blockstats = blockstats.drop('Unnamed: 0', 1)
blockstats = blockstats.set_index(blockstats.VTD)
blockstats.rename(columns={'population':'POP100'}, inplace=True)
totalpopulation = sum(blockstats.POP100)

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


foldername = "farster3/"
os.mkdir(foldername)
numwalkers = 5
walker=1
numsteps = 10
numplots = 100

starting_state = pd.read_csv('./farster3/state0.csv')
starting_state = pd.read_csv('./farster3/state5100.csv')
#starting_state = starting_state.drop('Unnamed: 0', 1)
#starting_state = contiguousStart2()
#starting_state.columns = ['key', 'value']

temp = dict(zip(starting_state.key, starting_state.value))
lowdists  = adjacencyFrame.low.replace(temp)
highdists = adjacencyFrame.high.replace(temp)
isSame = lowdists==highdists
adjacencyFrame['isSame'] = isSame

runningState = (starting_state.copy(), 1)
color_these_states(g, [runningState], foldername, 5100)
pd.DataFrame(runningState[0]).to_csv(foldername+"state0.csv", index = False)


numplots = 1000

runtimes = np.array([float(0)]*numplots)

for step in range(20, numplots):
    
    starttime = time.clock()
    
    runningState = MH(runningState[0], numsteps, ambitiousNeighbor, goodness, switchDistrict)
    #color_these_states(g, [runningState], foldername, step+1)
    pd.DataFrame(runningState[0]).to_csv(foldername+"state%d.csv"%(step+1), index = False)
    
    end = time.clock()
    
    runtimes[step] = end - starttime
    
    print("Written state %d of %d.  Runtime: %f"%(step, numplots, runtimes[step]))



os.mkdir("startingPoints")
for i in range(100):
    tempStart = contiguousStart2()
    tempStart.to_csv("./startingPoints/start%d.csv"%i, index = False)
    print("Finished writing %d.\n"%i)






