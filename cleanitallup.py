import os
import time
import random
import numpy as np
import pandas as pd

#os.chdir('/Users/marybarker/Documents/tarleton_misc/gerrymandering/Pennsylvania')
#os.chdir('/home/odin/Documents/gerrymandering/gerrymandering/Pennsylvania')
stateSHORT = 'PA'

blockstats = pd.read_csv("vtdstats.csv")
blockstats.rename(columns = {"POP100":"population"}, inplace = True)
blockstats = blockstats.drop('Unnamed: 0', 1)
blockstats = blockstats.set_index(blockstats.VTD)
totalpopulation = sum(blockstats.population)

cdtable = pd.read_csv('../cdbystate1.txt', '\t')
ndistricts = int(cdtable[cdtable['STATE']==stateSHORT].CD)
nvtd = len(blockstats.VTD)

adjacencyFrame = pd.read_csv('PRECINCTconnections.csv')
adjacencyFrame = adjacencyFrame.drop('Unnamed: 0', 1)
adjacencyFrame.columns = ['low', 'high', 'length']
metrics = {}

foldername = "slambp2/"
#os.mkdir(foldername)

numwalkers = 5
walker=1
numsteps = 100
numsaves = 1000
numplots = 10
startingPoint=0

starting_state = pd.read_csv('./startingPoints/start0.csv')
#starting_state = pd.read_csv('/Users/marybarker/Downloads/starting_states/start0.csv')
del starting_state['Unnamed: 0']

runningState = (starting_state.copy(), 1)

runtimes = np.array([float(0)]*numsaves)

for i in range(100, numsaves):
    
    #    starttime = time.clock()
    updateGlobals(runningState[0])
    runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
    runningState[0].to_csv(foldername+"state%d_save%d.csv"%(startingPoint, i + 1), index = False)
    #    end = time.clock()
    
    #    runtimes[i] = end - starttime
    
    print("Written state %d of %d."%(i+1, numsaves))

for i in range(100, numsaves):
    tempstate = pd.read_csv(foldername + "state%d_save%d.csv"%(startingPoint, i + 1))
    color_these_states(g, [(tempstate, 0)], foldername + "figures/", i+1)
    print("Colored %d of these states!"%(i+1))


##################################################################################
def MH(start, steps, neighbor, goodness, moveprob):
    #  object starting state   |         |
    #         integer steps to be taken for M-H algorithm.
    #                function returning a neighbor of current state
    #                          function for determining goodness.
    #                                    function which takes goodnesses and returns probabilities.
    global adjacencyFrame, metrics
    current = start.copy()
    best_state = start.copy()
    current_goodness = goodness(metrics)
    best_goodness = current_goodness
    better_hops = 0
    worse_hops = 0
    stays = 0
    for i in range(steps):
        possible = neighbor(current)
        possible_goodness = goodness(possible[2])
        if best_goodness < possible_goodness:
            best_state = possible[0].copy()
            best_goodness = possible_goodness
            best_metrics = possible[2].copy()
        if random.random() < moveprob(current_goodness, possible_goodness):
            if current_goodness < possible_goodness :
                better_hops += 1
            else:
                worse_hops += 1
            current = possible[0].copy()
            current_goodness = possible_goodness
            adjacencyFrame.update(possible[1])
            adjacencyFrame.lowdist  = adjacencyFrame.lowdist.astype(int)
            adjacencyFrame.highdist = adjacencyFrame.highdist.astype(int)
            metrics = possible[2].copy()
        else:
            stays += 1
    return((best_state, best_goodness, better_hops, worse_hops, stays))


def neighbor(state):
    
    #stConts = [contiguousness(runningState[0], i) for i in range(ndistricts)]
    #stPops  = [    population(runningState[0], i) for i in range(ndistricts)]
    #stBiz   = [   bizarreness(runningState[0], i) for i in range(ndistricts)]
    #stPerim = [     perimeter(runningState[0], i) for i in range(ndistricts)]
    #stArea  = [      distArea(runningState[0], i) for i in range(ndistricts)]
    newstate = state.copy()
    newmetrics = metrics.copy()

    missingdist = set.difference(set(range(ndistricts)), set(newstate['value']))
    #If we've blobbed out some districts, we wants to behave differently
    
    if len(missingdist) == 0:
        switchedge = np.random.choice(adjacencyFrame.index[-(adjacencyFrame.isSame == 1)])

        lownode      = adjacencyFrame.low[switchedge]
        highnode     = adjacencyFrame.high[switchedge]
        templowdist  = adjacencyFrame.lowdist[switchedge]
        temphighdist = adjacencyFrame.highdist[switchedge]
        #Randomly choose an adjacency.  Find the low node and high node for that adjacency.

        if random.random() < 0.5:
            
            #switch low node stuff to high node's district
            
            switchTo = (newstate[newstate.key == highnode].value).item()
            
            proposedChanges = adjacencyFrame[(adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)]
            #want proposedChanges to be a slice of adjacencyFrame where the values could be changing.
            
            newstate.value[newstate.key == lownode] = switchTo
            proposedChanges.lowdist[proposedChanges.low == lownode] = switchTo
            proposedChanges.highdist[proposedChanges.high == lownode] = switchTo
            proposedChanges.isSame = proposedChanges.lowdist == proposedChanges.highdist
            #change values in the state as well as the proposedChanges
            
            #update contiguousness
            newmetrics['contiguousness'][templowdist]  = contiguousness(newstate, templowdist)
            newmetrics['contiguousness'][temphighdist] = contiguousness(newstate, temphighdist)
            
            #change population
            popchange = blockstats.population[lownode]
            newmetrics['population'][templowdist]  -= popchange
            newmetrics['population'][temphighdist] += popchange
            
            #change bizarreness
            newmetrics['perimeter'][templowdist]  = perimeter(newstate, templowdist)
            newmetrics['perimeter'][temphighdist] = perimeter(newstate, temphighdist)
            
            areachange = blockstats.ALAND[lownode] + blockstats.AWATER[lownode]
            newmetrics['area'][templowdist] -= areachange
            newmetrics['area'][temphighdist] += areachange
            
            newmetrics['bizarreness'][templowdist] = bizarreness(newmetrics['area'][templowdist], \
                                                                  newmetrics['perimeter'][templowdist])
            newmetrics['bizarreness'][temphighdist] = bizarreness(newmetrics['area'][temphighdist], \
                                                                   newmetrics['perimeter'][temphighdist])
            
        else:
            
            #switch high node stuff to low node's district
            
            switchTo = (newstate[newstate.key == lownode].value).item()
            #switch to low node
            
            proposedChanges = adjacencyFrame[(adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)]
            #want proposedChanges to be a slice of adjacencyFrame where the values could be changing.
            
            newstate.value[newstate.key == highnode] = switchTo
            proposedChanges.lowdist[proposedChanges.low == highnode] = switchTo
            proposedChanges.highdist[proposedChanges.high == highnode] = switchTo
            proposedChanges.isSame = proposedChanges.lowdist == proposedChanges.highdist
            #change values in the state as well as the proposedChanges
            
            #update contiguousness
            newmetrics['contiguousness'][temphighdist]  = contiguousness(newstate, temphighdist)
            newmetrics['contiguousness'][templowdist] = contiguousness(newstate, templowdist)
            
            #change population
            popchange = blockstats.population[highnode]
            newmetrics['population'][temphighdist]  -= popchange
            newmetrics['population'][templowdist] += popchange
            
            #change bizarreness
            newmetrics['perimeter'][temphighdist]  = perimeter(newstate, temphighdist)
            newmetrics['perimeter'][templowdist] = perimeter(newstate, templowdist)
            
            areachange = blockstats.ALAND[highnode] + blockstats.AWATER[highnode]
            newmetrics['area'][temphighdist] -= areachange
            newmetrics['area'][templowdist] += areachange
            
            newmetrics['bizarreness'][temphighdist] = bizarreness(newmetrics['area'][temphighdist], \
                                                                  newmetrics['perimeter'][temphighdist])
            newmetrics['bizarreness'][templowdist] = bizarreness(newmetrics['area'][templowdist], \
                                                                   newmetrics['perimeter'][templowdist])
    else:
        #If there are some districts missing, 
        changenode = newstate.key.sample(1)
        olddist = newstate.value[changenode]
        newdist = list(missingdist)[0]
        newstate.value[newstate.key == changenode] = newdist
        #We want to select one randomly, and make it one of the missing districts
        proposedChanges = adjacencyFrame.loc[(adjacencyFrame.low == changenode) | \
                              (adjacencyFrame.high == changenode)]
        proposedChanges.lowdist[proposedChanges.low == changenode] = newdist
        proposedChanges.highdist[proposedChanges.high == changenode] = newdist
        proposedChanges.isSame = False
        # And none of its adjacencies match anymore.
        
        #change contiguousness
        newmetrics['contiguousness'][olddist] = contiguousness(newstate, olddist)
        
        #change population
        popchange = blockstats.population[changenode]
        newmetrics['population'][olddist] -= popchange
        newmetrics['population'][newdist] += popchange
        
        #change bizarreness
        newmetrics['perimeter'][olddist] = perimeter(newstate, olddist)
        newmetrics['perimeter'][newdist] = perimeter(newstate, newdist)
        
        areachange = blockstats.ALAND[changenode] + blockstats.AWATER[changenode]
        newmetrics['area'][olddist] -= areachange
        newmetrics['area'][newdist] += areachange
        
        newmetrics['bizarreness'][olddist] = bizarreness(newmetrics['area'][olddist], \
                                                              newmetrics['perimeter'][olddist])
        newmetrics['bizarreness'][newdist] = bizarreness(newmetrics['area'][newdist], \
                                                              newmetrics['perimeter'][newdist])
    
    
    return (newstate, proposedChanges, newmetrics)


def contiguousness(state, district):
    
    regions = 0
    regionlist = list(state.key[state.value == district])
    if len(regionlist) == 0:
        return 1
    
    subframe = adjacencyFrame.loc[adjacencyFrame.low.isin(regionlist) & adjacencyFrame.high.isin(regionlist)]
    subedges = subframe[subframe.length != 0][['low','high']]
    
    while len(regionlist) > 0:
        regions += 1
        currentregion = set()
        addons = {regionlist[0]}
        while len(addons) > 0:
            currentregion = currentregion.union(addons)
            subsubedges = subedges.loc[subedges.low.isin(currentregion) | subedges.high.isin(currentregion)]
            if(not subsubedges.empty):
                addons = set(subsubedges['low']).union(set(subsubedges['high'])) - currentregion
            else:
                addons = set()
        regionlist = [x for x in regionlist if x not in currentregion]
    return regions

def perimeter(state, district):
    regionlist = list(state.key[state.value == district])
    return sum(adjacencyFrame.length[adjacencyFrame.low.isin(regionlist) != adjacencyFrame.high.isin(regionlist)])

def interiorPerimeter(state, district):
    regionlist = list(state.key[state.value == district])
    return sum(adjacencyFrame.length[adjacencyFrame.low.isin(regionlist) & adjacencyFrame.high.isin(regionlist)])

def distArea(state, district):
    regionlist = list(state.key[state.value == district])
    return sum(blockstats.ALAND[blockstats.VTD.isin(regionlist)]) + \
           sum(blockstats.AWATER[blockstats.VTD.isin(regionlist)])

def population(state, district):
    return sum(blockstats.population[blockstats.VTD.isin(list(state.key[state.value == district]))])

def efficiency(state, district):
    #returns difference in percentage of votes wasted.  Negative values benefit R.
    subframe = blockstats.loc[blockstats.VTD.isin(list(state.key[state.value == district]))]
    rvotes = sum(subframe['repvotes'])
    dvotes = sum(subframe['demvotes'])
    
    if rvotes > dvotes:
        wastedR = max(rvotes, dvotes) - 0.5
        wastedD = min(rvotes,dvotes)
    else:
        wastedD = max(rvotes, dvotes) - 0.5
        wastedR = min(rvotes,dvotes)
    
    return wastedR-wastedD 

def bizarreness(A, p):
    return p/(2*np.sqrt(np.pi*A))   #Ratio of perimeter to circumference of circle with same area       


def updateGlobals(state):
    global metrics, adjacencyFrame
    temp = dict(zip(state.key, state.value))
    lowdists  = adjacencyFrame.low.replace(temp)
    highdists = adjacencyFrame.high.replace(temp)
    isSame = lowdists==highdists
    adjacencyFrame['isSame'] = isSame
    adjacencyFrame['lowdist'] = lowdists
    adjacencyFrame['highdist'] = highdists
    
    stConts = [contiguousness(state, i) for i in range(ndistricts)]
    stPops  = [    population(state, i) for i in range(ndistricts)]
    stPerim = [     perimeter(state, i) for i in range(ndistricts)]
    stArea  = [      distArea(state, i) for i in range(ndistricts)]
    
    stBiz   = [bizarreness(stArea[i], stPerim[i]) for i in range(ndistricts)]
    
    metrics = {'contiguousness': stConts,
               'population'    : stPops,
               'bizarreness'   : stBiz,
               'perimeter'     : stPerim,
               'area'          : stArea}


def goodness(metrics):
    #stConts = [contiguousness(runningState[0], i) for i in range(ndistricts)]
    #stPops  = [    population(runningState[0], i) for i in range(ndistricts)]
    #stBiz   = [   bizarreness(runningState[0], i) for i in range(ndistricts)]
    #stPerim = [     perimeter(runningState[0], i) for i in range(ndistricts)]
    #stArea  = [      distArea(runningState[0], i) for i in range(ndistricts)]
    
    tempStConts = metrics['contiguousness']
    tempStPops  = metrics['population']
    tempStBiz   = metrics['bizarreness']
    
    modTotalVar = sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in tempStPops])/(2*(1-float(1)/ndistricts))
    
    return -300*abs(sum(tempStConts) - ndistricts) - 300*modTotalVar - 100*np.nanmean(tempStBiz)

def switchDistrict(current_goodness, possible_goodness): # fix
    return float(1)/(1 + np.exp((current_goodness-possible_goodness)/100.0))




