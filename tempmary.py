import numpy as np
import random
import pandas as pd
import math
from osgeo import ogr
import os

###################################################
def MH(start, steps, neighbor, goodness, moveprob):
    #  object starting state   |         |
    #         integer steps to be taken for M-H algorithm.
    #                function returning a neighbor of current state
    #                          function for determining goodness.
    #                                    function which takes goodnesses and returns probabilities.
    
    current = start.copy()
    best_state = start.copy()
    current_goodness = goodness(current)
    best_goodness = current_goodness
    better_hops = 0
    worse_hops = 0
    stays = 0
    for i in range(steps):
        possible = neighbor(current)
        possible_goodness = goodness(possible)
        if best_goodness < possible_goodness:
            best_state = possible.copy()
            best_goodness = possible_goodness
        if random.random() < moveprob(current_goodness, possible_goodness):
            if current_goodness < possible_goodness :
                better_hops += 1
            else:
                worse_hops += 1
            current = possible.copy()
            current_goodness = possible_goodness
        else:
            stays += 1
    return((best_state, best_goodness, better_hops, worse_hops, stays))

#######################################################################

def MH_swarm(starts, steps, neighbor, goodness, moveprob):
    #        As in MH, but with multiple starts
    walkers = [MH(start, steps, neighbor, goodness, moveprob) for start in starts]
    return sorted(walkers, key = lambda x: x[1], reverse = True)


#######################################################################

def neighbor(state):
    
    newstate = state.copy()
    missingdist = set.difference(set(range(ndistricts)), set(state['value']))
    #If we've blobbed out some districts, we wants to behave differently
    
    if len(missingdist) == 0:
        switchedge = np.random.choice(adjacencyFrame.index[-(adjacencyFrame.isSame == 1)])

        lownode  = adjacencyFrame.low[switchedge]
        highnode = adjacencyFrame.high[switchedge]
        #Randomly choose an adjacency.  Find the low node and high node for that adjacency.

        if random.random() < 0.5:
            newstate.value[newstate.key ==  lownode] = (newstate[newstate.key == highnode].value).item()
            checks = adjacencyFrame.index[((adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)) & \
                                          (-(adjacencyFrame.isSame == 1))]
            adjacencyFrame.isSame[((adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)) & \
                                  adjacencyFrame.isSame] = False
            adjacencyFrame.isSame[checks] = [( newstate.value[newstate.key == adjacencyFrame.low[j]].item() == \
                                               newstate.value[newstate.key == adjacencyFrame.high[j]].item() ) for j in checks]
        else:
            newstate.value[newstate.key == highnode] = (newstate[newstate.key ==  lownode].value).item()
            checks = adjacencyFrame.index[((adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)) & \
                                          (-(adjacencyFrame.isSame == 1))]
            adjacencyFrame.isSame[((adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)) & \
                                  adjacencyFrame.isSame] = False
            adjacencyFrame.isSame[checks] = [( newstate.value[newstate.key == adjacencyFrame.low[j]].item() == \
                                               newstate.value[newstate.key == adjacencyFrame.high[j]].item() ) for j in checks]
        #We want to assign both nodes the same value, and there's a 50% chance for each value being chosen.
    else:
        #If there are some districts missing, 
        changenode = newstate.key.sample(1)
        newstate.value[newstate.key == changenode] = list(missingdist)[0]
        #We want to select one randomly, and make it one of the missing districts
        adjacencyFrame.isSame[(adjacencyFrame.low == changenode) | \
                              (adjacencyFrame.high == changenode)] = False
        # And none of its adjacencies match anymore.
    return newstate

def contiguousness_old(state, district):
    
    regions = 0
    regionlist = list(state.key[state.value == district])
    if len(regionlist) == 0:
        return 1
    
    subframe = adjacencyFrame[[(adjacencyFrame['low'][i] in regionlist) and (adjacencyFrame['high'][i] in regionlist) \
                               for i in range(adjacencyFrame.shape[0])]]
    subedges = subframe[subframe.length != 0][['low','high']]
    
    while len(regionlist) > 0:
        regions += 1
        currentregion = set()
        addons = {regionlist[0]}
        while len(addons) > 0:
            currentregion = currentregion.union(addons)
            subsubedges = subedges[[(subedges['low'][i] in currentregion) or (subedges['high'][i] in currentregion) \
                                    for i in subedges.index]]
            if(not subsubedges.empty):
                addons = set(subsubedges['low']).union(set(subsubedges['high'])) - currentregion
            else:
                addons = set()
        regionlist = [x for x in regionlist if x not in currentregion]
    return regions


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

def perimeter_old(state, district):
    regionlist = list(state.key[state.value == district])
    return sum(adjacencyFrame[[(adjacencyFrame['low'][i] in regionlist) != (adjacencyFrame['high'][i] in regionlist) \
                               for i in range(adjacencyFrame.shape[0])]]['length'])

def perimeter(state, district):
    regionlist = list(state.key[state.value == district])
    return sum(adjacencyFrame.length[adjacencyFrame.low.isin(regionlist) != adjacencyFrame.high.isin(regionlist)])

def interiorPerimeter(state, district):
    regionlist = list(state.key[state.value == district])
    return sum(adjacencyFrame.length[adjacencyFrame.low.isin(regionlist) & adjacencyFrame.high.isin(regionlist)])

def interiorPerimeter_old(state, district):
    regionlist = list(state.key[state.value == district])
    return sum(adjacencyFrame[[(adjacencyFrame['low'][i] in regionlist) and (adjacencyFrame['high'][i] in regionlist) \
                               for i in range(adjacencyFrame.shape[0])]]['length'])

def distArea(state, district):
    regionlist = list(state.key[state.value == district])
    return sum(blockstats.ALAND[blockstats.VTD.isin(regionlist)])

def population(state, district):
    return sum(blockstats.POP100[blockstats.VTD.isin(list(state.key[state.value == district]))])

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

def bizarreness_old(state, district):
    outer = perimeter(state, district)
    inner = interiorPerimeter(state, district)
    if inner + outer == 0:
        return np.nan
    return outer/(inner + outer)

def bizarreness(state, district):
    outer = perimeter(state, district)     #Perimeter of district
    area = distArea(state, district)       #Area of district
    return outer/(2*np.sqrt(np.pi*area))   #Ratio of perimeter to circumference of circle with same area       

def compactness1(state):
    return sum([perimeter(state, district) for district in range(ndistricts)])/2

def goodness(state):
    #Haves
        #contiguousness
        #evenness of population
        #efficiency
        #bizarreness
    #Needs
        #Compactness
    
    stconts = [contiguousness(state, i) for i in range(ndistricts)]
    stpops  = [population(state, i) for i in range(ndistricts)]
    #steffic = [efficiency(state, i) for i in range(ndistricts)]
    stbiz   = [bizarreness(state, i) for i in range(ndistricts)]
    
    modTotalVar = sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in stpops])/(2*(1-float(1)/ndistricts))
    
    #return -3000*abs(sum(stconts) - ndistricts) - 100*modTotalVar - 10*abs(sum(steffic)) -10*np.nansum(stbiz)
    return -300*abs(sum(stconts) - ndistricts) - 100*modTotalVar - 10*np.nansum(stbiz)

def switchDistrict(current_goodness, possible_goodness): # fix
    return float(1)/(1 + np.exp((current_goodness-possible_goodness)/1000.0))

###############################

class Configuration():

    # globals--for all class instances.
    precinctInfo = pd.DataFrame()
    adjacencyFrame = pd.DataFrame()
    nDistricts = 0
    nPrecincts = 0
    stepsTaken = 0


    def __init__(self, demographics_file, connectivities_file, stateName, ndistricts, startingState = 0):
        global nPrecincts, nDistricts, precinctInfo, adjacencyFrame
        blockstats = pd.read_csv(demographics_file)
        blockstats = blockstats.drop('Unnamed: 0', 1)
        blockstats = blockstats.set_index(blockstats.VTD)
        blockstats.rename(columns={'POP100':'population'}, inplace=True)
        precinctInfo = blockstats.copy()

        self.totalpopulation = sum(precinctInfo.population)
        nPrecincts = len(precinctInfo.VTD)
        nDistricts = ndistricts
        del blockstats

        adjacencyFrame = pd.read_csv(connectivities_file)
        adjacencyFrame = adjacencyFrame.drop('Unnamed: 0', 1)
        adjacencyFrame.columns = ['low', 'high', 'length']

        self.runningState = 0

    def contiguousStart(self):
        global nPrecincts, nDistricts, precinctInfo, adjacencyFrame

        state = pd.DataFrame([[precinctInfo.VTD[i], nDistricts] for i in range(nPrecincts)])
        state.columns = ['key', 'value']
        subAdj = adjacencyFrame.loc[adjacencyFrame.length != 0]
        
        missingdist = set(range(nDistricts))
        while len(list(missingdist)) > 0:
            state.value[random.randint(0,nPrecincts-1)] = list(missingdist)[0]
            missingdist = set.difference(set(range(nDistricts)), set(state['value']))
        #Above loop gives each district exactly one VTD.  The rest will be equal to ndistricts
        
        pops = [population(state,x) for x in range(nDistricts)]
        
        while nDistricts in set(state['value']):
            
            targdistr = pops.index(min(pops))
            
            subframe = state.loc[state.value!=nDistricts]
            detDists = set(subframe.key)
            tbdDists = set.difference(set(state.key), detDists)
            relevantAdjacencies = subAdj.loc[(subAdj.low.isin(detDists)) != (subAdj.high.isin(detDists))]
            #adjacencies where either low or high have a value that still has value of ndistricts, but the other doesn't
            curRegion = state.key[state.value == targdistr]
            relevantAdjacencies = subAdj.loc[((subAdj.low.isin(curRegion)) & (subAdj.high.isin(tbdDists))) |
                                             ((subAdj.high.isin(curRegion)) & (subAdj.low.isin(tbdDists)))]
            #Adjacencies where either low or high are in the region, but the other is unassigned
            
            if relevantAdjacencies.shape[0] == 0 :
                pops[targdistr] = float('inf')
            else :
                #choose entry in relevantAdjacencies and switch the value of the other node.
                temp = relevantAdjacencies.loc[relevantAdjacencies.index[random.randint(0,relevantAdjacencies.shape[0]-1)]]
                if temp.high in tbdDists:
                    state.value[state.key == temp.high] = state.value[state.key == temp.low].item()
                    pops[targdistr] = pops[targdistr] + precinctInfo.population[temp.high]
                else:
                    state.value[state.key == temp.low] = state.value[state.key == temp.high].item()
        self.runningState = state

        #return state

    def starting_state(self, startingState = 0):
        global nPrecincts, nDistricts, precinctInfo, adjacencyFrame

        if startingState != 0:
            self.runningState = pd.read_csv(startingState)
        else:
            self.runningState = self.contiguousStart()


        temp = dict(zip(self.runningState.key, self.runningState.value))

        adjacencyFrame['lowdist']  = adjacencyFrame.low.replace(temp)
        adjacencyFrame['highdist'] = adjacencyFrame.high.replace(temp)



    # initialize state space contiguously

thing = Configuration('/home/tsugrad/Documents/gerrymandering/Pennsylvania/vtdstats.csv', '/home/tsugrad/Documents/gerrymandering/Pennsylvania/PRECINCTconnections.csv', 'PA', 18)

thing.starting_state()


