import numpy as np
import random
import pandas as pd
import math
from osgeo import ogr
import os

###################################################
def MH_old(start, steps, neighbor, goodness, moveprob):
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


############################
def MH2_old(start, steps, neighbor, goodness, moveprob):
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
        possible_goodness = goodness(possible[0])
        if best_goodness < possible_goodness:
            best_state = possible[0].copy()
            best_goodness = possible_goodness
        if random.random() < moveprob(current_goodness, possible_goodness):
            if current_goodness < possible_goodness :
                better_hops += 1
            else:
                worse_hops += 1
            current = possible[0].copy()
            current_goodness = possible_goodness
            adjacencyFrame.update(possible[1])
        else:
            stays += 1
    return((best_state, best_goodness, better_hops, worse_hops, stays))


############################
def MH(start, steps, neighbor, goodness, moveprob):
    #  object starting state   |         |
    #         integer steps to be taken for M-H algorithm.
    #                function returning a neighbor of current state
    #                          function for determining goodness.
    #                                    function which takes goodnesses and returns probabilities.
    
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
            metrics = possible[2].copy()
        else:
            stays += 1
    return((best_state, best_goodness, better_hops, worse_hops, stays))


#######################################################################

def MH_swarm(starts, steps, neighbor, goodness, moveprob):
    #        As in MH, but with multiple starts
    walkers = [MH(start, steps, neighbor, goodness, moveprob) for start in starts]
    return sorted(walkers, key = lambda x: x[1], reverse = True)


#######################################################################

def neighbor_old(state):
    
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


def neighbor2_old(state):
    
    newstate = state.copy()

    missingdist = set.difference(set(range(ndistricts)), set(newstate['value']))
    #If we've blobbed out some districts, we wants to behave differently
    
    if len(missingdist) == 0:
        switchedge = np.random.choice(adjacencyFrame.index[-(adjacencyFrame.isSame == 1)])

        lownode  = adjacencyFrame.low[switchedge]
        highnode = adjacencyFrame.high[switchedge]
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
            
    else:
        #If there are some districts missing, 
        changenode = newstate.key.sample(1)
        newstate.value[newstate.key == changenode] = list(missingdist)[0]
        #We want to select one randomly, and make it one of the missing districts
        proposedChanges = adjacencyFrame.loc[(adjacencyFrame.low == changenode) | \
                              (adjacencyFrame.high == changenode)]
        proposedChanges.isSame = False
        # And none of its adjacencies match anymore.
    return (newstate, proposedChanges)

def neighbor(state):
    
    #stConts = [contiguousness(runningState[0], i) for i in range(nDistricts)]
    #stPops  = [    population(runningState[0], i) for i in range(nDistricts)]
    #stBiz   = [   bizarreness(runningState[0], i) for i in range(nDistricts)]
    #stPerim = [     perimeter(runningState[0], i) for i in range(nDistricts)]
    #stArea  = [      distArea(runningState[0], i) for i in range(nDistricts)]
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
            newmetrics['contiguousness'][templowdist]  = contiguoussness(newstate, templowdist)
            newmetrics['contiguousness'][temphighdist] = contiguoussness(newstate, temphighdist)
            
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
            
            newmeterics['bizarreness'][templowdist] = bizarreness(newmetrics['area'][templowdist], \
                                                                  newmetrics['perimeter'][templowdist])
            newmeterics['bizarreness'][temphighdist] = bizarreness(newmetrics['area'][temphighdist], \
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
            newmetrics['contiguousness'][temphighdist]  = contiguoussness(newstate, temphighdist)
            newmetrics['contiguousness'][templowdist] = contiguoussness(newstate, templowdist)
            
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
            
            newmeterics['bizarreness'][temphighdist] = bizarreness(newmetrics['area'][temphighdist], \
                                                                  newmetrics['perimeter'][temphighdist])
            newmeterics['bizarreness'][templowdist] = bizarreness(newmetrics['area'][templowdist], \
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
        
        newmeterics['bizarreness'][olddist] = bizarreness(newmetrics['area'][olddist], \
                                                              newmetrics['perimeter'][olddist])
        newmeterics['bizarreness'][newdist] = bizarreness(newmetrics['area'][newdist], \
                                                              newmetrics['perimeter'][newdist])
    
    
    return (newstate, proposedChanges, newmetrics)



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

def interiorPerimeter_old(state, district):
    regionlist = list(state.key[state.value == district])
    return sum(adjacencyFrame[[(adjacencyFrame['low'][i] in regionlist) and (adjacencyFrame['high'][i] in regionlist) \
                               for i in range(adjacencyFrame.shape[0])]]['length'])

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

def bizarreness_old(state, district):
    outer = perimeter(state, district)
    inner = interiorPerimeter(state, district)
    if inner + outer == 0:
        return np.nan
    return outer/(inner + outer)

def bizarreness2_old(state, district):
    outer = perimeter(state, district)     #Perimeter of district
    area = distArea(state, district)       #Area of district
    return outer/(2*np.sqrt(np.pi*area))   #Ratio of perimeter to circumference of circle with same area       

def bizarreness(A, p):
    return p/(2*np.sqrt(np.pi*A))   #Ratio of perimeter to circumference of circle with same area       

def compactness1(state):
    return sum([perimeter(state, district) for district in range(ndistricts)])/2

def goodness_old(state):
    #Haves
        #contiguousness
        #evenness of population
        #efficiency
        #compactness
    #Needs
        #Bizarreness
    
    stconts = [contiguousness(state, i) for i in range(ndistricts)]
    stpops  = [population(state, i) for i in range(ndistricts)]
    #steffic = [efficiency(state, i) for i in range(ndistricts)]
    stbiz   = [bizarreness(state, i) for i in range(ndistricts)]
    
    modTotalVar = sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in stpops])/(2*(1-float(1)/ndistricts))
    
    #return -3000*abs(sum(stconts) - ndistricts) - 100*modTotalVar - 10*abs(sum(steffic)) -10*np.nansum(stbiz)
    return -300*abs(sum(stconts) - ndistricts) - 100*modTotalVar - 10*np.nansum(stbiz)

def updateGlobals(state):
    
    temp = dict(zip(state.key, state.value))
    lowdists  = adjacencyFrame.low.replace(temp)
    highdists = adjacencyFrame.high.replace(temp)
    isSame = lowdists==highdists
    adjacencyFrame['isSame'] = isSame
    adjacencyFrame['lowdist'] = lowdists
    adjacencyFrame['highdist'] = highdists
    
    stConts = [contiguousness(state, i) for i in range(nDistricts)]
    stPops  = [    population(state, i) for i in range(nDistricts)]
    stBiz   = [   bizarreness(state, i) for i in range(nDistricts)]
    stPerim = [     perimeter(state, i) for i in range(nDistricts)]
    stArea  = [      distArea(state, i) for i in range(nDistricts)]
    
    metrics = {'contiguousness': stConts,
               'population'    : stPops,
               'bizarreness'   : stBiz,
               'perimeter'     : stPerim,
               'area'          : stArea}


def goodness(metrics):
    
    #stConts = [contiguousness(runningState[0], i) for i in range(nDistricts)]
    #stPops  = [    population(runningState[0], i) for i in range(nDistricts)]
    #stBiz   = [   bizarreness(runningState[0], i) for i in range(nDistricts)]
    #stPerim = [     perimeter(runningState[0], i) for i in range(nDistricts)]
    #stArea  = [      distArea(runningState[0], i) for i in range(nDistricts)]
    
    tempStConts = metrics['contiguousness']
    tempStPops  = metrics['population']
    tempStBiz   = metrics['bizarreness']
    
    modTotalVar = sum([abs(float(x)/totalpopulation - float(1)/nDistricts) for x in tempStPops])/(2*(1-float(1)/nDistricts))
    
    return -300*abs(sum(tempStConts) - nDistricts) - 100*modTotalVar - 10*np.nansum(tempStBiz)

def switchDistrict(current_goodness, possible_goodness): # fix
    return float(1)/(1 + np.exp((current_goodness-possible_goodness)/100.0))

def contiguousStart_old():
    state = pd.DataFrame([[blockstats.VTD[i], ndistricts] for i in range(0,nvtd)])
    state.columns = ['key', 'value']
    subAdj = adjacencyFrame.loc[adjacencyFrame.length != 0]

    missingdist = set(range(ndistricts))
    while len(list(missingdist)) > 0:
        state.value[random.randint(0,nvtd-1)] = list(missingdist)[0]
        missingdist = set.difference(set(range(ndistricts)), set(state['value']))
    #Above loop gives each district exactly one VTD.  The rest will be equal to ndistricts
    
    while ndistricts in set(state['value']):
         
        subframe = state.loc[state.value!=ndistricts]
        detDists = set(subframe.key)
        tbdDists = set.difference(set(state.key), detDists)
        relevantAdjacencies = subAdj.loc[(subAdj.low.isin(detDists)) != (subAdj.high.isin(detDists))]
        #adjacencies where either low or high have a value that still has value of ndistricts, but the other doesn't
        
        #choose entry in relevantAdjacencies and switch the value of the other node.
        temp = relevantAdjacencies.loc[relevantAdjacencies.index[random.randint(0,relevantAdjacencies.shape[0]-1)]]
        if temp.high in tbdDists:
            state.value[state.key == temp.high] = state.value[state.key == temp.low].item()
        else:
            state.value[state.key == temp.low] = state.value[state.key == temp.high].item()
            
    return state

def contiguousStart2_old():
    state = pd.DataFrame([[blockstats.VTD[i], ndistricts] for i in range(0,nvtd)])
    state.columns = ['key', 'value']
    subAdj = adjacencyFrame.loc[adjacencyFrame.length != 0]
    
    missingdist = set(range(ndistricts))
    while len(list(missingdist)) > 0:
        state.value[random.randint(0,nvtd-1)] = list(missingdist)[0]
        missingdist = set.difference(set(range(ndistricts)), set(state['value']))
    #Above loop gives each district exactly one VTD.  The rest will be equal to ndistricts
    
    pops = [population(state,x) for x in range(ndistricts)]
    
    while ndistricts in set(state['value']):
        
        targdistr = pops.index(min(pops))
        
        subframe = state.loc[state.value!=ndistricts]
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
                pops[targdistr] = pops[targdistr] + blockstats.POP100[temp.high]
            else:
                state.value[state.key == temp.low] = state.value[state.key == temp.high].item()
            
    return state



def contiguousStart3_old():
    state = pd.DataFrame([[blockstats.VTD[i], ndistricts] for i in range(0,nvtd)])
    state.columns = ['key', 'value']
    subAdj = adjacencyFrame.loc[adjacencyFrame.length != 0]
    
    missingdist = range(ndistricts)
    assignments = np.random.choice(state.key, ndistricts)
    state.value[state.key.isin(assignments)] = missingdist
    #Assign a single precinct to each CD.
    
    tbdDists = set(state,key)
    pops = [population(state,x) for x in range(ndistricts)]
    
    while ndistricts in set(state['value']):
        
        targdistr = pops.index(min(pops))
        
        subframe = state.loc[state.value!=ndistricts]
        
    #    relevantAdjacencies = subAdj.loc[(subAdj.low.isin(tbdDists)) != (subAdj.high.isin(tbdDists))]
        #adjacencies where either low or high have a value that still has value of ndistricts, but the other doesn't
        
        relevantAdjacencies = subAdj.loc[((subAdj.low.isin(state.key[state.value == targdistr])) & (subAdj.high.isin(tbdDists))) |
                                         ((subAdj.high.isin(state.key[state.value == targdistr])) & (subAdj.low.isin(tbdDists)))]
        #Adjacencies where either low or high are in the region, but the other is unassigned
        
        if relevantAdjacencies.shape[0] == 0 :
            pops[targdistr] = float('inf')
        else :
            #choose entry in relevantAdjacencies and switch the value of the other node.
            temp = relevantAdjacencies.loc[np.random.choice(relevantAdjacencies.index)]
            if temp.high in tbdDists:
                state.value[state.key == temp.high] = state.value[state.key == temp.low].item()
                pops[targdistr] = pops[targdistr] + blockstats.POP100[temp.high]
            else:
                state.value[state.key == temp.low] = state.value[state.key == temp.high].item()
                pops[targdistr] = pops[targdistr] + blockstats.POP100[temp.low]
    
    return state

def contiguousStart():
    
    #Begin with [ndistricts] different vtds to be the congressional districts.
    #Keep running list of series which are adjacent to the districts.
    #Using adjacencies, let the congressional districts grow by unioning with the remaining districts 
    
    precinctList = list(precinctStats.VTD)
    state = pd.DataFrame([[precinctList[i], nDistricts] for i in range(0,nPrecincts)])
    state.columns = ['key', 'value']
    subAdj = adjacencyFrame.loc[adjacencyFrame.length != 0]
    subAdj['lowdist'] = [nDistricts]*subAdj.shape[0]
    subAdj['highdist'] = [nDistricts]*subAdj.shape[0]
    
    missingdist = range(nDistricts)
    assignments = np.random.choice(precinctList, nDistricts, replace = False)
    state.value[state.key.isin(assignments)] = missingdist
    for i in range(nDistricts):
        subAdj.lowdist[subAdj.low   == assignments[i]] = i
        subAdj.highdist[subAdj.high == assignments[i]] = i
    #Assign a single precinct to each CD.
    
    pops = [population(state,x) for x in range(ndistricts)]
    
    while nDistricts in set(state.value):
        
        targdistr = pops.index(min(pops))
        
        relevantAdjacencies = subAdj.loc[((subAdj.lowdist == targdistr) & (subAdj.highdist == nDistricts)) |
                                         ((subAdj.highdist == targdistr) & (subAdj.lowdist == nDistricts))]
        #Adjacencies where either low or high are in the region, but the other is unassigned
        
        if relevantAdjacencies.shape[0] == 0 :
            pops[targdistr] = float('inf')
        else :
            #choose entry in relevantAdjacencies and switch the value of the other node.
            changes = set(relevantAdjacencies.low).union(\
                      set(relevantAdjacencies.high))
            state.value[state.key.isin(changes)] = targdistr
            pops[targdistr] = pops[targdistr] + sum(blockstats.population[changes])
            subAdj.lowdist[subAdj.low.isin(changes)] = targdistr
            subAdj.highdist[subAdj.high.isin(changes)] = targdistr
                
        print("%d districts left to assign."%(sum(state.value==nDistricts)))
    return state


def bizarreness(state, district, numlines = 1000):
    
    #uses shape file [precinctShapes]
    
    vtds = precinctInfo.VTD
    pvector = np.array([vtds.population], dtype=float)[0,:]/totalPopulation
    #randomly selects [numlines] pairs of precints based on the discrete distribution
    # induced by populations.
    
    districtShape = reduce(lambda x, y : x.Union(y), [f.geometry() for f in precinctShapes \
                           if f.GEOID10 + f.NAME10 in set(state.key[state.value == district])])
    distBounds = districtShape.Boundary()
    
    exits = 0
    
    for i in range(numlines):
        pair = np.random.choice(precinctShapes, 2, p = pvector)
        line = ogr.Geometry(ogr.wkbLineString)
        line.AddPoint(*pair[0].geometry().Centroid().GetPoint())
        line.AddPoint(*pair[1].geometry().Centroid().GetPoint())
        
        if line.Crosses(distBounds):
            exits++
    
    return float(exits)/numlines

def dataFrameEquiv(df1, df2):
    if df1.shape[0] != df2.shape[0]:
        return False
    if not all(df1.columns == df2.columns):
        return False
    
    return all([all(df1[col] == df2[col]) for col in df1.columns])

###############################

"""

#Lookup number of congressional districts state gets
cdtable = pd.read_csv('../../cdbystate1.txt', '\t')
ndistricts = int(cdtable[cdtable['STATE']=='NH'].CD)

#Lookup number of VTDs state has
ds = ogr.Open("./nh_final.shp")
nlay = ds.GetLayerCount()
lyr = ds.GetLayer(0)
nvtd = len(lyr)

#Read adjacency frame
adjacencyFrame = pd.read_csv('../HarvardData/VTDconnections.csv')
adjacencyFrame = adjacencyFrame.drop('Unnamed: 0', 1)
adjacencyFrame.columns = ['low', 'high', 'length']
adjacencyFrame.low  = [x[5:] for x in adjacencyFrame.low]
adjacencyFrame.high = [x[5:] for x in adjacencyFrame.high]

#Read blockstats
blockstats = pd.read_csv("../HarvardData/NHVTDstats.csv")
blockstats = blockstats.drop('Unnamed: 0', 1)
blockstats.set_index(blockstats.VTD)

totalpopulation = sum(blockstats.population)
"""




























