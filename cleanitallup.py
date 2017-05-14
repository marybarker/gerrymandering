import os
import time
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

#os.chdir('/Users/marybarker/Documents/tarleton_misc/gerrymandering/Pennsylvania')
#os.chdir('/home/odin/Documents/gerrymandering/gerrymandering/Pennsylvania')
stateSHORT = 'PA'

blockstats = pd.read_csv("noIslandsVTDStats.csv")
#blockstats.rename(columns = {"POP100":"population"}, inplace = True)
#blockstats = blockstats.drop('Unnamed: 0', 1)
blockstats = blockstats.set_index(blockstats.VTD)
totalpopulation = sum(blockstats.population)

cdtable = pd.read_csv('../cdbystate1.txt', '\t')
ndistricts = int(cdtable[cdtable['STATE']==stateSHORT].CD)
nvtd = len(blockstats.VTD)

adjacencyFrame = pd.read_csv('noIslandsPrecinctConnections.csv')
adjacencyFrame = adjacencyFrame.drop('Unnamed: 0', 1)
#adjacencyFrame.columns = ['low', 'high', 'length']
metrics = pd.DataFrame()

foldername = "bluebrangus1/"
foldername = "ffffffffffffffffffff1/"
foldername = "compare_annealing_vs_not/"
foldername = "million/"

#os.mkdir(foldername)

numstates= 1
numsteps = 400
numsaves = 2500
numplots = 1

for startingpoint in range(numstates):
    starting_state = pd.read_csv('./startingPoints2/start%d.csv'%startingpoint)
    starting_state = pd.read_csv(foldername+'anneal0_2499.csv')
    updateGlobals(starting_state)
    runningState = (starting_state.copy() ,1)
    #for i in range(numsaves):
    for i in range(numsaves, 2*numsaves):
        exploration = 1.0#10**(4*(1 - i*1.0/numsaves))
        runningState = MH(runningState[0], numsteps, neighbor, goodness, anneal)
        runningState[0].to_csv(foldername+'anneal%d_%d.csv'%(startingpoint, i))
        metrics.to_csv(foldername+'metrics_anneal%d_%d.csv'%(startingpoint, i))
        print 'finished with %d with annealing'%i

    updateGlobals(starting_state)
    runningState = (starting_state.copy() ,1)
    for i in range(numsaves):
        exploration = 10**(4*(1 - i*1.0/numsaves))
        runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
        runningState[0].to_csv(foldername+'plain%d_%d.csv'%(startingpoint, i))
        metrics.to_csv(foldername+'metrics_plain%d_%d.csv'%(startingpoint, i))
        print 'finished with %d in no-annealing'%i
    

maxBizArray1 = np.zeros((numstates,numsaves))
meanBizArray1 = np.zeros((numstates,numsaves))
totalVarArray1 = np.zeros((numstates,numsaves))
maxContArray1 = np.zeros((numstates,numsaves))
maxPopArray1 = np.zeros((numstates,numsaves))
popDiffArray1 = np.zeros((numstates,numsaves))

overallGoodnessArray1 = np.zeros((numstates,numsaves))

maxBizArray2 = np.zeros((numstates,numsaves))
meanBizArray2 = np.zeros((numstates,numsaves))
totalVarArray2 = np.zeros((numstates,numsaves))
maxContArray2 = np.zeros((numstates,numsaves))
maxPopArray2 = np.zeros((numstates,numsaves))
popDiffArray2 = np.zeros((numstates,numsaves))

overallGoodnessArray2 = np.zeros((numstates,numsaves))

for startingpoint in range(numstates):
    for j in range(numsaves):
        metrics_anneal = pd.read_csv(foldername+'metrics_anneal%d_%d.csv'%(startingpoint, j))
        #metrics_plain  = pd.read_csv(foldername+'metrics_plain%d_%d.csv'%(startingpoint, j))

        meanBizArray1[startingpoint,j]         = np.mean(metrics_anneal['bizarreness'])
        maxBizArray1[startingpoint,j]          = np.max(metrics_anneal['bizarreness'])
        maxContArray1[startingpoint,j]         = np.max(metrics_anneal['contiguousness'])
        maxPopArray1[startingpoint,j]          = np.max(metrics_anneal['population'])
        popDiffArray1[startingpoint,j]         = np.max(metrics_anneal['population']) - np.min(metrics_anneal['population'])
        totalVarArray1[startingpoint,j]        = np.sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in metrics_anneal['population']])/(2*(1-float(1)/ndistricts))
        overallGoodnessArray1[startingpoint,j] = -300*abs(sum(metrics_anneal['contiguousness']) - ndistricts) - 3000*totalVarArray1[startingpoint,j] - 1000*meanBizArray1[startingpoint,j] - \
            float(max(0, popDiffArray1[startingpoint,j] - 25000 )**2)/1000000

        meanBizArray2[startingpoint,j]         = np.mean(metrics_plain['bizarreness'])
        maxBizArray2[startingpoint,j]          = np.max(metrics_plain['bizarreness'])
        maxContArray2[startingpoint,j]         = np.max(metrics_plain['contiguousness'])
        maxPopArray2[startingpoint,j]          = np.max(metrics_plain['population'])
        popDiffArray2[startingpoint,j]         = np.max(metrics_plain['population']) - np.min(metrics_plain['population'])
        totalVarArray2[startingpoint,j]        = np.sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in metrics_plain['population']])/(2*(1-float(1)/ndistricts))
        overallGoodnessArray2[startingpoint,j] = -300*abs(sum(metrics_plain['contiguousness']) - ndistricts) - 3000*totalVarArray2[startingpoint,j] - 1000*meanBizArray2[startingpoint,j] - \
            float(max(0, popDiffArray2[startingpoint,j] - 25000 )**2)/1000000

for i in range(numstates):
    plt.title('max biz array')
    plt.plot(maxBizArray1[i,:], label='annealing')
    #plt.plot(maxBizArray2[i,:], label='no annealing')
    plt.legend()
    plt.show()
    plt.clf()
    
for i in range(numstates):
    plt.title('mean biz array')
    plt.plot(meanBizArray1[i,:], label='annealing')
    #plt.plot(meanBizArray2[i,:], label='no annealing')
    plt.legend()
    plt.show()
    plt.clf()

for i in range(numstates):
    plt.title('max cont array')
    plt.plot(maxContArray1[i,:], label='annealing')
    #plt.plot(maxContArray2[i,:], label='no annealing')
    plt.legend()
plt.show()

for i in range(numstates):
    plt.title('max pop array')
    plt.plot(maxPopArray1[i,:], label='annealing')
    #plt.plot(maxPopArray2[i,:], label='no annealing')
    plt.legend()
    plt.show()
    plt.clf()

for i in range(numstates):
    plt.title('mean pop array')
    plt.plot(popDiffArray1[i,:], label='annealing')
    #plt.plot(popDiffArray2[i,:], label='no annealing')
    plt.legend()
    plt.show()
    plt.clf()

for i in range(numstates):
    plt.title('overall goodness array')
    plt.plot(-overallGoodnessArray1[i,:], label='annealing')
    #plt.plot(-overallGoodnessArray2[i,:], label='no annealing')
    plt.gca().invert_yaxis()
    plt.gca().set_yscale('log')
    plt.tick_params(axis='y', which='minor')
    plt.legend()
    plt.show()
    plt.clf()






startingpoint = 14
i = 999
runningState = (pd.read_csv(foldername+"state%d_save%d.csv"%(startingpoint, i + 1)),0)
color_these_states(g, [runningState], 'two_million_', 0)

foldername = 'bluebrangus2/'
#os.mkdir(foldername)

for startingpoint in range(75, 76):
    
    starting_state = pd.read_csv('./startingPoints2/start%d.csv'%startingpoint)
    updateGlobals(starting_state)
    runningState = (starting_state.copy(), 1)
    for i in range(numsaves):
        exploration = 10**(4 - i/250.0)
        
        #oldState = runningState[0].copy()
        #oldAdjacencyFrame = adjacencyFrame.copy()
        #oldMetrics = metrics.copy()
        
        #def MH(start, steps, neighbor, goodness, moveprob):
        runningState = MH(runningState[0], numsteps, neighbor, goodness, anneal)
        #runningState[0].to_csv(foldername+"state%d_save%d.csv"%(startingpoint, i + 1), index = False)
        
        #updateGlobalsFromOld(oldState, runningState[0], oldAdjacencyFrame, oldMetrics)
        metrics.to_csv(foldername + 'metrics%d_save%d.csv'%(startingpoint, i+1), index = False)
        
        print("Written to state%d_save%d.csv"%(startingpoint, i + 1))
    runningState[0].to_csv(foldername+"state%d_save%d.csv"%(startingpoint, i + 1), index = False)



maxBizArray = np.zeros((numstates,numsaves))
meanBizArray = np.zeros((numstates,numsaves))
totalVarArray = np.zeros((numstates,numsaves))
maxContArray = np.zeros((numstates,numsaves))
maxPopArray = np.zeros((numstates,numsaves))
popDiffArray = np.zeros((numstates,numsaves))

overallGoodnessArray = np.zeros((numstates,numsaves))

for startingpoint in range(75):
    for j in range(numsaves):
        #tempstate = pd.read_csv(foldername + "state%d_save%d.csv"%(i, j+1))
        #updateGlobals(tempstate)
        #pd.DataFrame(metrics).to_csv(foldername + 'metrics%d_save%d.csv'%(1, i+50), index = False)
        thismetrics = pd.read_csv(foldername+'metrics%d_save%d.csv'%(startingpoint, j+1))
        
        meanBizArray[startingpoint,j]         = np.mean(thismetrics['bizarreness'])
        maxBizArray[startingpoint,j]          = np.max(thismetrics['bizarreness'])
        maxContArray[startingpoint,j]         = np.max(thismetrics['contiguousness'])
        maxPopArray[startingpoint,j]          = np.max(thismetrics['population'])
        popDiffArray[startingpoint,j]         = np.max(thismetrics['population']) - np.min(thismetrics['population'])
        totalVarArray[startingpoint,j]        = np.sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in thismetrics['population']])/(2*(1-float(1)/ndistricts))
        overallGoodnessArray[startingpoint,j] = -300*abs(sum(thismetrics['contiguousness']) - ndistricts) - 3000*totalVarArray[startingpoint,j] - 1000*meanBizArray[startingpoint,j] - \
            float(max(0, popDiffArray[startingpoint,j] - 25000 )**2)/1000000
        #metrics = {'contiguousness': metrics['contiguousness'],
        #           'population'    : stPops,
        #           'bizarreness'   : stBiz,
        #           'perimeter'     : stPerim,
        #           'area'          : stArea}
    
    print("Stored metrics for state %d"%(startingpoint))

num = len(meanBizArray)

startingpoint = 13

plt.plot(meanBizArray[startingpoint,:])
plt.title('mean Biz')
plt.show()
plt.clf()

plt.plot(maxBizArray[startingpoint,:])
plt.title('max Biz')
plt.show()
plt.clf()

plt.plot(maxContArray[startingpoint,:])
plt.title('max contig')
plt.show()
plt.clf()

plt.plot(maxPopArray[startingpoint,:])
plt.title('max pop')
plt.show()
plt.clf()

plt.plot(popDiffArray[startingpoint,:])
plt.title('pop diff')
plt.show()
plt.clf()

plt.plot(totalVarArray[startingpoint,:])
plt.title('mean Pop')
plt.show()
plt.clf()

plt.plot(overallGoodnessArray[startingpoint,:])
plt.title('goodness')
plt.show()

for i in range(numstates):
    plt.plot(meanBizArray[i,:])
plt.title('mean Biz')
plt.show()
plt.clf()
for i in range(numstates):
    plt.plot(maxBizArray[i,:])
plt.title('max Biz')
plt.show()
plt.clf()
for i in range(numstates):
    plt.plot(maxContArray[i,:])
plt.title('max contig')
plt.show()
plt.clf()
for i in range(numstates):
    plt.plot(maxPopArray[i,:])
plt.title('max pop')
plt.show()
plt.clf()
for i in range(numstates):
    plt.plot(popDiffArray[i,:])
plt.title('pop diff')
plt.show()
plt.clf()
for i in range(numstates):
    plt.plot(totalVarArray[i,:])
plt.title('mean Pop')
plt.show()
plt.clf()
for i in range(numstates):
    plt.plot(overallGoodnessArray[i,:])
plt.title('goodness')
plt.show()

color_these_states(g, [(tempstate, 0)], foldername+'theverylast_', 0)
tempstate = pd.read_csv(foldername + "state%d_save%d.csv"%(1, 1))
color_these_states(g, [(tempstate, 0)], foldername+'theveryfirst_', 0)


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
    best_adjacency = adjacencyFrame.copy()
    best_metrics = metrics.copy()
    
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
            best_adjacency = adjacencyFrame.copy()
            best_adjacency.update(possible[1])
            #best_adjacency.lowdist  = best_adjacency.lowdist.astype(int)
            #best_adjacency.highdist = best_adjacency.highdist.astype(int)
            
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
    
    adjacencyFrame.update(best_adjacency)
    adjacencyFrame.lowdist  = adjacencyFrame.lowdist.astype(int)
    adjacencyFrame.highdist = adjacencyFrame.highdist.astype(int)
    #Update adjacencyframe to the best that we ever had.
    
    metrics = best_metrics.copy()
    
    # " for metrics
    
    return((best_state, best_goodness, better_hops, worse_hops, stays))


def neighbor(state):
    
    #stConts = [contiguousness(runningState[0], i) for i in range(ndistricts)]
    #stPops  = [    population(runningState[0], i) for i in range(ndistricts)]
    #stBiz   = [   bizarreness(runningState[0], i) for i in range(ndistricts)]
    #stPerim = [     perimeter(runningState[0], i) for i in range(ndistricts)]
    #stArea  = [      distArea(runningState[0], i) for i in range(ndistricts)]
    global adjacencyFrame, metrics
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
            
            switchTo = temphighdist
            
            previousVersion = adjacencyFrame[(adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)]
            proposedChanges = previousVersion.copy()
            #want proposedChanges to be a slice of adjacencyFrame where the values could be changing.
    
            newstate.ix[newstate.key == lownode, 'value'] = switchTo
            proposedChanges.ix[proposedChanges.low == lownode, 'lowdist'] = switchTo
            proposedChanges.ix[proposedChanges.high == lownode, 'highdist'] = switchTo
            proposedChanges.isSame = proposedChanges.lowdist == proposedChanges.highdist
            #change values in the state as well as the proposedChanges
            
            #change population
            popchange = blockstats.population[lownode]
            newmetrics.ix[templowdist, 'population']  -= popchange
            newmetrics.ix[temphighdist, 'population'] += popchange
    
            #change bizarreness
            newmetrics.ix[templowdist,'perimeter']  += \
                (sum(previousVersion.length[ (previousVersion.isSame==1) & ((previousVersion.lowdist == templowdist) | (previousVersion.highdist == templowdist))]) -\
                 sum(previousVersion.length[-(previousVersion.isSame==1) & ((previousVersion.lowdist == templowdist) | (previousVersion.highdist == templowdist))]))
            newmetrics.ix[temphighdist, 'perimeter'] += \
                (sum(proposedChanges.length[-(proposedChanges.isSame==1) & ((proposedChanges.lowdist == templowdist) | (proposedChanges.highdist == templowdist))]) -\
                 sum(proposedChanges.length[ (proposedChanges.isSame==1) & ((proposedChanges.lowdist == templowdist) | (proposedChanges.highdist == templowdist))]))
    
            areachange = blockstats.ALAND[lownode] + blockstats.AWATER[lownode]
            newmetrics.ix[templowdist, 'area'] -= areachange
            newmetrics.ix[temphighdist,'area'] += areachange
            
            newmetrics.ix[templowdist, 'bizarreness']  = bizarreness(newmetrics['area'][templowdist], \
                                                                  newmetrics['perimeter'][templowdist])
            newmetrics.ix[temphighdist, 'bizarreness'] = bizarreness(newmetrics['area'][temphighdist], \
                                                                  newmetrics['perimeter'][temphighdist])

        else:
            
            #switch high node stuff to low node's district
            
            switchTo = templowdist
            #switch to low node
            
            previousVersion = adjacencyFrame[(adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)]
            proposedChanges = previousVersion.copy()
            #want proposedChanges to be a slice of adjacencyFrame where the values could be changing.
            
            newstate.ix[newstate.key == highnode, 'value'] = switchTo
            proposedChanges.ix[proposedChanges.low == highnode, 'lowdist'] = switchTo
            proposedChanges.ix[proposedChanges.high == highnode, 'highdist'] = switchTo
            proposedChanges.isSame = proposedChanges.lowdist == proposedChanges.highdist
            #change values in the state as well as the proposedChanges
            
            #change population
            popchange = blockstats.population[highnode]
            newmetrics.ix[temphighdist, 'population'] -= popchange
            newmetrics.ix[templowdist, 'population']  += popchange
            
            #change bizarreness
            newmetrics.ix[templowdist, 'perimeter']  += \
                (sum(proposedChanges.length[-(proposedChanges.isSame==1) & ((proposedChanges.lowdist == temphighdist) | (proposedChanges.highdist == temphighdist))]) -\
                 sum(proposedChanges.length[ (proposedChanges.isSame==1) & ((proposedChanges.lowdist == temphighdist) | (proposedChanges.highdist == temphighdist))]))
            newmetrics.ix[temphighdist, 'perimeter'] += \
                (sum(previousVersion.length[ (previousVersion.isSame==1) & ((previousVersion.lowdist == temphighdist) | (previousVersion.highdist == temphighdist))]) -\
                 sum(previousVersion.length[-(previousVersion.isSame==1) & ((previousVersion.lowdist == temphighdist) | (previousVersion.highdist == temphighdist))]))
            
            areachange = blockstats.ALAND[highnode] + blockstats.AWATER[highnode]
            newmetrics.ix[temphighdist, 'area'] -= areachange
            newmetrics.ix[templowdist, 'area']  += areachange
            
            newmetrics.ix[temphighdist, 'bizarreness'] = bizarreness(newmetrics['area'][temphighdist], \
                                                                     newmetrics['perimeter'][temphighdist])
            newmetrics.ix[templowdist, 'bizarreness'] = bizarreness(newmetrics['area'][templowdist], \
                                                                    newmetrics['perimeter'][templowdist])
        
        #update contiguousness
        neighborhood = set(proposedChanges.low).union(set(proposedChanges.high))
        oldContNeighborhoodLow  = contiguousness(   state.loc[   state.key.isin(neighborhood)], templowdist,  proposedChanges)
        oldContNeighborhoodHigh = contiguousness(   state.loc[   state.key.isin(neighborhood)], temphighdist, proposedChanges)
        newContNeighborhoodLow  = contiguousness(newstate.loc[newstate.key.isin(neighborhood)], templowdist,  proposedChanges)
        newContNeighborhoodHigh = contiguousness(newstate.loc[newstate.key.isin(neighborhood)], temphighdist, proposedChanges)
        
        if ((oldContNeighborhoodLow != newContNeighborhoodLow)|(oldContNeighborhoodHigh != newContNeighborhoodHigh)):
            
            tempframe = adjacencyFrame.copy()
            tempframe.update(proposedChanges)
            tempframe.lowdist  = tempframe.lowdist.astype(int)
            tempframe.highdist = tempframe.highdist.astype(int)
            
            if (oldContNeighborhoodLow != newContNeighborhoodLow):
                newmetrics.ix[templowdist, 'contiguousness']  = contiguousness(newstate, templowdist, tempframe)
            else:
                pass
            
            if (oldContNeighborhoodHigh != newContNeighborhoodHigh):
                newmetrics.ix[temphighdist, 'contiguousness'] = contiguousness(newstate, temphighdist, tempframe)
            else:
                pass
        
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


def contiguousness(state, district, subframe = "DEFAULT"):
    
    regions = 0
    regionlist = list(state.key[state.value == district])
    if len(regionlist) == 0:
        return 1
    
    if type(subframe) == str:
        subframe = adjacencyFrame.loc[(adjacencyFrame.lowdist == district) & (adjacencyFrame.highdist == district)]
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
    return sum(adjacencyFrame.length[(adjacencyFrame.lowdist == district) != (adjacencyFrame.highdist == district)])

def interiorPerimeter(state, district):
    return sum(adjacencyFrame.length[(adjacencyFrame.lowdist == district) & (adjacencyFrame.highdist == district)])

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
    
    metrics = pd.DataFrame({'contiguousness': stConts,
                            'population'    : stPops,
                            'bizarreness'   : stBiz,
                            'perimeter'     : stPerim,
                            'area'          : stArea}
                          )


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
    
    return -30000*abs(sum(tempStConts) - ndistricts) - 3000*modTotalVar - 300*np.nanmean(tempStBiz) - \
            float(max(0, (np.max(tempStPops) - np.min(tempStPops)) - 25000 )**2)/1000000000

def switchDistrict(current_goodness, possible_goodness): # fix
    return float(1)/(1 + np.exp((current_goodness-possible_goodness)/10.0))

def anneal(current_goodness, possible_goodness): # fix
    return 1/(1 + np.exp(float(current_goodness-possible_goodness)/exploration))

def updateGlobalsFromOld(state1, state2, oldAdjacencyFrame, oldMetrics):
    global metrics, adjacencyFrame
    
    adjacencyFrame = oldAdjacencyFrame
    
    substate1 = state1.loc[state1.value != state2.value]
    if substate1.shape[0] == 0:
        return
    
    substate2 = state2.loc[state1.value != state2.value]
    subadjacency1 = oldAdjacencyFrame.loc[oldAdjacencyFrame.low.isin(substate2.key) | oldAdjacencyFrame.high.isin(substate2.key)]
    subadjacency2 = subadjacency1.copy()
    temp = dict(zip(state2.key, state2.value))
    
    lowdists  = subadjacency2.low.replace(temp)
    highdists = subadjacency2.high.replace(temp)
    isSame = lowdists==highdists
    subadjacency2['isSame']   = isSame
    subadjacency2['lowdist']  = lowdists
    subadjacency2['highdist'] = highdists
    
    popdiff   = [population(substate2, i) - population(substate1, i) for i in range(ndistricts)]
    areadiff  = [  distArea(substate2, i) -   distArea(substate1, i) for i in range(ndistricts)]
    perimdiff = [sum(subadjacency2.length[~(subadjacency2.isSame==1) & ((subadjacency2.lowdist == i) | (subadjacency2.highdist == i))])-\
                 sum(subadjacency1.length[~(subadjacency1.isSame==1) & ((subadjacency1.lowdist == i) | (subadjacency1.highdist == i))]) for i in range(ndistricts)]
    
    stPops  = [metrics['population'][i] + popdiff[i]   for i in range(ndistricts)]
    stPerim = [metrics['perimeter'][i]  + perimdiff[i] for i in range(ndistricts)]
    stArea  = [metrics['area'][i]       + areadiff[i]  for i in range(ndistricts)]
    
    stBiz   = [bizarreness(stArea[i], stPerim[i]) for i in range(ndistricts)]
    
    metrics = {'contiguousness': metrics['contiguousness'],
               'population'    : stPops,
               'bizarreness'   : stBiz,
               'perimeter'     : stPerim,
               'area'          : stArea}
    
    adjacencyFrame.update(subadjacency2)
    adjacencyFrame.lowdist  = adjacencyFrame.lowdist.astype(int)
    adjacencyFrame.highdist = adjacencyFrame.highdist.astype(int)

def contiguousStart():
    
    #Begin with [ndistricts] different vtds to be the congressional districts.
    #Keep running list of series which are adjacent to the districts.
    #Using adjacencies, let the congressional districts grow by unioning with the remaining districts 

    state = pd.DataFrame({"key":blockstats.VTD.copy(), "value":ndistricts })
    subAdj = adjacencyFrame.ix[adjacencyFrame.length != 0, ['low','high']]
    subAdj['lowdist']  = ndistricts
    subAdj['highdist'] = ndistricts
    
    missingdist = range(ndistricts)
    assignments = np.random.choice(blockstats.VTD, ndistricts, replace = False)
    
    state.ix[state.key.isin(assignments), 'value'] = missingdist
    for i in range(ndistricts):
        subAdj.ix[subAdj.low  == assignments[i], 'lowdist' ] = i
        subAdj.ix[subAdj.high == assignments[i], 'highdist'] = i
    #Assign a single precinct to each CD.
    
    pops = [population(state,x) for x in range(ndistricts)]
    
    while ndistricts in set(state.value):
        
        targdistr = pops.index(min(pops))
        
        relevantAdjacencies = subAdj.ix[((subAdj.lowdist  == targdistr) & (subAdj.highdist == ndistricts)) |
                                        ((subAdj.highdist == targdistr) & (subAdj.lowdist  == ndistricts))]
        #Adjacencies where either low or high are in the region, but the other is unassigned
        
        if relevantAdjacencies.shape[0] == 0 :
            pops[targdistr] = float('inf')
        else :
            #choose entry in relevantAdjacencies and switch the value of the other node.
            changes = set(relevantAdjacencies.low).union(\
                      set(relevantAdjacencies.high))
            #changes = set(relevantAdjacencies.low.append(relevantAdjacencies.high))
            #changes = (relevantAdjacencies.low.append(relevantAdjacencies.high)).unique()
            state.ix[state.key.isin(changes), 'value'] = targdistr
            pops[targdistr] += sum(blockstats.population[changes])
            subAdj.ix[subAdj.low.isin(changes),  'lowdist' ] = targdistr
            subAdj.ix[subAdj.high.isin(changes), 'highdist'] = targdistr
        #print("%d districts left to assign."%(sum(state.value==ndistricts)))
    return state






