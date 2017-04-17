
time1 = time.clock()
dangle = neighbor(state)
time2 = time.clock()


time6 - time5
time5 - time4
time4 - time3
time3 - time2
time2 - time1


def neighbor(state):
    
    #stConts = [contiguousness(runningState[0], i) for i in range(ndistricts)]
    #stPops  = [    population(runningState[0], i) for i in range(ndistricts)]
    #stBiz   = [   bizarreness(runningState[0], i) for i in range(ndistricts)]
    #stPerim = [     perimeter(runningState[0], i) for i in range(ndistricts)]
    #stArea  = [      distArea(runningState[0], i) for i in range(ndistricts)]
    #global adjacencyFrame, metrics
    
    time1 = time.clock()
    newstate = state.copy()
    newmetrics = metrics.copy()

    missingdist = set.difference(set(range(ndistricts)), set(newstate['value']))
    #If we've blobbed out some districts, we wants to behave differently
    
    time2  = time.clock()
    if len(missingdist) == 0:
        switchedge = np.random.choice(adjacencyFrame.index[-(adjacencyFrame.isSame == 1)])

        lownode      = adjacencyFrame.low[switchedge]
        highnode     = adjacencyFrame.high[switchedge]
        templowdist  = adjacencyFrame.lowdist[switchedge]
        temphighdist = adjacencyFrame.highdist[switchedge]
        #Randomly choose an adjacency.  Find the low node and high node for that adjacency.
        time3  = time.clock()
        if random.random() < 0.5:
            
            #switch low node stuff to high node's district
            
            switchTo = temphighdist
            
            previousVersion = adjacencyFrame[(adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)]
            proposedChanges = previousVersion.copy()
            #want proposedChanges to be a slice of adjacencyFrame where the values could be changing.
            
            newstate.value[newstate.key == lownode] = switchTo
            proposedChanges.lowdist[proposedChanges.low == lownode] = switchTo
            proposedChanges.highdist[proposedChanges.high == lownode] = switchTo
            proposedChanges.isSame = proposedChanges.lowdist == proposedChanges.highdist
            #change values in the state as well as the proposedChanges
            
            #change population
            popchange = blockstats.population[lownode]
            newmetrics['population'][templowdist]  -= popchange
            newmetrics['population'][temphighdist] += popchange
            
            #change bizarreness
            newmetrics['perimeter'][templowdist]  += \
                (sum(previousVersion.length[ (previousVersion.isSame==1) & ((previousVersion.lowdist == templowdist) | (previousVersion.highdist == templowdist))]) -\
                 sum(previousVersion.length[-(previousVersion.isSame==1) & ((previousVersion.lowdist == templowdist) | (previousVersion.highdist == templowdist))]))
            newmetrics['perimeter'][temphighdist] += \
                (sum(proposedChanges.length[-(proposedChanges.isSame==1) & ((proposedChanges.lowdist == templowdist) | (proposedChanges.highdist == templowdist))]) -\
                 sum(proposedChanges.length[ (proposedChanges.isSame==1) & ((proposedChanges.lowdist == templowdist) | (proposedChanges.highdist == templowdist))]))
            
            areachange = blockstats.ALAND[lownode] + blockstats.AWATER[lownode]
            newmetrics['area'][templowdist]  -= areachange
            newmetrics['area'][temphighdist] += areachange
            
            newmetrics['bizarreness'][templowdist]  = bizarreness(newmetrics['area'][templowdist], \
                                                                  newmetrics['perimeter'][templowdist])
            newmetrics['bizarreness'][temphighdist] = bizarreness(newmetrics['area'][temphighdist], \
                                                                  newmetrics['perimeter'][temphighdist])
            
        else:
            
            #switch high node stuff to low node's district
            
            switchTo = templowdist
            #switch to low node
            
            previousVersion = adjacencyFrame[(adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)]
            proposedChanges = previousVersion.copy()
            #want proposedChanges to be a slice of adjacencyFrame where the values could be changing.
            
            newstate.value[newstate.key == highnode] = switchTo
            proposedChanges.lowdist[proposedChanges.low == highnode] = switchTo
            proposedChanges.highdist[proposedChanges.high == highnode] = switchTo
            proposedChanges.isSame = proposedChanges.lowdist == proposedChanges.highdist
            #change values in the state as well as the proposedChanges
            
            #change population
            popchange = blockstats.population[highnode]
            newmetrics['population'][temphighdist] -= popchange
            newmetrics['population'][templowdist]  += popchange
            
            #change bizarreness
            newmetrics['perimeter'][templowdist]  += \
                (sum(proposedChanges.length[-(proposedChanges.isSame==1) & ((proposedChanges.lowdist == temphighdist) | (proposedChanges.highdist == temphighdist))]) -\
                 sum(proposedChanges.length[ (proposedChanges.isSame==1) & ((proposedChanges.lowdist == temphighdist) | (proposedChanges.highdist == temphighdist))]))
            newmetrics['perimeter'][temphighdist] += \
                (sum(previousVersion.length[ (previousVersion.isSame==1) & ((previousVersion.lowdist == temphighdist) | (previousVersion.highdist == temphighdist))]) -\
                 sum(previousVersion.length[-(previousVersion.isSame==1) & ((previousVersion.lowdist == temphighdist) | (previousVersion.highdist == temphighdist))]))
            
            areachange = blockstats.ALAND[highnode] + blockstats.AWATER[highnode]
            newmetrics['area'][temphighdist] -= areachange
            newmetrics['area'][templowdist]  += areachange
            
            newmetrics['bizarreness'][temphighdist] = bizarreness(newmetrics['area'][temphighdist], \
                                                                  newmetrics['perimeter'][temphighdist])
            newmetrics['bizarreness'][templowdist] = bizarreness(newmetrics['area'][templowdist], \
                                                                   newmetrics['perimeter'][templowdist])
        time4  = time.clock()
        #update contiguousness
        neighborhood = set(proposedChanges.low).union(set(proposedChanges.high))
        oldContNeighborhoodLow  = contiguousness(   state.loc[   state.key.isin(neighborhood)], templowdist)
        oldContNeighborhoodHigh = contiguousness(   state.loc[   state.key.isin(neighborhood)], temphighdist)
        newContNeighborhoodLow  = contiguousness(newstate.loc[newstate.key.isin(neighborhood)], templowdist)
        newContNeighborhoodHigh = contiguousness(newstate.loc[newstate.key.isin(neighborhood)], temphighdist)
        time5 = time.clock()
        if (oldContNeighborhoodLow != newContNeighborhoodLow):
            newmetrics['contiguousness'][templowdist]  = contiguousness(newstate, templowdist)
        else:
            pass
        
        if (oldContNeighborhoodHigh != newContNeighborhoodHigh):
            newmetrics['contiguousness'][temphighdist] = contiguousness(newstate, temphighdist)
        else:
            pass
        time6 = time.clock()
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
