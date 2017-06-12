
time1 = time.time()
for i in range(100):
    dangle = neighbor(state)
time2 = time.time()
dungus  = MH(state, 100, neighbor, goodness, switchDistrict)
time3 = time.time()


t2-t1
t3-t2
t4-t3
t5-t4
t6-t5
t7-t6
t8-t7

t8-t2

##################################################################################

#oldContNeighborhoodHigh  = contiguousness(   state.loc[neighborhood], temphighdist,  previousVersion)
#statecopy = state.copy()
state = statecopy.loc[neighborhood]
district = temphighdist
subframe = previousVersion

#state = statecopy

def contiguousness(state, district, subframe = "DEFAULT"):
    #This function is going to count the numbr of disjoint, connected regeions of the district.
    #The arguments are a state of assignments of precincts to CDs, a district to evaluate, and
    #  a subframe, which is a subset of the adjacencies so we can check contiguousness on a relative
    #  topology.
    t1 = time.time()
    regions = 0
        #start with 0
        
    regionlist = list(state.key[state.value == district])
        #We're going to keep track of the precincts that have been used already.
        
    if len(regionlist) == 0:
        pass
        #If there's nothing in this district...
    #    return 1
            # ... we say that it's in one piece.
    
    #if type(subframe) == str:
        #If the subframe passed is the default, then use anything in the adjacencyframe that's in the district.
    subframe = adjacencyFrame.ix[(adjacencyFrame.lowdist == district) & (adjacencyFrame.highdist == district), :]
    #else:
    
    subedges = subframe[subframe.length != 0][['low','high']]
      
    t1a = np.empty(20)
    t2a = np.empty(20)
    t3a = np.empty(20)
    t4a = np.empty(20)
    
    t2 = time.time()
    while len(regionlist) > 0:
        loops = 0
        regions += 1
        currentregion = set()
        addons = {regionlist[0]}
        while len(addons) > 0:
            t1a[loops] = time.time()
            currentregion = currentregion.union(addons)
            t2a[loops] = time.time()
            subsubedges = subedges.loc[subedges.low.isin(currentregion) | subedges.high.isin(currentregion)]
            t3a[loops] = time.time()
            if(not subsubedges.empty):
                addons = set(subsubedges['low']).union(set(subsubedges['high'])) - currentregion
            else:
                addons = set()
            t4a[loops] = time.time()
            loops += 1
        regionlist = [x for x in regionlist if x not in currentregion]
        
    t1a = t1a[:loops]
    t2a = t2a[:loops]
    t3a = t3a[:loops]
    t4a = t4a[:loops]
    
    t3 = time.time()
    return regions

t4a-t3a
t3a-t2a
t2a-t1a

sum(t4a-t1a)
sum(t4a-t3a) # 0.002
sum(t3a-t2a) # 0.009
sum(t2a-t1a) # 3e-05


def contiguousness(state, district, subframe = "DEFAULT"):
    #This function is going to count the numbr of disjoint, connected regeions of the district.
    #The arguments are a state of assignments of precincts to CDs, a district to evaluate, and
    #  a subframe, which is a subset of the adjacencies so we can check contiguousness on a relative
    #  topology.
    #It is assumed that the precincts are uniquely identified by integers.
    t1 = time.time()
    regions = 0
        #start with 0
        
    regionlist = state.key[state.value == district].values
        #We're going to keep track of the precincts that have been used already.
        
            # ... we say that it's in one piece.
    
    #if type(subframe) == str:
        #If the subframe passed is the default, then use anything in the adjacencyframe that's in the district.
    subframe = adjacencyFrame.ix[(adjacencyFrame.lowdist == district) & (adjacencyFrame.highdist == district), :]
    
    subedges = subframe[subframe.length != 0][['low','high']]
    
    t1a = np.empty(20)
    t2a = np.empty(20)
    t3a = np.empty(20)
    t4a = np.empty(20)
    
    t2 = time.time()
    while len(regionlist) > 0:
        loops = 0
        regions += 1
        currentregion = np.array([], dtype=int)
        addons = np.array([regionlist[0]])
        while len(addons) > 0:
            t1a[loops] = time.time()
            currentregion = np.concatenate((currentregion, addons))
            t2a[loops] = time.time()
            subsubedges = subedges.ix[subedges.low.isin(addons) | subedges.high.isin(addons),:]
            t3a[loops] = time.time()
            if(not subsubedges.empty):
                addons = np.setdiff1d(np.concatenate((subsubedges['low'].values, subsubedges['high'].values)), currentregion)
            else:
                addons = np.array([], dtype=int)
            t4a[loops] = time.time()
            loops += 1
        regionlist = np.setdiff1d(regionlist, currentregion, assume_unique=True)
    
    t3 = time.time()
    
    t1a = t1a[:loops]
    t2a = t2a[:loops]
    t3a = t3a[:loops]
    t4a = t4a[:loops]
    
    return regions

t4a-t3a # 0.002
t3a-t2a # 0.01
t2a-t1a # 3e-05

sum(t4a-t1a)
sum(t4a-t3a)
sum(t3a-t2a)
sum(t2a-t1a)

def neighbor(state):
    
    #stConts = [contiguousness(runningState[0], i) for i in range(ndistricts)]
    #stPops  = [    population(runningState[0], i) for i in range(ndistricts)]
    #stBiz   = [   bizarreness(runningState[0], i) for i in range(ndistricts)]
    #stPerim = [     perimeter(runningState[0], i) for i in range(ndistricts)]
    #stArea  = [      distArea(runningState[0], i) for i in range(ndistricts)]
    global adjacencyFrame, metrics
    
    t1 = time.time()
    
    newstate = state.copy()
    newmetrics = metrics.copy()

    missingdist = set.difference(set(range(ndistricts)), set(newstate['value']))
    #If we've blobbed out some districts, we wants to behave differently
    t2 = time.time()
    if len(missingdist) == 0:
        switchedge = np.random.choice(adjacencyFrame.index[-(adjacencyFrame.isSame == 1)])
        
        lownode      = adjacencyFrame.low[switchedge]
        highnode     = adjacencyFrame.high[switchedge]
        templowdist  = adjacencyFrame.lowdist[switchedge]
        temphighdist = adjacencyFrame.highdist[switchedge]
        #Randomly choose an adjacency.  Find the low node and high node for that adjacency.
        t3 = time.time()
    
        
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
        
        t4 = time.time()
        #change population
        popchange = blockstats.population[lownode]
        newmetrics.ix[templowdist, 'population']  -= popchange
        newmetrics.ix[temphighdist, 'population'] += popchange
        
        t5 = time.time()
        #change bizarreness
        newmetrics.ix[templowdist,'perimeter']  += \
            (sum(previousVersion.length[ (previousVersion.isSame==1) & ((previousVersion.low == lownode) | (previousVersion.high == lownode))]) -\
             sum(previousVersion.length[~(previousVersion.isSame==1) & ((previousVersion.low == lownode) | (previousVersion.high == lownode))]))
        newmetrics.ix[temphighdist, 'perimeter'] += \
            (sum(proposedChanges.length[~(proposedChanges.isSame==1) & ((proposedChanges.low == lownode) | (proposedChanges.high == lownode))]) -\
             sum(proposedChanges.length[ (proposedChanges.isSame==1) & ((proposedChanges.low == lownode) | (proposedChanges.high == lownode))]))

        areachange = blockstats.ALAND[lownode] + blockstats.AWATER[lownode]
        newmetrics.ix[templowdist, 'area'] -= areachange
        newmetrics.ix[temphighdist,'area'] += areachange
        
        newmetrics.ix[templowdist, 'bizarreness']  = bizarreness(newmetrics['area'][templowdist], \
                                                              newmetrics['perimeter'][templowdist])
        newmetrics.ix[temphighdist, 'bizarreness'] = bizarreness(newmetrics['area'][temphighdist], \
                                                              newmetrics['perimeter'][temphighdist])
        
        t6 = time.time()
        #update boundary information
        newmetrics.ix[temphighdist,'sumAframDiff'] = newmetrics.ix[temphighdist,'sumAframDiff']\
                                                   + np.sum((-proposedChanges.isSame)*proposedChanges.aframdiff)\
                                                   - np.sum((-previousVersion.isSame)*previousVersion.aframdiff)
        newmetrics.ix[templowdist,'sumAframDiff'] = newmetrics.ix[templowdist,'sumAframDiff']\
                                                   - np.sum((-proposedChanges.isSame)*proposedChanges.aframdiff)\
                                                   + np.sum((-previousVersion.isSame)*previousVersion.aframdiff)
        
        newmetrics.ix[temphighdist,'sumHispDiff'] = newmetrics.ix[temphighdist,'sumHispDiff']\
                                                   + np.sum((-proposedChanges.isSame)*proposedChanges.hispdiff)\
                                                   - np.sum((-previousVersion.isSame)*previousVersion.hispdiff)
        newmetrics.ix[templowdist,'sumHispDiff'] = newmetrics.ix[templowdist,'sumHispDiff']\
                                                   - np.sum((-proposedChanges.isSame)*proposedChanges.hispdiff)\
                                                   + np.sum((-previousVersion.isSame)*previousVersion.hispdiff)
        
        newmetrics.ix[temphighdist,'numedges'] = newmetrics.ix[temphighdist,'numedges']\
                                               + np.sum(-(proposedChanges.isSame))\
                                               - np.sum(-(previousVersion.isSame))
        newmetrics.ix[templowdist,'numedges'] = newmetrics.ix[templowdist,'numedges']\
                                               - np.sum(-(proposedChanges.isSame))\
                                               + np.sum(-(previousVersion.isSame))
    
        t7 = time.time()
        #update contiguousness
        
        #t1 = time.time()
        
        neighborhood = set(proposedChanges.low).union(set(proposedChanges.high))

        oldContNeighborhoodLow  = contiguousness(   state.loc[neighborhood], templowdist,  previousVersion)
        oldContNeighborhoodHigh = contiguousness(   state.loc[neighborhood], temphighdist, previousVersion)
        newContNeighborhoodLow  = contiguousness(newstate.loc[neighborhood], templowdist,  proposedChanges)
        newContNeighborhoodHigh = contiguousness(newstate.loc[neighborhood], temphighdist, proposedChanges)
        
        #t2 = time.time()
        n1dict = {i: (0,1-float(state.value[i])/ndistricts, float(state.value[i])/ndistricts) for i in list(neighborhood)}
        n1dict[lownode] = tuple((x + 1.0)/2 for x in n1dict[lownode])
        n1dict[highnode] = tuple((x + 1.0)/2 for x in n1dict[highnode])
        color_by_rgb(g, n1dict, "tests/neighborhood.png", linewidth = 0.03, DPI = 1000)
        
        n2dict = {i: (0,1-float(newstate.value[i])/ndistricts, float(newstate.value[i])/ndistricts) for i in list(neighborhood)}
        n2dict[lownode] = tuple((x + 1.0)/2 for x in n2dict[lownode])
        n2dict[highnode] = tuple((x + 1.0)/2 for x in n2dict[highnode])
        color_by_rgb(g, n2dict, "tests/neighborhood2.png", linewidth = 0.03, DPI = 1000)
        
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
        t8 = time.time()
        
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
