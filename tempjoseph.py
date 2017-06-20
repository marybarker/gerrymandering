

#####
#Supplements to cleanitallup.py
#####

def contScore(metrics):
    if any([x!=1 for x in metrics['contiguousness']]):
        return float('-inf')
    return 1

def popDiffScore(metrics):
    return 1 - float(max(0, (np.max(metrics['population']) - np.min(metrics['population'])) - 25000 ))/ \
               ((totalpopulation - 25000))

def popVarScore(metrics):
    return 1 - sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in metrics['population']])/(2*(1-float(1)/ndistricts))

def bizMeanScore(metrics):
    return 1.0/np.nanmean(metrics['bizarreness'])

def bizMaxScore(metrics):
    return 1.0/max(metrics['bizarreness'])

def aframEdgeScore(metrics):
    mindists = metrics['mincon'].argsort()[-numMajMinDists:][::-1]
    return np.sum(metrics['sumAframDiff'][mindists]) / np.sum(metrics['numedges'])

def hispEdgeScore(metrics):
    mindists = metrics['mincon'].argsort()[-numMajMinDists:][::-1]
    return np.sum(metrics['sumHispDiff'][mindists]) / np.sum(metrics['numedges'])

"""
def goodness(metrics):
    
    #A break in contiguousness should be vetoed; it has a goodness of -infinity.
    if any([x!=1 for x in metrics['contiguousness']]):
        return float('-inf')
        
    #Scores should be scaled such that they are between zero and one,
    #    with one being "good".
    scores  = np.array([popVarScore(metrics),
                        bizMeanScore(metrics),
                        bizMaxScore(metrics),
                        aframEdgeScore(metrics),
                        hispEdgeScore(metrics)])
    
    #Because the scores are normalized, weights are more intuitive,
    #    and roughly correspond to scaling the slope of the metric
    #    in a given dimension.
    weights = np.array([500,
                        100,
                        100,
                        10,
                        10])
    
    #Multiply the scores by their weights, but then renormalize to between zero and one.
    return scores.dot(weights)/np.sum(weights)
"""

goodnessParams  = [contScore, popVarScore, bizMeanScore, bizMaxScore, aframEdgeScore, hispEdgeScore]
goodnessWeights = np.array([1, 500, 100, 100, 10, 10])

def goodness(metrics):
    return sum(x*f(metrics) for f,x in zip(goodnessParams, goodnessWeights))/sum(goodnessWeights)

def contiguousStart(stats = "DEFAULT"):

    #Begin with [ndistricts] different vtds to be the congressional districts.
    #Keep running list of series which are adjacent to the districts.
    #Using adjacencies, let the congressional districts grow by unioning with the remaining districts 

    if type(stats) == str:
        stats = blockstats
    
    state = pd.DataFrame({"key":stats.index.copy(), "value":ndistricts })
    subAdj = adjacencyFrame.ix[adjacencyFrame.low.isin(stats.index) & adjacencyFrame.high.isin(stats.index) & (adjacencyFrame.length != 0), ['low','high']]
    subAdj['lowdist']  = ndistricts
    subAdj['highdist'] = ndistricts
    
    missingdist = range(ndistricts)
    assignments = np.random.choice(stats.index, ndistricts, replace = False)
    
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
                      set(relevantAdjacencies.high)) - set(state.key[state.value == targdistr])
            #changes = set(relevantAdjacencies.low.append(relevantAdjacencies.high))
            #changes = (relevantAdjacencies.low.append(relevantAdjacencies.high)).unique()
            state.ix[state.key.isin(changes), 'value'] = targdistr
            pops[targdistr] += sum(stats.population[changes])
            subAdj.ix[subAdj.low.isin(changes),  'lowdist' ] = targdistr
            subAdj.ix[subAdj.high.isin(changes), 'highdist'] = targdistr
        print("%d districts left to assign."%(sum(state.value==ndistricts)))
    return state.set_index(state.key)


def distArea(state, district):
    regionlist = list(state.key[state.value == district])
    return sum(blockstats.Shape_area[blockstats.ID.isin(regionlist)])

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
            switchNode = lownode
            winnerDist = temphighdist
            loserDist  = templowdist
        else:
            #switch high node stuff to low node's district
            switchNode = highnode
            winnerDist = templowdist
            loserDist  = temphighdist
        
        #Keep track of the parts of adjacencyFrame which could be changing;
        #   Also keep track of previous version of adjacencyFrame in case we want to go back.
        previousVersion = adjacencyFrame[(adjacencyFrame.low == switchNode) | (adjacencyFrame.high == switchNode)]
        proposedChanges = previousVersion.copy()

        newstate.ix[newstate.key == switchNode, 'value'] = winnerDist
        proposedChanges.ix[proposedChanges.low == switchNode, 'lowdist'] = winnerDist
        proposedChanges.ix[proposedChanges.high == switchNode, 'highdist'] = winnerDist
        proposedChanges.isSame = proposedChanges.lowdist == proposedChanges.highdist
        #change values in the state as well as the proposedChanges
        
        #change population
        popChange = blockstats.population[switchNode]
        newmetrics.ix[loserDist, 'population']  -= popChange
        newmetrics.ix[winnerDist, 'population'] += popChange

        #change bizarreness
        newmetrics.ix[loserDist,'perimeter']  += \
            (sum(previousVersion.length[ (previousVersion.isSame==1) & ((previousVersion.low == switchNode) | (previousVersion.high == switchNode))]) -\
             sum(previousVersion.length[~(previousVersion.isSame==1) & ((previousVersion.low == switchNode) | (previousVersion.high == switchNode))]))
        newmetrics.ix[winnerDist, 'perimeter'] += \
            (sum(proposedChanges.length[~(proposedChanges.isSame==1) & ((proposedChanges.low == switchNode) | (proposedChanges.high == switchNode))]) -\
             sum(proposedChanges.length[ (proposedChanges.isSame==1) & ((proposedChanges.low == switchNode) | (proposedChanges.high == switchNode))]))

        areachange = blockstats.ALAND[switchNode] + blockstats.AWATER[switchNode]
        newmetrics.ix[loserDist, 'area'] -= areachange
        newmetrics.ix[winnerDist,'area'] += areachange
        
        newmetrics.ix[loserDist, 'bizarreness']  = bizarreness(newmetrics['area'][loserDist], \
                                                               newmetrics['perimeter'][loserDist])
        newmetrics.ix[winnerDist, 'bizarreness'] = bizarreness(newmetrics['area'][winnerDist], \
                                                               newmetrics['perimeter'][winnerDist])
        
        #update boundary information
        newmetrics.ix[winnerDist,'sumAframDiff'] = newmetrics.ix[winnerDist,'sumAframDiff']\
                                                   + np.sum((-proposedChanges.isSame)*proposedChanges.aframdiff)\
                                                   - np.sum((-previousVersion.isSame)*previousVersion.aframdiff)
        newmetrics.ix[loserDist,'sumAframDiff']  = newmetrics.ix[loserDist,'sumAframDiff']\
                                                   - np.sum((-proposedChanges.isSame)*proposedChanges.aframdiff)\
                                                   + np.sum((-previousVersion.isSame)*previousVersion.aframdiff)
        
        newmetrics.ix[winnerDist,'sumHispDiff']  = newmetrics.ix[winnerDist,'sumHispDiff']\
                                                   + np.sum((-proposedChanges.isSame)*proposedChanges.hispdiff)\
                                                   - np.sum((-previousVersion.isSame)*previousVersion.hispdiff)
        newmetrics.ix[loserDist,'sumHispDiff']   = newmetrics.ix[loserDist,'sumHispDiff']\
                                                   - np.sum((-proposedChanges.isSame)*proposedChanges.hispdiff)\
                                                   + np.sum((-previousVersion.isSame)*previousVersion.hispdiff)
        
        newmetrics.ix[winnerDist,'numedges']     = newmetrics.ix[winnerDist,'numedges']\
                                                   + np.sum(-(proposedChanges.isSame))\
                                                   - np.sum(-(previousVersion.isSame))
        newmetrics.ix[loserDist,'numedges']      = newmetrics.ix[loserDist,'numedges']\
                                                   - np.sum(-(proposedChanges.isSame))\
                                                   + np.sum(-(previousVersion.isSame))
        
        #update contiguousness
        #First check if our switch changes local contiguousness.
        neighborhood = set(proposedChanges.low).union(set(proposedChanges.high))
        
        nhadj = adjacencyFrame.ix[adjacencyFrame.low.isin(neighborhood) & adjacencyFrame.high.isin(neighborhood), ['low','high','length', 'lowdist', 'highdist']]
        oldContNeighborhood = contiguousness(   state.loc[neighborhood], loserDist, nhadj)
        
        nhadj.update(proposedChanges)
        newContNeighborhood = contiguousness(newstate.loc[neighborhood], loserDist, nhadj)
        
        #If local contiguousness changes, check the whole loserDist, since it could be an annulus.
        if (oldContNeighborhood != newContNeighborhood):
            tempframe = adjacencyFrame.copy()
            tempframe.update(proposedChanges)
            tempframe.lowdist  = tempframe.lowdist.astype(int)
            tempframe.highdist = tempframe.highdist.astype(int)
            tempframe.low      = tempframe.low.astype(int)
            tempframe.high     = tempframe.high.astype(int)
            
            newmetrics.ix[loserDist, 'contiguousness']  = contiguousness(newstate, loserDist, tempframe)
        
    else:
        #WARNING:
        #    WE HAVE NEVER GOTTEN TO THIS PIECE OF CODE;  IT HAS NEVER BEEN RUN.
        
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


#####
#Supplements to setup_stuff.py
#####

#None yet.

#####
#Supplements to setup.py
#####

blockstats["aframcon"] = blockstats.e_blak/blockstats.population
blockstats["hispcon"]  = blockstats.e_hsp/blockstats.population
blockstats["mincon"]   = blockstats.e_bh/blockstats.population
blockstats.ix[blockstats.aframcon.isnull(), "aframcon"] = 0
blockstats.ix[blockstats.hispcon.isnull(),  "hispcon" ] = 0
blockstats.ix[blockstats.mincon.isnull(),   "mincon"  ] = 0

adjacencyFrame['low']  = adjacencyFrame.merge(adjacencyFrame.merge(blockstats.ix[:, ['ID', 'VTD']], 
                             left_on='low',  right_on = 'VTD', sort=False)).ix[:, ['ID']]
adjacencyFrame['high'] = adjacencyFrame.merge(adjacencyFrame.merge(blockstats.ix[:, ['ID', 'VTD']], 
                             left_on='high', right_on = 'VTD', sort=False)).ix[:, ['ID']]

conlow = pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "aframcon"]], left_on = 'low', right_on = 'ID').aframcon
conhigh= pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "aframcon"]], left_on = 'high', right_on = 'ID').aframcon
adjacencyFrame["aframdiff"] = conhigh - conlow

conlow = pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "hispcon"]], left_on = 'low', right_on = 'ID').hispcon
conhigh= pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "hispcon"]], left_on = 'high', right_on = 'ID').hispcon
adjacencyFrame["hispdiff"] = conhigh - conlow
