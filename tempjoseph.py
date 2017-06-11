

#####
#Supplements to cleanitallup.py
#####

def popDiffScore(metrics):
    return 1 - float(max(0, (np.max(metrics['population']) - np.min(metrics['population'])) - 25000 ))/ \
               ((totalpopulation - 25000))

def popVarScore(metrics):
    return 1 - sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in metrics['population']])/(2*(1-float(1)/ndistricts))

def bizScore(metrics):
    return 1.0/np.nanmean(metrics['bizarreness'])

def aframEdgeScore(metrics):
    mindists = metrics['mincon'].argsort()[-numMajMinDists:][::-1]
    return np.sum(metrics['sumAframDiff'][mindists]) / np.sum(metrics['numedges'])

def hispEdgeScore(metrics):
    mindists = metrics['mincon'].argsort()[-numMajMinDists:][::-1]
    return np.sum(metrics['sumHispDiff'][mindists]) / np.sum(metrics['numedges'])

def goodness(metrics):
    
    #A break in contiguousness should be vetoed; it has a goodness of -infinity.
    if any([x!=1 for x in metrics['contiguousness']]):
        return float('-inf')
        
    #Scores should be scaled such that they are between zero and one,
    #    with one being "good".
    scores  = np.array([popVarScore(metrics),
                        bizScore(metrics),
                        aframEdgeScore(metrics),
                        hispEdgeScore(metrics)])
    
    #Because the scores are normalized, weights are more intuitive,
    #    and roughly correspond to scaling the slope of the metric
    #    in a given dimension.
    weights = np.array([1000,
                        10,
                        1,
                        1])
    
    #Multiply the scores by their weights, but then renormalize to between zero and one.
    return scores.dot(weights)/np.sum(weights)

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
