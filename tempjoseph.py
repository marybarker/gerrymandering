

#####
#Supplements to cleanitallup.py
#####

def popDiffScore(metrics):
    return float(max(0, (np.max(metrics['population']) - np.min(metrics['population'])) - 25000 )**2)/ \
           ((totalpopulation - 25000)**2)

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

    if any([x!=1 for x in metrics['contiguousness']]):
        return float('-inf')
    
    return (1     * popDiffScore(metrics)   + \
            1     * popVarScore(metrics)    + \
            1     * bizScore(metrics)       + \
            1     * aframEdgeScore(metrics) + \
            1     * hispEdgeScore(metrics)      )
    #functions should be written such that the numbers being scaled are between zero and one.

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




#####
#Supplements to setup.py
#####

blockstats["aframcon"] = blockstats.e_blak/blockstats.population
blockstats["hispcon"]  = blockstats.e_hsp/blockstats.population
blockstats["mincon"]   = blockstats.e_bh/blockstats.population
blockstats.ix[blockstats.aframcon.isnull(), "aframcon"] = 0
blockstats.ix[blockstats.hispcon.isnull(),  "hispcon" ] = 0
blockstats.ix[blockstats.mincon.isnull(),   "mincon"  ] = 0


conlow = pd.merge(adjacencyFrame, blockstats.ix[:, ["VTD", "aframcon"]], left_on = 'low', right_on = 'VTD').aframcon
conhigh= pd.merge(adjacencyFrame, blockstats.ix[:, ["VTD", "aframcon"]], left_on = 'high', right_on = 'VTD').aframcon
adjacencyFrame["aframdiff"] = conhigh - conlow

conlow = pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "hispcon"]], left_on = 'low', right_on = 'ID').hispcon
conhigh= pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "hispcon"]], left_on = 'high', right_on = 'ID').hispcon
adjacencyFrame["hispdiff"] = conhigh - conlow

adjacencyFrame[(-adjacencyFrame.high.isin(blockstats.VTD))|(-adjacencyFrame.low.isin(blockstats.VTD))]
