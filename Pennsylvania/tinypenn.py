
##############

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
            
            newmetrics = metrics.copy()
            updateGlobals(current)
            if not dataFrameEquiv(metrics, newmetrics):
                diffperim = metrics.perimeter - newmetrics.perimeter
                print(diffperim)
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
                (sum(previousVersion.length[ (previousVersion.isSame==1) & ((previousVersion.low == lownode) | (previousVersion.highdist == templowdist))]) -\
                 sum(previousVersion.length[-(previousVersion.isSame==1) & ((previousVersion.low == templowdist) | (previousVersion.highdist == templowdist))]))
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
            """
            def perimeter(state, district):
                return sum(adjacencyFrame.length[(adjacencyFrame.lowdist == district) != (adjacencyFrame.highdist == district)])

            """
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


def tiny_color_these_states(subset, geom_to_plot, list_of_states, foldername, number):
    colors = colorDict(ndistricts)
    #colors = {0:'yellow',1:'green'}

    #paths = geom_to_plot['paths']
    #names = geom_to_plot['names']

    thing = zip(geom_to_plot['paths'], geom_to_plot['names'])
    paths = [x[0] for x in thing if x[1] in subset]
    names = [x[1] for x in thing if x[1] in subset]

    for i in range(len(list_of_states)):#state in list_of_states:
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim([-77.1, -76.3])
        ax.set_ylim([39.5, 40.4])
        
        this_state = list_of_states[i]
        redistricting = this_state[0]
        #redistricting = redistricting.drop('Unnamed: 0', 1)
        #redistricting.columns = ['key', 'value']
        for p in range(len(paths)):
            path = paths[p]
            
            facecolor = redistricting.value[np.array(redistricting.key) == names[p]].item()
            patch = mpatches.PathPatch(path,facecolor=colors[facecolor],edgecolor='black')
            ax.add_patch(patch)
        ax.set_aspect(1.0)
        #plt.show()
        plt.savefig(foldername+'output%04d.png'%(number+i), dpi = 600)
        plt.clf()
        del fig



def updateGlobals(state):
    global metrics, adjacencyFrame
    
    lowdists  = pd.merge(adjacencyFrame, state, left_on = 'low' , right_on = 'key').value
    highdists = pd.merge(adjacencyFrame, state, left_on = 'high', right_on = 'key').value
    
    isSame = lowdists==highdists
    adjacencyFrame['isSame'] = isSame
    adjacencyFrame['lowdist'] = lowdists
    adjacencyFrame['highdist'] = highdists
    
    stConts  = [contiguousness(state, i) for i in range(ndistricts)]
    stPops   = [    population(state, i) for i in range(ndistricts)]
    stPerim  = [     perimeter(state, i) for i in range(ndistricts)]
    stArea   = [      distArea(state, i) for i in range(ndistricts)]
    
    stdAfram = [conDiffSum(state, i, 'aframdiff') for i in range(ndistricts)]
    stdHisp  = [conDiffSum(state, i,  'hispdiff') for i in range(ndistricts)]
    
    stMincon = [minorityConc(state, i, 'mincon') for i in range(ndistricts)]
    stBiz    = [bizarreness(stArea[i], stPerim[i]) for i in range(ndistricts)]
    
    metrics  = pd.DataFrame({'contiguousness': stConts,
                             'population'    : stPops,
                             'bizarreness'   : stBiz,
                             'perimeter'     : stPerim,
                             'area'          : stArea,
                             'mincon'        : stMincon,
                             'sumAframDiff'  : stdAfram,
                             'sumHispDiff'   : stdHisp
                             })

random.seed(100)
foldername = "tests/"
#os.mkdir(foldername)
#execfile('setup.py') #Stack overflow doesn't like this, for the record.
#tinyPenn = set(random.sample(blockstats.index, 1)); ndistricts = 3
degree = 2

for i in range(degree):
    tinyPenn = tinyPenn.union(
                                set(adjacencyFrame.high[adjacencyFrame.low.isin(tinyPenn)])
                             ).union(
                                set(adjacencyFrame.low[adjacencyFrame.high.isin(tinyPenn)])
                             )

tinyStats = blockstats.ix[list(tinyPenn), :]
adjacencyFrame = adjacencyFrame.ix[adjacencyFrame.low.isin(tinyPenn) & adjacencyFrame.high.isin(tinyPenn), :]

ndistricts = 4
numstates = 1
numsteps = 50
numsaves = 1000
numplots = 10
startingPoint=0


for j in range(numstates):
    
    tinystate = contiguousStart(tinyStats)
    runningState = (tinystate.copy(), 1)
    updateGlobals(runningState[0])

    for i in range(numsaves):
        
        #oldState = runningState[0].copy()
        #oldAdjacencyFrame = adjacencyFrame.copy()
        #oldMetrics = metrics.copy()
        
        #def MH(start, steps, neighbor, goodness, moveprob):
        runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
        runningState[0].to_csv(foldername+"state%d_save%d.csv"%(j, i + 1), index = False)
        
        #updateGlobalsFromOld(oldState, runningState[0], oldAdjacencyFrame, oldMetrics)
        pd.DataFrame(metrics).to_csv(foldername + 'metrics%d_save%d.csv'%(j, i+1), index = False)
        
        print("Written to state%d_save%d.csv"%(j, i + 1))


