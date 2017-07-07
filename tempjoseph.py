
def demoWastedVotes(state, district, demo, party1, party2):
    
    subframe = blockstats.ix[blockstats.ID.isin(list(state.key[state.value == district])), [demo, party1, party2]]
    
    p1Votes  = sum(subframe[party1])
    p2Votes  = sum(subframe[party2])
    numVotes = p1Votes + p2Votes
    
    if p1Votes > p2Votes:
        #p1 wins, waste sum(min/pop*p2)
        #         waste (p1Votes - 0.5*numVotes)*sum(min)/sum(pop)
        return float(sum(blockstats[demo]*blockstats[party2]))/numVotes + \
               (p1Votes - 0.5*numVotes)*sum(blockstats[demo])/numVotes
    else:
        #then p2 wins, do the opposite
        return float(sum(blockstats[demo]*blockstats[party1]))/numVotes + \
               (p2Votes - 0.5*numVotes)*sum(blockstats[demo])/numVotes

def demoEfficiency(state, demo, popcol, party1, party2):
    wasted = [0,0]
    for district in range(ndistricts):
        
        subframe = blockstats.ix[blockstats.ID.isin(list(state.key[state.value == district])), [demo, party1, party2]]
        
        p1Votes  = sum(subframe[party1])
        p2Votes  = sum(subframe[party2])
        numVotes = p1Votes + p2Votes
        
        if p1Votes > p2Votes:
            #p1 wins, waste sum(min/pop*p2)
            #         waste (p1Votes - 0.5*numVotes)*sum(min)/sum(pop)
            wasted[0] = wasted[0] + float(sum(blockstats[demo]*blockstats[party2]))/numVotes + \
                                    (p1Votes - 0.5*numVotes)*sum(blockstats[demo])/numVotes
            wasted[1] = wasted[1] + float(sum((blockstats[popcol] - blockstats[demo])*blockstats[party2]))/numVotes + \
                                    (p1Votes - 0.5*numVotes)*sum((blockstats[popcol] - blockstats[demo]))/numVotes
        else:
            #then p2 wins, do the opposite
            wasted[0] = wasted[0] + float(sum(blockstats[demo]*blockstats[party1]))/numVotes + \
                                    (p2Votes - 0.5*numVotes)*sum(blockstats[demo])/numVotes
            wasted[1] = wasted[1] + float(sum((blockstats[popcol] - blockstats[demo])*blockstats[party1]))/numVotes + \
                                    (p2Votes - 0.5*numVotes)*sum((blockstats[popcol] - blockstats[demo]))/numVotes
    return [float(wasted[0])/sum(blockstats[demo]), float(wasted[1])/sum(blockstats[popcol] - blockstats[demo])]

def singleDemoEfficiency(state, demo, popcol = "population"):
    wasted = [0,0]
    for district in range(ndistricts):
        subframe = blockstats.ix[blockstats.ID.isin(list(state.key[state.value == district])), [demo, popcol]]
        numVotes = sum(subframe[popcol])
        posVotes = sum(subframe[demo])
        negVotes = numVotes - posVotes
        
        if posVotes > negVotes:
            wasted[0] = wasted[0] + posVotes - 0.5*numVotes
            wasted[1] = wasted[1] + negVotes
        else:
            wasted[1] = wasted[1] + negVotes - 0.5*numVotes
            wasted[0] = wasted[0] + posVotes
    return [wasted[0]/sum(blockstats[demo]), wasted[1]/sum(blockstats[popcol])]


###
#Temporary run
###
republican = 'g2012_USH_rv'
democrat   = 'g2012_USH_dv'
hispcol    = 'hispPop'
aframcol   = 'aframPop'

statereads = 1
savereads  = numsaves

gridrange = [paramList[x] for x in distinctParam]

arrayDict = {"maxBiz"    : (np.zeros((len(gridrange)*statereads,numsaves)), "Maximum Bizarreness"                          ),
             "meanBiz"   : (np.zeros((len(gridrange)*statereads,numsaves)), "Mean Bizarreness"                             ),
             "totalVar"  : (np.zeros((len(gridrange)*statereads,numsaves)), "Total Population Variation"                   ),
             "maxCont"   : (np.zeros((len(gridrange)*statereads,numsaves)), "Maximum Contiguousness"                       ),
             "maxPop"    : (np.zeros((len(gridrange)*statereads,numsaves)), "Maximum Population"                           ),
             "popDiff"   : (np.zeros((len(gridrange)*statereads,numsaves)), "Maximum Population Difference"                ),
             "hispDiff"  : (np.zeros((len(gridrange)*statereads,numsaves)), "Hispanic Boundary Difference Measure"         ),
             "aframDiff" : (np.zeros((len(gridrange)*statereads,numsaves)), "African American Boundary Difference Measure" ),
             "baseEff"   : (np.zeros((len(gridrange)*statereads,numsaves)), "Baseline Vote Efficiency"                     ),
             "aframEff"  : (np.zeros((len(gridrange)*statereads,numsaves)), "African American Vote Efficiency"             ),
             "hispEff"   : (np.zeros((len(gridrange)*statereads,numsaves)), "Hispanic Vote Efficiency"                     ),
             "goodness"  : (np.zeros((len(gridrange)*statereads,numsaves)), "Goodness"                                     )}

for startingpoint in range(statereads):
    for point in range(len(gridrange)):
        weights = gridrange[point]
        subfoldername = "grid%04d.%04d.%04d/"%(weights[1],weights[3],weights[5])
        position = startingpoint * len(gridrange) + point
        
        for j in range(savereads):
            
            thismetrics = pd.read_csv(foldername+subfoldername + 'metrics%04d_save%04d.csv'%(startingpoint, j+1))
            
            mindists = thismetrics['mincon'].argsort()[-numMajMinDists:][::-1]
            
            arrayDict["maxBiz"   ][0][position,j] = np.max(thismetrics['bizarreness'])
            arrayDict["meanBiz"  ][0][position,j] = np.mean(thismetrics['bizarreness'])
            arrayDict["totalVar" ][0][position,j] = np.sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in thismetrics['population']])/(2*(1-float(1)/ndistricts))
            arrayDict["maxCont"  ][0][position,j] = np.max(thismetrics['contiguousness'])
            arrayDict["maxPop"   ][0][position,j] = np.max(thismetrics['population'])
            arrayDict["popDiff"  ][0][position,j] = np.max(thismetrics['population']) - np.min(thismetrics['population'])
            arrayDict["hispDiff" ][0][position,j] = np.sum(thismetrics['sumHispDiff'][mindists])
            arrayDict["aframDiff"][0][position,j] = np.sum(thismetrics['sumAframDiff'][mindists])
            arrayDict["goodness" ][0][position,j] = goodness(thismetrics)
        
        print("Loaded grid%04d.%04d.%04d state %04d"%(weights[1],weights[3],weights[5], startingpoint))

for startingpoint in range(statereads):
    for point in range(len(gridrange)):
        weights = gridrange[point]
        subfoldername = "grid%04d.%04d.%04d/"%(weights[1],weights[3],weights[5])
        position = startingpoint * len(gridrange) + point
        
        for j in range(savereads):
            
            thisstate = pd.read_csv(foldername+subfoldername + 'state%04d_save%04d.csv'%(startingpoint, j+1))
            temphisp  = demoEfficiency(thisstate, hispcol,      "population", democrat, republican)
            tempafram = demoEfficiency(thisstate, aframcol,     "population", democrat, republican)
            arrayDict["hispEff" ][0][position, j] = temphisp[0]*100
            arrayDict["aframEff"][0][position, j] = tempafram[0]*100
            arrayDict["baseEff" ][0][position, j] = temphisp[1]*50 + tempafram[1]*50
        print("Loaded grid%04d.%04d.%04d state %04d"%(weights[1],weights[3],weights[5], startingpoint))

#os.mkdir(foldername + "summaryFigures")
for arr in arrayDict.keys():
    for point in range(len(gridrange)):
        
        weights = gridrange[point]
        colorWeight = (0.8*np.log10(weights[1]/10)/2, 0.8*np.log10(weights[3])/2, 0.8*np.log10(weights[5])/2)
        
        for startingpoint in range(statereads):
            position = startingpoint * len(gridrange) + point
            plt.plot(arrayDict[arr][0][position,:], color = colorWeight)
    plt.title(arrayDict[arr][1])
    plt.savefig(foldername + "summaryFigures/" + arr + '.png')
    plt.clf()


#I want to view all of the trandlines on one chart.
for weights in paramList:
    #Assign a color based on these weights
    colorWeight = (0.8*np.log10(weights[1]/10)/2, 0.8*np.log10(weights[3])/2, 0.8*np.log10(weights[5])/2)
    #Load state zero for these weights
    subfoldername = "grid%04d.%04d.%04d/"%(weights[1],weights[3],weights[5])
    for i in range(numstates):
        if "maps_state%04d"%i not in os.listdir(foldername + subfoldername):
            os.mkdir(foldername + subfoldername+ "maps_state%04d"%i)
        for j in samplerate*np.arange(numreads/samplerate):
            thisstate = pd.read_csv(foldername + subfoldername + "state%04d_save%04d.csv"%(i, j+1))
            color_this_state(g, thisstate, foldername + subfoldername+ "maps_state%04d/save%04dmap.png"%(i, j), linewidth=0.3)
            print("Made map of state %d, save %d"%(i,j))




#####
#Population Evenness attempts
#####

def flatPopulationRun_Better(state, threshold = 25000, report = 10000):
    
    global adjacencyFrame, metrics
    #Prepare new state to change, and update globals
    idealpop = float(sum(blockstats.population))/ndistricts
    newstate = state.copy()
    updateGlobals(newstate)
    
    currentdiff = np.max(metrics['population']) - np.min(metrics['population'])
    freshreport = currentdiff + report
    
    while currentdiff > threshold:
        
        if currentdiff <= freshreport - report:
            print("Even-ing population.  Current range: %d"%currentdiff)
            freshreport = currentdiff
        
        #Make changes to newstate based on randomly selected district.  Extremes are more likely to be chosen.
        diffs = ((metrics.population - idealpop).abs() - float(threshold)/2).clip(0, np.inf)
        weight = diffs/sum(diffs)
        choicedist = np.random.choice(range(ndistricts), p = weight)
        #Randomly select edge on the border of maxDist
        switchedge = np.random.choice(adjacencyFrame.index[-(adjacencyFrame.isSame == 1) & (adjacencyFrame.length > 0) & ((adjacencyFrame.lowdist == choicedist) | (adjacencyFrame.highdist == choicedist))])
        if metrics.population[adjacencyFrame.lowdist[switchedge]] < metrics.population[adjacencyFrame.highdist[switchedge]]:
            smolnode     = adjacencyFrame.low[switchedge]
            biggnode     = adjacencyFrame.high[switchedge]
            tempsmoldist = adjacencyFrame.lowdist[switchedge]
            tempbiggdist = adjacencyFrame.highdist[switchedge]
        else:
            biggnode     = adjacencyFrame.low[switchedge]
            smolnode     = adjacencyFrame.high[switchedge]
            tempbiggdist = adjacencyFrame.lowdist[switchedge]
            tempsmoldist = adjacencyFrame.highdist[switchedge]
        
        biggadjacent = adjacencyFrame.ix[(adjacencyFrame.low == biggnode) | (adjacencyFrame.high == biggnode), ["high", "low", "highdist", "lowdist"]]
        proposedChanges = biggadjacent.copy()
        proposedChanges.ix[proposedChanges.low  == biggnode, "lowdist" ] = tempsmoldist
        proposedChanges.ix[proposedChanges.high == biggnode, "highdist"] = tempsmoldist
        proposedChanges.ix[:, "isSame"] = (proposedChanges.ix[:, "lowdist"] == proposedChanges.ix[:, "highdist"])
        
        #Check if this change would violate contiguousness
        neighborhood = set(biggadjacent.low).union(set(biggadjacent.high))
        proposedState = newstate.ix[neighborhood, :].copy()
        proposedState.ix[biggnode, "value"] = tempsmoldist
        
        nhadj = adjacencyFrame.ix[adjacencyFrame.low.isin(neighborhood) & adjacencyFrame.high.isin(neighborhood), ['low','high','length', 'lowdist', 'highdist']]
        oldContNeighborhood = contiguousness(newstate.loc[neighborhood], tempbiggdist, nhadj)
        
        nhadj.update(proposedChanges)
        newContNeighborhood = contiguousness(proposedState, tempbiggdist, nhadj)
        
        #If local contiguousness changes, check the whole loserDist, since it could be an annulus.
        if (oldContNeighborhood != newContNeighborhood):
            tempframe = adjacencyFrame.copy()
            tempframe.update(proposedChanges)
            tempframe.lowdist  = tempframe.lowdist.astype(int)
            tempframe.highdist = tempframe.highdist.astype(int)
            tempframe.low      = tempframe.low.astype(int)
            tempframe.high     = tempframe.high.astype(int)
            tempstate = newstate.copy()
            tempstate.value[biggnode] = tempsmoldist
            newCont = contiguousness(tempstate, tempbiggdist, tempframe)
        else:
            newCont = newContNeighborhood
        
        if newCont == 1:
            #Change everything for realz
            popchange = blockstats.population[biggnode]
            newstate.ix[biggnode, "value"] = tempsmoldist
            adjacencyFrame.ix[(adjacencyFrame.low  == biggnode), "lowdist"]  = tempsmoldist
            adjacencyFrame.ix[(adjacencyFrame.high == biggnode), "highdist"] = tempsmoldist
            adjacencyFrame.ix[:, "isSame"] = (adjacencyFrame.ix[:, "lowdist"] == adjacencyFrame.ix[:, "highdist"])
            metrics.ix[tempbiggdist, "population"] -= popchange
            metrics.ix[tempsmoldist, "population"] += popchange
        else:
            pass
            #Reject these changes, and hope for a better one on the next pass.
        
        currentdiff = np.max(metrics['population']) - np.min(metrics['population'])
    print("\n")
    #return once currentdiff is less than threshold
    
    updateGlobals(newstate)
    return newstate

begin = time.time()
final = flatPopulationRun_Better(start, threshold=20000)
end = time.time()


numstates= 1
numsteps = 100
numsaves = 500
temp = product(*([[1, 10, 100]]*3))
distinctParam = [0,1,2,3,4,5,6,7,8,9,10,11,12,15,18,19,20,21,24]

goodnessParams[1] = popDiffScore 
paramList = [[1, x[0], x[1], x[1], x[2], x[2]] for x in temp]
"""
    The way paramList works:
        - It is assumed that globalWeights is of a known length, in this case, 6.
        - I want a low, medium, and high importance option for the goodnesParams,
          but the 3rd and 4th (indices 2 and 3) I wanted to have the same weight,
          and similarly for the 5th and 6th.
        - After numstates have been created with a given system of weights,
          we can do some runs with the next set of weights.
        - This should give us some radically different maps, which we can analyze later.
"""
samplerate = 1
numreads = numsaves
#numreads = 1000

for startingpoint in range(numstates):
    for weights in [paramList[i] for i in distinctParam]:
        subfoldername = "grid%04d.%04d.%04d/"%(weights[1],weights[3],weights[5])
        if subfoldername[:-1] not in os.listdir(foldername):
            os.mkdir(foldername + subfoldername)
        
        goodnessWeights = np.array(weights)
        
        starting_state = contiguousStart()
        updateGlobals(starting_state)
        starting_state.to_csv(foldername+subfoldername + "contiguous_start%04d.csv"%startingpoint, index = False)
        metrics.to_csv(foldername + subfoldername + 'contiguous_metrics%04d.csv'%startingpoint, index = False)
        
        runningState = (flatPopulationRun_Better(starting_state, report = 25000), 1)
        updateGlobals(runningState[0])
        runningState[0].to_csv(foldername+subfoldername + "evenpop_start%04d.csv"%startingpoint, index = False)
        metrics.to_csv(foldername + subfoldername + 'evenpop_metrics%04d.csv'%startingpoint, index = False)
        
        for i in range(numsaves):
            
            runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
            
            runningState[0].to_csv(foldername+subfoldername + "state%04d_save%04d.csv"%(startingpoint, i + 1), index = False)
            
            metrics.to_csv(foldername + subfoldername + 'metrics%04d_save%04d.csv'%(startingpoint, i + 1), index = False)
            
            if i%(numsaves/20) == 0:
                print("run%04d.%04d.%04d.state%04d : %d\% complete."%(weights[1], weights[3], weights[5], startingpoint, 100*i/numsaves))
            
        print("\nFinished run%04d.%04d.%04d.state%04d\n"%(weights[1], weights[3], weights[5], startingpoint))



#####
#Supplements to cleanitallup.py
#####


#####
#Supplements to setup_stuff.py
#####


#####
#Supplements to setup.py
#####

