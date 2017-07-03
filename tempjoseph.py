
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
#Supplements to cleanitallup.py
#####

def flatPopulationRun_Bad(state, threshold = 25000):
    
    global adjacencyFrame, metrics
    #Prepare new state to change, and update globals
    newstate = state.copy()
    updateGlobals(newstate)
    
    currentdiff = float(max(0, (np.max(metrics['population']) - np.min(metrics['population']))))
    
    while currentdiff > threshold:
        
        print(currentdiff)
        
        #Make changes to newstate based on maximum and minimum populations
        maxDist = metrics.index[metrics.population == max(metrics.population)].values[0]
        minDist = metrics.index[metrics.population == min(metrics.population)].values[0]
        
        #Randomly select something on the border of maxDist
        switch = adjacencyFrame.loc[np.random.choice(adjacencyFrame[(adjacencyFrame.lowdist == maxDist) != (adjacencyFrame.highdist == maxDist)].index, 1)]
        if (switch.lowdist == maxDist).values[0]:
            maxnode   = switch.low.values[0]
            othernode = switch.high.values[0]
        else:
            maxnode   = switch.high.values[0]
            othernode = switch.low.values[0]
        #Switch the maxDist node to the other district
        popchange = blockstats.population[maxnode]
        newdist = state.value[othernode]
        newstate.ix[maxnode, "value"] = newdist
        adjacencyFrame.ix[(adjacencyFrame.low  == maxnode), "lowdist"]  == newdist
        adjacencyFrame.ix[(adjacencyFrame.high == maxnode), "highdist"] == newdist
        metrics.ix[maxDist, "population"] -= popchange
        metrics.ix[newdist, "population"] += popchange
        
        #Repeat for minDist, but different.
        switch = adjacencyFrame.loc[np.random.choice(adjacencyFrame[(adjacencyFrame.lowdist == minDist) != (adjacencyFrame.highdist == minDist)].index, 1)]
        if (switch.lowdist == maxDist).values[0]:
            minnode   = switch.low.values[0]
            othernode = switch.high.values[0]
        else:
            minnode   = switch.high.values[0]
            othernode = switch.low.values[0]
        #Switch the maxDist node to the other district
        popchange = blockstats.population[othernode]
        newdist = minDist
        olddist = state.value[othernode]
        newstate.ix[othernode, "value"] = newdist
        adjacencyFrame.ix[(adjacencyFrame.low  == othernode), "lowdist"]  == newdist
        adjacencyFrame.ix[(adjacencyFrame.high == othernode), "highdist"] == newdist
        metrics.ix[minDist, "population"] += popchange
        metrics.ix[olddist, "population"] -= popchange
        
        currentdiff = float(max(0, (np.max(metrics['population']) - np.min(metrics['population']))))
    #return once currentdiff is less than threshold
    return newstate

#####
#Supplements to setup_stuff.py
#####

#####
#Supplements to setup.py
#####
