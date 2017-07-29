
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

foldername = "gridyesflat/"

arrayDict = {"maxBiz"        : (np.zeros((len(gridrange)*statereads,numsaves)), "Maximum Bizarreness"                          ),
             "meanBiz"       : (np.zeros((len(gridrange)*statereads,numsaves)), "Mean Bizarreness"                             ),
             "totalVar"      : (np.zeros((len(gridrange)*statereads,numsaves)), "Total Population Variation"                   ),
             "maxCont"       : (np.zeros((len(gridrange)*statereads,numsaves)), "Maximum Contiguousness"                       ),
             "maxPop"        : (np.zeros((len(gridrange)*statereads,numsaves)), "Maximum Population"                           ),
             "popDiff"       : (np.zeros((len(gridrange)*statereads,numsaves)), "Maximum Population Difference"                ),
             "hispDiff"      : (np.zeros((len(gridrange)*statereads,numsaves)), "Hispanic Boundary Difference Measure"         ),
             "aframDiff"     : (np.zeros((len(gridrange)*statereads,numsaves)), "African American Boundary Difference Measure" ),
             "baseEff"       : (np.zeros((len(gridrange)*statereads,numsaves)), "Baseline Vote Efficiency"                     ),
             "aframEff"      : (np.zeros((len(gridrange)*statereads,numsaves)), "African American Vote Efficiency"             ),
             "hispEff"       : (np.zeros((len(gridrange)*statereads,numsaves)), "Hispanic Vote Efficiency"                     ),
             "aframEffRatio" : (np.zeros((len(gridrange)*statereads,numsaves)), "African American VE/Baseline"             ),
             "hispEffRatio"  : (np.zeros((len(gridrange)*statereads,numsaves)), "Hispanic VE/Baseline"                     ),
             "goodness"      : (np.zeros((len(gridrange)*statereads,numsaves)), "Goodness"                                     )}

for startingpoint in range(statereads):
    for point in range(len(gridrange)):
        weights = gridrange[point]
        goodnessWeights = np.array(weights)
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
            temphisp  = demoEfficiency(thisstate, hispcol,  "population", democrat, republican)
            tempafram = demoEfficiency(thisstate, aframcol, "population", democrat, republican)
            arrayDict["hispEff"      ][0][position, j] = temphisp[0]*100
            arrayDict["aframEff"     ][0][position, j] = tempafram[0]*100
            arrayDict["hispEffRatio" ][0][position, j] = temphisp[0]/temphisp[1]
            arrayDict["aframEffRatio"][0][position, j] = tempafram[0]/tempafram[1]
            arrayDict["baseEff"      ][0][position, j] = temphisp[1]*50 + tempafram[1]*50
        print("Loaded grid%04d.%04d.%04d state %04d"%(weights[1],weights[3],weights[5], startingpoint))

#os.mkdir(foldername + "summaryFigures")
for arr in arrayDict.keys():
    for point in range(len(gridrange)):
        
        weights = gridrange[point]
        colorWeight = (0.8*np.log10(weights[1])/2, 0.8*np.log10(weights[3])/2, 0.8*np.log10(weights[5])/2)
        
        for startingpoint in range(statereads):
            position = startingpoint * len(gridrange) + point
            plt.plot(arrayDict[arr][0][position,:], color = colorWeight)
    plt.title(arrayDict[arr][1])
    plt.savefig(foldername + "summaryFigures/" + arr + '.png')
    plt.clf()


#####
#Population Evenness attempts
#####

def flatPopulationRun(state, threshold = 25000, report = 10000):
    
    global adjacencyFrame, metrics, mutableBlockStats
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
        diffs = (metrics.population - idealpop).abs()
        weight = diffs/sum(diffs)
        choicedist = np.random.choice(range(ndistricts), p = weight)
        
        """
        SELECTION OF SMOLDIST AND BIGGNODE
        """
        
        #Look at district boundaries
        bounds = adjacencyFrame.index[-adjacencyFrame.isSame & (adjacencyFrame.length != 0) & \
                                      ((adjacencyFrame.lowdist == choicedist ) | (adjacencyFrame.highdist == choicedist))]
        #select other district based on population difference from choicedist
        choicediff = (metrics.population[set(adjacencyFrame.lowdist[bounds]).union(set(adjacencyFrame.highdist[bounds]))] - metrics.population[choicedist]).abs()
        choiceweight = choicediff/sum(choicediff)
        otherdist = np.random.choice(choiceweight.index, p = choiceweight)
        
        #Compare sizes:
        #    Set biggdist and smoldist
        if metrics.population[choicedist] < metrics.population[otherdist]:
            tempsmoldist = choicedist
            tempbiggdist = otherdist
        else:
            tempbiggdist = choicedist
            tempsmoldist = otherdist
        
        #Nodes in biggdist where the other node is in smoldist.
        templow  = adjacencyFrame.ix[((adjacencyFrame.lowdist == tempbiggdist) & (adjacencyFrame.highdist == tempsmoldist)), "low"]
        temphigh = adjacencyFrame.ix[((adjacencyFrame.lowdist == tempsmoldist) & (adjacencyFrame.highdist == tempbiggdist)), "high"]
        
        #For each of these, find the length of the internal boundary and the length of the boundary with smoldist
        borderLands   = set(templow).union(set(temphigh))
        borderLengths = pd.DataFrame({"inner": [sum(adjacencyFrame.length[adjacencyFrame.isSame & ((adjacencyFrame.low == i) | (adjacencyFrame.high == i))]) for i in borderLands], 
                                      "outer": [sum(adjacencyFrame.length[((adjacencyFrame.low == i) | (adjacencyFrame.high    == i)) &\
                                                                          ((adjacencyFrame.lowdist == tempsmoldist) | (adjacencyFrame.highdist == tempsmoldist))]) for i in borderLands]}, 
                                     index = borderLands)
        borderLengths["proportionDiff"] = borderLengths.outer/borderLengths.inner
        #Choose randomly from the ones with the smallest number of neighbors within.
        lengthWeight = borderLengths.proportionDiff/sum(borderLengths.proportionDiff)
        
        biggnode = np.random.choice(borderLengths.index, p = lengthWeight)
        
        """
        END SELECTION OF SMOLDIST AND BIGGNODE
        """
        
        #Check if this change would violate contiguousness
        biggadjacent = adjacencyFrame.ix[((adjacencyFrame.low == biggnode) | (adjacencyFrame.high == biggnode)) & (adjacencyFrame.length != 0),["low","high", "lowdist", "highdist","isSame"]]
        proposedChanges = biggadjacent.copy()
        proposedChanges.ix[proposedChanges.low  == biggnode,  "lowdist"] = tempsmoldist
        proposedChanges.ix[proposedChanges.high == biggnode, "highdist"] = tempsmoldist
        proposedChanges.ix[:, "isSame"] = proposedChanges.lowdist == proposedChanges.highdist
        neighborhood = set(biggadjacent.low).union(set(biggadjacent.high))
        proposedState = newstate.ix[neighborhood, :].copy()
        proposedState.ix[biggnode, "value"] = tempsmoldist
        
        nhadj = adjacencyFrame.ix[(adjacencyFrame.length != 0) & (adjacencyFrame.low.isin(neighborhood) & adjacencyFrame.high.isin(neighborhood)), ['low','high','length', 'lowdist', 'highdist']]
        oldContNeighborhood = contiguousness(newstate.loc[neighborhood], tempbiggdist, nhadj)
        
        nhadj.ix[nhadj.low  == biggnode,  "lowdist"] = tempsmoldist
        nhadj.ix[nhadj.high == biggnode, "highdist"] = tempsmoldist
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
            
            #Change numedge information
            biggNewEdges   = biggadjacent.index[  biggadjacent.isSame ] #Were the same
            biggLostEdges  = biggadjacent.index[-(biggadjacent.isSame)] #Were different
            biggadjacent = adjacencyFrame.ix[(adjacencyFrame.low == biggnode) | (adjacencyFrame.high == biggnode),["low","high", "isSame"]]
            smolNewEdges  = biggadjacent.index[-(biggadjacent.isSame)] #Are no longer the same
            smolLostEdges = biggadjacent.index[  biggadjacent.isSame ] #Are now the same

            metrics.ix[tempbiggdist,'numedges'] +=\
                len(biggNewEdges) - len(biggLostEdges)
            metrics.ix[tempsmoldist,'numedges'] +=\
                len(smolNewEdges) - len(smolLostEdges)
            
            """
            #For i in neighborhood, update mutableBlockStats to the correct value.
            subMute = pd.DataFrame({"boundAdjacent": \
                      [len(adjacencyFrame.index[-adjacencyFrame.isSame & ((adjacencyFrame.low == i) | \
                                                (adjacencyFrame.high == i))])\
                       for i in neighborhood]}, index = neighborhood)
            
            # NOTE FROM MARY: SHOULD THIS BE: 
            # subMute.index = [i for i in neighborhood]
            mutableBlockStats.ix[subMute.index, "boundAdjacent"] = subMute.boundAdjacent
            """
        else:
            pass
            #Reject these changes, and hope for a better one on the next pass.
        
        currentdiff = np.max(metrics['population']) - np.min(metrics['population'])
    print("\n")
    #return once currentdiff is less than threshold
    
    updateGlobals(newstate)
    return newstate

start = contiguousStart()
begin = time.time()
final = flatPopulationRun(start, threshold=50000)
end = time.time()
#color_this_state(g, start, "tests/start.png")
color_this_state(g, final, "tests/final.png")
#color_this_state(g, newstate, "tests/dangle.png")
numstates= 1
numsteps = 100
numsaves = 250
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

#flatpop run
foldername = "gridyesflat/"
#os.mkdir(foldername)
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
        
        runningState = (flatPopulationRun(starting_state, report = 25000), 1)
        updateGlobals(runningState[0])
        runningState[0].to_csv(foldername+subfoldername + "evenpop_start%04d.csv"%startingpoint, index = False)
        metrics.to_csv(foldername + subfoldername + 'evenpop_metrics%04d.csv'%startingpoint, index = False)
        
        for i in range(numsaves):
            
            runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
            
            runningState[0].to_csv(foldername+subfoldername + "state%04d_save%04d.csv"%(startingpoint, i + 1), index = False)
            
            metrics.to_csv(foldername + subfoldername + 'metrics%04d_save%04d.csv'%(startingpoint, i + 1), index = False)
            
            if i%(numsaves/20) == 0:
                print("run%04d.%04d.%04d.state%04d : %d%% complete."%(weights[1], weights[3], weights[5], startingpoint, 100*i/numsaves))
            
        print("\nFinished run%04d.%04d.%04d.state%04d\n"%(weights[1], weights[3], weights[5], startingpoint))

#flatpop run
foldername = "gridyesflat/"
#os.mkdir(foldername)
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
        
        runningState = (flatPopulationRun(starting_state, report = 25000), 1)
        updateGlobals(runningState[0])
        runningState[0].to_csv(foldername+subfoldername + "evenpop_start%04d.csv"%startingpoint, index = False)
        metrics.to_csv(foldername + subfoldername + 'evenpop_metrics%04d.csv'%startingpoint, index = False)
        
        for i in range(numsaves):
            
            runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
            
            runningState[0].to_csv(foldername+subfoldername + "state%04d_save%04d.csv"%(startingpoint, i + 1), index = False)
            
            metrics.to_csv(foldername + subfoldername + 'metrics%04d_save%04d.csv'%(startingpoint, i + 1), index = False)
            
            if i%(numsaves/20) == 0:
                print("run%04d.%04d.%04d.state%04d : %d%% complete."%(weights[1], weights[3], weights[5], startingpoint, 100*i/numsaves))
            
        print("\nFinished run%04d.%04d.%04d.state%04d\n"%(weights[1], weights[3], weights[5], startingpoint))

#not flatpop run
foldername = "gridnoflat/"
#os.mkdir(foldername)
for startingpoint in range(numstates):
    #for weights in [paramList[i] for i in distinctParam]:
    for weights in paramList:
        subfoldername = "grid%04d.%04d.%04d/"%(weights[1],weights[3],weights[5])
        if subfoldername[:-1] not in os.listdir(foldername):
            os.mkdir(foldername + subfoldername)
        
        goodnessWeights = np.array(weights)
        
        starting_state = contiguousStart()
        updateGlobals(starting_state)
        starting_state.to_csv(foldername+subfoldername + "contiguous_start%04d.csv"%startingpoint, index = False)
        metrics.to_csv(foldername + subfoldername + 'contiguous_metrics%04d.csv'%startingpoint, index = False)
        
        runningState = (starting_state.copy(), 1)
        for i in range(numsaves):
            
            runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
            
            runningState[0].to_csv(foldername+subfoldername + "state%04d_save%04d.csv"%(startingpoint, i + 1), index = False)
            
            metrics.to_csv(foldername + subfoldername + 'metrics%04d_save%04d.csv'%(startingpoint, i + 1), index = False)
            
            if i%(numsaves/20) == 0:
                print("run%04d.%04d.%04d.state%04d : %d%% complete."%(weights[1], weights[3], weights[5], startingpoint, 100*i/numsaves))
            
        print("\nFinished run%04d.%04d.%04d.state%04d\n"%(weights[1], weights[3], weights[5], startingpoint))



#####
#Supplements to cleanitallup.py
#####


#####
#Supplements to setup_stuff.py
#####

def adjacencyGeom(connectivitydf, boundaries):
    edgeGeoms = list()
    totallength = np.shape(connectivitydf)[0]
    for i in range(totallength): 
        lo = connectivitydf.low[i]
        hi = connectivitydf.high[i]
        b1 = boundaries[lo]
        b2 = boundaries[hi]
        frontier = list()
        for b11 in b1: 
            b = len(b11)
            if b > 1:
                for b22 in b2:
                    pointsInCommon = [ [b11[a], b11[(a+1)%b] ] for a in range(b) if ((b11[a] in b22) and (b11[(a+1)%b] in b22)) ]
                    if len(pointsInCommon) > 0:
                        addThis = [x[0] for x in pointsInCommon] + [pointsInCommon[-1][1]]
                        frontier = frontier + [addThis]
        edgeGeoms.append(frontier)
        
        if (i%(totallength/50) == 0):
            print("%d%% completed with boundaries."%((i*100)/totallength))
        
    connectivitydf['geom'] = edgeGeoms
    return connectivitydf

def district_bounds_to_geom():
    global adjacencyFrame, adjacencyGeomFrame, g, perimeterNC
    whereBound = adjacencyFrame.index[-(adjacencyFrame.isSame)]
    return {"names" : [x for x in whereBound] + ["outside"],
            "paths" : [x for x in adjacencyGeomFrame.geom[whereBound]] + [[[(x[1], x[0], 0) for x in zip(perimeterNC.latitude,perimeterNC.longitude)]]],
            "xlim"  : g['xlim'],
            "ylim"  : g['ylim']}

def patch_edges(geom_to_plot, color = "black", linewidth = 1):
    #            In setup.py for each state, there should be a line like the following
    #                g = package_vtds("./VTDS_of_Interest.shp")
    #            g contains the geometries and names of the VTDS for your state space.
    #                          |
    #                          vtds_rgb_dict should be a dictionary of VTD names with RGB triples to color them.
    
    paths = geom_to_plot['paths']
    names = geom_to_plot['names']
    
    for p in range(len(paths)):
        path = paths[p]
        for subpath in path:
            if len(subpath) > 1:
                temp = [x for x in zip(*subpath)]
                plt.plot(temp[0], temp[1], color = color, linewidth = linewidth)

def color_by_rgb(geom_to_plot, vtds_rgb_dict, filename, linewidth = 1, DPI = 300, district_boundary = True, bound_args = [district_bounds_to_geom()]):
    #            In setup.py for each state, there should be a line like the following
    #                g = package_vtds("./VTDS_of_Interest.shp")
    #            g contains the geometries and names of the VTDS for your state space.
    #                          |
    #                          vtds_rgb_dict should be a dictionary of VTD names with RGB triples to color them.
    
    subset = vtds_rgb_dict.keys()
    thing = zip(geom_to_plot['paths'], geom_to_plot['names'])
    paths = [x[0] for x in thing if x[1] in subset]
    names = [x[1] for x in thing if x[1] in subset]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(geom_to_plot['xlim'])
    ax.set_ylim(geom_to_plot['ylim'])
    
    for p in range(len(paths)):
        path = paths[p]
        facecolor = vtds_rgb_dict[names[p]]
        patch = mpatches.PathPatch(path,facecolor=facecolor, edgecolor='black', linewidth=linewidth)
        ax.add_patch(patch)
    
    if district_boundary:
        patch_edges(*bound_args)
    
    ax.set_aspect(1.0)
    #plt.show()
    plt.savefig(filename, dpi=DPI)
    plt.clf()
    del fig


vtdfile = 'precinct/precinct.shp'
ds = ogr.Open(vtdfile)
lyr = ds.GetLayer(0)
vtds = features(lyr)
vtd_boundaries = boundaries(vtds, ['GEOID10'])

connectivitydf = pd.read_csv("PRECINCTconnections.csv").ix[:, ['low', 'high']]
adjacencyGeomFrame = adjacencyGeom(connectivitydf, vtd_boundaries)

district_geom = district_bounds_to_geom()

draw_these_edges(district_geom, "tests/dangle.png")

highestblack = max(blockstats.aframcon)
black = {blockstats.index[i] : (1 - blockstats.aframcon[i]/highestblack, 1 - 0.3*blockstats.aframcon[i]/highestblack, 1 - blockstats.aframcon[i]/highestblack) for i in blockstats.index}

color_by_rgb(g, black, "tests/aframmapnoedge.png", 0, district_boundary = False)
color_by_rgb(g, black, "tests/aframmapedge.png", 0)


#####
#Supplements to setup.py
#####

