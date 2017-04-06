from matplotlib import pyplot as plt
import time
dangus = [(node, adjacencyFrame[adjacencyFrame.low == node].shape[0] + adjacencyFrame[adjacencyFrame.high == node].shape[0])\
          for node in blockstats.VTD]


dangus.count(6)

dump = [dangus.count(i) for i in range(max(dangus) + 1)]

plt.bar(range(1, len(dump) + 1), dump)

outliers = [node for node in dangus if node[1] > 15]
sum(x[1] for x in dangus)

runtimes = np.array([float(0)]*10)

for i in range(10):
    starttime = time.clock()
    
    relevantAdjacencies = subAdj.loc[(subAdj.low.isin(subframe.key)) != (subAdj.high.isin(subframe.key))]
    temp = relevantAdjacencies.loc[relevantAdjacencies.index[random.randint(0,relevantAdjacencies.shape[0]-1)]]
    
    end = time.clock()
    runtimes[i] = end - starttime
    print(runtimes[i])

for i in range(10):
    starttime = time.clock()
    
    relevantAdjacencies = subAdj.loc[((subAdj.low.isin(state.key[state.value == targdistr])) & (subAdj.high.isin(tbdDists))) |
                                     ((subAdj.high.isin(state.key[state.value == targdistr])) & (subAdj.low.isin(tbdDists)))]

    temp = relevantAdjacencies.loc[relevantAdjacencies.index[random.randint(0,relevantAdjacencies.shape[0]-1)]]
    
    end = time.clock()
    runtimes[i] = end - starttime
    print(runtimes[i])

for kappa in range(500):
    targdistr = pops.index(min(pops))
    
    subframe = state.loc[state.value!=ndistricts]
    detDists = set(subframe.key)
    tbdDists = set.difference(set(state.key), detDists)
    relevantAdjacencies = subAdj.loc[(subAdj.low.isin(subframe.key)) != (subAdj.high.isin(subframe.key))]
    #adjacencies where either low or high have a value that still has value of ndistricts, but the other doesn't
    
    relevantAdjacencies = subAdj.loc[((subAdj.low.isin(state.key[state.value == targdistr])) & (-subAdj.high.isin(subframe.key))) |
                                     ((subAdj.high.isin(state.key[state.value == targdistr])) & (-subAdj.low.isin(subframe.key)))]
    #Adjacencies where either low or high are in the region, but the other is unassigned
    
    if relevantAdjacencies.shape[0] == 0 :
        pops[targdistr] = float('inf')
    else :
        #choose entry in relevantAdjacencies and switch the value of the other node.
        temp = relevantAdjacencies.loc[relevantAdjacencies.index[random.randint(0,relevantAdjacencies.shape[0]-1)]]
        if temp.high in tbdDists:
            state.value[state.key == temp.high] = state.value[state.key == temp.low].item()
            pops[targdistr] = pops[targdistr] + blockstats.POP100[temp.high]
        else:
            state.value[state.key == temp.low] = state.value[state.key == temp.high].item()

Class Configuration():
    ###################################################
    def MH(self, steps):
        #Takes a number of steps from the current state.
        
        current = self.runningState.copy()
        best_state = self.runningState.copy()
        current_goodness = self.goodness(current)
        best_goodness = current_goodness
        
        for i in range(steps):
            possible = self.neighbor(current)
            possible_goodness = self.goodness(possible)
            if best_goodness < possible_goodness:
                best_state = possible.copy()
                best_goodness = possible_goodness
            if random.random() < self.moveprob(current_goodness, possible_goodness):
                if current_goodness < possible_goodness:
                    self.better_hops += 1
                else:
                    self.worse_hops += 1
                current = possible.copy()
                current_goodness = possible_goodness
            else:
                self.stays += 1
            self.stepsTaken += 1
        self.runningState = best_state.copy()
        return 
    
    #######################################################################
    
    def neighbor(state):
        
        newstate = state.copy()
        missingdist = set.difference(set(range(self.ndistricts)), set(state['value']))
        #If we've blobbed out some districts, we wants to behave differently
        
        if len(missingdist) == 0:
            switchedge = np.random.choice(self.adjacencyFrame.index[-(self.adjacencyFrame.isSame == 1)])
            
            lownode  = self.adjacencyFrame.low[switchedge]
            highnode = self.adjacencyFrame.high[switchedge]
            #Randomly choose an adjacency.  Find the low node and high node for that adjacency.
            
            if random.random() < 0.5:
                newstate.value[newstate.key ==  lownode] = (newstate[newstate.key == highnode].value).item()
                checks = self.adjacencyFrame.index[((self.adjacencyFrame.low == lownode) | (self.adjacencyFrame.high == lownode)) & \
                                                   (-(self.adjacencyFrame.isSame == 1))]
                self.adjacencyFrame.isSame[((self.adjacencyFrame.low == lownode) | (self.adjacencyFrame.high == lownode)) & \
                                           self.adjacencyFrame.isSame] = False
                self.adjacencyFrame.isSame[checks] = [( newstate.value[newstate.key == self.adjacencyFrame.low[j]].item() == \
                                                        newstate.value[newstate.key == self.adjacencyFrame.high[j]].item() ) for j in checks]
            else:
                newstate.value[newstate.key == highnode] = (newstate[newstate.key ==  lownode].value).item()
                checks = self.adjacencyFrame.index[((self.adjacencyFrame.low == highnode) | (self.adjacencyFrame.high == highnode)) & \
                                                   (-(self.adjacencyFrame.isSame == 1))]
                self.adjacencyFrame.isSame[((self.adjacencyFrame.low == highnode) | (self.adjacencyFrame.high == highnode)) & \
                                           self.adjacencyFrame.isSame] = False
                self.adjacencyFrame.isSame[checks] = [( newstate.value[newstate.key == self.adjacencyFrame.low[j]].item() == \
                                                        newstate.value[newstate.key == self.adjacencyFrame.high[j]].item() ) for j in checks]
            #We want to assign both nodes the same value, and there's a 50% chance for each value being chosen.
        else:
            #If there are some districts missing, 
            changenode = newstate.key.sample(1)
            newstate.value[newstate.key == changenode] = list(missingdist)[0]
            #We want to select one randomly, and make it one of the missing districts
            self.adjacencyFrame.isSame[(self.adjacencyFrame.low == changenode) | \
                                       (self.adjacencyFrame.high == changenode)] = False
            # And none of its adjacencies match anymore.
        return newstate
    
    def contiguousness(state, district):
    
        regions = 0
        regionlist = list(state.key[state.value == district])
        if len(regionlist) == 0:
            return 1
        
        subframe = self.adjacencyFrame.loc[self.adjacencyFrame.low.isin(regionlist) & self.adjacencyFrame.high.isin(regionlist)]
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
        regionlist = list(state.key[state.value == district])
        return sum(self.adjacencyFrame.length[self.adjacencyFrame.low.isin(regionlist) != self.adjacencyFrame.high.isin(regionlist)])
    
    def interiorPerimeter(state, district):
        regionlist = list(state.key[state.value == district])
        return sum(self.adjacencyFrame.length[self.adjacencyFrame.low.isin(regionlist) & self.adjacencyFrame.high.isin(regionlist)])
    
    def distArea(state, district):
        regionlist = list(state.key[state.value == district])
        return sum(precinctInfo.ALAND[precinctInfo.VTD.isin(regionlist)])
    
    def population(state, district):
        return sum(precinctInfo.POP100[precinctInfo.VTD.isin(list(state.key[state.value == district]))])
    
    def efficiency(state, district):
        #returns difference in percentage of votes wasted.  Negative values benefit R.
        subframe = precinctInfo.loc[precinctInfo.VTD.isin(list(state.key[state.value == district]))]
        rvotes = sum(subframe['repvotes'])
        dvotes = sum(subframe['demvotes'])
        
        if rvotes > dvotes:
            wastedR = max(rvotes, dvotes) - 0.5
            wastedD = min(rvotes,dvotes)
        else:
            wastedD = max(rvotes, dvotes) - 0.5
            wastedR = min(rvotes,dvotes)
        
        return wastedR-wastedD             
    
    def compactness(state, district):
        outer = perimeter(state, district)     #Perimeter of district
        area = distArea(state, district)       #Area of district
        return (2*np.sqrt(np.pi*area))/outer   #Ratio of circumference of circle with same area to perimeter
    
    def compactness2(state):
        return sum([perimeter(state, district) for district in range(nDistricts)])/2
    
    def goodness(state):
        #Haves
            #contiguousness
            #evenness of population
            #efficiency
            #bizarreness
        #Needs
            #Compactness
        
        modTotalVar = sum([abs(float(x)/totalPopulation - float(1)/nDistricts) for x in self.statePops])/(2*(1-float(1)/ndistricts))
        
        #return -3000*abs(sum(stconts) - ndistricts) - 100*modTotalVar - 10*abs(sum(steffic)) -10*np.nansum(stbiz)
        return -300*abs(sum(self.stateConts) - self.nDistricts) - 100*modTotalVar - 10*np.nansum(self.stateComp)
    
    def switchDistrict(current_goodness, possible_goodness): # fix
        return float(1)/(1 + np.exp((current_goodness-possible_goodness)/self.exploration))
    
    def contiguousStart():
        state = pd.DataFrame([[precinctInfo.VTD[i], nDistricts] for i in range(0,self.nPrecincts)])
        state.columns = ['key', 'value']
        subAdj = self.adjacencyFrame.loc[self.adjacencyFrame.length != 0]
        
        missingdist = range(nDistricts)
        assignments = np.random.choice(state.key, nDistricts)
        state.value[state.key.isin(assignments)] = missingdist
        #Assign a single precinct to each CD.
        
        tbdDists = set(state.key)
        pops = [population(state,x) for x in range(nDistricts)]
        
        while nDistricts in set(state['value']):
            
            targdistr = pops.index(min(pops))
            
            subframe = state.loc[state.value!=nDistricts]
            
        #    relevantAdjacencies = subAdj.loc[(subAdj.low.isin(tbdDists)) != (subAdj.high.isin(tbdDists))]
            #adjacencies where either low or high have a value that still has value of ndistricts, but the other doesn't
            
            relevantAdjacencies = subAdj.loc[((subAdj.low.isin(state.key[state.value == targdistr])) & (subAdj.high.isin(tbdDists))) |
                                             ((subAdj.high.isin(state.key[state.value == targdistr])) & (subAdj.low.isin(tbdDists)))]
            #Adjacencies where either low or high are in the region, but the other is unassigned
            
            if relevantAdjacencies.shape[0] == 0 :
                pops[targdistr] = float('inf')
            else :
                #choose entry in relevantAdjacencies and switch the value of the other node.
                temp = relevantAdjacencies.loc[np.random.choice(relevantAdjacencies.index)]
                if temp.high in tbdDists:
                    state.value[state.key == temp.high] = state.value[state.key == temp.low].item()
                    pops[targdistr] = pops[targdistr] + precinctInfo.POP100[temp.high]
                    tbdDists.remove(temp.high)
                else:
                    state.value[state.key == temp.low] = state.value[state.key == temp.high].item()
                    pops[targdistr] = pops[targdistr] + precinctInfo.POP100[temp.low]
                    tbdDists.remove(temp.low)
        
        return state


















numstates = 32
numsteps  = 100
numsaves  = 25

meanBizArray  = np.zeros((numstates, numsaves))
maxBizArray   = np.zeros((numstates, numsaves))
totalVarArray = np.zeros((numstates, numsaves))
popDiffArray  = np.zeros((numstates, numsaves))

for startingpoint in range(numstates):
    
    starting_state = pd.read_csv('./startingPoints/start%d.csv'%startingpoint)
    
    temp = dict(zip(starting_state.key, starting_state.value))
    lowdists  = adjacencyFrame.low.replace(temp)
    highdists = adjacencyFrame.high.replace(temp)
    isSame = lowdists==highdists
    adjacencyFrame['isSame'] = isSame
    runningState = (starting_state.copy(), 1)
    
    color_these_states(g, [runningState], 'farsterplot/', 0)
    
    for step in range(numsaves):
        
        runningState = MH(runningState[0], numsteps, ambitiousNeighbor_old, goodness, switchDistrict)
        runningState[0].to_csv(foldername+"state%d_save%d.csv"%(startingpoint + 1, step + 1), index = False)
        
        stpops  = [population(runningState[0], i) for i in range(ndistricts)]
        stbiz   = [bizarreness(runningState[0], i) for i in range(ndistricts)]
        
        meanBizArray[startingpoint,step]  = np.mean(stbiz)
        maxBizArray[startingpoint,step]   = np.max(stbiz)
        totalVarArray[startingpoint,step] = np.sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in stpops])/(2*(1-float(1)/ndistricts))
        popDiffArray[startingpoint,step]  = np.max(stpops) - min(stpops)
        
        print("Saved " + foldername+"state%d_save%d.csv"%(startingpoint + 1, step + 1))
    
    color_these_states(g, [runningState], 'farsterplot/', 2500)

