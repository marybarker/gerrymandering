import numpy as np
import random
import pandas as pd
import math
from osgeo import ogr
import os

###################################################

class Configuration():

    def __init__(self, demographics_file, connectivities_file, stateName, ndistricts, startingState = 0):
        global nPrecincts, nDistricts, precinctInfo, adjacencyFrame
        blockstats = pd.read_csv(demographics_file)
        blockstats = blockstats.drop('Unnamed: 0', 1)
        blockstats = blockstats.set_index(blockstats.VTD)
        blockstats.rename(columns={'POP100':'population'}, inplace=True)
        precinctInfo = blockstats.copy()

        self.totalPopulation = sum(precinctInfo.population)
        self.nPrecincts = len(precinctInfo.VTD)
        self.nDistricts = ndistricts
        del blockstats

        self.adjacencyFrame = pd.read_csv(connectivities_file)
        self.adjacencyFrame = self.adjacencyFrame.drop('Unnamed: 0', 1)
        #self.adjacencyFrame.rename(columns={'lo':'low'},inplace=True)
        self.adjacencyFrame.columns = ['low', 'high', 'length']

        self.runningState = 0
        self.exploration = 1000.0
        self.stateConts = []
        self.stateComps = []
        self.statePops = []

        self.precinctInfo = pd.read_csv(demographics_file)
        self.stepsTaken = 0


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
    
    def neighbor(self, state):
        
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
    
    def contiguousness(self, state, district):
    
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
    
    def perimeter(self, state, district):
        regionlist = list(state.key[state.value == district])
        return sum(self.adjacencyFrame.length[self.adjacencyFrame.low.isin(regionlist) != self.adjacencyFrame.high.isin(regionlist)])
    
    def interiorPerimeter(self, state, district):
        regionlist = list(state.key[state.value == district])
        return sum(self.adjacencyFrame.length[self.adjacencyFrame.low.isin(regionlist) & self.adjacencyFrame.high.isin(regionlist)])
    
    def distArea(self, state, district):
        regionlist = list(state.key[state.value == district])
        return sum(self.precinctInfo.ALAND[self.precinctInfo.VTD.isin(regionlist)])
    
    def population(self, state, district):
        return sum(self.precinctInfo.population[self.precinctInfo.VTD.isin(list(state.key[state.value == district]))])
    
    def efficiency(self, state, district):
        #returns difference in percentage of votes wasted.  Negative values benefit R.
        subframe = self.precinctInfo.loc[self.precinctInfo.VTD.isin(list(state.key[state.value == district]))]
        rvotes = sum(subframe['repvotes'])
        dvotes = sum(subframe['demvotes'])
        
        if rvotes > dvotes:
            wastedR = max(rvotes, dvotes) - 0.5
            wastedD = min(rvotes,dvotes)
        else:
            wastedD = max(rvotes, dvotes) - 0.5
            wastedR = min(rvotes,dvotes)
        
        return wastedR-wastedD             
    
    def compactness(self, state, district):
        outer = perimeter(state, district)     #Perimeter of district
        area = distArea(state, district)       #Area of district
        return (2*np.sqrt(np.pi*area))/outer   #Ratio of circumference of circle with same area to perimeter
    
    def compactness2(self, state):
        return sum([perimeter(state, district) for district in range(nDistricts)])/2
    
    def goodness(self, state):
        #Haves
            #contiguousness
            #evenness of population
            #efficiency
            #bizarreness
        #Needs
            #Compactness
        
        modTotalVar = sum([abs(float(x)/self.totalPopulation - float(1)/nDistricts) for x in self.statePops])/(2*(1-float(1)/ndistricts))
        
        #return -3000*abs(sum(stconts) - ndistricts) - 100*modTotalVar - 10*abs(sum(steffic)) -10*np.nansum(stbiz)
        return -300*abs(sum(self.stateConts) - self.nDistricts) - 100*modTotalVar - 10*np.nansum(self.stateComp)
    
    def switchDistrict(self, current_goodness, possible_goodness): # fix
        return float(1)/(1 + np.exp((current_goodness-possible_goodness)/self.exploration))


    def contiguousStart(self):

        state = pd.DataFrame([[self.precinctInfo.VTD[i], self.nDistricts] for i in range(self.nPrecincts)])
        state.columns = ['key', 'value']
        subAdj = self.adjacencyFrame.loc[self.adjacencyFrame.length != 0]
        
        missingdist = set(range(self.nDistricts))
        while len(list(missingdist)) > 0:
            state.value[random.randint(0,self.nPrecincts-1)] = list(missingdist)[0]
            missingdist = set.difference(set(range(self.nDistricts)), set(state['value']))
        #Above loop gives each district exactly one VTD.  The rest will be equal to ndistricts
        
        pops = [self.population(state,x) for x in range(self.nDistricts)]
        
        while self.nDistricts in set(state['value']):
            
            targdistr = pops.index(min(pops))
            
            subframe = state.loc[state.value!=self.nDistricts]
            detDists = set(subframe.key)
            tbdDists = set.difference(set(state.key), detDists)
            relevantAdjacencies = subAdj.loc[(subAdj.low.isin(detDists)) != (subAdj.high.isin(detDists))]
            #adjacencies where either low or high have a value that still has value of ndistricts, but the other doesn't
            curRegion = state.key[state.value == targdistr]
            relevantAdjacencies = subAdj.loc[((subAdj.low.isin(curRegion)) & (subAdj.high.isin(tbdDists))) |
                                             ((subAdj.high.isin(curRegion)) & (subAdj.low.isin(tbdDists)))]
            #Adjacencies where either low or high are in the region, but the other is unassigned
            
            if relevantAdjacencies.shape[0] == 0 :
                pops[targdistr] = float('inf')
            else :
                #choose entry in relevantAdjacencies and switch the value of the other node.
                temp = relevantAdjacencies.loc[relevantAdjacencies.index[random.randint(0,relevantAdjacencies.shape[0]-1)]]
                if temp.high in tbdDists:
                    state.value[state.key == temp.high] = state.value[state.key == temp.low].item()
                    pops[targdistr] = pops[targdistr] + self.precinctInfo.population[temp.high]
                else:
                    state.value[state.key == temp.low] = state.value[state.key == temp.high].item()
        #self.runningState = state
        return state

    def startingState(self, starting_state = 0):
        if starting_state != 0:
            self.runningState = pd.read_csv(starting_state)
        else:
            self.runningState = self.contiguousStart()

        temp = dict(zip(self.runningState.key, self.runningState.value))

        self.adjacencyFrame['lowdist']  = self.adjacencyFrame.low.replace(temp)
        self.adjacencyFrame['highdist'] = self.adjacencyFrame.high.replace(temp)
        self.adjacencyFrame['isSame'] = (self.adjacencyFrame.lowdist == self.adjacencyFrame.highdist)

        self.stateConts = [self.contiguousness(self.runningState, i) for i in range(self.nDistricts)]
        self.statePops  = [self.population(self.runningState, i) for i in range(self.nDistricts)]
        self.stateComps   = [self.compactness2(self.runningState, i) for i in range(self.nDistricts)]





    # initialize state space contiguously

#thing = Configuration('/home/tsugrad/Documents/gerrymandering/Pennsylvania/vtdstats.csv', '/home/tsugrad/Documents/gerrymandering/Pennsylvania/PRECINCTconnections.csv', 'PA', 18)
file1 = '/home/thisisme/Documents/NewHampshire/HarvardData/NHVTDstats.csv'
file2 = '/home/thisisme/Documents/NewHampshire/HarvardData/VTDconnections.csv'
thing = Configuration(file1, file2, 'NH', 18 )
thing.startingState()

thing.runningState



