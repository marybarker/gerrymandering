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
    
    color_these_states(g, [runningState], foldername, step+1)
    
    end = time.clock()
    runtimes[i] = end - starttime
    print(runtimes[i])

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

