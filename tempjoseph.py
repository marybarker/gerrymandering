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



