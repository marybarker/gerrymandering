from matplotlib import pyplot as plt
import time
dangus = [(node, adjacencyFrame[adjacencyFrame.low == node].shape[0] + adjacencyFrame[adjacencyFrame.high == node].shape[0])\
          for node in blockstats.VTD]


dangus.count(6)

dump = [dangus.count(i) for i in range(max(dangus) + 1)]

plt.bar(range(1, len(dump) + 1), dump)

outliers = [node for node in dangus if node[1] > 15]
sum(x[1] for x in dangus)

starttime = time.clock()

goodness(starting_state)

end = time.clock()
runtime = end - starttime
