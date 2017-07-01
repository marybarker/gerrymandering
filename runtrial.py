import os
import time
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from itertools import product

#os.chdir('/Users/marybarker/Documents/tarleton_misc/gerrymandering/Pennsylvania')
#os.chdir('/Users/marybarker/Documents/tarleton_misc/gerrymandering/NorthCarolina/')
#os.chdir('/home/odin/Documents/gerrymandering/gerrymandering/Pennsylvania')
#os.chdir('/home/odin/Documents/gerrymandering/gerrymandering/Texas')
#os.chdir('/home/odin/Documents/gerrymandering/gerrymandering/NorthCarolina')

execfile('../cleanitallup.py')
execfile('../setup_stuff.py')
execfile('setup.py') #Stack overflow doesn't like this, for the record.
#execfile('../tempjoseph.py')
#%run -i simulator #Supposedly Stack overflow is okay with this, maybe?

metrics = pd.DataFrame()

foldername = "fffffff2/"
foldername = "slambp3ALLOFTHESTATES/"
foldername = "muffle/" # even when global metrics are incorrectly updated, we keep the incorrect version
foldername = "huffle/" # reset global metrics after every MH call
foldername = "buffle/" # low to high or high to low
foldername = "boundarydangle/"
foldername = "awnw/"
foldername = "gridruns/"
#os.mkdir(foldername)

numstates= 1
numsteps = 100
numsaves = 1000
samplerate = 1
numreads = numsaves
#numreads = 734

#########
#Run numstates instances from scratch, without annealing
#########

for startingpoint in range(numstates):
    
    starting_state = contiguousStart()
    runningState = (starting_state.copy(), 1)
    updateGlobals(runningState[0])
    for i in range(numsaves):
        
        runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
        
        runningState[0].to_csv(foldername+"state%d_save%d.csv"%(startingpoint, i + 1), index = False)
        
        metrics.to_csv(foldername + 'metrics%d_save%d.csv'%(startingpoint, i+1), index = False)
        
        print("Written to state%d_save%d.csv"%(startingpoint, i + 1))


#########
#Run numstates instances FOR EACH set of parameters from scratch, without annealing
#########

numstates= 5
numsteps = 100
numsaves = 1000
temp = product(*([[1, 10, 100]]*3))
distinctParam = [0,1,2,3,4,5,6,7,8,9,10,11,12,15,18,19,20,21,24]
paramList = [[1, 10* x[0], x[1], x[1], x[2], x[2]] for x in temp]
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

for weights in [paramList[i] for i in distinctParam]:
    subfoldername = "grid%04d.%04d.%04d/"(weights[1],weights[3],weights[5])
    if subfoldername not in os.listdir(foldername):
        os.mkdir(foldername + subfoldername)
    
    goodnessWeights = np.array(weights)
    
    for startingpoint in range(numstates):
        
        starting_state = contiguousStart()
        runningState = (starting_state.copy(), 1)
        updateGlobals(runningState[0])
        for i in range(numsaves):
            
            runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
            
            runningState[0].to_csv(foldername+subfoldername + "state%04d_save%04d.csv"%(startingpoint, i + 1), index = False)
            
            metrics.to_csv(foldername + subfoldername + 'metrics%04d_save%04d.csv'%(startingpoint, i + 1), index = False)
            
            print(str(i))
            
        print("\nFinished run%04d.%04d.%04d.state%04d\n"%(weights[1], weights[3], weights[5], startingpoint))

"""
tempi = i
tempstartingpoint = startingpoint+1
tempweights = weights

tempstartingpoint = 2
tempweights = [1, 10, 10, 10, 100, 100]
"""

"""
#For restarting inner loop from specified save range.
for i in range(tempi, numsaves):
    
    runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
    
    runningState[0].to_csv(foldername+"run%04d.%04d.%04d.state%04d_save%04d.csv"%(weights[1], weights[3], weights[5], startingpoint, i + 1), index = False)
    
    metrics.to_csv(foldername + 'run%04d.%04d.%04d.metrics%04d_save%04d.csv'%(weights[1], weights[3], weights[5], startingpoint, i + 1), index = False)
    
    print(str(i))
    
print("\nFinished run%04d.%04d.%04d.state%04d\n"%(weights[1], weights[3], weights[5], startingpoint))
"""

"""
#For restarting middle loop from specified state number
for startingpoint in range(tempstartingpoint, numstates):
    
    starting_state = contiguousStart()
    runningState = (starting_state.copy(), 1)
    updateGlobals(runningState[0])
    for i in range(numsaves):
        
        runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
        
        runningState[0].to_csv(foldername+"run%04d.%04d.%04d.state%04d_save%04d.csv"%(weights[1], weights[3], weights[5], startingpoint, i + 1), index = False)
        
        metrics.to_csv(foldername + 'run%04d.%04d.%04d.metrics%04d_save%04d.csv'%(weights[1], weights[3], weights[5], startingpoint, i + 1), index = False)
        
        print(str(i))
        
    print("\nFinished run%04d.%04d.%04d.state%04d\n"%(weights[1], weights[3], weights[5], startingpoint))

"""

"""
#For restarting outer loop from specified weight
for weights in paramList[14:]:
    
    goodnessWeights = np.array(weights)
    
    for startingpoint in range(numstates):
        
        starting_state = contiguousStart()
        runningState = (starting_state.copy(), 1)
        updateGlobals(runningState[0])
        for i in range(numsaves):
            
            runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
            
            runningState[0].to_csv(foldername+"run%04d.%04d.%04d.state%04d_save%04d.csv"%(weights[1], weights[3], weights[5], startingpoint, i + 1), index = False)
            
            metrics.to_csv(foldername + 'run%04d.%04d.%04d.metrics%04d_save%04d.csv'%(weights[1], weights[3], weights[5], startingpoint, i + 1), index = False)
            
            print(str(i))
            
        print("\nFinished run%04d.%04d.%04d.state%04d\n"%(weights[1], weights[3], weights[5], startingpoint))

"""
#########
#Read metrics and make graphs of stats for the run.
#########
"""
tempweights = paramList[5]
extendedname = foldername + "run%04d.%04d.%04d."%(tempweights[1], tempweights[3], tempweights[5])
plotMetricsByState(createMetricsArrays(extendedname, numstates, numreads, samplerate, pad = True))
"""
plotMetricsByState(createMetricsArrays(foldername, numstates, numreads, samplerate))

#########
#Keep track of metrics
#########
"""
maxBizArray = np.zeros((numstates,numreads))
meanBizArray = np.zeros((numstates,numreads))
totalVarArray = np.zeros((numstates,numreads))
maxContArray = np.zeros((numstates,numreads))
maxPopArray = np.zeros((numstates,numreads))
popDiffArray = np.zeros((numstates,numreads))
hispDiffArray = np.zeros((numstates,numreads))
aframDiffArray = np.zeros((numstates,numreads))

overallGoodnessArray = np.zeros((numstates,numreads))

#States as they are being created

# ...

#States in folder
for startingpoint in range(numstates):
    for j in range(numreads):
        #tempstate = pd.read_csv(foldername + "state%d_save%d.csv"%(i, j+1))
        #updateGlobals(tempstate)
        #pd.DataFrame(metrics).to_csv(foldername + 'metrics%d_save%d.csv'%(1, i+50), index = False)
        thismetrics = pd.read_csv(foldername+'metrics%d_save%d.csv'%(startingpoint, j+1))
        
        mindists = thismetrics['mincon'].argsort()[-numMajMinDists:][::-1]

        meanBizArray[startingpoint,j]         = np.mean(thismetrics['bizarreness'])
        maxBizArray[startingpoint,j]          = np.max(thismetrics['bizarreness'])
        maxContArray[startingpoint,j]         = np.max(thismetrics['contiguousness'])
        maxPopArray[startingpoint,j]          = np.max(thismetrics['population'])
        popDiffArray[startingpoint,j]         = np.max(thismetrics['population']) - np.min(thismetrics['population'])
        hispDiffArray[startingpoint,j]        = np.sum(thismetrics['sumHispDiff'][mindists])
        aframDiffArray[startingpoint,j]       = np.sum(thismetrics['sumAframDiff'][mindists])
        totalVarArray[startingpoint,j]        = np.sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in thismetrics['population']])/(2*(1-float(1)/ndistricts))
        overallGoodnessArray[startingpoint,j] = goodness(thismetrics)

        #metrics = {'contiguousness': metrics['contiguousness'],
        #           'population'    : stPops,
        #           'bizarreness'   : stBiz,
        #           'perimeter'     : stPerim,
        #           'area'          : stArea}
    
    print("Stored metrics for state %d"%(startingpoint))
"""
#####
#Show Plots (Not Save)
#####
"""
startingpoint = 0

plt.plot(meanBizArray[startingpoint,:])
plt.title('mean Biz')
plt.show()
plt.clf()

plt.plot(maxBizArray[startingpoint,:])
plt.title('max Biz')
plt.show()
plt.clf()

plt.plot(maxContArray[startingpoint,:])
plt.title('max contig')
plt.show()
plt.clf()

plt.plot(maxPopArray[startingpoint,:])
plt.title('max pop')
plt.show()
plt.clf()

plt.plot(popDiffArray[startingpoint,:])
plt.title('pop diff')
plt.show()
plt.clf()

plt.plot(totalVarArray[startingpoint,:])
plt.title('mean Pop')
plt.show()
plt.clf()

plt.plot(aframDiffArray[startingpoint,:])
plt.title('Afram Diff')
plt.show()
plt.clf()

plt.plot(hispDiffArray[startingpoint,:])
plt.title('Hisp Diff')
plt.show()
plt.clf()

plt.plot(overallGoodnessArray[startingpoint,:])
plt.title('goodness')
plt.show()
plt.clf()
"""
#####
#Save Plots
#####

"""
startingpoint = 0

plt.plot(meanBizArray[startingpoint,:])
plt.title('mean Biz')
plt.savefig(foldername+'meanBiz%04d.png'%(startingpoint))
plt.clf()

plt.plot(maxBizArray[startingpoint,:])
plt.title('max Biz')
plt.savefig(foldername+'maxBiz%04d.png'%(startingpoint))
plt.clf()

plt.plot(maxContArray[startingpoint,:])
plt.title('max contig')
plt.savefig(foldername+'maxCont%04d.png'%(startingpoint))
plt.clf()

plt.plot(maxPopArray[startingpoint,:])
plt.title('max pop')
plt.savefig(foldername+'maxPop%04d.png'%(startingpoint))
plt.clf()

plt.plot(popDiffArray[startingpoint,:])
plt.title('pop diff')
plt.savefig(foldername+'popDiff%04d.png'%(startingpoint))
plt.clf()

plt.plot(totalVarArray[startingpoint,:])
plt.title('pop variation')
plt.savefig(foldername+'popVar%04d.png'%(startingpoint))
plt.clf()

plt.plot(aframDiffArray[startingpoint,:])
plt.title('Afram Diff')
plt.savefig(foldername+'aframDiff%04d.png'%(startingpoint))
plt.clf()

plt.plot(hispDiffArray[startingpoint,:])
plt.title('Hisp Diff')
plt.savefig(foldername+'hispDiff%04d.png'%(startingpoint))
plt.clf()

plt.plot(overallGoodnessArray[startingpoint,:])
plt.title('goodness')
plt.savefig(foldername+'goodness%04d.png'%(startingpoint))
plt.clf()
"""

#####
#Make maps of states
#####

samplerate = numsaves/250
"""
for i in range(numstates):
    if "maps_state%04d"%i not in os.listdir(foldername):
        os.mkdir(foldername + "maps_state%04d"%i)
    for j in samplerate*np.arange(numreads/samplerate):
        thisstate = pd.read_csv(foldername + "state%d_save%d.csv"%(i, j+1))
        color_this_state(g, thisstate, foldername + "maps_state%04d/save%04dmap.png"%(i, j), linewidth=0.3)
        print("Made map of state %d, save %d"%(i,j))
"""
for weights in [paramList[i] for i in distinctParam]:
    subfoldername = "grid%04d.%04d.%04d/"%(weights[1],weights[3],weights[5])
    for i in range(numstates):
        if "maps_state%04d"%i not in os.listdir(foldername + subfoldername):
            os.mkdir(foldername + subfoldername+ "maps_state%04d"%i)
        for j in samplerate*np.arange(numreads/samplerate):
            thisstate = pd.read_csv(foldername + subfoldername + "state%04d_save%04d.csv"%(i, j+1))
            color_this_state(g, thisstate, foldername + subfoldername+ "maps_state%04d/save%04dmap.png"%(i, j), linewidth=0.3)
            print("Made map of state %d, save %d"%(i,j))
#####

#Plots of all states
#####

for i in range(numstates):
    plt.plot(meanBizArray[i,:])
plt.title('mean Biz')
plt.show()
plt.clf()
for i in range(numstates):
    plt.plot(maxBizArray[i,:])
plt.title('max Biz')
plt.show()
plt.clf()
for i in range(numstates):
    plt.plot(maxContArray[i,:])
plt.title('max contig')
plt.show()
plt.clf()
for i in range(numstates):
    plt.plot(maxPopArray[i,:])
plt.title('max pop')
plt.show()
plt.clf()
for i in range(numstates):
    plt.plot(popDiffArray[i,:])
plt.title('pop diff')
plt.show()
plt.clf()
for i in range(numstates):
    plt.plot(totalVarArray[i,:])
plt.title('mean Pop')
plt.show()
plt.clf()
for i in range(numstates):
    plt.plot(overallGoodnessArray[i,:])
plt.title('goodness')
plt.show()

#########
#Determine efficiency gaps of states.
#########
demo1 = "DEM_C"
demo2 = "REP_C"

#States as they are being created
efficiencyGapArray = np.zeros(numstates)
gapArray = np.zeros((numstates, ndistricts))
popArray = np.zeros((numstates, ndistricts))
for i in range(numstates):
    state = contiguousStart()
    state.to_csv(foldername + "state%d_start.csv"%(i), index = False)
    temp = [demoEfficiency(state, dist, demo1, demo2) for dist in range(ndistricts)]
    gapArray[i,:] = [x[0] - x[1] for x in temp]
    efficiencyGapArray[i] = np.sum(gapArray[i,:])
    print(i)

#States in folder
for i in range(numstates):
    state = pd.read_csv(foldername + "state%d_start.csv"%(i))
    temp = [demoEfficiency(state, dist, demo1, demo2) for dist in range(ndistricts)]
    gapArray[i,:] = [x[0] - x[1] for x in temp]
    efficiencyGapArray[i] = np.sum(gapArray[i,:])
    
    popArray[i, :] = [population(state, dist) for dist in range(ndistricts)]

plt.hist(efficiencyGapArray)

