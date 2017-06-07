import os
import time
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

#os.chdir('/Users/marybarker/Documents/tarleton_misc/gerrymandering/Pennsylvania')
#os.chdir('/home/odin/Documents/gerrymandering/gerrymandering/Pennsylvania')
#os.chdir('/home/odin/Documents/gerrymandering/gerrymandering/Texas')
#os.chdir('/home/odin/Documents/gerrymandering/gerrymandering/NorthCarolina')

execfile('../cleanitallup.py')
execfile('../setup_stuff.py')
execfile('setup.py') #Stack overflow doesn't like this, for the record.
execfile('../tempjoseph.py')
#%run -i simulator #Supposedly Stack overflow is okay with this, maybe?

metrics = pd.DataFrame()

foldername = "fffffff2/"
foldername = "slambp3ALLOFTHESTATES/"
foldername = "muffle/" # even when global metrics are incorrectly updated, we keep the incorrect version
foldername = "huffle/" # reset global metrics after every MH call
foldername = "buffle/" # low to high or high to low
foldername = "boundarydangle/"
foldername = "awnw/"
#os.mkdir(foldername)

numstates= 1
numsteps = 200
numsaves = 10

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
#Keep track of metrics
#########

maxBizArray = np.zeros((numstates,numsaves))
meanBizArray = np.zeros((numstates,numsaves))
totalVarArray = np.zeros((numstates,numsaves))
maxContArray = np.zeros((numstates,numsaves))
maxPopArray = np.zeros((numstates,numsaves))
popDiffArray = np.zeros((numstates,numsaves))
hispDiffArray = np.zeros((numstates,numsaves))
aframDiffArray = np.zeros((numstates,numsaves))

overallGoodnessArray = np.zeros((numstates,numsaves))

#States as they are being created

# ...

#States in folder
for startingpoint in range(numstates):
    for j in range(numsaves):
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

#####
#Make maps of states
#####

samplerate = 1

for i in range(numstates):
    if "maps_state%04d"%i not in os.listdir(foldername):
        os.mkdir(foldername + "maps_state%04d"%i)
    for j in samplerate*np.arange(numsaves/samplerate):
        thisstate = pd.read_csv(foldername + "state%d_save%d.csv"%(i, j+1) ,dtype={'key':str})
        color_this_state(g, thisstate, foldername + "maps_state%04d/save%dmap.png"%(i, j), linewidth=0.3)

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

