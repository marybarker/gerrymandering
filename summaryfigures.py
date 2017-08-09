import pandas as pd
import os
import numpy as np
import sys

"""
Turn this into a script that accepts arguments at some point.  For now, just run it.

"""

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

if "summaryFigures" not in os.listdir(foldername):
    os.mkdir(foldername + "summaryFigures")

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
