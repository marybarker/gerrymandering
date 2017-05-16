from osgeo import ogr
import numpy as np
import os
import pandas as pd
import math

""" * * * * * * * * * * * * * * * * * """
"""   make a  list of all features    """
""" * * * * * * * * * * * * * * * * * """
def features(layer):
    features = []
    for feat in layer: 
        features.append( feat )
    return features

blockFile   = './tabblock2010_42_pophu/tabblock2010_42_pophu.shp' 
vtdFile     = './precinct/precinct.shp'
#groupFile   = '/home/joseph/Dropbox/gerrymandering/PA_BlockGroupData.csv'
groupFile   = './PA_BlockGroupData.csv'
vtdInfoFile = './vtdstats.csv'
assignmentfile = './blockframe.csv'
#blockFrame = pd.read_csv(assignmentfile)

ds = ogr.Open(blockFile)
lyr = ds.GetLayer(0)
blocks = features(lyr)

ds = ogr.Open(vtdFile)
lyr = ds.GetLayer(0)
vtds = features(lyr)

groupStats = pd.read_csv(groupFile)
del groupStats['Unnamed: 0']

vtdInfo = pd.read_csv(vtdInfoFile)
del vtdInfo['Unnamed: 0']

def coordDist(c1, c2):
    return np.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2)

def closestkvtds(p1, k=100):
    dists = [coordDist(p1,vtdCentroids[i]) for i in range(len(vtdCentroids))]
    return np.argpartition(dists, 10)[:10]


blockCentroids = [block.geometry().Centroid() for block in blocks]
blockCentroidPoints = [x.GetPoint() for x in blockCentroids]
vtdCentroids   = [  vtd.geometry().Centroid().GetPoint() for   vtd in   vtds]

blockVTD = []
for i in range(len(blocks)):
    assignedVTD = False
    for j in range(len(vtds)):
        if blockCentroids[i].Intersects(vtds[j].geometry()):
            blockVTD.append((blocks[i].BLOCKID10, vtds[j].GEOID10 + vtds[j].NAME10, i))
            assignedVTD = True
            break
    if not assignedVTD:
        blockVTD.append((blocks[i].BLOCKID10, 'NO VTD COVERS CENTROID', i))
    if i%1000 ==0:
        print("We've gone through %d blocks!"%i)

blockFrame = pd.DataFrame([(block.BLOCKID10, block.BLOCKID10[:-3], block.POP10) for block in blocks], columns = ('BLOCKID', 'GROUPID', 'population'))
blockFrame['VTDID'] = [x[1] for x in blockVTD]
blockFrame.to_csv(assignmentFile, index=False)

problems = blockFrame.loc[blockFrame.VTDID == "NO VTD COVERS CENTROID"].copy()
closestVTD = [] #This currently looks at the nearest centroid.
                #What we actually want is the VTD that most intersects the block
                #Check all VTDs for whether they intersect.  Then, among those, 
                # find the area of the intersection.  Assign accordingly.
for ID in problems.BLOCKID:
    blockShape = [block for block in blocks if block.BLOCKID10 == ID][0]
    closest = vtds[closestkvtds(blockShape.geometry().Centroid().GetPoint(), 1)[0]]
    closestVTD.append(closest.GEOID10 + closest.NAME10)


precinctConnectionsFile = 'PRECINCTconnections.csv'
g = package_vtds()

ds = ogr.Open(precinctBoundaryFile)
lyr = ds.GetLayer(0)
precincts = features(lyr)

#VTD block group pairings.
# VTD | block group | % of block group

#blockFrame = pd.read_csv(assignmentfile)
#groupStats = pd.read_csv(groupFile)

vtdgrouppairing = { "VTD"      : [],
                    "GROUP"    : [],
                    "PCTGROUP" : [] }

i = 0
numgroups = groupStats.shape[0]
grouplist = list(set(blockFrame.GROUPID))

for group in grouplist:
    temp = blockFrame.ix[blockFrame.GROUPID == group, ["population", "VTDID"]]
    grouppop = sum(temp.population)

    if grouppop != 0:
        for vtd in temp.VTDID.unique():
            vtdgrouppairing['VTD'].append(vtd)
            vtdgrouppairing['GROUP'].append(group)
            vtdgrouppairing['PCTGROUP'].append(float(sum(temp.population[temp.VTDID == vtd]))/grouppop)
    else:
        for vtd in temp.VTDID.unique():
            vtdgrouppairing['VTD'].append(vtd)
            vtdgrouppairing['GROUP'].append(group)
            vtdgrouppairing['PCTGROUP'].append(float(1)/len(temp.VTDID.unique()))
    i = i + 1
    print("Group %d of %d"%(i, numgroups))

#vpgtlist = list(set([x for x in zip(vtdgrouppairing['VTD'], vtdgrouppairing['GROUP'], vtdgrouppairing['PCTGROUP'])]))
#vgpFrame = pd.DataFrame(vpgtlist, columns=["VTD", "GROUP", "PCTGROUP"])
vgpFrame = pd.DataFrame(vtdgrouppairing)

problems = [vgpFrame.GROUP[i] for i in vgpFrame.index if \
            abs(sum(vgpFrame.PCTGROUP[vgpFrame.GROUP == vgpFrame.GROUP[i]]) - 1) > 0.001]

vgpFrame.ix[ vgpFrame.GROUP == problems[0] ]



def applyDemographicInfo(vtds, listOfStatsToCompute, VTDToBlockGroup, blockGroupValues, outputfilename):
    blockGroupValues = blockGroupValues.set_index("GEOID10")
    thewholething = pd.DataFrame(columns=listOfStatsToCompute)

    for vtd in vtds:
        tempRow = pd.DataFrame(dict(zip(listOfStatsToCompute, [0 for x in listOfStatsToCompute])), index=[1])
        stufftomultiply = VTDToBlockGroup.ix[VTDToBlockGroup.VTD == vtd, ['PCTGROUP', 'GROUP']]
        
        for group in stufftomultiply.GROUP:
            tempRow += blockGroupValues.ix[group, listOfStatsToCompute].values*stufftomultiply.ix[stufftomultiply.GROUP == group, 'PCTGROUP'].values
        thewholething = thewholething.append(tempRow)

    thewholething['VTD'] = vtds
    thewholething.to_csv(outputfilename)

applyDemographicInfo(startingState, ['B03002m1', 'C17002m1'], vgpFrame, groupStats, "all_the_apportionments.csv")
columns = [
            'B02001e1',  # 
            'B02008e1',  # 
            'B02009e1',  # 
            'B02010e1',  # 
            'B02011e1',  # 
            'B02012e1',  # 
            'B02013e1',  # 
            'B03002e1',  # 
            'B03003e1',  # 
            'B17017e1',  # 
            'B17021e1',  # 
            'B19301e1',  # 
            'C17002e1',  # 
          ]
thing = blockstats.VTD.copy()
applyDemographicInfo(thing, columns, vgpFrame, groupStats, "all_the_apportionments_try_2.csv")



