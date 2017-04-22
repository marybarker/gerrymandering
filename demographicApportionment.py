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
groupFile   = '/home/joseph/Dropbox/gerrymandering/PA_BlockGroupData.csv'
vtdInfoFile = './vtdstats.csv'

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
    for j in range(len(vtds)):
        if blockCentroids[i].Intersects(vtds[j].geometry()):
            blockVTD.append(vtds[j].GEOID10 + vtds[j].NAME10)
            break
    if i%1000 ==0:
        print("We've gone through %d blocks!"%i)

blockFrame = pd.DataFrame([(block.BLOCKID10, block.BLOCKID10[:-3], block.POP10) for block in blocks], columns = ('BLOCKID', 'GROUPID', 'population'))

precinctConnectionsFile = 'PRECINCTconnections.csv'
g = package_vtds()

ds = ogr.Open(precinctBoundaryFile)
lyr = ds.GetLayer(0)
precincts = features(lyr)




