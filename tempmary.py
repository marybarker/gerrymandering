import os
import time
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import math
from osgeo import ogr
os.chdir('/Users/marybarker/Documents/tarleton_misc/gerrymandering/Pennsylvania')
foldername = 'attempt1/'


# create new shapefile with merged vtds
blockstats = pd.read_csv("vtdstats.csv")
adjacencyFrame = pd.read_csv("PRECINCTconnections.csv")
del adjacencyFrame['Unnamed: 0']
adjacencyFrame.columns = ['low','high','length']

# first merge all doughnut-ed vtds into single vtd
singles = [vtd for vtd in blockstats.VTD if adjacencyFrame.loc[(adjacencyFrame.low == vtd) | (adjacencyFrame.high == vtd)].shape[0] == 1]

toglom = []
for single in singles: 
    if any(adjacencyFrame.high == single):
        toglom.append(adjacencyFrame.low[adjacencyFrame.high == single].item())
    else:
        toglom.append(adjacencyFrame.high[adjacencyFrame.low == single].item())

"""
        mynewdf = blockstats.copy()
        thingstoadd = list(mynewdf.columns.copy())
        thingstoadd.remove('VTD')
        thingstoadd.remove('PERIM')
        
        for i in range(len(singles)):
            mynewdf.ix[toglom[i], thingstoadd] += mynewdf.ix[singles[i], thingstoadd]
        mynewdf = mynewdf.ix[~(mynewdf.index.isin(singles))]
        
        apportionmentdata = pd.read_csv('all_the_apportionments_try_3.csv')
        del apportionmentdata["Unnamed: 0"]
        apportionmentdata = apportionmentdata.set_index("VTD")
        thingstoadd = list(apportionmentdata.columns.copy())
        
        for i in range(len(singles)):
            apportionmentdata.ix[toglom[i], thingstoadd ] += apportionmentdata.ix[singles[i], thingstoadd]
        
        blockstats = pd.read_csv("noIslandsVTDStats.csv")
        blockstats = blockstats.drop("Unnamed: 0", 1)
        apportionmentdata = apportionmentdata.ix[~(apportionmentdata.index.isin(singles))]
        apportionmentdata['VTD'] = apportionmentdata.index.copy()
        apportionmentdata.to_csv("noIslandsApportionmentDataTry2.csv")
        newdata = apportionmentdata.merge(blockstats, how = 'outer', on='VTD')
"""

# now read the shapefile and merge together desired ones. 

allthevtds = ogr.Open('precinct/precinct.shp')
lyr = allthevtds.GetLayer(0)
vtds = [feat for feat in lyr]

glommers = []
glommees = []
indices = []
sindices = []
for count in range(len(vtds)):
    vtd = vtds[count]
    name = vtd.GEOID10+vtd.NAME10
    if name in toglom:
        indices.append(count)
        glommers.append(vtd)
    if name in singles:
        sindices.append(count)
        glommees.append(vtd)

allimportantbs = boundaries( [vtds[x] for x in indices])
singlesbds     = boundaries( [vtds[x] for x in sindices])

bettervtds = [vtds[i] for i in range(len(vtds)) if (i not in sindices and i not in indices)]

thingstoaddup = ['ALAND10', 'AWATER10', 'POP100']

for thing in glommers:
    g = thing.geometry()
    name = thing.GEOID10+thing.NAME10
    othernames = [singles[x] for x in range(len(toglom)) if toglom[x] == name]
    for other in glommees:
        if other.GEOID10+other.NAME10 in othernames:
            g = g.Union(other.geometry())
            for key in thingstoaddup:
                thing[key] += other[key]
    thing.SetGeometry(g)

    bettervtds.append(thing)

outShapefile = str(os.getcwd()) + '/'+foldername+"with_empty_polys.shp"
outDriver = ogr.GetDriverByName("ESRI Shapefile")
if os.path.exists(outShapefile):
    outDriver.DeleteDataSource(outShapefile)
# Create the output shapefile
outDataSource = outDriver.CreateDataSource(outShapefile)
outLayer = outDataSource.CreateLayer('vtds', geom_type=ogr.wkbPolygon)

for key in vtds[0].keys():
    thisfield = ogr.FieldDefn(key, ogr.OFTInteger)
    outLayer.CreateField(thisfield)

for vtd in bettervtds:
    featureDefn = outLayer.GetLayerDefn()
    g = vtd.geometry()
    geoid = vtd.GEOID10
    name = vtd.NAME10
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(g)

    for key in vtd.keys():
        feature.SetField(key, vtd[key])
    outLayer.CreateFeature(feature)
    feature = None

# Close DataSource
outDataSource.Destroy()








#centroids
centroids = []
names = []
for vtd in vtds:

    g = vtd.geometry()
    if (vtd.GEOID10+vtd.NAME10) in toglom:
        other = vtd.geometry()
        for others in vtds:
            if (others.GEOID10+others.NAME10) in singles:
                if toglom[singles.index(others.GEOID10+others.NAME10)] == (vtd.GEOID10+vtd.NAME10):
                    other = others.geometry()
                    g = g.Union(other)
        centroid = g.Centroid().GetPoint()
    elif (vtd.GEOID10+vtd.NAME10) in singles:
        pass
    else:
        centroid = g.Centroid().GetPoint()

    names.append(vtd.GEOID10+vtd.NAME10)
    centroids.append(centroid)

pd.DataFrame({"VTD":names, "centroid":centroids}).to_csv("noIslandsCentroids.csv")


############################################################################################################
#                                Contested districts: 23, 26, 27, 35 
# (https://www.dallasnews.com/news/politics/2017/03/10/report-us-court-voids-texas-congressional-districts)
############################################################################################################
os.chdir("/Users/marybarker/Downloads/Texas/")
os.chdir("/home/thisisme/Documents/gerrymandering/Texas/")
CDsToLookAt = [15, 16, 20, 21, 23, 27, 28, 34, 35]

# read Congressional District shapefile first
ds = ogr.Open("TX_CURRENT_CD/CDS.shp")
lyr = ds.GetLayer(0)
CDS = [feat for feat in lyr] 
CDS = [CDS[x - 1] for x in CDsToLookAt]

# now get all vtds 
ds = ogr.Open("precinct/precinct.shp")
lyr = ds.GetLayer(0)
vtds = [feat for feat in lyr]

lookup = []
#counter = 0
#for vtd in vtds:
for counter in range(len(vtds)):
    vtd = vtds[counter]
    vtdg = vtd.geometry()
    mycd = -1
    myintarea = 0.0
    for cd in CDS:
        cdg = cd.geometry()
        if vtdg.Intersects(cdg):
            newintarea = vtdg.Intersection(cdg).Area()
            if newintarea > myintarea:
                mycd = cd.DISTRICT
                myintarea = newintarea
    if mycd != -1:
        lookup.append( (counter, vtd.CNTYVTD, vtd.VTDKEY, mycd) )

thing = pd.DataFrame(lookup, columns=['number', 'GEOID10', 'NAME10', 'CD'])
thing.to_csv("VTD_to_CD.csv")


# now build connectivity frame for vtds
adjacencyFrame = adjacencies([vtds[i] for i in myindices])
adjacencyFrame.to_csv("temporary_edges.csv")
vtdboundaries = boundaries([vtds[i] for i in myindices])
adjacencyFrame = adjacentEdgeLengths(adjacencyFrame, vtdboundaries)
adjacencyFrame.to_csv("PRECINCTconnections.csv")


#now write relevant vtd statistics to file
myindices = zip(*lookup)[0]
allthestats = pd.DataFrame()
for key in vtds[0].keys():
    allthestats[key] = [vtds[i][key] for i in myindices] 
allthestats['PERIM'] = [sum(adjacencyFrame.ix[(adjacencyFrame.low.astype(str) == (str(vtds[i].CNTYVTD)+str(vtds[i].VTDKEY))) | (adjacencyFrame.high.astype(str) == (str(vtds[i].CNTYVTD)+str(vtds[i].VTDKEY))), 'length']) for i in myindices]
allthestats['ALAND'] = [vtds[i].geometry().Area() for i in myindices]
allthestats['AWATER'] = 0
allthestats['CNTYVTD']  = allthestats['GEOID10']
allthestats[ 'VTDKEY'] = allthestats['NAME10']
allthestats.rename(columns={'e_total':'population'}, inplace=True)
allthestats['VTD'] = pd.Series([str(x) for x in allthestats.GEOID10.values]).str.cat([str(y) for y in allthestats.NAME10.values], sep='')
allthestats.to_csv("alt_vtdstats_1.csv")


outShapefile = str(os.getcwd()) + '/VTDS_of_Interest.shp'
outDriver = ogr.GetDriverByName("ESRI Shapefile")
if os.path.exists(outShapefile):
    outDriver.DeleteDataSource(outShapefile)
# Create the output shapefile
outDataSource = outDriver.CreateDataSource(outShapefile)
outLayer = outDataSource.CreateLayer('vtds', geom_type=ogr.wkbPolygon)

for key in vtds[0].keys():
    thisfield = ogr.FieldDefn(key, ogr.OFTInteger)
    outLayer.CreateField(thisfield)
outLayer.CreateField(ogr.FieldDefn('GEOID10', ogr.OFTInteger))
outLayer.CreateField(ogr.FieldDefn('NAME10', ogr.OFTInteger))

for vtd in [vtds[i] for i in myindices]:
    featureDefn = outLayer.GetLayerDefn()
    g = vtd.geometry()
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(g)

    for key in vtd.keys():
        feature.SetField(key, vtd[key])
    feature.SetField('GEOID10', vtd['CNTYVTD'])
    feature.SetField('NAME10', vtd['VTDKEY'])
    outLayer.CreateFeature(feature)
    feature = None

# Close DataSource
outDataSource.Destroy()


############################################################################################################
ds = ogr.Open("07_11_block_group_Texas.gdb")
lyr1 = ds.GetLayer(0)
lyr2 = ds.GetLayer(1)
stuff1 = [feat for feat in lyr1]
stuff2 = [feat for feat in lyr2]

# first find out all of the keys that we want and write their explanations to a csv
lookupFrame = [(a['Short_Name'], a['Full_Name']) for a in stuff1 if a['Short_Name'].endswith('e1')]
lookupFrame = zip(*lookupFrame)
pd.DataFrame({'ShortName':lookupFrame[0], 'LongName':lookupFrame[1]}).to_csv("meaningOfKeys.csv")


allTheColumnsIWant1 = list(stuff2[0].keys()[:16])
allTheColumnsIWant2 = list(lookupFrame[0])
allTheColumnsIWant = allTheColumnsIWant1+allTheColumnsIWant2

# Now we want to extract the important data for each block and write to csv 
bigDataCsv = pd.DataFrame({key:[stuffthing[key] for stuffthing in stuff2] for key in allTheColumnsIWant})
bigDataCsv.to_csv("blockGroupData.csv")

# want to find the proportion of block to block group. 
ds = ogr.Open("/home/thisisme/Downloads/tabblock2010_48_pophu.shp")
lyr = ds.GetLayer(0)
newblocks = [feat for feat in lyr]

# just the blocks that are in our subset of Texas for efficiency's sake. 
blockToCD = pd.read_csv('BlockAssign_ST48_TX_CD.txt')
blockToCD.BLOCKID = blockToCD.BLOCKID.astype(str)
blockToCD.set_index("BLOCKID", inplace=True)
blockToVTD = pd.read_csv("BlockAssign_ST48_TX_VTD.txt")
blockToVTD.BLOCKID = blockToVTD.BLOCKID.astype(str)
blockToVTD.set_index(blockToVTD.BLOCKID, inplace=True)
# write out population for each block
newblocks = [b for b in newblocks if int(blockToCD.DISTRICT[str(b.BLOCKID10)]) in CDsToLookAt]
pd.DataFrame({'id':[b.BLOCKID10 for b in newblocks], 'pop':[b.POP10 for b in newblocks]}).to_csv("blockPops.csv")
# cut short our blockToCD and blockToVTD values just to the subset of Texas we care about
blockToCD = blockToCD.ix[blockToCD.DISTRICT.isin(CDsToLookAt)]
blockToVTD = blockToVTD.ix[blockToCD.index]

blockstats=pd.read_csv("vtdstats.csv")
blockstats.GEOID10 = blockstats.GEOID10.astype(str)
blockstats.set_index(blockstats.GEOID10, inplace=True)
vtds = blockstats.CNTY.astype(str).str.zfill(3).str.cat(blockstats.CNTYVTD.astype(str).str.slice(start=3))
blockstats['compVTDS'] = vtds

indices = blockToVTD.ix[ blockToVTD.COUNTYFP.astype(str).str.zfill(3).str.cat(blockToVTD.DISTRICT.astype(str).str.zfill(4)).isin(blockstats.compVTDS.astype(str))]
blockToVTD = indices

blockpops = pd.read_csv("blockPops.csv")
blockpops['id'] = blockpops['id'].astype(str)
blockpops.set_index(blockpops['id'], inplace=True)

blockVTD = [ (blockstats.GEOID10[blockstats.compVTDS == (str(blockToVTD.COUNTYFP[block]).zfill(3)+str(blockToVTD.DISTRICT[block]).zfill(4) ) ], block)  for block in blockToVTD.index]
pd.DataFrame(blockVTD, columns=['vtd', 'block']).to_csv("block_to_VTD_in_id_form.csv")
blockVTD['population'] = [blockpops.population[x] for x in blockVTD.blockid]
blockVTD.to_csv("block_to_VTD_in_id_form.csv")


len(set(blockVTD.geoid.values))
