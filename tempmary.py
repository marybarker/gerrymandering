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
ds = ogr.Open("cds/cds.shp")
lyr = ds.GetLayer(0)
CDS = [feat for feat in lyr] 
CDS = [CDS[x - 1] for x in CDsToLookAt]

# now get all vtds 
ds = ogr.Open("precinct/precinct.shp")
ds = ogr.Open("baseline_VTDS.shp")
lyr = ds.GetLayer(0)
vtds = [feat for feat in lyr]

funname = 'current_'
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
        #lookup.append( (counter, vtd.GEOID10, vtd.NAME10, mycd) )
        lookup.append( (counter, str(vtd.CNTYVTD), str(vtd.VTDKEY), mycd) )

thing = pd.DataFrame(lookup, columns=['number', 'GEOID10', 'NAME10', 'CD'])
thing.to_csv(funname+"VTD_to_CD.csv")


# now build connectivity frame for vtds
adjacencyFrame = adjacencies([vtds[i] for i in myindices])
adjacencyFrame.to_csv(funname+"temporary_edges.csv")
vtdboundaries = boundaries([vtds[i] for i in myindices])
adjacencyFrame = adjacentEdgeLengths(adjacencyFrame, vtdboundaries)
adjacencyFrame.to_csv(funname+"PRECINCTconnections.csv")


#now write relevant vtd statistics to file
myindices = thing['number'].values#zip(*lookup)[0]
allthestats = pd.DataFrame()
for key in vtds[0].keys():
    allthestats[key] = [vtds[i][key] for i in myindices] 

allthestats['PERIM'] = [sum(adjacencyFrame.ix[(adjacencyFrame.low.astype(str) == (str(vtds[i].CNTYVTD)+str(vtds[i].VTDKEY))) | (adjacencyFrame.high.astype(str) == (str(vtds[i].CNTYVTD)+str(vtds[i].VTDKEY))), 'length']) for i in myindices]
allthestats['ALAND'] = [vtds[i].geometry().Area() for i in myindices]
allthestats['AWATER'] = 0
allthestats['GEOID10'] = allthestats['CNTYVTD'].astype(str)
allthestats['NAME10'] = allthestats[ 'VTDKEY'].astype(str)
#allthestats['CNTYVTD']  = allthestats['GEOID10']
#allthestats[ 'VTDKEY'] = allthestats['NAME10']
allthestats.rename(columns={'e_total':'population'}, inplace=True)
allthestats['VTD'] = allthestats.GEOID10.astype(str).str.cat(allthestats.NAME10.astype(str))
allthestats.to_csv(funname+"vtdstats.csv")


outShapefile = str(os.getcwd()) + '/'+ 'best_ever_' + 'block_groups.shp'
outDriver = ogr.GetDriverByName("ESRI Shapefile")
if os.path.exists(outShapefile):
    outDriver.DeleteDataSource(outShapefile)
# Create the output shapefile
outDataSource = outDriver.CreateDataSource(outShapefile)
outLayer = outDataSource.CreateLayer('vtds', geom_type=ogr.wkbPolygon)

for key in vtds[0].keys():
    thisfield = ogr.FieldDefn(key, ogr.OFTString)
    outLayer.CreateField(thisfield)
#outLayer.CreateField(ogr.FieldDefn('GEOID10', ogr.OFTString))
#outLayer.CreateField(ogr.FieldDefn('NAME10', ogr.OFTString))

for vtd in vtds: #[vtds[i] for i in myindices]:
    featureDefn = outLayer.GetLayerDefn()
    g = vtd.geometry()
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(g)

    for key in vtd.keys():
        feature.SetField(key, str(vtd[key]))
    #feature.SetField('GEOID10', str(vtd.CNTYVTD))
    #feature.SetField('NAME10', str(vtd.VTDKEY))
    outLayer.CreateFeature(feature)
    feature = None

# Close DataSource
outDataSource.Destroy()



os.chdir("/Users/marybarker/Dropbox/temp/")
ds = ogr.Open("current_bgs_of_Interest.shp")
lyr = ds.GetLayer()
bgs = [feat for feat in lyr]

thing1 = list(bgs[0].keys()[:15])
thing2 = list(pd.read_csv("meaningOfKeys.csv").ShortName.values)
thingstokeep = thing1 + thing2

outputfilename = 'blockGroupData.csv'
outputstuff = pd.DataFrame({reference:[str(x[str(reference)]) for x in bgs] for reference in thingstokeep})
outputstuff.to_csv(outputfilename)


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


############################################################################################################
# North Carolina

os.chdir("/Users/marybarker/Documents/tarleton_misc/gerrymandering/NorthCarolina/")
popdata = pd.read_excel("VTDTotalPopRaceAndEthnicity.xlsx", skiprows=1)
cntyToFips = pd.read_csv("county_to_fips.csv", header=None)
cntyToFips.columns=['state', 'stateFIPS', 'countyFIPS', 'countyName', 'fipsClass']
cntyToFips.countyName = [x.replace(' County', '') for x in cntyToFips.countyName]
countyFIPS = cntyToFips.countyFIPS.copy()
cntyToFips.set_index(countyFIPS, inplace=True)
blockstats = pd.read_csv("vtdstats.csv")


lookupvtd = blockstats.VTDST10
lookupcnty = [cntyToFips.countyName[x] for x in blockstats.COUNTYFP10]

blockstats.VTDST10  == popdata.VTD
cntyToFips[ blockstats.COUNTYFP10, 'countyName'] == popdata.County
lookup_table = pd.DataFrame({'GEOID10':blockstats.GEOID10, 'excelCountyName':lookupcnty, 'excelVTDId':lookupvtd })
lookup_table.to_csv("blockstats_to_excel_lookup.csv")
lookup_table.set_index('GEOID10', inplace=True)
population = [popdata.ix[((popdata.County == lookup_table.excelCountyName[x]) & (popdata.VTD == lookup_table.excelVTDId[x])), 'Total'] for x in blockstats.GEOID10]
population = [sum(x.values) for x in population]

blockstats['population'] = population

allTheColumnsIWant1 = list(stuff2[0].keys()[:16])
allTheColumnsIWant2 = list(lookupFrame[0])
allTheColumnsIWant = allTheColumnsIWant1+allTheColumnsIWant2

# Now we want to extract the important data for each block and write to csv 
bigDataCsv = pd.DataFrame({key:[stuffthing[key] for stuffthing in stuff2] for key in allTheColumnsIWant})
bigDataCsv.to_csv("blockGroupData.csv")

