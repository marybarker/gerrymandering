# create new shapefile with merged vtds

# first merge all doughnut-ed vtds into single vtd
singles = [vtd for vtd in blockstats.VTD if adjacencyFrame.loc[(adjacencyFrame.low == vtd) | (adjacencyFrame.high == vtd)].shape[0] == 1]
toglom = []
for single in singles: 
    if any(adjacencyFrame.high == single):
        toglom.append(adjacencyFrame.low[adjacencyFrame.high == single].item())
    else:
        toglom.append(adjacencyFrame.high[adjacencyFrame.low == single].item())

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

# now read the shapefile and merge together desired ones. 
from osgeo import ogr

allthevtds = ogr.Open('precinct/precinct.shp')
lyr = allthevtds.GetLayer(0)
vtds = [feat for feat in lyr]

inShapefile = "/Users/marybarker/Documents/tarleton_misc/gerrymandering/Pennsylvania/precinct/precinct.shp"
outShapefile = "/Users/marybarker/Documents/tarleton_misc/gerrymandering/Pennsylvania/noIslandsprecinct1.shp"
outDriver = ogr.GetDriverByName("ESRI Shapefile")
if os.path.exists(outShapefile):
    outDriver.DeleteDataSource(outShapefile)
# Create the output shapefile
#outDataSource = outDriver.CreateDataSource(outShapefile)
# get all keys from previous shapefile
outDataSource = outDriver.CopyDataSource(ogr.Open(inShapefile), outShapefile)
outLayer = outDataSource.CreateLayer("vtds", geom_type=ogr.wkbPolygon)

for vtd in vtds:
    featureDefn = outLayer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)

    if (vtd.GEOID10+vtd.NAME10) in toglom:
        g = vtd.geometry()
        other = vtd.geometry()
        for others in vtds:
            if (others.GEOID10+others.NAME10) in singles:
                if toglom[singles.index(others.GEOID10+others.NAME10)] == (vtd.GEOID10+vtd.NAME10):
                    other = others.geometry()
                    g = g.Union(other)
        feature.SetGeometry(g)
    elif (vtd.GEOID10+vtd.NAME10) in singles:
        pass
    else:
        feature.SetGeometry(vtd.geometry())
    outLayer.CreateFeature(feature)
    feature = None

# Close DataSource
outDataSource.Destroy()




#centroids
centroid = []
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

