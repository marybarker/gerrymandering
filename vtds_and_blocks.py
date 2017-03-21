from osgeo import ogr
import numpy as np
import os
import pandas as pd
import math
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
"""   make a  list of all voting tabulation districts   """
""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
def features(layer):
    features = []
    for feat in layer: 
        features.append( feat )
    return features


""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
"""   total distance for a list of points in lat/long   """
""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
def ToFeet(listofpoints):
    R = 3961 * 5280
    x = np.asarray([point[0] for point in listofpoints])
    y = np.asarray([point[1] for point in listofpoints])

    dlon = (y[1:] - y[:-1])
    dlat = (x[1:] - x[:-1])
    a = (np.sin(math.pi * 0.5*dlat/180.0))**2 +  np.cos(math.pi*x[1:]/180.0) *  np.cos(math.pi*x[:-1]/180.0) * (np.sin(math.pi*0.5*dlon/180.0))**2
    c = 2 * np.arctan2( np.sqrt(a), np.sqrt(1-a) )
    d = R * sum( c )
    return d


""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
"""  create a lookup table of adjacencies between vtds  """
""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
def adjacencies(mylistoffeatures):
    l1 = list()
    l2 = list()
    for count in range(len(mylistoffeatures)):
        f1 = mylistoffeatures[count]
        name = f1['GEOID10']
        g1 = f1.geometry()
        for f2 in mylistoffeatures[count+1:]:
            g2 = f2.geometry()
            if g1.Touches(g2):
                l1.append(name)
                l2.append(f2['GEOID10'])
    newthing = pd.DataFrame(np.column_stack((np.array(l1), np.array(l2))))
    newthing.columns=['lo','hi']
    return newthing


""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
"""  create a lookup dict boundaries as lat/long lists  """
""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
def boundaries(mylistoffeatures):
    boundaries = {}
    for feat in mylistoffeatures:
        geom = feat.geometry()
        gtype = geom.GetGeometryType()

        if gtype == 6: 
            x = []
            y = []
            for i in xrange(geom.GetGeometryCount()):
                g = geom.GetGeometryRef(i)
                for ring in g:
                    for j in xrange(ring.GetPointCount()):
                        point = ring.GetPoint(j)
                        x.append(point[0])
                        y.append(point[1])
            boundaries[feat['GEOID10']] = zip(x, y)
        elif gtype == 3: # polygon
            x = []
            y = []
            for ring in geom:
                for i in xrange(ring.GetPointCount()):
                    point = ring.GetPoint(i)
                    x.append(point[0])
                    y.append(point[1])
            boundaries[feat['GEOID10']] = zip(x, y)
        else:
            b = geom.GetBoundary()
            boundaries[feat['GEOID10']] = b.GetPoints()
    return boundaries


""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
"""          get lengths of each connectivity           """
""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
def adjancentEdgeLengths(connectivitydf, boundaries):
    edgelengths = list()
    for i in range(np.shape(connectivitydf)[0]): 
        low = connectivitydf.lo[i]
        hi = connectivitydf.hi[i]
        b1 = boundaries[low]
        b2 = boundaries[hi]
        pointsInCommon = [point for point in b1 if point in b2]
        l = ToFeet(pointsInCommon)
        edgelengths.append(l)
    connectivitydf['length'] = edgelengths
    return connectivitydf

def package_vtds(filetouse):
    this_geom = {}

    # get extents of the geometry first of all 
    ds = ogr.Open(filetouse)
    nlay = ds.GetLayerCount()
    lyr = ds.GetLayer(0)
    ext = lyr.GetExtent()
    xoffset = (ext[1] - ext[0])/50
    yoffset = (ext[3] - ext[2])/50
    this_geom['xlim'] = [ext[0]-xoffset,ext[1]+xoffset]
    this_geom['ylim'] = [ext[2]-yoffset,ext[3]+yoffset]

    lyr.ResetReading()
    names = []
    paths = []

    for vtd in lyr:
        names.append(vtd['VTDST10_1'])
        geom = vtd.geometry()
        gtype = geom.GetGeometryType()
        
        codes = []
        x = []
        y = []

        # extract boundary points
        if gtype == 6: # who knows? 
            for i in range(geom.GetGeometryCount()):
                g = geom.GetGeometryRef(i)
                for ring in g:
                    for j in xrange(ring.GetPointCount()):
                        point = ring.GetPoint(j)
                        x.append(point[0])
                        y.append(point[1])
        elif gtype == 3: # polygon
            for ring in geom:
                for i in xrange(ring.GetPointCount()):
                    point = ring.GetPoint(i)
                    x.append(point[0])
                    y.append(point[1])
        else:
            all_x = []; all_y = []
            for i in range(geom.GetGeometryCount()):
                r = geom.GetGeometryRef(i)
                all_x = [r.GetX(j) for j in range(r.GetPointCount())]
                all_y = [r.GetY(j) for j in range(r.GetPointCount())]
                x = x + all_x
                y = y + all_y

        # create closed paths for each feature
        x = np.asarray(x)
        y = np.asarray(y)
        thing = [mpath.Path.MOVETO] + (len(x) - 1)*[mpath.Path.LINETO]
        codes = codes + thing
        A = np.column_stack((x, y))
        path = mpath.Path(vertices=A, codes=np.asarray(codes[:len(x)]))
        paths.append(path)

    this_geom['paths'] = paths
    this_geom['names'] = names
    return this_geom


""" NOW USE ALL OF THE FUNCTIONS:"""

blockfile = 'block/block.shp'
blocktotabdistrict = pd.read_csv('pa_block_to_vtd.txt')
state='NewHampshire/'
state='Pennsylvania/'
blockfile = 'NewHampshireBlockShapefiles/tl_2016_33_tabblock10.shp'
vtdfile = ("HarvardData/nh_final.shp")

# current working directory
os.chdir("/Users/marybarker/Dropbox/Tarleton/forfun/gerrymandering/"+state)

# census block file
ds = ogr.Open(blockfile)
lyr = ds.GetLayer(0)
censusblocks = features(lyr)
#block_connectivities = create_adjacencies(censusblocks)

# voting tabulation district 
vtdfile = 'precinct/precinct.shp'
ds = ogr.Open(vtdfile)
lyr = ds.GetLayer(0)
vtds = features(lyr)
vtd_boundaries = boundaries(vtds)
vtd_connectivities = adjacencies(vtds)
vtd_connectivities = adjancentEdgeLengths(vtd_connectivities, vtd_boundaries)
vtd_connectivities.to_csv('PRECINCTconnections.csv')



g = package_vtds(vtdfile)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(g['xlim'])
ax.set_ylim(g['ylim'])
paths = g['paths']

for p in range(len(paths)):
    path = paths[p]
    patch = mpatches.PathPatch(path,facecolor="green",edgecolor='black')
    ax.add_patch(patch)
ax.set_aspect(1.0)
plt.show()

lyr.ResetReading()
vtds = features(lyr)
"""

# figuring out which ones are neither polygons nor point files
"""
for b, val in vtd_boundaries.iteritems():
    if isinstance(val, list):
        pass
    else: 
        print val, b

faulty_names = ['33015ZZZZZZ', '33017SOME03']
faulty_feats = [f for f in vtds if f['GEOID10'] in faulty_names]
f1 = faulty_feats[0]
f2 = faulty_feats[1]
g1 = f1.geometry()
g2 = f2.geometry()
"""


