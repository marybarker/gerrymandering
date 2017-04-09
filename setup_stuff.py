from osgeo import ogr
import numpy as np
import os
import pandas as pd
import math
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import colorsys
plt.rcParams['agg.path.chunksize'] = 1000
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
        name = str(f1['GEOID10']) + f1['NAME10']
        g1 = f1.geometry()
        for f2 in mylistoffeatures[count+1:]:
            g2 = f2.geometry()
            if g1.Touches(g2):
                l1.append(name)
                l2.append(str(f2['GEOID10']) + f2['NAME10'])
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
            boundaries[str(feat['GEOID10']) + feat['NAME10']] = zip(x, y)
        elif gtype == 3: # polygon
            x = []
            y = []
            for ring in geom:
                for i in xrange(ring.GetPointCount()):
                    point = ring.GetPoint(i)
                    x.append(point[0])
                    y.append(point[1])
            boundaries[str(feat['GEOID10']) + feat['NAME10']] = zip(x, y)
        else:
            b = geom.GetBoundary()
            boundaries[str(feat['GEOID10']) + feat['NAME10']] = b.GetPoints()
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
        name = str(vtd['GEOID10']) + vtd['NAME10']
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
        names.append(name)

    this_geom['paths'] = paths
    this_geom['names'] = names
    return this_geom

def colorDict(n):
    return {i: colorsys.hsv_to_rgb(float(i)/n, 1, 1) for i in range(n)}

def color_these_states(geom_to_plot, list_of_states, foldername, number):
    colors = colorDict(ndistricts)
    #colors = {0:'yellow',1:'green'}
    #ax.set_xlim([-71.8, -71.2])
    #ax.set_ylim([42.6, 43.2])

    paths = geom_to_plot['paths']
    names = geom_to_plot['names']

    for i in range(len(list_of_states)):#state in list_of_states:
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim(geom_to_plot['xlim'])
        ax.set_ylim(geom_to_plot['ylim'])
        
        this_state = list_of_states[i]
        redistricting = this_state[0]
        redistricting = redistricting.drop('Unnamed: 0', 1)
        redistricting.columns = ['key', 'value']
        for p in range(len(paths)):
            path = paths[p]
            
            facecolor = redistricting.value[np.array(redistricting.key) == names[p]].item()
            patch = mpatches.PathPatch(path,facecolor=colors[facecolor],edgecolor='black')
            ax.add_patch(patch)
        ax.set_aspect(1.0)
        plt.show()
        #plt.savefig(foldername+'output%04d.png'%(number+i))
        #plt.clf()
        #del fig
    

precinctBoundaryFile =  'precinct/precinct.shp'
precinctStatsFile = 'vtdstats.csv'
precinctConnectionsFile = 'PRECINCTconnections.csv'

ds = ogr.Open(precinctBoundaryFile)
lyr = ds.GetLayer(0)
precincts = features(lyr)
precinctBoundaries = boundaries(precincts)
precinctStats = pd.read_csv(precinctStatsFile)
precinctConns = pd.read_csv(precinctConnectionsFile)

# need: 
#  - contiguous start w r t shapefiles (globbed out as before, or using an overlaid shape? )
#  - way of associating precinct stats wtih the features of the shapefile





""" NOW USE ALL OF THE FUNCTIONS:"""

"""
# voting tabulation district 
vtdfile = 'precinct/precinct.shp'
ds = ogr.Open(vtdfile)
lyr = ds.GetLayer(0)
vtds = features(lyr)

vtd_boundaries = boundaries(vtds)
vtd_connectivities = adjacencies(vtds)
vtd_connectivities = adjancentEdgeLengths(vtd_connectivities, vtd_boundaries)
vtd_connectivities.to_csv('NEWVTDconnections.csv')






# plot neighbors of each VTD
names = g['names']
count = 0
for connection in names:
    count = count + 1
    lo_connected = vtd_connectivities.hi[vtd_connectivities.lo == connection]
    hi_connected = vtd_connectivities.lo[vtd_connectivities.hi == connection]
    colors = [0 for x in range(len(names))]
    l1 = [x for x in range(len(names)) if names[x] in set(lo_connected)]
    l2 = [x for x in range(len(names)) if names[x] in set(hi_connected)]
    for x in l1+l2:
        colors[x] = 1
    newframe = pd.DataFrame({'key':names, 'value':colors})
    color_these_states(g, [(newframe, 0)], foldername, count)

"""

