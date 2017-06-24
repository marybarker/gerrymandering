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
    features = [feat for feat in layer]
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
def adjacencies(mylistoffeatures, name_of_keys=['GEOID10', 'NAME10']):
    l1 = list()
    l2 = list()
    for count in range(len(mylistoffeatures)):
        f1 = mylistoffeatures[count]
        name =  ''.join([ str(f1[x]) for x in name_of_keys])
        g1 = f1.geometry()
        for f2 in mylistoffeatures[count+1:]:
            g2 = f2.geometry()
            if g1.Touches(g2):
                l1.append(name)
                l2.append( ''.join([ str(f2[x]) for x in name_of_keys]) )
    newthing = pd.DataFrame(np.column_stack((np.array(l1), np.array(l2))))
    newthing.columns=['low','high']
    return newthing


""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
"""  create a lookup dict boundaries as lat/long lists  """
""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
def boundaries(mylistoffeatures, name_of_keys = ['GEOID10', 'NAME10']):
    boundaries = {}
    for feat in mylistoffeatures:
        geom = feat.geometry()
        gtype = geom.GetGeometryType()

        if gtype == 6: 
            allxy = []
            for i in xrange(geom.GetGeometryCount()):
                g = geom.GetGeometryRef(i)
                for ring in g:
                    xy = []
                    for j in xrange(ring.GetPointCount()):
                        point = ring.GetPoint(j)
                        xy.append(point)
                    allxy.append(xy)
            boundaries[ ''.join([ str(feat[x]) for x in name_of_keys]) ] = allxy
        elif gtype == 3: # polygon
            allxy = []
            for ring in geom:
                xy = []
                for i in xrange(ring.GetPointCount()):
                    point = ring.GetPoint(i)
                    xy.append(point)
                allxy.append(xy)
            boundaries[ ''.join([ str(feat[x]) for x in name_of_keys]) ] = allxy
        else:
            b = geom.GetBoundary()
            boundaries[ ''.join([ str(feat[x]) for x in name_of_keys]) ] = [b.GetPoints()]
    return boundaries


""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
"""          get lengths of each connectivity           """
""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
def adjacentEdgeLengths(connectivitydf, boundaries):
    edgelengths = list()
    for i in range(np.shape(connectivitydf)[0]): 
        lo = connectivitydf.low[i]
        hi = connectivitydf.high[i]
        b1 = boundaries[lo]
        b2 = boundaries[hi]
        l = 0.0
        for b11 in b1: 
            b = len(b11)
            if b > 1:
                for b22 in b2:
                    pointsInCommon = [ [b11[a], b11[(a+1)%b] ] for a in range(b) if ((b11[a] in b22) and (b11[(a+1)%b] in b22)) ]
                    l += sum( [ToFeet(a) for a in pointsInCommon] )
        edgelengths.append(l)
    connectivitydf['length'] = edgelengths
    return connectivitydf

def package_vtds(shapefile_to_use, id_to_number_lookup_file, name_of_keys=['GEOID10', 'NAME10']):
    this_geom = {}

    # get extents of the geometry first of all 
    ds = ogr.Open(shapefile_to_use)
    nlay = ds.GetLayerCount()
    lyr = ds.GetLayer(0)
    ext = lyr.GetExtent()
    xoffset = (ext[1] - ext[0])/50
    yoffset = (ext[3] - ext[2])/50
    this_geom['xlim'] = [ext[0]-xoffset,ext[1]+xoffset]
    this_geom['ylim'] = [ext[2]-yoffset,ext[3]+yoffset]

    lookup = pd.read_csv(id_to_number_lookup_file)
    lookup = dict(zip(lookup.GEOID, lookup.IDNUM))

    lyr.ResetReading()
    names = []
    paths = []

    for vtd in lyr:
        name = lookup[ ''.join([str(vtd[keyname]) for keyname in name_of_keys]) ]
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
    districtColors = {i: colorsys.hsv_to_rgb(float(i)/n, 1, 1) for i in range(n)}
    districtColors[n] = colorsys.hsv_to_rgb(0, 0, 0.5)
    return districtColors

def color_these_states(geom_to_plot, list_of_states, foldername, number, linewidth = 1, DPI = 300):
    colors = colorDict(ndistricts)

    paths = geom_to_plot['paths']
    names = geom_to_plot['names']

    for i in range(len(list_of_states)):#state in list_of_states:
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim(geom_to_plot['xlim'])
        ax.set_ylim(geom_to_plot['ylim'])
        
        this_state = list_of_states[i]
        redistricting = this_state[0]
        for p in range(len(paths)):
            path = paths[p]
            if names[p] in redistricting.key.values:
                facecolor = redistricting.value[np.array(redistricting.key) == names[p]].item()
                patch = mpatches.PathPatch(path,facecolor=colors[facecolor],edgecolor='black', linewidth = linewidth)
                ax.add_patch(patch)
        ax.set_aspect(1.0)
        #plt.show()
        plt.savefig(foldername+'output%04d.png'%(number+i), dpi=DPI)
        plt.clf()
        del fig

def color_this_state(geom_to_plot, state, filename, linewidth = 1, DPI = 300):
    #                In setup.py for each state, there should be a line like the following
    #                    g = package_vtds("./VTDS_of_Interest.shp")
    #                g contains the geometries and names of the VTDS for your state space.
    #                              |
    #                              state should be a pandas dataframe that has columns "key" and "value",
    #                                  with index set to be key.  ( state = state.set_index(state.key) )
    #                              'key' is a column of vtds to color, and 'value' are integers corresponding
    #                                  to district number.
    
    colors = colorDict(ndistricts)
    paths = geom_to_plot['paths']
    names = geom_to_plot['names']

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(geom_to_plot['xlim'])
    ax.set_ylim(geom_to_plot['ylim'])
    
    #redistricting = redistricting.drop('Unnamed: 0', 1)
    #redistricting.columns = ['key', 'value']
    for p in range(len(paths)):
        path = paths[p]
        if names[p] in state.key.values:
            facecolor = state.value[np.array(state.key) == names[p]].item()
            patch = mpatches.PathPatch(path,facecolor=colors[facecolor],edgecolor='black', linewidth = linewidth)#colors[facecolor])#'black')
            ax.add_patch(patch)
    ax.set_aspect(1.0)
    #plt.show()
    plt.savefig(filename, dpi=DPI)
    plt.clf()
    del fig


def color_by_rgb(geom_to_plot, vtds_rgb_dict, filename, linewidth = 1, DPI = 300):
    #            In setup.py for each state, there should be a line like the following
    #                g = package_vtds("./VTDS_of_Interest.shp")
    #            g contains the geometries and names of the VTDS for your state space.
    #                          |
    #                          vtds_rgb_dict should be a dictionary of VTD names with RGB triples to color them.
    
    subset = vtds_rgb_dict.keys()
    thing = zip(geom_to_plot['paths'], geom_to_plot['names'])
    paths = [x[0] for x in thing if x[1] in subset]
    names = [x[1] for x in thing if x[1] in subset]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(geom_to_plot['xlim'])
    ax.set_ylim(geom_to_plot['ylim'])
    
    for p in range(len(paths)):
        path = paths[p]
        facecolor = vtds_rgb_dict[names[p]]
        patch = mpatches.PathPatch(path,facecolor=facecolor, edgecolor='black', linewidth=linewidth)
        ax.add_patch(patch)
    ax.set_aspect(1.0)
    #plt.show()
    plt.savefig(filename, dpi=DPI)
    plt.clf()
    del fig


#precinctBoundaryFile =  'precinct/precinct.shp'
#precinctStatsFile = 'vtdstats.csv'
#precinctConnectionsFile = 'PRECINCTconnections.csv'

#ds = ogr.Open(precinctBoundaryFile)
#lyr = ds.GetLayer(0)
#precincts = features(lyr)
#precinctBoundaries = boundaries(precincts)
#precinctStats = pd.read_csv(precinctStatsFile)
#precinctConns = pd.read_csv(precinctConnectionsFile)

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

vtd_boundaries = boundaries(vtds, ['GEOID10'])
vtd_connectivities = adjacencies(vtds, ['GEOID10'])
vtd_connectivities = adjacentEdgeLengths(vtd_connectivities, vtd_boundaries)
vtd_connectivities.to_csv('PRECINCTconnections.csv')






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














