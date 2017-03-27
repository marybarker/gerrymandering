from osgeo import ogr
import numpy as np
import os
import pandas as pd
import math
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import colorsys

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
        if stateSHORT == 'PA':
            name = str(f1['GEOID10']) + f1['NAME10']
        elif stateSHORT == 'NH':
            name = str(f1['GEOID10'])
        g1 = f1.geometry()
        for f2 in mylistoffeatures[count+1:]:
            g2 = f2.geometry()
            if g1.Touches(g2):
                l1.append(name)
                if stateSHORT == 'PA':
                    l2.append(str(f2['GEOID10']) + f2['NAME10'])
                elif stateSHORT == 'NH':
                    l2.append(str(f2['GEOID10']))
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
            if stateSHORT == 'PA':
                boundaries[str(feat['GEOID10']) + feat['NAME10']] = zip(x, y)
            elif stateSHORT == 'NH':
                boundaries[str(feat['GEOID10'])] = zip(x, y)
        elif gtype == 3: # polygon
            x = []
            y = []
            for ring in geom:
                for i in xrange(ring.GetPointCount()):
                    point = ring.GetPoint(i)
                    x.append(point[0])
                    y.append(point[1])
            if stateSHORT == 'PA':
                boundaries[str(feat['GEOID10']) + feat['NAME10']] = zip(x, y)
            elif stateSHORT == 'NH':
                boundaries[str(feat['GEOID10'])] = zip(x, y)
        else:
            b = geom.GetBoundary()
            if stateSHORT == 'PA':
                boundaries[str(feat['GEOID10']) + feat['NAME10']] = b.GetPoints()
            elif stateSHORT == 'NH':
                boundaries[str(feat['GEOID10'])] = b.GetPoints()
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
        #names.append(str(vtd['GEOID10']) + vtd['NAME10'])
        names.append(str(vtd['GEOID10']))
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

def colorDict(n):
    return {i: colorsys.hsv_to_rgb(float(i)/n, 1, 1) for i in range(n)}

def color_these_states(geom_to_plot, list_of_states, foldername, number):
    colors = colorDict(ndistricts)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(geom_to_plot['xlim'])
    ax.set_ylim(geom_to_plot['ylim'])

    paths = geom_to_plot['paths']
    names = geom_to_plot['names']

    for i in range(len(list_of_states)):#state in list_of_states:
        this_state = list_of_states[i]
        redistricting = this_state[0]
        for p in range(len(paths)):
            path = paths[p]
            facecolor = redistricting.value[np.array(redistricting.key) == names[p][5:]].item()
            patch = mpatches.PathPatch(path,facecolor=colors[facecolor],edgecolor='black')
            ax.add_patch(patch)
        ax.set_aspect(1.0)
        plt.savefig(foldername+'/output%04d.png'%(number+i))
        #plt.show()

""" NOW USE ALL OF THE FUNCTIONS:"""

"""
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


names = [str(v['GEOID10'])+v['NAME10'] for v in vtds]
aland = [v['ALAND10'] for v in vtds]
awater = [v['AWATER10'] for v in vtds]
perim = [ToFeet(vtd_boundaries[name]) for name in names]
population = [v['POP100'] for v in vtds]
pres_repvotes = [v['USPRV2000'] for v in vtds]
pres_demvotes = [v['USPDV2000'] for v in vtds]

cong_repvotes = [v['USCRV2000'] for v in vtds]
cong_demvotes = [v['USCDV2000'] for v in vtds]

sen_repvotes = [v['USSRV2000'] for v in vtds]
sen_demvotes = [v['USSDV2000'] for v in vtds]

#Office Code Table
#------------------------------------------------------------
# 1 USP President of the United States
# 2 USS United States Senator
# 3 GOV Governor
# 4 LTG Lieutenant Governor
# 5 ATT Attorney General
# 6 AUD Auditor General
# 7 TRE State Treasurer
# 8 SPM Justice of the Supreme Court
# 9 SPR Judge of the Superior Court
#10 CCJ Judge of the Commonwealth Court
#11 USC Representative in Congress
#12 STS Senator in the General Assembly
#13 STH Representative in the General Assembly

VTDcsv = pd.DataFrame({
    'VTD':names, 
    'ALAND':aland, 
    'AWATER':awater, 
    'PERIM':perim, 
    'POP100':population, 
    'REP_C':cong_repvotes,
    'DEM_C':cong_demvotes,
    'REP_S':sen_repvotes,
    'DEM_S':sen_demvotes,
    'REP_P':pres_repvotes,
    'DEM_P':pres_demvotes
})
VTDcsv.to_csv('vtdstats.csv')
v.keys()
v['USPDV2000']
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

# figuring out which ones are neither polygons nor point files
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








