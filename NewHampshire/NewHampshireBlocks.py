from osgeo import ogr
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import os
import pandas as pd
import math

def ToFeet(listofpoints):
    R = 3961 * 5280
    x = np.asarray([point[0] for point in listofpoints])
    y = np.asarray([point[1] for point in listofpoints])

    dlon = (y[1:] - y[:-1])
    dlat = (x[1:] - x[:-1])
    a = (np.sin(math.pi * 0.5*dlat/180.0))**2 + np.cos(math.pi*x[1:]/180.0) * np.cos(math.pi*x[:-1]/180.0) * (np.sin(math.pi*0.5*dlon/180.0))**2
    c = 2 * np.arctan2( np.sqrt(a), np.sqrt(1-a) )
    d = R * sum( c )
    return d

os.chdir("/Users/marybarker/Dropbox/Tarleton/forfun/gerrymandering/NewHampshire/NewHampshireBlockShapefiles/")
ds = ogr.Open("tl_2016_33_tabblock10.shp")
nlay = ds.GetLayerCount()
lyr = ds.GetLayer(0)



""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
ext = lyr.GetExtent()
xoff = (ext[1] - ext[0])/50
yoff = (ext[3] - ext[2])/50

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(ext[0]-xoff,ext[1]+xoff)
ax.set_ylim(ext[2]-yoff,ext[3]+yoff)

paths = []
lyr.ResetReading()

for feat in lyr: 
    geom = feat.geometry()
    codes = []
    all_x = []
    all_y = []
    for i in range(geom.GetGeometryCount()):
        r = geom.GetGeometryRef(i)
        x = [r.GetX(j) for j in range(r.GetPointCount())]
        y = [r.GetY(j) for j in range(r.GetPointCount())]
        thing = [mpath.Path.MOVETO] + \
                (len(x) - 1)*[mpath.Path.LINETO]
        codes = codes + thing
        all_x = all_x + x
        all_y = all_y + y

    all_x = np.asarray(all_x)
    all_y = np.asarray(all_y)

    A = np.column_stack((all_x, all_y))
    path = mpath.Path(vertices=A, codes=np.asarray(codes[:len(all_x)]))
    paths.append(path)

for path in paths: 
    patch = mpatches.PathPatch(path, \
                facecolor="green", edgecolor="yellow")
    ax.add_patch(patch)

ax.set_aspect(1.0)
plt.show()
""" * * * * * * * * * * * * * * * * * * * * * * * * * * """


""" * * * * * * * * * * * * * * * * * * * * * * * * * * """
geoid = []
area  = []
warea = []
perim = []
lyr.ResetReading()
for feat in lyr: 

    geoid.append(feat['GEOID10'])
    area.append(feat['ALAND10'])
    warea.append(feat['AWATER10'])
    geom = feat.geometry()
    l = 0
    for i in range(geom.GetGeometryCount()):
        r = geom.GetGeometryRef(i)
        x = np.asarray([r.GetX(j) for j in range(r.GetPointCount())])
        y = np.asarray([r.GetY(j) for j in range(r.GetPointCount())])
        d = ToFeet(zip(x,y))

        """
        dlon = (y[1:] - y[:-1])
        dlat = (x[1:] - x[:-1])
        a = (np.sin(0.5*dlat))**2 + np.cos(x[1:]) * np.cos(x[:-1]) * (np.sin(0.5*dlon))**2
        c = 2 * np.arctan2( np.sqrt(a), np.sqrt(1-a) )
        d = R * sum( c )
        """
        l = l + d
    perim.append(l)

geoid = np.array(geoid)
area = np.array(area)
warea = np.array(warea)
perim = np.array(perim)
mydata = pd.DataFrame(np.column_stack((geoid, area, warea)))

NewHampshire = pd.DataFrame(np.column_stack((geoid, area, warea, perim)))
NewHampshire.columns = ['GEOID', 'AREA', 'WaterAREA', 'PERIMETER']
filename = 'NewHampshire.csv'
NewHampshire.to_csv(filename)
""" * * * * * * * * * * * * * * * * * * * * * * * * * * """


lyr.ResetReading()
# will be in lat/long, but that's ok because it's compared with other boundaries.
boundaries = {}
for feat in lyr:
    geom = feat.geometry()
    b = geom.GetBoundary()
    boundaries[feat['BLOCKCE10']] = b.GetPoints()


lyr.ResetReading()
features = []
for feat in lyr: 
    features.append(feat)

"""
myconnectivitydict = {}
for f1 in features[:10]: 
    name = f1['GEOID10']
    g1 = f1.geometry()
    listofconns = list()
    for f2 in features: 
        if f2['GEOID10'] != name: 
            g2 = f2.geometry()
            if g1.Touches(g2): 
                listofconns.append(f2['GEOID10'])
                print name, f2['GEOID10']

    myconnectivitydict[name] = listofconns
"""

things_to_save = ['GEOID10', 'ALAND10','AWATER10',
                  'POP100','POPMALE','POPFEMALE',
                  'POP1RACE','POPWHITE','POPBLACK',
                  'POPAMERIND','POPASIAN','POPPACIFIC',
                  'POPOTHER','POPMULTRAC','POPHISP']

mydata = pd.DataFrame(columns=things_to_save+['PERIM'])
for feat in features: 
    print feat.keys()
    things = [feat[t] for t in things_to_save]
    geom = feat.geometry()
    l = 0
    for i in range(geom.GetGeometryCount()):
        r = geom.GetGeometryRef(i)
        x = np.asarray([r.GetX(j) for j in range(r.GetPointCount())])
        y = np.asarray([r.GetY(j) for j in range(r.GetPointCount())])
        d = ToFeet(zip(x,y))
        l = l + d
    things.append(l)
    things = pd.DataFrame(np.array([things]), columns=things_to_save+['PERIM'])
    mydata = pd.concat([mydata,things])
mydata.to_csv('NHBlocks.csv')


# now get median income for each blockgroup and accumulate to blocks. 
os.chdir('/Users/marybarker/Dropbox/Tarleton/forfun/gerrymandering/NewHampshire/NewHampshireBlockShapefiles/')
blockgroups = ogr.Open('../BlockGroupInfo/2011_ACS_5YR_BG_33_NEW_HAMPSHIRE.gdb')
lyr1 = blockgroups.GetLayer(0)
lyr1.ResetReading()
features1 = []
for feat in lyr1: 
    features1.append(feat)

# for all codes, see bottom of file. 
income_by_blockgroup = np.matrix([[str(feat['GEOID']), feat['B19013e1'], feat['B19013m1']] for feat in features1])
income_by_blockgroup = pd.DataFrame(income_by_blockgroup)
income_by_blockgroup.columns=['GEOID','expected','median']
blocks_and_blockgroups = pd.read_csv('../NewHampshireBlockAssignment/BlockAssign_ST33_NH_VTD.txt')
mydata = pd.read_csv('../HarvardData/NHBlocks.csv')

income = list()
for block in mydata['GEOID10']:
    shortened = str(block)[:12]
    value = income_by_blockgroup.expected[np.array(income_by_blockgroup.GEOID) == shortened]
    print value
    income.append(value)
income


"""
    bg = income_by_blockgroup['B19013e1']
    ids = list(income_by_blockgroup['GEOID'])
    #try: 
    whereto = [f for f in ids if str(block) in str(f)]
    print whereto
    #    #amount = bg[np.array(income_by_blockgroup.GEOID) == whereto[0]]#block in list(income_by_blockgroup['GEOID'])]
    #    #income.append(amount)
    #except: 
    #    pass
    #    #print str(block)
"""


"""
Number|Title
B01001|Sex By Age
B01002|Median Age By Sex
B02001|Race
B02008|White Alone Or In Combination With One Or More Other Races
B02009|Black Or African American Alone Or In Combination With One Or More Other Races
B02010|American Indian And Alaska Native  Alone Or In Combination With One Or More Other Races
B02011|Asian Alone Or In Combination With One Or More Other Races
B02012|Native   Hawaiian   And   Other   Pacific   Islander   Alone   Or   In Combination With One Or More Other Races
B02013|Some Other Race Alone Or In Combination With One Or More Other Races
B03002|Hispanic Or Latino Origin By Race
B03003|Hispanic Or Latino Origin
B08301|Means Of Transportation To Work
B08303|Travel Time To Work
B09002|Own Children Under 18 Years By Family Type And Age
B09017|Relationship By Household Type (Including Living Alone) For The Population 65 Years And Over
B11001|Household Type (Including Living Alone)
B11016|Household Type By Household Size
B12001|Sex By Marital Status For The Population 15 Years And Over
B15002|Sex By Educational Attainment For The Population 25 Years And Over
B16004|Age By Language Spoken At Home By Ability To Speak English For The Population 5 Years And Over
B17017|Poverty Status In The Past 12 Months By Household Type By Age Of Householder
B17021|Poverty Status Of Individuals In The Past 12 Months By Living Arrangement
B19001|Household Income In The Past 12 Months (In 2011 Inflation- Adjusted Dollars)
B19013|Median  Household  Income  In  The  Past  12  Months  (In  2011 Inflation-Adjusted Dollars)
B19051|Earnings In The Past 12 Months For Households
B19052|Wage Or Salary Income In The Past 12 Months For Households
B19053|Self-Employment Income In The Past 12 Months For Households
B19054|"Interest, Dividends, Or Net Rental Income In The Past 12 Months For Households"
B19055|Social Security Income In The Past 12 Months For Households
B19056|Supplemental Security Income (Ssi) In The Past 12 Months For Households
B19057|Public Assistance Income In The Past 12 Months For Households
B19059|Retirement Income In The Past 12 Months For Households
B19060|Other Types Of Income In The Past 12 Months For Households
B19101|Family Income In The Past 12 Months (In 2011 Inflation-Adjusted Dollars)
B19113|Median Family Income In The Past 12 Months (In 2011 Inflation- Adjusted Dollars)
B19201|Nonfamily Household Income In The Past 12 Months (In 2011 Inflation-Adjusted Dollars)
B19202|Median Nonfamily Household Income In The Past 12 Months (In 2011 Inflation-Adjusted Dollars)
B19301|Per Capita Income In The Past 12 Months (In 2011 Inflation-Adjusted Dollars)
B21001|Sex By Age By Veteran Status For The Civilian Population 18 Years And Over
B21002|Period Of Military Service For Civilian Veterans 18 Years And Over
B25001|Housing Units
B25002|Occupancy Status
B25003|Tenure
B25004|Vacancy Status
B25006|Race Of Householder
B25010|Average Household Size Of Occupied Housing Units By Tenure
B25014|Tenure By Occupants Per Room
B25017|Rooms
B25018|Median Number Of Rooms
B25024|Units In Structure
B25034|Year Structure Built
B25035|Median Year Structure Built
B25040|House Heating Fuel
B25041|Bedrooms
B25044|Tenure By Vehicles Available
B25047|Plumbing Facilities For All Housing Units
B25051|Kitchen Facilities For All Housing Units
B25056|Contract Rent
B25061|Rent Asked
B25063|Gross Rent
B25064|Median Gross Rent (Dollars)
B25075|Value
B25077|Median Value (Dollars)
B25081|Mortgage Status
B25085|Price Asked
B25087|Mortgage Status And Selected Monthly Owner Costs
C17002|Ratio Of Income To Poverty Level In The Past 12 Months
C24010|Sex By Occupation For The Civilian Employed Population 16 Years And Over
C25095|Household  Income  By  Selected  Monthly  Owner  Costs  As  A Percentage Of Household Income In The Past 12 Months
"""



