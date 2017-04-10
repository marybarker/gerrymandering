from osgeo import ogr
import numpy as np
import pandas as pd
import os

ds = ogr.Open("PA_CD_2011/PA_CD.shp")
lyr = ds.GetLayer(0)
mydata = features(lyr)

vtdfile = 'precinct/precinct.shp'
ds = ogr.Open(vtdfile)
lyr = ds.GetLayer(0)
vtds = features(lyr)


vtdNAME = []
vtdDIST = []
for vtd in vtds:
    vtdNAME.append(vtd.GEOID10+vtd.NAME10)
    vtdDist = ''
    current_area = 0.0
    for dist in mydata:
        if vtd.geometry().Within(dist.geometry()):
            vtdDist = dist.OBJECTID

    if vtdDist == '':
        for dist in mydata:
            if vtd.geometry().Intersects(dist.geometry()):
                new_area = dist.geometry().Intersection(vtd.geometry()).Area()
                if new_area >= current_area:
                    vtdDist = dist.OBJECTID
                    current_area = new_area
        
    vtdDIST.append(vtdDist)

Pennsylvanias_current = pd.DataFrame({'key':vtdNAME, 'value':vtdDIST})
Pennsylvanias_current.to_csv("current_PA_CDs.csv")

startfolder = 'starting_states/'
starting_states = [f for f in os.listdir(startfolder) if f.endswith('.csv')]

metrics = {}
updateGlobals(Pennsylvanias_current)
actualGoodness = goodness(metrics)

startGoodness = []
for file in starting_states:
    state = pd.read_csv(startfolder+file)
    updateGlobals(state)
    startGoodness.append(goodness(metrics))

print actualGoodness 
print startGoodness

