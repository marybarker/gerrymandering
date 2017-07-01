# Things I now know: 

# Percent black  21.5
# Percent white  68.5
# Percent Hisp    8.4

# Percent registered as democrat: 46
# Percent registered as republican: 31 
# Percent registered unafilliated: 23.5
# Percent registered libertarian: 0.1

# average registered voter age: 47.2 w standard deviation ~3


"""
allthingy = pd.merge(pd.read_csv("TotalPopRaceAndEthnicity.csv"), pd.read_csv("RegPartyRace.csv"), on='GEOID10')
allthingy = allthingy.merge(pd.read_csv("RegGenderAgeEthnicity.csv"), on='GEOID10')

plt.scatter(allthingy['VR: All Reps'] / allthingy['VR Total'], allthingy['% White Non Hisp'])
plt.xlabel("percent republican")
plt.ylabel("percent white non-hispanic")
plt.title('White Non-Hispanic Republican')

plt.scatter(allthingy['VR: All Reps'] / allthingy['VR Total'], 1.0 - allthingy.loc[:, ['% Total Black', '%  Hisp ']].sum(axis=1))
plt.xlabel("percent republican")
plt.ylabel("percent non-minority")
plt.title("Non-minority republican")

plt.scatter(allthingy['VR: All Dems'] / allthingy['VR Total'], allthingy['% Total Black'])
plt.xlabel("percent democratic")
plt.ylabel("percent black")
plt.title("black democrat")

plt.scatter(allthingy['Voter Registration by Gender Male %'], allthingy['VR: All Dems'] / allthingy['VR Total'])

plt.scatter(allthingy['Voter Registration by Gender Female %'], allthingy['% Total Black'])

x = len(allthingy['VR: All Reps'])
plt.scatter(allthingy['% Total Black'], 
            allthingy['VR: All Dems']/allthingy['VR Total'], 
            c=np.column_stack((np.ones(x), 
                               np.zeros(x), 
                               np.zeros(x), 
                               allthingy['VR Total'].astype(float) / np.maximum( allthingy['Total'].astype(float), 
                               allthingy['VR Total'].astype(float)))) )


meanAge = allthingy.loc[:, ['Voter Registration by Age 18-25 %', \
                            'Voter Registration by Age 26-40 %', \
                            'Voter Registration by Age 41-65 %', \
                            'Voter Registration by Age 66+ %']].apply(lambda x: sum(x * np.array([21.5, 33, 53, 66])), axis=1)

# allthingy['Voter Registration by Age 18-25 %'] stats: 
#  average is 0.10648 with standard dev of 0.0602

np.average(allthingy.loc[allthingy['Voter Registration by Age 18-25 %'] >= 0.22, '% White '] - np.average(allthingy['% White ']))
#-0.10491356745559763 is the value, meaning that when we are 2 standard deviations away from the average for 
# young voter registration, then we are on average 10 percent lower white (i.e. young voters != white)


allthingy.loc[allthingy['Voter Registration by Age 26-40 %'] >= 0.3, 'VR Total']


plt.scatter(allthingy['% Total Black'], 
            allthingy['VR: All Dems']/allthingy['VR Total'], 
            #c=np.column_stack((np.zeros(x), allthingy[''], 1.0 - allthingy[''], np.ones(x)*0.5 )))
            c=np.column_stack((np.zeros(x), 
                               allthingy['VR Total'].astype(float) / np.maximum( allthingy['Total'].astype(float), 
                               allthingy['VR Total'].astype(float)), 
                               1.0 - allthingy['VR Total'].astype(float) / np.maximum( allthingy['Total'].astype(float), 
                               allthingy['VR Total'].astype(float)), np.ones(x)*0.25 )))

plt.scatter(allthingy['VR: All Reps']/allthingy['VR Total'], 
            allthingy['VR: All Dems']/allthingy['VR Total'], 
            c=np.column_stack((np.zeros(x), meanAge / max(meanAge), 1. - meanAge / max(meanAge), np.ones(x)*0.125 )))
"""

goodnessWeights=[1, 100, 100, 100, 0, 0]
_exploration = 1
currentNCstate = pd.read_csv("VTD_to_CD.csv").loc[:, ['GEOID10', 'CD']].rename(columns={"GEOID10":"key", 'CD':'value'})
currentNCstate.key   = [lookup.get(x) for x in currentNCstate.key]
currentNCstate.value = [int(x) - 3701 for x in currentNCstate.value]
compare_current_state_to_possible_perturbations(currentNCstate, 100, 1000000, foldername, "initial look at NC demographics")

states = [(pd.read_csv(foldername+'data/state%d.csv'%i), 0) for i in range(10)]


def compare_current_state_to_possible_perturbations(current_state, num_perturbations, steps_per_perturbation, foldername, misc_data =''):
    import datetime, inspect
    global goodnessWeights, goodnessParams
    if not os.path.isdir(foldername):
        os.mkdir(foldername)
    os.chdir(foldername)
    
    metadata = 30*'+'+\
                '\nDate: %s'%str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                '\nSteps per run: %d'%(steps_per_perturbation)+\
                '\ngoodness function: %s'%inspect.getsource(goodness)+\
                '\ngoodness Weights: %s'%' '.join([str(x) for x in goodnessWeights])+\
                '\ngoodness Params: %s'%' '.join([str(x) for x in goodnessParams])+\
                '\nmisc: %s'%misc_data
    with open("Info.txt", "w") as metafile:
        metafile.write(metadata)
    
    if not os.path.isdir('data/'):
        os.mkdir('data/')
    if not os.path.isdir('pictures/'):
        os.mkdir('pictures/')
    updateGlobals(current_state)
    
    allStates = []
    allMetrics = []
    
    current_state.to_csv("data/initialState.csv", index=False)
    metrics.to_csv("data/initialMetrics.csv", index=False)
    initmet = metrics.copy()
    
    for i in range(num_perturbations):
        updateGlobals(current_state)
        runningState = MH(current_state, steps_per_perturbation, neighbor, goodness, switchDistrict)
        runningState[0].to_csv("data/state%d.csv"%i, index=False)
        metrics.to_csv("data/metrics%d.csv"%i, index=False)
        allStates.append(runningState[0])
        allMetrics.append(metrics)
        print 'finished with step %d of %d. Saved in %sdata/'%(i+1, num_perturbations, foldername)
    
    tempgoodnessParams  = [x for x in goodnessParams]
    tempgoodnessWeights = [x for x in goodnessWeights]

    goodnessParams  = [contScore, popVarScore]
    goodnessWeights = np.array([1, 500])

    allGoodnesses = [goodness(allMetrics[i]) for i in range(num_perturbations)]
    plt.plot( allGoodnesses )
    plt.savefig("pictures/goodness.png")
    plt.clf()
    for met in allMetrics:
        plt.scatter(range(ndistricts), met.mincon)
    plt.scatter(range(ndistricts), initmet.mincon, c='black',s=10)
    plt.scatter(range(ndistricts), initmet.mincon, c='black',s=250, alpha=0.125)
    plt.savefig("pictures/minorityConc.png", dpi=600)
    plt.clf()
    
    goodnessParams  = tempgoodnessParams
    goodnessWeights = tempgoodnessWeights
    os.chdir("../")

#def compare_current_state_to_possible_other_states()


if False:
    realStateOfTexas = pd.read_csv("current_state_of_texas.csv")
    del realStateOfTexas['GEOID']
    rep = {'CD':'CD', 'ID':'key'}
    realStateOfTexas.columns = [rep.get(x) for x in realStateOfTexas.columns]
    
    thelist = list(set(realStateOfTexas.CD))
    realStateOfTexas['value'] = 0
    for x in thelist:
        realStateOfTexas.loc[realStateOfTexas.CD == x, 'value'] = thelist.index(x)
    plt.errorbar(range(-1, num_perturbations+2), [allGoodnesses[0]]+allGoodnesses+[allGoodnesses[-1]], yerr = np.zeros(num_perturbations+3), label='unconverged goodness')
    plt.errorbar(range(-1, num_perturbations+2), [each_mean[0]]+each_mean+[each_mean[-1]], yerr = [0]+each_std+[0], label='average converged goodness')
    plt.legend()
    
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
    
    GEE = pd.read_excel("VTDRegGenderAgeEthnicity.xlsx", skiprows=1, skipfooter=1)
    PRE = pd.read_excel("VTDTotalPopRaceAndEthnicity.xlsx")
    VAP = pd.read_excel("VTDVotingAgePopulation.xlsx", skipfooter=1, header=[0,1])
    RPR = pd.read_excel("VTDRegPartyRace.xlsx", skiprows=1, header=[0, 1], skipfooter=1)
    
    RPR.columns = [' '.join([ x for x in col if 'Unnamed' not in x]) for col in RPR.columns.values]
    GEE.columns = [' '.join([a for a in x if 'Unnamed' not in a]) for x in GEE.columns.values]


