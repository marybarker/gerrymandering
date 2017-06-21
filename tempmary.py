def neighbor(state):
    
    global adjacencyFrame, metrics
    newstate = state.copy()
    newmetrics = metrics.copy()
    
    missingdist = set.difference(set(range(ndistricts)), set(newstate['value']))
    #If we've blobbed out some districts, we wants to behave differently
    
    if len(missingdist) == 0:
        switchedge = np.random.choice(adjacencyFrame.index[-(adjacencyFrame.isSame == 1)])
        
        if random.random() < 0.5:
            lownode      = adjacencyFrame.low[switchedge]
            highnode     = adjacencyFrame.high[switchedge]
            templowdist  = adjacencyFrame.lowdist[switchedge]
            temphighdist = adjacencyFrame.highdist[switchedge]
        else:
            lownode      = adjacencyFrame.high[switchedge]
            highnode     = adjacencyFrame.low[switchedge]
            templowdist  = adjacencyFrame.highdist[switchedge]
            temphighdist = adjacencyFrame.lowdist[switchedge]
            
        switchTo = temphighdist
        
        previousVersion = adjacencyFrame[(adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)]
        proposedChanges = previousVersion.copy()
        #want proposedChanges to be a slice of adjacencyFrame where the values could be changing.
        
        newstate.ix[newstate.key == lownode, 'value'] = switchTo
        proposedChanges.ix[proposedChanges.low == lownode, 'lowdist'] = switchTo
        proposedChanges.ix[proposedChanges.high == lownode, 'highdist'] = switchTo
        proposedChanges.isSame = proposedChanges.lowdist == proposedChanges.highdist
        #change values in the state as well as the proposedChanges
        
        #change minority concentration while population is unchanged. 
        newmetrics.ix[templowdist,  'mincon'] -= float(blockstats.ix[lownode, 'population'] * blockstats.ix[lownode, 'mincon']) / np.nansum(newmetrics.ix[templowdist,  'population'])
        newmetrics.ix[temphighdist, 'mincon'] += float(blockstats.ix[lownode, 'population'] * blockstats.ix[lownode, 'mincon']) / np.nansum(newmetrics.ix[temphighdist, 'population'])
        
        #change population
        popchange = blockstats.population[lownode]
        newmetrics.ix[templowdist, 'population']  -= popchange
        newmetrics.ix[temphighdist, 'population'] += popchange
        
        #change bizarreness
        newmetrics.ix[templowdist,'perimeter']  += \
            (sum(previousVersion.length[ (previousVersion.isSame==1) & ((previousVersion.low == lownode) | (previousVersion.high == lownode))]) -\
             sum(previousVersion.length[~(previousVersion.isSame==1) & ((previousVersion.low == lownode) | (previousVersion.high == lownode))]))
        newmetrics.ix[temphighdist, 'perimeter'] += \
            (sum(proposedChanges.length[~(proposedChanges.isSame==1) & ((proposedChanges.low == lownode) | (proposedChanges.high == lownode))]) -\
             sum(proposedChanges.length[ (proposedChanges.isSame==1) & ((proposedChanges.low == lownode) | (proposedChanges.high == lownode))]))
        
        areachange = blockstats.ALAND[lownode] + blockstats.AWATER[lownode]
        newmetrics.ix[templowdist, 'area'] -= areachange
        newmetrics.ix[temphighdist,'area'] += areachange
        
        newmetrics.ix[templowdist, 'bizarreness']  = bizarreness(newmetrics['area'][templowdist], \
                                                              newmetrics['perimeter'][templowdist])
        newmetrics.ix[temphighdist, 'bizarreness'] = bizarreness(newmetrics['area'][temphighdist], \
                                                              newmetrics['perimeter'][temphighdist])
        
        #update boundary information
        newmetrics.ix[temphighdist,'sumAframDiff'] = newmetrics.ix[temphighdist,'sumAframDiff']\
                                                     + np.sum((-proposedChanges.isSame)*proposedChanges.aframdiff)\
                                                     - np.sum((-previousVersion.isSame)*previousVersion.aframdiff)
        newmetrics.ix[templowdist,'sumAframDiff'] = newmetrics.ix[templowdist,'sumAframDiff']\
                                                    - np.sum((-proposedChanges.isSame)*proposedChanges.aframdiff)\
                                                    + np.sum((-previousVersion.isSame)*previousVersion.aframdiff)
        
        newmetrics.ix[temphighdist,'sumHispDiff'] = newmetrics.ix[temphighdist,'sumHispDiff']\
                                                    + np.sum((-proposedChanges.isSame)*proposedChanges.hispdiff)\
                                                    - np.sum((-previousVersion.isSame)*previousVersion.hispdiff)
        newmetrics.ix[templowdist,'sumHispDiff']  = newmetrics.ix[templowdist,'sumHispDiff']\
                                                    - np.sum((-proposedChanges.isSame)*proposedChanges.hispdiff)\
                                                    + np.sum((-previousVersion.isSame)*previousVersion.hispdiff)
        
        newmetrics.ix[temphighdist,'numedges'] = newmetrics.ix[temphighdist,'numedges']\
                                                 + np.sum(-(proposedChanges.isSame))\
                                                 - np.sum(-(previousVersion.isSame))
        newmetrics.ix[templowdist,'numedges']  = newmetrics.ix[templowdist,'numedges']\
                                                 - np.sum(-(proposedChanges.isSame))\
                                                 + np.sum(-(previousVersion.isSame))
        #update contiguousness
        neighborhood = set(proposedChanges.low).union(set(proposedChanges.high))
        nhadj = adjacencyFrame.ix[adjacencyFrame.low.isin(neighborhood) & adjacencyFrame.high.isin(neighborhood), ['low','high','length', 'lowdist', 'highdist']]
        oldContNeighborhoodLow  = contiguousness(   state.loc[neighborhood], templowdist,  nhadj)
        oldContNeighborhoodHigh = contiguousness(   state.loc[neighborhood], temphighdist, nhadj)
        
        nhadj.update(proposedChanges)
        newContNeighborhoodLow  = contiguousness(newstate.loc[neighborhood], templowdist,  nhadj)
        newContNeighborhoodHigh = contiguousness(newstate.loc[neighborhood], temphighdist, nhadj)
        
        if ((oldContNeighborhoodLow != newContNeighborhoodLow)|(oldContNeighborhoodHigh != newContNeighborhoodHigh)):
            tempframe = adjacencyFrame.copy()
            tempframe.update(proposedChanges)
            tempframe.lowdist  = tempframe.lowdist.astype(int)
            tempframe.highdist = tempframe.highdist.astype(int)
            
            if (oldContNeighborhoodLow != newContNeighborhoodLow):
                newmetrics.ix[templowdist, 'contiguousness']  = contiguousness(newstate, templowdist, tempframe)
            else:
                pass
            
            if (oldContNeighborhoodHigh != newContNeighborhoodHigh):
                newmetrics.ix[temphighdist, 'contiguousness'] = contiguousness(newstate, temphighdist, tempframe)
            else:
                pass
        
    else:
        #If there are some districts missing, 
        changenode = newstate.key.sample(1)
        olddist = newstate.value[changenode]
        newdist = list(missingdist)[0]
        newstate.value[newstate.key == changenode] = newdist
        #We want to select one randomly, and make it one of the missing districts
        proposedChanges = adjacencyFrame.loc[(adjacencyFrame.low == changenode) | \
                              (adjacencyFrame.high == changenode)]
        proposedChanges.lowdist[proposedChanges.low == changenode] = newdist
        proposedChanges.highdist[proposedChanges.high == changenode] = newdist
        proposedChanges.isSame = False
        # And none of its adjacencies match anymore.
        
        #change contiguousness
        newmetrics['contiguousness'][olddist] = contiguousness(newstate, olddist)
        
        #change population
        popchange = blockstats.population[changenode]
        newmetrics['population'][olddist] -= popchange
        newmetrics['population'][newdist] += popchange
        
        #change bizarreness
        newmetrics['perimeter'][olddist] = perimeter(newstate, olddist)
        newmetrics['perimeter'][newdist] = perimeter(newstate, newdist)
        
        areachange = blockstats.ALAND[changenode] + blockstats.AWATER[changenode]
        newmetrics['area'][olddist] -= areachange
        newmetrics['area'][newdist] += areachange
        
        newmetrics['bizarreness'][olddist] = bizarreness(newmetrics['area'][olddist], \
                                                              newmetrics['perimeter'][olddist])
        newmetrics['bizarreness'][newdist] = bizarreness(newmetrics['area'][newdist], \
                                                              newmetrics['perimeter'][newdist])
        newmetrics['mincon'][olddist] = minorityConc(newstate, olddist, 'mincon')
        newmetrics['mincon'][newdist] = blockstats.ix[changenode, 'mincon']

    return (newstate, proposedChanges, newmetrics)


def goodness(metrics):
    
    tempStConts  = metrics['contiguousness']
    
    if any([x!=1 for x in tempStConts]):
        return float('-inf')
    
    tempStPops   = metrics['population']
    tempStBiz    = metrics['bizarreness']
    
    modTotalVar = sum([abs(float(x)/totalpopulation - float(1)/ndistricts) for x in tempStPops])/(2*(1-float(1)/ndistricts))
    
    return -3000*modTotalVar - 300*np.nanmean(tempStBiz) - \
            float(max(0, (np.max(tempStPops) - np.min(tempStPops)) - 25000 )**2)/1000000 
    #functions should be written such that the numbers being scaled are between zero and one.


def compare_current_state_to_possible_perturbations(current_state, num_perturbations, steps_per_perturbation, foldername, misc_data =''):
    import datetime
    if not os.path.isdir(foldername):
        os.mkdir(foldername)
    os.chdir(foldername)
    
    metadata = 80*'+'+"\nDate: %s"%str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+\
                '\nSteps per run: %d'%(steps_per_perturbation)+\
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
    


