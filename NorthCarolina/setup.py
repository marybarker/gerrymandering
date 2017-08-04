import ast

stateSHORT = 'NC'
ndistricts = 13

blockstats = pd.read_csv('./vtdstats.csv', dtype={"VTD":str}).rename(columns={"ALAND10":"ALAND", "AWATER10":"AWATER"})
blockstats = blockstats.set_index(blockstats.ID)
#GEE = pd.read_csv("RegGenderAgeEthnicity.csv")
PRE = pd.read_csv("TotalPopRaceAndEthnicity.csv")
#VAP = pd.read_csv("VotingAgePopulation.csv")
RPR = pd.read_csv("RegPartyRace.csv")

RPR = RPR.loc[:, ["GEOID10", "VR Total", \
                   "VR: All Dems", "VR:White Dems", "VR: Black Dems", "VR: Asian Dems", "VR: Mult Race Ds", \
                   "VR: All Reps", "VR: White Reps", "VR: Black Reps", "VR: Asian Reps", "VR: Mult Race Rs", \
                   "VR: All Libs", "VR: White Libs", "VR: Black Libs", "VR: Asian Libs", "VR: Mult Race Ls", \
                   "VR:All Unaf.", "VR: White Unafil.", "VR: Black Unafil.", "VR: Asian Unafil.", "VR: Other Unafil."
                  ]]
#totalVR = RPR.loc[:, ['VR: All Dems', 'VR: All Reps', 'VR: All Libs', 'VR:All Unaf.']].sum(axis=1)
totalVR = RPR['VR Total']
RPR['percentDem'] = RPR['VR: All Dems'] / totalVR
RPR['percentRep'] = RPR['VR: All Reps'] / totalVR
RPR['percentOther'] = 1.0 - RPR.loc[:, ['percentDem', 'percentRep']].sum(axis=1)

RPR['DemPercentWhite'] = RPR['VR:White Dems'] / RPR['VR: All Dems']
RPR['DemPercentBlack'] = RPR['VR: Black Dems'] / RPR['VR: All Dems']
RPR['DemPercentAsian'] = RPR['VR: Asian Dems'] / RPR['VR: All Dems']
RPR['DemPercentOther'] = 1.0 - RPR.loc[:, ['DemPercentWhite', 'DemPercentBlack', 'DemPercentAsian']].sum(axis=1)
RPR['RepPercentWhite'] = RPR['VR: White Reps'] / RPR['VR: All Reps']
RPR['RepPercentBlack'] = RPR['VR: Black Reps'] / RPR['VR: All Reps']
RPR['RepPercentAsian'] = RPR['VR: Asian Reps'] / RPR['VR: All Reps']
RPR['RepPercentOther'] = 1.0 - RPR.loc[:, ['RepPercentWhite', 'RepPercentBlack', 'RepPercentAsian']].sum(axis=1)
RPR = RPR.loc[:, ['GEOID10', 'percentDem', 'percentRep', 'percentOther', \
                  'DemPercentWhite', 'DemPercentBlack', 'DemPercentAsian', 'DemPercentOther', \
                  'RepPercentWhite', 'RepPercentBlack', 'RepPercentAsian', 'RepPercentOther'
                  ]]

PRE = PRE.loc[:, ['Total', 'GEOID10', '% Total Black', '%  Hisp ']].rename(columns={'Total':'population', '% Total Black':'aframcon', '%  Hisp ':'hispcon'})
PRE["mincon"] = PRE.loc[:, ['aframcon', 'hispcon']].sum(axis=1)
PRE.mincon[PRE.mincon > 1.0] = 1.0

votingResultsData = pd.read_csv("NC_2012_with_GEOID.csv").loc[:, ['GEOID10', 'g2012_USH_dv', 'g2012_USH_rv', 'g2012_USH_tv']]
blockstats = pd.merge(votingResultsData, pd.merge(blockstats, pd.merge(PRE, RPR, on='GEOID10'), on='GEOID10'), on='GEOID10')
blockstats = blockstats.set_index(blockstats.ID).sort_index()

totalpopulation = sum(blockstats.population)
nvtd = len(blockstats.VTD)

stateconcentration = np.nansum(blockstats.mincon * blockstats.population) * 1.0 / np.nansum(blockstats.population)
numMajMinDists = int(ndistricts*stateconcentration)


lookup = dict(zip(blockstats.GEOID10, blockstats.ID))
adjacencyFrame = pd.read_csv("PRECINCTconnections.csv")
adjacencyFrame.low = [lookup.get(x) for x in adjacencyFrame.low]
adjacencyFrame.high = [lookup.get(x) for x in adjacencyFrame.high]

conlow = pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "aframcon"]], left_on = 'low', right_on = 'ID').aframcon
conhigh= pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "aframcon"]], left_on = 'high', right_on = 'ID').aframcon
adjacencyFrame["aframdiff"] = conhigh - conlow

conlow = pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "hispcon"]], left_on = 'low', right_on = 'ID').hispcon
conhigh= pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "hispcon"]], left_on = 'high', right_on = 'ID').hispcon
adjacencyFrame["hispdiff"] = conhigh - conlow

blockstats["repPop"  ] = blockstats.percentRep*blockstats.population
blockstats["demPop"  ] = blockstats.percentDem*blockstats.population
blockstats["hispPop" ] = blockstats.hispcon*blockstats.population
blockstats["aframPop"] = blockstats.aframcon*blockstats.population

blockstats["numAdjacent"] = [len(adjacencyFrame.index[ (adjacencyFrame.low == i) | (adjacencyFrame.high == i)]) for i in blockstats.index]
blockstats["totalBorder"] = [sum(adjacencyFrame.length[(adjacencyFrame.low == i) | (adjacencyFrame.high == i)]) for i in blockstats.index]

mutableBlockStats = {}

g = package_vtds("./precinct/precinct.shp", "GEOIDToIDNUM.csv", ['GEOID10'])

metrics = pd.DataFrame()

perimeterNC = pd.read_csv("./currentNCBoundary.csv")
adjacencyGeomFrame = pd.read_csv("./adjacencyGeoms.csv")

adjacencyGeomFrame.ix[:, "geom"] = [ast.literal_eval(x) for x in adjacencyGeomFrame.geom]



