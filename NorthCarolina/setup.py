stateSHORT = 'NC'
blockstats = pd.read_csv('./vtdstats.csv', dtype={"VTD":str})

#GEE = pd.read_csv("RegGenderAgeEthnicity.csv")
PRE = pd.read_csv("TotalPopRaceAndEthnicity.csv")
#VAP = pd.read_csv("VotingAgePopulation.csv")
RPR = pd.read_csv("RegPartyRace.csv")

RPR = RPR.loc[:, ["GEOID10", \
                   "VR: All Dems", "VR: White Dems", "VR: Black Dems", "VR: Asian Dems", "VR: Mult Race Ds", \
                   "VR: All Reps", "VR: White Reps", "VR: Black Reps", "VR: Asian Reps", "VR: Mult Race Rs", \
                   "VR: All Libs", "VR: White Libs", "VR: Black Libs", "VR: Asian Libs", "VR: Other Libs", \
                   "VR:All Unaf.", "VR: White Unafil.", "VR: Black Unafil.", "VR: Asian Unafil.", "VR: Other Unafil."
                  ]]
totalVR = RPR.loc[:, ['VR: All Dems', 'VR: All Reps', 'VR: All Libs', 'VR:All Unaf.']].sum(axis=1)
RPR['percentDem'] = RPR['VR: All Dems'] / totalVR
RPR['percentRep'] = RPR['VR: All Reps'] / totalVR
RPR['percentOther'] = 1.0 - RPR.loc[:, ['percentDem', 'percentRep']].sum(axis=1)

RPR['DemPercentWhite'] = RPR['VR: White Dems'] / RPR['VR: All Dems']
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

PRE = PRE.loc[:, ['GEOID10', '% Total Black', '%  Hisp ']].rename(columns={'% Total Black':'aframcon', '%  Hisp ':'hispcon'})
PRE["mincon"]   = PRE.loc[:, ['aframcon', 'hispcon']].sum(axis=1)

blockstats = pd.merge(blockstats, pd.merge(PRE, RPR, on='GEOID10'), on='GEOID10')

# want to look at minority concentration in each county 
# cross-referenced with the voting habits of both the county and of the minority pop

adjacencyFrame = pd.read_csv("PRECINCTconnections.csv")
g = package_vtds("./precinct/precinct.shp", "GEOIDToIDNUM.csv", ['GEOID10'])

