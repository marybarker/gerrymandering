stateSHORT = 'NC'

blockstats = pd.read_csv('./vtdstats.csv', dtype={"VTD":str})

GEE = pd.read_csv("RegGenderAgeEthnicity.csv")
PRE = pd.read_csv("TotalPopRaceAndEthnicity.csv")
VAP = pd.read_csv("VotingAgePopulation.csv")
RPR = pd.read_csv("RegPartyRace.csv")

blockstats = pd.merge(blockstats, GEE, on='GEOID10')

adjacencyFrame = pd.read_csv("PRECINCTconnections.csv", index=False)
