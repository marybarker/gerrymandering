stateSHORT = 'PA'

#blockstats = pd.read_csv('./noIslandsVTDStats.csv').merge(pd.read_csv('GEOIDToIDNUM.csv'),left_on="VTD",right_on="GEOID")
#blockstats.rename(columns={"IDNUM":"ID"},inplace=True)
#blockstats=blockstats.merge(pd.read_csv('noIslandsApportionmentData.csv'),on="ID")
#blockstats['VTD']=blockstats['VTD_x'].copy()
blockstats = pd.read_csv('./precinctstats.csv').rename(columns={"Unnamed: 0":"ID"}).merge(pd.read_csv('./noIslandsApportionmentData.csv'), on="ID")
blockstats['VTD']=blockstats['VTD.1'].copy()
for key in blockstats.keys():
    if (key[:4] == "ID.") or (key[:9] == 'Unnamed: '):
        blockstats = blockstats.drop(key, 1)
    if (key == "POP100"):
        blockstats.rename(columns={key:"population"},inplace=True)
    if (key[:3]=='VTD') and (len(key)>4):
        blockstats=blockstats.drop(key,1)

blockstats = blockstats.set_index(blockstats.ID).sort_index()
totalpopulation = sum(blockstats.population)

conccolumn = 'hispcon'
blockstats['hispcon'] = 0
blockstats['aframcon'] = blockstats.B02009e1.values/blockstats.B02001e1.values
stateconcentration = np.nansum(blockstats.B02009e1.values)/np.nansum(blockstats.B02001e1.values)
numMajMinDists = min(totalpopulation/1975932, int(18*stateconcentration)) 
#number found from blobbing out majorityminority districts once.  This result is suspect.

blockstats['mincon'] = blockstats.loc[:,['aframcon','hispcon']].sum(axis=1)
blockstats.mincon[blockstats.mincon > 1.0] = 1.0

cdtable = pd.read_csv('../cdbystate1.txt', '\t')
ndistricts = int(cdtable[cdtable['STATE']==stateSHORT].CD)
nvtd = len(blockstats.VTD)

adjacencyFrame = pd.read_csv('noIslandsPrecinctConnections.csv')

lookup=dict(zip(blockstats.VTD,blockstats.ID))
adjacencyFrame['low']=[lookup.get(x) for x in adjacencyFrame.low]
adjacencyFrame['high']=[lookup.get(x) for x in adjacencyFrame.high]

conlow = pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "aframcon"]], left_on = 'low', right_on = 'ID').aframcon
conhigh= pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "aframcon"]], left_on = 'high', right_on = 'ID').aframcon
adjacencyFrame["aframdiff"] = conhigh - conlow

conlow = pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "hispcon"]], left_on = 'low', right_on = 'ID').hispcon
conhigh= pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "hispcon"]], left_on = 'high', right_on = 'ID').hispcon
adjacencyFrame["hispdiff"] = conhigh - conlow

mutableBlockStats = {}

g = package_vtds("precinct/precinct.shp", "GEOIDToIDNUM.csv")



