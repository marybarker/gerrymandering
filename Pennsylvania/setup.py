stateSHORT = 'PA'

blockstats = pd.read_csv('./noIslandsVTDStats.csv').merge(pd.read_csv('./noIslandsApportionmentData.csv'), on="VTD")
for key in blockstats.keys():
    if (key[:4] == "VTD.") or (key[:9] == 'Unnamed: '):
        blockstats = blockstats.drop(key, 1)

blockstats = blockstats.set_index(blockstats.VTD)
totalpopulation = sum(blockstats.population)

conccolumn = 'aframcon'
blockstats['aframcon'] = blockstats.B02009e1.values/blockstats.B02001e1.values
stateconcentration = np.nansum(blockstats.B02009e1.values)/np.nansum(blockstats.B02001e1.values)
numMajMinDists = min(totalpopulation/1975932, int(18*stateconcentration)) 
#number found from blobbing out majorityminority districts once.  This result is suspect.

cdtable = pd.read_csv('../cdbystate1.txt', '\t')
ndistricts = int(cdtable[cdtable['STATE']==stateSHORT].CD)
nvtd = len(blockstats.VTD)

adjacencyFrame = pd.read_csv('noIslandsPrecinctConnections.csv')
adjacencyFrame = adjacencyFrame.drop('Unnamed: 0', 1)
#adjacencyFrame.columns = ['low', 'high', 'length']

g = package_vtds("precinct/precinct.shp")