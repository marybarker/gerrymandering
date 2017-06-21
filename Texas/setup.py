stateSHORT = 'TX'

blockstats = pd.read_csv('./vtdstats.csv', dtype={"VTD":str}).merge(pd.read_csv('./noIslandsApportionmentData.csv'), on="ID")
for key in blockstats.keys():
    if (key[:4] == "ID.") or (key[:9] == 'Unnamed: '):
        blockstats = blockstats.drop(key, 1)
    if key == "POP100":
        blockstats = blockstats.rename(columns={"POP100":"population"})
#newindex = [str(blockstats.NAME10[i]) + str(blockstats.VTD[i]) for i in blockstats.index]
#blockstats = blockstats.set_index(pd.Series(newindex))
#blockstats.VTD = blockstats.VTD.astype(str)
blockstats = blockstats.set_index(blockstats.ID)  #ID is now a unique identifier
totalpopulation = sum(blockstats.population)

#conccolumn = 'aframcon'
#blockstats['aframcon'] = blockstats.B02009e1.values/blockstats.B02001e1.values
stateconcentration = np.nansum(blockstats.e_bh.values)/np.nansum(blockstats.e_bh.values)
numMajMinDists = min(totalpopulation/1975932, int(18*stateconcentration))
#number found from blobbing out majorityminority districts once.  This result is suspect.

cdtable = pd.read_csv('../cdbystate1.txt', '\t')
#ndistricts = int(cdtable[cdtable['STATE']==stateSHORT].CD)
ndistricts = 9 # Because we're looking at a subset

nvtd = len(blockstats.VTD)

adjacencyFrame = pd.read_csv('PRECINCTconnections.csv', dtype = {'lo' : str, 'hi' : str, 'low' : str, 'high' : str})
for key in adjacencyFrame.keys():
    if key[:9] == 'Unnamed: ':
        adjacencyFrame = adjacencyFrame.drop(key, 1)
        
replacements = {'lo': 'low', 'hi': 'high'}
adjacencyFrame.columns = [replacements.get(x,x) for x in adjacencyFrame.columns]
    #This setup is superior to our previous ones because we might have other information in adjacencyFrame.
    #We just add more things to replacements as they become relevant.
lookup = pd.read_csv("GEOIDToIDNUM.csv")
lookup = dict(zip(lookup.GEOID, lookup.IDNUM))
adjacencyFrame.low = [int(lookup.get(x)) for x in adjacencyFrame.low.values]
adjacencyFrame.high = [int(lookup.get(x)) for x in adjacencyFrame.high.values]
g = package_vtds("./precinct/precinct.shp", "GEOIDToIDNUM.csv")

###
#Demographic differences across VTD boundaries
###

blockstats["aframcon"] = blockstats.e_blak/blockstats.population
blockstats["hispcon"]  = blockstats.e_hsp/blockstats.population
blockstats["mincon"]   = blockstats.e_bh/blockstats.population
blockstats.ix[blockstats.aframcon.isnull(), "aframcon"] = 0
blockstats.ix[blockstats.hispcon.isnull(),  "hispcon" ] = 0
blockstats.ix[blockstats.mincon.isnull(),   "mincon"  ] = 0

conlow = pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "aframcon"]], left_on = 'low', right_on = 'ID').aframcon
conhigh= pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "aframcon"]], left_on = 'high', right_on = 'ID').aframcon
adjacencyFrame["aframdiff"] = conhigh - conlow

conlow = pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "hispcon"]], left_on = 'low', right_on = 'ID').hispcon
conhigh= pd.merge(adjacencyFrame, blockstats.ix[:, ["ID", "hispcon"]], left_on = 'high', right_on = 'ID').hispcon
adjacencyFrame["hispdiff"] = conhigh - conlow

###
#Overwritten functions because of naming differences
###

def distArea(state, district):
    regionlist = list(state.key[state.value == district])
    return sum(blockstats.Shape_area[blockstats.ID.isin(regionlist)])
