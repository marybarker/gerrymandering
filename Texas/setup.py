stateSHORT = 'TX'

blockstats = pd.read_csv('./vtdstats.csv')#.merge(pd.read_csv('./noIslandsApportionmentData.csv'), on="VTD")
for key in blockstats.keys():
    if (key[:4] == "VTD.") or (key[:9] == 'Unnamed: '):
        blockstats = blockstats.drop(key, 1)
    if key == "POP100":
        blockstats = blockstats.rename(columns={"POP100":"population"})
#newindex = [str(blockstats.NAME10[i]) + str(blockstats.VTD[i]) for i in blockstats.index]
#blockstats = blockstats.set_index(pd.Series(newindex))
blockstats = blockstats.set_index(blockstats.VTD)  #VTD is now a unique identifier
totalpopulation = sum(blockstats.population)

#conccolumn = 'aframcon'
#blockstats['aframcon'] = blockstats.B02009e1.values/blockstats.B02001e1.values
#stateconcentration = np.nansum(blockstats.B02009e1.values)/np.nansum(blockstats.B02001e1.values)
#numMajMinDists = min(totalpopulation/1975932, int(18*stateconcentration)) 
#number found from blobbing out majorityminority districts once.  This result is suspect.

cdtable = pd.read_csv('../cdbystate1.txt', '\t')
#ndistricts = int(cdtable[cdtable['STATE']==stateSHORT].CD)
ndistricts = 7 # Because we're looking at a subset

nvtd = len(blockstats.VTD)

adjacencyFrame = pd.read_csv('PRECINCTconnections.csv', dtype = {'lo' : str, 'hi' : str, 'low' : str, 'high' : str})
for key in adjacencyFrame.keys():
    if key[:9] == 'Unnamed: ':
        adjacencyFrame = adjacencyFrame.drop(key, 1)
        
replacements = {'lo': 'low', 'hi': 'high'}
adjacencyFrame.columns = [replacements.get(x,x) for x in adjacencyFrame.columns]
    #This setup is superior to our previous ones because we might have other information in adjacencyFrame.
    #We just add more things to replacements as they become relevant.

g = package_vtds("./VTDS_of_Interest.shp")

