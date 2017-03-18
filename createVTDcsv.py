import pandas as pd
import os

#Read lookup table for block to VTD
lookuptable = pd.read_csv('../HarvardData/BlockAssign_ST33_NH_VTD.txt')

#Read data from blockshapefiles table, 
blockinfo = pd.read_csv('../NewHampshireBlockShapefiles/NewHampshire.csv')

#Read voter history
voteinfo = pd.read_csv('../HarvardData/nh_vote_merge.tab','\t')

vtdpop = []
vtdwhite = []
vtdrep = []
vtddem = []

#and add the things together
vtdname = list(set(lookuptable.DISTRICT))
for i in range(len(vtdname)):
    
    dist = vtdname[i]
    subframe1 = blockinfo.loc[lookuptable.BLOCKID[lookuptable.DISTRICT == dist].index]
    subframe2 = voteinfo.loc[voteinfo.vtdst10 == dist]
    
    vtdpop.append(np.nansum(subframe1.POP100))
    vtdwhite.append(np.nanmean(subframe1.POPWHITE))
    
    try:
        repvotes = int(subframe2.gov06_coburn_rep.all())
    except:
        repvotes = 0
        
    try:
        demvotes = int(subframe2.gov06_lynch_dem.all())
    except:
        demvotes = 0
    
    if (demvotes == 0) & (repvotes == 0):
        vtdrep.append(0)
        vtddem.append(0)
    else:
        vtdrep.append(float(repvotes)/(repvotes + demvotes))
        vtddem.append(float(demvotes)/(repvotes + demvotes))

    

vtdstats = pd.DataFrame({
    'VTD': vtdname,
    'population': vtdpop,
    'whiteperc' : vtdwhite,
    'repvotes' : vtdrep,
    'demvotes' : vtddem})
    

#write that to a new csv file for reading
vtdstats.to_csv('../HarvardData/NHVTDstats.csv')
