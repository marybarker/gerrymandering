

#####
#Population Evenness attempts
#####



start = contiguousStart()
begin = time.time()
final = flatPopulationRun(start, threshold=50000)
end = time.time()
#color_this_state(g, start, "tests/start.png")
color_this_state(g, final, "tests/final.png")
#color_this_state(g, newstate, "tests/dangle.png")
numstates= 1
numsteps = 100
numsaves = 250
temp = product(*([[1, 10, 100]]*3))
distinctParam = [0,1,2,3,4,5,6,7,8,9,10,11,12,15,18,19,20,21,24]

goodnessParams[1] = popDiffScore 
paramList = [[1, x[0], x[1], x[1], x[2], x[2]] for x in temp]
"""
    The way paramList works:
        - It is assumed that globalWeights is of a known length, in this case, 6.
        - I want a low, medium, and high importance option for the goodnesParams,
          but the 3rd and 4th (indices 2 and 3) I wanted to have the same weight,
          and similarly for the 5th and 6th.
        - After numstates have been created with a given system of weights,
          we can do some runs with the next set of weights.
        - This should give us some radically different maps, which we can analyze later.
"""
samplerate = 1
numreads = numsaves
statereads = numstates
#numreads = 1000

#flatpop run
foldername = "gridyesflat/"
#os.mkdir(foldername)
for startingpoint in range(numstates):
    for weights in [paramList[i] for i in distinctParam]:
        subfoldername = "grid%04d.%04d.%04d/"%(weights[1],weights[3],weights[5])
        if subfoldername[:-1] not in os.listdir(foldername):
            os.mkdir(foldername + subfoldername)
        
        goodnessWeights = np.array(weights)
        
        starting_state = contiguousStart()
        updateGlobals(starting_state)
        starting_state.to_csv(foldername+subfoldername + "contiguous_start%04d.csv"%startingpoint, index = False)
        metrics.to_csv(foldername + subfoldername + 'contiguous_metrics%04d.csv'%startingpoint, index = False)
        
        runningState = (flatPopulationRun(starting_state, report = 25000), 1)
        updateGlobals(runningState[0])
        runningState[0].to_csv(foldername+subfoldername + "evenpop_start%04d.csv"%startingpoint, index = False)
        metrics.to_csv(foldername + subfoldername + 'evenpop_metrics%04d.csv'%startingpoint, index = False)
        
        for i in range(numsaves):
            
            runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
            
            runningState[0].to_csv(foldername+subfoldername + "state%04d_save%04d.csv"%(startingpoint, i + 1), index = False)
            
            metrics.to_csv(foldername + subfoldername + 'metrics%04d_save%04d.csv'%(startingpoint, i + 1), index = False)
            
            if i%(numsaves/20) == 0:
                print("run%04d.%04d.%04d.state%04d : %d%% complete."%(weights[1], weights[3], weights[5], startingpoint, 100*i/numsaves))
            
        print("\nFinished run%04d.%04d.%04d.state%04d\n"%(weights[1], weights[3], weights[5], startingpoint))

#flatpop run
foldername = "gridyesflat/"
#os.mkdir(foldername)
for startingpoint in range(numstates):
    for weights in [paramList[i] for i in distinctParam]:
        subfoldername = "grid%04d.%04d.%04d/"%(weights[1],weights[3],weights[5])
        if subfoldername[:-1] not in os.listdir(foldername):
            os.mkdir(foldername + subfoldername)
        
        goodnessWeights = np.array(weights)
        
        starting_state = contiguousStart()
        updateGlobals(starting_state)
        starting_state.to_csv(foldername+subfoldername + "contiguous_start%04d.csv"%startingpoint, index = False)
        metrics.to_csv(foldername + subfoldername + 'contiguous_metrics%04d.csv'%startingpoint, index = False)
        
        runningState = (flatPopulationRun(starting_state, report = 25000), 1)
        updateGlobals(runningState[0])
        runningState[0].to_csv(foldername+subfoldername + "evenpop_start%04d.csv"%startingpoint, index = False)
        metrics.to_csv(foldername + subfoldername + 'evenpop_metrics%04d.csv'%startingpoint, index = False)
        
        for i in range(numsaves):
            
            runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
            
            runningState[0].to_csv(foldername+subfoldername + "state%04d_save%04d.csv"%(startingpoint, i + 1), index = False)
            
            metrics.to_csv(foldername + subfoldername + 'metrics%04d_save%04d.csv'%(startingpoint, i + 1), index = False)
            
            if i%(numsaves/20) == 0:
                print("run%04d.%04d.%04d.state%04d : %d%% complete."%(weights[1], weights[3], weights[5], startingpoint, 100*i/numsaves))
            
        print("\nFinished run%04d.%04d.%04d.state%04d\n"%(weights[1], weights[3], weights[5], startingpoint))

#not flatpop run
foldername = "gridnoflat/"
#os.mkdir(foldername)
for startingpoint in range(numstates):
    #for weights in [paramList[i] for i in distinctParam]:
    for weights in paramList:
        subfoldername = "grid%04d.%04d.%04d/"%(weights[1],weights[3],weights[5])
        if subfoldername[:-1] not in os.listdir(foldername):
            os.mkdir(foldername + subfoldername)
        
        goodnessWeights = np.array(weights)
        
        starting_state = contiguousStart()
        updateGlobals(starting_state)
        starting_state.to_csv(foldername+subfoldername + "contiguous_start%04d.csv"%startingpoint, index = False)
        metrics.to_csv(foldername + subfoldername + 'contiguous_metrics%04d.csv'%startingpoint, index = False)
        
        runningState = (starting_state.copy(), 1)
        for i in range(numsaves):
            
            runningState = MH(runningState[0], numsteps, neighbor, goodness, switchDistrict)
            
            runningState[0].to_csv(foldername+subfoldername + "state%04d_save%04d.csv"%(startingpoint, i + 1), index = False)
            
            metrics.to_csv(foldername + subfoldername + 'metrics%04d_save%04d.csv'%(startingpoint, i + 1), index = False)
            
            if i%(numsaves/20) == 0:
                print("run%04d.%04d.%04d.state%04d : %d%% complete."%(weights[1], weights[3], weights[5], startingpoint, 100*i/numsaves))
            
        print("\nFinished run%04d.%04d.%04d.state%04d\n"%(weights[1], weights[3], weights[5], startingpoint))



#####
#Supplements to cleanitallup.py
#####


#####
#Supplements to setup_stuff.py
#####

highestblack = max(blockstats.aframcon)
black = {blockstats.index[i] : (1 - blockstats.aframcon[i]/highestblack, 1 - 0.3*blockstats.aframcon[i]/highestblack, 1 - blockstats.aframcon[i]/highestblack) for i in blockstats.index}

color_by_rgb(g, black, "tests/aframmapnoedge.png", 0, district_boundary = False)
color_by_rgb(g, black, "tests/aframmapedge.png", 0)


#####
#Supplements to setup.py
#####


#####
#Create gifs for the maps as they develop over time
#####


temp = product(*([[1, 10, 100]]*3))
distinctParam = [0,1,2,3,4,5,6,7,8,9,10,11,12,15,18,19,20,21,24]

goodnessParams[1] = popDiffScore 
paramList = [[1, x[0], x[1], x[1], x[2], x[2]] for x in temp]

gridrange = [paramList[x] for x in distinctParam]
samplerate = 10
statereads = 1
numsaves = 250
numreads = numsaves
foldername = "gridyesflat/"

if "mapgifs" not in os.listdir(foldername):
    os.mkdir(foldername + "mapgifs")

if "CD_color" not in os.listdir(foldername + "mapgifs"):
    os.mkdir(foldername + "mapgifs/CD_color")

if "CD_afram" not in os.listdir(foldername + "mapgifs"):
    os.mkdir(foldername + "mapgifs/CD_afram")

if "CD_hisp" not in os.listdir(foldername + "mapgifs"):
    os.mkdir(foldername + "mapgifs/CD_hisp")

colors = colorDict(ndistricts)

import PIL
from scipy.misc import imread

aframimg = imread("demographicMaps/aframcon300DPI.png")
hispimg  = imread("demographicMaps/hispcon300DPI.png")

for startingpoint in range(statereads):
    for point in range(len(gridrange)):
        weights = gridrange[point]
        goodnessWeights = np.array(weights)
        subfoldername = "grid%04d.%04d.%04d/"%(weights[1],weights[3],weights[5])
        
        if subfoldername[:-1] not in os.listdir(foldername + "mapgifs/CD_color"):
            os.mkdir(foldername + "mapgifs/CD_color/"+ subfoldername)
        if "state%04d"%startingpoint not in os.listdir(foldername + "mapgifs/CD_color/"+ subfoldername):
            os.mkdir(foldername + "mapgifs/CD_color/"+ subfoldername + "state%04d"%startingpoint)
        if subfoldername[:-1] not in os.listdir(foldername + "mapgifs/CD_afram"):
            os.mkdir(foldername + "mapgifs/CD_afram/"+ subfoldername)
        if "state%04d"%startingpoint not in os.listdir(foldername + "mapgifs/CD_afram/"+ subfoldername):
            os.mkdir(foldername + "mapgifs/CD_afram/"+ subfoldername + "state%04d"%startingpoint)
        if subfoldername[:-1] not in os.listdir(foldername + "mapgifs/CD_hisp"):
            os.mkdir(foldername + "mapgifs/CD_hisp/"+ subfoldername)
        if "state%04d"%startingpoint not in os.listdir(foldername + "mapgifs/CD_hisp/"+ subfoldername):
            os.mkdir(foldername + "mapgifs/CD_hisp/"+ subfoldername + "state%04d"%startingpoint)
        
        for j in samplerate*np.arange(numreads/samplerate):
            thisstate = pd.read_csv(foldername+subfoldername + 'state%04d_save%04d.csv'%(startingpoint, j+1))
            updateGlobals(thisstate)
            thisgeom = district_bounds_to_geom()
            
            color_by_rgb(g, {vtd : colors[thisstate.value[vtd]] for vtd in thisstate.index}, 
                         foldername + "mapgifs/CD_color/" + subfoldername + "state%04d/fig%04d.png"%(startingpoint,j),
                         district_boundary=False)
            
            #Make map of outline with African American concentrations visible.
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlim(g['xlim'])
            ax.set_ylim(g['ylim'])
            plt.imshow(aframimg, zorder=0, extent=[g['xlim'][0], g['xlim'][1], g['ylim'][0], g['ylim'][1]])
            
            patch_edges(thisgeom)
            
            ax.set_aspect(1.0)
            
            plt.savefig(foldername + "mapgifs/CD_afram/" + subfoldername + "state%04d/fig%04d.png"%(startingpoint,j+1),
                        dpi=DPI)
            plt.clf()
            del fig
            
            #Make map of outline with Hispanic concentrations visible.
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlim(g['xlim'])
            ax.set_ylim(g['ylim'])
            plt.imshow(aframimg, zorder=0, extent=[g['xlim'][0], g['xlim'][1], g['ylim'][0], g['ylim'][1]])
            
            patch_edges(thisgeom)
            
            ax.set_aspect(1.0)
            
            plt.savefig(foldername + "mapgifs/CD_hisp/" + subfoldername + "state%04d/fig%04d.png"%(startingpoint,j+1),
                        dpi=DPI)
            plt.clf()
            del fig
            
        print("Loaded grid%04d.%04d.%04d state %04d"%(weights[1],weights[3],weights[5], startingpoint))

"""

highestblack = max(blockstats.aframcon)
black = {blockstats.index[i] : (1 - blockstats.aframcon[i]/highestblack, 1 - 0.6*blockstats.aframcon[i]/highestblack, 1 - 0.6*blockstats.aframcon[i]/highestblack) for i in blockstats.index}

highesthisp = max(blockstats.hispcon)
hisp = {blockstats.index[i] : (1 - 0.5*blockstats.hispcon[i]/highesthisp, 1 - blockstats.hispcon[i]/highesthisp, 1 - 0.5*blockstats.hispcon[i]/highesthisp) for i in blockstats.index}

vtds_rgb_dict = black; filename = "./demographicMaps/aframcon300DPI.png";
vtds_rgb_dict = hisp;  filename = "./demographicMaps/hispcon300DPI.png";

#Make frameless plot
xSize = g['xlim'][1] - g['xlim'][0]
ySize = g['ylim'][1] - g['ylim'][0]
fig = plt.figure(figsize=(xSize,ySize))
ax = fig.add_axes([0,0,1,1])

plt.xlim(geom_to_plot['xlim'])
plt.ylim(geom_to_plot['ylim'])

ax.axis('off')
ax.set_autoscale_on(False)

for p in range(len(paths)):
    path = paths[p]
    facecolor = vtds_rgb_dict[names[p]]
    patch = mpatches.PathPatch(path,facecolor=facecolor, edgecolor='black', linewidth=0)
    ax.add_patch(patch)

ax.set_aspect(1.0)


plt.savefig(filename, dpi=DPI)

plt.clf()
del fig
"""



