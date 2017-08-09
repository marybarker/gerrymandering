from scipy.misc import imread
import os
import time
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from itertools import product
from subprocess import call
#import animate
#  (Not on hard drive.  Nicer than ffmpeg)
os.chdir('/home/odin/Documents/gerrymandering/gerrymandering/NorthCarolina')

execfile('../cleanitallup.py')
execfile('../setup_stuff.py')
execfile('setup.py')

temp = product(*([[1, 10, 100]]*3))
distinctParam = [0,1,2,3,4,5,6,7,8,9,10,11,12,15,18,19,20,21,24]

goodnessParams[1] = popDiffScore 
paramList = [[1, x[0], x[1], x[1], x[2], x[2]] for x in temp]

gridrange = [paramList[x] for x in distinctParam]
#gridrange = [paramList[x] for x in [2,6,15]]
samplerate = 10
statereads = 1
numsaves = 250
numreads = numsaves
foldername = "gridyesflat/"
DPI = 300

if "mapgifs" not in os.listdir(foldername):
    os.mkdir(foldername + "mapgifs")

if "CD_color" not in os.listdir(foldername + "mapgifs"):
    os.mkdir(foldername + "mapgifs/CD_color")

if "CD_afram" not in os.listdir(foldername + "mapgifs"):
    os.mkdir(foldername + "mapgifs/CD_afram")

if "CD_hisp" not in os.listdir(foldername + "mapgifs"):
    os.mkdir(foldername + "mapgifs/CD_hisp")

aframimg = imread("demographicMaps/aframcon300DPI.png")
hispimg  = imread("demographicMaps/hispcon300DPI.png")
colors = colorDict(ndistricts, light = 0.85)

"""
The following section fills up folders with pngs of the appropriate maps.

"""

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
                         district_boundary=[thisgeom], linewidth = 0.2)
            
            #color_by_rgb(g, {vtd : colors[thisstate.value[vtd]] for vtd in thisstate.index}, 
            #             foldername + "mapgifs/CD_color/" + subfoldername + "state%04d/fig%04d.png"%(startingpoint,j),
            #             district_boundary=False)
            
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
            plt.close(fig)
            
            #Make map of outline with Hispanic concentrations visible.
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlim(g['xlim'])
            ax.set_ylim(g['ylim'])
            plt.imshow(hispimg, zorder=0, extent=[g['xlim'][0], g['xlim'][1], g['ylim'][0], g['ylim'][1]])
            
            patch_edges(thisgeom)
            
            ax.set_aspect(1.0)
            
            plt.savefig(foldername + "mapgifs/CD_hisp/" + subfoldername + "state%04d/fig%04d.png"%(startingpoint,j+1),
                        dpi=DPI)
            plt.clf()
            plt.close(fig)
            
            
            print("...")
        print("Created maps for grid%04d.%04d.%04d state %04d"%(weights[1],weights[3],weights[5], startingpoint))


"""
The following section makes repeated calls for ffmpeg to stitch the pngs created above into mp4s.
The original plan was to make gifs, but they were ugly  ^_^

"""
for startingpoint in range(statereads):
    for point in range(len(gridrange)):
        weights = gridrange[point]
        goodnessWeights = np.array(weights)
        subfoldername = "grid%04d.%04d.%04d/"%(weights[1],weights[3],weights[5])
        call(("ffmpeg -framerate 5 -i ./gridyesflat/mapgifs/CD_afram/grid%04d.%04d.%04d/state%04d/fig%%03d1.png ./gridyesflat/mapgifs/afram%04d.%04d.%04d_state%04d.mp4"%(weights[1],weights[3],weights[5], startingpoint, weights[1],weights[3],weights[5], startingpoint)).split())
        call(("ffmpeg -framerate 5 -i ./gridyesflat/mapgifs/CD_hisp/grid%04d.%04d.%04d/state%04d/fig%%03d1.png ./gridyesflat/mapgifs/hisp%04d.%04d.%04d_state%04d.mp4"%(weights[1],weights[3],weights[5], startingpoint, weights[1],weights[3],weights[5], startingpoint)).split())
        call(("ffmpeg -framerate 5 -i ./gridyesflat/mapgifs/CD_color/grid%04d.%04d.%04d/state%04d/fig%%03d0.png ./gridyesflat/mapgifs/color%04d.%04d.%04d_state%04d.mp4"%(weights[1],weights[3],weights[5], startingpoint, weights[1],weights[3],weights[5], startingpoint)).split())

#If we really want gifs, we can make higher quality gifs by first creating a palette
#    ffmpeg -framerate 3 -i ./fig%03d0.png -vf fps=10,scale=320:-1:flags=lanczos,palettegen palette.png
#And then we can create a gif with that palette
#    ffmpeg -framerate 4 -i ./fig%03d0.png -i palette.png -filter_complex "scale=1080:-1:flags=lanczos[x];[x][1:v]paletteuse" output.gif

