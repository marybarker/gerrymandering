import os
os.chdir('/home/odin/Documents/gerrymandering/gerrymandering/NewHampshire')
from redistricting import *
import plot_stuff

#print(os.getcwd())
os.chdir('./HarvardData')
g = plot_stuff.package_vtds()

steps = 100
for j in range(10):
    starting_state = contiguousStart()
    starting_state.columns = ['key', 'value']
    foldername = 'rangledangle' + str(j+10) + '/'
    runningState = (starting_state.copy(), 1)
    plot_stuff.color_these_states(g, [runningState], foldername, 0)

    for i in range(10):
        runningState = MH(runningState[0], steps, neighbor, goodness, switchDistrict)
        plot_stuff.color_these_states(g, [runningState], foldername, i+1)



