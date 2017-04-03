ambition = 10

def ambitiousNeighbor_old(state):
    
    newstate = state.copy()
    
    switchedges = np.random.choice(adjacencyFrame.index[-(adjacencyFrame.isSame == 1)], ambition, replace = False)
    for switchedge in switchedges:
        
        missingdist = set.difference(set(range(ndistricts)), set(newstate['value']))
        #If we've blobbed out some districts, we wants to behave differently
        
        if len(missingdist) == 0:
            lownode  = adjacencyFrame.low[switchedge]
            highnode = adjacencyFrame.high[switchedge]
            #Randomly choose an adjacency.  Find the low node and high node for that adjacency.
    
            if random.random() < 0.5:
                newstate.value[newstate.key ==  lownode] = (newstate[newstate.key == highnode].value).item()
                checks = adjacencyFrame.index[((adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)) & \
                                              (-(adjacencyFrame.isSame == 1))]
                adjacencyFrame.isSame[((adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)) & \
                                      adjacencyFrame.isSame] = False
                adjacencyFrame.isSame[checks] = [( newstate.value[newstate.key == adjacencyFrame.low[j]].item() == \
                                                   newstate.value[newstate.key == adjacencyFrame.high[j]].item() ) for j in checks]
            else:
                newstate.value[newstate.key == highnode] = (newstate[newstate.key ==  lownode].value).item()
                checks = adjacencyFrame.index[((adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)) & \
                                              (-(adjacencyFrame.isSame == 1))]
                adjacencyFrame.isSame[((adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)) & \
                                      adjacencyFrame.isSame] = False
                adjacencyFrame.isSame[checks] = [( newstate.value[newstate.key == adjacencyFrame.low[j]].item() == \
                                                   newstate.value[newstate.key == adjacencyFrame.high[j]].item() ) for j in checks]
            #We want to assign both nodes the same value, and there's a 50% chance for each value being chosen.
        else:
            #If there are some districts missing, 
            changenode = newstate.key.sample(1)
            newstate.value[newstate.key == changenode] = list(missingdist)[0]
            #We want to select one randomly, and make it one of the missing districts
            adjacencyFrame.isSame[(adjacencyFrame.low == changenode) | \
                                  (adjacencyFrame.high == changenode)] = False
            # And none of its adjacencies match anymore.
    return newstate




def ambitiousNeighbor_old(state):
    
    newstate = state.copy()
    
    switchedges = np.random.choice(adjacencyFrame.index[-(adjacencyFrame.isSame == 1)], ambition, replace = False)
    for switchedge in switchedges:
        
        missingdist = set.difference(set(range(ndistricts)), set(newstate['value']))
        #If we've blobbed out some districts, we wants to behave differently
        
        if len(missingdist) == 0:
            lownode  = adjacencyFrame.low[switchedge]
            highnode = adjacencyFrame.high[switchedge]
            #Randomly choose an adjacency.  Find the low node and high node for that adjacency.
    
            if random.random() < 0.5:
                newstate.value[newstate.key ==  lownode] = (newstate[newstate.key == highnode].value).item()
                checks = adjacencyFrame.index[((adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)) & \
                                              (-(adjacencyFrame.isSame == 1))]
                adjacencyFrame.isSame[((adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)) & \
                                      adjacencyFrame.isSame] = False
                adjacencyFrame.isSame[checks] = [( newstate.value[newstate.key == adjacencyFrame.low[j]].item() == \
                                                   newstate.value[newstate.key == adjacencyFrame.high[j]].item() ) for j in checks]
            else:
                newstate.value[newstate.key == highnode] = (newstate[newstate.key ==  lownode].value).item()
                checks = adjacencyFrame.index[((adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)) & \
                                              (-(adjacencyFrame.isSame == 1))]
                adjacencyFrame.isSame[((adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)) & \
                                      adjacencyFrame.isSame] = False
                adjacencyFrame.isSame[checks] = [( newstate.value[newstate.key == adjacencyFrame.low[j]].item() == \
                                                   newstate.value[newstate.key == adjacencyFrame.high[j]].item() ) for j in checks]
            #We want to assign both nodes the same value, and there's a 50% chance for each value being chosen.
        else:
            #If there are some districts missing, 
            changenode = newstate.key.sample(1)
            newstate.value[newstate.key == changenode] = list(missingdist)[0]
            #We want to select one randomly, and make it one of the missing districts
            adjacencyFrame.isSame[(adjacencyFrame.low == changenode) | \
                                  (adjacencyFrame.high == changenode)] = False
            # And none of its adjacencies match anymore.
    return newstate

def ambitiousNeighbor_old(state):
    
    newstate = state.copy()
    
    switchedges = np.random.choice(adjacencyFrame.index[-(adjacencyFrame.isSame == 1)], ambition, replace = False)
    for switchedge in switchedges:
        
        missingdist = set.difference(set(range(ndistricts)), set(newstate['value']))
        #If we've blobbed out some districts, we wants to behave differently
        
        if len(missingdist) == 0:
            lownode  = adjacencyFrame.low[switchedge]
            highnode = adjacencyFrame.high[switchedge]
            #Randomly choose an adjacency.  Find the low node and high node for that adjacency.
    
            if random.random() < 0.5:
                newstate.value[newstate.key ==  lownode] = (newstate[newstate.key == highnode].value).item()
                checks = adjacencyFrame.index[((adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)) & \
                                              (-(adjacencyFrame.isSame == 1))]
                adjacencyFrame.isSame[((adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)) & \
                                      adjacencyFrame.isSame] = False
                adjacencyFrame.isSame[checks] = [( newstate.value[newstate.key == adjacencyFrame.low[j]].item() == \
                                                   newstate.value[newstate.key == adjacencyFrame.high[j]].item() ) for j in checks]
            else:
                newstate.value[newstate.key == highnode] = (newstate[newstate.key ==  lownode].value).item()
                checks = adjacencyFrame.index[((adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)) & \
                                              (-(adjacencyFrame.isSame == 1))]
                adjacencyFrame.isSame[((adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)) & \
                                      adjacencyFrame.isSame] = False
                adjacencyFrame.isSame[checks] = [( newstate.value[newstate.key == adjacencyFrame.low[j]].item() == \
                                                   newstate.value[newstate.key == adjacencyFrame.high[j]].item() ) for j in checks]
            #We want to assign both nodes the same value, and there's a 50% chance for each value being chosen.
        else:
            #If there are some districts missing, 
            changenode = newstate.key.sample(1)
            newstate.value[newstate.key == changenode] = list(missingdist)[0]
            #We want to select one randomly, and make it one of the missing districts
            adjacencyFrame.isSame[(adjacencyFrame.low == changenode) | \
                                  (adjacencyFrame.high == changenode)] = False
            # And none of its adjacencies match anymore.
    return newstate


def ambitiousNeighbor(state):
    
    newstate = state.copy()
    
    if (not adjacencyFrame.isSame.dtype == bool):
        adjacencyFrame.isSame = adjacencyFrame.isSame == 1
    
    switchedges = np.random.choice(adjacencyFrame.index[-adjacencyFrame.isSame], ambition, replace = False)
    for switchedge in switchedges:
        
        missingdist = set.difference(set(range(ndistricts)), set(newstate['value']))
        #If we've blobbed out some districts, we wants to behave differently
        
        if len(missingdist) == 0:
            lownode  = adjacencyFrame.low[switchedge]
            highnode = adjacencyFrame.high[switchedge]
            #Randomly choose an adjacency.  Find the low node and high node for that adjacency.
    
            if random.random() < 0.5:
                newstate.value[newstate.key ==  lownode] = (newstate[newstate.key == highnode].value).item()
                checks = adjacencyFrame.index[((adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)) & \
                                              (-adjacencyFrame.isSame]
                adjacencyFrame.isSame[((adjacencyFrame.low == lownode) | (adjacencyFrame.high == lownode)) & \
                                      adjacencyFrame.isSame] = False
                adjacencyFrame.isSame[checks] = [( newstate.value[newstate.key == adjacencyFrame.low[j]].item() == \
                                                   newstate.value[newstate.key == adjacencyFrame.high[j]].item() ) for j in checks]
            else:
                newstate.value[newstate.key == highnode] = (newstate[newstate.key ==  lownode].value).item()
                checks = adjacencyFrame.index[((adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)) & \
                                              (-adjacencyFrame.isSame)]
                adjacencyFrame.isSame[((adjacencyFrame.low == highnode) | (adjacencyFrame.high == highnode)) & \
                                      adjacencyFrame.isSame] = False
                adjacencyFrame.isSame[checks] = [( newstate.value[newstate.key == adjacencyFrame.low[j]].item() == \
                                                   newstate.value[newstate.key == adjacencyFrame.high[j]].item() ) for j in checks]
            #We want to assign both nodes the same value, and there's a 50% chance for each value being chosen.
        else:
            #If there are some districts missing, 
            changenode = newstate.key.sample(1)
            newstate.value[newstate.key == changenode] = list(missingdist)[0]
            #We want to select one randomly, and make it one of the missing districts
            adjacencyFrame.isSame[(adjacencyFrame.low == changenode) | \
                                  (adjacencyFrame.high == changenode)] = False
            # And none of its adjacencies match anymore.
    return newstate












