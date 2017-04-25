singles = [vtd for vtd in blockstats.VTD if adjacencyFrame.loc[(adjacencyFrame.low == vtd) | (adjacencyFrame.high == vtd)].shape[0] == 1]

dungus = singles + \
         [adjacencyFrame.high[adjacencyFrame.low == x].item() for x in singles if adjacencyFrame.high[adjacencyFrame.low == x].shape[0] != 0] + \
         [adjacencyFrame.low[adjacencyFrame.high == x].item() for x in singles if adjacencyFrame.low[adjacencyFrame.high == x].shape[0] != 0]

stupidstate = pd.DataFrame({'key':blockstats.VTD, 'value': 1*blockstats.VTD.isin(singles) +1*blockstats.VTD.isin(dungus)})

color_these_states(g, [(stupidstate,0)], './', 0)


