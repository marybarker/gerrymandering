df1 = pd.DataFrame({'c1':['a','a','b','b'], 'c2':['x','y','x','y'], 'val':0})
df1

df2 = pd.DataFrame({'c1':['a','a','b'], 'c2':['x','y','y'], 'val':[12,31,14]})
df2

df1.merge(df2, 'outer')

df1.update(df2)
df1

numsteps = 10
starting_state = state

temp = dict(zip(starting_state.key, starting_state.value))
lowdists  = adjacencyFrame.low.replace(temp)
highdists = adjacencyFrame.high.replace(temp)
isSame = lowdists==highdists
adjacencyFrame['isSame'] = isSame
adjacencyFrame['lowdist'] = lowdists
adjacencyFrame['highdist'] = highdists

runningState = MH(starting_state, numsteps, neighbor, goodness, switchDistrict)


runningState = MH(starting_state, numsteps, neighbor, goodness, switchDistrict)
