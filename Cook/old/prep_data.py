for (idx,e) in enumerate(edges):
    if(len(e) != 2):
        raise Exception('edge %d: %s is not length 2'%(idx, e))
    if((e[0]%1 != 0) | (e[1]%1 != 0)):
        raise Exception('edge %d: %s contains a non-integer'%(idx, e))
    if(e[0] < e[1]):
        edges[idx] = (int(e[0]), int(e[1]))
    elif(e[0] == e[1]):
        edges[idx] = (0, 0)
    elif(e[0] > e[1]):
        edges[idx] = (int(e[1]), int(e[0]))

edges = set(edges)
edges.discard((0, 0))
edges = np.array(sorted(list(edges))).astype('uint16')
#display(edges)

if num_vtds != edges.max()+1:
    raise Exception('length of population vector != max vtd in edge list.  Is there an isolated vtd?')

vtd0 = edges[:,0]
vtd1 = edges[:,1]
num_edges = len(vtd0)

nbrs0 = [[] for v in range(num_vtds)]
for e in edges:
    nbrs0[e[0]].append(e[1])
    nbrs0[e[1]].append(e[0])
nbrs = [np.array(list(n)).astype('uint16') for n in nbrs0]
del nbrs0

degree = np.array([len(n) for n in nbrs])
nbr_rng = np.insert(np.cumsum(degree),0,0)

def print_data():
    #print("edges\n%s"%edges)
    print("num_edges = %s"%num_edges)
    print("vtd0 = %s"%vtd0)
    print("vtd1 = %s"%vtd1)
    print("num_vtds = %s"%num_vtds)
    print("neighbors = %s"%nbrs)
    print("degree = %s"%degree)
    print("neighbor index range = %s"%nbr_rng)
    print(DIV2)
    print(DIV2)
    
#print_data()