#%%time
#%%capture
from helper.setup import *
from helper.pycuda_setup import *

gridDims = (1,1,1)
blockDims = (50,1,1)
max_steps = 100
num_sims = np.product([gridDims,blockDims])
anneal_rate = 2


num_clrs = 8
num_vtds = 50
vtd_pop = np.random.randint(50,size=num_vtds).astype('uint16')
num_vtds = len(vtd_pop)

#If given neighbors list
#edges = [(v,w) for (v,N) in enumerate nbrs for w in N if v<w]

num_edges = int((num_vtds*(num_vtds-1)/2)*.50)
edges = np.random.choice(num_vtds, size=[num_edges,2]).tolist()
#edges = [(0,1), (1,2), (2,3), (3,4), (4,0), (0,5), (0,6)]

exec(open("prep_data.py").read())

vtd_clr_arr = num_clrs * np.ones([num_sims,num_vtds]).astype('uint16')
num_cmps = np.zeros(num_sims).astype('uint16')
TV = np.zeros(num_sims).astype('float32')
good = np.zeros(num_sims).astype('float32')


exec(open("MH_parallel.py").read())
