#%%time
#%%capture
from helper.setup import *
from helper.pycuda_setup import *

gridDims = (1,1,1)
blockDims = (1,1,1)
max_steps = 2
num_sims = np.product([gridDims,blockDims])


num_clrs = 3
num_vtds = 10
vtd_pop = np.random.randint(20,size=num_vtds).astype('uint16')
num_vtds = len(vtd_pop)

#If given neighbors list
#edges = [(v,w) for (v,N) in enumerate nbrs for w in N if v<w]

num_edges = int((num_vtds*(num_vtds-1)/2)*.75)
edges = np.random.choice(num_vtds, size=[num_edges,2]).tolist()

exec(open("prep_data_new.py").read())
#exec(open("connected_start_gpu.py").read())
exec(open("connected_start_array_new.py").read())
# print(ctd_clr_arr.shape)
print(vtd_clr_arr)
exec(open("main_new.py").read())