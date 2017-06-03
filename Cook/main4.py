import numpy as np
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import pycuda.gpuarray as gpuarray

kernel_code = """
    #include <stdio.h>
    __global__ void main_kernel(uint *A)
    {   
        __shared__ ushort A_shared[%(A_LEN)s];
        if(threadIdx.x < A_LEN){
            A_shared[threadIdx.x] = A[threadIdx.x];
            printf("I'm thread %d.  I see that A[%d] = %d\\n",threadIdx.x,threadIdx.x,A_shared[threadIdx.x]);
        }
    }
"""
n = 40
A = np.arange(n)
#A_gpu = togpu(A,'uint16')
A = togpu(A)

rep = {'A_LEN':n}
kernel_code = kernel_code % rep

#for key, val in macros.items():
#    kernel_code = kernel_code.replace(str(key),str(val))



mod = SourceModule(kernel_code);
main = mod.get_function("main_kernel")

blockDims = (num_vtds,1,1)
gridDims = (1,1,1)

main(A_gpu, block=blockDims, grid=gridDims, shared=0)



num_colors = 3

num_edges = 20
vtd0 = np.zeros(num_edges)
vtd1 = np.zeros(num_edges)
num_edges = len(vtd0)
if len(vtd1) != num_edges:
    raise Exception("Mismatch: length of vtd0 = %d, length of vtd1 = %d"%(len(vtd0), len(vtd1)))

nbrs = [[1]
       ,[0,2]
       ,[1]
       ]
num_vtds = len(nbrs)
degree = np.array([len(n) for n in nbrs])
nbr_rng = np.insert(degree.cumsum(),0,0)







num_vtds_gpu = np.uint16(num_vtds)
num_edges_gpu = np.uint32(num_edges)
num_colors_gpu = np.uint16(num_colors)
vtd0_gpu = togpu(vtd0,'uint16')
vtd1_gpu = togpu(vtd1,'uint16')

nbrs_gpu = togpu(np.concatenate(nbrs),'uint16')
nbr_rng_gpu = togpu(nbr_rng,'uint16')


