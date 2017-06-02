from .setup import *
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import pycuda.gpuarray as gpuarray
import pycuda.curandom as curandom

#Lessons learned - it is best to "ravel" a numpy array BEFORE sending to GPU.
#The automatic raveling done by pycuda seems to be unpredictable.
#I have hit frustrating and difficult to debug issues when relying on it.
#To be fair, I was also learning how to use pycuda, so this may have been
#my fault.  But I find it safest to ravel first to be safe.

def togpu(obj,dtype=None):
    try:
        if dtype is None:
            return gpuarray.to_gpu(np.asarray(obj).ravel())
        else:
            return gpuarray.to_gpu(np.asarray(obj,dtype=dtype).ravel())
    except:
        raise Exception('to_gpu failed: object must be converable to numpy array')