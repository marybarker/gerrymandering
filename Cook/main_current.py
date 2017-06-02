from helper.setup import *
from helper.pycuda_setup import *



kernel_code = """
    #include <stdio.h>    
    __device__ void to_shared(ushort *A, ushort *A_shr, ushort size, ushort threadIdLoc) {
        const ushort threads_per_block = blockDim.z * blockDim.y * blockDim.x;
        ushort idx;
        ushort len = size / threads_per_block;
        ushort idx_start = len*threadIdLoc;
        
        for(idx = idx_start; idx < idx_start + len; idx++) {
            A_shr[idx] = A[idx];
        }
        idx = len*threads_per_block + threadIdLoc;
        if(idx < size) {
            A_shr[idx] = A[idx];
        }
        __syncthreads();
    }


    __device__ void print_shared(ushort *A_shr, ushort size) {
        uint idx;        
        for(idx = 0; idx < size; idx++) {
            printf("A_shr[%d] = %d\\n",idx,A_shr[idx]);
        }
    }

    __device__ void get_components_local(ushort *nbrs, ushort *nbr_rng, ushort *color, ushort *component, ushort *num_comps, ushort w) { 
        printf("num vtds = %d, w = %d, num_comps = %d\\n",NUM_VTDS, w, (*num_comps));
        ushort u, v, idx, orig_comp, old_comp, new_comp;
        
        orig_comp = component[w];
        new_comp = (*num_comps);
        component[w] = new_comp;
        printf("num vtds = %d, w = %d, num_comps = %d\\n",NUM_VTDS, w, (*num_comps));
        printf("\\n\\nTrying to assign vtd %d to a connected component.  Previously in component %d\\n", w, orig_comp);
        printf("There are currently %d connected components\\n",(*num_comps));
        for(idx = nbr_rng[w]; idx < nbr_rng[w+1]; idx++) {
            v = nbrs[idx];
            printf("\\nvtd %d, color %d, component %d\\n", w, color[w], component[w]);
            printf("Checking neighbor vtd %d, color %d, component %d\\n", v, color[v], component[v]);
            if(color[v] == color[w]) {
                printf("Colors match\\n");
                if(component[v] < new_comp) {
                    new_comp = component[v];
                    printf("Smaller component found!!\\nvtd %d: component %d -> %d\\n", w, component[w], new_comp);
                    component[w] = new_comp;
                }
                else {
                    printf("Larger component found.  Not reassigning.");
                }
            }
            if(new_comp == 0){
                break;
            }
        }
        component[w] = new_comp;
        if(new_comp == (*num_comps)) {
            printf("Incrementing num_comps\\n");
            (*num_comps)++;
        }        


        char found = 0;
        if(new_comp != orig_comp) {
            printf("\\nvtd %d switched from component %d to component %d.  Let's check if somebody else is still in component %d\\n",w, orig_comp, new_comp, orig_comp);
            for(v = 0; v < NUM_VTDS; v++) {
                if(component[v] == orig_comp) {
                    printf("Yes, vtd %d is still in component %d\\n", v, component[v]);
                    found = 1;
                    break;
                }
            }
            if(found == 0) {
                (*num_comps)--;                
                printf("No, component %d is now empty.  I will relabel every vtd in component %d as component %d\\n", orig_comp, (*num_comps), orig_comp);
                printf("There are currently %d connected components\\n", (*num_comps));
                for(v = 0; v < NUM_VTDS; v++) {
                    if(component[v] == (*num_comps)){
                        printf("vtd %d: component %d -> %d\\n", v, component[v], old_comp);
                        component[v] = orig_comp;
                    }
                }
            }
        }

        printf("\\nNow we check if this fused >=2 previously distinct components into 1 big component\\n");
        for(idx = nbr_rng[w]; idx < nbr_rng[w+1]; idx++) {
            v = nbrs[idx];
            printf("Checking neighbor vtd %d, color %d, component %d\\n", v, color[v], component[v]);
            if(color[v] == color[w]) {
                if(component[v] > component[w]) {
                    (*num_comps)--;
                    old_comp = component[v];
                    new_comp = component[w];
                    
                    printf("Yes, component %d fused with component %d.\\nRelabeling all vtds in component %d as %d and all vtds in component %d as %d.\\n", old_comp, new_comp, old_comp, new_comp, (*num_comps), old_comp);
                    
                    for(u = 0; u < NUM_VTDS; u++) {
                        if(component[u] == old_comp) {
                            printf("vtd %d: component %d -> %d\\n", u, component[u], new_comp);
                            component[u] = new_comp;                            
                        } else if(component[u] == (*num_comps)) {
                            printf("vtd %d: component %d -> %d", u, component[u], old_comp);
                            component[u] = old_comp;                            
                        }
                    }
                }
            }
        }
        printf("\\nThere are currently %d connected components\\n", (*num_comps));
        for(v = 0; v < NUM_VTDS; v++) {
            printf("%d:%d  ", v, component[v]);
        }
        
        printf("w = %d",w);
        printf("\\nThere are currently %d connected components\\n", (*num_comps));
        printf("END OF ASSIGNMENT PROCESS FOR vtd %d\\n\\n\\n", w);        
    }



    __global__ void gerrymander_kernel(ushort *vtd0, ushort *vtd1, ushort *nbrs, ushort *nbr_rng) {
        const ushort threads_per_block = blockDim.z * blockDim.y * blockDim.x;
        const ushort blockId = blockIdx.z*(gridDim.y*gridDim.x) + blockIdx.y*(gridDim.x) + blockIdx.x;
        const ushort threadIdLoc = threadIdx.z*(blockDim.y*blockDim.x) + threadIdx.y*(blockDim.x) + threadIdx.x;
        const uint threadIdGlob = blockId*threads_per_block + threadIdLoc;
        
        ushort w;
        
        
        __shared__ ushort nbrs_shr[NUM_NBRS];
        to_shared(nbrs, nbrs_shr, NUM_NBRS, threadIdLoc);
        if(threadIdLoc == 0) {
            printf("Printing A_shr = nbrs_shr\\n");
            print_shared(nbrs_shr, NUM_NBRS);
            printf("\\n");
        }

        __shared__ ushort nbr_rng_shr[NUM_VTDS+1];
        to_shared(nbr_rng, nbr_rng_shr, NUM_VTDS+1, threadIdLoc);
        if(threadIdLoc == 0) {
            printf("Printing A_shr = nbr_rng_shr\\n");
            print_shared(nbr_rng_shr, NUM_VTDS+1);
            printf("\\n");
        }
        
        ushort color[NUM_VTDS];
        color[0] = 0;
        color[1] = 0;
        color[2] = 0;
        
        ushort component[NUM_VTDS];
        ushort num_comps = NUM_VTDS;
        for(w = 0; w < NUM_VTDS; w++) {
            component[w] = w;
        }
        
        for(w = 0; w < NUM_VTDS; w++) {
            printf("There are currently %d connected components\\n",num_comps);
            get_components_local(nbrs_shr, nbr_rng_shr, color, component, &num_comps, w);
            printf("There are currently %d connected components\\n",num_comps);
            printf("\\n\\n");
        }
    }
"""


blockDims = (1,1,1)
gridDims = (1,1,1)

num_colors = 2

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

vtd0_gpu = togpu(vtd0,'uint16')
vtd1_gpu = togpu(vtd1,'uint16')
nbrs_gpu = togpu(np.concatenate(nbrs),'uint16')
nbr_rng_gpu = togpu(nbr_rng,'uint16')

macros = {'DIV1':"*" * 80 + "\\n"
         ,'DIV2':"#" * 80 + "\\n"
         ,'NUM_VTDS':num_vtds
         ,'NUM_NBRS':len(nbrs_gpu)
         ,'NUM_COLORS':num_colors
         ,'NUM_EDGES':num_edges
         }
for key, val in macros.items():
    kernel_code = kernel_code.replace(str(key),str(val))

mod = SourceModule(kernel_code);
gerrymander = mod.get_function("gerrymander_kernel")

gerrymander(vtd0_gpu, vtd1_gpu, nbrs_gpu, nbr_rng_gpu, block=blockDims, grid=gridDims, shared=0)