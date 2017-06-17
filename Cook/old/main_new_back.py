from helper.setup import *
from helper.pycuda_setup import *



kernel_code = """
    #include <stdio.h>
    #include <curand.h>
    #include <curand_kernel.h>

extern "C"
{
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
        ushort idx;        
        for(idx = 0; idx < size; idx++) {
            printf("A_shr[%u] = %u\\n", idx, A_shr[idx]);
        }
    }

    __device__ void print_vtd(ushort *vtd_clr, ushort *vtd_cmp, ushort v) {
        printf("vtd%u:clr%u:cmp%u  ", v, vtd_clr[v], vtd_cmp[v]);
    }


    __device__ void print_all(ushort *vtd_clr, ushort *vtd_cmp, ushort num_cmps) {
        ushort c, v;
        for(c = 0; c < num_cmps; c++) {
            printf("\\nComp %u - ",c);
            for(v = 0; v < NUM_VTDS; v++) {
                if(vtd_cmp[v] ==c) {
                    print_vtd(vtd_clr, vtd_cmp, v);
                }
            }
        }
        printf("\\n");        
    }

    __device__ void get_cmps_local(ushort *nbrs, ushort *nbr_rng, ushort *vtd_clr, ushort *vtd_cmp, ushort *num_cmps, ushort w) {         
        printf("Trying to assign vtd%u to a connected cmp.\\n", w);
        print_all(vtd_clr, vtd_cmp, (*num_cmps));
        ushort u, v, idx, orig_cmp, old_cmp, new_cmp;
        
        orig_cmp = vtd_cmp[w];
        new_cmp = (*num_cmps);
        vtd_cmp[w] = new_cmp;                
        printf("vtd%u:clr%u:cmp%u -> cmp%u\\n", w, vtd_clr[w], orig_cmp, new_cmp);
        for(idx = nbr_rng[w]; idx < nbr_rng[w+1]; idx++) {
            v = nbrs[idx];
            printf("Checking neighbor vtd%u:clr%u:cmp %u\\n", v, vtd_clr[v], vtd_cmp[v]);
            if(vtd_clr[v] == vtd_clr[w]) {
                printf("Clrs match\\n");
                if(vtd_cmp[v] < new_cmp) {
                    new_cmp = vtd_cmp[v];
                    printf("Smaller cmp found!!\\n");
                    printf("vtd%u:clr%u:cmp%u -> cmp%u\\n", w, vtd_clr[w], vtd_cmp[w], new_cmp);
                    vtd_cmp[w] = new_cmp;
                }
                else {
                    printf("Larger cmp found.  Not reassigning.\\n");
                }
            }
            if(new_cmp == 0){
                break;
            }
        }
        vtd_cmp[w] = new_cmp;
        if(new_cmp == (*num_cmps)) {
            printf("Incrementing num_cmps\\n");
            (*num_cmps)++;
        }        


        char found = 0;
        if(new_cmp != orig_cmp) {
            printf("\\nvtd%u switched from cmp%u to cmp%u.  Let's check if somebody else is still in cmp%u\\n",w, orig_cmp, new_cmp, orig_cmp);
            for(v = 0; v < NUM_VTDS; v++) {
                if(vtd_cmp[v] == orig_cmp) {
                    printf("Yes, vtd%u:clr%u:cmp%u\\n", v, vtd_clr[v], vtd_cmp[v]);
                    found = 1;
                    break;
                }
            }
            if(found == 0) {
                (*num_cmps)--;                
                printf("No, cmp%u is now empty.  I will relabel every vtd in cmp%u -> cmp%u\\n", orig_cmp, (*num_cmps), orig_cmp);
                for(v = 0; v < NUM_VTDS; v++) {
                    if(vtd_cmp[v] == (*num_cmps)){
                        printf("vtd%u:clr%u:cmp%u -> cmp%u\\n", v, vtd_clr[v], vtd_cmp[v], orig_cmp);
                        vtd_cmp[v] = orig_cmp;
                    }
                }
            }
        }

        printf("\\nNow we check if this fused >=2 previously distinct cmps into 1 big cmp\\n");
        for(idx = nbr_rng[w]; idx < nbr_rng[w+1]; idx++) {
            v = nbrs[idx];
            printf("Checking neighbor vtd%u:clr %u:cmp%u\\n", v, vtd_clr[v], vtd_cmp[v]);
            if(vtd_clr[v] == vtd_clr[w]) {
                if(vtd_cmp[v] > vtd_cmp[w]) {
                    (*num_cmps)--;
                    old_cmp = vtd_cmp[v];
                    new_cmp = vtd_cmp[w];
                    
                    printf("Yes, cmp%u fused with cmp%u.\\nRelabeling all vtds in cmp%u -> cmp%u and all vtds in cmp%u -> cmp%u.\\n", old_cmp, new_cmp, old_cmp, new_cmp, (*num_cmps), old_cmp);
                    
                    for(u = 0; u < NUM_VTDS; u++) {
                        if(vtd_cmp[u] == old_cmp) {
                            printf("vtd %u: cmp %u -> %u\\n", u, vtd_cmp[u], new_cmp);
                            vtd_cmp[u] = new_cmp;                            
                        } else if(vtd_cmp[u] == (*num_cmps)) {
                            printf("vtd %u: cmp %u -> %u\\n", u, vtd_cmp[u], old_cmp);
                            vtd_cmp[u] = old_cmp;                            
                        }
                    }
                }
            }
        }
        
        print_all(vtd_clr, vtd_cmp, (*num_cmps));
        printf("END OF ASSIGNMENT PROCESS FOR vtd %u\\n", w);
    }


    __device__ void get_TV_local(ushort *clr_pop, float pop_tot, float unif_dens, float *TV_vec, ushort c) {
        TV_vec[c] = fabs(clr_pop[c] / pop_tot - unif_dens);
    }


    __global__ void gerrymander_kernel(ushort *vtd_edge0, ushort *vtd_edge1, ushort *nbrs, ushort *nbr_rng, ushort *vtd_pop, ushort *vtd_clr_arr, uint *rand_seeds) {    
        const ushort threads_per_block = blockDim.z * blockDim.y * blockDim.x;
        const ushort blockId = blockIdx.z*(gridDim.y*gridDim.x) + blockIdx.y*(gridDim.x) + blockIdx.x;
        const ushort threadIdLoc = threadIdx.z*(blockDim.y*blockDim.x) + threadIdx.y*(blockDim.x) + threadIdx.x;
        const uint threadIdGlob = blockId*threads_per_block + threadIdLoc;
        
        curandState rand_state;
        curand_init(rand_seeds[threadIdGlob], threadIdGlob, 0, &rand_state);
        
        ushort u, v, w, c;
        
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

        __shared__ ushort vtd_pop_shr[NUM_VTDS];
        to_shared(vtd_pop, vtd_pop_shr, NUM_VTDS, threadIdLoc);
        if(threadIdLoc == 0) {
            printf("Printing A_shr = vtd_pop_shr\\n");
            print_shared(vtd_pop_shr, NUM_VTDS);
            printf("\\n");
        }
        
               
        ushort vtd_clr[NUM_VTDS];
        uint idx_start = threadIdGlob*NUM_VTDS;
        for(v = 0; v < NUM_VTDS; v++) {            
            vtd_clr[v] = vtd_clr_arr[idx_start+v];
        }
        
       

        ushort vtd_cmp[NUM_VTDS];
        ushort num_cmps = NUM_VTDS;
        for(w = 0; w < NUM_VTDS; w++) {
            vtd_cmp[w] = w;
        }

        for(w = 0; w < NUM_VTDS; w++) {
            printf("DIV1");
            get_cmps_local(nbrs_shr, nbr_rng_shr, vtd_clr, vtd_cmp, &num_cmps, w);
            printf("There are currently %u connected cmps\\nDIV1",num_cmps);
        }
        printf("All vtds assigned to a component!\\n");
        printf("DIV2DIV2");
        ushort clr_pop[NUM_CLRS];
        for(c = 0; c < NUM_CLRS; c++) {
            clr_pop[c] = 0;
        }
        
        float pop_tot = 0;
        for(w = 0; w < NUM_VTDS; w++) {
            clr_pop[vtd_clr[w]] += vtd_pop_shr[w];
            pop_tot += vtd_pop_shr[w];
        }
        
        float unif_dens = 1 / (float)NUM_CLRS;
        float TV_vec[NUM_CLRS];
        float TV = 0;
        for(c = 0; c < NUM_CLRS; c++) {
            get_TV_local(clr_pop, pop_tot, unif_dens, TV_vec, c);
            TV += TV_vec[c];
        }
        
        printf("unif_dens=%f\\n",unif_dens);
        for(c = 0; c < NUM_CLRS; c++) {
            printf("clr %u: pop=%u, dens=%f, TV_local=%f\\n", c, clr_pop[c], clr_pop[c]/pop_tot, TV_vec[c]);            
        }
        printf("TV=%f\\n",TV);
        printf("DIV2DIV2");

        ushort step = 0;
        ushort draw = 0;
        ushort const max_draws = 10;
        ushort edge = 0;
        float r;
        print_all(vtd_clr, vtd_cmp, num_cmps);
        for(step = 0; step < MAX_STEPS; step++) {
            printf("STARTING STEP %u\\n",step);
            for(draw = 0; draw < max_draws; draw++) {
                edge = (ushort)(curand_uniform(&rand_state) * NUM_EDGES);
                v = vtd_edge0[edge];
                w = vtd_edge1[edge];
                printf("I choose edge %u between vtd%u:clr%u:cmp%u and vtd%u:clr%u:cmp%u,\\n", edge, v, vtd_clr[v], vtd_cmp[v], w, vtd_clr[w], vtd_cmp[w]);
                if(vtd_clr[v] != vtd_clr[w]){
                    printf("clrs differ!! Moving on\\n");
                    break;
                }
            }
            if(draw >= max_draws) {
                printf("ERROR - could not find edge with different clrs\\n");
                break;;
            }
            r = curand_uniform(&rand_state); 
            if(r <= 0.5) {
                printf("I drew r = %f <= 0.5, so vtd%u:clr%u -> clr%u\\n\\n", r, w, vtd_clr[w], vtd_clr[v]);                
            } 
            else {
                //Swap v & w
                u = v;
                v = w;
                w = u;
                printf("I drew r = %f > 0.5, so vtd%u:clr%u -> clr%u\\n\\n", r, w, vtd_clr[w], vtd_clr[v]);                
            }
            for(c = 0; c < NUM_CLRS; c++) {
                printf("clr%u:pop%u:TV%f  ",c, clr_pop[c], TV_vec[c]);
            }
            printf("  TV=%f\\n\\n",TV);
            printf("Pop%u=%u  clr%u -> clr%u\\n\\n", w, vtd_pop[w], vtd_clr[w], vtd_clr[v]);
            clr_pop[vtd_clr[v]] += vtd_pop[w];
            clr_pop[vtd_clr[w]] -= vtd_pop[w];

            TV -= TV_vec[v];
            TV -= TV_vec[w];
            get_TV_local(clr_pop, pop_tot, unif_dens, TV_vec, vtd_clr[v]);
            get_TV_local(clr_pop, pop_tot, unif_dens, TV_vec, vtd_clr[w]);
            TV += TV_vec[v];
            TV += TV_vec[w];            
            vtd_clr[w] = vtd_clr[v];
            for(c = 0; c < NUM_CLRS; c++) {
                printf("clr%u:pop%u:TV%f  ",c, clr_pop[c], TV_vec[c]);
            }
            printf("  TV=%f\\n",TV);
            printf("DIV1");
            get_cmps_local(nbrs_shr, nbr_rng_shr, vtd_clr, vtd_cmp, &num_cmps, w);
            
            
            printf("DIV1After step %u, num components = %u and TV = %f\\n", step, num_cmps, TV);
            print_all(vtd_clr, vtd_cmp, num_cmps);
            printf("DIV2DIV2");
        }
        
        idx_start = threadIdGlob*NUM_VTDS;
        for(v = 0; v < NUM_VTDS; v++) {            
            vtd_clr_arr[idx_start+v] = vtd_clr[v];
        }

    }    
}
"""

rand_seeds = np.arange(num_sims)
rand_seeds_gpu = togpu(rand_seeds,'uint32')
vtd_edge0_gpu = togpu(vtd_edge0,'uint16')
vtd_edge1_gpu = togpu(vtd_edge1,'uint16')
nbrs_gpu = togpu(np.concatenate(nbrs),'uint16')
nbr_rng_gpu = togpu(nbr_rng,'uint16')
vtd_pop_gpu = togpu(vtd_pop,'uint16')
#vtd_clr_arr = num_clrs * np.ones([num_sims,num_vtds]).astype('uint16')
vtd_clr_arr_gpu = togpu(vtd_clr_arr,'uint16')


macros = {'DIV1':"*" * 80 + "\\n"
         ,'DIV2':"#" * 80 + "\\n"
         ,'DIV3':"!" * 80 + "\\n"
         ,'NUM_VTDS':num_vtds
         ,'NUM_NBRS':len(nbrs_gpu)
         ,'NUM_CLRS':num_clrs
         ,'NUM_EDGES':num_edges
         ,'MAX_STEPS':max_steps
         }
for key, val in macros.items():
    kernel_code = kernel_code.replace(str(key),str(val))

mod = SourceModule(kernel_code, no_extern_c=True);
gerrymander = mod.get_function("gerrymander_kernel")

gerrymander(vtd_edge0_gpu, vtd_edge1_gpu, nbrs_gpu, nbr_rng_gpu, vtd_pop_gpu, vtd_clr_arr_gpu, rand_seeds_gpu, block=blockDims, grid=gridDims, shared=0)