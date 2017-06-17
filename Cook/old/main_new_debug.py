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


    __device__ void print_vtd(ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp, ushort v) {
        printf("vtd%u:clr%u:cmp%u:pop%u", v, vtd_clr[v], vtd_cmp[v], vtd_pop[v]);
    }


    __device__ void print_clr(ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp, ushort c) {
        ushort v;
        uint pop = 0;
        printf("Clr%u - ",c);
        for(v = 0; v < NUM_VTDS; v++) {
            if(vtd_clr[v] ==c) {
                print_vtd(vtd_pop, vtd_clr, vtd_cmp, v);
                printf("  ");
                pop += vtd_pop[v];
            }
        }
        printf("|||  TOTAL POP:%u\\n",pop);
    }



    __device__ void print_clrs(ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp) {
        ushort c;
        for(c = 0; c < NUM_CLRS; c++) {
            print_clr(vtd_pop, vtd_clr, vtd_cmp, c);
        }


    __device__ void print_cmp(ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp, ushort c) {
        ushort v;
        uint pop = 0;
        printf("Cmp%u - ", c);
        for(v = 0; v < NUM_VTDS; v++) {
            if(vtd_cmp[v] ==c) {
                print_vtd(vtd_pop, vtd_clr, vtd_cmp, v);
                printf("  ");
                pop += vtd_pop[v];
            }
        }
        printf("|||  TOTAL POP:%u\\n", pop);
    }



    __device__ void print_cmps(ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp, ushort num_cmps) {
        ushort c;
        for(c = 0; c < num_cmps; c++) {
            print_cmp(vtd_pop, vtd_clr, vtd_cmp, c);
        }


    }

    
    __device__ void argsrt(ushort *val, ushort *idx, ushort N) {
        ushort i, j, a, b;
/*
        for(i = 0; i < N; i++) {
            printf("%u  ", i);
        }
        printf("\\n");
        for(i = 0; i < N; i++) {
            printf("%u  ", val[i]);
        }
        printf("\\n\\n");
*/
        for(i = 0; i < N; i++) {
            idx[i] = N;
        }
        for(j = 0; j < N; j++) {
            a = j;
            for(i = 0; i < N; i++) {
                if( (idx[i] >= N) || (val[a] < val[idx[i]]) || ( (val[a] == val[idx[i]]) && (a < idx[i]) ) ) {
                    b = idx[i];
                    idx[i] = a;
                    a = b;
                }
            }
        }
/*
        for(i = 0; i < N; i++) {
            printf("%u  ", idx[i]);
        }
        printf("\\n");
        for(i = 0; i < N; i++) {
            printf("%u  ", val[idx[i]]);
        }
        printf("\\n\\n");
*/
    }




    __device__ void reclr_vtd(ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp, ushort *clr_count, ushort *clr_pop, ushort v, ushort c) {
        print_vtd(vtd_pop, vtd_clr, vtd_cmp, v);
        if(vtd_clr[v] < NUM_CLRS) {
            clr_count[vtd_clr[v]] -= 1;
            clr_pop[vtd_clr[v]] -= vtd_pop[v];
        }
        vtd_clr[v] = c;
        clr_count[vtd_clr[v]] += 1;
        clr_pop[vtd_clr[v]] += vtd_pop[v];
        printf("->clr%u", vtd_clr[v]);
    }



    __device__ void recmp_vtd(ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp, ushort v, ushort c) {
        print_vtd(vtd_pop, vtd_clr, vtd_cmp, v);
        vtd_cmp[v] = c;
        printf("->clr%u", vtd_clr[v]);
    }



    __device__ void get_cmp(ushort *nbrs, ushort *nbr_rng, ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp, ushort *num_cmps, ushort w) {         
        printf("Trying to assign vtd%u to a connected cmp.\\n", w);        
        ushort u, v, idx, orig_cmp, old_cmp, new_cmp;
        
        orig_cmp = vtd_cmp[w];
        new_cmp = (*num_cmps);
        //vtd_cmp[w] = new_cmp;
        //print_vtd(vtd_pop, vtd_clr, vtd_cmp, w);        
        //printf("vtd%u:clr%u:cmp%u->cmp%u\\n", w, vtd_clr[w], orig_cmp, new_cmp);
        recmp_vtd(vtd_pop, vtd_clr, vtd_cmp, w, new_cmp);
        for(idx = nbr_rng[w]; idx < nbr_rng[w+1]; idx++) {
            v = nbrs[idx];
            printf("Checking neighbor ");// vtd%u:clr%u:cmp%u\\n", v, vtd_clr[v], vtd_cmp[v]);
            print_vtd(vtd_pop, vtd_clr, vtd_cmp, v);
            printf("\\n");
            if(vtd_clr[v] == vtd_clr[w]) {
                printf("Clrs match\\n");
                if(vtd_cmp[v] < new_cmp) {
                    new_cmp = vtd_cmp[v];
                    printf("Smaller cmp found!!\\n");
                    print_vtd(vtd_pop, vtd_clr, vtd_cmp, v);
                    vtd_cmp[w] = new_cmp;
                    printf("->cmp%u\\n", vtd_cmp[w]);
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
            printf("Incrementing num_cmps.\\n");
            (*num_cmps)++;
        }

        char found = 0;
        if(new_cmp != orig_cmp) {
            printf("\\nvtd%u switched cmp%u->cmp%u.  Let's check if somebody else is still in cmp%u.\\n",w, orig_cmp, new_cmp, orig_cmp);
            for(v = 0; v < NUM_VTDS; v++) {
                if(vtd_cmp[v] == orig_cmp) {
                    printf("Yes");
                    print_vtd(vtd_pop, vtd_clr, vtd_cmp, v);
                    found = 1;
                    break;
                }
            }
            if(found == 0) {
                (*num_cmps)--;                
                printf("No, cmp%u is now empty.  I will relabel every vtd in cmp%u->cmp%u\\n", orig_cmp, (*num_cmps), orig_cmp);
                for(v = 0; v < NUM_VTDS; v++) {
                    if(vtd_cmp[v] == (*num_cmps)){
                            print_vtd(vtd_pop, vtd_clr, vtd_cmp, v);
                            vtd_cmp[v] = orig_cmp;
                            printf("->%u", vtd_cmp[v]);

                    }
                }
                printf("\\n");
            }
        }

        printf("Now we check if this fused >=2 previously distinct cmps into 1 big cmp\\n");
        for(idx = nbr_rng[w]; idx < nbr_rng[w+1]; idx++) {
            v = nbrs[idx];
            print_vtd(vtd_pop, vtd_clr, vtd_cmp, v);
            printf("  ");
            if(vtd_clr[v] == vtd_clr[w]) {
                if(vtd_cmp[v] > vtd_cmp[w]) {
                    (*num_cmps)--;
                    old_cmp = vtd_cmp[v];
                    new_cmp = vtd_cmp[w];
                    
                    printf("\\nYes, cmp%u fused with cmp%u.  Relabeling all vtds in cmp%u->cmp%u and all vtds in cmp%u->cmp%u.\\n", old_cmp, new_cmp, old_cmp, new_cmp, (*num_cmps), old_cmp);
                    
                    for(u = 0; u < NUM_VTDS; u++) {
                        if(vtd_cmp[u] == old_cmp) {
                            print_vtd(vtd_pop, vtd_clr, vtd_cmp, u);
                            vtd_cmp[u] = new_cmp;
                            printf("->%u", vtd_cmp[u]);
                        } else if(vtd_cmp[u] == (*num_cmps)) {
                            print_vtd(vtd_pop, vtd_clr, vtd_cmp, u);
                            vtd_cmp[u] = old_cmp;
                            printf("->%u", vtd_cmp[u]);
                        }
                    }
                }
            }
        }
        printf("\\n");
        print_cmps(vtd_pop, vtd_clr, vtd_cmp, (*num_cmps));
        printf("END OF ASSIGNMENT PROCESS FOR vtd%u\\n", w);
    }


    __device__ void get_TV(ushort *clr_pop, float tot_pop, float unif_dens, float *TV_vec, ushort c) {
        TV_vec[c] = fabs(clr_pop[c] / tot_pop - unif_dens);
    }




    __global__ void gerrymander_kernel(ushort *vtd_edge0, ushort *vtd_edge1, ushort *nbrs_glob, ushort *nbr_rng_glob, ushort *vtd_pop_glob, ushort *vtd_clr_glob, uint *rand_seeds) {
    
        const ushort threads_per_block = blockDim.z * blockDim.y * blockDim.x;
        const ushort blockId = blockIdx.z*(gridDim.y*gridDim.x) + blockIdx.y*(gridDim.x) + blockIdx.x;
        const ushort threadIdLoc = threadIdx.z*(blockDim.y*blockDim.x) + threadIdx.y*(blockDim.x) + threadIdx.x;
        const uint threadIdGlob = blockId*threads_per_block + threadIdLoc;
 
        curandState rand_state;
        curand_init(rand_seeds[threadIdGlob], threadIdGlob, 0, &rand_state);
        
        ushort u, v, w, c, i, j, k, step, idx;
        
        printf("DIV1DIV2DIV3");
        printf("Start reading data.\\n");
        __shared__ ushort nbrs[NUM_NBRS];
        to_shared(nbrs_glob, nbrs, NUM_NBRS, threadIdLoc);
        if(threadIdLoc == 0) {
            printf("Printing A_shr = nbrs\\n");
            print_shared(nbrs, NUM_NBRS);
            printf("\\n");
        }

        __shared__ ushort nbr_rng[NUM_VTDS+1];
        to_shared(nbr_rng_glob, nbr_rng, NUM_VTDS+1, threadIdLoc);
        if(threadIdLoc == 0) {
            printf("Printing A_shr = nbr_rng\\n");
            print_shared(nbr_rng, NUM_VTDS+1);
            printf("\\n");
        }

        __shared__ ushort vtd_pop[NUM_VTDS];        
        to_shared(vtd_pop_glob, vtd_pop, NUM_VTDS, threadIdLoc);
        if(threadIdLoc == 0) {
            printf("Printing A_shr = vtd_pop\\n");
            print_shared(vtd_pop, NUM_VTDS);
            printf("\\n");
        }

        ushort vtd_degree[NUM_VTDS];
        float tot_pop = 0;
        ushort vtd_clr[NUM_VTDS];
        ushort vtd_cmp[NUM_VTDS];
        ushort num_cmps = NUM_VTDS;
        uint idx_start = threadIdGlob*NUM_VTDS;
        for(v = 0; v < NUM_VTDS; v++) {
            vtd_degree[v] = nbr_rng[v+1] - nbr_rng[v];
            tot_pop += (float)vtd_pop[v];
            vtd_clr[v] = vtd_clr_glob[idx_start+v];
            vtd_cmp[v] = v;            
        }        
        print_cmps(vtd_pop, vtd_clr, vtd_cmp, num_cmps);
        printf("Finished reading data.\\n");
        printf("DIV1DIV2DIV3");




        printf("Start initial color assignment.\\n");
        ushort vtd_degree_argsrt[NUM_VTDS];
        argsrt(vtd_degree, vtd_degree_argsrt, NUM_VTDS);

        ushort clr_list[NUM_CLRS];
        for(c = 0; c < NUM_CLRS; c++) {
            clr_list[c] = c;
        }
        ushort clr_count[NUM_CLRS];
        ushort clr_pop[NUM_CLRS];
        ushort clr_pop_argsrt[NUM_CLRS];        
        
        //Start coloring with the vtds of smallest degree to minimize chance of islands
        ushort clrs_left = NUM_CLRS;
        for(i = 0; i < NUM_CLRS; i++) {
            v = vtd_degree_argsrt[i];
            //Select randomly from unused colors
            j = (ushort)(curand_uniform(&rand_state) * clrs_left);
            c = clr_list[j];
            
            vtd_clr[v] = c;
            clr_count[c] = 1;
            clr_pop[c] = vtd_pop[v];
            
            clrs_left -= 1;
            //Keep all unused colors before used colors in clr_list            
            clr_list[j] = clr_list[clrs_left];
            clr_list[clrs_left] = c;
            
            printf("DIV2");
            printf("DIV2");
            printf("STEP %u\\n",i);
            printf("Assigning clr%u to vtd%u with degree%u.\\n", c, v, vtd_degree[v]);
            print_clrs(vtd_pop, vtd_clr, vtd_cmp, num_cmps);
        }

        ushort deg;
        ushort skip;
        ushort v_offset, w_offset;
        char found = 0;
        for(step = NUM_CLRS; step < NUM_VTDS; step++) {
            printf("DIV2");
            printf("DIV2");
            printf("STEP %u\\n",step);

            found = 0;
            argsrt(clr_pop, clr_pop_argsrt, NUM_CLRS);
            for(i = 0; i < NUM_CLRS; i++) {
                c = clr_pop_argsrt[i];
                printf("DIV2");
                printf("Looking for an uncolored neighbor of ");
                
                
                print_vtd(vtd_pop, vtd_clr, vtd_cmp, v);
                //, vtd_cmp
                
                //clr%u:pop%u.  The vtds with clr%u are:.\\n", c, clr_pop[c], c);
                for(v = 0; v < NUM_VTDS; v++) {
                    if(vtd_clr[v] == c) {
                        printf("%u  ", v);
                    }
                }
                printf("\\n");
                //Randomize by starting with the "skip"^th vtd with clr c
                //The "skip"^th vtd with color c is in slot "offset"
                skip = (ushort)(curand_uniform(&rand_state) * clr_count[c]);
                v_offset = 0;
                j = 0;
                while(j <= skip) {
                    if(vtd_clr[v_offset] == c) {
                        j++;
                    }
                    v_offset++;
                }
                v_offset--;

                for(j = 0; j < NUM_VTDS; j++) {
                    v = (v_offset + j) % NUM_VTDS;
                    if(vtd_clr[v] == c) {
                        printf("DIV1");
                        printf("Looking for uncolored neighbor of vtd%u.  The neighbors are:.\\n", v);
                        for(k = nbr_rng[v]; k < nbr_rng[v+1]; k++) {
                            printf("%u ",nbrs[k]);
                        }
                        printf("\\n");
                        deg = vtd_degree[v];
                        //Randomize by starting with the "skip"^th neighbor of v
                        w_offset = (ushort)(curand_uniform(&rand_state) * deg);
                        for(k = 0; k < deg; k++) {
                            idx = nbr_rng[v] + (w_offset + k) % deg;
                            w = nbrs[idx];
                            printf("Trying neighbor vtd%u.\\n",w);
                            if(vtd_clr[w] == NUM_CLRS) {
                                printf("HURRAY, it is not colored.  Assigning clr%u to vtd%u.\\n", c, w);
                                reclr_vtd(vtd_pop, vtd_clr, vtd_cmp, clr_count, clr_pop, w, c);
                                found = 1;
                                break;                            
                            }
                            else {
                                printf("BOO, vtd%u already has clr%u.\\n", w, vtd_clr[w]);
                            }                            
                        }
                        if(found == 0) {
                            printf("BOO, all neighbors of vtd%u are already colored.\\n", v);
                        }
                        else {
                            break;
                        }
                    }
                }
                if(found == 0) {
                    printf("BOO, all neighbors of CLR%u are already colored.\\n", c);            
                }
                else {
                    break;
                }
            }
            if(found == 0) {
                c = clr_pop_argsrt[0];
                printf("DIV3");
                printf("DIV3");
                printf("DIV3");
                printf("BOO, all neighbors of ALL CLRS are already colored.\\n");
                for(j = 0; j < NUM_VTDS; j++) {
                    v = vtd_degree_argsrt[j];
                    if(vtd_clr[v] == NUM_CLRS) {
                        printf("Assigning clr%u to vtd%u with degree%u.\\n", c, v, vtd_degree[v]);
                        reclr_vtd(vtd_pop, vtd_clr, vtd_cmp, clr_count, clr_pop, v, c);
                        found = 1;
                        break;
                    }
                }
            }
            printf("DIV1");
            print_clrs(vtd_pop, vtd_clr, vtd_cmp, num_cmps);
        }
        printf("Finished initial color assignment.\\n");
        printf("DIV1DIV2DIV3");



        printf("Start initial component detection.\\n");
        print_cmps(vtd_pop, vtd_clr, vtd_cmp, num_cmps);
        for(v = 0; v < NUM_VTDS; v++) {
            printf("DIV1");
            get_cmp(nbrs, nbr_rng, vtd_pop, vtd_clr, vtd_cmp, &num_cmps, v);
            printf("There are currently %u connected cmps\\nDIV1",num_cmps);
        }
        printf("Finish initial component detection.\\n");
        printf("DIV1DIV2DIV3");
        


        printf("Start initial TV calcuation.\\n.");
        float unif_dens = 1 / (float)NUM_CLRS;
        float TV_vec[NUM_CLRS];
        float TV = 0;
        for(c = 0; c < NUM_CLRS; c++) {
            get_TV(clr_pop, tot_pop, unif_dens, TV_vec, c);
            TV += TV_vec[c];
        }
        
        printf("tot_pop=%f,  unif_dens=%f\\n", tot_pop, unif_dens);
        for(c = 0; c < NUM_CLRS; c++) {
            printf("clr %u: pop=%u, dens=%f, TV_local=%f\\n", c, clr_pop[c], clr_pop[c]/tot_pop, TV_vec[c]);            
        }
        printf("TV=%f\\n",TV);
        printf("Finish initial TV calcuation.\\n");        
        printf("DIV1DIV2DIV3");
        
        
        
        
// /*
        printf("Start Markov chain evolution.\\n");
        ushort draw = 0;
        ushort const max_draws = 10;
        ushort edge = 0;
        float r;
        print_cmps(vtd_pop, vtd_clr, vtd_cmp, num_cmps);
        printf("DIV2");
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
                break;
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
            get_TV(clr_pop, tot_pop, unif_dens, TV_vec, vtd_clr[v]);
            get_TV(clr_pop, tot_pop, unif_dens, TV_vec, vtd_clr[w]);
            TV += TV_vec[v];
            TV += TV_vec[w];            
            vtd_clr[w] = vtd_clr[v];
            for(c = 0; c < NUM_CLRS; c++) {
                printf("clr%u:pop%u:TV%f  ",c, clr_pop[c], TV_vec[c]);
            }
            printf("  TV=%f\\n",TV);
            printf("DIV1");
            get_cmp(nbrs, nbr_rng, vtd_pop, vtd_clr, vtd_cmp, &num_cmps, v);
            
            
            printf("DIV1After step %u, num components = %u and TV = %f\\n", step, num_cmps, TV);
            print_cmps(vtd_pop, vtd_clr, vtd_cmp, num_cmps);
            printf("DIV2");
        }
        printf("Finish Markov chain evolution.\\n");
        printf("DIV1DIV2DIV3");
        

        printf("Writing data to global memory for CPU.\\n");
        idx_start = threadIdGlob*NUM_VTDS;
        for(v = 0; v < NUM_VTDS; v++) {            
            vtd_clr_glob[idx_start+v] = vtd_clr[v];
        }
        printf("I used newer\\n\\n");
// */
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
vtd_clr_arr_gpu = togpu(vtd_clr_arr,'uint16')


DIV1 = "*" * 80 + "\\n"
DIV2 = "#" * 80 + "\\n"

macros = {'DIV1': DIV1
         ,'DIV2': DIV2
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