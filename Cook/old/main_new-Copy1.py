from helper.setup import *
from helper.pycuda_setup import *



kernel_code = """
    #include <stdio.h>
    #include <math.h>
    #include <curand.h>
    #include <curand_kernel.h>

extern "C"
{
    __device__ void goodness(ushort num_cmps, float TV, float *good) {
        float good_cmp = (float)NUM_CLRS / num_cmps;
        float TV_max = 2*(1 - 1/(float)NUM_CLRS);
        float good_TV = 1 - (TV / TV_max);
        float coef_cmp = 0.7;
        float coef_TV = 1-coef_cmp;
        (*good) = (coef_cmp*good_cmp) + (coef_TV*good_TV);
        //printf("num_cmps=%u  good_cmp=%f;  TV=%f  TV_max=%f  good_TV=%f;  good=%f\\n", num_cmps, good_cmp, TV, TV_max, good_TV, (*good));
    }


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


    __device__ void print_vtd(ushort w, ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp) {
        printf("vtd%u:clr%u:cmp%u:pop%u", w, vtd_clr[w], vtd_cmp[w], vtd_pop[w]);
    }


    __device__ void print_clr(ushort c, ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp) {
        ushort w;
        ushort count = 0;
        uint pop = 0;
        //printf("Clr%u - ",c);
        for(w = 0; w < NUM_VTDS; w++) {
            if(vtd_clr[w] ==c) {
                //print_vtd(w, vtd_pop, vtd_clr, vtd_cmp);
                //printf("  ");
                count += 1;
                pop += vtd_pop[w];
            }
        }
        //printf("|||  VTD COUNT:%u  |||  TOTAL POP:%u\\n", count, pop);
    }
    
    
    __device__ void print_clrs(ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp) {
        ushort c;
        for(c = 0; c < NUM_CLRS; c++) {
            print_clr(c, vtd_pop, vtd_clr, vtd_cmp);
        }
    }


    __device__ void get_TV(ushort c, ushort *clr_pop, float tot_pop, float unif_dens, float *TV_vec) {
        TV_vec[c] = fabs(clr_pop[c] / tot_pop - unif_dens);
    }


    __device__ void reclr_vtd(ushort w, ushort c_new, ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp, ushort *clr_count, ushort *clr_pop, float tot_pop, float unif_dens, float *TV_vec, float *TV) {
        ushort c_old = vtd_clr[w];
        print_vtd(w, vtd_pop, vtd_clr, vtd_cmp);
        if(c_old < NUM_CLRS) {
            clr_count[c_old] -= 1;
            clr_pop[c_old] -= vtd_pop[w];
            (*TV) -= TV_vec[c_old];
        }
        (*TV) -= TV_vec[c_new];
        vtd_clr[w] = c_new;
        clr_count[c_new] += 1;
        clr_pop[c_new] += vtd_pop[w];
        get_TV(c_old, clr_pop, tot_pop, unif_dens, TV_vec);
        get_TV(c_new, clr_pop, tot_pop, unif_dens, TV_vec);
        (*TV) += TV_vec[c_old];
        (*TV) += TV_vec[c_new];
        
        //printf("->clr%u", vtd_clr[w]);
    }



    __device__ void print_cmp(ushort c, ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp) {
        ushort w;
        ushort count = 0;
        uint pop = 0;
        //printf("Cmp%u - ", c);
        for(w = 0; w < NUM_VTDS; w++) {
            if(vtd_cmp[w] ==c) {
                print_vtd(w, vtd_pop, vtd_clr, vtd_cmp);
                //printf("  ");
                count += 1;
                pop += vtd_pop[w];
            }
        }
        //printf("|||  VTD COUNT:%u  |||  TOTAL POP:%u\\n", count, pop);
    }





    __device__ void print_cmps(ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp, ushort num_cmps) {
        ushort c;
        for(c = 0; c < num_cmps; c++) {
            print_cmp(c, vtd_pop, vtd_clr, vtd_cmp);
        }
    }



    __device__ void recmp_vtd(ushort w, ushort c, ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp) {
        print_vtd(w, vtd_pop, vtd_clr, vtd_cmp);
        vtd_cmp[w] = c;
        //printf("->cmp%u", vtd_cmp[w]);
    }

    __device__ void print_TV(ushort c, ushort *clr_count, ushort *clr_pop, float tot_pop, float unif_dens, float *TV_vec) {
        printf("clr%u:#vtds%u:pop%u:dens%f:unif:%f:TV%f", c, clr_count[c], clr_pop[c], clr_pop[c]/tot_pop, unif_dens, TV_vec[c]);
        printf("  ");
    }
    

    __device__ void print_TVs(ushort *clr_count, ushort *clr_pop, float tot_pop, float unif_dens, float *TV_vec) {
        ushort c;
        float TV = 0;
        for(c = 0; c < NUM_CLRS; c++) {
            print_TV(c, clr_count, clr_pop, tot_pop, unif_dens, TV_vec);
            printf("\\n");
            TV += TV_vec[c];
        }
        printf("TV=%f\\n",TV);
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




    __device__ void get_cmp(ushort w, ushort *nbrs, ushort *nbr_rng, ushort *vtd_pop, ushort *vtd_clr, ushort *vtd_cmp, ushort *num_cmps) {
        ushort u, v, idx, orig_cmp;
        //printf("START OF CMP ASSIGNMENT PROCESS FOR VTD%u\\n", w);
        //print_vtd(w, vtd_pop, vtd_clr, vtd_cmp);
        //printf(" to a connected component.\\n");
        orig_cmp = vtd_cmp[w];
        recmp_vtd(w, (*num_cmps), vtd_pop, vtd_clr, vtd_cmp);
        //printf("\\n");
        for(idx = nbr_rng[w]; idx < nbr_rng[w+1]; idx++) {
            v = nbrs[idx];
            //printf("Checking neighbor  ");
            //print_vtd(v, vtd_pop, vtd_clr, vtd_cmp);            
            if(vtd_clr[v] != vtd_clr[w]) {
                //printf("  Clrs do not match.");
            }
            else {
                //printf("  Clrs match.");
                if(vtd_cmp[v] >= vtd_cmp[w]) {
                    //printf("  Larger cmp.  Not reassigning.");
                }
                else {
                    //printf("  Smaller cmp.  Reassigning.  ");
                    recmp_vtd(w, vtd_cmp[v], vtd_pop, vtd_clr, vtd_cmp);                    
                }
            }
            //printf("\\n");            
            if(vtd_cmp[w] == 0) {
                //printf("vtd%u assigned to cmp0.  Stop.\\n", w);
                break;
            }
        }        
        //print_vtd(w, vtd_pop, vtd_clr, vtd_cmp);
        
        if(vtd_cmp[w] == (*num_cmps)) {
            //printf("\\nIncrementing num_cmps.");
            (*num_cmps)++;
        }
        //printf("\\nDIV1");
        
        
        
        //printf("Now we check if this disconnected cmp%u.\\n", orig_cmp);
        char count = 0;
        for(idx = nbr_rng[w]; idx < nbr_rng[w+1]; idx++) {
            v = nbrs[idx];
            //printf("Checking neighbor  ");
            //print_vtd(v, vtd_pop, vtd_clr, vtd_cmp);
            if(vtd_clr[v] == vtd_clr[w]) {
                //printf("  Clrs match.  No potential problem.\\n");
            }
            else {
                //printf("  Clrs do not match.  ");
                if(vtd_cmp[v] != orig_cmp) {
                    //printf("  But vtd%u is NOT in cmp%u.  No potential problem.\\n", v, orig_cmp);
                }
                else {
                    count += 1;
                    //printf("  And vtd%u IS in cmp%u.  ", v, orig_cmp);
                    if(count <= 1) {
                        //printf("  If we find another, there could be a problem.\\n");
                    }
                    else {
                        //printf("  That's >=2 nbrs of vtd%u in cmp%u.  So this step might have split cmp%u into multiple pieces.  We must check.\\n", w, orig_cmp, orig_cmp);
                        //printf("DIV3DIV3DIV3");
                        break;
                    }
                }
            }
        }
        if(count > 1) {
            ushort cutoff = (*num_cmps)-1;
            //printf("First, I will reassign every vtd in currently cmp%u to its own unique cmp", orig_cmp);
            if(orig_cmp < cutoff) {
                //printf(" and every vtd in cmp%u to cmp%u", cutoff, orig_cmp);
            }
            //printf(".\\n");
            (*num_cmps)--;
            for(v = 0; v < NUM_VTDS; v++) {
                if(vtd_cmp[v] == orig_cmp) {
                    (*num_cmps)++;
                    recmp_vtd(v, (*num_cmps)-1, vtd_pop, vtd_clr, vtd_cmp);
                    //printf("  ");
                    
                }
                else if(vtd_cmp[v] == cutoff) {
                    recmp_vtd(v, orig_cmp, vtd_pop, vtd_clr, vtd_cmp);
                    //printf("  ");
                }
            }
            //printf("\\nNow, I will run get_cmp on each of those vtd's separately.\\n");
            printf("DIV1DIV1");
            for(v = 0; v < NUM_VTDS; v++) {                
                if(vtd_cmp[v] >= cutoff) {
                    get_cmp(v, nbrs, nbr_rng, vtd_pop, vtd_clr, vtd_cmp, num_cmps);
                    printf("DIV1DIV1");
                }
            }
        }
        printf("DIV1");
        
        
        

        char found = 0;
        if(vtd_cmp[w] == orig_cmp) {
            printf("vtd%u did not change cmp.\\n", w);
        }
        else {
            printf("vtd%u switched cmp%u->cmp%u.  Let's check if somebody else is still in cmp%u.\\n",w, orig_cmp, vtd_cmp[w], orig_cmp);
            for(v = 0; v < NUM_VTDS; v++) {
                if(vtd_cmp[v] == orig_cmp) {
                    printf("Yes  ");
                    print_vtd(v, vtd_pop, vtd_clr, vtd_cmp);
                    found = 1;
                    break;
                }
            }
            if(found == 0) {
                (*num_cmps)--;
                printf("No, cmp%u is now empty.  I will relabel every vtd in cmp%u->cmp%u.\\n", orig_cmp, (*num_cmps), orig_cmp);
                for(v = 0; v < NUM_VTDS; v++) {
                    if(vtd_cmp[v] == (*num_cmps)){
                            printf("   ");
                            recmp_vtd(v, orig_cmp, vtd_pop, vtd_clr, vtd_cmp);
                    }
                }
            }
            printf("\\nDIV1");
            
            
            


            printf("Now we check if this fused >=2 previously distinct cmps into 1 big cmp.\\n");
            for(idx = nbr_rng[w]; idx < nbr_rng[w+1]; idx++) {
                v = nbrs[idx];
                printf("Checking neighbor.  ");
                print_vtd(v, vtd_pop, vtd_clr, vtd_cmp);

                if(vtd_clr[v] != vtd_clr[w]) {
                    printf("  Clrs do not match.");
                }
                else {
                    printf("  Clrs match.  ");
                    if(vtd_cmp[v] <= vtd_cmp[w]) {
                        printf("  Smaller (or equal) cmp.  Not reassigning.");
                    }
                    else {
                        (*num_cmps)--;
                        orig_cmp = vtd_cmp[v];

                        printf("  Larger cmp.  So cmp%u fused with cmp%u.  Reassigning all vtds in cmp%u->cmp%u", orig_cmp, vtd_cmp[w], orig_cmp, vtd_cmp[w]);
                        if(orig_cmp < (*num_cmps)) {
                            printf(" and all vtds in cmp%u->cmp%u", (*num_cmps), orig_cmp);
                        }
                        printf(".\\n");
                        for(u = 0; u < NUM_VTDS; u++) {
                            if(vtd_cmp[u] == orig_cmp) {
                                printf("   ");
                                recmp_vtd(u, vtd_cmp[w], vtd_pop, vtd_clr, vtd_cmp);
                            } else if(vtd_cmp[u] == (*num_cmps)) {
                                printf("   ");
                                recmp_vtd(u, orig_cmp, vtd_pop, vtd_clr, vtd_cmp);
                            }
                        }
                    }
                }
                printf("\\n");
            }
        }
        printf("DIV1");
        print_cmps(vtd_pop, vtd_clr, vtd_cmp, (*num_cmps));
        printf("END OF CMP ASSIGNMENT PROCESS FOR VTD%u\\n", w);
    }






    __global__ void gerrymander_kernel(ushort *vtd_edge0, ushort *vtd_edge1, ushort *nbrs_glob, ushort *nbr_rng_glob, ushort *vtd_pop_glob, ushort *vtd_clr_glob, uint *rand_seeds, float *good_glob, float anneal_rate) {
    
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
/*
        if(threadIdLoc == 0) {
            printf("Printing A_shr = nbrs\\n");
            print_shared(nbrs, NUM_NBRS);
            printf("\\n");
        }
*/

        __shared__ ushort nbr_rng[NUM_VTDS+1];
        to_shared(nbr_rng_glob, nbr_rng, NUM_VTDS+1, threadIdLoc);
/*
        if(threadIdLoc == 0) {
           printf("Printing A_shr = nbr_rng\\n");
            print_shared(nbr_rng, NUM_VTDS+1);
            printf("\\n");
        }
*/
        __shared__ ushort vtd_pop[NUM_VTDS];        
        to_shared(vtd_pop_glob, vtd_pop, NUM_VTDS, threadIdLoc);
/*
        if(threadIdLoc == 0) {
            printf("Printing A_shr = vtd_pop\\n");
            print_shared(vtd_pop, NUM_VTDS);
            printf("\\n");
        }
*/
        ushort vtd_degree[NUM_VTDS];
        float tot_pop = 0;
        ushort vtd_clr[NUM_VTDS];
        ushort vtd_cmp[NUM_VTDS];
        ushort num_cmps = NUM_VTDS;
        uint idx_start = threadIdGlob*NUM_VTDS;
        for(w = 0; w < NUM_VTDS; w++) {
            vtd_degree[w] = nbr_rng[w+1] - nbr_rng[w];
            tot_pop += (float)vtd_pop[w];
            vtd_clr[w] = vtd_clr_glob[idx_start+w];
            vtd_cmp[w] = w;            
        }        
        print_cmps(vtd_pop, vtd_clr, vtd_cmp, num_cmps);
        printf("Finished reading data.\\n");
        printf("DIV1DIV2DIV3DIV1DIV2DIV3");




        printf("Start initial color assignment.\\n");
        printf("DIV2DIV2");
        ushort vtd_degree_argsrt[NUM_VTDS];
        argsrt(vtd_degree, vtd_degree_argsrt, NUM_VTDS);

        ushort clr_list[NUM_CLRS];
        ushort clr_count[NUM_CLRS];
        ushort clr_pop[NUM_CLRS];
        ushort clr_pop_argsrt[NUM_CLRS];
        float unif_dens = 1 / (float)NUM_CLRS;
        float TV_vec[NUM_CLRS];
        float TV = 0;

        for(c = 0; c < NUM_CLRS; c++) {
            clr_list[c] = c;
            clr_count[c] = 0;
            clr_pop[c] = 0;
            clr_pop_argsrt[c] = c;
            TV_vec[c] = 0;
        }


        //Start coloring with the vtds of smallest degree to minimize chance of islands
        ushort clrs_left = NUM_CLRS;
        for(step = 0; step < NUM_CLRS; step++) {
            printf("STEP %u\\n", step);
            printf("DIV2");
            w = vtd_degree_argsrt[step];
            //Select randomly from unused colors
            j = (ushort)(curand_uniform(&rand_state) * clrs_left);
            c = clr_list[j];
            
            printf("Assigning vtd%u with degree%u to clr%u.\\n", w, vtd_degree[w], c);
            reclr_vtd(w, c, vtd_pop, vtd_clr, vtd_cmp, clr_count, clr_pop, tot_pop, unif_dens, TV_vec, &TV);
            printf("\\n");

            clrs_left -= 1;
            //Keep all unused colors before used colors in clr_list
            clr_list[j] = clr_list[clrs_left];
            clr_list[clrs_left] = c;

            print_clrs(vtd_pop, vtd_clr, vtd_cmp);
            printf("DIV2DIV2");
        }

        ushort deg;
        ushort skip;
        ushort v_offset, w_offset;
        char found = 0;
        for(step = NUM_CLRS; step < NUM_VTDS; step++) {
            printf("STEP %u\\n", step);
            found = 0;
            argsrt(clr_pop, clr_pop_argsrt, NUM_CLRS);
            for(i = 0; i < NUM_CLRS; i++) {
                c = clr_pop_argsrt[i];
                printf("DIV2");
                printf("Looking for an uncolored neighbor of ");
                print_clr(c, vtd_pop, vtd_clr, vtd_cmp);

                //Randomize by starting with the "skip"^th vtd with clr c
                //The "skip"^th vtd with color c is in slot "offset"
                skip = (ushort)(curand_uniform(&rand_state) * clr_count[c]);
                printf("rand:%u   max:%u\\n", skip, clr_count[c]-1);
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
                        printf("Looking for uncolored neighbor of ");
                        print_vtd(v, vtd_pop, vtd_clr, vtd_cmp);
                        printf(".  Neighbors are:\\n");
                        for(k = nbr_rng[v]; k < nbr_rng[v+1]; k++) {
                            w = nbrs[k];
                            print_vtd(w, vtd_pop, vtd_clr, vtd_cmp);
                            printf("  ");
                        }
                        printf("\\n");
                        deg = vtd_degree[v];
                        //Randomize by starting with the "skip"^th neighbor of v
                        w_offset = (ushort)(curand_uniform(&rand_state) * deg);
                        printf("rand:%u   max:%u\\n", w_offset, deg-1);
                        for(k = 0; k < deg; k++) {
                            idx = nbr_rng[v] + (w_offset + k) % deg;
                            w = nbrs[idx];
                            printf("Checking neighbor  ");
                            print_vtd(w, vtd_pop, vtd_clr, vtd_cmp);
                            if(vtd_clr[w] == NUM_CLRS) {
                                printf("  HURRAY, it is not colored.  Assigning vtd%u to clr%u.\\n", w, c);
                                reclr_vtd(w, c, vtd_pop, vtd_clr, vtd_cmp, clr_count, clr_pop, tot_pop, unif_dens, TV_vec, &TV);

                                printf("\\n");
                                found = 1;
                                break;
                            }
                            else {
                                printf("  BOO, it is already colored.\\n");
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
                    w = vtd_degree_argsrt[j];
                    if(vtd_clr[w] == NUM_CLRS) {
                        printf("Assigning vtd%u with degree%u to clr%u.\\n", w, vtd_degree[w], c);
                        reclr_vtd(w, c, vtd_pop, vtd_clr, vtd_cmp, clr_count, clr_pop, tot_pop, unif_dens, TV_vec, &TV);

                        printf("\\n");
                        found = 1;
                        break;
                    }
                }
            }
            printf("DIV1");
            print_clrs(vtd_pop, vtd_clr, vtd_cmp);
            printf("DIV2DIV2");
        }
        printf("Finished initial color assignment.\\n");
        printf("DIV1DIV2DIV3DIV1DIV2DIV3");



        printf("Start initial component detection.\\n");
        print_cmps(vtd_pop, vtd_clr, vtd_cmp, num_cmps);
        printf("DIV2DIV2");
        for(w = 0; w < NUM_VTDS; w++) {            
            get_cmp(w, nbrs, nbr_rng, vtd_pop, vtd_clr, vtd_cmp, &num_cmps);
            printf("DIV2DIV2");
        }
        printf("Finish initial component detection.\\n");
        printf("DIV1DIV2DIV3DIV1DIV2DIV3");      


        printf("Start initial TV calculation.\\n");
        TV = 0;
        for(c = 0; c < NUM_CLRS; c++) {
            get_TV(c, clr_pop, tot_pop, unif_dens, TV_vec);
            TV += TV_vec[c];
        }
        
        print_TVs(clr_count, clr_pop, tot_pop, unif_dens, TV_vec);
        printf("Finish initial TV calculation.\\n");        
        printf("DIV1DIV2DIV3DIV1DIV2DIV3");
        
        
        printf("START MARKOV CHAIN EVOLUTION.\\n");
        ushort draw = 0;
        ushort const max_draws = 100;
        ushort edge = 0;
        float r;
        char reject;
        float accept_ratio;
        
        ushort c_old;
        float TV_old = 2;
        float TV_best = 2;
        float good_old = 0;
        float good = 0;
        float good_diff = 0;
        float good_best = 0;
        ushort vtd_clr_best[NUM_VTDS];
        
        print_cmps(vtd_pop, vtd_clr, vtd_cmp, num_cmps);
        printf("DIV2DIV2");
        for(step = 0; step < MAX_STEPS; step++) {
            printf("STARTING STEP %u\\n",step);
            for(draw = 0; draw < max_draws; draw++) {
                edge = (ushort)(curand_uniform(&rand_state) * NUM_EDGES);
                v = vtd_edge0[edge];
                w = vtd_edge1[edge];
                printf("I choose edge %u between ", edge);
                print_vtd(v, vtd_pop, vtd_clr, vtd_cmp);
                printf(" and ");
                print_vtd(w, vtd_pop, vtd_clr, vtd_cmp);
                if(vtd_clr[v] == vtd_clr[w]){
                    printf(".  Same clr.  Try again.\\n");
                }
                else {                
                    printf(".  Different clrs!!\\n");
                    break;
                }                
            }
            if(draw >= max_draws) {
                printf("ERROR - could not find edge with different clrs\\n");
                break;
            }
            r = curand_uniform(&rand_state);            
            if(r <= 0.5) {
                printf("I drew r = %f <= 0.5  ", r);
            } 
            else {
                printf("I drew r = %f > 0.5  ", r);
                //Swap v & w
                u = v;
                v = w;
                w = u;
            }
            c_old = vtd_clr[w];
            TV_old = TV;
            good_old = good;
            
            reclr_vtd(w, vtd_clr[v], vtd_pop, vtd_clr, vtd_cmp, clr_count, clr_pop, tot_pop, unif_dens, TV_vec, &TV);
            printf("\\nDIV2");
            get_cmp(w, nbrs, nbr_rng, vtd_pop, vtd_clr, vtd_cmp, &num_cmps);
            printf("DIV2");
            print_TVs(clr_count, clr_pop, tot_pop, unif_dens, TV_vec);

            reject = 0;
            for(c = 0; c < NUM_CLRS; c++) {
                if(clr_count[c] <= 0) {
                    printf("Clr%u has no vtds.  Invalid coloring.\\n", c);
                    reject = 1;
                    break;
                }
            }

            if(reject == 0) {
                goodness(num_cmps, TV, &good);
                good_diff = good - good_old;
                if(good_diff >= 0) {
                    printf("good_new:%f >= good_old:%f.  Accept.\\n", good, good_old);
                    reject = 0;
                    if(good > good_best) {
                        printf("This is the best coloring so far - good:%f > good_best:%f!!\\n", good, good_best);
                        TV_best = TV;
                        good_best = good;
                        for(w = 0; w < NUM_VTDS; w++) {
                            vtd_clr_best[w] = vtd_clr[w];
                        }
                    }
                }
                else {                    
                    accept_ratio = expf(anneal_rate * good_diff / good_old);
                    r = curand_uniform(&rand_state);
                    //printf("good_new:%f < good_old:%f  ", good, good_old);
                    if(r < accept_ratio) {
                        reject = 0;
                        //printf("  accepting anyway because r:%f < a:%f\\n", r, accept_ratio);
                    }
                    else {
                        reject = 1;
                        //printf("  rejecting because r:%f >= a:%f\\n", r, accept_ratio);
                    }
                }
            }
                    
                    
            if(reject == 1) {
                //printf("Rejecting\\n");
                reclr_vtd(w, c_old, vtd_pop, vtd_clr, vtd_cmp, clr_count, clr_pop, tot_pop, unif_dens, TV_vec, &TV);
                //printf("\\nDIV2");
                get_cmp(w, nbrs, nbr_rng, vtd_pop, vtd_clr, vtd_cmp, &num_cmps);
                //printf("DIV2");
                print_TVs(clr_count, clr_pop, tot_pop, unif_dens, TV_vec);
                goodness(num_cmps, TV, &good);
                //printf("TV_old = %f,  TV = %f, good_old = %f, good = %f\\n", TV_old, TV, good_old, good);
            }
                
                    
                    
            
                
            
            //printf("END OF STEP%u\\n", step);
            //printf("DIV2DIV2");
             
        }
        //printf("FINISHED MARKOV CHAIN EVOLUTION\\n");
        //printf("DIV1DIV2DIV3DIV1DIV2DIV3");
        

        //printf("Writing data to global memory for CPU.  DONE.\\n");
        idx_start = threadIdGlob*NUM_VTDS;
        for(v = 0; v < NUM_VTDS; v++) {            
            vtd_clr_glob[idx_start+v] = vtd_clr_best[v];
        }
        good_glob[threadIdGlob] = good_best;
    }    
}
"""

good = np.zeros(num_sims)
rand_seeds = np.arange(num_sims)

rand_seeds_gpu = togpu(rand_seeds,'uint32')
vtd_edge0_gpu = togpu(vtd_edge0,'uint16')
vtd_edge1_gpu = togpu(vtd_edge1,'uint16')
nbrs_gpu = togpu(np.concatenate(nbrs),'uint16')
nbr_rng_gpu = togpu(nbr_rng,'uint16')
vtd_pop_gpu = togpu(vtd_pop,'uint16')
vtd_clr_arr_gpu = togpu(vtd_clr_arr,'uint16')
good_gpu = togpu(good,'float32')

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
         #,'printf':'//printf'
         }
for key, val in macros.items():
    kernel_code = kernel_code.replace(str(key),str(val))

mod = SourceModule(kernel_code, no_extern_c=True);
gerrymander = mod.get_function("gerrymander_kernel")


gerrymander(vtd_edge0_gpu, vtd_edge1_gpu, nbrs_gpu, nbr_rng_gpu, vtd_pop_gpu, vtd_clr_arr_gpu, rand_seeds_gpu, good_gpu, np.float32(anneal_rate), block=blockDims, grid=gridDims, shared=0)


    


#vtd_clr_arr = vtd_clr_arr_gpu.get()
#good = good_gpu.get()

#print(vtd_clr_arr_gpu)
#print(good_gpu)