def print_stats():
    print("vtd colors")
    print(clr)
    print("#vtds per color")
    print(clr_count)
    print("population per color")
    print(clr_pop)
    print("number of uncolored vtds = %d"%num_uclr)

def recolor(w, c):
    global clr_count, clr_pop, num_uclr
    clr_count[clr[w]] -= 1
    clr_pop[clr[w]] -= vtd_pop[w]
    clr[w] = c
    clr_count[clr[w]] += 1
    clr_pop[clr[w]] += vtd_pop[w]
    num_uclr -= 1

clr_arr = num_clrs * np.ones([num_sims,num_vtds]).astype('uint16')
small_deg_vtds = np.argpartition(degree,num_clrs)[:num_clrs]
for s in range(num_sims):
    #print("Map %d"%s)
    #print("Initial coloring - place colors on vtds with smallest degree")    
    clr = clr_arr[s,:]   
    
    init_clr = np.random.permutation(num_clrs)
    clr[small_deg_vtds] = init_clr    
    clr_count = np.bincount(clr).astype('uint16')
    clr_pop = np.bincount(clr,weights=vtd_pop).astype('uint16')    
    num_uclr = clr_count[-1]   
    #print_stats()

    Steps = num_uclr
    for step in range(Steps):
        clr_pop_srt = np.argsort(clr_pop[:-1])
        #print(DIV2)
        #print(DIV2)
        #print("STEP %d"%step)        
        #print_stats()
        #print(clr_pop_srt)
        #print(DIV1)
        for c in clr_pop_srt:
            c_idx = np.where(clr==c)[0]
            np.random.shuffle(c_idx)
            #print("Trying to add vtd to color %d.  Here are vtds with that color (in random order):"%c)
            #print(c_idx)
            found = False
            for v in c_idx:
                M = nbrs[v]
                N = M[clr[M] == num_clrs]
                #print("Trying to color a neighbor of %d.  Here are the uncolored neighbors."%v)
                #print(N)
                if(len(N) == 0):
                    #print("All neighbors of %d are already colored.  Trying next vtd with color %d."%(v, c))
                    pass
                else:
                    w = np.random.choice(N,1)
                    #print("vtd %d is NOT already colored.  Coloring it %d."%(w,c))
                    recolor(w, c)
                    found = True
                    break
            if(found == False):
                #print("Unable to add ANY neighboring vtds to color %d because all neighbors are already colored.  Trying next color."%c)
                #print(DIV1)
                #print(DIV1)
                pass
            if(found == True):            
                break
        if(found == False):
            c = clr_pop_srt[0]
            nbr = np.where(clr==num_clrs)[0]
            w = nbr[np.argmin(degree[nbr])]
            #print("Can't add ANY neighboring vtds to ANY color.  Putting color %d on a vtd %d - it has smallest degree of the uncolored vtds."%(c,w))
            recolor(w, c)
        #print(DIV2)
        #print_stats()
    #print(DIV2)
    #print(DIV1)
    #print(DIV2)
    #print(DIV1)
    #clr_arr[s,:] = clr
#print(clr_arr)