from helper.setup import *
from helper.pycuda_setup import *


#neighbors(u) = int vector of the neighbors of u
#neighbors = concat all neighbor(u)
#degrees = int of degress of size num_vtd
#index = int of start/end neighbors size num_vtd+1 (ie the cum sum of degrees)
#vtd0 = int vector of size num_edges
#vtd1 = int vector of size num_edges
#population = int vector of size num_vtd


num_vtds = 12
#num_vtds = 3
num_edges = 11
num_colors = 2

vtd0 = np.zeros(num_edges).astype('uint16')
vtd0_gpu = togpu(vtd0)

vtd1 = np.zeros(num_edges).astype('uint16')
vtd1_gpu = togpu(vtd1)

#color = np.random.randint(num_colors,size=num_vtds).astype('uint16')
color = np.array([0,0,1,1,0,1,1,0,1,1,1,1]).astype('uint16')
#color = np.array([0,1,0]).astype('uint16')
color_gpu = togpu(color)

component = np.arange(num_vtds).astype('uint16')
component_gpu = togpu(component)

num_components = np.uint16(num_vtds)
#num_components_gpu = togpu(num_components)

#adj = np.random.randint(2,size=(num_vtd,num_vtd)).astype('uint8')
#adj = np.ones((num_vtds,num_vtds)).astype('uint8')
#adj = np.array([[0,1,0,0,0,0,0,0,0,0,0,0],
#                [1,0,0,0,1,0,0,0,0,0,1,1],
#                [0,0,0,0,0,1,0,0,0,1,1,0],
#                [0,0,0,0,0,0,1,0,1,1,0,1],
#                [0,1,0,0,0,0,0,1,0,0,1,1],
#                [0,0,1,0,0,0,0,1,0,1,1,0],
#                [0,0,0,1,0,0,0,1,0,1,0,1],
#                [0,0,0,0,1,1,1,0,0,1,1,1],
#                [0,0,0,1,0,0,0,0,0,0,0,0],
#                [0,0,1,1,0,1,1,1,0,0,0,0], 
#                [0,1,1,0,1,1,0,1,0,0,0,0],
#                [0,1,0,1,1,0,1,1,0,0,0,0]]).astype('uint8')


#adj = np.triu(adj,k=1)
#adj = adj + adj.T
#print(np.all(adj == adj.T))
#adj_gpu = togpu(adj)

neighbors



mod = SourceModule("""
    #include <stdio.h>
    __global__ void compute_components_gpu(ushort *component, ushort num_components, ushort *color, ushort num_vtds, char *adj)
    {
        const uint blockId = blockIdx.z*(gridDim.y*gridDim.x) + blockIdx.y*(gridDim.x) + blockIdx.x;
        const uint threadIdLoc = threadIdx.z*(blockDim.y*blockDim.x) + threadIdx.y*(blockDim.x) + threadIdx.x;
        const uint threadIdGlob = blockId*(blockDim.z*blockDim.y*blockDim.x) + threadIdLoc;

        const char *div1 = "***********************************************************************************************\\n";
        const char *div2 = "###############################################################################################\\n";

        /*
        printf("%sHello from block (%d,%d,%d) and thread (%d,%d,%d)\\n%s",
                div1, blockIdx.x, blockIdx.y, blockIdx.x, threadIdx.x, threadIdx.y, threadIdx.x, div2);
        */

        
        int u, v, w;
        num_components = 0;
        //printf("%d\\n",num_vtds);
        //printf("%d\\n",num_components);
        for(w = 0; w < num_vtds; w++)
        {
            printf("\\n\\nLooking at vtd %d\\n\\n",w);
            component[w] = num_components;
            //printf("vtd w=%d is in component %d\\n", w,component[w]);
            for(v = 0; v < w; v++)
            {
                if(color[v] == color[w])
                {
                    printf("Matched color on %d and %d.\\n", v, w);
                    printf("vtd v=%d is in component %d\\n", v,component[v]);
                    if(component[v] < component[w])
                    {
                        if(adj[v*num_vtds + w] == 1)
                        {
                            printf("%d and %d were adjacent, with comp[%d] < comp[%d].\\n", v, w, v, w);
                            component[w] = component[v];
                        }
                    }
                }
            }
            if(component[w] == num_components)
            {
                //printf("Incrementing num_components to %d",num_components);
                num_components++;
            }
            
            for(u = 0; u < w; u++)
            {
                if(color[u] == color[w])
                {
                    if(component[u] > component[w])
                    {
                        if(adj[u*num_vtds + w] == 1)
                        {
                            for(v = 0; v < w; v++)
                            {
                                if(v != u)
                                {
                                    if(component[v] == component[u])
                                    {
                                        printf("Updating component of %d from %d to %d\\n", v, component[v], component[w]);
                                        component[v] = component[w];
                                    }
                                }
                            }
                            component[u] = component[w];
                        }
                    }
                }
            }
        }
        printf("MADE IT\\n");
    }
""")
compute_components = mod.get_function("compute_components_gpu")

blockDims = (1,1,1)
gridDims = (1,1,1)
print(component_gpu)
print(color_gpu)
compute_components(component_gpu, num_components, color_gpu, np.uint16(num_vtds), adj_gpu, block=blockDims, grid=gridDims, shared=0)
print(color)
print(component_gpu)
print(num_components)