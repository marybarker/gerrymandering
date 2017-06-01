from helper.setup import *
from helper.pycuda_setup import *

num_vtds = 8
num_edges = 20
num_colors = 3

vtd0 = np.zeros(num_edges).astype('uint16')
vtd0_gpu = togpu(vtd0)

vtd1 = np.zeros(num_edges).astype('uint16')
vtd1_gpu = togpu(vtd1)

color = np.random.randint(num_colors,size=num_vtds).astype('uint16')
color_gpu = togpu(color)

component = np.arange(num_vtds).astype('uint16')
component_gpu = togpu(component)

num_components = np.uint16(num_vtds)
#num_components_gpu = togpu(num_components)

#adj = np.random.randint(2,size=(num_vtd,num_vtd)).astype('uint8')
adj = np.ones((num_vtds,num_vtds)).astype('uint8')
adj = np.triu(adj,k=1)
adj = adj + adj.T
adj_gpu = togpu(adj)



mod = SourceModule("""
    #include <stdio.h>
    __global__ void compute_components_gpu(ushort *component, ushort num_components, ushort *color, ushort num_vtd, char *adj)
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

        
        num_components = 0;
        for(v = 0; v < num_vtd; v++){
            component[v] = num_components;
            num_components++;
            for(u = 0; u < v; u++){
                if(color[u] == color[v]){
                    if(adj[u*num_vtd + v] == 1){
                        component[v] = component[u];
                        num_components--;
                    }
                }
            }
            
            for(w = 0; w < num_vtd; w++){
                if(v != w){
                    if(color[v] == color[w]){
                        if(component[v] != component[w]){
                            if(adj[v*num_vtd + w] == 1){
                                component[w] = component[v];
                
            
        
            ##if(component[v] == v){
                component[v] = num_components;
                num_components++;
                for(w = v; w < num_vtd; w++){
                    if(component[w] == w){
                        if(color[v] == color[w]){
                            if(adj[v*num_vtd + w] == 1){
                                component[w] = component[v];






        int u, v, w;
        num_components = 0;
        for(w = 0; w < num_vtd; w++){
            component[w] = num_components;
            for(v = 0; v < w; v++){
                if(color[v] == color[w]){
                    if(component[v] < component[w])
                        if(adj[v*num_vtd + w] == 1){
                            component[w] = component[v];
                        }
                    }
                }
            }
            if(component[w] == num_components){
                num_components++;
            
            for(u = 0; u < w; u++){
                if(color[u] == color[w]){
                    if(component[u] > component[w]){
                        if(adj[u*num_vtd + w] == 1){
                            old_component = component[u];
                            new_component = component[w];
                            for(v = 0; v < w; v++);
                                if(component[v] == old_component);
                                    component[v] = new_component;
                                }

            
            




                    if(adj[v*num_vtd + w] == 1){
                        if(
                        num_components--;
                        component[w] = component[v];
                        break;
            
            for(u = 0; u < num_vtd; u++){
                if(u != w){                    
                    if(color[u] == color[w]){
                        if(adj[u*num_vtd + w] == 1){
                        
                        
                        
                        old_component = component[w]
                        new_component = component[v]                        
                        for(u = 0; u <= w; u++){
                            if(component[u] == old_component){
                                component[u] = new_component
                            }                        
                        }
                    }
                }
            }
        }




int v, w;
        num_components = 0;
        for(v = 0; v < num_vtd; v++){
            component[v] = num_components;
            for(w = v; v < num_vtd; v++){
                if(color[v] == color[w]){
                    if(adj[v*num_vtd + w] == 1){
                        component[w] = component[v];
                        num_components--;
                        break;
                    }
                }
            }
            num_components++;
        }


}
""")
compute_components = mod.get_function("compute_components_gpu")

blockDims = (1,1,1)
gridDims = (1,1,1)
compute_components(component_gpu, np.uint16(0), color_gpu, np.uint16(num_vtd), adj_gpu, block=blockDims, grid=gridDims, shared=0)
print(color)
print(component_gpu)
print(num_components_gpu)