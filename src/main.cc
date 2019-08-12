#include <set>
#include <iostream>
#include <vector>
#include "mesh.h"
#include "fem.h"
#include "utils.h"

extern void gpu_fem(float *u, Mesh &M);
extern void gpu_femses(float *u, Mesh &M);

int main(int argc, char** argv){
    int nr[2];
    float x[2], y[2];
    float *u_gpu, *u_gpu_femses;
    int order;
    
    parse_arguments(argc, argv);

    nr[0] = n, nr[1] = m;
    x[0] = a; x[1] = a + dr;

    if(annulus)
        y[0] = 0.0, y[1] = 1.0;
    else
        y[0] = a, y[1] = a + dr;

    order = (nr[0]+1)*(nr[0]+1);
    
    Mesh M(nr,x,y);
    M.deform(annulus_seg_map);
    
    if(cpu){ 
        FEM F(M);
        F.solve();
        F.output("output_cpu.csv");
    }
    
    if(gpu_f){ 
        u_gpu = new float[order];
        gpu_fem(u_gpu, M);
        output_csv("output_gpu.csv", M, u_gpu, order);
    }
    
    if(gpu_fs){
        u_gpu_femses = new float[order]; 
        gpu_femses(u_gpu, M);
        output_csv("output_femses.csv", M, u_gpu_femses, order);
    }

    delete[] u_gpu; delete[] u_gpu_femses;

    return 0;
}
