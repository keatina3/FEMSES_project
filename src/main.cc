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
    float a[2], b[2];
    float *u_gpu, *u_gpu_femses;
    int order;
    
    parse_arguments(argc, argv);

    nr[0] = n, nr[1] = m;
    a[0] = a; b[0] = b;

    if(annulus)
        a[1] = 0; b[1] = 1;
    else
        a[1] = a; b[1] = b;

    order = (nr[0]+1)*(nr[0]+1);

    u_gpu = new float[order];
    u_gpu_femses = new float[order];
    
    Mesh M(nr,a,b);
    M.deform(annulus_seg_map);

    FEM F(M);
    // F.solve();
    // F.output("output_cpu.csv");
    
    gpu_fem(u_gpu, M);
    output_csv("output_gpu.csv", M, u_gpu, order);

    gpu_femses(u_gpu, M);
    output_csv("output_femses.csv", M, u_gpu_femses, order);
    
    delete[] u_gpu; delete[] u_gpu_femses;

    return 0;
}
