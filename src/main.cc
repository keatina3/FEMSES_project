#include <set>
#include <iostream>
#include <vector>
#include "mesh.h"
#include "fem.h"
#include "utils.h"

extern void gpu_fem(float *u, Mesh &M);

int main(int argc, char** argv){
    int nr[2];
    float a[2], b[2];
    float *u_gpu;
    int order;
    
    nr[0] = 2, nr[1] = 2;
    a[0] = 3; a[1] = 0;
    b[0] = 10; b[1] = 1;

    order = (nr[0]+1)*(nr[0]+1);

    u_gpu = new float[order];
    
    Mesh M(nr,a,b);
    // M.deform(annulus_seg_map);

    FEM F(M);
    // F.solve();
    // F.output("output_cpu.csv");
    
    gpu_fem(u_gpu, M);
    output_csv("output_gpu.csv", M, u_gpu, order);

    delete[] u_gpu;

    return 0;
}
