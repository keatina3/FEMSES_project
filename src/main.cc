#include <ctime>
#include <cmath>
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
    float *u, *u_gpu, *u_gpu_femses;
    int order;
    tau tau_cpu, tau_gpu_f, tau_gpu_fs;

    parse_arguments(argc, argv);

    nr[0] = n, nr[1] = m;
    x[0] = a; x[1] = a + dr;
    
    //if(annulus)
    //    y[0] = 0.0, y[1] = 1.0;
    //else
        y[0] = a, y[1] = a + 0.25*dr;

    order = (nr[0]+1)*(nr[0]+1);
    
    Mesh M(nr,x,y);
    //M.deform(annulus_seg_map, 1.0);
    M.deform(annulus_seg_map, -M_PI/6);
     
    float xy[2];
    for(int v=0; v<nr[0]*nr[1]*2; v++){
        M.get_xy(xy, v);
        std::cout << xy[0] << " " << xy[1] << std::endl;
    }

    u = new float[order];
    analytical(u, M, x[0], x[1], order);

    if(cpu){ 
        FEM F(M);
        F.solve();
        F.output("output_cpu.csv", u);
    }
    
    if(gpu_f){
        u_gpu = new float[order]();
        gpu_fem(u_gpu, M);
        output_csv("output_gpu.csv", M, u_gpu, u, order);
    }
    
    if(gpu_fs){
        u_gpu_femses = new float[order](); 
        gpu_femses(u_gpu_femses, M);
        output_csv("output_femses.csv", M, u_gpu_femses, u, order);
    }

    delete[] u; delete[] u_gpu; delete[] u_gpu_femses;

    return 0;
}
