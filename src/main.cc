// ============================================================================= //
// Name:            fem_solver                                                   //
// Author:          Alex Keating                                                 //
// Version:         02/09/19                                                     //
// Description:     Performs the finite element method on the poisson equation,  //
//                  utilising serial, standard GPU and FEMSES GPU approaches.    //
//                                                                               //
// WRITTEN FOR COMPLETION OF MSC. HIGH PERFORMANCE COMPUTING RESEARCH PROJECT    //
// ============================================================================= //

#include <cassert>
#include <ctime>
#include <cmath>
#include <set>
#include <iostream>
#include <vector>
#include <chrono>
#include "mesh.h"
#include "utils.h"
#include "fem.h"

extern void dummy(float *dat, int n);
extern void gpu_fem(float *u, Mesh &M, Tau &t);
extern void gpu_femses(float *u, Mesh &M, Tau &t);

int main(int argc, char** argv){
    int nr[2];
    float x[2], y[2];
    float *u, *u_gpu, *u_gpu_femses;
    float sse_cpu, sse_gpuf, sse_gpufs;
    int order;
    Tau tau_cpu = tau_default, tau_gpu_f = tau_default, tau_gpu_fs = tau_default;

    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    parse_arguments(argc, argv);

    if(verbose) init_screen();

    nr[0] = n, nr[1] = m;
    x[0] = a; x[1] = a + dr;
    
    if(annulus)     y[0] = 0.0, y[1] = 1.0;
    else            y[0] = a, y[1] = a + dr;

    order = (n+1)*(m+1);
    
    // FIXME // 
    /*
    if(order >= 5E4 && dense){
        std::cout << "Too many unknowns to create dense matrix. Changing to sparse solver\n";
        dense = false;
    }
    */
    if(order >= MAX_UNKNOWNS){
        std::cerr << "Problem too large. Exiting.\n";
        std::exit(1);
    }

    
    Mesh M(nr,x,y);
    if(annulus)     M.deform(annulus_seg_map, 1.0);
    //M.deform(annulus_seg_map, -M_PI/6);

    u = new float[order];
    analytical(u, M, x[0], x[1], order);
   
    //////////////      CPU       ///////////////////
 
    if(cpu){
        start = std::chrono::high_resolution_clock::now();
        
        FEM F(M, tau_cpu);
        F.solve(tau_cpu);
        
        end = std::chrono::high_resolution_clock::now(); 
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        tau_cpu.tot = duration.count();
        
        F.output(u);
        sse_cpu = F.sse_fem(u);
    }
    
    /////////////////////////////////////////////////


    if(gpu_f || gpu_fs)     dummy(u, order);        // dummy kernel to prevent slowdown


    //////////////      GPU       ///////////////////
    
    if(gpu_f){
        start = std::chrono::high_resolution_clock::now();
        
        u_gpu = new float[order]();
        gpu_fem(u_gpu, M, tau_gpu_f);
        
        end = std::chrono::high_resolution_clock::now(); 
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        tau_gpu_f.tot = duration.count();
        
        output_results(M, u, u_gpu, order, 1, sse_gpuf);
        sse_gpuf = sse(u, u_gpu, order);
    }
        
    /////////////////////////////////////////////////
    
    
    /////////////      FEMSES       /////////////////
    
    if(gpu_fs){
        start = std::chrono::high_resolution_clock::now();
        
        u_gpu_femses = new float[order](); 
        gpu_femses(u_gpu_femses, M, tau_gpu_fs);
    
        end = std::chrono::high_resolution_clock::now(); 
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        tau_gpu_fs.tot = duration.count();
        
        output_results(M, u, u_gpu_femses, order, 2, sse_gpufs);
        sse_gpufs = sse(u, u_gpu_femses, order);
    }
    
    /////////////////////////////////////////////////
     
     
    if(verbose) output(tau_cpu, tau_gpu_f, tau_gpu_fs, sse_cpu, sse_gpuf, sse_gpufs);

    if(timing){
        if(cpu)     output_times(tau_cpu, 0);
        if(gpu_f)   output_times(tau_gpu_f, 1);
        if(gpu_fs)  output_times(tau_gpu_fs, 2);
    }

    delete[] u;
    if(gpu_f)   delete[] u_gpu; 
    if(gpu_fs)  delete[] u_gpu_femses;

    return 0;
}
