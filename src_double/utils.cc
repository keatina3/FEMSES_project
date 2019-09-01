#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <mkl.h>
#include "mesh.h"
#include "utils.h"

bool verbose = false, timing = false, cpu = true, gpu_f = true, gpu_fs = true;
bool annulus = false, dense = false, dnsspr = false, debug =  false, mem_config = true;
int n = 2, m = 2, k = 1, block_size_X = 32; 
double a = 3.0, dr = 7.0, ui = 2.0, uo = 6.0;
const struct Tau tau_default = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

////////////////// Parse command line arguments //////////////////////////
int parse_arguments(int argc, char **argv){
    int opt;

    while((opt = getopt(argc, argv, "vhtcgfdskDCMn:m:a:r:i:o:b:")) != -1) {
        switch(opt){
            case 'v': 
                verbose =  true; break;
            case 'h':
                print_usage(); exit(0); break;
            case 't':
                timing = true; break;
            case 'c':
                cpu = false; break;
            case 'g':
                gpu_f = false; break;
            case 'f':
                gpu_fs = false; break;
            case 'd':
                dense = true; break;
            case 's':
                dnsspr = true, dense = true; break;
            case 'k':
                k = 0; break;
            case 'D':
                debug = true; break;
            case 'C':
                annulus = true; break;
            case 'M':
                mem_config = false; break;
            case 'n':
                n = atoi(optarg); break;
            case 'm':
                m = atoi(optarg); break;
            case 'a':
                a = atof(optarg); break;
            case 'r':
                dr = atof(optarg); break;
            case 'i':
                ui = atof(optarg); break;
            case 'o':
                uo = atof(optarg); break;
            case 'b':
                block_size_X = atoi(optarg); break;
            default:
                fprintf(stderr, "Invalid option given\n");
                print_usage();
                exit(1);
        }
    }
    return 0;
}
///////


////////////////////////// Prints usage information /////////////////////////
void print_usage(){
    printf("\n==============================================================================\n");
    printf("\nFinite element method - single element solution GPU program\n");
    printf("This program numerically calculates the solution to a PDE using FEM\n");
    printf("Usage:\n");
    printf("fem_solver [options]\n");
    printf("    -v          : will activate verbose mode for output (default: no)\n");
    printf("    -h          : will show this usage message\n");
    printf("    -t          : will activate timing output (default: no)\n");
    printf("    -c          : will skip CPU testing\n");
    printf("    -g          : will skip standard FEM GPU test\n");
    printf("    -f          : will skip FEMSES GPU test\n");
    printf("    -d          : will use dense linear solvers\\matrix assembly\n");
    printf("    -s          : will use dense assembly & conversion to CSR, sparse solver\n");
    printf("    -k          : turns off Tesla K40, changes GPU to RTX2080 Super\n");
    printf("    -D          : turns on debugging mode\n");
    printf("    -C          : turns on mesh deformation from rectangle to annulus\n");
    printf("    -M          : turns off GPU shared/register memory reconfiguration\n");
    printf("    -n          : number of rectangles in x-axis (default: 2)\n");
    printf("    -m          : number of rectangles in y-axis (deafult: 2)\n");
    printf("    -a          : radius of inner circle in annulus/left corner of rectangle (default: 3)\n");
    printf("    -r          : difference between inner radius and outer radius (default: 7)\n");
    printf("    -i         : inside/left, dirichlet boundary condition (default: 2.0)\n");
    printf("    -o         : outside/right, dirichlet boundary condition (default: 6.0)\n");
    printf("    -b         : sets block_size_X (default: 32)\n");
    printf("\n==============================================================================\n\n");
}
//////


//////////////////// Prints initial screen for program run ////////////////////
void init_screen(){
    printf("\n==============================================================================\n\n");
    printf(BLUE "fem_solver - Finite Element Method GPU Implementation\n");
    printf("Alex Keating\n" RESET);
    printf("\n================================================\n");
    printf(GREEN "Intial setup...\n" RESET);
    printf("                Dirichlet BCs   : %0.00lf, %0.00lf\n", ui, uo);
    printf("                Range           : %0.00lf, %0.00lf\n", a, a+dr);    
    printf("                CPU enabled     : %d\n", cpu);
    printf("                GPU enabled     : %d\n", gpu_f);
    printf("                FEMSES enabled  : %d\n", gpu_fs);   
    if(dnsspr)      printf("                Solver          : Dense-Sparse conv\n");
    else if(dense)  printf("                Solver          : Dense\n");
    else            printf("                Solver          : Sparse\n");
    printf("                Memory reconfig : %d\n", mem_config);
    printf("                block_size_X    : %d\n", block_size_X);
    printf("                # of unknowns   : %d\n", (n+1)*(m+1));   
    printf("================================================\n\n");
}
//////

//////////////////// Prints program run info to stderr file ////////////////////
void error_log(){
    fprintf(stderr, "\n==============================================================================\n\n");
    fprintf(stderr, BLUE "fem_solver - Finite Element Method GPU Implementation\n");
    fprintf(stderr, "Error Log\n" RESET);
    fprintf(stderr, "\n================================================\n");
    fprintf(stderr, GREEN "Intial setup...\n" RESET);
    fprintf(stderr, "                Dirichlet BCs   : %0.00f, %0.00f\n", ui, uo);
    fprintf(stderr, "                Range           : %0.00f, %0.00f\n", a, a+dr);    
    fprintf(stderr, "                CPU enabled     : %d\n", cpu);
    fprintf(stderr, "                GPU enabled     : %d\n", gpu_f);
    fprintf(stderr, "                FEMSES enabled  : %d\n", gpu_fs);   
    if(dnsspr)      fprintf(stderr, "                Solver          : Dense-Sparse conv\n");
    else if(dense)  fprintf(stderr, "                Solver          : Dense\n");
    else            fprintf(stderr, "                Solver          : Sparse\n");
    fprintf(stderr, "                Memory reconfig : %d\n", mem_config);
    fprintf(stderr, "                block_size_X    : %d\n", block_size_X);
    fprintf(stderr, "                # of unknowns   : %d\n", (n+1)*(m+1));   
    fprintf(stderr, "================================================\n\n");

}
///////


/////////////////////// Prints output from program run //////////////////////
void output(tau &t_cpu, tau &t_gpu, tau &t_gpufs, double sse_cpu, double sse_gpu, double sse_gpufs){
    printf("\n================================================\n");
    printf(GREEN "Performance results...\n" RESET);
    printf("               CPU         GPU         FEMSES\n");
    printf("SSE" RED "            %lf" BLUE "    %lf" MAG "    %lf\n" RESET, 
                        sse_cpu, sse_gpu, sse_gpufs);
    if(timing){
    printf("Total(ms)   " RED "   %lf" BLUE "    %lf" MAG "    %lf\n" RESET, 
                        t_cpu.tot, t_gpu.tot, t_gpufs.tot);
    printf("Alloc(ms)   " RED"   %lf" BLUE "    %lf" MAG "    %lf\n" RESET, 
                        t_cpu.alloc, t_gpu.alloc, t_gpufs.alloc);
    printf("Trans(ms)   " RED "   %lf" BLUE "    %lf" MAG "    %lf\n" RESET, 
                        t_cpu.transfer, t_gpu.transfer, t_gpufs.transfer);
    printf("Elem Mats(ms)" RED "  %lf" BLUE "    %lf" MAG "    %lf\n" RESET, 
                        t_cpu.elem_mats, t_gpu.elem_mats, t_gpufs.elem_mats);
    printf("Assembly(ms)" RED "   %lf" BLUE "    %lf" MAG "    %lf\n" RESET, 
                        t_cpu.assembly, t_gpu.assembly, t_gpufs.assembly);
    printf("Solver(ms)  " RED "   %lf" BLUE "    %lf" MAG "    %lf\n" RESET, 
                        t_cpu.solve, t_gpu.solve, t_gpufs.solve);
    printf("Convert(ms) " RED "   %lf" BLUE "    %lf" MAG "    %lf\n" RESET, 
                        t_cpu.convert, t_gpu.convert, t_gpufs.convert);
    printf("Sparsity(ms)" RED "   %lf" BLUE "    %lf" MAG "    %lf\n" RESET, 
                        t_cpu.sparsity_scan, t_gpu.sparsity_scan, t_gpufs.sparsity_scan);
    }
    printf("================================================\n");
    printf("\n==============================================================================\n\n");
}
//////


//////////////////// Returns SSE of two vectors a,b ////////////////////////
double sse(double *u, double *u_hat, int dim){
    const MKL_INT n = dim;
    const MKL_INT incx = 1, incy = 1;
    const double alpha = -1.0;
    double sse = 0.0;
    
    cblas_daxpy(n, alpha, u, incx, u_hat, incy);

    sse = cblas_dnrm2(n, u_hat, incx);

    return sse;
}
//////


//////////////// Analytical solution to annulus problem ////////////////////
void analytical(double *u, Mesh &M, int a, int b, int order){
    double xy[2];
    double r;
    double m;

    m = (uo - ui)/(b-a);
    for(int v=0; v<order; v++){
        if(annulus) {
            r = sqrt(xy[0]*xy[0] + xy[1]*xy[1]);
            u[v] = (uo*log(a) - ui*log(b) + ui*log(r) - uo*log(r)) /(log(a) - log(b));
        } else {
            M.get_xy(xy, v);
            u[v] = m*(xy[0]-a) + ui;
        }
    }
}
//////


//////////////////////// Outputs solution to CSV file ////////////////////
void output_results(Mesh &M, double *u, double *u_hat, int order, int routine){
    FILE *fptr;
    double xy[2];
    std::string fname = "results/output_";
    
    if(routine==0)      fname.append("cpu_");
    else if(routine==1) fname.append("gpu_");
    else                fname.append("femses_");

    if(routine != 0){
        if(k==0)    fname.append("GTX2080_");
        else        fname.append("Tesla_");
    }
    
    if(dnsspr && routine != 0)  fname.append("dnsspr");
    else if(dense)              fname.append("dense");
    else                        fname.append("sparse");
    
    fname.append("_results.csv");
    
    fptr = fopen(&fname[0],"w");
    if(!fptr)
        printf("Couldn't open file %s\n",&fname[0]);

    fprintf(fptr, "x, y, u, u_analytical\n");
    for(int v=0; v<order; v++){
        M.get_xy(xy, v);
        fprintf(fptr,"%lf, %lf, %lf, %lf\n", xy[0], xy[1], u_hat[v], u[v]);
    }

    fclose(fptr);
    std::cout << "      Successfully written to file.\n";
}
////////


////////////////// Checks if existing file is empty //////////////////
int is_empty(FILE *file){
    size_t size;

    fseek(file, 0, SEEK_END);
    size=ftell(file);

    return size ? 0 : 1;
}
///////


/////////////////////// Output timings to file ////////////////////////
void output_times(Tau &t, int routine, double sse, int iters, int reconfig){
    FILE *fptr;
    std::string fname = "timings/";
    
    if(routine==0)      fname.append("cpu_");
    else if(routine==1) fname.append("gpu_");
    else                fname.append("femses_");

    if(routine != 0){
        if(k==0)    fname.append("GTX2080_");
        else        fname.append("Tesla_");
    }
    
    if(dnsspr && routine != 0)  fname.append("dnsspr");
    else if(dense)              fname.append("dense");
    else                        fname.append("sparse");

    fname.append("_times.csv");

    fptr = fopen(&fname[0],"a");
    if(!fptr)
        printf("Couldn't open file %s\n",&fname[0]);

    if(is_empty(fptr))
        fprintf(fptr, "n, m, block_size_X, reconfig, total, allocation, transfer, elem_mats, assembly, solve, convert, sparsity scan, sse, iterations\n");

    fprintf(fptr, "%d, %d, %d, %d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %d\n", 
            n, m, block_size_X, reconfig, t.tot, t.alloc, t.transfer, t.elem_mats, t.assembly, t.solve, t.convert, t.sparsity_scan, sse, iters);

    fclose(fptr);
}
////////


/////////////// Allocates and assigns 2d pointers ////////////////////
void assign_ptrs(double*** arr_ptr, double** arr, int n, int m){
    double **arr_tmp, *tmp_vals;

    *arr_ptr = (double**)malloc(n*sizeof(double*));
    *arr = (double*)calloc(n*m,sizeof(double));
    
    arr_tmp = *arr_ptr;
    tmp_vals = *arr;

    for(int i=0; i<n; i++)
        arr_tmp[i] = &tmp_vals[i*m];
}

void assign_ptrs(int*** arr_ptr, int** arr, int n, int m){
    int **arr_tmp, *tmp_vals;
    
    *arr_ptr = (int**)malloc(n*sizeof(int*));
    *arr = (int*)calloc(n*m,sizeof(int));

    arr_tmp = *arr_ptr;
    tmp_vals = *arr;

    for(int i=0; i<n; i++)
        arr_tmp[i] = &tmp_vals[i*m];
}
///////


////////////////// Prints out CSR matrix /////////////////////////////
void print_csr(int m, const double *csrValA, const int *csrRowPtrA, const int *csrColIndA){

    for(int row = 0; row < m; row++){
        const int start = csrRowPtrA[row  ];
        const int end   = csrRowPtrA[row+1];
        for(int colidx = start ; colidx < end ; colidx++){
            const int col = csrColIndA[colidx];
            const double Areg = csrValA[colidx];
            printf("(%d,%d) = %lf\n", row, col, Areg);
        }
    }
}
///////

