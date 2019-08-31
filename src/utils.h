// =========================================================================== //
// General utility functions for use with the FEM program.                     //
// =========================================================================== //

#ifndef _UTILS_H_
#define _UTILS_H_

#define MAX_UNKNOWNS 5E5
#define EPS 1.0E-04
#define MAX_ITERS 1E6

#define RED  "\x1B[91m"     // for colour setting in print
#define GREEN  "\x1B[92m"
#define BLUE  "\x1B[94m"
#define MAG  "\x1B[95m"
#define RESET "\x1B[0m"

// command line arguments - see usage function for explanations //
extern bool verbose, timing, cpu, gpu_f, gpu_fs;
extern bool annulus, dense, dnsspr, debug, mem_config;
extern int n, m, k, block_size_X; 
extern float a, dr, ui, uo;

// struct containing all the various timings    //
typedef struct Tau {
    float sparsity_scan;    // sparse vs dense //
    float tot;
    float alloc;            // sparse vs dense && cpu vs gpu (both versions) //
    float transfer;         // csr vs dense, forward & back //
    float elem_mats;        // serial vs embarassingly parallel (pot try other blck sizes) //
    float assembly;         // csr vs dense, serial vs parallel //
    float assem_p_elem;       // time for elem_mats + assembly
    float solve;            // csr vs dense vs dnsspr //
    float convert;
} tau, *tau_ptr;

extern const struct Tau tau_default;

int parse_arguments(int argc, char **argv);
void print_usage();
void init_screen();
void error_log();
void output(Tau &t_cpu, Tau &t_gpu, Tau &t_gpufs, float sse_cpu, float sse_gpu, float sse_gpufs);
float sse(float *a, float *b, int n);
void analytical(float *u, Mesh &M, int a, int b, int order);
void output_results(Mesh &M, float *u, float *u_hat, int order, int routine);
void output_times(Tau &t, int routine, float sse, int iters, int reconfig);
int is_empty(FILE *fptr);
void assign_ptrs(float*** arr_ptr, float** arr, int n, int m);
void assign_ptrs(int*** arr_ptr, int** arr, int n, int m);
void print_csr(int m, const float *csrValA, const int *csrRowPtrA, 
                            const int *csrColIndA);

#endif
