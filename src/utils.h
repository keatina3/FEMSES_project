#ifndef _UTILS_H_
#define _UTILS_H_

extern bool verbose, timing, cpu, gpu_f, gpu_fs;
extern bool annulus, dense, dnsspr, debug;
extern int n, m; 
extern float a, dr;
//extern int block_size_X, block_size_Y;

typedef struct Tau {
    float sparsity_scan;    // sparse vs dense //
    float tot;
    float alloc;            // sparse vs dense && cpu vs gpu (both versions) //
    float transfer;         // csr vs dense, forward & back //
    float elem_mats;        // serial vs embarassingly parallel (pot try other blck sizes) //
    float assembly;         // csr vs dense, serial vs parallel //
    float solve;            // csr vs dense vs dnsspr //
    float convert;
} tau, *tau_ptr;

// need new output functions //
// need to fix fname string parse //
int parse_arguments(int argc, char **argv);
void print_usage();

float sse(float *a, float *b, int n);
void output_csv(char *fname, Mesh &M, float *u, int order);
void assign_ptrs(float*** arr_ptr, float** arr, int n, int m);
void assign_ptrs(int*** arr_ptr, int** arr, int n, int m);
void print_csr(int m, const float *csrValA, const int *csrRowPtrA, 
                            const int *csrColIndA);

#endif
