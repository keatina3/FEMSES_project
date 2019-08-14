#include <cmath>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include "mesh.h"
#include "utils.h"

bool verbose = false, timing = false, cpu = true, gpu_f = true, gpu_fs = true;
bool annulus = true, dense = false, dnsspr = false, debug =  false;
int n = 2, m = 2; 
float a = 3.0, dr = 7.0, ui = 2.0, uo = 6.0;

int parse_arguments(int argc, char **argv){
    int opt;

    while((opt = getopt(argc, argv, "vhtcgfdsDCn:m:a:r:i:o:")) != -1) {
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
            case 'D':
                debug = true; break;
            case 'C':
                annulus = false; break;
            case 'n':
                n = atoi(optarg); break;
            case 'm':
                m = atoi(optarg); break;
            case 'a':
                a = atoi(optarg); break;
            case 'r':
                dr = atoi(optarg); break;
            case 'i':
                ui = atoi(optarg); break;
            case 'o':
                uo = atoi(optarg); break;
            default:
                fprintf(stderr, "Invalid option given\n");
                print_usage();
                return -1;
        }
    }
    return 0;
}

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
    printf("    -D          : turns onn debugging mode\n");
    printf("    -C          : turns off mesh deformation from rectangle to annulus\n");
    printf("    -n          : number of rectangles in x-axis (default: 2)\n");
    printf("    -m          : number of rectangles in y-axis (deafult: 2)\n");
    printf("    -a          : radius of inner circle in annulus/left corner of rectangle (default: 3)\n");
    printf("    -r          : difference between inner radius and outer radius (default: 7)\n");
    printf("    -ui         : inside/left, dirichlet boundary condition (default: 2.0)\n");
    printf("    -uo         : outside/right, dirichlet boundary condition (default: 6.0)\n");
    printf("\n==============================================================================\n\n");
}

//////////////////// Returns SSE of two vectors a,b ////////////////////////
float sse(float *a, float *b, int n){
    float sse = 0.0;

    for(int i=0; i<n; i++)
        sse += (a[i] -b[i]) * (a[i] - b[i]);

    return sse;
}
//////


//////////////// Analytical solution to annulus problem ////////////////////
void analytical(float *u, Mesh &M, int a, int b, int order){
    float xy[2];
    float r;

    for(int v=0; v<order; v++){
        M.get_xy(xy, v);
        r = sqrt(xy[0]*xy[0] + xy[1]*xy[1]);
        u[v] = (uo*log(a) - ui*log(b) + ui*log(r) - uo*log(r)) /(log(a) - log(b));
    }
}
//////


///////////////////////// FIXME //////
void output_csv(char *fname, Mesh &M, float *u, float *u_an, int order){
    FILE* fptr;
    float xy[2];

    fptr = fopen(fname,"w");
    if(!fptr)
        printf("Couldn't open file %s\n",fname);

    fprintf(fptr, "x, y, u(x,y)\n");
    for(int v=0; v<order; v++){
        M.get_xy(xy, v);
        fprintf(fptr,"%f, %f, %f, %f\n", xy[0], xy[1], u[v], u_an[v]);
    }

    fclose(fptr);
}
////////


/////////////// Allocates and assigns 2d pointers ////////////////////
void assign_ptrs(float*** arr_ptr, float** arr, int n, int m){
    float **arr_tmp, *tmp_vals;

    *arr_ptr = (float**)malloc(n*sizeof(float*));
    *arr = (float*)calloc(n*m,sizeof(float));
    
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
void print_csr(int m, const float *csrValA, const int *csrRowPtrA, const int *csrColIndA){

    for(int row = 0; row < m; row++){
        const int start = csrRowPtrA[row  ];
        const int end   = csrRowPtrA[row+1];
        for(int colidx = start ; colidx < end ; colidx++){
            const int col = csrColIndA[colidx];
            const float Areg = csrValA[colidx];
            printf("(%d,%d) = %f\n", row, col, Areg);
        }
    }
}
///////

