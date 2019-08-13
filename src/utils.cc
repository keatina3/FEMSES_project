#include <vector>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include "mesh.h"
#include "utils.h"

bool verbose = false, timing = false, cpu = true, gpu_f = true, gpu_fs = true;
bool annulus = true, dense = false, dnsspr = false, debug =  false;
int n = 2, m = 2;
float a = 3.0, dr = 7.0;

int parse_arguments(int argc, char **argv){
    int opt;

    while((opt = getopt(argc, argv, "vhtcgfdsDCn:m:a:r:")) != -1) {
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
            default:
                fprintf(stderr, "Invalid option given\n");
                print_usage();
                return -1;
        }
    }
    return 0;
}

void print_usage(){
    printf("Finite element method - single element solution GPU program\n");
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
    printf("    -n          : number of rectangles in x-axis\n");
    printf("    -m          : number of rectangles in y-axis\n");
    printf("    -a          : radius of inner circle in annulus/left corner of rectangle\n");
    printf("    -r          : difference between inner radius and outer radius\n");
    printf("        \n");
}

float sse(float *a, float *b, int n){
    float sse = 0.0;

    for(int i=0; i<n; i++)
        sse += (a[i] -b[i]) * (a[i] - b[i]);

    return sse;
}

//// FIX THIS //////
void output_csv(char *fname, Mesh &M, float *u, int order){
    FILE* fptr;
    float xy[2];

    fptr = fopen(fname,"w");
    if(!fptr)
        printf("Couldn't open file %s\n",fname);

    fprintf(fptr, "x, y, u(x,y)\n");
    for(int v=0; v<order; v++){
        M.get_xy(xy, v);
        fprintf(fptr,"%f, %f, %f\n", xy[0], xy[1], u[v]);
    }

    fclose(fptr);
}

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
