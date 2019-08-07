#include <vector>
#include <cstdlib>
#include <cstdio>
#include "mesh.h"
#include "utils.h"

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

void printCsr(int m, const float *csrValA, const int *csrRowPtrA, const int *csrColIndA){

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
