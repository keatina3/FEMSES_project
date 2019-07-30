#include <vector>
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
