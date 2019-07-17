#include <iostream>
#include <vector>
#include "mesh.h"
#include "fem.h"

int main(int argc, char** argv){
    int nr[2];
    float a[2], b[2];

    nr[0] = 8, nr[1] = 8;
    a[0] = 3; a[1] = 0;
    b[0] = 10; b[1] = 1;

    Mesh M(nr,a,b);
    M.deform(annulus_seg_map);

    FEM F(&M);
    F.solve();
    F.output("output.csv");
     
    return 0;
}
