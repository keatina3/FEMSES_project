#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <cstdio>
#include "mkl.h"
#include "mesh.h"
#include "fem.h"

// this also needs parameters //
// change a,b to x, y maybe?? //
FEM::FEM(Mesh *M){
    // need to incorporate DOF here also //
    // dim L = (#Nodes + #DOF/node)^2
    // FIX DEGREES OF FREEDOM AND NUM DIMENSIONS //
    int nr[2];
    int num_nodes, dim;
    
    M->get_recs(nr);
    
    num_nodes = (nr[0]+1)*(nr[1]+1);
    dim = 2+1;
    order = num_nodes + 0;
    num_cells = 2*nr[0]*nr[1];
    
    L_vals = new float[order*order]();
    L = new float*[order];

    for(int i=0;i<order;i++)
        L[i] = &L_vals[i*order];

    this->b = new float[order]();

    Le.resize(num_cells, std::vector<std::vector <float> >(dim, std::vector<float>(dim,0.0)));
    be.resize(num_cells, std::vector<float>(dim,0.0));

    //M = new Mesh(nr, a, b);
    this->M = M;
}

FEM::~FEM(){
    delete[] L;
    delete[] L_vals;
    delete[] b;
}

void FEM::assemble(){
    int dof_r, dof_s;
    
    // change all these iterations to .size() for genralisation later //
    for(unsigned int e=0; e<Le.size(); e++){
        // can this be put in here ?? //
        //elem_mat(e);
        for(unsigned int r=0; r<Le[e].size(); r++){
            dof_r = M->dof_map(e,r);
            for(unsigned int s=0; s<Le[e][r].size(); s++){
                dof_s = M->dof_map(e,s);
                L[dof_r][dof_s] += Le[e][r][s];
                //std::cout << e << " " << dof_r << " " << dof_s << " " << std::endl;
            }
            b[dof_r] += be[e][r];
        }
    }
}

/*
float FEM::phi_P1(const float* x, const float del) const {

}
*/

// change this to CSC format in time //
void FEM::solve(){
    // add LAPACK library call here //
    MKL_INT n = order, nrhs = 1, lda = order, ldb = 1, info;
    MKL_INT ipiv[order];

    for(unsigned int e=0; e<Le.size(); e++)
        elem_mat(e);
     
    /*
    for(int e=0; e<num_cells; e++)
        std::cout << M->dof_map(e,0) << " " << M->dof_map(e,1) << " " << M->dof_map(e,2) << std::endl;
        //std::cout << M->get_vertex(e,0) << " " << M->get_vertex(e,1) << " " << M->get_vertex(e,2) << std::endl;

    for(int i=0;i<Le[0].size();i++){
        for(int j=0;j<Le[0][0].size();j++){
            std::cout << Le[1][i][j] << " ";
        }
        std::cout << "be[i] = " << be[1][i] << 
            " Boundary = " << M->get_bound(M->get_vertex(1,i)) << std::endl;
    }
    */

    assemble();    
   
     
    // testing symmetry // 
    for(int i=0; i<order; i++){
        for(int j=0; j<order; j++){
            if(L[i][j] != L[j][i])
                std::cout << L[i][j] << std::endl;
            //std::cout << L[i][j] << " ";
        }
        //std::cout << std::endl;
    }
    
    info = LAPACKE_ssysv(LAPACK_ROW_MAJOR, 'L', n, nrhs, L_vals, lda, ipiv, b, ldb);
}

// THIS FUNCTION IS VERY INEFFICIENT // 
// NEED TO REDUCE UNECESSARY FN CALLS //
void FEM::elem_mat(const int e) {
    // these need to be expanded for more DOF //
    float alpha[3], beta[3], gamma[3];
    float del, bound;
    float xi[3][3];
    int v;
    bool is_bound;

    for(int i=0;i<3;i++){
        xi[i][0] = 1.0;
        M->get_xy(&xi[i][1], M->get_vertex(e,i));
    }

    // check this //
    del = e%2 == 0 ? area(xi) : (-1)*area(xi);
    //del = area(xi);
    
    // THIS WONT WORK FOR P2, P3 etc //
    for(unsigned int i=0; i<Le[e].size(); i++){
        // alpha = xi[(i+1)%3][1] * xi[(i+2)%3][2] - xi[(i+2)%3][1] * xi[(i+1)%3][2];
        beta[i] = xi[(i+1)%3][2] - xi[(i+2)%3][2];
        gamma[i] = xi[(i+1)%3][1] - xi[(i+2)%3][1];
        
        // need fn for 
        // also fn ptr here for Int(fv) //
        // be[e][i] = INT(fv) //
        //be[e][i] = M->get_bound(M->get_vertex(e,i));
        be[e][i] = 0.0;;


        // remove some operations from this // 
        for(unsigned int j=0; j<=i; j++){
            // possibly change this to fn ptr //
            Le[e][i][j] = 0.25 * del*del*del * (beta[i]*beta[j] + gamma[i]*gamma[j]);
            Le[e][j][i] = Le[e][i][j];
        }
    }

    for(unsigned int i=0; i<Le[e].size(); i++){
        v = M->get_vertex(e,i);
        is_bound = M->is_bound(v);
        
        if(is_bound){
            bound = M->get_bound(v);
            std::cout << v << " Boundary present\n";
            for(unsigned int j=0; j<Le[e][i].size(); j++){
                if(i != j){
                    // WILL THIS WORK AFTER Le = 0.0 FOR 2ND BOUNDARY ??? //
                    be[e][j] -= Le[e][j][i] * bound;
                    //std::cout << "i = " << i << " j = " << j << std::endl;
                    Le[e][i][j] = 0.0;
                    Le[e][j][i] = 0.0;
                }
            }
            Le[e][i][i] = 1.0;
            be[e][i] = bound;
        }
        //if(e==0){
        //    std::cout << be[e][i] << std::endl;
        //}
    }
}

// potentially expand this to numerical integration instead //
float FEM::area(float xi[3][3]) const {
    float tmp = 0.0;

    tmp += xi[0][0] * (xi[1][1]*xi[2][2] - xi[1][2]*xi[2][1]);
    tmp -= xi[0][1] * (xi[1][0]*xi[2][2] - xi[1][2]*xi[2][0]);
    tmp += xi[0][2] * (xi[1][0]*xi[2][1] - xi[1][1]*xi[2][0]);

    return 0.5*tmp;
}

void FEM::output(char* fname) const {
    FILE* fptr;
    float xy[2];

    fptr = fopen(fname, "w");
    if(!fptr)
        printf("Couldn't open file %s\n", fname);
    
    // ORDER OR NUM_NODES ?? //
    fprintf(fptr, "x, y, z\n");
    for(int v=0; v<order; v++){
        M->get_xy(xy, v);
        fprintf(fptr, "%f, %f, %f\n", xy[0], xy[1], b[v]);
    }

    fclose(fptr);
}
