#include <vector>
#include "mesh.h"
#include "fem.h"

// this also needs parameters //
// change a,b to x, y maybe?? //
FEM::FEM(const int* nr, const int* a, const int* b){
    // need to incorporate DOF here also //
    // dim L = (#Nodes + #DOF/node)^2
    // FIX DEGREES OF FREEDOM AND NUM DIMENSIONS //

    int num_nodes = (nr[0]+1)*(nr[1]+1);
    int dim = num_nodes + 0;
    
    L.resize(dim, std::vector<float>(dim,0.0)); 
    u.resize(dim, 0.0); 
    this->b.resize(dim, 0.0); 
    
    Le.resize(dim, std::vector<std::vector <float> >(3, std::vector<float>(3,0.0)));
    be.resize(dim, std::vector<float>(3,0.0));

    M = new Mesh(nr, a, b);
}

void FEM::assemble(){
    int dof_r, dof_s;
    
    // change all these iterations to .size() for genralisation later //
    for(unsigned int e=0; e<Le.size(); e++){
        elem_mat(e);
        for(unsigned int r=0; r<Le[e].size(); r++){
            dof_r = M->dof_map(e,r);
            for(unsigned int s=0; s<Le[e][r].size(); s++){
                dof_s = M->dof_map(e,s);
                L[dof_r][dof_s] += Le[e][r][s];
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
    // add BLAS library call here //
}

// THIS FUNCTION IS VERY INEFFICIENT // 
// NEED TO REDUCE UNECESSARY FN CALLS //
void FEM::elem_mat(const int e) {
    // these need to be expanded for more DOF //
    float alpha[3], beta[3], gamma[3];
    float del;
    float xi[3][3];
    int v;
    
    for(int i=0;i<3;i++){
        xi[i][0] = 1.0;
        M->get_xy(&xi[i][1], M->get_vertex(e,i));
    }

    // check this //
    del = e%2 == 0 ? area(xi) : (-1)*area(xi);

    for(unsigned int i=0; i<Le[e].size(); i++){
        // alpha = xi[(i+1)%3][1] * xi[(i+2)%3][2] - xi[(i+2)%3][1] * xi[(i+1)%3][2];
        beta[i] = xi[(i+1)%3][2] - xi[(i+2)%3][2];
        gamma[i] = xi[(i+1)%3][1] - xi[(i+2)%3][1];
        
        v = M->get_vertex(e,i);
        // need fn for 
        // also fn ptr here for Int(fv) //
        be[e][i] = M->get_bound(v);
       
        // remove some operations from this // 
        for(unsigned int j=0; j<Le[e].size(); j++){
            // possibly change this to fn ptr //
            Le[e][i][j] = (1/4.0) * del*del*del * (beta[i]*beta[j] + gamma[i]*gamma[j]);
            if(M->get_bound(v) && i != j){
                be[e][j] -= Le[e][j][i] * M->get_bound(v);
                Le[e][i][j] = 0.0;
                Le[e][j][i] = 0.0;
            }
        }
        
        if(M->get_bound(v)){
            be[e][i] = M->get_bound(M->get_vertex(e,i));
            Le[e][i][i] = 1.0;
        }
    }
}

// potentially expand this to numerical integration instead //
float FEM::area(float xi[3][3]) const {
    float tmp = 0.0;

    tmp += xi[0][0] * (xi[1][1]*xi[2][2] - xi[1][2]*xi[2][1]);
    tmp -= xi[0][1] * (xi[1][0]*xi[2][2] - xi[1][2]*xi[2][0]);
    tmp += xi[0][2] * (xi[1][0]*xi[2][1] - xi[1][1]*xi[2][0]);

    return 2*tmp;
}
