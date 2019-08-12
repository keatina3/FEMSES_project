#include <cstring>
#include <set>
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <cstdio>
// #include <gsl/gsl_linalg.h>
#include <mkl.h>
#include "mkl_types.h"
#include "mesh.h"
#include "utils.h"
#include "fem.h"

// this also needs parameters //
// change a,b to x, y maybe?? //
FEM::FEM(Mesh &M){
    // need to incorporate DOF here also //
    // dim L = (#Nodes + #DOF/node)^2
    // FIX DEGREES OF FREEDOM AND NUM DIMENSIONS //
    int nr[2];
    int num_nodes, dim;
    
    M.get_recs(nr);
    
    num_nodes = (nr[0]+1)*(nr[1]+1);
    dim = 2+1;
    order = num_nodes + 0;
    num_cells = 2*nr[0]*nr[1];
    
    this->M = &M;
    
    if(dense){
        L = new float*[order];
        L_vals = new float[order*order]();

        for(int i=0;i<order;i++)
            L[i] = &L_vals[i*order];

    } else {
        // this possibly shouldn't return a value, change to passing by reference maybe //
        this->nnz = M.sparsity_pass_half(valsL, rowPtrL, colIndL);
    }
 
    this->b = new float[order]();
    Le.resize(num_cells, std::vector<std::vector <float> >(dim, std::vector<float>(dim,0.0)));
    be.resize(num_cells, std::vector<float>(dim,0.0));
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
        
    for(unsigned int i=0; i<9; i++){
        for(unsigned int j=0; j<9; j++){
            std::cout << L[i][j] << " ";
        }
        std::cout << "bi = " << b[i] << std::endl;
    }     
}

void FEM::assemble_csr(){
    int dof_r, dof_s;
    int off = 0;
    int *tmp;
    
    for(unsigned int e=0; e<Le.size(); e++){
        // can this be put in here ?? //
        //elem_mat(e);
        for(unsigned int r=0; r<Le[e].size(); r++){
            dof_r = M->dof_map(e,r);
            for(unsigned int s=0; s<Le[e][r].size(); s++){
                dof_s = M->dof_map(e,s);
                if(dof_s < dof_r)   continue;       // dealing with half matrix //
                
                tmp = &colIndL[rowPtrL[dof_r]];
                while(*tmp != dof_s){
                    off++;
                    tmp++;
                }
                valsL[rowPtrL[dof_r] + off] += Le[e][r][s];
                off = 0;
            }
            b[dof_r] += be[e][r];
        }
    }
    std::cout << "TEST\n";
    print_csr(order, &valsL[0], &rowPtrL[0], &colIndL[0]);
}

// change this to CSC format in time //
void FEM::solve(){
    // add LAPACK library call here //
    MKL_INT n = order, nrhs = 1, lda = order, ldb = 1, info;
    // MKL_INT ipiv[order];
    
    ///////// GENERATING ELEMENT MATRICES //////////////
    for(unsigned int e=0; e<Le.size(); e++)
        elem_mat(e);
    
    /////////////////////////////////////////////////////

    /////////// ASSEMBLING GLOBAL STIFFNESS MATRIX //////
    if(dense) 
        assemble();
    else 
        assemble_csr();

    /////////////////////////////////////////////////////
    
    //////////// SOLVING LINEAR SYSTEM //////////////////
    if(dense){
        info = LAPACKE_sposv(LAPACK_ROW_MAJOR, 'L', n, nrhs, L_vals, lda, b, ldb);
        assert(!info);
    } else
        MKL_solve();

    /////////////////////////////////////////////////////
}

void FEM::MKL_solve(){
    MKL_INT status = MKL_DSS_SUCCESS;
    const MKL_INT nRows = order;
    const MKL_INT nCols = order;
    const MKL_INT nNonZeros = nnz;
    const MKL_INT nRhs = 1;
    const _INTEGER_t *rowInd = &rowPtrL[0];
    const _INTEGER_t *columns = &colIndL[0];
    const _REAL_t *values = &valsL[0];

    const _REAL_t *rhs = &b[0];
    _REAL_t *solVals = new float[order]();
    _MKL_DSS_HANDLE_t handle;
    MKL_INT opt = MKL_DSS_DEFAULTS;
    MKL_INT sym = MKL_DSS_SYMMETRIC;
    MKL_INT type = MKL_DSS_POSITIVE_DEFINITE;
    
    std::cout << "order = " << order << std::endl;
    // need option here for double precision versions // 
    opt = MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR 
                        + MKL_DSS_SINGLE_PRECISION + MKL_DSS_ZERO_BASED_INDEXING;
    status = dss_create(handle, opt);
    assert(status == MKL_DSS_SUCCESS);
    
    std::cout << "MKL" << std::endl; 
    status = dss_define_structure(handle, sym, rowInd, nRows, nCols, columns, nNonZeros);
    assert(status == MKL_DSS_SUCCESS);

    std::cout << "MKL3" << std::endl; 
    opt = MKL_DSS_DEFAULTS;
    status = dss_reorder(handle, opt, 0); 
    assert(status == MKL_DSS_SUCCESS);
    
    std::cout << "MKL" << std::endl; 
    status = dss_factor_real(handle, type, values);
    assert(status == MKL_DSS_SUCCESS);
    
    // will this work for same b?? //
    std::cout << "MKL" << std::endl; 
    status = dss_solve_real(handle, opt, rhs, nRhs, solVals);
    assert(status == MKL_DSS_SUCCESS);
    
    std::cout << "MKL" << std::endl; 
    status = dss_delete(handle, opt);
    assert(status == MKL_DSS_SUCCESS);

    memcpy(b, solVals, order*sizeof(float));
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
    // CHANGE ORDERING OF NODES //
    //del = e%2 == 0 ? area(xi) : (-1)*area(xi);
    del = area(xi);
    
    // THIS WONT WORK FOR P2, P3 etc //
    for(unsigned int i=0; i<Le[e].size(); i++){
        // alpha = xi[(i+1)%3][1] * xi[(i+2)%3][2] - xi[(i+2)%3][1] * xi[(i+1)%3][2];
        beta[i] = xi[(i+1)%3][2] - xi[(i+2)%3][2];
        gamma[i] = xi[(i+2)%3][1] - xi[(i+1)%3][1];
        
        // be[e][i] = INT(fv) //
        be[e][i] = 0.0;

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
            //std::cout << v << " Boundary present\n";
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
    
    output_csv(fname, *M, b, order);
}
