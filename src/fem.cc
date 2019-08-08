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
    

    L_vals = new float[order*order]();
    L = new float*[order];

    for(int i=0;i<order;i++)
        L[i] = &L_vals[i*order];

    this->b = new float[order]();
    
    // add if statements for sparse/dense //
    Le.resize(num_cells, std::vector<std::vector <float> >(dim, std::vector<float>(dim,0.0)));
    be.resize(num_cells, std::vector<float>(dim,0.0));

    //M = new Mesh(nr, a, b);
    this->M = &M;

    // this possibly shouldn't return a value, change to passing by reference maybe //
    this->nnz = M.sparsity_pass(valsL, rowPtrL, colPtrL);

    std::cout <<"test i\n";
    
    /* 
    for(int i=0; i<=order; i++){
        std::cout << rowPtrL[i] << std::endl;
    }
    
    for(int i=0; i<nnz; i++){
        std::cout << colPtrL[i] << std::endl;
    }
    */
}

FEM::~FEM(){
    // delete[] L;
    // delete[] L_vals;
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
                //std::cout << e << " " << dof_r << " " << dof_s << " " << std::endl;
                tmp = &colPtrL[rowPtrL[dof_r]];
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
    print_csr(order, &valsL[0], &rowPtrL[0], &colPtrL[0]);
}

// change this to CSC format in time //
void FEM::solve(){
    // add LAPACK library call here //
    MKL_INT n = order, nrhs = 1, lda = order, ldb = 1, info;
    // MKL_INT ipiv[order];

    for(unsigned int e=0; e<Le.size(); e++)
        elem_mat(e);
     
    /*
    for(int e=0; e<num_cells; e++)
        std::cout << M->dof_map(e,0) << " " << M->dof_map(e,1) << " " << M->dof_map(e,2) << std::endl;
        //std::cout << M->get_vertex(e,0) << " " << M->get_vertex(e,1) << " " << M->get_vertex(e,2) << std::endl;
    
    for(int i=0;i<Le[0].size();i++){
        for(int j=0;j<Le[0][0].size();j++){
            std::cout << Le[0][i][j] << " ";
        }
        std::cout << "be[i] = " << be[0][i] << 
            " Boundary = " << M->get_bound(M->get_vertex(1,i)) << std::endl;
    }
    */

    assemble();
    assemble_csr(); 
   /* 
    for(unsigned int e=0; e<be.size();e++){
        for(int i=0; i<3; i++){
            std::cout << be[e][i] << " ";
        }
        std::cout << std::endl;
    }

    // testing symmetry // 
    for(int i=0; i<order; i++){
        for(int j=0; j<order; j++){
            if(L[i][j] != L[j][i])
                std::cout << L[i][j] << std::endl;
            //std::cout << L[i][j] << " ";
        }
        //std::cout << std::endl;
    }
    */
    
    //const gsl_matrix_float_view L_gsl = gsl_matrix_float_view_array(L_vals, order, order);
    //gsl_vector_float_view b_gsl = gsl_vector_float_view_array(b, order);

    // gsl_linalg_cholesky_svx(&L_gsl.matrix, &b_gsl.vector);

    // info = LAPACKE_ssysv(LAPACK_ROW_MAJOR, 'L', n, nrhs, L_vals, lda, ipiv, b, ldb);
    // info = LAPACKE_sgesv(LAPACK_ROW_MAJOR, n, nrhs, L_vals, lda, ipiv, b, ldb);
    
    //info = LAPACKE_sposv(LAPACK_ROW_MAJOR, 'L', n, nrhs, L_vals, lda, b, ldb);
    //assert(!info);

    MKL_solve();
}

void FEM::MKL_solve(){
    MKL_INT status = MKL_DSS_SUCCESS;
    _MKL_DSS_HANDLE_t handle;
    MKL_INT opt;
    MKL_INT* rowInd = &rowPtrL[0];
    MKL_INT* columns = &colPtrL[0];
    MKL_INT nRows = order;
    MKL_INT nCols = order;
    MKL_INT nNonZeros = nnz;
    MKL_INT nRhs = 1;
    
    std::cout << "MKL" << std::endl;
    // need option here for double precision versions // 
    opt = MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR 
                        + MKL_DSS_SINGLE_PRECISION + MKL_DSS_ZERO_BASED_INDEXING;
    status = dss_create(handle, opt);
    assert(status == MKL_DSS_SUCCESS);
    
    std::cout << "MKL" << std::endl; 
    opt = MKL_DSS_SYMMETRIC;
    status = dss_define_structure(handle, opt, rowInd, nRows, nCols, columns, nNonZeros);
    assert(status == MKL_DSS_SUCCESS);

    // printCsr(order, order, nnz, &valsL[0], &rowPtrL[0], &colPtrL[0]);

    std::cout << "MKL3" << std::endl; 
    opt = MKL_DSS_AUTO_ORDER;
    status = dss_reorder(handle, opt, 0); 
    assert(status == MKL_DSS_SUCCESS);
    
    std::cout << "MKL" << std::endl; 
    opt = MKL_DSS_POSITIVE_DEFINITE;
    status = dss_factor_real(handle, opt, &valsL[0]);
    assert(status == MKL_DSS_SUCCESS);
    
    // will this work for same b?? //
    std::cout << "MKL" << std::endl; 
    opt = MKL_DSS_DEFAULTS;
    status = dss_solve_real(handle, opt, b, nRhs, b);
    assert(status == MKL_DSS_SUCCESS);
    
    std::cout << "MKL" << std::endl; 
    status = dss_delete(handle, opt);
    assert(status == MKL_DSS_SUCCESS);
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
        
        // need fn for 
        // also fn ptr here for Int(fv) //
        // be[e][i] = INT(fv) //
        //be[e][i] = M->get_bound(M->get_vertex(e,i));
        be[e][i] = 0.0;


        // remove some operations from this // 
        for(unsigned int j=0; j<=i; j++){
            // possibly change this to fn ptr //
            Le[e][i][j] = 0.25 * del*del*del * (beta[i]*beta[j] + gamma[i]*gamma[j]);
            Le[e][j][i] = Le[e][i][j];
        }
    }
    
    /*
    std::cout << "cell = " << e << " Area = " << del << " Nodes = " << 
        M->get_vertex(e,0) << " " << M->get_vertex(e,1) << " " << M->get_vertex(e,2)
        << " gamma[i] = " << gamma[0] << " " << gamma[1] << " " << gamma[2]
        << " beta[i] = " << beta[0] << " " << beta[1] << " " << beta[2]
        << " Lii = " << Le[e][0][0] << " " << Le[e][1][1] << " " << Le[e][2][2] << std::endl;
    
    int tmp = 7;
    if(e==tmp){
    for(unsigned int i=0;i<Le[0].size();i++){
        for(unsigned int j=0;j<Le[tmp][0].size();j++){
            std::cout << Le[tmp][i][j] << " ";
        }
        std::cout << "be[i] = " << be[tmp][i] << 
            " Boundary = " << M->get_bound(M->get_vertex(tmp,i)) << std::endl;
    }
    }
    */
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
        //if(e==0){
        //    std::cout << be[e][i] << std::endl;
        //}
    }
    
   /* 
    if(e==tmp){
    for(unsigned int i=0;i<Le[0].size();i++){
        for(unsigned int j=0;j<Le[tmp][0].size();j++){
            std::cout << Le[tmp][i][j] << " ";
        }
        std::cout << "be[i] = " << be[tmp][i] << 
            " Boundary = " << M->get_bound(M->get_vertex(tmp,i)) << std::endl;
    }
    }
    */
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

    /*
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
    */
}
/*
int FEM::sparsity_pass(float **valsL, float **rowPtr, float **colPtrL){
    std::vector<std::set<int> > sparsity;
    std::vector<int> colTmp;
    int v1, v2;
    int n = 0;
    int *rowTmp;
    int count=0;

    sparsity.resize(order);
    colTmp.resize(order*order, 0);
    this->rowPtrL.resize(order+1,0);

    rowTmp = &rowPtrL[0];

    for(int e=0; e<num_cells; e++){
        for(int i=0; i<3; i++){
            v1 = M->get_vertex(e, i);
            for(int j=0; j<3; j++){
                v2 = M->get_vertex(e, j);
                sparsity[v1].insert(v2);
                sparsity[v2].insert(v1);
            }
        }
    }
    
    rowTmp++;
    for(std::vector<std::set<int> >::iterator v = sparsity.begin(); v != sparsity.end(); ++v){
        for(std::set<int>::iterator vi = (*v).begin(); vi != (*v).end(); ++vi){
            colTmp[count] = (*vi);
            count++;
        }
        n += (*v).size();
        *(rowTmp++) = n;
    }
    
    this->valsL.resize(n,0.0);
    this->colPtrL.resize(n,0);
    
    std::memcpy(&colPtrL[0], &colTmp[0], n*sizeof(int));
    std::cout << "sparsity test\n";    

    return n;
}
*/
