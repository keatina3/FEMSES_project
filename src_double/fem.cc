#include <cstring>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <mkl.h>
#include "mesh.h"
#include "utils.h"
#include "fem.h"

FEM::FEM(Mesh &M, Tau &t){
    int nr[2];
    
    M.get_recs(nr);
    
    order = (nr[0]+1)*(nr[1]+1);
    num_cells = 2*nr[0]*nr[1];
    
    this->M = &M;

    // assigning memory whether dense or CSR format //
    auto start = std::chrono::high_resolution_clock::now();
    if(dense){
        try {
            L = new double*[order];
            L_vals = new double[order*order];
        } catch(std::bad_alloc const &err) {
            error_log();
            std::cerr << "Bad allocation of stiffness matrix\n";
            std::cerr << err.what() << std::endl;
            std::exit(1);
        }
        for(int i=0;i<order;i++)
            L[i] = &L_vals[i*order];
    }
    try {
        this->b = new double[order]();
        Le.resize(num_cells, std::vector<std::vector <double> >(3, std::vector<double>(3,0.0)));
        be.resize(num_cells, std::vector<double>(3,0.0));
    } catch(std::bad_alloc const &err) {
        error_log();
        std::cerr << "Bad allocation of stress vector matrix\n";
        std::cerr << err.what() << std::endl;
        std::exit(1);
    }
    auto end = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    t.alloc = duration.count();
}

FEM::~FEM(){
    if(dense) {
        delete[] L;
        delete[] L_vals;
    }
    delete[] b;
}

////////////////// Assembles global stiffness matrix and stress vector ///////////
// Assembles L in dense format from all the element matrices in the
// vector Le. Same applies for stress vector b and be element vectors 
void FEM::assemble(){
    int dof_r, dof_s;       // global node numbering 
                            // for local nodes r,s
    
    for(unsigned int e=0; e<Le.size(); e++){
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
/////////


////////////////// Assembles global stiffness matrix and stress vector ///////////
// Assembles L in CSR format from all the element matrices in the
// vector Le. Same applies for stress vector b and be element vectors 
// Only assembles upper half of matrix
// Uses the already found sparsity pattern, populated in colInd and rowPtr
void FEM::assemble_csr(){
    int dof_r, dof_s;
    int off = 0;            // offset from start of row
    int *tmp;
    
    for(unsigned int e=0; e<Le.size(); e++){
        for(unsigned int r=0; r<Le[e].size(); r++){
            dof_r = M->dof_map(e,r);
            for(unsigned int s=0; s<Le[e][r].size(); s++){
                dof_s = M->dof_map(e,s);
                
                if(dof_s < dof_r)   continue;       // no lower half values //
                
                tmp = &colIndL[rowPtrL[dof_r]];
                // increments until column dof_s is found in row //
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
}
///////


////////////////////// Overall solver function //////////////////////////////
// Takes FEM object, gets element matrices, assembles L, 
// solves linear system
void FEM::solve(Tau &t){
    MKL_INT n = order, nrhs = 1, lda = order, ldb = 1, info;    
    double alloc_sprs, tau_sprs;

    std::cout << GREEN "\nCPU Solver...\n" RESET;
    
    // See mesh.h for function info //
    if(!dense){  
        std::cout << "      Sparsity pass...\n";
        M->sparsity_pass_half(valsL, rowPtrL, colIndL, this->nnz, alloc_sprs, tau_sprs);
    }
    t.alloc += alloc_sprs;
    t.sparsity_scan = tau_sprs;


    ///////// GENERATING ELEMENT MATRICES //////////////
   
    std::cout << "      Getting element matrices...\n";
    auto start = std::chrono::high_resolution_clock::now();
    
    for(unsigned int e=0; e<Le.size(); e++)
        elem_mat(e);
    
    auto end = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    t.elem_mats = duration.count();
    
    /////////////////////////////////////////////////////


    /////////// ASSEMBLING GLOBAL STIFFNESS MATRIX //////
    
    std::cout << "      Assembling stiffness matrix...\n";
    start = std::chrono::high_resolution_clock::now();
    
    if(dense)   assemble();
    else        assemble_csr();
    
    end = std::chrono::high_resolution_clock::now(); 
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    t.assembly = duration.count();
    
    int nr[2];
    
    M->get_recs(nr);
    
    order = (nr[0]+1)*(nr[1]+1);
    // print_csr(order, &valsL[0], &rowPtrL[0], &colIndL[0]);

    /////////////////////////////////////////////////////
    

    //////////// SOLVING LINEAR SYSTEM //////////////////
    
    std::cout << "      Solving linear system...\n";
    if(!dense) std::cout << "      (nnz = " << nnz << ")\n";
    
    start = std::chrono::high_resolution_clock::now();
    
    if(dense){
        info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'L', n, nrhs, L_vals, lda, b, ldb);
        assert(!info);
    } else
        MKL_solve();
    
    end = std::chrono::high_resolution_clock::now(); 
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    t.solve = duration.count();
    
    for(int i=0;  i<order; i++){
        std::cout << b[i] << std::endl;
    }
    /////////////////////////////////////////////////////
}
//////


/////////////// Linear solver for CSR matrices using MKL's DSS lib //////////////
// Takes stifness matrix in CSR and stress vector b and solves system
// Value written back over b
void FEM::MKL_solve(){
    const MKL_INT nRows = order;
    const MKL_INT nCols = order;
    const MKL_INT nNonZeros = nnz;
    const MKL_INT nRhs = 1;
    const _INTEGER_t *rowInd = &rowPtrL[0];
    const _INTEGER_t *columns = &colIndL[0];
    const _DOUBLE_PRECISION_t *values = &valsL[0];
    const _DOUBLE_PRECISION_t *rhs = &b[0];

    MKL_INT status = MKL_DSS_SUCCESS;
    _DOUBLE_PRECISION_t *solVals = new double[order]();
    _MKL_DSS_HANDLE_t handle;
    MKL_INT opt = MKL_DSS_DEFAULTS;
    MKL_INT sym = MKL_DSS_SYMMETRIC;
    MKL_INT type = MKL_DSS_POSITIVE_DEFINITE;
    
    // setting initial options for solver  and creating handle // 
    opt = MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR + MKL_DSS_ZERO_BASED_INDEXING;
    status = dss_create(handle, opt);
    assert(status == MKL_DSS_SUCCESS);
    
    // defining structure of matrix as symmetric //
    status = dss_define_structure(handle, sym, rowInd, nRows, nCols, columns, nNonZeros);
    assert(status == MKL_DSS_SUCCESS);

    // auto-reordering matrix //
    opt = MKL_DSS_DEFAULTS;
    status = dss_reorder(handle, opt, 0); 
    assert(status == MKL_DSS_SUCCESS);
    
    // factoring matrix using Cholesky decomposition //
    status = dss_factor_real(handle, type, values);
    assert(status == MKL_DSS_SUCCESS);
    
    // solving linear system //
    status = dss_solve_real(handle, opt, rhs, nRhs, solVals);
    assert(status == MKL_DSS_SUCCESS);
    
    // freeing handle //
    status = dss_delete(handle, opt);
    assert(status == MKL_DSS_SUCCESS);

    // copying memory back over b //
    memcpy(b, solVals, order*sizeof(double));

    delete[] solVals;
}
////////


/////////////// Calculates element matrix and element vector for cell e ////////////////
// Uses some analytical formulas based on the fact that 
// each element is a triangle and of order P1
void FEM::elem_mat(const int e){
    double beta[3], gamma[3];    // constants needed for formula
    double del, bound;           // del = area of triangle
    double xi[3][3];             // matrix needed for shoelace equation
    int v;
    bool is_bound;

    // gets x,y coordinates of each of the 3 nodes and fills them into matrix //
    for(int i=0;i<3;i++){
        xi[i][0] = 1.0;
        M->get_xy(&xi[i][1], M->get_vertex(e,i));
    }

    // uses shoelace formula/determinant of matrix xi to get area //
    del = area(xi);
    
    // hand calculates analytical formulae for integral on LHS //
    for(unsigned int i=0; i<Le[e].size(); i++){
        beta[i] = xi[(i+1)%3][2] - xi[(i+2)%3][2];
        gamma[i] = xi[(i+2)%3][1] - xi[(i+1)%3][1];
        
        be[e][i] = 0.0;     // Poisson -> stress vector = 0 unless boundary //

        for(unsigned int j=0; j<=i; j++){
            Le[e][i][j] = 0.25 * del*del*del * (beta[i]*beta[j] + gamma[i]*gamma[j]);
            Le[e][j][i] = Le[e][i][j];      // symmetrical matrix
        }
    }
    
    // applying boundary conditions to element matrix and stress vector //
    for(unsigned int i=0; i<Le[e].size(); i++){
        v = M->get_vertex(e,i);
        is_bound = M->is_bound(v);

        // checks if v is boundary point //
        if(is_bound){
            bound = M->get_bound(v);

            ////////// Applying below algorithm /////////////
            // suppose B.C. = D @ k
            // bi -= Lik * D        (taking BC from RHS of interior points)
            // Lek = Lke = 0        (setting entire row/col @ k = 0)
            // Lkk = 1.0            (val on diagonal = 1.0)
            // bk D                 (RHS = B.C. @ k)
            for(unsigned int j=0; j<Le[e][i].size(); j++){
                if(i != j){
                    be[e][j] -= Le[e][j][i] * bound;
                    Le[e][i][j] = 0.0;
                    Le[e][j][i] = 0.0;
                }
            }
            Le[e][i][i] = 1.0;
            be[e][i] = bound;
        }
    }
}
////////


////////////// Evaluates area of triangle, given glob coords ////////////////
double FEM::area(double xi[3][3]) const {
    double tmp = 0.0;

    tmp += xi[0][0] * (xi[1][1]*xi[2][2] - xi[1][2]*xi[2][1]);
    tmp -= xi[0][1] * (xi[1][0]*xi[2][2] - xi[1][2]*xi[2][0]);
    tmp += xi[0][2] * (xi[1][0]*xi[2][1] - xi[1][1]*xi[2][0]);

    return 0.5*tmp;
}

// need an updated version of this //
void FEM::output(double *u_an) const {

    output_results(*M, u_an, b, order, 0);
}

double FEM::sse_fem(double *u){
    this->SSE = sse(u,b,order);
    
    return this->SSE;
}