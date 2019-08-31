// =========================================================================== //
// Class for completing the standard Finite Element Method on a given          //
// mesh M.                                                                     //
// PDE in question is the Poisson equation defined over the mesh.              //
// Contains vectors to store stiffness matrix in both CSR and dense storage    //
// formats.                                                                    //
// Uses intel MKL library to solve linear system.                               //
// =========================================================================== //

#ifndef _FEM_H_
#define _FEM_H_

class FEM {
private:
    float **L;                      // 2d array for storage of dense stiffness matrix
    float *L_vals;                  // array to store values of stiffness matrix
    float *b;                       // array to store RHS/stress vector
    
    std::vector<float> valsL;       // stores values in CSR format of stiffness matrix
    std::vector<int> rowPtrL;       // storing rowPtrs of stiffness matrix
    std::vector<int> colIndL;       // storing column indices of stiffness matrix

    // Note: reason for array/vector mix was due to 
    // needing contiguous 2d data for cuda transfer

    std::vector<std::vector<std::vector<float> > > Le;  // vector of all element matrices
    std::vector<std::vector<float> > be;                // vector of all element vectors
    
    Mesh* M;                        // mesh passed through in constructor
    int order;                      // order = num_nodes since P1
    int num_cells;
    int nnz;                        // number of non-zeros in CSR of stifness matrix
    float SSE;

public:
    FEM(Mesh &M, Tau &t);
    ~FEM();
    
    void assemble();
    void assemble_csr();
    void solve(Tau &t);
    void MKL_solve();
    void elem_mat(const int e);
    float area(float xi[3][3]) const;
    void output(float *u_an) const;
    float sse_fem(float *u);
};

#endif
