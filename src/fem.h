#ifndef _FEM_H_
#define _FEM_H_

#define ERR 1.0E-08

class FEM {
private:
    //float **L;
    //float *L_vals;
    float *b;
    
    std::vector<float> valsL;
    std::vector<int> rowPtrL;
    std::vector<int> colPtrL;

    std::vector<std::vector<std::vector<float> > > Le;
    std::vector<std::vector<float> > be;
    
    // instert parameters here, or change default //
    Mesh* M;
    int order;
    int num_cells;
    int nnz;

public:
    FEM(Mesh &M);
    ~FEM();
    
    void assemble();
    void assemble_csr();
    void solve();
    void MKL_solve();
    void elem_mat(const int e);
    float area(float xi[3][3]) const;
    void output(char* fname) const;
    int sparsity_pass();
};

#endif
