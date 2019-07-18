#ifndef _FEM_H_
#define _FEM_H_

#define ERR 1.0E-08

class FEM {
private:
    float **L;
    float *L_vals;
    float *b;
    
    std::vector<std::vector<std::vector<float> > > Le;
    std::vector<std::vector<float> > be;
    
    // instert parameters here, or change default //
    Mesh* M;
    int order;
    int num_cells;

public:
    FEM(Mesh *M);
    ~FEM();
    
    void assemble();
    void solve();
    void elem_mat(const int e);
    float area(float xi[3][3]) const;
    void output(char* fname) const;
};

#endif
