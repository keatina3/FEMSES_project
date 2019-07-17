#ifndef _FEM_H_
#define _FEM_H_

class FEM {
private:
    std::vector<std::vector<float> > L;
    std::vector<float> u;
    std::vector<float> b;
    
    std::vector<std::vector<std::vector<float> > > Le;
    std::vector<std::vector<float> > be;
    
    // instert parameters here, or change default //
    Mesh* M;

public:
    FEM(const int* nr, const int* a, const int* b);
    ~FEM(){ delete M;};
    
    void assemble();
    void solve();
    void elem_mat(const int e);
    float area(float xi[3][3]) const;
};

#endif
