// =========================================================================== //
// Class for handling the mesh invoked in the FEM.                             //
// Contains arrays storing the global x,y coordinates of each node,            //
// global node index and dirichlet boundary conditions.                        //
// =========================================================================== //

#ifndef _MESH_H_
#define _MESH_H_

class Mesh {
private:

    float **vertices, *vert_vals;       // glob x,y coordinates
    int **cells, *cells_vals;           // glob node numbering
    int **dof, *dof_vals;               // same as above for P1
    int *boundary;                      // bool for boundary/interior node
    float *bdry_vals;                   // boundary value if true
    
    int nr[2];                          // # rectangles in mesh in each axis
    float x[2], y[2];                   // x0,x1,y0,y1 

public:
    Mesh(const int* nr, const float* a, const float* b);
    ~Mesh();

    void deform(void (*map)(float*, float*, float*, float, int), float theta);
    void get_xy(float *xy, const int v) const;
    int get_vertex(const int e, const int i) const;
    int dof_map(const int e, const int r) const;
    float get_bound(const int v) const;
    int is_bound(const int v) const;
    void get_recs(int* nrecs) const;
    void get_arrays(float **vertices, int **cells, int **dof, int **is_bound, float **bdry_vals);
    void sparsity_pass(std::vector<float> &valsL, std::vector<int> &rowPtr, 
                        std::vector<int> &colPtrL, int &nnz, float &alloc, float &tau);
    void sparsity_pass_half(std::vector<float> &valsL, std::vector<int> &rowPtrL,
                        std::vector<int> &colPtrL, int &nnz, float &alloc, float &tau);
};

void annulus_seg_map(float *vertex, float *a, float *b, float theta, int s);

#endif
