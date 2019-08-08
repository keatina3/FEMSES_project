#ifndef _MESH_H_
#define _MESH_H_

class Mesh {
private:

    float **vertices, *vert_vals;
    int **cells, *cells_vals;
    int **dof, *dof_vals;
    int *boundary;
    float *bdry_vals;
    
    /*
    std::vector<std::vector<float> > vertices;
    std::vector<std::vector<int> > cells;
    std::vector<std::vector<int> > dof;
    std::vector<int> boundary;
    std::vector<float> bdr_val;
    */

    int nr[2];
    float x[2], y[2];

public:
    Mesh(const int* nr, const float* a, const float* b);
    ~Mesh();

    void deform(void (*map)(float*, float*, float*, float, int));
    void get_xy(float *xy, const int v) const;
    int get_vertex(const int e, const int i) const;
    int dof_map(const int e, const int r) const;
    float get_bound(const int v) const;
    int is_bound(const int v) const;
    void get_recs(int* nrecs) const;
    void get_arrays(float **vertices, int **cells, int **dof, int **is_bound, float **bdry_vals);
    int sparsity_pass(std::vector<float> &valsL, std::vector<int> &rowPtr, 
                        std::vector<int> &colPtrL);
};

void annulus_seg_map(float *vertex, float *a, float *b, float theta, int s);

#endif
