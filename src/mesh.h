#ifndef _MESH_H_
#define _MESH_H_

class Mesh {
private:
    std::vector<std::vector<float> > vertices;
    std::vector<std::vector<int> > cells;
    std::vector<std::vector<int> > dof;
    std::vector<bool> boundary;
    std::vector<float> bdr_val;

    int nr[2];
    float a[2], b[2];

public:
    Mesh(const int* nr, const float* a, const float* b);
    ~Mesh(){};

    void deform(void (*map)(std::vector<float>&, float*, float*, float, int));
    void get_xy(float *xy, const int v) const;
    int get_vertex(const int e, const int i) const;
    int dof_map(const int e, const int r) const;
    float get_bound(const int v) const;
    bool is_bound(const int v) const;
    void get_recs(int* nrecs) const;
};

void annulus_seg_map(std::vector<float> &vertex, float *a, float *b, float theta, int s);

#endif
