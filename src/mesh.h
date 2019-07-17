#ifndef _MESH_H_
#define _MESH_H_

class Mesh {
private:
    std::vector<std::vector<float> > vertices;
    std::vector<std::vector<int> > cells;
    std::vector<std::vector<int> > dof;
    std::vector<float> boundary;

    int nr[2], a[2], b[2];

public:
    Mesh(const int* nr, const int* a, const int* b);
    ~Mesh(){};

    void deform(void (*map)());
    void get_xy(float *xy, const int v) const;
    int get_vertex(const int e, const int i) const;
    int dof_map(const int e, const int r) const;
    float get_bound(const int v) const;
};

#endif
