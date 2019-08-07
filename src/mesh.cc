#include <set>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include "mesh.h"
#include "utils.h"

// put in BC parameters //
Mesh::Mesh(const int* nr, const float* a, const float *b){
    // allow for more DOFs
    //int dof = 3;
    int dim = 2;
    int dofs = 3;
    float dx, dy;
    int count=0;
    int order = (nr[0]+1)*(nr[1]+1);
    int num_cells = 2*nr[0]*nr[1];

    for(int i=0; i<2; i++){
        this->nr[i] = nr[i];
        this->a[i] = a[i];
        this->b[i] = b[i];
    }

    dx = (this->b[0]-this->a[0])/nr[0], dy = (this->b[1]-this->a[1])/nr[1];
    //std::cout << dx << std::endl;
    
    /*
    vertices.resize((nr[0]+1)*(nr[1]+1), std::vector<float>(dim, 0.0));
    cells.resize(2*nr[0]*nr[1], std::vector<int>(3,0));
    dof.resize(2*nr[0]*nr[1], std::vector<int>(3,0));
    boundary.resize((nr[0]+1)*(nr[1]+1), false);
    bdr_val.resize((nr[0]+1)*(nr[1]+1), 0.0);
    */
    std::cout << "test test 1\n";

    assign_ptrs(&vertices, &vert_vals, order, dim);
    assign_ptrs(&cells, &cells_vals, num_cells, dofs);
    assign_ptrs(&dof, &dof_vals, num_cells, dofs);
    boundary = (int*)calloc(order,sizeof(int));
    bdry_vals = (float*)calloc(order,sizeof(int));
    
    std::cout << "test test 2\n";
    
    // (x1,y1) - (x2,y2) //
    //
    for(int i=0; i<=nr[1]; i++){
        for(int j=0; j<=nr[0]; j++){
            // fix this bit //
            vertices[count][1] = a[1] + i*dy;
            vertices[count][0] = a[0] + j*dx;

            if(j==0){
                boundary[count] = true;
                bdry_vals[count] = 2.0;
            } else if(j==nr[0]){
                boundary[count] = true;
                bdry_vals[count] = 6.0;
            } else {
                boundary[count] = false;
            }
            count++;
        }
    }
    
    std::cout << "testing \n";

    for(int i=0; i<this->nr[1]; i++){
        for(int j=0; j<this->nr[0]; j++){
            cells[2* (j + (i*nr[1]) ) ][0] = j + i*(nr[0]+1); 
            cells[2* (j + (i*nr[1]) )][1] = j + i*(nr[0]+1) + 1;
            cells[2* (j + (i*nr[1]) )][2] = j + (i+1)*(nr[0]+1); 
            
            cells[2* (j + (i*nr[1]) ) + 1][0] = j + i*(nr[0]+1) + 1; 
            cells[2* (j + (i*nr[1]) ) + 1][2] = j + (i+1)*(nr[0]+1);
            cells[2* (j + (i*nr[1]) ) + 1][1] = j + (i+1)*(nr[0]+1) + 1;

            dof[2* (j + (i*nr[1]) ) ][0] = j + i*(nr[0]+1); 
            dof[2* (j + (i*nr[1]) )][1] = j + i*(nr[0]+1) + 1;
            dof[2* (j + (i*nr[1]) )][2] = j + (i+1)*(nr[0]+1); 
            
            dof[2* (j + (i*nr[1]) ) + 1][0] = j + i*(nr[0]+1) + 1; 
            dof[2* (j + (i*nr[1]) ) + 1][2] = j + (i+1)*(nr[0]+1);
            dof[2* (j + (i*nr[1]) ) + 1][1] = j + (i+1)*(nr[0]+1) + 1;
        }
    }
    
    std::cout << "testing \n";
}

Mesh::~Mesh(){
    delete[] vertices; delete[] vert_vals;
    delete[] cells; delete[] cells_vals;
    delete[] dof; delete[] dof_vals;
    delete[] boundary; delete[] bdry_vals;
}

void Mesh::deform(void (*map)(float*, float*, float*, float, int)){
    float **v = vertices;
    int order = (nr[0]+1)*(nr[1]+1);

    for(int i=0; i<order; i++)
        map(v[i], a, b, M_PI/2.0, 2);

    std::cout << "testing \n";
    /*
    // check this works properly //
    for(std::vector<std::vector<float> >::iterator v=vertices.begin(); v!=vertices.end(); ++v){ 
        map(*v, a, b, M_PI/2.0, 2);
    }
    */
}

void Mesh::get_xy(float *xy, const int v) const {
    xy[0] = vertices[v][0];
    xy[1] = vertices[v][1];
}

int Mesh::get_vertex(const int e, const int i) const {
    int v;
    v = cells[e][i];

    return v;
}

int Mesh::dof_map(const int e, const int r) const {
    int dof;
    dof = this->dof[e][r];

    return dof;
}

float Mesh::get_bound(const int v) const {
    float boundary;
    boundary = bdry_vals[v];

    return boundary;
}

int Mesh::is_bound(const int v) const {
    int is_bound;
    is_bound = boundary[v];

    return is_bound;
}

// expand this to be N-dimensional //
void Mesh::get_recs(int* nrecs) const {
    nrecs[0] = nr[0];
    nrecs[1] = nr[1];
}

void Mesh::get_arrays(float **vertices, int **cells, int **dof, int **is_bound, float **bdry_vals){
    *vertices = this->vert_vals;
    *cells = this->cells_vals;
    *dof = this->dof_vals;
    *is_bound = this->boundary;
    *bdry_vals = this->bdry_vals;
}

int Mesh::sparsity_pass(std::vector<float> &valsL, std::vector<int> &rowPtrL, std::vector<int> &colPtrL){
    std::vector<std::set<int> > sparsity;
    std::vector<int> colTmp;
    int v1, v2;
    int n = 0;
    int *rowTmp;
    int count=0;
    int num_cells = nr[0]*nr[1]*2;
    int order = (nr[0]+1)*(nr[1]+1);

    sparsity.resize(order);
    colTmp.resize(order*order, 0);
    rowPtrL.resize(order+1,0);

    rowTmp = &rowPtrL[0];

    for(int e=0; e<num_cells; e++){
        for(int i=0; i<3; i++){
            v1 = get_vertex(e, i);
            for(int j=0; j<3; j++){
                v2 = get_vertex(e, j);
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
    
    valsL.resize(n,0.0);
    colPtrL.resize(n,0);
    
    std::memcpy(&colPtrL[0], &colTmp[0], n*sizeof(int));
    std::cout << "sparsity test\n";    

    return n;
}

void annulus_seg_map(float *vertex, float *a, float *b, float theta, int s){
    float x_hat, y_hat;
    
    x_hat = a[0] + (b[0]-a[0]) * pow( (vertex[0]-a[0]) / (b[0]-a[0]), s);
    y_hat = vertex[1];

    vertex[0] = x_hat*cos(theta*y_hat);
    vertex[1] = x_hat*sin(theta*y_hat);
}
