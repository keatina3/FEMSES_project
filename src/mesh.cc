#include <iostream>
#include <vector>
#include <cmath>
#include "mesh.h"

/*
Mesh::Mesh(){
     
}
*/ 
// put in BC parameters //
Mesh::Mesh(const int* nr, const float* a, const float *b){
    // allow for more DOFs
    //int dof = 3;
    int dim = 2;
    float dx, dy;
    int count=0;
    
    for(int i=0; i<2; i++){
        this->nr[i] = nr[i];
        this->a[i] = a[i];
        this->b[i] = b[i];
    }

    dx = (this->b[0]-this->a[0])/nr[0], dy = (this->b[1]-this->a[1])/nr[1];
    //std::cout << dx << std::endl;

    vertices.resize((nr[0]+1)*(nr[1]+1), std::vector<float>(dim, 0.0));
    cells.resize(2*nr[0]*nr[1], std::vector<int>(3,0));
    dof.resize(2*nr[0]*nr[1], std::vector<int>(3,0));
    boundary.resize((nr[0]+1)*(nr[1]+1), false);
    bdr_val.resize((nr[0]+1)*(nr[1]+1), 0.0);
    
    // (x1,y1) - (x2,y2) //
    //
    for(int i=0; i<=nr[1]; i++){
        for(int j=0; j<=nr[0]; j++){
            // fix this bit //
            vertices[count][1] = a[1] + i*dy;
            vertices[count][0] = a[0] + j*dx;

            if(j==0){
                boundary[count] = true;
                bdr_val[count] = 2.0;
            } else if(j==nr[0]){
                boundary[count] = true;
                bdr_val[count] = 6.0;
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
}

void Mesh::deform(void (*map)(std::vector<float>&, float*, float*, float, int)){
    
    // check this works properly //
    for(std::vector<std::vector<float> >::iterator v=vertices.begin(); v!=vertices.end(); ++v){ 
        map(*v, a, b, M_PI/2.0, 2);
    }
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
    boundary = bdr_val[v];

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
    *vertices = &this->vertices[0][0];
    *cells = &this->cells[0][0];
    *dof = &this->dof[0][0];
    *is_bound = &this->boundary[0];
    *bdry_vals = &this->bdr_val[0];
}

void annulus_seg_map(std::vector<float> &vertex, float *a, float *b, float theta, int s){
    float x_hat, y_hat;
    
    x_hat = a[0] + (b[0]-a[0]) * pow( (vertex[0]-a[0]) / (b[0]-a[0]), s);
    y_hat = vertex[1];

    vertex[0] = x_hat*cos(theta*y_hat);
    vertex[1] = x_hat*sin(theta*y_hat);
}
