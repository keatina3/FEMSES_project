#include <vector>
#include <cmath>
#include "mesh.h"

/*
Mesh::Mesh(){
     
}
*/ 

Mesh::Mesh(const int* nr, const int* a, const int *b){
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

    dx = this->b[0]-this->a[0], dy = this->b[1]-this->a[1];
    
    vertices.resize((nr[0]+1)*(nr[1]+1), std::vector<float>(dim, 0.0));
    cells.resize(2*nr[0]*nr[1], std::vector<int>(3,0));
    dof.resize(2*nr[0]*nr[1], std::vector<int>(3,0));
    boundary.resize((nr[0]+1)*(nr[1]+1), 0.0);

    // CHECK IF a[1] = y-axis or x-axis ??? //
    for(int i=this->a[1]; i<=this->b[1]; i+=dy){
        for(int j=this->a[0]; j<=this->b[0]; j+=dx){
            vertices[count][1] = i;
            vertices[count][0] = j;

            // DO MORE GENERAL BCs IMPLEMENTATION //
            if(j==a[0]){
                boundary[count] = a[0];
            } else if(j==b[0]){
                boundary[count] = b[0];
            } else {
                boundary[count] = 0.0;
            }
            count++;
        }
    }
     
    for(int i=0; i<this->nr[1]; i++){
        for(int j=0; j<this->nr[0]; j++){
            cells[(2*j) + (i*nr[1])][0] = j + i*(nr[0]+1); 
            cells[(2*j) + (i*nr[1])][1] = j + i*(nr[0]+1) + 1;
            cells[(2*j) + (i*nr[1])][2] = j + (i+1)*(nr[0]+1); 
            
            cells[(2*j) + (i*nr[1]) + 1][0] = j + i*(nr[0]+1) + 1; 
            cells[(2*j) + (i*nr[1]) + 1][1] = j + (i+1)*(nr[0]+1);
            cells[(2*j) + (i*nr[1]) + 1][2] = j + (i+1)*(nr[0]+1) + 1;

            dof[(2*j) + (i*nr[1])][0] = j + i*(nr[0]+1); 
            dof[(2*j) + (i*nr[1])][1] = j + i*(nr[0]+1) + 1;
            dof[(2*j) + (i*nr[1])][2] = j + (i+1)*(nr[0]+1); 
            
            dof[(2*j) + (i*nr[1]) + 1][0] = j + i*(nr[0]+1) + 1; 
            dof[(2*j) + (i*nr[1]) + 1][1] = j + (i+1)*(nr[0]+1);
            dof[(2*j) + (i*nr[1]) + 1][2] = j + (i+1)*(nr[0]+1) + 1;
        }
    }

    // deform();
}

void Mesh::deform(void (*map)()){
    
    // check this works properly //
    for(std::vector<std::vector<float> >::iterator v=vertices.begin(); v!=vertices.end(); ++v){ 
        //map(v);
    }
}

void Mesh::get_xy(float *xy, const int v) const{
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
    boundary = this->boundary[v];

    return boundary;
}

void annulus_portion(std::vector<int> &vertex, int *a, int *b, float theta){
    float x_hat, y_hat;
    
    x_hat = a[0] + (b[0]-a[0]) * pow( (vertex[0]-a[0]) / (b[0]-a[0]), 2);
    y_hat = vertex[1];

    vertex[0] = x_hat*cos(theta*y_hat/2);
    vertex[1] = x_hat*sin(theta*y_hat/2);
}
