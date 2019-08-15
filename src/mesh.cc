#include <set>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "mesh.h"
#include "utils.h"

#include <iostream>

Mesh::Mesh(const int* nr, const float* x, const float *y){
    float dx, dy;
    int count=0;
    int order = (nr[0]+1)*(nr[1]+1);
    int num_cells = 2*nr[0]*nr[1];

    // assigning dimensions of mesh //
    for(int i=0; i<2; i++){
        this->nr[i] = nr[i];
        this->x[i] = x[i];
        this->y[i] = y[i];
    }
    dx = (this->x[1]-this->x[0])/(float)nr[0], dy = (this->y[1]-this->y[0])/(float)nr[1];
    
    /////////////// Allocating memory for mesh //////////////
    assign_ptrs(&vertices, &vert_vals, order, 2);
    assign_ptrs(&cells, &cells_vals, num_cells, 3);
    assign_ptrs(&dof, &dof_vals, num_cells, 3);
    boundary = (int*)calloc(order,sizeof(int));
    bdry_vals = (float*)calloc(order,sizeof(int));
    /////////////////////////////////////////////////////////////
    

    //////// Setting mesh vertices & boundary values ////////////
    for(int i=0; i<=nr[1]; i++){
        for(int j=0; j<=nr[0]; j++){
            // fix this bit //
            vertices[count][1] = y[0] + i*dy;
            vertices[count][0] = x[0] + j*dx;

            if(j==0){
                boundary[count] = true;
                bdry_vals[count] = ui;
            } else if(j==nr[0]){
                boundary[count] = true;
                bdry_vals[count] = uo;
            } else {
                boundary[count] = false;
            }
            count++;
        }
    }
    ///////////////////////////////////////////////////////////////

    std::cout << "Test mesh\n";    
    
    ///////// Setting global vertex values for each cell //////////
    for(int i=0; i<this->nr[1]; i++){
        for(int j=0; j<this->nr[0]; j++){
            cells[2* (j + (i*nr[1]) )][0] = j + i*(nr[0]+1); 
            cells[2* (j + (i*nr[1]) )][1] = j + i*(nr[0]+1) + 1;
            cells[2* (j + (i*nr[1]) )][2] = j + (i+1)*(nr[0]+1); 
            
            cells[2* (j + (i*nr[1]) ) + 1][0] = j + i*(nr[0]+1) + 1; 
            cells[2* (j + (i*nr[1]) ) + 1][2] = j + (i+1)*(nr[0]+1);
            cells[2* (j + (i*nr[1]) ) + 1][1] = j + (i+1)*(nr[0]+1) + 1;

            dof[2* (j + (i*nr[1]) )][0] = j + i*(nr[0]+1); 
            dof[2* (j + (i*nr[1]) )][1] = j + i*(nr[0]+1) + 1;
            dof[2* (j + (i*nr[1]) )][2] = j + (i+1)*(nr[0]+1); 
            
            dof[2* (j + (i*nr[1]) ) + 1][0] = j + i*(nr[0]+1) + 1; 
            dof[2* (j + (i*nr[1]) ) + 1][2] = j + (i+1)*(nr[0]+1);
            dof[2* (j + (i*nr[1]) ) + 1][1] = j + (i+1)*(nr[0]+1) + 1;
        }
    }
    //////////////////////////////////////////////////////////////
    std::cout << "Test mesh\n";    
}

Mesh::~Mesh(){
    delete[] vertices;  delete[] vert_vals;
    delete[] cells;     delete[] cells_vals;
    delete[] dof;       delete[] dof_vals;
    delete[] boundary;  delete[] bdry_vals;
}


//////////// Pass function "map" as parameter to deform mesh from rectangle //////
void Mesh::deform(void (*map)(float*, float*, float*, float, int), float theta){
    float **v = vertices;
    int order = (nr[0]+1)*(nr[1]+1);

    for(int i=0; i<order; i++)
        map(v[i], x, y, theta, 1);
}
/////////


//////////////// Gives x,y global coordinates of node v /////////////////////
void Mesh::get_xy(float *xy, const int v) const {
    xy[0] = vertices[v][0];
    xy[1] = vertices[v][1];
}
//////


////////////// Gives global vertex number of node i in cell e //////////////
int Mesh::get_vertex(const int e, const int i) const {
    int v;
    v = cells[e][i];

    return v;
}


////////////// Gives mapping of node i in cell e to global dof /////////////
// NOTE: For P1 problems this is the same as get_vertex();
int Mesh::dof_map(const int e, const int r) const {
    int dof;
    dof = this->dof[e][r];

    return dof;
}
//////


////////////// Gets boundary value of global node v ////////////////////////
float Mesh::get_bound(const int v) const {
    float boundary;
    boundary = bdry_vals[v];

    return boundary;
}
//////



////////////// Checks if global node v is interior or dirichlet boundary ////////
int Mesh::is_bound(const int v) const {
    int is_bound;
    is_bound = boundary[v];

    return is_bound;
}
//////


/////////////// Gets number of rectangles in x,y direction of mesh ////////////
void Mesh::get_recs(int* nrecs) const {
    nrecs[0] = nr[0];
    nrecs[1] = nr[1];
}
//////


///////////// Returns pointers to all mesh arrays (for passing vals to GPU) ////////////
void Mesh::get_arrays(float **vertices, int **cells, int **dof, int **is_bound, float **bdry_vals){
    *vertices = this->vert_vals;
    *cells = this->cells_vals;
    *dof = this->dof_vals;
    *is_bound = this->boundary;
    *bdry_vals = this->bdry_vals;
}
//////


///////////////////// Returns sparsity pattern of matrix /////////////////////////////////
// Passes through the matrix, ammending any nodes which connect into a vector of sets
// The vector of sets is then converted into the rowPtr and colInd vectors of a CSR matrix
void Mesh::sparsity_pass(
            std::vector<float> &valsL,
            std::vector<int> &rowPtrL,
            std::vector<int> &colIndL,
            int &nnz)
{
    std::vector<std::set<int> > sparsity;       // one set of connected nodes for each node //
    std::vector<int> colTmp;                    // tmp vector for colInd array //
    int v1, v2;
    int n = 0;
    int *rowTmp;
    int count=0;
    int num_cells = nr[0]*nr[1]*2;
    int order = (nr[0]+1)*(nr[1]+1);

    sparsity.resize(order);                     // order = num_nodes        //
    colTmp.resize(order*order, 0);              // max size => 0 non-zeros  //
    rowPtrL.resize(order+1,0);                  // nodes+1 rows 

    rowTmp = &rowPtrL[0];

    /*  Incrementing through each cell 
    adding connected nodes into each set  */
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
    

    /*  Incrementing through vector of sets (connected nodes)
    For each connected node, a column entry is added to that row
    Count is incremented to keep track of the valsL index
    rowPtr is incremented after each set has been iterated
    */ 
    rowTmp++;
    for(std::vector<std::set<int> >::iterator v = sparsity.begin(); v != sparsity.end(); ++v){
        for(std::set<int>::iterator vi = (*v).begin(); vi != (*v).end(); ++vi){
            colTmp[count] = (*vi);
            count++;
        }
        n += (*v).size();
        *(rowTmp++) = n;
    }
    
    //// resized when number non-zeros known ////
    valsL.resize(n,0.0);
    colIndL.resize(n,0);
    
    //// populated tmp column index vector copied to main column vector ////
    std::memcpy(&colIndL[0], &colTmp[0], n*sizeof(int));

    nnz = n;
}
///////


////////////// See above function - Only populates upper half of matrix //////////
// All comments same as above //
void Mesh::sparsity_pass_half(
            std::vector<float> &valsL,
            std::vector<int> &rowPtrL,
            std::vector<int> &colPtrL,
            int &nnz)
{
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
    int row = 0;
    for(std::vector<std::set<int> >::iterator v = sparsity.begin(); v != sparsity.end(); ++v){
        for(std::set<int>::iterator vi = (*v).begin(); vi != (*v).end(); ++vi){
            if((*vi) >= row){           // if upper half of matrix //
                colTmp[count] = (*vi);
                count++;
            } else
                continue;
        }
        row++;
        *(rowTmp++) = count;
    }
    n = count;

    valsL.resize(n,0.0);
    colPtrL.resize(n,0);
    
    std::memcpy(&colPtrL[0], &colTmp[0], n*sizeof(int));

    nnz = n;
}
//////////


/////////////////// Deformation function to pass to M.deform() ////////////////////
// Maps a rectangle to an annulus portion
void annulus_seg_map(float *vertex, float *x, float *y, float theta, int s){
    float x_hat, y_hat;
    float x_hat2;

    x_hat2 = x[0] + (x[1]-x[0]) * pow( (vertex[0]-x[0]) / (x[1]-x[0]), s);
    //y_hat = vertex[1];
    /*
    vertex[0] = x_hat*cos(theta*y_hat);
    vertex[1] = x_hat*sin(theta*y_hat);
    */
    
    x_hat = x_hat2*cos(theta) - vertex[1]*sin(theta);
    y_hat = x_hat2*sin(theta) + vertex[1]*cos(theta);

    vertex[0] = x_hat;
    vertex[1] = y_hat;
}
///////

