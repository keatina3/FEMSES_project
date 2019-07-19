#include "mesh.h"
#include "gpu_fem.h"

__device__ float area(float *xi){
    float tmp = 0.0;

    tmp += xi[0]*(xi[4]*xi[8] - xi[5]*xi[7]);
    tmp -= xi[1]*(xi[3]*xi[8] - xi[5]*xi[6]);
    tmp += xi[2]*(xi[3]*xi[7] - xi[4]*xi[6]);

    return tmp;
}

__device__ void elem_mat_gpu(float *vertices, float *cells,  float *is_bound, float *bdry_vals, float *temp1, int idx, int idy){
    float *Le, *be, *xi, *consts;
    int v;
    float bound;

    Le = temp1;
    be = &temp1[9];
    xi = &temp1[12]
    consts = &temp1[21]
    
    v = cells[(idx*3) + idy];

    xi[3*idy] = 1.0;
    xi[(3*idy) + 1] = vertices[2*v]; 
    xi[(3*idy) + 2] = vertices[(2*v)+1];

    if(idy==0)
        consts[9] = idx%2 ? area(xi) : (-1)*area(xi);
    __syncthreads();

    // consts[(3*idy)] = xi[(idy+1)%3][1] * xi[(i+2)%3][2] - xi[(i+2)%3][1] * xi[(i+1)%3][2];
    consts[(3*idy)+1] = xi[ 3*((idy+1)%3) +2] - xi[ 3*((idy+2)%3) + 2]
    consts[(3*idy)+2] = xi[ 3*((idy+1)%3) +1] - xi[ 3*((idy+2)%3) + 1]
    
    be[idy] = 0.0;      // or Int(fv) //
     
    for(int i=0; i<idy; i++){
        Le[(3*idy)+i] = 0.25*consts[9]*consts[9]*consts[9] * (consts[(3*idy)+1]*consts[(3*i)+1] 
                                    + consts[(3*idy)+2]*consts[(3*i)+2]);
        Le[(3*i)+idy] = Le[(3*idy)+i];
    }

    if(is_bound[v]){
        bound = bdry_vals[v];
        // change this appropriately for more DOF if necessary //
        for(int j=0; j<3; j++){
            if(idy != j){
                be[j] -= Le[(3*j) + idy]*bound;
                Le[(3*j) + idy] = 0.0;
                Le[(3*idy) + j] = 0.0;
            }
        }
        Le[(3*idy)+idy] = 1.0;
        be[idy] = bound;
    }                            
}

__global__ void assemble_gpu(int num_cells, float *Le, float *be, float *vertices, float *cells, float *is_bound, float *bdry_vals){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;
    extern __shared__ float temp1;

    if(idx < num_cells && idy < 3){

    }
} 

extern void gpu_fem(Mesh &M){
    int nr[2];
    int num_nodes, dim;
    int block_size_X, block_size_Y;
    float *vertices_gpu, *vertices;
    float *cells_gpu, *cells;
    float *dof_gpu, *dof;
    bool *is_bound_gpu, *is_bound;
    float *bdry_vals_gpu, *bdry_vals;
    float *L, *Le *b, *be;

    M->get_recs(nr);

    num_nodes = (nr[0]+1)*(nr[1]+1);
    dim = 2+1;      // needs expansion here //
    order = num_nodes + 0;
    num_cells = 2*nr[0]*nr[1];
    
    // Sorting Mesh // 
    get_arrays(&vertices, &cells, &dof, &is_bound, &bdry_vals);

    cudaMalloc( (void**)&vertices_gpu, 2*num_nodes*sizeof(float));
    cudaMalloc( (void**)&cells_gpu, num_cells*sizeof(float));
    cudaMalloc( (void**)&dof_gpu, num_cells*sizeof(float));
    cudaMalloc( (void**)&is_bound_gpu, num_nodes*sizeof(bool));
    cudaMalloc( (void**)&bdry_vals_gpu, num_nodes*sizeof(float));
    
    cudaMemcpy(vertices_gpu, vertices, 2*num_nodes(sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(cells_gpu, cells, num_cells*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dof_gpu, dof, num_cells*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(is_bound_gpu, is_bound, num_nodes*sizeof(bool), cudaMemcpyHostToDevice);
    cudaMemcpy(bdry_vals_gpu, bdry_vals, num_nodes*sizeof(float), cudaMemcpyHostToDevice);
    //////////////////

    cudaMalloc( (void**)&L, order*order*sizeof(float));
    cudaMalloc( (void**)&b, order*sizeof(float));
    //cudaMalloc(b);
    //cudaMalloc(be);
    
    block_size_X = 1, block_size_Y = 3;
    dim3 dimBlock(block_size_X, block_size_Y);
    dim3 dimGrid((n/dimBlock.x)+(!(n%dimBlock.x)?0:1),
                (numSamples/dimBlock.y)+(!(numSamples%dimBlock.y)?0:1));
    
    // this assumes 1 dof per triangle //
    assemble_gpu<<<dimGrid, dimBlock, 31*sizeof(float)>>>(Le, be, vertices_gpu, cells_gpu, is_bound_gpu, bdry_vals_gpu);

    // solve<<<
}
