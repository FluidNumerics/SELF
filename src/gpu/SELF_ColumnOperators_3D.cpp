#include "SELF_GPU_Macros.h"

// Device kernels of the column operators on extruded meshes (SELF_ColumnOperators_3D_t).
// Element e3 = e2*nLayers + k (0-based) is layer k, counted from the bottom, of column e2.
// A volume field of nVar variables is indexed i + (N+1)*(j + (N+1)*(l + (N+1)*(e3 + nEl*ivar)))
// and a surface field i + (N+1)*(j + (N+1)*(e2 + nColumns*ivar)), as SC_3D_INDEX and SC_2D_INDEX.
// One thread owns one vertical grid line (i,j) of one column and one variable, and visits
// the layers of the column in order. Consecutive threads take consecutive i, so every load
// and store of a warp is contiguous.

__global__ void ColumnIntegral_3D_gpukernel(real *q, real *dxds, real *qWeights, real *qbar,
                                            int N, int nLayers, int nColumns){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ivar = blockIdx.y;
  uint32_t nq = (N+1)*(N+1);
  uint32_t ndof = nq*nColumns;

  if( idof < ndof ){
    uint32_t ij = idof % nq;
    uint32_t e2 = idof/nq;
    uint32_t nEl = nColumns*nLayers;
    uint32_t ndof3 = nq*(N+1)*nEl;
    real s = 0.0;
    for(int k = 0; k < nLayers; k++){
      uint32_t iel = e2*nLayers + k;
      for(int l = 0; l < N+1; l++){
        uint32_t node = ij + nq*(l + (N+1)*iel);
        // dxds(i,j,l,iel,1,3,3) is Z_zeta; TE_3D_INDEX with nVar = 1, row = col = 2
        s += qWeights[l]*dxds[node + 8*ndof3]*q[node + ndof3*ivar];
      }
    }
    qbar[idof + ndof*ivar] = s;
  }

}

__global__ void ColumnBroadcast_3D_gpukernel(real *qbar, real *q, int N, int nLayers, int nColumns){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ivar = blockIdx.y;
  uint32_t nq = (N+1)*(N+1);
  uint32_t ndof3 = nq*(N+1)*nColumns*nLayers;

  if( idof < ndof3 ){
    uint32_t ij = idof % nq;
    uint32_t iel = idof/(nq*(N+1));
    uint32_t e2 = iel/nLayers;
    q[idof + ndof3*ivar] = qbar[ij + nq*(e2 + nColumns*ivar)];
  }

}

__global__ void VerticalSweep_3D_gpukernel(real *aInv, real *r, real *inflow, real *q,
                                           int N, int nLayers, int nColumns, int upward){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ivar = blockIdx.y;
  uint32_t nq = (N+1)*(N+1);
  uint32_t ndof = nq*nColumns;

  if( idof < ndof ){
    uint32_t ij = idof % nq;
    uint32_t e2 = idof/nq;
    uint32_t nEl = nColumns*nLayers;
    uint32_t ndof3 = nq*(N+1)*nEl;
    // The node of each layer that passes the inflow on: the top node going up, the bottom one going down
    int lout = upward ? N : 0;
    real qin = inflow[idof + ndof*ivar];
    for(int kk = 0; kk < nLayers; kk++){
      int k = upward ? kk : nLayers-1-kk;
      uint32_t iel = e2*nLayers + k;
      real qout = qin;
      for(int l = 0; l < N+1; l++){
        // q_l = q_in + sum_m aInv(m,l) r_m, aInv in Fortran (column major) order
        real s = 0.0;
        for(int m = 0; m < N+1; m++){
          s += aInv[m + (N+1)*l]*r[ij + nq*(m + (N+1)*iel) + ndof3*ivar];
        }
        real ql = qin + s;
        q[ij + nq*(l + (N+1)*iel) + ndof3*ivar] = ql;
        if( l == lout ) qout = ql;
      }
      qin = qout;
    }
  }

}

extern "C"
{
  void ColumnIntegral_3D_gpu(real *q, real *dxds, real *qWeights, real *qbar, int N, int nLayers, int nColumns, int nvar)
  {
    int ndof = (N+1)*(N+1)*nColumns;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nvar,1);
    dim3 nthreads(threads_per_block,1,1);

    ColumnIntegral_3D_gpukernel<<<nblocks, nthreads, 0, 0>>>(q, dxds, qWeights, qbar, N, nLayers, nColumns);
  }
}

extern "C"
{
  void ColumnBroadcast_3D_gpu(real *qbar, real *q, int N, int nLayers, int nColumns, int nvar)
  {
    int ndof3 = (N+1)*(N+1)*(N+1)*nColumns*nLayers;
    int threads_per_block = 256;
    int nblocks_x = ndof3/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nvar,1);
    dim3 nthreads(threads_per_block,1,1);

    ColumnBroadcast_3D_gpukernel<<<nblocks, nthreads, 0, 0>>>(qbar, q, N, nLayers, nColumns);
  }
}

extern "C"
{
  void VerticalSweep_3D_gpu(real *aInv, real *r, real *inflow, real *q, int N, int nLayers, int nColumns, int nvar, int upward)
  {
    int ndof = (N+1)*(N+1)*nColumns;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;

    dim3 nblocks(nblocks_x,nvar,1);
    dim3 nthreads(threads_per_block,1,1);

    VerticalSweep_3D_gpukernel<<<nblocks, nthreads, 0, 0>>>(aInv, r, inflow, q, N, nLayers, nColumns, upward);
  }
}
