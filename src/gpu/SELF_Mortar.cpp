#include "SELF_GPU_Macros.h"

/*
  2:1 nonconforming (mortar) interface kernels.

  Edge traces are staged in a buffer with Fortran layout buff(1:N+1, 1:4, 1:nMortars, 1:nl),
  in the big side's edge orientation :
    slot 0, 1 - big-side trace (one copy per sub-edge so MPI receives never alias)
    slot 2, 3 - small-side traces on sub-edges 1 and 2
  nl counts the (variable, direction) lines; boundary-type arrays are indexed with the
  SCB_2D layout with the line index in place of the variable index, which also covers
  vector data (VEB_2D layout) when nl = 2*nVar.

  mortarInfo has Fortran layout (1:8, 1:nMortars); rows (0-based here) :
    0 - big element (global) ; 1 - big local side
    2 - small element on sub-edge 1 ; 3 - 10*(small local side) + flip
    4 - small element on sub-edge 2 ; 5 - 10*(small local side) + flip
    6, 7 - global side ids of the sub-edges (MPI tags; unused in the kernels)

  All element ids are global; elemToRank/offset translate to rank-local addressing,
  exactly as in SideExchange_2D.
*/

#define MORTARBUFF_INDEX(i,slot,m,l,N,nMortars) i + (N+1)*(slot + 4*(m + nMortars*l))

__global__ void MortarGather_2D(real *buff, real *boundary, int *mortarInfo, int *elemToRank,
                                int rankId, int offset, int N, int nMortars, int nEl){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*nMortars;
  uint32_t i = idof % (N+1);
  uint32_t m = idof/(N+1);
  uint32_t l = blockIdx.y;

  if(idof < ndof){

    int eB = mortarInfo[0 + 8*m];
    int sB = mortarInfo[1 + 8*m];

    if(elemToRank[eB-1] == rankId){
      real fb = boundary[SCB_2D_INDEX(i,sB-1,eB-1-offset,l,N,nEl)];
      buff[MORTARBUFF_INDEX(i,0,m,l,N,nMortars)] = fb;
      buff[MORTARBUFF_INDEX(i,1,m,l,N,nMortars)] = fb;
    }

    for(int k = 0; k < 2; k++){
      int eS = mortarInfo[2+2*k + 8*m];
      if(elemToRank[eS-1] == rankId){
        int sS = mortarInfo[3+2*k + 8*m]/10;
        int flip = mortarInfo[3+2*k + 8*m] - 10*sS;
        int i1 = (flip == 0) ? i : (N-i);
        buff[MORTARBUFF_INDEX(i,2+k,m,l,N,nMortars)] =
            boundary[SCB_2D_INDEX(i1,sS-1,eS-1-offset,l,N,nEl)];
      }
    }
  }
}

extern "C"
{
  void MortarGather_2D_gpu(real *buff, real *boundary, int *mortarInfo, int *elemToRank,
                           int rankId, int offset, int N, int nl, int nMortars, int nEl)
  {
    int ndof = (N+1)*nMortars;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    MortarGather_2D<<<dim3(nblocks_x,nl,1), dim3(threads_per_block,1,1), 0, 0>>>(
        buff, boundary, mortarInfo, elemToRank, rankId, offset, N, nMortars, nEl);
  }
}

__global__ void MortarFlip_2D(real *buff, int *mortarInfo, int *elemToRank,
                              int rankId, int N, int nMortars){

  // Reorients small-side traces received over MPI (on the big side's rank) into the
  // big side's edge orientation. Locally gathered traces were already reoriented in
  // MortarGather_2D. One thread per swap pair.
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t nhalf = (N+1)/2;
  uint32_t ndof = nhalf*nMortars;
  uint32_t i = idof % nhalf;
  uint32_t m = idof/nhalf;
  uint32_t l = blockIdx.y;

  if(idof < ndof){

    int eB = mortarInfo[0 + 8*m];
    if(elemToRank[eB-1] == rankId){
      for(int k = 0; k < 2; k++){
        int eS = mortarInfo[2+2*k + 8*m];
        int sS = mortarInfo[3+2*k + 8*m]/10;
        int flip = mortarInfo[3+2*k + 8*m] - 10*sS;
        if(elemToRank[eS-1] != rankId && flip == 1){
          real tmp = buff[MORTARBUFF_INDEX(i,2+k,m,l,N,nMortars)];
          buff[MORTARBUFF_INDEX(i,2+k,m,l,N,nMortars)] =
              buff[MORTARBUFF_INDEX(N-i,2+k,m,l,N,nMortars)];
          buff[MORTARBUFF_INDEX(N-i,2+k,m,l,N,nMortars)] = tmp;
        }
      }
    }
  }
}

extern "C"
{
  void MortarFlip_2D_gpu(real *buff, int *mortarInfo, int *elemToRank,
                         int rankId, int N, int nl, int nMortars)
  {
    int ndof = ((N+1)/2)*nMortars;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    MortarFlip_2D<<<dim3(nblocks_x,nl,1), dim3(threads_per_block,1,1), 0, 0>>>(
        buff, mortarInfo, elemToRank, rankId, N, nMortars);
  }
}

__global__ void MortarScatter_2D(real *extBoundary, real *buff, real *mortarR, real *mortarP,
                                 int *mortarInfo, int *elemToRank,
                                 int rankId, int offset, int N, int nMortars, int nEl){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*nMortars;
  uint32_t i = idof % (N+1);
  uint32_t m = idof/(N+1);
  uint32_t l = blockIdx.y;

  if(idof < ndof){

    // Small sides : restrict the big-side trace to each sub-edge (exact)
    for(int k = 0; k < 2; k++){
      int eS = mortarInfo[2+2*k + 8*m];
      if(elemToRank[eS-1] == rankId){
        int sS = mortarInfo[3+2*k + 8*m]/10;
        int flip = mortarInfo[3+2*k + 8*m] - 10*sS;
        real fm = 0.0;
        for(int ii = 0; ii < N+1; ii++){
          fm += mortarR[ii + (N+1)*(i + (N+1)*k)]*
                buff[MORTARBUFF_INDEX(ii,k,m,l,N,nMortars)];
        }
        int iout = (flip == 0) ? i : (N-i);
        extBoundary[SCB_2D_INDEX(iout,sS-1,eS-1-offset,l,N,nEl)] = fm;
      }
    }

    // Big side : L2 projection of the small-side traces
    int eB = mortarInfo[0 + 8*m];
    if(elemToRank[eB-1] == rankId){
      int sB = mortarInfo[1 + 8*m];
      real fm = 0.0;
      for(int k = 0; k < 2; k++){
        for(int ii = 0; ii < N+1; ii++){
          fm += mortarP[ii + (N+1)*(i + (N+1)*k)]*
                buff[MORTARBUFF_INDEX(ii,2+k,m,l,N,nMortars)];
        }
      }
      extBoundary[SCB_2D_INDEX(i,sB-1,eB-1-offset,l,N,nEl)] = fm;
    }
  }
}

extern "C"
{
  void MortarScatter_2D_gpu(real *extBoundary, real *buff, real *mortarR, real *mortarP,
                            int *mortarInfo, int *elemToRank,
                            int rankId, int offset, int N, int nl, int nMortars, int nEl)
  {
    int ndof = (N+1)*nMortars;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    MortarScatter_2D<<<dim3(nblocks_x,nl,1), dim3(threads_per_block,1,1), 0, 0>>>(
        extBoundary, buff, mortarR, mortarP, mortarInfo, elemToRank,
        rankId, offset, N, nMortars, nEl);
  }
}

__global__ void MortarFluxScatter_2D(real *boundaryNormal, real *buff, real *mortarP,
                                     int *mortarInfo, int *elemToRank,
                                     int rankId, int offset, int N, int nMortars, int nEl){

  // Replaces the big side's surface-flux integrand with -2 * sum_k P_k g_k, where the
  // g_k are the small sides' integrands staged in buffer slots 2 and 3. The factor of
  // two converts the solution-space projection into the integrand-space projection and
  // the sign accounts for the opposing outward normals (see MortarFluxCollect).
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*nMortars;
  uint32_t i = idof % (N+1);
  uint32_t m = idof/(N+1);
  uint32_t l = blockIdx.y;

  if(idof < ndof){

    int eB = mortarInfo[0 + 8*m];
    if(elemToRank[eB-1] == rankId){
      int sB = mortarInfo[1 + 8*m];
      real fm = 0.0;
      for(int k = 0; k < 2; k++){
        for(int ii = 0; ii < N+1; ii++){
          fm += mortarP[ii + (N+1)*(i + (N+1)*k)]*
                buff[MORTARBUFF_INDEX(ii,2+k,m,l,N,nMortars)];
        }
      }
      boundaryNormal[SCB_2D_INDEX(i,sB-1,eB-1-offset,l,N,nEl)] = -2.0*fm;
    }
  }
}

extern "C"
{
  void MortarFluxScatter_2D_gpu(real *boundaryNormal, real *buff, real *mortarP,
                                int *mortarInfo, int *elemToRank,
                                int rankId, int offset, int N, int nl, int nMortars, int nEl)
  {
    int ndof = (N+1)*nMortars;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    MortarFluxScatter_2D<<<dim3(nblocks_x,nl,1), dim3(threads_per_block,1,1), 0, 0>>>(
        boundaryNormal, buff, mortarP, mortarInfo, elemToRank,
        rankId, offset, N, nMortars, nEl);
  }
}
