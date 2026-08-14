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

/*
  3D face-mortar kernels.

  Face traces are staged in a buffer with Fortran layout
  buff(1:N+1, 1:N+1, 1:8, 1:nMortars, 1:nl), in the big face's coordinates :
    slots 0..3 - big-face trace (one copy per sub-face so MPI receives never alias)
    slots 4..7 - small-face traces on sub-faces (quadrants) 1..4
  nl counts the (variable, direction) lines; boundary-type arrays are indexed with the
  SCB_3D layout with the line index in place of the variable index, which also covers
  vector data (VEB_3D layout) when nl = 3*nVar.

  mortarInfo has Fortran layout (1:14, 1:nMortars); rows (0-based here) :
    0 - big element (global) ; 1 - big local face
    2+2*q - small element on sub-face q ; 3+2*q - 10*(small local face) + flip
    10..13 - global side ids of the sub-faces (MPI tags; unused in the kernels)

  The flip maps big-face indices to the small face's own indices (MortarFaceMap in
  SELF_Mesh_3D_t); sub-face q = kx + 2*ky covers the big-face half-intervals
  (kx, ky) with kx, ky in {0,1}, and the 1-D operator pair (mortarR/mortarP) is
  applied per direction.

  All element ids are global; elemToRank/offset translate to rank-local addressing,
  exactly as in SideExchange_3D.
*/

#define MORTARBUFF3D_INDEX(i,j,slot,m,l,N,nMortars) i + (N+1)*(j + (N+1)*(slot + 8*(m + nMortars*l)))

// Compile-time bound on N+1 for the shared-memory face staging in MortarFlip_3D;
// enforced at mortar-buffer allocation in the gpu mapped-data classes.
#define MORTAR3D_MAXNP 16

__device__ __forceinline__ void MortarFaceMap_3d(int i, int j, int N, int flip,
                                                 int *i2, int *j2){
  switch(flip){
    case 0: *i2 = i;   *j2 = j;   break;
    case 1: *i2 = N-i; *j2 = j;   break;
    case 2: *i2 = N-i; *j2 = N-j; break;
    case 3: *i2 = i;   *j2 = N-j; break;
    case 4: *i2 = j;   *j2 = i;   break;
    case 5: *i2 = N-j; *j2 = i;   break;
    case 6: *i2 = N-j; *j2 = N-i; break;
    default: *i2 = j;  *j2 = N-i; break; // flip == 7
  }
}

__global__ void MortarGather_3D(real *buff, real *boundary, int *mortarInfo, int *elemToRank,
                                int rankId, int offset, int N, int nMortars, int nEl){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t np = (N+1)*(N+1);
  uint32_t ndof = np*nMortars;
  uint32_t i = idof % (N+1);
  uint32_t j = (idof/(N+1)) % (N+1);
  uint32_t m = idof/np;
  uint32_t l = blockIdx.y;

  if(idof < ndof){

    int eB = mortarInfo[0 + 14*m];
    int sB = mortarInfo[1 + 14*m];

    if(elemToRank[eB-1] == rankId){
      real fb = boundary[SCB_3D_INDEX(i,j,sB-1,eB-1-offset,l,N,nEl)];
      for(int q = 0; q < 4; q++){
        buff[MORTARBUFF3D_INDEX(i,j,q,m,l,N,nMortars)] = fb;
      }
    }

    for(int q = 0; q < 4; q++){
      int eS = mortarInfo[2+2*q + 14*m];
      if(elemToRank[eS-1] == rankId){
        int sS = mortarInfo[3+2*q + 14*m]/10;
        int flip = mortarInfo[3+2*q + 14*m] - 10*sS;
        int i1, j1;
        MortarFaceMap_3d(i,j,N,flip,&i1,&j1);
        buff[MORTARBUFF3D_INDEX(i,j,4+q,m,l,N,nMortars)] =
            boundary[SCB_3D_INDEX(i1,j1,sS-1,eS-1-offset,l,N,nEl)];
      }
    }
  }
}

extern "C"
{
  void MortarGather_3D_gpu(real *buff, real *boundary, int *mortarInfo, int *elemToRank,
                           int rankId, int offset, int N, int nl, int nMortars, int nEl)
  {
    int ndof = (N+1)*(N+1)*nMortars;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    MortarGather_3D<<<dim3(nblocks_x,nl,1), dim3(threads_per_block,1,1), 0, 0>>>(
        buff, boundary, mortarInfo, elemToRank, rankId, offset, N, nMortars, nEl);
  }
}

__global__ void MortarFlip_3D(real *buff, int *mortarInfo, int *elemToRank,
                              int rankId, int N, int nMortars){

  // Reorients small-face traces received over MPI (on the big face's rank) into the
  // big face's coordinates. Locally gathered traces were already reoriented in
  // MortarGather_3D. One block per (mortar, line); the face is staged through shared
  // memory because the general flips include transposes, which cannot be applied
  // in place with independent swap pairs.
  __shared__ real sbuf[(MORTAR3D_MAXNP)*(MORTAR3D_MAXNP)];
  uint32_t i = threadIdx.x % (N+1);
  uint32_t j = threadIdx.x/(N+1);
  uint32_t m = blockIdx.x;
  uint32_t l = blockIdx.y;
  bool active = (threadIdx.x < (N+1)*(N+1));

  int eB = mortarInfo[0 + 14*m];
  bool bigLocal = (elemToRank[eB-1] == rankId);

  for(int q = 0; q < 4; q++){
    int eS = mortarInfo[2+2*q + 14*m];
    int sS = mortarInfo[3+2*q + 14*m]/10;
    int flip = mortarInfo[3+2*q + 14*m] - 10*sS;
    bool reorient = bigLocal && (elemToRank[eS-1] != rankId) && (flip != 0);
    if(reorient && active){
      sbuf[i + (N+1)*j] = buff[MORTARBUFF3D_INDEX(i,j,4+q,m,l,N,nMortars)];
    }
    __syncthreads();
    if(reorient && active){
      int i1, j1;
      MortarFaceMap_3d(i,j,N,flip,&i1,&j1);
      buff[MORTARBUFF3D_INDEX(i,j,4+q,m,l,N,nMortars)] = sbuf[i1 + (N+1)*j1];
    }
    __syncthreads();
  }
}

extern "C"
{
  void MortarFlip_3D_gpu(real *buff, int *mortarInfo, int *elemToRank,
                         int rankId, int N, int nl, int nMortars)
  {
    int threads_per_block = (N+1)*(N+1);
    MortarFlip_3D<<<dim3(nMortars,nl,1), dim3(threads_per_block,1,1), 0, 0>>>(
        buff, mortarInfo, elemToRank, rankId, N, nMortars);
  }
}

__global__ void MortarScatter_3D(real *extBoundary, real *buff, real *mortarR, real *mortarP,
                                 int *mortarInfo, int *elemToRank,
                                 int rankId, int offset, int N, int nMortars, int nEl){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t np = (N+1)*(N+1);
  uint32_t ndof = np*nMortars;
  uint32_t i = idof % (N+1);
  uint32_t j = (idof/(N+1)) % (N+1);
  uint32_t m = idof/np;
  uint32_t l = blockIdx.y;

  if(idof < ndof){

    // Small faces : restrict the big-face trace to each sub-face (exact). The
    // tensor-product restriction is applied as a double sum; kx, ky select the
    // half-interval operator per direction.
    for(int q = 0; q < 4; q++){
      int eS = mortarInfo[2+2*q + 14*m];
      if(elemToRank[eS-1] == rankId){
        int sS = mortarInfo[3+2*q + 14*m]/10;
        int flip = mortarInfo[3+2*q + 14*m] - 10*sS;
        int kx = q % 2;
        int ky = q/2;
        real fm = 0.0;
        for(int jj = 0; jj < N+1; jj++){
          real rowsum = 0.0;
          for(int ii = 0; ii < N+1; ii++){
            rowsum += mortarR[ii + (N+1)*(i + (N+1)*kx)]*
                      buff[MORTARBUFF3D_INDEX(ii,jj,q,m,l,N,nMortars)];
          }
          fm += mortarR[jj + (N+1)*(j + (N+1)*ky)]*rowsum;
        }
        int i1, j1;
        MortarFaceMap_3d(i,j,N,flip,&i1,&j1);
        extBoundary[SCB_3D_INDEX(i1,j1,sS-1,eS-1-offset,l,N,nEl)] = fm;
      }
    }

    // Big face : L2 projection of the small-face traces
    int eB = mortarInfo[0 + 14*m];
    if(elemToRank[eB-1] == rankId){
      int sB = mortarInfo[1 + 14*m];
      real fm = 0.0;
      for(int q = 0; q < 4; q++){
        int kx = q % 2;
        int ky = q/2;
        for(int jj = 0; jj < N+1; jj++){
          real rowsum = 0.0;
          for(int ii = 0; ii < N+1; ii++){
            rowsum += mortarP[ii + (N+1)*(i + (N+1)*kx)]*
                      buff[MORTARBUFF3D_INDEX(ii,jj,4+q,m,l,N,nMortars)];
          }
          fm += mortarP[jj + (N+1)*(j + (N+1)*ky)]*rowsum;
        }
      }
      extBoundary[SCB_3D_INDEX(i,j,sB-1,eB-1-offset,l,N,nEl)] = fm;
    }
  }
}

extern "C"
{
  void MortarScatter_3D_gpu(real *extBoundary, real *buff, real *mortarR, real *mortarP,
                            int *mortarInfo, int *elemToRank,
                            int rankId, int offset, int N, int nl, int nMortars, int nEl)
  {
    int ndof = (N+1)*(N+1)*nMortars;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    MortarScatter_3D<<<dim3(nblocks_x,nl,1), dim3(threads_per_block,1,1), 0, 0>>>(
        extBoundary, buff, mortarR, mortarP, mortarInfo, elemToRank,
        rankId, offset, N, nMortars, nEl);
  }
}

__global__ void MortarFluxScatter_3D(real *boundaryNormal, real *buff, real *mortarP,
                                     int *mortarInfo, int *elemToRank,
                                     int rankId, int offset, int N, int nMortars, int nEl){

  // Replaces the big face's surface-flux integrand with -4 * sum_q (P_kx x P_ky) g_q,
  // where the g_q are the small faces' integrands staged in buffer slots 4..7. The
  // factor of four converts the solution-space projection into the integrand-space
  // projection and the sign accounts for the opposing outward normals (see
  // MortarFluxCollect).
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t np = (N+1)*(N+1);
  uint32_t ndof = np*nMortars;
  uint32_t i = idof % (N+1);
  uint32_t j = (idof/(N+1)) % (N+1);
  uint32_t m = idof/np;
  uint32_t l = blockIdx.y;

  if(idof < ndof){

    int eB = mortarInfo[0 + 14*m];
    if(elemToRank[eB-1] == rankId){
      int sB = mortarInfo[1 + 14*m];
      real fm = 0.0;
      for(int q = 0; q < 4; q++){
        int kx = q % 2;
        int ky = q/2;
        for(int jj = 0; jj < N+1; jj++){
          real rowsum = 0.0;
          for(int ii = 0; ii < N+1; ii++){
            rowsum += mortarP[ii + (N+1)*(i + (N+1)*kx)]*
                      buff[MORTARBUFF3D_INDEX(ii,jj,4+q,m,l,N,nMortars)];
          }
          fm += mortarP[jj + (N+1)*(j + (N+1)*ky)]*rowsum;
        }
      }
      boundaryNormal[SCB_3D_INDEX(i,j,sB-1,eB-1-offset,l,N,nEl)] = -4.0*fm;
    }
  }
}

extern "C"
{
  void MortarFluxScatter_3D_gpu(real *boundaryNormal, real *buff, real *mortarP,
                                int *mortarInfo, int *elemToRank,
                                int rankId, int offset, int N, int nl, int nMortars, int nEl)
  {
    int ndof = (N+1)*(N+1)*nMortars;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block + 1;
    MortarFluxScatter_3D<<<dim3(nblocks_x,nl,1), dim3(threads_per_block,1,1), 0, 0>>>(
        boundaryNormal, buff, mortarP, mortarInfo, elemToRank,
        rankId, offset, N, nMortars, nEl);
  }
}
