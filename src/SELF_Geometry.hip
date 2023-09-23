#include <hip/hip_runtime.h>
#include <math.h>
#include "SELF_HIP_Macros.h"

__global__ void CalculateContravariantBasis_SEMQuad_gpu(real *dxds, real *dsdx, int N){

  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    // Ja1
    dsdx[TE_2D_INDEX(1,1,i,j,0,iEl,N,1)] = dxds[TE_2D_INDEX(2,2,i,j,0,iEl,N,1)];
    dsdx[TE_2D_INDEX(2,1,i,j,0,iEl,N,1)] = -dxds[TE_2D_INDEX(1,2,i,j,0,iEl,N,1)];

    // Ja2
    dsdx[TE_2D_INDEX(1,2,i,j,0,iEl,N,1)] = -dxds[TE_2D_INDEX(2,1,i,j,0,iEl,N,1)];
    dsdx[TE_2D_INDEX(2,2,i,j,0,iEl,N,1)] = dxds[TE_2D_INDEX(1,1,i,j,0,iEl,N,1)];

}
extern "C"
{
  void CalculateContravariantBasis_SEMQuad_gpu_wrapper(real **dxds, real **dsdx, int N, int nEl)
  { 
	  CalculateContravariantBasis_SEMQuad_gpu<<<dim3(nEl,1,1), dim3(N+1,N+1), 0, 0>>>(*dxds, *dsdx, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

__global__ void AdjustBoundaryContravariantBasis_SEMQuad_gpu(real *dsdx, real *J, int N){

  size_t iSide = blockIdx.x+1;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

    real fac = fabs(J[SCB_2D_INDEX(i,iSide,0,iEl,N,1)])/J[SCB_2D_INDEX(i,iSide,0,iEl,N,1)];
    if (iSide == 4 || iSide == 1 ){
      fac = -fac;
    } 

    for( int row = 1; row <= 2; row++ ){
      for( int col = 1; col <= 2; col++ ){
        dsdx[TEB_2D_INDEX(row,col,i,0,iSide,iEl,N,1)] = fac*dsdx[TEB_2D_INDEX(row,col,i,0,iSide,iEl,N,1)]; 
      }
    }

}

extern "C"
{
  void AdjustBoundaryContravariantBasis_SEMQuad_gpu_wrapper(real **dsdx, real **J, int N, int nEl)
  { 
	  AdjustBoundaryContravariantBasis_SEMQuad_gpu<<<dim3(4,nEl,1), dim3(N+1,1,1), 0, 0>>>(*dsdx, *J, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}
__global__ void CalculateContravariantBasis_SEMHex_gpu(real *dxds, real *dsdx, int N){

  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    // Ja1
    dsdx[TE_3D_INDEX(1,1,i,j,k,0,iEl,N,1)] = 
	      dxds[TE_3D_INDEX(2,2,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(3,3,i,j,k,0,iEl,N,1)]- 
              dxds[TE_3D_INDEX(3,2,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(2,3,i,j,k,0,iEl,N,1)];

    dsdx[TE_3D_INDEX(2,1,i,j,k,0,iEl,N,1)] =
              dxds[TE_3D_INDEX(1,3,i,j,k,0,iEl,N,1)]*
              dxds[TE_3D_INDEX(3,2,i,j,k,0,iEl,N,1)]-
              dxds[TE_3D_INDEX(3,3,i,j,k,0,iEl,N,1)]*
              dxds[TE_3D_INDEX(1,2,i,j,k,0,iEl,N,1)];

    dsdx[TE_3D_INDEX(3,1,i,j,k,0,iEl,N,1)] = 
              dxds[TE_3D_INDEX(1,2,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(2,3,i,j,k,0,iEl,N,1)]- 
              dxds[TE_3D_INDEX(2,2,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(1,3,i,j,k,0,iEl,N,1)];

    // Ja2
    dsdx[TE_3D_INDEX(1,2,i,j,k,0,iEl,N,1)] = 
              dxds[TE_3D_INDEX(2,3,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(3,1,i,j,k,0,iEl,N,1)]- 
              dxds[TE_3D_INDEX(3,3,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(2,1,i,j,k,0,iEl,N,1)];

    dsdx[TE_3D_INDEX(2,2,i,j,k,0,iEl,N,1)] = 
              dxds[TE_3D_INDEX(1,1,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(3,3,i,j,k,0,iEl,N,1)]- 
              dxds[TE_3D_INDEX(3,1,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(1,3,i,j,k,0,iEl,N,1)];

    dsdx[TE_3D_INDEX(3,2,i,j,k,0,iEl,N,1)] = 
              dxds[TE_3D_INDEX(1,3,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(2,1,i,j,k,0,iEl,N,1)]- 
              dxds[TE_3D_INDEX(2,3,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(1,1,i,j,k,0,iEl,N,1)];

    // Ja3
    dsdx[TE_3D_INDEX(1,3,i,j,k,0,iEl,N,1)] = 
              dxds[TE_3D_INDEX(2,1,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(3,2,i,j,k,0,iEl,N,1)]- 
              dxds[TE_3D_INDEX(3,1,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(2,2,i,j,k,0,iEl,N,1)];

    dsdx[TE_3D_INDEX(2,3,i,j,k,0,iEl,N,1)] = 
              dxds[TE_3D_INDEX(1,2,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(3,1,i,j,k,0,iEl,N,1)]- 
              dxds[TE_3D_INDEX(3,2,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(1,1,i,j,k,0,iEl,N,1)];

    dsdx[TE_3D_INDEX(3,3,i,j,k,0,iEl,N,1)] = 
              dxds[TE_3D_INDEX(1,1,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(2,2,i,j,k,0,iEl,N,1)]- 
              dxds[TE_3D_INDEX(2,1,i,j,k,0,iEl,N,1)]* 
              dxds[TE_3D_INDEX(1,2,i,j,k,0,iEl,N,1)];
}
extern "C"
{
  void CalculateContravariantBasis_SEMHex_gpu_wrapper(real **dxds, real **dsdx, int N, int nEl)
  { 
	  CalculateContravariantBasis_SEMHex_gpu<<<dim3(nEl,1,1), dim3(N+1,N+1,N+1), 0, 0>>>(*dxds, *dsdx, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

__global__ void AdjustBoundaryContravariantBasis_SEMHex_gpu(real *dsdx, real *J, int N){

  size_t iSide = blockIdx.x+1;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    real fac = fabs(J[SCB_3D_INDEX(i,j,iSide,0,iEl,N,1)])/J[SCB_3D_INDEX(i,j,iSide,0,iEl,N,1)];
    if (iSide == 5 || iSide == 1 || iSide == 2){
      fac = -fac;
    } 

    for( int row = 1; row <= 3; row++ ){
      for( int col = 1; col <= 3; col++ ){
        dsdx[TEB_3D_INDEX(row,col,i,j,0,iSide,iEl,N,1)] = fac*dsdx[TEB_3D_INDEX(row,col,i,j,0,iSide,iEl,N,1)]; 
      }
    }

}

extern "C"
{
  void AdjustBoundaryContravariantBasis_SEMHex_gpu_wrapper(real **dsdx, real **J, int N, int nEl)
  { 
	  AdjustBoundaryContravariantBasis_SEMHex_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*dsdx, *J, N);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}
