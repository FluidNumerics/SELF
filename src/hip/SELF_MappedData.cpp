#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

/*
// Template
__global__ void Template_{1D|2D|3D}_gpu( , int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;


  // How to access scalars 
  f[SC_1D_INDEX(i,iVar,iEl,N,nVar)];
  f[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];
  f[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

  // How to access vectors (2d)
  f[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)];
  f[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)];


  // How to access vectors (3d)
  f[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)];
  f[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)];
  f[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)];

  // How to access tensors (2d)
  f[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)];
  f[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)];
  f[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)];
  f[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)];


  // How to access tensors (3d)
  f[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)];
  f[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,N,nVar)];

}

extern "C"
{
  void Template_{1D|2D|3D}_gpu_wrapper(int N, int nVar)
  {
	  hipLaunchKernelGGL((Template_{1D|2D|3D}_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0,  N, M, nVar);
  } 
}

*/

// JacobianWeight_MappedScalar1D_gpu
__global__ void JacobianWeight_MappedScalar1D_gpu(real *scalar, real *dxds, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;

    scalar[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = scalar[SC_1D_INDEX(i,iVar,iEl,N,nVar)]/
                                             dxds[SC_1D_INDEX(i,0,iEl,N,1)];
}

extern "C"
{
  void JacobianWeight_MappedScalar1D_gpu_wrapper(real **scalar, real **dxds, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedScalar1D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0>>>(*scalar, *dxds, N, nVar);
  }
}

// JacobianWeight_MappedScalar2D_gpu
__global__ void JacobianWeight_MappedScalar2D_gpu(real *scalar, real *jacobian, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]/
                                               jacobian[SC_2D_INDEX(i,j,0,iEl,N,1)];
}

extern "C"
{
  void JacobianWeight_MappedScalar2D_gpu_wrapper(real **scalar, real **jacobian, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedScalar2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*scalar, *jacobian, N, nVar);
  }
}

// JacobianWeight_MappedScalar3D_gpu
__global__ void JacobianWeight_MappedScalar3D_gpu(real *scalar, real *jacobian, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)] = scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];
}

extern "C"
{
  void JacobianWeight_MappedScalar3D_gpu_wrapper(real **scalar, real **jacobian, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedScalar3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*scalar, *jacobian, N, nVar);
  }
}

// ContravariantWeight_MappedScalar2D_gpu
__global__ void ContravariantWeight_MappedScalar2D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    workTensor[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(1,1,i,j,0,iEl,N,1)]*
                                                       scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]; 
  
    workTensor[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(1,2,i,j,0,iEl,N,1)]*
                                                       scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]; 
  
    workTensor[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(2,1,i,j,0,iEl,N,1)]*
                                                       scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]; 
  
    workTensor[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(2,2,i,j,0,iEl,N,1)]*
                                                       scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]; 
}

extern "C"
{
  void ContravariantWeight_MappedScalar2D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantWeight_MappedScalar2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*scalar, *workTensor, *dsdx, N, nVar);
  } 
}

// ContravariantWeightBoundary_MappedScalar2D_gpu
__global__ void ContravariantWeightBoundary_MappedScalar2D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nVar){

  size_t iSide = blockIdx.x+1;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t iVar = threadIdx.y;

    workTensor[TEB_2D_INDEX(1,1,i,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_2D_INDEX(1,1,i,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)]; 
  
    workTensor[TEB_2D_INDEX(2,1,i,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_2D_INDEX(1,2,i,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)]; 
  
    workTensor[TEB_2D_INDEX(1,2,i,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_2D_INDEX(2,1,i,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)]; 
  
    workTensor[TEB_2D_INDEX(2,2,i,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_2D_INDEX(2,2,i,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)]; 
}

extern "C"
{
  void ContravariantWeightBoundary_MappedScalar2D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantWeightBoundary_MappedScalar2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*scalar, *workTensor, *dsdx, N, nVar);
  } 
}

// ContravariantWeight_MappedScalar3D_gpu
__global__ void ContravariantWeight_MappedScalar3D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  workTensor[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,1,i,j,k,0,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,2,i,j,k,0,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,3,i,j,k,0,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(2,1,i,j,k,0,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(2,2,i,j,k,0,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(2,3,i,j,k,0,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(3,1,i,j,k,0,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(3,2,i,j,k,0,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(3,3,i,j,k,0,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

}

extern "C"
{
  void ContravariantWeight_MappedScalar3D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantWeight_MappedScalar3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*scalar, *workTensor, *dsdx, N, nVar);
  } 
}

// ContravariantWeightBoundary_MappedScalar3D_gpu
__global__ void ContravariantWeightBoundary_MappedScalar3D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nVar){

  size_t iSide = blockIdx.x+1;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t iVar = threadIdx.z;

  workTensor[TEB_3D_INDEX(1,1,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(1,1,i,j,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(2,1,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(1,2,i,j,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(3,1,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(1,3,i,j,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(1,2,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(2,1,i,j,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(2,2,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(2,2,i,j,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(3,2,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(2,3,i,j,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(1,3,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(3,1,i,j,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(2,3,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(3,2,i,j,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(3,3,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(3,3,i,j,0,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

}

extern "C"
{
  void ContravariantWeightBoundary_MappedScalar3D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantWeightBoundary_MappedScalar3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nVar), 0, 0>>>(*scalar, *workTensor, *dsdx, N, nVar);
  } 
}

// ContravariantProjection_MappedVector2D_gpu
__global__ void ContravariantProjection_MappedVector2D_gpu(real *physVector, real *compVector, real *dsdx, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;


    compVector[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(1,1,i,j,0,iEl,N,1)]*
                                                     physVector[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)]+ 
						     dsdx[TE_2D_INDEX(2,1,i,j,0,iEl,N,1)]*
                                                     physVector[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)];

    compVector[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(1,2,i,j,0,iEl,N,1)]*
                                                     physVector[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)]+ 
						     dsdx[TE_2D_INDEX(2,2,i,j,0,iEl,N,1)]*
                                                     physVector[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)];
  
}

extern "C"
{
  void ContravariantProjection_MappedVector2D_gpu_wrapper(real **physVector, real **compVector, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantProjection_MappedVector2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*physVector, *compVector, *dsdx, N, nVar);
  } 
}

__global__ void ContravariantProjectionBoundary_MappedVector2D_gpu(real *physVector, real *boundaryNormal, real *nhat, int N, int nVar){

  size_t iSide = blockIdx.x+1;
  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.z;
  size_t i = threadIdx.x;

    boundaryNormal[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)] = nhat[VEB_2D_INDEX(1,i,0,iSide,iEl,N,1)]*
                                                     physVector[VEB_2D_INDEX(1,i,iVar,iSide,iEl,N,nVar)]+ 
						     nhat[VEB_2D_INDEX(1,i,0,iSide,iEl,N,1)]*
                                                     physVector[VEB_2D_INDEX(2,i,iVar,iSide,iEl,N,nVar)];
}

extern "C"
{
  void ContravariantProjectionBoundary_MappedVector2D_gpu_wrapper(real **physVector, real **boundaryNormal, real **nhat, int N, int nVar, int nEl)
  {
    ContravariantProjectionBoundary_MappedVector2D_gpu<<<dim3(4,nVar,nEl), dim3(N+1,1,1), 0, 0>>>(*physVector, *boundaryNormal, *nhat, N, nVar);
  } 
}

// JacobianWeight_MappedVector2D_gpu
__global__ void JacobianWeight_MappedVector2D_gpu(real *vector, real *jacobian, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    vector[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] = vector[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_2D_INDEX(i,j,0,iEl,N,1)];

    vector[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] = vector[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_2D_INDEX(i,j,0,iEl,N,1)];
}

extern "C"
{
  void JacobianWeight_MappedVector2D_gpu_wrapper(real **vector, real **jacobian, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedVector2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*vector, *jacobian, N, nVar);
  }
}

// ContravariantProjection_MappedVector3D_gpu
__global__ void ContravariantProjection_MappedVector3D_gpu(real *physVector, real *compVector, real *dsdx, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  compVector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,1,i,j,k,0,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)]+ 
                                                     dsdx[TE_3D_INDEX(2,1,i,j,k,0,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)]+
                                                     dsdx[TE_3D_INDEX(3,1,i,j,k,0,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)];

  compVector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,2,i,j,k,0,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)]+ 
                                                     dsdx[TE_3D_INDEX(2,2,i,j,k,0,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)]+
                                                     dsdx[TE_3D_INDEX(3,2,i,j,k,0,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)];

  compVector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,3,i,j,k,0,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)]+ 
                                                     dsdx[TE_3D_INDEX(2,3,i,j,k,0,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)]+
                                                     dsdx[TE_3D_INDEX(3,3,i,j,k,0,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)];

}

extern "C"
{
  void ContravariantProjection_MappedVector3D_gpu_wrapper(real **physVector, real **compVector, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantProjection_MappedVector3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*physVector, *compVector, *dsdx, N, nVar);
  } 
}

__global__ void ContravariantProjectionBoundary_MappedVector3D_gpu(real *physVector, real *boundaryNormal, real *nHat, int N, int nVar){

  size_t iSide = blockIdx.x+1;
  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.z;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  boundaryNormal[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)] = nHat[VEB_3D_INDEX(1,i,j,0,iSide,iEl,N,1)]*
                                                     physVector[VEB_3D_INDEX(1,i,j,iVar,iSide,iEl,N,nVar)]+ 
                                                     nHat[VEB_3D_INDEX(2,i,j,0,iSide,iEl,N,1)]*
                                                     physVector[VEB_3D_INDEX(2,i,j,iVar,iSide,iEl,N,nVar)]+
                                                     nHat[VEB_3D_INDEX(3,i,j,0,iSide,iEl,N,1)]*
                                                     physVector[VEB_3D_INDEX(3,i,j,iVar,iSide,iEl,N,nVar)];

}

extern "C"
{
  void ContravariantProjectionBoundary_MappedVector3D_gpu_wrapper(real **physVector, real **boundaryNormal, real **nHat, int N, int nVar, int nEl)
  {
    ContravariantProjectionBoundary_MappedVector3D_gpu<<<dim3(6,nVar,nEl), dim3(N+1,N+1,1), 0, 0>>>(*physVector, *boundaryNormal, *nHat, N, nVar);
  } 
}

// JacobianWeight_MappedVector3D_gpu
__global__ void JacobianWeight_MappedVector3D_gpu(real *vector, real *jacobian, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    vector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];

    vector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];

    vector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];
}

extern "C"
{
  void JacobianWeight_MappedVector3D_gpu_wrapper(real **vector, real **jacobian, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedVector3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*vector, *jacobian, N, nVar);
  }
}

// JacobianWeight_MappedTensor2D_gpu
__global__ void JacobianWeight_MappedTensor2D_gpu(real *tensor, real *jacobian, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    tensor[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)] = tensor[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_2D_INDEX(i,j,0,iEl,N,1)];

    tensor[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)] = tensor[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_2D_INDEX(i,j,0,iEl,N,1)];

    tensor[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)] = tensor[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_2D_INDEX(i,j,0,iEl,N,1)];

    tensor[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)] = tensor[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_2D_INDEX(i,j,0,iEl,N,1)];
}

extern "C"
{
  void JacobianWeight_MappedTensor2D_gpu_wrapper(real **tensor, real **jacobian, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedTensor2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*tensor, *jacobian, N, nVar);
  }
}

// JacobianWeight_MappedTensor3D_gpu
__global__ void JacobianWeight_MappedTensor3D_gpu(real *tensor, real *jacobian, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    tensor[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];

    tensor[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];

    tensor[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];

    tensor[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];

    tensor[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];

    tensor[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];

    tensor[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];

    tensor[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];

    tensor[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];
}

extern "C"
{
  void JacobianWeight_MappedTensor3D_gpu_wrapper(real **tensor, real **jacobian, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedTensor3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*tensor, *jacobian, N, nVar);
  }
}

// CalculateCurl_MappedTensor2D_gpu
__global__ void CalculateCurl_MappedTensor2D_gpu(real *dfdx, real *curlf, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    curlf[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = dfdx[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)]-
                                              dfdx[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)];
}

extern "C"
{
  void CalculateCurl_MappedTensor2D_gpu_wrapper(real **dfdx, real **curlf, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedTensor2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*dfdx, *curlf, N, nVar);
  }
}

// CalculateCurl_MappedTensor3D_gpu
__global__ void CalculateCurl_MappedTensor3D_gpu(real *dfdx, real *curlf, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    curlf[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] = dfdx[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)]-
                                                  dfdx[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)];

    curlf[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] = dfdx[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)]-
                                                  dfdx[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,N,nVar)];

    curlf[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] = dfdx[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,N,nVar)]-
                                                  dfdx[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)];
}

extern "C"
{
  void CalculateCurl_MappedTensor3D_gpu_wrapper(real **dfdx, real **curlf, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedTensor3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*dfdx, *curlf, N, nVar);
  }
}

// MapTo support routines
__global__ void MapToScalar_MappedVector2D_gpu(real* scalar, real* vector, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    size_t jVar = 2*(iVar);
    scalar[SC_2D_INDEX(i,j,jVar,iEl,N,nVar*2)] = vector[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)];

    jVar += 1;
    scalar[SC_2D_INDEX(i,j,jVar,iEl,N,nVar*2)] = vector[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)];

}

extern "C"
{
  void MapToScalar_MappedVector2D_gpu_wrapper(real **scalar, real **vector, int N, int nVar, int nEl)
  {
    MapToScalar_MappedVector2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*scalar, *vector, N, nVar);
  }
}

__global__ void MapToScalarBoundary_MappedVector2D_gpu(real* scalar, real* vector, int N, int nVar){

  size_t iSide = blockIdx.x+1;
  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.z;
  size_t j = threadIdx.x;

    size_t jVar = 2*(iVar);
    scalar[SCB_2D_INDEX(j,jVar,iSide,iEl,N,nVar*2)] = vector[VEB_2D_INDEX(1,j,iVar,iSide,iEl,N,nVar)];

    jVar += 1;
    scalar[SCB_2D_INDEX(j,jVar,iSide,iEl,N,nVar*2)] = vector[VEB_2D_INDEX(2,j,iVar,iSide,iEl,N,nVar)];

}

extern "C"
{
  void MapToScalarBoundary_MappedVector2D_gpu_wrapper(real **scalar, real **vector, int N, int nVar, int nEl)
  {
    MapToScalarBoundary_MappedVector2D_gpu<<<dim3(4,nVar,nEl), dim3(N+1,1,1), 0, 0>>>(*scalar, *vector, N, nVar);
  }
}

__global__ void MapToTensor_MappedVector2D_gpu(real* tensor, real* vector, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    size_t jVar = 2*(iVar);
    tensor[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)] = vector[VE_2D_INDEX(1,i,j,jVar,iEl,N,nVar*2)];

    tensor[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)] = vector[VE_2D_INDEX(2,i,j,jVar,iEl,N,nVar*2)];

    jVar +=1 ;
    tensor[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)] = vector[VE_2D_INDEX(1,i,j,jVar,iEl,N,nVar*2)];

    tensor[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)] = vector[VE_2D_INDEX(2,i,j,jVar,iEl,N,nVar*2)];

}

extern "C"
{
  void MapToTensor_MappedVector2D_gpu_wrapper(real **tensor, real **vector, int N, int nVar, int nEl)
  {
    MapToTensor_MappedVector2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*tensor, *vector, N, nVar);
  }
}

__global__ void MapToTensorBoundary_MappedVector2D_gpu(real* tensor, real* vector, int N, int nVar){

  size_t iSide = blockIdx.x+1;
  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.z;
  size_t j = threadIdx.x;

    size_t jVar =2*(iVar);
    tensor[TEB_2D_INDEX(1,1,j,iVar,iSide,iEl,N,nVar)] = vector[VEB_2D_INDEX(1,j,jVar,iSide,iEl,N,nVar*2)];

    tensor[TEB_2D_INDEX(1,2,j,iVar,iSide,iEl,N,nVar)] = vector[VEB_2D_INDEX(2,j,jVar,iSide,iEl,N,nVar*2)];

    jVar += 1;
    tensor[TEB_2D_INDEX(2,1,j,iVar,iSide,iEl,N,nVar)] = vector[VEB_2D_INDEX(1,j,jVar,iSide,iEl,N,nVar*2)];

    tensor[TEB_2D_INDEX(2,2,j,iVar,iSide,iEl,N,nVar)] = vector[VEB_2D_INDEX(2,j,jVar,iSide,iEl,N,nVar*2)];

}

extern "C"
{
  void MapToTensorBoundary_MappedVector2D_gpu_wrapper(real **tensor, real **vector, int N, int nVar, int nEl)
  {
    MapToTensorBoundary_MappedVector2D_gpu<<<dim3(4,nVar,nEl), dim3(N+1,1,1), 0, 0>>>(*tensor, *vector, N, nVar);
  }
}

__global__ void MapToScalar_MappedVector3D_gpu(real* scalar, real* vector, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    size_t jVar = 3*(iVar);
    scalar[SC_3D_INDEX(i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)];

    jVar += 1;
    scalar[SC_3D_INDEX(i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)];

    jVar += 1;
    scalar[SC_3D_INDEX(i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)];

}

extern "C"
{
  void MapToScalar_MappedVector3D_gpu_wrapper(real **scalar, real **vector, int N, int nVar, int nEl)
  {
    MapToScalar_MappedVector3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*scalar, *vector, N, nVar);
  }
}

__global__ void MapToScalarBoundary_MappedVector3D_gpu(real* scalar, real* vector, int N, int nVar){

  size_t iSide = blockIdx.x+1;
  size_t iEl = blockIdx.y;
  size_t j = threadIdx.x;
  size_t k = threadIdx.y;
  size_t iVar = threadIdx.z;

    size_t jVar = 3*(iVar);
    scalar[SCB_3D_INDEX(j,k,jVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(1,j,k,iVar,iSide,iEl,N,nVar)];

    jVar += 1;
    scalar[SCB_3D_INDEX(j,k,jVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(2,j,k,iVar,iSide,iEl,N,nVar)];

    jVar += 1;
    scalar[SCB_3D_INDEX(j,k,jVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(3,j,k,iVar,iSide,iEl,N,nVar)];

}

extern "C"
{
  void MapToScalarBoundary_MappedVector3D_gpu_wrapper(real **scalar, real **vector, int N, int nVar, int nEl)
  {
    MapToScalarBoundary_MappedVector3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nVar), 0, 0>>>(*scalar, *vector, N, nVar);
  }
}

__global__ void MapToTensor_MappedVector3D_gpu(real* tensor, real* vector, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    size_t jVar = 3*(iVar);
    tensor[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(1,i,j,k,jVar,iEl,N,nVar*3)];

    tensor[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(2,i,j,k,jVar,iEl,N,nVar*3)];

    tensor[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(3,i,j,k,jVar,iEl,N,nVar*3)];

    jVar += 1;
    tensor[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(1,i,j,k,jVar,iEl,N,nVar*3)];

    tensor[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(2,i,j,k,jVar,iEl,N,nVar*3)];

    tensor[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(3,i,j,k,jVar,iEl,N,nVar*3)];

    jVar += 1;
    tensor[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(1,i,j,k,jVar,iEl,N,nVar*3)];

    tensor[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(2,i,j,k,jVar,iEl,N,nVar*3)];

    tensor[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(3,i,j,k,jVar,iEl,N,nVar*3)];

}

extern "C"
{
  void MapToTensor_MappedVector3D_gpu_wrapper(real **tensor, real **vector, int N, int nVar, int nEl)
  {
    MapToTensor_MappedVector3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*tensor, *vector, N, nVar);
  }
}

__global__ void MapToTensorBoundary_MappedVector3D_gpu(real* tensor, real* vector, int N, int nVar){

  size_t iSide = blockIdx.x+1;
  size_t iEl = blockIdx.y;
  size_t j = threadIdx.x;
  size_t k = threadIdx.x;
  size_t iVar = threadIdx.z;

    size_t jVar = 3*(iVar);
    tensor[TEB_3D_INDEX(1,1,j,k,iVar,iSide,iEl,N,nVar)] = vector[VEB_3D_INDEX(1,j,k,jVar,iSide,iEl,N,nVar*3)];

    tensor[TEB_3D_INDEX(1,2,j,k,iVar,iSide,iEl,N,nVar)] = vector[VEB_3D_INDEX(2,j,k,jVar,iSide,iEl,N,nVar*3)];

    tensor[TEB_3D_INDEX(1,3,j,k,iVar,iSide,iEl,N,nVar)] = vector[VEB_3D_INDEX(3,j,k,jVar,iSide,iEl,N,nVar*3)];

    jVar += 1;
    tensor[TEB_3D_INDEX(2,1,j,k,iVar,iSide,iEl,N,nVar)] = vector[VEB_3D_INDEX(1,j,k,jVar,iSide,iEl,N,nVar*3)];

    tensor[TEB_3D_INDEX(2,2,j,k,iVar,iSide,iEl,N,nVar)] = vector[VEB_3D_INDEX(2,j,k,jVar,iSide,iEl,N,nVar*3)];

    tensor[TEB_3D_INDEX(2,3,j,k,iVar,iSide,iEl,N,nVar)] = vector[VEB_3D_INDEX(3,j,k,jVar,iSide,iEl,N,nVar*3)];

    jVar += 1;
    tensor[TEB_3D_INDEX(3,1,j,k,iVar,iSide,iEl,N,nVar)] = vector[VEB_3D_INDEX(1,j,k,jVar,iSide,iEl,N,nVar*3)];

    tensor[TEB_3D_INDEX(3,2,j,k,iVar,iSide,iEl,N,nVar)] = vector[VEB_3D_INDEX(2,j,k,jVar,iSide,iEl,N,nVar*3)];

    tensor[TEB_3D_INDEX(3,3,j,k,iVar,iSide,iEl,N,nVar)] = vector[VEB_3D_INDEX(3,j,k,jVar,iSide,iEl,N,nVar*3)];

}

extern "C"
{
  void MapToTensorBoundary_MappedVector3D_gpu_wrapper(real **tensor, real **vector, int N, int nVar, int nEl)
  {
    MapToTensorBoundary_MappedVector3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nVar), 0, 0>>>(*tensor, *vector, N, nVar);
  }
}

__global__ void SideExchange_MappedScalar2D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t ivar = threadIdx.y;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)]-1;
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank == rankId ){
      if(flip == 0){
        extBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nVar)] = boundary[SCB_2D_INDEX(i1,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 1){
        int i2 = N-i1;
        extBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nVar)] = boundary[SCB_2D_INDEX(i2,ivar,s2,e2,N,nVar)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedScalar2D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    SideExchange_MappedScalar2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, N, nVar);
  }

}

__global__ void SideExchange_MappedVector2D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t ivar = threadIdx.y;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)]-1;
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];
  int i2 = N-i1;

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank == rankId ){
      if(flip == 0){
        extBoundary[VEB_2D_INDEX(1,i1,ivar,s1,e1,N,nVar)] = boundary[VEB_2D_INDEX(1,i1,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_2D_INDEX(2,i1,ivar,s1,e1,N,nVar)] = boundary[VEB_2D_INDEX(2,i1,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 1){
        extBoundary[VEB_2D_INDEX(1,i1,ivar,s1,e1,N,nVar)] = boundary[VEB_2D_INDEX(1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_2D_INDEX(2,i1,ivar,s1,e1,N,nVar)] = boundary[VEB_2D_INDEX(2,i2,ivar,s2,e2,N,nVar)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedVector2D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    SideExchange_MappedVector2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, N, nVar);
  }

}

__global__ void SideExchange_MappedTensor2D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t s1 = blockIdx.x;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t ivar = threadIdx.y;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)]-1;
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];
  int i2 = N-i1;

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank == rankId ){
      if(flip == 0){
        extBoundary[TEB_2D_INDEX(1,1,i1,ivar,s1,e1,N,nVar)] = boundary[TEB_2D_INDEX(1,1,i1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_2D_INDEX(2,1,i1,ivar,s1,e1,N,nVar)] = boundary[TEB_2D_INDEX(2,1,i1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_2D_INDEX(1,2,i1,ivar,s1,e1,N,nVar)] = boundary[TEB_2D_INDEX(1,2,i1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_2D_INDEX(2,2,i1,ivar,s1,e1,N,nVar)] = boundary[TEB_2D_INDEX(2,2,i1,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 1){
        extBoundary[TEB_2D_INDEX(1,1,i1,ivar,s1,e1,N,nVar)] = boundary[TEB_2D_INDEX(1,1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_2D_INDEX(2,1,i1,ivar,s1,e1,N,nVar)] = boundary[TEB_2D_INDEX(2,1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_2D_INDEX(1,2,i1,ivar,s1,e1,N,nVar)] = boundary[TEB_2D_INDEX(1,2,i2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_2D_INDEX(2,2,i1,ivar,s1,e1,N,nVar)] = boundary[TEB_2D_INDEX(2,2,i2,ivar,s2,e2,N,nVar)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedTensor2D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    SideExchange_MappedTensor2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, N, nVar);
  }

}

__global__ void SideExchange_MappedScalar3D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  size_t ivar = threadIdx.z;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)]-1;
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];
  int i2 = N-i1;
  int j2 = N-j1;

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank == rankId ){
      if(flip == 0){
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = boundary[SCB_3D_INDEX(i1,j1,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 1){
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = boundary[SCB_3D_INDEX(j2,i1,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 2){
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = boundary[SCB_3D_INDEX(i2,j2,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 3){
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = boundary[SCB_3D_INDEX(j1,i2,ivar,s2,e2,N,nVar)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedScalar3D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    SideExchange_MappedScalar3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nVar), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, N, nVar);
  }

}

__global__ void SideExchange_MappedVector3D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  size_t ivar = threadIdx.z;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)]-1;
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];
  int i2 = N-i1;
  int j2 = N-j1;

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank == rankId ){
      if(flip == 0){
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(1,i1,j1,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(2,i1,j1,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(3,i1,j1,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 1){
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(1,j2,i1,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(2,j2,i1,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(3,j2,i1,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 2){
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(1,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(2,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(3,i2,j2,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 3){
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(1,j1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(2,j1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(3,j1,i2,ivar,s2,e2,N,nVar)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedVector3D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    SideExchange_MappedVector3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nVar), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, N, nVar);
  }

}

__global__ void SideExchange_MappedTensor3D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  size_t ivar = threadIdx.z;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)]-1;
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];
  int i2 = N-i1;
  int j2 = N-j1;

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank == rankId ){
      if(flip == 0){
        extBoundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nVar)];
        extBoundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nVar)];
        extBoundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nVar)];
        extBoundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 1){
        extBoundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(1,1,j2,i1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(2,1,j2,i1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(3,1,j2,i1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(1,2,j2,i1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(2,2,j2,i1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(3,2,j2,i1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(1,3,j2,i1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(2,3,j2,i1,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(3,3,j2,i1,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 2){
        extBoundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(1,1,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(2,1,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(3,1,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(1,2,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(2,2,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(3,2,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(1,3,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(2,3,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(3,3,i2,j2,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 3){
        extBoundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(1,1,j1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(2,1,j1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(3,1,j1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(1,2,j1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(2,2,j1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(3,2,j1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(1,3,j1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(2,3,j1,i2,ivar,s2,e2,N,nVar)];
        extBoundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[TEB_3D_INDEX(3,3,j1,i2,ivar,s2,e2,N,nVar)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedTensor3D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    SideExchange_MappedTensor3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nVar), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, N, nVar);
  }

}

__global__ void BassiRebaySides_MappedScalar2D_gpu(real* avgBoundary, real *boundary, real *extBoundary, int N, int nVar){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t ivar = threadIdx.y;
  
  avgBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nVar)] =0.5*(extBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nVar)]+
		                                     boundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nVar)]);
  
}

extern "C"
{
  void BassiRebaySides_MappedScalar2D_gpu_wrapper(real **avgBoundary, real **boundary, real **extBoundary, int N, int nVar, int nEl)
  {
    BassiRebaySides_MappedScalar2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*avgBoundary, *boundary, *extBoundary, N, nVar);
  }

}

__global__ void BassiRebaySides_MappedVector2D_gpu(real *extBoundary, real *boundary, int N, int nVar){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t ivar = threadIdx.y;
  
  boundary[VEB_2D_INDEX(1,i1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[VEB_2D_INDEX(1,i1,ivar,s1,e1,N,nVar)]+
                                                   boundary[VEB_2D_INDEX(1,i1,ivar,s1,e1,N,nVar)]);
  boundary[VEB_2D_INDEX(2,i1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[VEB_2D_INDEX(2,i1,ivar,s1,e1,N,nVar)]+
                                                   boundary[VEB_2D_INDEX(2,i1,ivar,s1,e1,N,nVar)]);
  
}

extern "C"
{
  void BassiRebaySides_MappedVector2D_gpu_wrapper(real **extBoundary, real **boundary, int N, int nVar, int nEl)
  {
    BassiRebaySides_MappedVector2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*extBoundary, *boundary, N, nVar);
  }

}

__global__ void BassiRebaySides_MappedTensor2D_gpu(real *extBoundary, real *boundary, int N, int nVar){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t ivar = threadIdx.y;
  
  boundary[TEB_2D_INDEX(1,1,i1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_2D_INDEX(1,1,i1,ivar,s1,e1,N,nVar)]+
                                                        boundary[TEB_2D_INDEX(1,1,i1,ivar,s1,e1,N,nVar)]);
  boundary[TEB_2D_INDEX(2,1,i1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_2D_INDEX(2,1,i1,ivar,s1,e1,N,nVar)]+
                                                        boundary[TEB_2D_INDEX(2,1,i1,ivar,s1,e1,N,nVar)]);
  boundary[TEB_2D_INDEX(1,2,i1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_2D_INDEX(1,2,i1,ivar,s1,e1,N,nVar)]+
                                                        boundary[TEB_2D_INDEX(1,2,i1,ivar,s1,e1,N,nVar)]);
  boundary[TEB_2D_INDEX(2,2,i1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_2D_INDEX(2,2,i1,ivar,s1,e1,N,nVar)]+
                                                        boundary[TEB_2D_INDEX(2,2,i1,ivar,s1,e1,N,nVar)]);
  
}

extern "C"
{
  void BassiRebaySides_MappedTensor2D_gpu_wrapper(real **extBoundary, real **boundary, int N, int nVar, int nEl)
  {
    BassiRebaySides_MappedTensor2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*extBoundary, *boundary, N, nVar);
  }

}

__global__ void BassiRebaySides_MappedScalar3D_gpu(real *avgBoundary, real *boundary, real *extBoundary, int N, int nVar){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  size_t ivar = threadIdx.z;
  
  avgBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)]+
                                                         boundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)]);
  
}

extern "C"
{
  void BassiRebaySides_MappedScalar3D_gpu_wrapper(real **avgBoundary, real **boundary, real **extBoundary, int N, int nVar, int nEl)
  {
    BassiRebaySides_MappedScalar3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nVar), 0, 0>>>(*avgBoundary, *boundary, *extBoundary, N, nVar);
  }

}

__global__ void BassiRebaySides_MappedVector3D_gpu(real *extBoundary, real *boundary, int N, int nVar){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  size_t ivar = threadIdx.z;
  
  boundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)]+
                                                      boundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)]);
  boundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nVar)]+
                                                      boundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)]);
  boundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nVar)]+
                                                      boundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)]);
  
}

extern "C"
{
  void BassiRebaySides_MappedVector3D_gpu_wrapper(real **extBoundary, real **boundary, int N, int nVar, int nEl)
  {
    BassiRebaySides_MappedVector3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nVar), 0, 0>>>(*extBoundary, *boundary, N, nVar);
  }

}

__global__ void BassiRebaySides_MappedTensor3D_gpu(real *extBoundary, real *boundary, int N, int nVar){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  size_t ivar = threadIdx.z;
  
      boundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nVar)]+
                                                            boundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nVar)]);
      boundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nVar)]+
                                                            boundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nVar)]);
      boundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nVar)]+
                                                            boundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nVar)]);
      boundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nVar)]+
                                                            boundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nVar)]);
      boundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nVar)]+
                                                            boundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nVar)]);
      boundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nVar)]+
                                                            boundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nVar)]);
      boundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nVar)]+
                                                            boundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nVar)]);
      boundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nVar)]+
                                                            boundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nVar)]);
      boundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nVar)] = 0.5*(extBoundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nVar)]+
                                                            boundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nVar)]);
  
}

extern "C"
{
  void BassiRebaySides_MappedTensor3D_gpu_wrapper(real **extBoundary, real **boundary, int N, int nVar, int nEl)
  {
    BassiRebaySides_MappedTensor3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nVar), 0, 0>>>(*extBoundary, *boundary, N, nVar);
  }

}

__global__ void ApplyFlip_MappedScalar2D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y;
  size_t e1 = blockIdx.z;
  size_t i1 = threadIdx.x;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1,e1,5,4)];
  int i2 = N-i1;

  __shared__ real extBuff[16];

  extBuff[i1] = extBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank /= rankId ){
      if(flip == 1){
        extBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nVar)] = extBuff[i2];
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedScalar2D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedScalar2D_gpu<<<dim3(nVar,4,nEl), dim3(N+1,1,1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nVar);
  }

}

__global__ void ApplyFlip_MappedVector2D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y;
  size_t e1 = blockIdx.z;
  size_t dir = threadIdx.x;
  size_t i1 = threadIdx.y;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1,e1,5,4)];
  int i2 = N-i1;

  __shared__ real extBuff[32];

  extBuff[dir+2*i1] = extBoundary[VEB_2D_INDEX(dir+1,i1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank /= rankId ){
      if(flip == 1){
        extBoundary[VEB_2D_INDEX(dir+1,i1,ivar,s1,e1,N,nVar)] = extBuff[dir+2*i2];
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedVector2D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedVector2D_gpu<<<dim3(nVar,4,nEl), dim3(2,N+1,1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nVar);
  }

}

__global__ void ApplyFlip_MappedTensor2D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y;
  size_t e1 = blockIdx.z;
  size_t row = threadIdx.x;
  size_t col = threadIdx.y;
  size_t i1 = threadIdx.z;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1,e1,5,4)];
  int i2 = N-i1;

  __shared__ real extBuff[64];

  extBuff[row+2*(col+2*i1)] = extBoundary[TEB_2D_INDEX(row+1,col+1,i1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank /= rankId ){
      if(flip == 1){
        extBoundary[TEB_2D_INDEX(row+1,col+1,i1,ivar,s1,e1,N,nVar)] = extBuff[row+2*(col+2*i2)];
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedTensor2D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedTensor2D_gpu<<<dim3(nVar,4,nEl), dim3(2,2,N+1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nVar);
  }

}

__global__ void ApplyFlip_MappedScalar3D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y;
  size_t e1 = blockIdx.z;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.x;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1,e1,5,4)];

  __shared__ real extBuff[256];

  extBuff[i1+(N+1)*j1] = extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank /= rankId ){
      if(flip == 2){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 3){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = extBuff[i2+(N+1)*j2];
      }
    }
  }
}

extern "C"
{
  void ApplyFlip_MappedScalar3D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedScalar3D_gpu<<<dim3(nVar,6,nEl), dim3(N+1,N+1,1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nVar);
  }

}

__global__ void ApplyFlip_MappedVector3D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y;
  size_t e1 = blockIdx.z;
  size_t dir = threadIdx.x;
  size_t i1 = threadIdx.y;
  size_t j1 = threadIdx.z;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1,e1,5,4)];

  __shared__ real extBuff[768];

  extBuff[dir+3*(i1+(N+1)*j1)] = extBoundary[VEB_3D_INDEX(dir+1,i1,j1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank /= rankId ){
      if(flip == 2){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[dir+3*(i2+(N+1)*j2)];
      }
      else if(flip == 3){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[dir+3*(i2+(N+1)*j2)];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[dir+3*(i2+(N+1)*j2)];
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedVector3D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedVector3D_gpu<<<dim3(nVar,6,nEl), dim3(3,N+1,N+1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nVar);
  }

}

__global__ void ApplyFlip_MappedTensor3D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y;
  size_t e1 = blockIdx.z;
  size_t dir = threadIdx.x;
  size_t i1 = threadIdx.y;
  size_t j1 = threadIdx.z;
  size_t row = dir/3;
  size_t col = dir - dir*row;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1,e1,5,4)];

  __shared__ real extBuff[2304];

  extBuff[row+3*(col+3*(i1+(N+1)*j1))] = extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank /= rankId ){
      if(flip == 2){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[row+3*(col+3*(i2+(N+1)*j2))];
      }
      else if(flip == 3){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[row+3*(col+3*(i2+(N+1)*j2))];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[row+3*(col+3*(i2+(N+1)*j2))];
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedTensor3D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedTensor3D_gpu<<<dim3(nVar,6,nEl), dim3(9,N+1,N+1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nVar);
  }

}
