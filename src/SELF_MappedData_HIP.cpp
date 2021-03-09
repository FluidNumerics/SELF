#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

/*
// Template
__global__ void Template_{1D|2D|3D}_gpu( , int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t k = hipThreadIdx_z;


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

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;

    scalar[SC_1D_INDEX(i,iVar,iEl,N,nVar)] = scalar[SC_1D_INDEX(i,iVar,iEl,N,nVar)]/
                                             dxds[SC_1D_INDEX(i,iVar,iEl,N,nVar)];
}

extern "C"
{
  void JacobianWeight_MappedScalar1D_gpu_wrapper(real **scalar, real **dxds, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((JacobianWeight_MappedScalar1D_gpu), dim3(nVar,nEl,1), dim3(N+1,1,1), 0, 0, *scalar, *dxds, N, nVar);
  }
}

// JacobianWeight_MappedScalar2D_gpu
__global__ void JacobianWeight_MappedScalar2D_gpu(real *scalar, real *jacobian, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;

    scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]/
                                               jacobian[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];
}

extern "C"
{
  void JacobianWeight_MappedScalar2D_gpu_wrapper(real **scalar, real **jacobian, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((JacobianWeight_MappedScalar2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *scalar, *jacobian, N, nVar);
  }
}

// JacobianWeight_MappedScalar3D_gpu
__global__ void JacobianWeight_MappedScalar3D_gpu(real *scalar, real *jacobian, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t k = hipThreadIdx_z;

    scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)] = scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];
}

extern "C"
{
  void JacobianWeight_MappedScalar3D_gpu_wrapper(real **scalar, real **jacobian, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((JacobianWeight_MappedScalar3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *scalar, *jacobian, N, nVar);
  }
}

// ContravariantWeight_MappedScalar2D_gpu
__global__ void ContravariantWeight_MappedScalar2D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;

    workTensor[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(1,1,i,j,1,iEl,N,1)]*
                                                       scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]; 
  
    workTensor[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(1,2,i,j,1,iEl,N,1)]*
                                                       scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]; 
  
    workTensor[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(2,1,i,j,1,iEl,N,1)]*
                                                       scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]; 
  
    workTensor[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(2,2,i,j,1,iEl,N,1)]*
                                                       scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]; 
}

extern "C"
{
  void ContravariantWeight_MappedScalar2D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((ContravariantWeight_MappedScalar2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *scalar, *workTensor, *dsdx, N, nVar);
  } 
}

// ContravariantWeightBoundary_MappedScalar2D_gpu
__global__ void ContravariantWeightBounday_MappedScalar2D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nVar){

  size_t iSide = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t iVar = hipThreadIdx_y;

#define TEB_3D_INDEX(row,col,i,iVar,iSide,iel,N,nVar) row-1 + 3*(col-1 + 3*(i + (N+1)*(j + (N+1)*(iVar + nVar*(iSide-1 + 6*iel)))))
    workTensor[TEB_2D_INDEX(1,1,i,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_2D_INDEX(1,1,i,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)]; 
  
    workTensor[TEB_2D_INDEX(2,1,i,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_2D_INDEX(1,2,i,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)]; 
  
    workTensor[TEB_2D_INDEX(1,2,i,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_2D_INDEX(2,1,i,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)]; 
  
    workTensor[TEB_2D_INDEX(2,2,i,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_2D_INDEX(2,2,i,j,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)]; 
}

extern "C"
{
  void ContravariantWeightBoundary_MappedScalar2D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((ContravariantWeightBoundary_MappedScalar2D_gpu), dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0, *scalar, *workTensor, *dsdx, N, nVar);
  } 
}

// ContravariantWeight_MappedScalar3D_gpu
__global__ void ContravariantWeight_MappedScalar3D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t k = hipThreadIdx_z;

  workTensor[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,1,i,j,k,1,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,2,i,j,k,1,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,3,i,j,k,1,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(2,1,i,j,k,1,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(2,2,i,j,k,1,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(2,3,i,j,k,1,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(3,1,i,j,k,1,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(3,2,i,j,k,1,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

  workTensor[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(3,3,i,j,k,1,iEl,N,1)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)]; 

}

extern "C"
{
  void ContravariantWeight_MappedScalar3D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((ContravariantWeight_MappedScalar3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *scalar, *workTensor, *dsdx, N, nVar);
  } 
}

// ContravariantWeightBoundary_MappedScalar3D_gpu
__global__ void ContravariantWeightBoundary_MappedScalar3D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nVar){

  size_t iSide = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t iVar = hipThreadIdx_z;

  workTensor[TEB_3D_INDEX(1,1,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(1,1,i,j,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(2,1,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(1,2,i,j,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(3,1,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(1,3,i,j,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(1,2,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(2,1,i,j,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(2,2,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(2,2,i,j,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(3,2,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(2,3,i,j,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(1,3,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(3,1,i,j,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(2,3,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(3,2,i,j,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

  workTensor[TEB_3D_INDEX(3,3,i,j,iVar,iSide,iEl,N,nVar)] = dsdx[TEB_3D_INDEX(3,3,i,j,1,iSide,iEl,N,1)]*
                                                       scalar[SCB_3D_INDEX(i,j,iVar,iSide,iEl,N,nVar)]; 

}

extern "C"
{
  void ContravariantWeightBoundary_MappedScalar3D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((ContravariantWeightBoundary_MappedScalar3D_gpu), dim3(6,nEl,1), dim3(N+1,N+1,nVar), 0, 0, *scalar, *workTensor, *dsdx, N, nVar);
  } 
}

// ContravariantProjection_MappedVector2D_gpu
__global__ void ContravariantProjection_MappedVector2D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;

    workTensor[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(1,1,i,j,1,iEl,N,1)]*
                                                       scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]; 
  
    workTensor[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(1,2,i,j,1,iEl,N,1)]*
                                                       scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]; 
  
    workTensor[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(2,1,i,j,1,iEl,N,1)]*
                                                       scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]; 
  
    workTensor[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(2,2,i,j,1,iEl,N,1)]*
                                                       scalar[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]; 
}

extern "C"
{
  void ContravariantProjection_MappedVector2D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((ContravariantProjection_MappedVector2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *scalar, *workTensor, *dsdx, N, nVar);
  } 
}

// JacobianWeight_MappedVector2D_gpu
__global__ void JacobianWeight_MappedVector2D_gpu(real *vector, real *jacobian, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;

    vector[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] = vector[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

    vector[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] = vector[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];
}

extern "C"
{
  void JacobianWeight_MappedVector2D_gpu_wrapper(real **vector, real **jacobian, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((JacobianWeight_MappedVector2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *vector, *jacobian, N, nVar);
  }
}

// ContravariantProjection_MappedVector3D_gpu
__global__ void ContravariantProjection_MappedVector3D_gpu(real *physVector, real *compVector, real *dsdx, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t k = hipThreadIdx_z;

  compVector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,1,i,j,k,1,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)]+ 
                                                     dsdx[TE_3D_INDEX(2,1,i,j,k,1,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)]+
                                                     dsdx[TE_3D_INDEX(3,1,i,j,k,1,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)];

  compVector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,2,i,j,k,1,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)]+ 
                                                     dsdx[TE_3D_INDEX(2,2,i,j,k,1,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)]+
                                                     dsdx[TE_3D_INDEX(3,2,i,j,k,1,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)];

  compVector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,3,i,j,k,1,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)]+ 
                                                     dsdx[TE_3D_INDEX(2,3,i,j,k,1,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)]+
                                                     dsdx[TE_3D_INDEX(3,3,i,j,k,1,iEl,N,1)]*
                                                     physVector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)];

}

extern "C"
{
  void ContravariantProjection_MappedVector3D_gpu_wrapper(real **physVector, real **compVector, real **dsdx, int N, int nVar, int nEl)
  {
	  hipLaunchKernelGGL((ContravariantProjection_MappedVector3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *physVector, *compVector, *dsdx, N, nVar);
  } 
}

// JacobianWeight_MappedVector3D_gpu
__global__ void JacobianWeight_MappedVector3D_gpu(real *vector, real *jacobian, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t k = hipThreadIdx_z;

    vector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

    vector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

    vector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] = vector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];
}

extern "C"
{
  void JacobianWeight_MappedVector3D_gpu_wrapper(real **vector, real **jacobian, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((JacobianWeight_MappedVector3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *vector, *jacobian, N, nVar);
  }
}

// JacobianWeight_MappedTensor2D_gpu
__global__ void JacobianWeight_MappedTensor2D_gpu(real *tensor, real *jacobian, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;

    tensor[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)] = tensor[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

    tensor[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)] = tensor[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

    tensor[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)] = tensor[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

    tensor[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)] = tensor[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)]/
                                                 jacobian[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];
}

extern "C"
{
  void JacobianWeight_MappedTensor2D_gpu_wrapper(real **tensor, real **jacobian, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((JacobianWeight_MappedTensor2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *tensor, *jacobian, N, nVar);
  }
}

// JacobianWeight_MappedTensor3D_gpu
__global__ void JacobianWeight_MappedTensor3D_gpu(real *tensor, real *jacobian, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t k = hipThreadIdx_z;

    tensor[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

    tensor[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

    tensor[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

    tensor[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

    tensor[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

    tensor[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

    tensor[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

    tensor[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];

    tensor[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,N,nVar)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)];
}

extern "C"
{
  void JacobianWeight_MappedTensor3D_gpu_wrapper(real **tensor, real **jacobian, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((JacobianWeight_MappedTensor3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *tensor, *jacobian, N, nVar);
  }
}

// CalculateCurl_MappedTensor2D_gpu
__global__ void CalculateCurl_MappedTensor2D_gpu(real *dfdx, real *curlf, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;

    curlf[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = dfdx[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)]-
                                              dfdx[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)];
}

extern "C"
{
  void CalculateCurl_MappedTensor2D_gpu_wrapper(real **dfdx, real **curlf, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((JacobianWeight_MappedTensor2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *dfdx, *curlf, N, nVar);
  }
}

// CalculateCurl_MappedTensor3D_gpu
__global__ void CalculateCurl_MappedTensor3D_gpu(real *dfdx, real *curlf, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t k = hipThreadIdx_z;

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
    hipLaunchKernelGGL((JacobianWeight_MappedTensor3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *dfdx, *curlf, N, nVar);
  }
}

// MapTo support routines
__global__ void MapToScalar_MappedVector2D_gpu(real* scalar, real* vector, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;

    jVar = 1 + 2*(iVar-1);
    scalar[SC_2D_INDEX(i,j,jVar,iEl,N,nVar*2)] = vector[VE_2D_INDEX(1,i,j,ivar,iel,N,nVar)];

    jVar = 2 + 2*(iVar-1);
    scalar[SC_2D_INDEX(i,j,jVar,iEl,N,nVar*2)] = vector[VE_2D_INDEX(2,i,j,ivar,iel,N,nVar)];

}

extern "C"
{
  void MapToScalar_MappedVector2D_gpu_wrapper(real **scalar, real **vector, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((MapToScalar_MappedVector2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *scalar, *vector, N, nVar);
  }
}

__global__ void MapToScalarBoundary_MappedVector2D_gpu(real* scalar, real* vector, int N, int nVar){

  size_t iSide = hipBlockIdx_x;
  size_t iVar = hipBlockIdx_y;
  size_t iEl = hipBlockIdx_z;
  size_t j = hipThreadIdx_x;

    jVar = 1 + 2*(iVar-1);
    scalar[SCB_2D_INDEX(j,jVar,iSide,iEl,N,nVar*2)] = vector[VEB_2D_INDEX(1,j,ivar,iSide,iel,N,nVar)];

    jVar = 2 + 2*(iVar-1);
    scalar[SCB_2D_INDEX(j,jVar,iSide,iEl,N,nVar*2)] = vector[VEB_2D_INDEX(2,j,ivar,iSide,iel,N,nVar)];

}

extern "C"
{
  void MapToScalarBoundary_MappedVector2D_gpu_wrapper(real **scalar, real **vector, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((MapToScalarBoundary_MappedVector2D_gpu), dim3(4,nVar,nEl), dim3(N+1,1,1), 0, 0, *scalar, *vector, N, nVar);
  }
}

__global__ void MapToTensor_MappedVector2D_gpu(real* tensor, real* vector, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;

    jVar = 1 + 2*(iVar-1);
    tensor[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar*2)] = vector[VE_2D_INDEX(1,i,j,jvar,iel,N,nVar)];

    tensor[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar*2)] = vector[VE_2D_INDEX(2,i,j,jvar,iel,N,nVar)];

    jVar = 2 + 2*(iVar-1);
    tensor[TE_2D_INDEX(2,1,i,j,jVar,iEl,N,nVar*2)] = vector[VE_2D_INDEX(1,i,j,jvar,iel,N,nVar)];

    tensor[TE_2D_INDEX(2,2,i,j,jVar,iEl,N,nVar*2)] = vector[VE_2D_INDEX(2,i,j,jvar,iel,N,nVar)];

}

extern "C"
{
  void MapToTensor_MappedVector2D_gpu_wrapper(real **tensor, real **vector, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((MapToTensor_MappedVector2D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0, *tensor, *vector, N, nVar);
  }
}

__global__ void MapToTensorBoundary_MappedVector2D_gpu(real* tensor, real* vector, int N, int nVar){

  size_t iSide = hipBlockIdx_x;
  size_t iVar = hipBlockIdx_y;
  size_t iEl = hipBlockIdx_z;
  size_t j = hipThreadIdx_x;

    jVar = 1 + 2*(iVar-1);
    tensor[TEB_2D_INDEX(1,1,j,iVar,iSide,iEl,N,nVar*2)] = vector[VEB_2D_INDEX(1,j,jvar,iSide,iel,N,nVar)];

    tensor[TEB_2D_INDEX(1,2,j,iVar,iSide,iEl,N,nVar*2)] = vector[VEB_2D_INDEX(2,j,jvar,iSide,iel,N,nVar)];

    jVar = 2 + 2*(iVar-1);
    tensor[TEB_2D_INDEX(2,1,j,iVar,iSide,iEl,N,nVar*2)] = vector[VEB_2D_INDEX(1,j,jvar,iSide,iel,N,nVar)];

    tensor[TEB_2D_INDEX(2,2,j,iVar,iSide,iEl,N,nVar*2)] = vector[VEB_2D_INDEX(2,j,jvar,iSide,iel,N,nVar)];

}

extern "C"
{
  void MapToTensorBoundary_MappedVector2D_gpu_wrapper(real **tensor, real **vector, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((MapToTensorBoundary_MappedVector2D_gpu), dim3(4,nVar,nEl), dim3(N+1,1,1), 0, 0, *tensor, *vector, N, nVar);
  }
}

__global__ void MapToScalar_MappedVector3D_gpu(real* scalar, real* vector, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t k = hipThreadIdx_z;

    jVar = 1 + 3*(iVar-1);
    scalar[SC_3D_INDEX(i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(1,i,j,k,ivar,iel,N,nVar)];

    jVar = 2 + 3*(iVar-1);
    scalar[SC_3D_INDEX(i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(2,i,j,k,ivar,iel,N,nVar)];

    jVar = 3 + 3*(iVar-1);
    scalar[SC_3D_INDEX(i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(3,i,j,k,ivar,iel,N,nVar)];

}

extern "C"
{
  void MapToScalar_MappedVector3D_gpu_wrapper(real **scalar, real **vector, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((MapToScalar_MappedVector3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *scalar, *vector, N, nVar);
  }
}

__global__ void MapToScalarBoundary_MappedVector3D_gpu(real* scalar, real* vector, int N, int nVar){

  size_t iSide = hipBlockIdx_x;
  size_t iVar = hipBlockIdx_y;
  size_t iEl = hipBlockIdx_z;
  size_t j = hipThreadIdx_x;
  size_t k = hipThreadIdx_y;

    jVar = 1 + 3*(iVar-1);
    scalar[SCB_3D_INDEX(j,k,jVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(1,j,k,ivar,iSide,iel,N,nVar)];

    jVar = 2 + 3*(iVar-1);
    scalar[SCB_3D_INDEX(j,k,jVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(2,j,k,ivar,iSide,iel,N,nVar)];

    jVar = 3 + 3*(iVar-1);
    scalar[SCB_3D_INDEX(j,k,jVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(3,j,k,ivar,iSide,iel,N,nVar)];

}

extern "C"
{
  void MapToScalarBoundary_MappedVector3D_gpu_wrapper(real **scalar, real **vector, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((MapToScalarBoundary_MappedVector3D_gpu), dim3(6,nVar,nEl), dim3(N+1,N+1,1), 0, 0, *scalar, *vector, N, nVar);
  }
}

__global__ void MapToTensor_MappedVector3D_gpu(real* tensor, real* vector, int N, int nVar){

  size_t iVar = hipBlockIdx_x;
  size_t iEl = hipBlockIdx_y;
  size_t i = hipThreadIdx_x;
  size_t j = hipThreadIdx_y;
  size_t k = hipThreadIdx_z;

    jVar = 1 + 3*(iVar-1);
    tensor[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(1,i,j,k,jvar,iel,N,nVar)];

    tensor[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(2,i,j,k,jvar,iel,N,nVar)];

    tensor[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(3,i,j,k,jvar,iel,N,nVar)];

    jVar = 2 + 3*(iVar-1);
    tensor[TE_3D_INDEX(2,1,i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(1,i,j,k,jvar,iel,N,nVar)];

    tensor[TE_3D_INDEX(2,2,i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(2,i,j,k,jvar,iel,N,nVar)];

    tensor[TE_3D_INDEX(2,3,i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(3,i,j,k,jvar,iel,N,nVar)];

    jVar = 3 + 3*(iVar-1);
    tensor[TE_3D_INDEX(3,1,i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(1,i,j,k,jvar,iel,N,nVar)];

    tensor[TE_3D_INDEX(3,2,i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(2,i,j,k,jvar,iel,N,nVar)];

    tensor[TE_3D_INDEX(3,3,i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(3,i,j,k,jvar,iel,N,nVar)];

}

extern "C"
{
  void MapToTensor_MappedVector3D_gpu_wrapper(real **tensor, real **vector, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((MapToTensor_MappedVector3D_gpu), dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0, *tensor, *vector, N, nVar);
  }
}

__global__ void MapToTensorBoundary_MappedVector3D_gpu(real* tensor, real* vector, int N, int nVar){

  size_t iSide = hipBlockIdx_x;
  size_t iVar = hipBlockIdx_y;
  size_t iEl = hipBlockIdx_z;
  size_t j = hipThreadIdx_x;
  size_t k = hipThreadIdx_x;

    jVar = 1 + 3*(iVar-1);
    tensor[TEB_3D_INDEX(1,1,j,k,iVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(1,j,k,jvar,iSide,iel,N,nVar)];

    tensor[TEB_3D_INDEX(1,2,j,k,iVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(2,j,k,jvar,iSide,iel,N,nVar)];

    tensor[TEB_3D_INDEX(1,3,j,k,iVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(3,j,k,jvar,iSide,iel,N,nVar)];

    jVar = 2 + 3*(iVar-1);
    tensor[TEB_3D_INDEX(2,1,j,k,iVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(1,j,k,jvar,iSide,iel,N,nVar)];

    tensor[TEB_3D_INDEX(2,2,j,k,iVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(2,j,k,jvar,iSide,iel,N,nVar)];

    tensor[TEB_3D_INDEX(2,3,j,k,iVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(3,j,k,jvar,iSide,iel,N,nVar)];

    jVar = 3 + 3*(iVar-1);
    tensor[TEB_3D_INDEX(3,1,j,k,iVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(1,j,k,jvar,iSide,iel,N,nVar)];

    tensor[TEB_3D_INDEX(3,2,j,k,iVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(2,j,k,jvar,iSide,iel,N,nVar)];

    tensor[TEB_3D_INDEX(3,3,j,k,iVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(3,j,k,jvar,iSide,iel,N,nVar)];

}

extern "C"
{
  void MapToTensorBoundary_MappedVector3D_gpu_wrapper(real **tensor, real **vector, int N, int nVar, int nEl)
  {
    hipLaunchKernelGGL((MapToTensorBoundary_MappedVector3D_gpu), dim3(6,nVar,nEl), dim3(N+1,N+1,1), 0, 0, *tensor, *vector, N, nVar);
  }
}
