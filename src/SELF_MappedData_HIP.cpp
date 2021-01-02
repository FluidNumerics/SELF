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
