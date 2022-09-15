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

__global__ void Determinant_Tensor2D_gpu(real *tensor, real *detTensor, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    detTensor[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] = tensor[TE_2D_INDEX(1,1,i,j,iVar,iEl,N,nVar)]*
	                                          tensor[TE_2D_INDEX(2,2,i,j,iVar,iEl,N,nVar)]-
						  tensor[TE_2D_INDEX(1,2,i,j,iVar,iEl,N,nVar)]*
						  tensor[TE_2D_INDEX(2,1,i,j,iVar,iEl,N,nVar)];
}
extern "C"
{
  void Determinant_Tensor2D_gpu_wrapper(real **tensor, real **detTensor, int N, int nVar, int nEl)
  {
	  Determinant_Tensor2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*tensor, *detTensor, N, nVar);
  } 
}

__global__ void Determinant_Tensor3D_gpu(real *tensor, real *detTensor, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    detTensor[SC_3D_INDEX(i,j,k,iVar,iEl,N,nVar)] = tensor[TE_3D_INDEX(1,1,i,j,k,iVar,iEl,N,nVar)]*
	                                            (tensor[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,N,nVar)]*
						     tensor[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,N,nVar)]-
						     tensor[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)]*
						     tensor[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)])-
                                                    tensor[TE_3D_INDEX(2,1,i,j,k,iVar,iEl,N,nVar)]*
	                                            (tensor[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)]*
						     tensor[TE_3D_INDEX(3,3,i,j,k,iVar,iEl,N,nVar)]-
						     tensor[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)]*
						     tensor[TE_3D_INDEX(3,2,i,j,k,iVar,iEl,N,nVar)])+
                                                    tensor[TE_3D_INDEX(3,1,i,j,k,iVar,iEl,N,nVar)]*
	                                            (tensor[TE_3D_INDEX(1,2,i,j,k,iVar,iEl,N,nVar)]*
						     tensor[TE_3D_INDEX(2,3,i,j,k,iVar,iEl,N,nVar)]-
						     tensor[TE_3D_INDEX(1,3,i,j,k,iVar,iEl,N,nVar)]*
						     tensor[TE_3D_INDEX(2,2,i,j,k,iVar,iEl,N,nVar)]);
}
extern "C"
{
  void Determinant_Tensor3D_gpu_wrapper(real **tensor, real **detTensor, int N, int nVar, int nEl)
  {
	  Determinant_Tensor3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*tensor, *detTensor, N, nVar);
  } 
}
