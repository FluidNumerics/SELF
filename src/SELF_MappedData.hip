#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

// JacobianWeight_MappedScalar1D_gpu
__global__ void JacobianWeight_MappedScalar1D_gpu(real *scalar, real *dxds, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;

    scalar[SC_1D_INDEX(i,iEl,iVar,N,nEl)] = scalar[SC_1D_INDEX(i,iEl,iVar,N,nEl)]/
                                             dxds[SC_1D_INDEX(i,iEl,0,N,nEl)];
}

extern "C"
{
  void JacobianWeight_MappedScalar1D_gpu_wrapper(real **scalar, real **dxds, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedScalar1D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,1,1), 0, 0>>>(*scalar, *dxds, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// GradientBR_MappedScalar2D_gpu
__global__ void GradientBR_MappedScalar2D_gpu(real *scalar, real *avgBoundary, real *dsdx, real *jacobian, real* nHat, real *nScale, real *gradF, real *dgMatrix, real *bMatrix, real *qWeights, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real gfx = 0.0;
  real gfy = 0.0;
  real f1, f2;

  for( int ii=0; ii < N+1; ii++ ){
	  
    f1 = scalar[SC_2D_INDEX(ii,j,iEl,iVar,N,nEl)]*
	 dsdx[TE_2D_INDEX(1,1,ii,j,iEl,0,N,nEl)];
    
    f2 = scalar[SC_2D_INDEX(i,ii,iEl,iVar,N,nEl)]*
	 dsdx[TE_2D_INDEX(1,2,i,ii,iEl,0,N,nEl)];

    gfx += dgMatrix[ii+i*(N+1)]*f1+
           dgMatrix[ii+j*(N+1)]*f2;

    f1 = scalar[SC_2D_INDEX(ii,j,iEl,iVar,N,nEl)]*
	 dsdx[TE_2D_INDEX(2,1,ii,j,iEl,0,N,nEl)];
    
    f2 = scalar[SC_2D_INDEX(i,ii,iEl,iVar,N,nEl)]*
	 dsdx[TE_2D_INDEX(2,2,i,ii,iEl,0,N,nEl)];

    gfy += dgMatrix[ii+i*(N+1)]*f1+
           dgMatrix[ii+j*(N+1)]*f2;

  }

  // East
  f1 = avgBoundary[SCB_2D_INDEX(j,iVar,2,iEl,N,nEl)]*
	  nHat[VEB_2D_INDEX(1,j,0,2,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(j,0,2,iEl,N,1)];
  // West
  f2 = avgBoundary[SCB_2D_INDEX(j,iVar,4,iEl,N,nEl)]*
	  nHat[VEB_2D_INDEX(1,j,0,4,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(j,0,4,iEl,N,1)];

  gfx += (f1*bMatrix[i+N+1] + f2*bMatrix[i])/qWeights[i];

  // North
  f1 = avgBoundary[SCB_2D_INDEX(i,iVar,3,iEl,N,nEl)]*
	  nHat[VEB_2D_INDEX(1,i,0,3,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(i,0,3,iEl,N,1)];
  // South
  f2 = avgBoundary[SCB_2D_INDEX(i,iVar,1,iEl,N,nEl)]*
	  nHat[VEB_2D_INDEX(1,i,0,1,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(i,0,1,iEl,N,1)];

  gfx += (f1*bMatrix[j+N+1] + f2*bMatrix[j])/qWeights[j];

  // East
  f1 = avgBoundary[SCB_2D_INDEX(j,iVar,2,iEl,N,nEl)]*
	  nHat[VEB_2D_INDEX(2,j,0,2,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(j,0,2,iEl,N,1)];
  // West
  f2 = avgBoundary[SCB_2D_INDEX(j,iVar,4,iEl,N,nEl)]*
	  nHat[VEB_2D_INDEX(2,j,0,4,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(j,0,4,iEl,N,1)];

  gfy += (f1*bMatrix[i+N+1] + f2*bMatrix[i])/qWeights[i];

  // North
  f1 = avgBoundary[SCB_2D_INDEX(i,iVar,3,iEl,N,nEl)]*
	  nHat[VEB_2D_INDEX(2,i,0,3,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(i,0,3,iEl,N,1)];
  // South
  f2 = avgBoundary[SCB_2D_INDEX(i,iVar,1,iEl,N,nEl)]*
	  nHat[VEB_2D_INDEX(2,i,0,1,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(i,0,1,iEl,N,1)];

  gfy += (f1*bMatrix[j+N+1] + f2*bMatrix[j])/qWeights[j];

  gradF[VE_2D_INDEX(1,i,j,iEl,iVar,N,nEl)] = gfx/jacobian[SC_2D_INDEX(i,j,iEl,0,N,nEl)];
  gradF[VE_2D_INDEX(2,i,j,iEl,iVar,N,nEl)] = gfy/jacobian[SC_2D_INDEX(i,j,iEl,0,N,nEl)];

}

extern "C"
{
  void GradientBR_MappedScalar2D_gpu_wrapper(real **scalar, real **avgBoundary, real **dsdx, real **jacobian, real **nHat, real **nScale, real **gradF, real **dgMatrix, real **bMatrix, real **qWeights, int N, int nVar, int nEl)
  {
    GradientBR_MappedScalar2D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(*scalar, *avgBoundary, *dsdx, *jacobian, *nHat, *nScale, *gradF, *dgMatrix, *bMatrix, *qWeights, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// GradientSF_MappedScalar2D_gpu
__global__ void GradientSF_MappedScalar2D_gpu(real *scalar, real *dsdx, real *jacobian, real *gradF, real *dMatrix, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real gf[2] = {0.0};

  for( int ii=0; ii < N+1; ii++ ){
    gf[0] += dMatrix[ii+i*(N+1)]*
	       scalar[SC_2D_INDEX(ii,j,iEl,iVar,N,nEl)]*
	       dsdx[TE_2D_INDEX(1,1,ii,j,iEl,0,N,nEl)]+
             dMatrix[ii+j*(N+1)]*
               scalar[SC_2D_INDEX(i,ii,iEl,iVar,N,nEl)]*   	       
               dsdx[TE_2D_INDEX(1,2,i,ii,iEl,0,N,nEl)];

    gf[1] += dMatrix[ii+i*(N+1)]*
	       scalar[SC_2D_INDEX(ii,j,iEl,iVar,N,nEl)]*
	       dsdx[TE_2D_INDEX(2,1,ii,j,iEl,0,N,nEl)]+
             dMatrix[ii+j*(N+1)]*
               scalar[SC_2D_INDEX(i,ii,iEl,iVar,N,nEl)]*   	       
               dsdx[TE_2D_INDEX(2,2,i,ii,iEl,0,N,nEl)];

  }
  gradF[VE_2D_INDEX(1,i,j,iEl,iVar,N,nEl)] = gf[0]/jacobian[SC_2D_INDEX(i,j,iEl,0,N,nEl)];
  gradF[VE_2D_INDEX(2,i,j,iEl,iVar,N,nEl)] = gf[1]/jacobian[SC_2D_INDEX(i,j,iEl,0,N,nEl)];

}

extern "C"
{
  void GradientSF_MappedScalar2D_gpu_wrapper(real **scalar, real **dsdx, real **jacobian, real **gradF, real **dMatrix, int N, int nVar, int nEl)
  {
    GradientSF_MappedScalar2D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(*scalar, *dsdx, *jacobian, *gradF, *dMatrix, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// JacobianWeight_MappedScalar2D_gpu
__global__ void JacobianWeight_MappedScalar2D_gpu(real *scalar, real *jacobian, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    scalar[SC_2D_INDEX(i,j,iEl,iVar,N,nEl)] = scalar[SC_2D_INDEX(i,j,iEl,iVar,N,nEl)]/
                                               jacobian[SC_2D_INDEX(i,j,iEl,0,N,nEl)];
}

extern "C"
{
  void JacobianWeight_MappedScalar2D_gpu_wrapper(real **scalar, real **jacobian, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedScalar2D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(*scalar, *jacobian, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// JacobianWeight_MappedScalar3D_gpu
__global__ void JacobianWeight_MappedScalar3D_gpu(real *scalar, real *jacobian, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    scalar[SC_3D_INDEX(i,j,k,iEl,iVar,N,nEl)] = scalar[SC_3D_INDEX(i,j,k,iEl,iVar,N,nEl)]/
                                                 jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];
}

extern "C"
{
  void JacobianWeight_MappedScalar3D_gpu_wrapper(real **scalar, real **jacobian, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedScalar3D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,N+1), 0, 0>>>(*scalar, *jacobian, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// ContravariantWeight_MappedScalar2D_gpu
__global__ void ContravariantWeight_MappedScalar2D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    workTensor[TE_2D_INDEX(1,1,i,j,iEl,iVar,N,nEl)] = dsdx[TE_2D_INDEX(1,1,i,j,iEl,0,N,nEl)]*
                                                       scalar[SC_2D_INDEX(i,j,iEl,iVar,N,nEl)]; 
  
    workTensor[TE_2D_INDEX(2,1,i,j,iEl,iVar,N,nEl)] = dsdx[TE_2D_INDEX(1,2,i,j,iEl,0,N,nEl)]*
                                                       scalar[SC_2D_INDEX(i,j,iEl,iVar,N,nEl)]; 
  
    workTensor[TE_2D_INDEX(1,2,i,j,iEl,iVar,N,nEl)] = dsdx[TE_2D_INDEX(2,1,i,j,iEl,0,N,nEl)]*
                                                       scalar[SC_2D_INDEX(i,j,iEl,iVar,N,nEl)]; 
  
    workTensor[TE_2D_INDEX(2,2,i,j,iEl,iVar,N,nEl)] = dsdx[TE_2D_INDEX(2,2,i,j,iEl,0,N,nEl)]*
                                                       scalar[SC_2D_INDEX(i,j,iEl,iVar,N,nEl)]; 
}

extern "C"
{
  void ContravariantWeight_MappedScalar2D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantWeight_MappedScalar2D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(*scalar, *workTensor, *dsdx, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

// ContravariantWeightBoundary_MappedScalar2D_gpu
__global__ void ContravariantWeightBoundary_MappedScalar2D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nEl){

  size_t iSide = blockIdx.x+1;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t iVar = threadIdx.y;

    workTensor[TEB_2D_INDEX(1,1,i,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_2D_INDEX(1,1,i,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_2D_INDEX(i,iSide,iEl,iVar,N,nEl)]; 
  
    workTensor[TEB_2D_INDEX(2,1,i,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_2D_INDEX(1,2,i,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_2D_INDEX(i,iSide,iEl,iVar,N,nEl)]; 
  
    workTensor[TEB_2D_INDEX(1,2,i,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_2D_INDEX(2,1,i,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_2D_INDEX(i,iSide,iEl,iVar,N,nEl)]; 
  
    workTensor[TEB_2D_INDEX(2,2,i,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_2D_INDEX(2,2,i,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_2D_INDEX(i,iSide,iEl,iVar,N,nEl)]; 
}

extern "C"
{
  void ContravariantWeightBoundary_MappedScalar2D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantWeightBoundary_MappedScalar2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*scalar, *workTensor, *dsdx, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

// ContravariantWeight_MappedScalar3D_gpu
__global__ void ContravariantWeight_MappedScalar3D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  workTensor[TE_3D_INDEX(1,1,i,j,k,iEl,iVar,N,nEl)] = dsdx[TE_3D_INDEX(1,1,i,j,k,iEl,0,N,nEl)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iEl,iVar,N,nEl)]; 

  workTensor[TE_3D_INDEX(2,1,i,j,k,iEl,iVar,N,nEl)] = dsdx[TE_3D_INDEX(1,2,i,j,k,iEl,0,N,nEl)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iEl,iVar,N,nEl)]; 

  workTensor[TE_3D_INDEX(3,1,i,j,k,iEl,iVar,N,nEl)] = dsdx[TE_3D_INDEX(1,3,i,j,k,iEl,0,N,nEl)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iEl,iVar,N,nEl)]; 

  workTensor[TE_3D_INDEX(1,2,i,j,k,iEl,iVar,N,nEl)] = dsdx[TE_3D_INDEX(2,1,i,j,k,iEl,0,N,nEl)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iEl,iVar,N,nEl)]; 

  workTensor[TE_3D_INDEX(2,2,i,j,k,iEl,iVar,N,nEl)] = dsdx[TE_3D_INDEX(2,2,i,j,k,iEl,0,N,nEl)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iEl,iVar,N,nEl)]; 

  workTensor[TE_3D_INDEX(3,2,i,j,k,iEl,iVar,N,nEl)] = dsdx[TE_3D_INDEX(2,3,i,j,k,iEl,0,N,nEl)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iEl,iVar,N,nEl)]; 

  workTensor[TE_3D_INDEX(1,3,i,j,k,iEl,iVar,N,nEl)] = dsdx[TE_3D_INDEX(3,1,i,j,k,iEl,0,N,nEl)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iEl,iVar,N,nEl)]; 

  workTensor[TE_3D_INDEX(2,3,i,j,k,iEl,iVar,N,nEl)] = dsdx[TE_3D_INDEX(3,2,i,j,k,iEl,0,N,nEl)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iEl,iVar,N,nEl)]; 

  workTensor[TE_3D_INDEX(3,3,i,j,k,iEl,iVar,N,nEl)] = dsdx[TE_3D_INDEX(3,3,i,j,k,iEl,0,N,nEl)]*
                                                       scalar[SC_3D_INDEX(i,j,k,iEl,iVar,N,nEl)]; 

}

extern "C"
{
  void ContravariantWeight_MappedScalar3D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantWeight_MappedScalar3D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,N+1), 0, 0>>>(*scalar, *workTensor, *dsdx, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

// ContravariantWeightBoundary_MappedScalar3D_gpu
__global__ void ContravariantWeightBoundary_MappedScalar3D_gpu(real *scalar, real *workTensor, real *dsdx, int N, int nEl){

  size_t iSide = blockIdx.x+1;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t iVar = threadIdx.z;

  workTensor[TEB_3D_INDEX(1,1,i,j,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_3D_INDEX(1,1,i,j,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_3D_INDEX(i,j,iSide,iEl,iVar,N,nEl)]; 

  workTensor[TEB_3D_INDEX(2,1,i,j,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_3D_INDEX(1,2,i,j,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_3D_INDEX(i,j,iSide,iEl,iVar,N,nEl)]; 

  workTensor[TEB_3D_INDEX(3,1,i,j,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_3D_INDEX(1,3,i,j,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_3D_INDEX(i,j,iSide,iEl,iVar,N,nEl)]; 

  workTensor[TEB_3D_INDEX(1,2,i,j,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_3D_INDEX(2,1,i,j,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_3D_INDEX(i,j,iSide,iEl,iVar,N,nEl)]; 

  workTensor[TEB_3D_INDEX(2,2,i,j,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_3D_INDEX(2,2,i,j,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_3D_INDEX(i,j,iSide,iEl,iVar,N,nEl)]; 

  workTensor[TEB_3D_INDEX(3,2,i,j,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_3D_INDEX(2,3,i,j,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_3D_INDEX(i,j,iSide,iEl,iVar,N,nEl)]; 

  workTensor[TEB_3D_INDEX(1,3,i,j,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_3D_INDEX(3,1,i,j,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_3D_INDEX(i,j,iSide,iEl,iVar,N,nEl)]; 

  workTensor[TEB_3D_INDEX(2,3,i,j,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_3D_INDEX(3,2,i,j,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_3D_INDEX(i,j,iSide,iEl,iVar,N,nEl)]; 

  workTensor[TEB_3D_INDEX(3,3,i,j,iSide,iEl,iVar,N,nEl)] = dsdx[TEB_3D_INDEX(3,3,i,j,iSide,iEl,0,N,nEl)]*
                                                       scalar[SCB_3D_INDEX(i,j,iSide,iEl,iVar,N,nEl)]; 

}

extern "C"
{
  void ContravariantWeightBoundary_MappedScalar3D_gpu_wrapper(real **scalar, real **workTensor, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantWeightBoundary_MappedScalar3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nEl), 0, 0>>>(*scalar, *workTensor, *dsdx, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

// ContravariantProjection_MappedVector2D_gpu
__global__ void ContravariantProjection_MappedVector2D_gpu(real *vector, real *dsdx, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;


    real Fx = vector[VE_2D_INDEX(1,i,j,iEl,iVar,N,nEl)];
    real Fy = vector[VE_2D_INDEX(2,i,j,iEl,iVar,N,nEl)];

    vector[VE_2D_INDEX(1,i,j,iEl,iVar,N,nEl)] = dsdx[TE_2D_INDEX(1,1,i,j,iEl,0,N,nEl)]*Fx+
                			         dsdx[TE_2D_INDEX(2,1,i,j,iEl,0,N,nEl)]*Fy;

    vector[VE_2D_INDEX(2,i,j,iEl,iVar,N,nEl)] = dsdx[TE_2D_INDEX(1,2,i,j,iEl,0,N,nEl)]*Fx+
						 dsdx[TE_2D_INDEX(2,2,i,j,iEl,0,N,nEl)]*Fy;
  
}

extern "C"
{
  void ContravariantProjection_MappedVector2D_gpu_wrapper(real **vector, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantProjection_MappedVector2D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(*vector, *dsdx, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

// ContravariantProjection_MappedP2Vector2D_gpu
__global__ void ContravariantProjection_MappedP2Vector2D_gpu(real *vector, real *physical, real *dsdx, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;


  for( int n = 0; n < N+1; n++){

    real Fx = physical[P2PVE_2D_INDEX(1,1,n,i,j,iEl,iVar,N,nEl)];
    real Fy = physical[P2PVE_2D_INDEX(2,1,n,i,j,iEl,iVar,N,nEl)];

    vector[P2VE_2D_INDEX(1,n,i,j,iEl,iVar,N,nEl)] = 
	    0.5*( dsdx[TE_2D_INDEX(1,1,i,j,iEl,0,N,nEl)] +
	          dsdx[TE_2D_INDEX(1,1,n,j,iEl,0,N,nEl)] )*Fx+
            0.5*( dsdx[TE_2D_INDEX(2,1,i,j,iEl,0,N,nEl)] +
                  dsdx[TE_2D_INDEX(2,1,n,j,iEl,0,N,nEl)] )*Fy;

    Fx = physical[P2PVE_2D_INDEX(1,2,n,i,j,iEl,iVar,N,nEl)];
    Fy = physical[P2PVE_2D_INDEX(2,2,n,i,j,iEl,iVar,N,nEl)];

    vector[P2VE_2D_INDEX(2,n,i,j,iEl,iVar,N,nEl)] = 
	    0.5*( dsdx[TE_2D_INDEX(1,2,i,j,iEl,0,N,nEl)] +
	          dsdx[TE_2D_INDEX(1,2,i,n,iEl,0,N,nEl)] )*Fx+
	    0.5*( dsdx[TE_2D_INDEX(2,2,i,j,iEl,0,N,nEl)] +
	          dsdx[TE_2D_INDEX(2,2,i,n,iEl,0,N,nEl)] )*Fy;
  }
  
}

extern "C"
{
  void ContravariantProjection_MappedP2Vector2D_gpu_wrapper(real **vector, real **physical, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantProjection_MappedP2Vector2D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(*vector, *physical, *dsdx, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

__global__ void ContravariantProjectionBoundary_MappedVector2D_gpu(real *physVector, real *boundaryNormal, real *nhat, int N, int nEl){

  size_t iSide = blockIdx.x+1;
  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.z;
  size_t i = threadIdx.x;

    boundaryNormal[SCB_2D_INDEX(i,iSide,iEl,iVar,N,nEl)] = nhat[VEB_2D_INDEX(1,i,iSide,iEl,0,N,nEl)]*
                                                     physVector[VEB_2D_INDEX(1,i,iSide,iEl,iVar,N,nEl)]+ 
						     nhat[VEB_2D_INDEX(1,i,iSide,iEl,0,N,nEl)]*
                                                     physVector[VEB_2D_INDEX(2,i,iSide,iEl,iVar,N,nEl)];
}

extern "C"
{
  void ContravariantProjectionBoundary_MappedVector2D_gpu_wrapper(real **physVector, real **boundaryNormal, real **nhat, int N, int nVar, int nEl)
  {
    ContravariantProjectionBoundary_MappedVector2D_gpu<<<dim3(4,nEl,nVar), dim3(N+1,1,1), 0, 0>>>(*physVector, *boundaryNormal, *nhat, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

// JacobianWeight_MappedVector2D_gpu
__global__ void JacobianWeight_MappedVector2D_gpu(real *vector, real *jacobian, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    vector[VE_2D_INDEX(1,i,j,iEl,iVar,N,nEl)] = vector[VE_2D_INDEX(1,i,j,iEl,iVar,N,nEl)]/
                                                 jacobian[SC_2D_INDEX(i,j,iEl,0,N,nEl)];

    vector[VE_2D_INDEX(2,i,j,iEl,iVar,N,nEl)] = vector[VE_2D_INDEX(2,i,j,iEl,iVar,N,nEl)]/
                                                 jacobian[SC_2D_INDEX(i,j,iEl,0,N,nEl)];
}

extern "C"
{
  void JacobianWeight_MappedVector2D_gpu_wrapper(real **vector, real **jacobian, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedVector2D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(*vector, *jacobian, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// ContravariantProjection_MappedVector3D_gpu
__global__ void ContravariantProjection_MappedVector3D_gpu(real *physVector, real *compVector, real *dsdx, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  real Fx = physVector[VE_3D_INDEX(1,i,j,k,iEl,iVar,N,nEl)];
  real Fy = physVector[VE_3D_INDEX(2,i,j,k,iEl,iVar,N,nEl)];
  real Fz = physVector[VE_3D_INDEX(3,i,j,k,iEl,iVar,N,nEl)];

  compVector[VE_3D_INDEX(1,i,j,k,iEl,iVar,N,nEl)] = dsdx[TE_3D_INDEX(1,1,i,j,k,iEl,0,N,nEl)]*Fx+
                                                     dsdx[TE_3D_INDEX(2,1,i,j,k,iEl,0,N,nEl)]*Fy+
                                                     dsdx[TE_3D_INDEX(3,1,i,j,k,iEl,0,N,nEl)]*Fz;

  compVector[VE_3D_INDEX(2,i,j,k,iEl,iVar,N,nEl)] = dsdx[TE_3D_INDEX(1,2,i,j,k,iEl,0,N,nEl)]*Fx+
                                                     dsdx[TE_3D_INDEX(2,2,i,j,k,iEl,0,N,nEl)]*Fy+
                                                     dsdx[TE_3D_INDEX(3,2,i,j,k,iEl,0,N,nEl)]*Fz;

  compVector[VE_3D_INDEX(3,i,j,k,iEl,iVar,N,nEl)] = dsdx[TE_3D_INDEX(1,3,i,j,k,iEl,0,N,nEl)]*Fx+
                                                     dsdx[TE_3D_INDEX(2,3,i,j,k,iEl,0,N,nEl)]*Fy+
                                                     dsdx[TE_3D_INDEX(3,3,i,j,k,iEl,0,N,nEl)]*Fz;

}

extern "C"
{
  void ContravariantProjection_MappedVector3D_gpu_wrapper(real **physVector, real **compVector, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantProjection_MappedVector3D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,N+1), 0, 0>>>(*physVector, *compVector, *dsdx, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

__global__ void ContravariantProjectionBoundary_MappedVector3D_gpu(real *physVector, real *boundaryNormal, real *nHat, int N, int nEl){

  size_t iSide = blockIdx.x+1;
  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.z;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  boundaryNormal[SCB_3D_INDEX(i,j,iSide,iEl,iVar,N,nEl)] = nHat[VEB_3D_INDEX(1,i,j,iSide,iEl,0,N,nEl)]*
                                                     physVector[VEB_3D_INDEX(1,i,j,iSide,iEl,iVar,N,nEl)]+ 
                                                     nHat[VEB_3D_INDEX(2,i,j,iSide,iEl,0,N,nEl)]*
                                                     physVector[VEB_3D_INDEX(2,i,j,iSide,iEl,iVar,N,nEl)]+
                                                     nHat[VEB_3D_INDEX(3,i,j,iSide,iEl,0,N,nEl)]*
                                                     physVector[VEB_3D_INDEX(3,i,j,iSide,iEl,iVar,N,nEl)];

}

extern "C"
{
  void ContravariantProjectionBoundary_MappedVector3D_gpu_wrapper(real **physVector, real **boundaryNormal, real **nHat, int N, int nVar, int nEl)
  {
    ContravariantProjectionBoundary_MappedVector3D_gpu<<<dim3(6,nEl,nVar), dim3(N+1,N+1,1), 0, 0>>>(*physVector, *boundaryNormal, *nHat, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

// JacobianWeight_MappedVector3D_gpu
__global__ void JacobianWeight_MappedVector3D_gpu(real *vector, real *jacobian, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    vector[VE_3D_INDEX(1,i,j,k,iEl,iVar,N,nEl)] = vector[VE_3D_INDEX(1,i,j,k,iEl,iVar,N,nEl)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];

    vector[VE_3D_INDEX(2,i,j,k,iEl,iVar,N,nEl)] = vector[VE_3D_INDEX(2,i,j,k,iEl,iVar,N,nEl)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];

    vector[VE_3D_INDEX(3,i,j,k,iEl,iVar,N,nEl)] = vector[VE_3D_INDEX(3,i,j,k,iEl,iVar,N,nEl)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];
}

extern "C"
{
  void JacobianWeight_MappedVector3D_gpu_wrapper(real **vector, real **jacobian, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedVector3D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,N+1), 0, 0>>>(*vector, *jacobian, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// JacobianWeight_MappedTensor2D_gpu
__global__ void JacobianWeight_MappedTensor2D_gpu(real *tensor, real *jacobian, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    tensor[TE_2D_INDEX(1,1,i,j,iEl,iVar,N,nEl)] = tensor[TE_2D_INDEX(1,1,i,j,iEl,iVar,N,nEl)]/
                                                 jacobian[SC_2D_INDEX(i,j,iEl,0,N,nEl)];

    tensor[TE_2D_INDEX(2,1,i,j,iEl,iVar,N,nEl)] = tensor[TE_2D_INDEX(2,1,i,j,iEl,iVar,N,nEl)]/
                                                 jacobian[SC_2D_INDEX(i,j,iEl,0,N,nEl)];

    tensor[TE_2D_INDEX(1,2,i,j,iEl,iVar,N,nEl)] = tensor[TE_2D_INDEX(1,2,i,j,iEl,iVar,N,nEl)]/
                                                 jacobian[SC_2D_INDEX(i,j,iEl,0,N,nEl)];

    tensor[TE_2D_INDEX(2,2,i,j,iEl,iVar,N,nEl)] = tensor[TE_2D_INDEX(2,2,i,j,iEl,iVar,N,nEl)]/
                                                 jacobian[SC_2D_INDEX(i,j,iEl,0,N,nEl)];
}

extern "C"
{
  void JacobianWeight_MappedTensor2D_gpu_wrapper(real **tensor, real **jacobian, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedTensor2D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(*tensor, *jacobian, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// JacobianWeight_MappedTensor3D_gpu
__global__ void JacobianWeight_MappedTensor3D_gpu(real *tensor, real *jacobian, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    tensor[TE_3D_INDEX(1,1,i,j,k,iEl,iVar,N,nEl)] = tensor[TE_3D_INDEX(1,1,i,j,k,iEl,iVar,N,nEl)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];

    tensor[TE_3D_INDEX(2,1,i,j,k,iEl,iVar,N,nEl)] = tensor[TE_3D_INDEX(2,1,i,j,k,iEl,iVar,N,nEl)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];

    tensor[TE_3D_INDEX(3,1,i,j,k,iEl,iVar,N,nEl)] = tensor[TE_3D_INDEX(3,1,i,j,k,iEl,iVar,N,nEl)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];

    tensor[TE_3D_INDEX(1,2,i,j,k,iEl,iVar,N,nEl)] = tensor[TE_3D_INDEX(1,2,i,j,k,iEl,iVar,N,nEl)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];

    tensor[TE_3D_INDEX(2,2,i,j,k,iEl,iVar,N,nEl)] = tensor[TE_3D_INDEX(2,2,i,j,k,iEl,iVar,N,nEl)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];

    tensor[TE_3D_INDEX(3,2,i,j,k,iEl,iVar,N,nEl)] = tensor[TE_3D_INDEX(3,2,i,j,k,iEl,iVar,N,nEl)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];

    tensor[TE_3D_INDEX(1,3,i,j,k,iEl,iVar,N,nEl)] = tensor[TE_3D_INDEX(1,3,i,j,k,iEl,iVar,N,nEl)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];

    tensor[TE_3D_INDEX(2,3,i,j,k,iEl,iVar,N,nEl)] = tensor[TE_3D_INDEX(2,3,i,j,k,iEl,iVar,N,nEl)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];

    tensor[TE_3D_INDEX(3,3,i,j,k,iEl,iVar,N,nEl)] = tensor[TE_3D_INDEX(3,3,i,j,k,iEl,iVar,N,nEl)]/
                                                   jacobian[SC_3D_INDEX(i,j,k,iEl,0,N,nEl)];
}

extern "C"
{
  void JacobianWeight_MappedTensor3D_gpu_wrapper(real **tensor, real **jacobian, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedTensor3D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,N+1), 0, 0>>>(*tensor, *jacobian, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// CalculateCurl_MappedTensor2D_gpu
__global__ void CalculateCurl_MappedTensor2D_gpu(real *dfdx, real *curlf, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    curlf[SC_2D_INDEX(i,j,iEl,iVar,N,nEl)] = dfdx[TE_2D_INDEX(2,1,i,j,iEl,iVar,N,nEl)]-
                                              dfdx[TE_2D_INDEX(1,2,i,j,iEl,iVar,N,nEl)];
}

extern "C"
{
  void CalculateCurl_MappedTensor2D_gpu_wrapper(real **dfdx, real **curlf, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedTensor2D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(*dfdx, *curlf, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// CalculateCurl_MappedTensor3D_gpu
__global__ void CalculateCurl_MappedTensor3D_gpu(real *dfdx, real *curlf, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    curlf[VE_3D_INDEX(1,i,j,k,iEl,iVar,N,nEl)] = dfdx[TE_3D_INDEX(3,2,i,j,k,iEl,iVar,N,nEl)]-
                                                  dfdx[TE_3D_INDEX(2,3,i,j,k,iEl,iVar,N,nEl)];

    curlf[VE_3D_INDEX(2,i,j,k,iEl,iVar,N,nEl)] = dfdx[TE_3D_INDEX(1,3,i,j,k,iEl,iVar,N,nEl)]-
                                                  dfdx[TE_3D_INDEX(3,1,i,j,k,iEl,iVar,N,nEl)];

    curlf[VE_3D_INDEX(3,i,j,k,iEl,iVar,N,nEl)] = dfdx[TE_3D_INDEX(2,1,i,j,k,iEl,iVar,N,nEl)]-
                                                  dfdx[TE_3D_INDEX(1,2,i,j,k,iEl,iVar,N,nEl)];
}

extern "C"
{
  void CalculateCurl_MappedTensor3D_gpu_wrapper(real **dfdx, real **curlf, int N, int nVar, int nEl)
  {
    JacobianWeight_MappedTensor3D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,N+1), 0, 0>>>(*dfdx, *curlf, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// MapTo support routines
__global__ void MapToScalar_MappedVector2D_gpu(real* scalar, real* vector, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    size_t jVar = 2*(iVar);
    scalar[SC_2D_INDEX(i,j,jVar,iEl,N,nVar*2)] = vector[VE_2D_INDEX(1,i,j,iEl,iVar,N,nEl)];

    jVar += 1;
    scalar[SC_2D_INDEX(i,j,jVar,iEl,N,nVar*2)] = vector[VE_2D_INDEX(2,i,j,iEl,iVar,N,nEl)];

}

extern "C"
{
  void MapToScalar_MappedVector2D_gpu_wrapper(real **scalar, real **vector, int N, int nVar, int nEl)
  {
    MapToScalar_MappedVector2D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(*scalar, *vector, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void MapToScalarBoundary_MappedVector2D_gpu(real* scalar, real* vector, int N, int nEl){

  size_t iSide = blockIdx.x+1;
  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.z;
  size_t j = threadIdx.x;

    size_t jVar = 2*(iVar);
    scalar[SCB_2D_INDEX(j,jVar,iSide,iEl,N,nVar*2)] = vector[VEB_2D_INDEX(1,j,iSide,iEl,iVar,N,nEl)];

    jVar += 1;
    scalar[SCB_2D_INDEX(j,jVar,iSide,iEl,N,nVar*2)] = vector[VEB_2D_INDEX(2,j,iSide,iEl,iVar,N,nEl)];

}

extern "C"
{
  void MapToScalarBoundary_MappedVector2D_gpu_wrapper(real **scalar, real **vector, int N, int nVar, int nEl)
  {
    MapToScalarBoundary_MappedVector2D_gpu<<<dim3(4,nEl,nVar), dim3(N+1,1,1), 0, 0>>>(*scalar, *vector, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void MapToTensor_MappedVector2D_gpu(real* tensor, real* vector, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    size_t jVar = 2*(iVar);
    tensor[TE_2D_INDEX(1,1,i,j,iEl,iVar,N,nEl)] = vector[VE_2D_INDEX(1,i,j,jVar,iEl,N,nVar*2)];

    tensor[TE_2D_INDEX(1,2,i,j,iEl,iVar,N,nEl)] = vector[VE_2D_INDEX(2,i,j,jVar,iEl,N,nVar*2)];

    jVar +=1 ;
    tensor[TE_2D_INDEX(2,1,i,j,iEl,iVar,N,nEl)] = vector[VE_2D_INDEX(1,i,j,jVar,iEl,N,nVar*2)];

    tensor[TE_2D_INDEX(2,2,i,j,iEl,iVar,N,nEl)] = vector[VE_2D_INDEX(2,i,j,jVar,iEl,N,nVar*2)];

}

extern "C"
{
  void MapToTensor_MappedVector2D_gpu_wrapper(real **tensor, real **vector, int N, int nVar, int nEl)
  {
    MapToTensor_MappedVector2D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,1), 0, 0>>>(*tensor, *vector, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void MapToTensorBoundary_MappedVector2D_gpu(real* tensor, real* vector, int N, int nEl){

  size_t iSide = blockIdx.x+1;
  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.z;
  size_t j = threadIdx.x;

    size_t jVar =2*(iVar);
    tensor[TEB_2D_INDEX(1,1,j,iSide,iEl,iVar,N,nEl)] = vector[VEB_2D_INDEX(1,j,jVar,iSide,iEl,N,nVar*2)];

    tensor[TEB_2D_INDEX(1,2,j,iSide,iEl,iVar,N,nEl)] = vector[VEB_2D_INDEX(2,j,jVar,iSide,iEl,N,nVar*2)];

    jVar += 1;
    tensor[TEB_2D_INDEX(2,1,j,iSide,iEl,iVar,N,nEl)] = vector[VEB_2D_INDEX(1,j,jVar,iSide,iEl,N,nVar*2)];

    tensor[TEB_2D_INDEX(2,2,j,iSide,iEl,iVar,N,nEl)] = vector[VEB_2D_INDEX(2,j,jVar,iSide,iEl,N,nVar*2)];

}

extern "C"
{
  void MapToTensorBoundary_MappedVector2D_gpu_wrapper(real **tensor, real **vector, int N, int nVar, int nEl)
  {
    MapToTensorBoundary_MappedVector2D_gpu<<<dim3(4,nEl,nVar), dim3(N+1,1,1), 0, 0>>>(*tensor, *vector, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void MapToScalar_MappedVector3D_gpu(real* scalar, real* vector, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    size_t jVar = 3*(iVar);
    scalar[SC_3D_INDEX(i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(1,i,j,k,iEl,iVar,N,nEl)];

    jVar += 1;
    scalar[SC_3D_INDEX(i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(2,i,j,k,iEl,iVar,N,nEl)];

    jVar += 1;
    scalar[SC_3D_INDEX(i,j,k,jVar,iEl,N,nVar*3)] = vector[VE_3D_INDEX(3,i,j,k,iEl,iVar,N,nEl)];

}

extern "C"
{
  void MapToScalar_MappedVector3D_gpu_wrapper(real **scalar, real **vector, int N, int nVar, int nEl)
  {
    MapToScalar_MappedVector3D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,N+1), 0, 0>>>(*scalar, *vector, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void MapToScalarBoundary_MappedVector3D_gpu(real* scalar, real* vector, int N, int nEl){

  size_t iSide = blockIdx.x+1;
  size_t iEl = blockIdx.x;
  size_t j = threadIdx.x;
  size_t k = threadIdx.y;
  size_t iVar = threadIdx.z;

    size_t jVar = 3*(iVar);
    scalar[SCB_3D_INDEX(j,k,jVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(1,j,k,iSide,iEl,iVar,N,nEl)];

    jVar += 1;
    scalar[SCB_3D_INDEX(j,k,jVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(2,j,k,iSide,iEl,iVar,N,nEl)];

    jVar += 1;
    scalar[SCB_3D_INDEX(j,k,jVar,iSide,iEl,N,nVar*3)] = vector[VEB_3D_INDEX(3,j,k,iSide,iEl,iVar,N,nEl)];

}

extern "C"
{
  void MapToScalarBoundary_MappedVector3D_gpu_wrapper(real **scalar, real **vector, int N, int nVar, int nEl)
  {
    MapToScalarBoundary_MappedVector3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nEl), 0, 0>>>(*scalar, *vector, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void MapToTensor_MappedVector3D_gpu(real* tensor, real* vector, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t iEl = blockIdx.x;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

    size_t jVar = 3*(iVar);
    tensor[TE_3D_INDEX(1,1,i,j,k,iEl,iVar,N,nEl)] = vector[VE_3D_INDEX(1,i,j,k,jVar,iEl,N,nVar*3)];

    tensor[TE_3D_INDEX(1,2,i,j,k,iEl,iVar,N,nEl)] = vector[VE_3D_INDEX(2,i,j,k,jVar,iEl,N,nVar*3)];

    tensor[TE_3D_INDEX(1,3,i,j,k,iEl,iVar,N,nEl)] = vector[VE_3D_INDEX(3,i,j,k,jVar,iEl,N,nVar*3)];

    jVar += 1;
    tensor[TE_3D_INDEX(2,1,i,j,k,iEl,iVar,N,nEl)] = vector[VE_3D_INDEX(1,i,j,k,jVar,iEl,N,nVar*3)];

    tensor[TE_3D_INDEX(2,2,i,j,k,iEl,iVar,N,nEl)] = vector[VE_3D_INDEX(2,i,j,k,jVar,iEl,N,nVar*3)];

    tensor[TE_3D_INDEX(2,3,i,j,k,iEl,iVar,N,nEl)] = vector[VE_3D_INDEX(3,i,j,k,jVar,iEl,N,nVar*3)];

    jVar += 1;
    tensor[TE_3D_INDEX(3,1,i,j,k,iEl,iVar,N,nEl)] = vector[VE_3D_INDEX(1,i,j,k,jVar,iEl,N,nVar*3)];

    tensor[TE_3D_INDEX(3,2,i,j,k,iEl,iVar,N,nEl)] = vector[VE_3D_INDEX(2,i,j,k,jVar,iEl,N,nVar*3)];

    tensor[TE_3D_INDEX(3,3,i,j,k,iEl,iVar,N,nEl)] = vector[VE_3D_INDEX(3,i,j,k,jVar,iEl,N,nVar*3)];

}

extern "C"
{
  void MapToTensor_MappedVector3D_gpu_wrapper(real **tensor, real **vector, int N, int nVar, int nEl)
  {
    MapToTensor_MappedVector3D_gpu<<<dim3(nEl,nVar,1), dim3(N+1,N+1,N+1), 0, 0>>>(*tensor, *vector, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void MapToTensorBoundary_MappedVector3D_gpu(real* tensor, real* vector, int N, int nEl){

  size_t iSide = blockIdx.x+1;
  size_t iEl = blockIdx.x;
  size_t j = threadIdx.x;
  size_t k = threadIdx.x;
  size_t iVar = threadIdx.z;

    size_t jVar = 3*(iVar);
    tensor[TEB_3D_INDEX(1,1,j,k,iSide,iEl,iVar,N,nEl)] = vector[VEB_3D_INDEX(1,j,k,jVar,iSide,iEl,N,nVar*3)];

    tensor[TEB_3D_INDEX(1,2,j,k,iSide,iEl,iVar,N,nEl)] = vector[VEB_3D_INDEX(2,j,k,jVar,iSide,iEl,N,nVar*3)];

    tensor[TEB_3D_INDEX(1,3,j,k,iSide,iEl,iVar,N,nEl)] = vector[VEB_3D_INDEX(3,j,k,jVar,iSide,iEl,N,nVar*3)];

    jVar += 1;
    tensor[TEB_3D_INDEX(2,1,j,k,iSide,iEl,iVar,N,nEl)] = vector[VEB_3D_INDEX(1,j,k,jVar,iSide,iEl,N,nVar*3)];

    tensor[TEB_3D_INDEX(2,2,j,k,iSide,iEl,iVar,N,nEl)] = vector[VEB_3D_INDEX(2,j,k,jVar,iSide,iEl,N,nVar*3)];

    tensor[TEB_3D_INDEX(2,3,j,k,iSide,iEl,iVar,N,nEl)] = vector[VEB_3D_INDEX(3,j,k,jVar,iSide,iEl,N,nVar*3)];

    jVar += 1;
    tensor[TEB_3D_INDEX(3,1,j,k,iSide,iEl,iVar,N,nEl)] = vector[VEB_3D_INDEX(1,j,k,jVar,iSide,iEl,N,nVar*3)];

    tensor[TEB_3D_INDEX(3,2,j,k,iSide,iEl,iVar,N,nEl)] = vector[VEB_3D_INDEX(2,j,k,jVar,iSide,iEl,N,nVar*3)];

    tensor[TEB_3D_INDEX(3,3,j,k,iSide,iEl,iVar,N,nEl)] = vector[VEB_3D_INDEX(3,j,k,jVar,iSide,iEl,N,nVar*3)];

}

extern "C"
{
  void MapToTensorBoundary_MappedVector3D_gpu_wrapper(real **tensor, real **vector, int N, int nVar, int nEl)
  {
    MapToTensorBoundary_MappedVector3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nEl), 0, 0>>>(*tensor, *vector, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void SideExchange_MappedScalar2D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nEl){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t ivar = threadIdx.y;
  
  int e2Global = sideInfo[INDEX3(2,s1-1,e1,5,4)]-1;
  int e2 = e2Global - offset;
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];

  if(bcid == 0){
    int neighborRank = elemToRank[e2Global];
    if( neighborRank == rankId ){
      if(flip == 0){
        extBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nEl)] = boundary[SCB_2D_INDEX(i1,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 1){
        int i2 = N-i1;
        extBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nEl)] = boundary[SCB_2D_INDEX(i2,ivar,s2,e2,N,nEl)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedScalar2D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    dim3 nblocks(4,nEl,1);
    dim3 nthreads(N+1,nVar,1);
    SideExchange_MappedScalar2D_gpu<<<nblocks,nthreads>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, offset, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void SideExchange_MappedVector2D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nEl){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t ivar = threadIdx.y;
  
  int e2Global = sideInfo[INDEX3(2,s1-1,e1,5,4)]-1;
  int e2 = e2Global - offset;
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];
  int i2 = N-i1;

  if(bcid == 0){
    int neighborRank = elemToRank[e2Global];
    if( neighborRank == rankId ){
      if(flip == 0){
        extBoundary[VEB_2D_INDEX(1,i1,ivar,s1,e1,N,nEl)] = boundary[VEB_2D_INDEX(1,i1,ivar,s2,e2,N,nEl)];
        extBoundary[VEB_2D_INDEX(2,i1,ivar,s1,e1,N,nEl)] = boundary[VEB_2D_INDEX(2,i1,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 1){
        extBoundary[VEB_2D_INDEX(1,i1,ivar,s1,e1,N,nEl)] = boundary[VEB_2D_INDEX(1,i2,ivar,s2,e2,N,nEl)];
        extBoundary[VEB_2D_INDEX(2,i1,ivar,s1,e1,N,nEl)] = boundary[VEB_2D_INDEX(2,i2,ivar,s2,e2,N,nEl)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedVector2D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    SideExchange_MappedVector2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, offset, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void SideExchange_MappedTensor2D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nEl){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t ivar = threadIdx.y;
  
  int e2Global = sideInfo[INDEX3(2,s1-1,e1,5,4)]-1;
  int e2 = e2Global - offset;
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];
  int i2 = N-i1;

  if(bcid == 0){
    int neighborRank = elemToRank[e2Global];
    if( neighborRank == rankId ){
      if(flip == 0){
        extBoundary[TEB_2D_INDEX(1,1,i1,ivar,s1,e1,N,nEl)] = boundary[TEB_2D_INDEX(1,1,i1,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_2D_INDEX(2,1,i1,ivar,s1,e1,N,nEl)] = boundary[TEB_2D_INDEX(2,1,i1,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_2D_INDEX(1,2,i1,ivar,s1,e1,N,nEl)] = boundary[TEB_2D_INDEX(1,2,i1,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_2D_INDEX(2,2,i1,ivar,s1,e1,N,nEl)] = boundary[TEB_2D_INDEX(2,2,i1,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 1){
        extBoundary[TEB_2D_INDEX(1,1,i1,ivar,s1,e1,N,nEl)] = boundary[TEB_2D_INDEX(1,1,i2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_2D_INDEX(2,1,i1,ivar,s1,e1,N,nEl)] = boundary[TEB_2D_INDEX(2,1,i2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_2D_INDEX(1,2,i1,ivar,s1,e1,N,nEl)] = boundary[TEB_2D_INDEX(1,2,i2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_2D_INDEX(2,2,i1,ivar,s1,e1,N,nEl)] = boundary[TEB_2D_INDEX(2,2,i2,ivar,s2,e2,N,nEl)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedTensor2D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    SideExchange_MappedTensor2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, offset, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void SideExchange_MappedScalar3D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nEl){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  size_t ivar = threadIdx.z;
  
  int e2Global = sideInfo[INDEX3(2,s1-1,e1,5,6)]-1;
  int e2 = e2Global - offset;
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,6)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,6)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,6)];

  if(bcid == 0){
    int neighborRank = elemToRank[e2Global];
    if( neighborRank == rankId ){
      if(flip == 0){
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)] = boundary[SCB_3D_INDEX(i1,j1,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 1){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)] = boundary[SCB_3D_INDEX(i2,j2,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 2){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)] = boundary[SCB_3D_INDEX(i2,j2,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 3){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)] = boundary[SCB_3D_INDEX(i2,j2,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)] = boundary[SCB_3D_INDEX(i2,j2,ivar,s2,e2,N,nEl)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedScalar3D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    SideExchange_MappedScalar3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nEl), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, offset, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void SideExchange_MappedVector3D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nEl){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  size_t ivar = threadIdx.z;
  
  int e2Global = sideInfo[INDEX3(2,s1-1,e1,5,6)]-1;
  int e2 = e2Global - offset;
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,6)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,6)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,6)];

  if(bcid == 0){
    int neighborRank = elemToRank[e2Global];
    if( neighborRank == rankId ){
      if(flip == 0){
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(1,i1,j1,ivar,s2,e2,N,nEl)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(2,i1,j1,ivar,s2,e2,N,nEl)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(3,i1,j1,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 1){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(3,i2,j2,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 2){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(3,i2,j2,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 3){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(3,i2,j2,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = i1;
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[VEB_3D_INDEX(3,i2,j2,ivar,s2,e2,N,nEl)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedVector3D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    SideExchange_MappedVector3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nEl), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, offset, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void SideExchange_MappedTensor3D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nEl){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  size_t ivar = threadIdx.z;
  
  int e2Global = sideInfo[INDEX3(2,s1-1,e1,5,6)]-1;
  int e2 = e2Global - offset;
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,6)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,6)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,6)];

  if(bcid == 0){
    int neighborRank = elemToRank[e2Global];
    if( neighborRank == rankId ){
      if(flip == 0){
        extBoundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nEl)];
        extBoundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 1){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,3,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,3,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,3,i2,j2,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 2){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,3,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,3,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,3,i2,j2,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 3){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,3,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,3,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,3,i2,j2,ivar,s2,e2,N,nEl)];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = i1;
        extBoundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,1,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(1,2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(2,2,i2,j2,ivar,s2,e2,N,nEl)];
        extBoundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nEl)] = boundary[TEB_3D_INDEX(3,2,i2,j2,ivar,s2,e2,N,nEl)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedTensor3D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    SideExchange_MappedTensor3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nEl), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, offset, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void BassiRebaySides_MappedScalar2D_gpu(real* avgBoundary, real *boundary, real *extBoundary, int N, int nEl){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t ivar = threadIdx.y;
  
  avgBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nEl)] =0.5*(extBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nEl)]+
		                                     boundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nEl)]);
  
}

extern "C"
{
  void BassiRebaySides_MappedScalar2D_gpu_wrapper(real **avgBoundary, real **boundary, real **extBoundary, int N, int nVar, int nEl)
  {
    BassiRebaySides_MappedScalar2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*avgBoundary, *boundary, *extBoundary, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void BassiRebaySides_MappedVector2D_gpu(real *extBoundary, real *boundary, int N, int nEl){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t ivar = threadIdx.y;
  
  boundary[VEB_2D_INDEX(1,i1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[VEB_2D_INDEX(1,i1,ivar,s1,e1,N,nEl)]+
                                                   boundary[VEB_2D_INDEX(1,i1,ivar,s1,e1,N,nEl)]);
  boundary[VEB_2D_INDEX(2,i1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[VEB_2D_INDEX(2,i1,ivar,s1,e1,N,nEl)]+
                                                   boundary[VEB_2D_INDEX(2,i1,ivar,s1,e1,N,nEl)]);
  
}

extern "C"
{
  void BassiRebaySides_MappedVector2D_gpu_wrapper(real **extBoundary, real **boundary, int N, int nVar, int nEl)
  {
    BassiRebaySides_MappedVector2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*extBoundary, *boundary, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void BassiRebaySides_MappedTensor2D_gpu(real *extBoundary, real *boundary, int N, int nEl){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t ivar = threadIdx.y;
  
  boundary[TEB_2D_INDEX(1,1,i1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_2D_INDEX(1,1,i1,ivar,s1,e1,N,nEl)]+
                                                        boundary[TEB_2D_INDEX(1,1,i1,ivar,s1,e1,N,nEl)]);
  boundary[TEB_2D_INDEX(2,1,i1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_2D_INDEX(2,1,i1,ivar,s1,e1,N,nEl)]+
                                                        boundary[TEB_2D_INDEX(2,1,i1,ivar,s1,e1,N,nEl)]);
  boundary[TEB_2D_INDEX(1,2,i1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_2D_INDEX(1,2,i1,ivar,s1,e1,N,nEl)]+
                                                        boundary[TEB_2D_INDEX(1,2,i1,ivar,s1,e1,N,nEl)]);
  boundary[TEB_2D_INDEX(2,2,i1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_2D_INDEX(2,2,i1,ivar,s1,e1,N,nEl)]+
                                                        boundary[TEB_2D_INDEX(2,2,i1,ivar,s1,e1,N,nEl)]);
  
}

extern "C"
{
  void BassiRebaySides_MappedTensor2D_gpu_wrapper(real **extBoundary, real **boundary, int N, int nVar, int nEl)
  {
    BassiRebaySides_MappedTensor2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*extBoundary, *boundary, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void BassiRebaySides_MappedScalar3D_gpu(real *avgBoundary, real *boundary, real *extBoundary, int N, int nEl){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  size_t ivar = threadIdx.z;
  
  avgBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)]+
                                                         boundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)]);
  
}

extern "C"
{
  void BassiRebaySides_MappedScalar3D_gpu_wrapper(real **avgBoundary, real **boundary, real **extBoundary, int N, int nVar, int nEl)
  {
    BassiRebaySides_MappedScalar3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nEl), 0, 0>>>(*avgBoundary, *boundary, *extBoundary, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void BassiRebaySides_MappedVector3D_gpu(real *extBoundary, real *boundary, int N, int nEl){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  size_t ivar = threadIdx.z;
  
  boundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nEl)]+
                                                      boundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nEl)]);
  boundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nEl)]+
                                                      boundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nEl)]);
  boundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nEl)]+
                                                      boundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nEl)]);
  
}

extern "C"
{
  void BassiRebaySides_MappedVector3D_gpu_wrapper(real **extBoundary, real **boundary, int N, int nVar, int nEl)
  {
    BassiRebaySides_MappedVector3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nEl), 0, 0>>>(*extBoundary, *boundary, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void BassiRebaySides_MappedTensor3D_gpu(real *extBoundary, real *boundary, int N, int nEl){

  size_t s1 = blockIdx.x+1;
  size_t e1 = blockIdx.y;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.y;
  size_t ivar = threadIdx.z;
  
      boundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nEl)]+
                                                            boundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nEl)]);
      boundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nEl)]+
                                                            boundary[TEB_3D_INDEX(2,1,i1,j1,ivar,s1,e1,N,nEl)]);
      boundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nEl)]+
                                                            boundary[TEB_3D_INDEX(3,1,i1,j1,ivar,s1,e1,N,nEl)]);
      boundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_3D_INDEX(1,2,i1,j1,ivar,s1,e1,N,nEl)]+
                                                            boundary[TEB_3D_INDEX(1,1,i1,j1,ivar,s1,e1,N,nEl)]);
      boundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nEl)]+
                                                            boundary[TEB_3D_INDEX(2,2,i1,j1,ivar,s1,e1,N,nEl)]);
      boundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nEl)]+
                                                            boundary[TEB_3D_INDEX(3,2,i1,j1,ivar,s1,e1,N,nEl)]);
      boundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nEl)]+
                                                            boundary[TEB_3D_INDEX(1,3,i1,j1,ivar,s1,e1,N,nEl)]);
      boundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nEl)]+
                                                            boundary[TEB_3D_INDEX(2,3,i1,j1,ivar,s1,e1,N,nEl)]);
      boundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nEl)] = 0.5*(extBoundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nEl)]+
                                                            boundary[TEB_3D_INDEX(3,3,i1,j1,ivar,s1,e1,N,nEl)]);
  
}

extern "C"
{
  void BassiRebaySides_MappedTensor3D_gpu_wrapper(real **extBoundary, real **boundary, int N, int nVar, int nEl)
  {
    BassiRebaySides_MappedTensor3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nEl), 0, 0>>>(*extBoundary, *boundary, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void ApplyFlip_MappedScalar2D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t s1 = blockIdx.y+1;
  size_t e1 = blockIdx.z;
  size_t i1 = threadIdx.x;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];
  int i2 = N-i1;
  int neighborRank = elemToRank[e2];

  if(bcid == 0){ // Interior Element
    if( neighborRank /= rankId ){ // Side shared with another rank
      if(flip == 1){ // Neighboring elements have different orientation

        __shared__ real extBuff[16];
        extBuff[i1] = extBoundary[SCB_2D_INDEX(i2,ivar,s1,e1,N,nEl)];
        __syncthreads();

        extBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nEl)] = extBuff[i1];
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedScalar2D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedScalar2D_gpu<<<dim3(nVar,4,nEl), dim3(N+1,1,1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void ApplyFlip_MappedVector2D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t s1 = blockIdx.y+1;
  size_t e1 = blockIdx.z;
  size_t i1 = threadIdx.x;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];
  int i2 = N-i1;
  int neighborRank = elemToRank[e2];


  if(bcid == 0){
    if( neighborRank /= rankId ){
      if(flip == 1){

        __shared__ real extBuff[16];
        extBuff[i1] = extBoundary[VEB_2D_INDEX(1,i2,ivar,s1,e1,N,nEl)];
        __syncthreads();

        extBoundary[VEB_2D_INDEX(1,i1,ivar,s1,e1,N,nEl)] = extBuff[i1];

        extBuff[i1] = extBoundary[VEB_2D_INDEX(2,i2,ivar,s1,e1,N,nEl)];
        __syncthreads();

        extBoundary[VEB_2D_INDEX(2,i1,ivar,s1,e1,N,nEl)] = extBuff[i2];

      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedVector2D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedVector2D_gpu<<<dim3(nVar,4,nEl), dim3(N+1,1,1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void ApplyFlip_MappedTensor2D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t s1 = blockIdx.y+1;
  size_t e1 = blockIdx.z;
  size_t row = threadIdx.x;
  size_t col = threadIdx.y;
  size_t i1 = threadIdx.z;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];
  int i2 = N-i1;

  __shared__ real extBuff[64];

  extBuff[row+2*(col+2*i1)] = extBoundary[TEB_2D_INDEX(row+1,col+1,i1,ivar,s1,e1,N,nEl)];

  __syncthreads();

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank /= rankId ){
      if(flip == 1){
        extBoundary[TEB_2D_INDEX(row+1,col+1,i1,ivar,s1,e1,N,nEl)] = extBuff[row+2*(col+2*i2)];
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedTensor2D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedTensor2D_gpu<<<dim3(nVar,4,nEl), dim3(2,2,N+1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void ApplyFlip_MappedScalar3D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t s1 = blockIdx.y+1;
  size_t e1 = blockIdx.z;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.x;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];

  __shared__ real extBuff[256];

  extBuff[i1+(N+1)*j1] = extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)];

  __syncthreads();

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank /= rankId ){
      if(flip == 1){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 2){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 3){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nEl)] = extBuff[i2+(N+1)*j2];
      }
    }
  }
}

extern "C"
{
  void ApplyFlip_MappedScalar3D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedScalar3D_gpu<<<dim3(nVar,6,nEl), dim3(N+1,N+1,1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void ApplyFlip_MappedVector3D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t s1 = blockIdx.y+1;
  size_t e1 = blockIdx.z;
  size_t dir = threadIdx.x;
  size_t i1 = threadIdx.y;
  size_t j1 = threadIdx.z;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];

  __shared__ real extBuff[768];

  extBuff[dir+3*(i1+(N+1)*j1)] = extBoundary[VEB_3D_INDEX(dir+1,i1,j1,ivar,s1,e1,N,nEl)];

  __syncthreads();

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank /= rankId ){
      if(flip == 1){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nEl)] = extBuff[dir+3*(i2+(N+1)*j2)];
      }
      else if(flip == 2){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nEl)] = extBuff[dir+3*(i2+(N+1)*j2)];
      }
      else if(flip == 3){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nEl)] = extBuff[dir+3*(i2+(N+1)*j2)];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = i1;
        extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nEl)] = extBuff[dir+3*(i2+(N+1)*j2)];
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedVector3D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedVector3D_gpu<<<dim3(nVar,6,nEl), dim3(3,N+1,N+1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void ApplyFlip_MappedTensor3D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nEl){

  size_t iVar = blockIdx.y;
  size_t s1 = blockIdx.y+1;
  size_t e1 = blockIdx.z;
  size_t dir = threadIdx.x;
  size_t i1 = threadIdx.y;
  size_t j1 = threadIdx.z;
  size_t row = dir/3;
  size_t col = dir - dir*row;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];

  __shared__ real extBuff[2304];

  extBuff[row+3*(col+3*(i1+(N+1)*j1))] = extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nEl)];

  __syncthreads();

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank /= rankId ){
      if(flip == 1){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nEl)] = extBuff[row+3*(col+3*(i2+(N+1)*j2))];
      }
      else if(flip == 2){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nEl)] = extBuff[row+3*(col+3*(i2+(N+1)*j2))];
      }
      else if(flip == 3){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nEl)] = extBuff[row+3*(col+3*(i2+(N+1)*j2))];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = i1;
        extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nEl)] = extBuff[row+3*(col+3*(i2+(N+1)*j2))];
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedTensor3D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedTensor3D_gpu<<<dim3(nVar,6,nEl), dim3(9,N+1,N+1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nEl);
    HIP_SAFE_CALL(hipGetLastError());
  }

}
