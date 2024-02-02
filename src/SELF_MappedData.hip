#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

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
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// GradientInterior_MappedScalar2D_gpu
// Performs the interior operations for a scalar gradient calculation
// Use a derivative matrix for strong form;
// use a DG derivative matrix for weak form and follow up with the boundary contributions
__global__ void GradientInterior_MappedScalar2D_gpu(real *scalar, real *dsdx, real *jacobian, real *gradF, real *dMatrix, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real gf[2] = {0.0};

  for( int ii=0; ii < N+1; ii++ ){
    gf[0] += dMatrix[ii+i*(N+1)]*
	       scalar[SC_2D_INDEX(ii,j,iVar,iEl,N,nVar)]*
	       dsdx[TE_2D_INDEX(1,1,ii,j,0,iEl,N,1)]+
             dMatrix[ii+j*(N+1)]*
               scalar[SC_2D_INDEX(i,ii,iVar,iEl,N,nVar)]*   	       
               dsdx[TE_2D_INDEX(1,2,i,ii,0,iEl,N,1)];

    gf[1] += dMatrix[ii+i*(N+1)]*
	       scalar[SC_2D_INDEX(ii,j,iVar,iEl,N,nVar)]*
	       dsdx[TE_2D_INDEX(2,1,ii,j,0,iEl,N,1)]+
             dMatrix[ii+j*(N+1)]*
               scalar[SC_2D_INDEX(i,ii,iVar,iEl,N,nVar)]*   	       
               dsdx[TE_2D_INDEX(2,2,i,ii,0,iEl,N,1)];

  }
  gradF[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] = gf[0]/jacobian[SC_2D_INDEX(i,j,0,iEl,N,1)];
  gradF[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] = gf[1]/jacobian[SC_2D_INDEX(i,j,0,iEl,N,1)];

}

// GradientBR_MappedScalar2D_gpu
__global__ void DGGradientBoundary_MappedScalar2D_gpu(real *scalar, real *boundary, real *dsdx, real *jacobian, real* nHat, real *nScale, real *gradF, real *bMatrix, real *qWeights, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

  real gfx = 0.0;
  real gfy = 0.0;
  real f1, f2;
/*
  for( int ii=0; ii < N+1; ii++ ){
	  
    f1 = scalar[SC_2D_INDEX(ii,j,iVar,iEl,N,nVar)]*
	 dsdx[TE_2D_INDEX(1,1,ii,j,0,iEl,N,1)];
    
    f2 = scalar[SC_2D_INDEX(i,ii,iVar,iEl,N,nVar)]*
	 dsdx[TE_2D_INDEX(1,2,i,ii,0,iEl,N,1)];

    gfx += dgMatrix[ii+i*(N+1)]*f1+
           dgMatrix[ii+j*(N+1)]*f2;

    f1 = scalar[SC_2D_INDEX(ii,j,iVar,iEl,N,nVar)]*
	 dsdx[TE_2D_INDEX(2,1,ii,j,0,iEl,N,1)];
    
    f2 = scalar[SC_2D_INDEX(i,ii,iVar,iEl,N,nVar)]*
	 dsdx[TE_2D_INDEX(2,2,i,ii,0,iEl,N,1)];

    gfy += dgMatrix[ii+i*(N+1)]*f1+
           dgMatrix[ii+j*(N+1)]*f2;

  }
*/
  // East
  f1 = boundary[SCB_2D_INDEX(j,iVar,2,iEl,N,nVar)]*
	  nHat[VEB_2D_INDEX(1,j,0,2,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(j,0,2,iEl,N,1)];
  // West
  f2 = boundary[SCB_2D_INDEX(j,iVar,4,iEl,N,nVar)]*
	  nHat[VEB_2D_INDEX(1,j,0,4,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(j,0,4,iEl,N,1)];

  gfx = (f1*bMatrix[i+N+1] + f2*bMatrix[i])/qWeights[i];

  // North
  f1 = boundary[SCB_2D_INDEX(i,iVar,3,iEl,N,nVar)]*
	  nHat[VEB_2D_INDEX(1,i,0,3,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(i,0,3,iEl,N,1)];
  // South
  f2 = boundary[SCB_2D_INDEX(i,iVar,1,iEl,N,nVar)]*
	  nHat[VEB_2D_INDEX(1,i,0,1,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(i,0,1,iEl,N,1)];

  gfx += (f1*bMatrix[j+N+1] + f2*bMatrix[j])/qWeights[j];

  // East
  f1 = boundary[SCB_2D_INDEX(j,iVar,2,iEl,N,nVar)]*
	  nHat[VEB_2D_INDEX(2,j,0,2,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(j,0,2,iEl,N,1)];
  // West
  f2 = boundary[SCB_2D_INDEX(j,iVar,4,iEl,N,nVar)]*
	  nHat[VEB_2D_INDEX(2,j,0,4,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(j,0,4,iEl,N,1)];

  gfy = (f1*bMatrix[i+N+1] + f2*bMatrix[i])/qWeights[i];

  // North
  f1 = boundary[SCB_2D_INDEX(i,iVar,3,iEl,N,nVar)]*
	  nHat[VEB_2D_INDEX(2,i,0,3,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(i,0,3,iEl,N,1)];
  // South
  f2 = boundary[SCB_2D_INDEX(i,iVar,1,iEl,N,nVar)]*
	  nHat[VEB_2D_INDEX(2,i,0,1,iEl,N,1)]*
	  nScale[SCB_2D_INDEX(i,0,1,iEl,N,1)];

  gfy += (f1*bMatrix[j+N+1] + f2*bMatrix[j])/qWeights[j];

  gradF[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] += gfx/jacobian[SC_2D_INDEX(i,j,0,iEl,N,1)];
  gradF[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] += gfy/jacobian[SC_2D_INDEX(i,j,0,iEl,N,1)];

}

extern "C"
{
  void GradientBR_MappedScalar2D_gpu_wrapper(real **scalar, real **avgBoundary, real **dsdx, real **jacobian, real **nHat, real **nScale, real **gradF, real **dgMatrix, real **bMatrix, real **qWeights, int N, int nVar, int nEl)
  {

    GradientInterior_MappedScalar2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1)>>>(*scalar, *dsdx, *jacobian, *gradF,*dgMatrix, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());

    DGGradientBoundary_MappedScalar2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*scalar, *avgBoundary, *dsdx, *jacobian, *nHat, *nScale, *gradF, *bMatrix, *qWeights, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }
}


extern "C"
{
  void GradientSF_MappedScalar2D_gpu_wrapper(real **scalar, real **dsdx, real **jacobian, real **gradF, real **dMatrix, int N, int nVar, int nEl)
  {
    GradientInterior_MappedScalar2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*scalar, *dsdx, *jacobian, *gradF, *dMatrix, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// GradientInterior_MappedScalar3D_gpu
// Performs the interior operations for a scalar gradient calculation
// Use a derivative matrix for strong form;
// use a DG derivative matrix for weak form and follow up with the boundary contributions
__global__ void GradientInterior_MappedScalar3D_gpu(real *scalar, real *dsdx, real *jacobian, real *gradF, real *dMatrix, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  real gf[3] = {0.0};

  for( int ii=0; ii < N+1; ii++ ){
    gf[0] += dMatrix[ii+i*(N+1)]*
	       scalar[SC_3D_INDEX(ii,j,k,iVar,iEl,N,nVar)]*
	       dsdx[TE_3D_INDEX(1,1,ii,j,k,0,iEl,N,1)]+
             dMatrix[ii+j*(N+1)]*
               scalar[SC_3D_INDEX(i,ii,k,iVar,iEl,N,nVar)]*   	       
               dsdx[TE_3D_INDEX(1,2,i,ii,k,0,iEl,N,1)]+
             dMatrix[ii+k*(N+1)]*
               scalar[SC_3D_INDEX(i,j,ii,iVar,iEl,N,nVar)]*   	       
               dsdx[TE_3D_INDEX(1,3,i,j,ii,0,iEl,N,1)];

    gf[1] += dMatrix[ii+i*(N+1)]*
	       scalar[SC_3D_INDEX(ii,j,k,iVar,iEl,N,nVar)]*
	       dsdx[TE_3D_INDEX(2,1,ii,j,k,0,iEl,N,1)]+
             dMatrix[ii+j*(N+1)]*
               scalar[SC_3D_INDEX(i,ii,k,iVar,iEl,N,nVar)]*   	       
               dsdx[TE_3D_INDEX(2,2,i,ii,k,0,iEl,N,1)]+
             dMatrix[ii+k*(N+1)]*
               scalar[SC_3D_INDEX(i,j,ii,iVar,iEl,N,nVar)]*   	       
               dsdx[TE_3D_INDEX(2,3,i,j,ii,0,iEl,N,1)];

    gf[2] += dMatrix[ii+i*(N+1)]*
	       scalar[SC_3D_INDEX(ii,j,k,iVar,iEl,N,nVar)]*
	       dsdx[TE_3D_INDEX(3,1,ii,j,k,0,iEl,N,1)]+
             dMatrix[ii+j*(N+1)]*
               scalar[SC_3D_INDEX(i,ii,k,iVar,iEl,N,nVar)]*   	       
               dsdx[TE_3D_INDEX(3,2,i,ii,k,0,iEl,N,1)]+
             dMatrix[ii+k*(N+1)]*
               scalar[SC_3D_INDEX(i,j,ii,iVar,iEl,N,nVar)]*   	       
               dsdx[TE_3D_INDEX(3,3,i,j,ii,0,iEl,N,1)];

  }
  gradF[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] = gf[0]/jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];
  gradF[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] = gf[1]/jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];
  gradF[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] = gf[2]/jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];

}

// GradientBR_MappedScalar3D_gpu
__global__ void DGGradientBoundary_MappedScalar3D_gpu(real *scalar, real *boundary, real *dsdx, real *jacobian, real* nHat, real *nScale, real *gradF, real *bMatrix, real *qWeights, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.y;

  real gfx = 0.0;
  real gfy = 0.0;
  real gfz = 0.0;
  real f1, f2;

  // East
  f1 = boundary[SCB_3D_INDEX(j,k,iVar,2,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(1,j,k,0,2,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(j,k,0,2,iEl,N,1)];
  // West
  f2 = boundary[SCB_3D_INDEX(j,k,iVar,4,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(1,j,k,0,4,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(j,k,0,4,iEl,N,1)];

  gfx = (f1*bMatrix[i+N+1] + f2*bMatrix[i])/qWeights[i];

  // North
  f1 = boundary[SCB_3D_INDEX(i,k,iVar,3,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(1,i,k,0,3,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(i,k,0,3,iEl,N,1)];
  // South
  f2 = boundary[SCB_3D_INDEX(i,k,iVar,1,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(1,i,k,0,1,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(i,k,0,1,iEl,N,1)];

  gfx += (f1*bMatrix[j+N+1] + f2*bMatrix[j])/qWeights[j];

  // Top
  f1 = boundary[SCB_3D_INDEX(i,j,iVar,6,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(1,i,j,0,6,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(i,j,0,6,iEl,N,1)];
  // Bottom
  f2 = boundary[SCB_3D_INDEX(i,j,iVar,5,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(1,i,j,0,5,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(i,j,0,5,iEl,N,1)];

  gfx += (f1*bMatrix[k+N+1] + f2*bMatrix[k])/qWeights[k];

  // East
  f1 = boundary[SCB_3D_INDEX(j,k,iVar,2,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(2,j,k,0,2,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(j,k,0,2,iEl,N,1)];
  // West
  f2 = boundary[SCB_3D_INDEX(j,k,iVar,4,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(2,j,k,0,4,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(j,k,0,4,iEl,N,1)];

  gfy = (f1*bMatrix[i+N+1] + f2*bMatrix[i])/qWeights[i];

  // North
  f1 = boundary[SCB_3D_INDEX(i,k,iVar,3,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(2,i,k,0,3,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(i,k,0,3,iEl,N,1)];
  // South
  f2 = boundary[SCB_3D_INDEX(i,k,iVar,1,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(2,i,k,0,1,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(i,k,0,1,iEl,N,1)];

  gfy += (f1*bMatrix[j+N+1] + f2*bMatrix[j])/qWeights[j];

  // Top
  f1 = boundary[SCB_3D_INDEX(i,j,iVar,6,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(2,i,j,0,6,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(i,j,0,6,iEl,N,1)];
  // Bottom
  f2 = boundary[SCB_3D_INDEX(i,j,iVar,5,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(2,i,j,0,5,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(i,j,0,5,iEl,N,1)];

  gfy += (f1*bMatrix[k+N+1] + f2*bMatrix[k])/qWeights[k];

  // East
  f1 = boundary[SCB_3D_INDEX(j,k,iVar,2,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(3,j,k,0,2,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(j,k,0,2,iEl,N,1)];
  // West
  f2 = boundary[SCB_3D_INDEX(j,k,iVar,4,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(3,j,k,0,4,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(j,k,0,4,iEl,N,1)];

  gfz = (f1*bMatrix[i+N+1] + f2*bMatrix[i])/qWeights[i];

  // North
  f1 = boundary[SCB_3D_INDEX(i,k,iVar,3,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(3,i,k,0,3,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(i,k,0,3,iEl,N,1)];
  // South
  f2 = boundary[SCB_3D_INDEX(i,k,iVar,1,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(3,i,k,0,1,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(i,k,0,1,iEl,N,1)];

  gfz += (f1*bMatrix[j+N+1] + f2*bMatrix[j])/qWeights[j];

  // Top
  f1 = boundary[SCB_3D_INDEX(i,j,iVar,6,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(3,i,j,0,6,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(i,j,0,6,iEl,N,1)];
  // Bottom
  f2 = boundary[SCB_3D_INDEX(i,j,iVar,5,iEl,N,nVar)]*
	  nHat[VEB_3D_INDEX(3,i,j,0,5,iEl,N,1)]*
	  nScale[SCB_3D_INDEX(i,j,0,5,iEl,N,1)];

  gfz += (f1*bMatrix[k+N+1] + f2*bMatrix[k])/qWeights[k];

  gradF[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] += gfx/jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];
  gradF[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] += gfy/jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];
  gradF[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] += gfz/jacobian[SC_3D_INDEX(i,j,k,0,iEl,N,1)];

}

extern "C"
{
  void GradientBR_MappedScalar3D_gpu_wrapper(real **scalar, real **avgBoundary, real **dsdx, real **jacobian, real **nHat, real **nScale, real **gradF, real **dgMatrix, real **bMatrix, real **qWeights, int N, int nVar, int nEl)
  {

    GradientInterior_MappedScalar3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1)>>>(*scalar, *dsdx, *jacobian, *gradF,*dgMatrix, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());

    DGGradientBoundary_MappedScalar3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*scalar, *avgBoundary, *dsdx, *jacobian, *nHat, *nScale, *gradF, *bMatrix, *qWeights, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }
}


extern "C"
{
  void GradientSF_MappedScalar3D_gpu_wrapper(real **scalar, real **dsdx, real **jacobian, real **gradF, real **dMatrix, int N, int nVar, int nEl)
  {
    GradientInterior_MappedScalar3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*scalar, *dsdx, *jacobian, *gradF, *dMatrix, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

// ContravariantProjection_MappedVector2D_gpu
__global__ void ContravariantProjection_MappedVector2D_gpu(real *vector, real *dsdx, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;


    real Fx = vector[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)];
    real Fy = vector[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)];

    vector[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(1,1,i,j,0,iEl,N,1)]*Fx+
                			         dsdx[TE_2D_INDEX(2,1,i,j,0,iEl,N,1)]*Fy;

    vector[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] = dsdx[TE_2D_INDEX(1,2,i,j,0,iEl,N,1)]*Fx+
						 dsdx[TE_2D_INDEX(2,2,i,j,0,iEl,N,1)]*Fy;
  
}

extern "C"
{
  void ContravariantProjection_MappedVector2D_gpu_wrapper(real **vector, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantProjection_MappedVector2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*vector, *dsdx, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  } 
}

// ContravariantProjection_MappedP2Vector2D_gpu
__global__ void ContravariantProjection_MappedP2Vector2D_gpu(real *vector, real *physical, real *dsdx, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;


  for( int n = 0; n < N+1; n++){

    real Fx = physical[P2PVE_2D_INDEX(1,1,n,i,j,iVar,iEl,N,nVar)];
    real Fy = physical[P2PVE_2D_INDEX(2,1,n,i,j,iVar,iEl,N,nVar)];

    vector[P2VE_2D_INDEX(1,n,i,j,iVar,iEl,N,nVar)] = 
	    0.5*( dsdx[TE_2D_INDEX(1,1,i,j,0,iEl,N,1)] +
	          dsdx[TE_2D_INDEX(1,1,n,j,0,iEl,N,1)] )*Fx+
            0.5*( dsdx[TE_2D_INDEX(2,1,i,j,0,iEl,N,1)] +
                  dsdx[TE_2D_INDEX(2,1,n,j,0,iEl,N,1)] )*Fy;

    Fx = physical[P2PVE_2D_INDEX(1,2,n,i,j,iVar,iEl,N,nVar)];
    Fy = physical[P2PVE_2D_INDEX(2,2,n,i,j,iVar,iEl,N,nVar)];

    vector[P2VE_2D_INDEX(2,n,i,j,iVar,iEl,N,nVar)] = 
	    0.5*( dsdx[TE_2D_INDEX(1,2,i,j,0,iEl,N,1)] +
	          dsdx[TE_2D_INDEX(1,2,i,n,0,iEl,N,1)] )*Fx+
	    0.5*( dsdx[TE_2D_INDEX(2,2,i,j,0,iEl,N,1)] +
	          dsdx[TE_2D_INDEX(2,2,i,n,0,iEl,N,1)] )*Fy;
  }
  
}

extern "C"
{
  void ContravariantProjection_MappedP2Vector2D_gpu_wrapper(real **vector, real **physical, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantProjection_MappedP2Vector2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*vector, *physical, *dsdx, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
  }
}

// ContravariantProjection_MappedVector3D_gpu
__global__ void ContravariantProjection_MappedVector3D_gpu(real *vector, real *dsdx, int N, int nVar){

  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;
  size_t k = threadIdx.z;

  real Fx = vector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)];
  real Fy = vector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)];
  real Fz = vector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)];

  vector[VE_3D_INDEX(1,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,1,i,j,k,0,iEl,N,1)]*Fx+
                                                     dsdx[TE_3D_INDEX(2,1,i,j,k,0,iEl,N,1)]*Fy+
                                                     dsdx[TE_3D_INDEX(3,1,i,j,k,0,iEl,N,1)]*Fz;

  vector[VE_3D_INDEX(2,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,2,i,j,k,0,iEl,N,1)]*Fx+
                                                     dsdx[TE_3D_INDEX(2,2,i,j,k,0,iEl,N,1)]*Fy+
                                                     dsdx[TE_3D_INDEX(3,2,i,j,k,0,iEl,N,1)]*Fz;

  vector[VE_3D_INDEX(3,i,j,k,iVar,iEl,N,nVar)] = dsdx[TE_3D_INDEX(1,3,i,j,k,0,iEl,N,1)]*Fx+
                                                     dsdx[TE_3D_INDEX(2,3,i,j,k,0,iEl,N,1)]*Fy+
                                                     dsdx[TE_3D_INDEX(3,3,i,j,k,0,iEl,N,1)]*Fz;

}

extern "C"
{
  void ContravariantProjection_MappedVector3D_gpu_wrapper(real **vector, real **dsdx, int N, int nVar, int nEl)
  {
    ContravariantProjection_MappedVector3D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,N+1), 0, 0>>>(*vector, *dsdx, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
  }
}

__global__ void SideExchange_MappedScalar2D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nVar){

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
  void SideExchange_MappedScalar2D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    dim3 nblocks(4,nEl,1);
    dim3 nthreads(N+1,nVar,1);
    SideExchange_MappedScalar2D_gpu<<<nblocks,nthreads>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, offset, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void SideExchange_MappedVector2D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nVar){

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
  void SideExchange_MappedVector2D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    SideExchange_MappedVector2D_gpu<<<dim3(4,nEl,1), dim3(N+1,nVar,1), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, offset, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void SideExchange_MappedScalar3D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nVar){

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
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = boundary[SCB_3D_INDEX(i1,j1,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 1){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = boundary[SCB_3D_INDEX(i2,j2,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 2){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = boundary[SCB_3D_INDEX(i2,j2,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 3){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = boundary[SCB_3D_INDEX(i2,j2,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = boundary[SCB_3D_INDEX(i2,j2,ivar,s2,e2,N,nVar)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedScalar3D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    SideExchange_MappedScalar3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nVar), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, offset, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void SideExchange_MappedVector3D_gpu(real *extBoundary, real *boundary, int *sideInfo, int *elemToRank, int rankId, int offset, int N, int nVar){

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
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(1,i1,j1,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(2,i1,j1,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(3,i1,j1,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 1){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(1,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(2,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(3,i2,j2,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 2){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(1,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(2,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(3,i2,j2,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 3){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(1,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(2,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(3,i2,j2,ivar,s2,e2,N,nVar)];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = i1;
        extBoundary[VEB_3D_INDEX(1,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(1,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(2,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(2,i2,j2,ivar,s2,e2,N,nVar)];
        extBoundary[VEB_3D_INDEX(3,i1,j1,ivar,s1,e1,N,nVar)] = boundary[VEB_3D_INDEX(3,i2,j2,ivar,s2,e2,N,nVar)];
      }
    }
  }
  
}

extern "C"
{
  void SideExchange_MappedVector3D_gpu_wrapper(real **extBoundary, real **boundary, int **sideInfo, int **elemToRank, int rankId, int offset, int N, int nVar, int nEl)
  {
    SideExchange_MappedVector3D_gpu<<<dim3(6,nEl,1), dim3(N+1,N+1,nVar), 0, 0>>>(*extBoundary, *boundary, *sideInfo, *elemToRank, rankId, offset, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
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
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void ApplyFlip_MappedScalar2D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t ivar = blockIdx.x;
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
        extBuff[i1] = extBoundary[SCB_2D_INDEX(i2,ivar,s1,e1,N,nVar)];
        __syncthreads();

        extBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nVar)] = extBuff[i1];
      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedScalar2D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedScalar2D_gpu<<<dim3(nVar,4,nEl), dim3(N+1,1,1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void ApplyFlip_MappedVector2D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t ivar = blockIdx.x;
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
        extBuff[i1] = extBoundary[VEB_2D_INDEX(1,i2,ivar,s1,e1,N,nVar)];
        __syncthreads();

        extBoundary[VEB_2D_INDEX(1,i1,ivar,s1,e1,N,nVar)] = extBuff[i1];

        extBuff[i1] = extBoundary[VEB_2D_INDEX(2,i2,ivar,s1,e1,N,nVar)];
        __syncthreads();

        extBoundary[VEB_2D_INDEX(2,i1,ivar,s1,e1,N,nVar)] = extBuff[i2];

      }
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedVector2D_gpu_wrapper(real **extBoundary, int **sideInfo, int **elemToRank, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedVector2D_gpu<<<dim3(nVar,4,nEl), dim3(N+1,1,1), 0, 0>>>(*extBoundary, *sideInfo, *elemToRank, rankId, N, nVar);
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void ApplyFlip_MappedScalar3D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y+1;
  size_t e1 = blockIdx.z;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.x;
  
  int e2 = sideInfo[INDEX3(2,s1-1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1-1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1-1,e1,5,4)]-s2*10;
  int bcid = sideInfo[INDEX3(4,s1-1,e1,5,4)];

  __shared__ real extBuff[256];

  extBuff[i1+(N+1)*j1] = extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank /= rankId ){
      if(flip == 1){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 2){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 3){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = extBuff[i2+(N+1)*j2];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = i1;
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
    HIP_SAFE_CALL(hipGetLastError());
  }

}

__global__ void ApplyFlip_MappedVector3D_gpu(real *extBoundary, int *sideInfo, int *elemToRank, int rankId, int N, int nVar){

  size_t ivar = blockIdx.x;
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

  extBuff[dir+3*(i1+(N+1)*j1)] = extBoundary[VEB_3D_INDEX(dir+1,i1,j1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0){
    int neighborRank = elemToRank[e2];
    if( neighborRank /= rankId ){
      if(flip == 1){
        int i2 = j1;
        int j2 = N-i1;
        extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[dir+3*(i2+(N+1)*j2)];
      }
      else if(flip == 2){
        int i2 = N-i1;
        int j2 = N-j1;
        extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[dir+3*(i2+(N+1)*j2)];
      }
      else if(flip == 3){
        int i2 = N-j1;
        int j2 = i1;
        extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[dir+3*(i2+(N+1)*j2)];
      }
      else if(flip == 4){
        int i2 = j1;
        int j2 = i1;
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
    HIP_SAFE_CALL(hipGetLastError());
  }

}
