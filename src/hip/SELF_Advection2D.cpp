#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"
#include <cstdio>

__global__ void UpdateGRK3_Advection2D_gpu(real *gRK3, real *solution, real *dSdt, real rk3A, real rk3G, real dt, int N, int nVar){

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    gRK3[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] =  rk3A*gRK3[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]+
	    dSdt[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

    solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] += rk3G*dt*gRK3[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

}

extern "C"
{
  void UpdateGRK3_Advection2D_gpu_wrapper(real **gRK3, real **solution, real **dSdt, real rk3A, real rk3G, real dt, int N, int nVar, int nEl)
  {
    UpdateGRK3_Advection2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*gRK3, *solution, *dSdt, rk3A, rk3G, dt, N, nVar);
  }
}

__global__ void InternalFlux_Advection2D_gpu(real *flux, real *solution, real *velocity, real *dsdx, int N, int nVar){

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    real Fx = velocity[VE_2D_INDEX(1,i,j,0,iEl,N,1)]*
              solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

    real Fy = velocity[VE_2D_INDEX(2,i,j,0,iEl,N,1)]*
              solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

    // x-component of the advective flux
    flux[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] = 
                dsdx[TE_2D_INDEX(1,1,i,j,0,iEl,N,1)]*Fx+
                dsdx[TE_2D_INDEX(2,1,i,j,0,iEl,N,1)]*Fy; 

    // y-component of the advective flux
    flux[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] = 
                dsdx[TE_2D_INDEX(1,2,i,j,0,iEl,N,1)]*Fx+
                dsdx[TE_2D_INDEX(2,2,i,j,0,iEl,N,1)]*Fy;

}

extern "C"
{
  void InternalFlux_Advection2D_gpu_wrapper(real **flux, real **solution, real **velocity, real **dsdx, int N, int nVar, int nEl)
  {

    // Block size is set to match the size of the element exactly
    // Grid size is set to ( number of tracers X number of elements )
    // DGSEM is beautiful
    InternalFlux_Advection2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*flux, *solution, *velocity, *dsdx, N, nVar);
  }
}

__global__ void InternalDiffusiveFlux_Advection2D_gpu(real *flux, real *solutionGradient, real *dsdx, real diffusivity, int N, int nVar){

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    real Fx = diffusivity*
              solutionGradient[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)];

    real Fy = diffusivity*
              solutionGradient[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)];

    // xi-component of the advective flux
    flux[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] += 
                -(dsdx[TE_2D_INDEX(1,1,i,j,0,iEl,N,1)]*Fx+
                  dsdx[TE_2D_INDEX(2,1,i,j,0,iEl,N,1)]*Fy); 

    // eta-component of the advective flux
    flux[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] += 
                -(dsdx[TE_2D_INDEX(1,2,i,j,0,iEl,N,1)]*Fx+
                  dsdx[TE_2D_INDEX(2,2,i,j,0,iEl,N,1)]*Fy);

}

extern "C"
{
  void InternalDiffusiveFlux_Advection2D_gpu_wrapper(real **flux, real **solutionGradient, real **dsdx, real diffusivity, int N, int nVar, int nEl)
  {

    // Block size is set to match the size of the element exactly
    // Grid size is set to ( number of tracers X number of elements )
    // DGSEM is beautiful
    InternalDiffusiveFlux_Advection2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*flux, *solutionGradient,*dsdx, diffusivity, N, nVar);
  }
}

__global__ void SideFlux_Advection2D_gpu(real *flux, real *boundarySol, real *extSol, real *velocity, real *nhat, real* nscale, int N, int nVar){

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iSide = blockIdx.y+1;
  size_t iEl = blockIdx.z;
  size_t i = threadIdx.x;

  // Calculate the normal velocity at the cell sides/edges
  float un = velocity[VEB_2D_INDEX(1,i,0,iSide,iEl,N,1)]*
	     nhat[VEB_2D_INDEX(1,i,0,iSide,iEl,N,1)] +
             velocity[VEB_2D_INDEX(2,i,0,iSide,iEl,N,1)]*
	     nhat[VEB_2D_INDEX(2,i,0,iSide,iEl,N,1)];

  // Get the external and internal states for the Riemann solver
  float extState = extSol[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)];
  float intState = boundarySol[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)];
  float nmag = nscale[SCB_2D_INDEX(i,0,iSide,iEl,N,1)];

  // Calculate the normal flux at the cell sides/edges
  flux[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)] =
                    0.5*( un*(intState + extState) - fabs(un)*(extState - intState) )*nmag;


}

extern "C"
{
  void SideFlux_Advection2D_gpu_wrapper(real **flux, real **boundarySol, real **extSol, real **velocity, real **nhat, real **nscale, int N, int nVar, int nEl)
  {

    // Block size is set to match the size of the element exactly
    // Grid size is set to ( number of tracers X number of elements )
    // DGSEM is beautiful
    SideFlux_Advection2D_gpu<<<dim3(nVar,4,nEl), dim3(N+1,1,1), 0, 0>>>(*flux, *boundarySol, *extSol, *velocity, *nhat, *nscale, N, nVar);
  }
}

__global__ void SideDiffusiveFlux_Advection2D_gpu(real *flux, real *boundarySolGrad, real *extSolGrad, real *nhat, real* nscale, real diffusivity, int N, int nVar){

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iSide = blockIdx.y+1;
  size_t iEl = blockIdx.z;
  size_t i = threadIdx.x;


  float extState = extSolGrad[VEB_2D_INDEX(1,i,iVar,iSide,iEl,N,nVar)]*
	           nhat[VEB_2D_INDEX(1,i,0,iSide,iEl,N,1)] +
                   extSolGrad[VEB_2D_INDEX(2,i,iVar,iSide,iEl,N,nVar)]*
	           nhat[VEB_2D_INDEX(2,i,0,iSide,iEl,N,1)];

  float intState = boundarySolGrad[VEB_2D_INDEX(1,i,iVar,iSide,iEl,N,nVar)]*
	           nhat[VEB_2D_INDEX(1,i,0,iSide,iEl,N,1)] +
                   boundarySolGrad[VEB_2D_INDEX(2,i,iVar,iSide,iEl,N,nVar)]*
	           nhat[VEB_2D_INDEX(2,i,0,iSide,iEl,N,1)];

  float nmag = nscale[SCB_2D_INDEX(i,0,iSide,iEl,N,1)];

  // Calculate the normal flux at the cell sides/edges
  flux[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)] +=
                    -0.5*diffusivity*(intState + extState)*nmag;


}

extern "C"
{
  void SideDiffusiveFlux_Advection2D_gpu_wrapper(real **flux, real **boundarySolGrad, real **extSolGrad, real **nhat, real **nscale, real diffusivity, int N, int nVar, int nEl)
  {

    SideDiffusiveFlux_Advection2D_gpu<<<dim3(nVar,4,nEl), dim3(N+1,1,1), 0, 0>>>(*flux, *boundarySolGrad, *extSolGrad, *nhat, *nscale, diffusivity, N, nVar);

  }
}
