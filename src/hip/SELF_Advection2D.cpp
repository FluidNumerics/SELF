#include <hip/hip_runtime.h>
#include "SELF_HIP_Macros.h"

__global__ void InitializeGRK3_Advection2D_gpu(real *gRK3, int N, int nVar){

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    gRK3[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] =  0.0; 

}

extern "C"
{
  void InitializeGRK3_Advection2D_gpu_wrapper(real **gRK3, int N, int nVar, int nEl)
  {

    // Block size is set to match the size of the element exactly
    // Grid size is set to ( number of tracers X number of elements )
    // DGSEM is beautiful
    InitializeGRK3_Advection2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*gRK3, N, nVar);
  }
}

__global__ void UpdateGRK3_Advection2D_gpu(real *gRK3, real *solution, real *dSdt, int rk3_a, int rk3_g, real dt, int N, int nVar){

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    gRK3[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] =  rk3_a*gRK3[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]+
	    dSdt[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

    solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)] =  solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)]+
	    rk3_g*dt*gRK3[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

}

extern "C"
{
  void UpdateGRK3_Advection2D_gpu_wrapper(real **gRK3, real **solution, real **dSdt, int rk3_a, int rk3_g, real dt, int N, int nVar, int nEl)
  {

    // Block size is set to match the size of the element exactly
    // Grid size is set to ( number of tracers X number of elements )
    // DGSEM is beautiful
    UpdateGRK3_Advection2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*gRK3, *solution, *dSdt, rk3_a, rk3_g, dt, N, nVar);
  }
}

__global__ void InternalFlux_Advection2D_gpu(real *flux, real *solution, real *velocity, int N, int nVar){

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iEl = blockIdx.y;
  size_t i = threadIdx.x;
  size_t j = threadIdx.y;

    // x-component of the advective flux
    flux[VE_2D_INDEX(1,i,j,iVar,iEl,N,nVar)] = 
	    velocity[VE_2D_INDEX(1,i,j,1,iEl,N,nVar)]*
            solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

    // y-component of the advective flux
    flux[VE_2D_INDEX(2,i,j,iVar,iEl,N,nVar)] = 
	    velocity[VE_2D_INDEX(2,i,j,1,iEl,N,nVar)]*
            solution[SC_2D_INDEX(i,j,iVar,iEl,N,nVar)];

}

extern "C"
{
  void InternalFlux_Advection2D_gpu_wrapper(real **flux, real **solution, real **velocity, int N, int nVar, int nEl)
  {

    // Block size is set to match the size of the element exactly
    // Grid size is set to ( number of tracers X number of elements )
    // DGSEM is beautiful
    InternalFlux_Advection2D_gpu<<<dim3(nVar,nEl,1), dim3(N+1,N+1,1), 0, 0>>>(*flux, *solution, *velocity, N, nVar);
  }
}

__global__ void SideFlux_Advection2D_gpu(real *flux, real *boundarySol, real *extSol, real *velocity, real *nhat, real* nscale, int N, int nVar){

  // Get the array indices from the GPU thread IDs
  size_t iVar = blockIdx.x;
  size_t iSide = blockIdx.y;
  size_t iEl = blockIdx.z;
  size_t i = threadIdx.x;

  // Calculate the normal velocity at the cell sides/edges
  float un = velocity[VEB_2D_INDEX(1,i,1,iSide,iEl,N,nVar)]*nhat[0] +
             velocity[VEB_2D_INDEX(2,i,1,iSide,iEl,N,nVar)]*nhat[1];

  // Get the external and internal states for the Riemann solver
  float extState = extSol[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)];
  float intState = boundarySol[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)];
  float nmag = nscale[SCB_2D_INDEX(i,iVar,iSide,iEl,N,nVar)];

  // Calculate the normal flux at the cell sides/edges
  // using a Lax-Friedrich's upwind Riemann solver
  // Right.... the flux is a vector.. gotta hack this..
  // The normal flux at the cell boundaries will be stored in the
  // first vector component..I won't forget this..
  flux[VEB_2D_INDEX(1,i,iVar,iSide,iEl,N,nVar)] =
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
