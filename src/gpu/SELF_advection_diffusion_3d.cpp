#include "SELF_GPU_Macros.h"


__global__ void setboundarycondition_advection_diffusion_3d_gpukernel(real *extBoundary, real *boundary, int *sideInfo, int N, int nel, int nvar){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*(N+1)*6*nel;

  if(idof < ndof){
    uint32_t i = idof % (N+1);
    uint32_t j = (idof/(N+1)) % (N+1);
    uint32_t s1 = (idof/(N+1)/(N+1)) % 6;
    uint32_t e1 = idof/(N+1)/(N+1)/6;
    uint32_t e2 = sideInfo[INDEX3(2,s1,e1,5,6)];
    if( e2 == 0){
        uint32_t ivar = blockIdx.y;
        extBoundary[SCB_3D_INDEX(i,j,s1,e1,ivar,N,nel)] = 0.0;
    }
  }
}

extern "C" 
{
  void setboundarycondition_advection_diffusion_3d_gpu(real *extBoundary, real *boundary, int *sideInfo, int N, int nel, int nvar){
    int threads_per_block = 256;
    int ndof = (N+1)*(N+1)*6*nel;
    int nblocks_x = ndof/threads_per_block +1;

    dim3 nblocks(nblocks_x,nvar,1);
    dim3 nthreads(threads_per_block,1,1);

	setboundarycondition_advection_diffusion_3d_gpukernel<<<nblocks,nthreads, 0, 0>>>(extBoundary,boundary,sideInfo,N,nel,nvar);
  }
}

__global__ void setgradientboundarycondition_advection_diffusion_3d_gpukernel(real *extBoundary, real *boundary, int *sideInfo, int N, int nel, int nvar){

  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*(N+1)*6*nel;

  if(idof < ndof){
    uint32_t i = idof % (N+1);
    uint32_t j = (idof/(N+1)) % (N+1);
    uint32_t s1 = (idof/(N+1)/(N+1)) % 6;
    uint32_t e1 = idof/(N+1)/(N+1)/6;
    uint32_t e2 = sideInfo[INDEX3(2,s1,e1,5,6)];
    if( e2 == 0){
        uint32_t ivar = blockIdx.y;
        uint32_t idir = blockIdx.z;
        extBoundary[VEB_3D_INDEX(i,j,s1,e1,ivar,idir,N,nel,nvar)] = boundary[VEB_3D_INDEX(i,j,s1,e1,ivar,idir,N,nel,nvar)];
    }
  }
}

extern "C" 
{
  void setgradientboundarycondition_advection_diffusion_3d_gpu(real *extBoundary, real *boundary, int *sideInfo, int N, int nel, int nvar){
    int threads_per_block = 256;
    int ndof = (N+1)*(N+1)*6*nel;
    int nblocks_x = ndof/threads_per_block +1;

    dim3 nblocks(nblocks_x,nvar,3);
    dim3 nthreads(threads_per_block,1,1);

	setgradientboundarycondition_advection_diffusion_3d_gpukernel<<<nblocks,nthreads, 0, 0>>>(extBoundary,boundary,sideInfo,N,nel,nvar);
  }
}

__global__ void fluxmethod_advection_diffusion_3d_gpukernel(real *solution, real *solutiongradient, real *flux, real u, real v, real w, real nu, int ndof){
  uint32_t i = threadIdx.x + blockIdx.x*blockDim.x;

  if( i < ndof ){
    flux[i] = u*solution[i] - nu*solutiongradient[i];
    flux[i+ndof] = v*solution[i] - nu*solutiongradient[i+ndof];
    flux[i+2*ndof] = w*solution[i] - nu*solutiongradient[i+2*ndof];
  }

}
extern "C"
{
  void fluxmethod_advection_diffusion_3d_gpu(real *solution, real *solutiongradient, real *flux, real u, real v, real w, real nu, int N, int nel, int nvar){
    int ndof = (N+1)*(N+1)*(N+1)*nel*nvar;
    int threads_per_block = 256;
    int nblocks_x = ndof/threads_per_block +1;
    fluxmethod_advection_diffusion_3d_gpukernel<<<dim3(nblocks_x,1,1), dim3(threads_per_block,1,1), 0, 0>>>(solution,solutiongradient,flux,u,v,w,nu,ndof);
  }

}

__global__ void boundaryflux_advection_diffusion_3d_gpukernel(real *fb, real *fextb, real *dfavg, real *nhat, real *nscale, real *flux, real u, real v, real w, real nu, int N, int nel, int nvar){
  uint32_t idof = threadIdx.x + blockIdx.x*blockDim.x;
  uint32_t ndof = (N+1)*(N+1)*6*nel;

  if( idof < ndof ){

    uint32_t i = idof % (N+1);
    uint32_t j = (idof/(N+1)) % (N+1);
    uint32_t k = (idof/(N+1)/(N+1)) % 6;
    uint32_t iel = idof/(N+1)/(N+1)/6;
    uint32_t ivar = blockIdx.y;

    real nx = nhat[VEB_3D_INDEX(i,j,k,iel,0,0,N,nel,1)];
    real ny = nhat[VEB_3D_INDEX(i,j,k,iel,0,1,N,nel,1)];
    real nz = nhat[VEB_3D_INDEX(i,j,k,iel,0,2,N,nel,1)];

    real un = u*nx+v*ny+w*nz;

    real dfdn = dfavg[VEB_3D_INDEX(i,j,k,iel,ivar,0,N,nel,nvar)]*nx+
                dfavg[VEB_3D_INDEX(i,j,k,iel,ivar,1,N,nel,nvar)]*ny+
                dfavg[VEB_3D_INDEX(i,j,k,iel,ivar,2,N,nel,nvar)]*nz;

    real nmag = nscale[SCB_3D_INDEX(i,j,k,iel,0,N,nel)];

    flux[idof+ivar*ndof] = (0.5*(un*(fb[idof+ivar*ndof]+fextb[idof+ivar*ndof])+
                            fabsf(un)*(fb[idof+ivar*ndof]-fextb[idof+ivar*ndof]))-
                            nu*dfdn)*nmag;
  }

}
extern "C"
{
  void boundaryflux_advection_diffusion_3d_gpu(real *fb, real *fextb, real *dfavg, real *nhat, real *nscale, real *flux, real u, real v, real w, real nu, int N, int nel, int nvar){
    int threads_per_block = 256;
    uint32_t ndof = (N+1)*(N+1)*6*nel;
    int nblocks_x = ndof/threads_per_block +1;

    dim3 nblocks(nblocks_x,nvar,1);
    dim3 nthreads(threads_per_block,1,1);

    boundaryflux_advection_diffusion_3d_gpukernel<<<nblocks,nthreads>>>(fb,fextb,dfavg,nhat,nscale,flux,u,v,w,nu,N,nel,nvar);
  }

}

