#ifndef SELF_GPU_H
#define SELF_GPU_H
#endif

#include <float.h>

#ifdef DOUBLE_PRECISION
  typedef double real;
#else
  typedef float real;
#endif


#include <cassert>
#include <climits>
#include <cstdio>

#ifdef HAVE_HIP

#include <hip/hip_runtime.h>

#ifdef __HIPCC__
#define __CUDACC__
#endif

__attribute__((unused))
static void check(const hipError_t err, const char *const file, const int line)
{
  if (err == hipSuccess) return;
  fprintf(stderr,"HIP ERROR AT LINE %d OF FILE '%s': %s %s\n",line,file,hipGetErrorName(err),hipGetErrorString(err));
  fflush(stderr);
  exit(err);
}

#elif HAVE_CUDA

#include <cuda_runtime.h>

#define hipLaunchKernelGGL(F,G,B,M,S,...) F<<<G,B,M,S>>>(__VA_ARGS__)

__attribute__((unused))
static void check(const cudaError_t err, const char *const file, const int line)
{
  if (err == cudaSuccess) return;
  fprintf(stderr,"CUDA ERROR AT LINE %d OF FILE '%s': %s %s\n",line,file,cudaGetErrorName(err),cudaGetErrorString(err));
  fflush(stderr);
  exit(err);
}

#endif

#define CHECK(X) check(X,__FILE__,__LINE__)


#define INDEX(i,j,N) i+j*(N+1)
#define INDEX3(i,j,k,Ni,Nj) i+Ni*(j+Nj*k)

#define SC_1D_INDEX(i,iel,iVar,N,nEl) i+(N+1)*(iEl + nEl*iVar)
#define SCB_1D_INDEX(j,iel,iVar,N,nEl) j + 2*(iel + nEl*iVar)

#define SC_2D_INDEX(i,j,iel,iVar,N,nEl) i+(N+1)*(j + (N+1)*(iel + nEl*iVar)) 
#define SCB_2D_INDEX(i,j,iel,iVar,N,nEl) i+(N+1)*(j + 4*(iel + nEl*iVar))

#define SC_3D_INDEX(i,j,k,iel,iVar,N,nEl) i+(N+1)*(j + (N+1)*(k + (N+1)*(iel + nEl*iVar)))
#define SCB_3D_INDEX(i,j,k,iel,iVar,N,nEl) i+(N+1)*(j + (N+1)*(k + 6*(iel + nEl*iVar)))

#define VE_2D_INDEX(i,j,iel,iVar,idir,N,nEl,nVar) i + (N+1)*(j + (N+1)*(iel + nEl*(iVar + nVar*idir)))
#define VEB_2D_INDEX(i,j,iel,iVar,idir,N,nEl,nVar) i + (N+1)*(j + 4*(iel + nEl*(iVar + nVar*idir)))

#define VE_3D_INDEX(i,j,k,iel,iVar,idir,N,nEl,nVar) i + (N+1)*(j + (N+1)*(k + (N+1)*(iel + nEl*(iVar + nVar*idir))))
#define VEB_3D_INDEX(i,j,k,iel,iVar,idir,N,nEl,nVar) i + (N+1)*(j + (N+1)*(k + 6*(iel + nEl*(iVar + nVar*idir))))

#define TE_2D_INDEX(i,j,iel,iVar,row,col,N,nEl,nVar) i + (N+1)*(j + (N+1)*(iel + nEl*(iVar + nVar*(row + 2*col))))
#define TEB_2D_INDEX(i,j,iel,iVar,row,col,N,nEl,nVar) i + (N+1)*(j + 4*(iel + nEl*(iVar + nVar*(row + 2*col))))

#define TE_3D_INDEX(i,j,k,iel,iVar,row,col,N,nEl,nVar) i + (N+1)*(j + (N+1)*(k + (N+1)*(iel + nEl*(iVar + nVar*(row + 3*col)))))
#define TEB_3D_INDEX(i,j,k,iel,iVar,row,col,N,nEl,nVar) i + (N+1)*(j + (N+1)*(k + 6*(iel + nEl*(iVar + nVar*(row + 3*col)))))


// Boundary condition flags //
//
//  Conditions on the solution
#define SELF_BC_PRESCRIBED 100
#define SELF_BC_RADIATION 101
#define SELF_BC_NONORMALFLOW 102

// Conditions on the solution gradients
#define SELF_BC_PRESCRIBED_STRESS 200
#define SELF_BC_NOSTRESS 201