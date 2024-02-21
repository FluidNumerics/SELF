// ******************************************************************** //
// 
// Copyright 2020 Fluid Numerics LLC
// Author : Joseph Schoonover (joe@fluidnumerics.com)
// Support : self-fluids@fluidnumerics.com
//
// ******************************************************************** //
#ifndef SELF_Macros_H
#define SELF_Macros_H
#endif

#include <float.h>

#ifdef DOUBLE_PRECISION
  typedef double real;
#else
  typedef float real;
#endif

#define HIP_SAFE_CALL( call) {                                    \
    hipError_t err = call;                                                    \
    if( hipSuccess != err) {                                                \
        fprintf(stderr, "HIP error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, hipGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }


#define INDEX(i,j,N) i+j*(N+1)
#define INDEX3(i,j,k,Ni,Nj) i+Ni*(j+Nj*k)

#define SC_1D_INDEX(i,iel,iVar,N,nEl) i+(N+1)*(iEl + nEl*iVar)
#define SCB_1D_INDEX(iSide,iel,iVar,N,nEl) iSide + 2*(iel + nEl*iVar)

#define SC_2D_INDEX(i,j,iel,iVar,N,nEl) i+(N+1)*(j + (N+1)*(iel + nEl*iVar)) 
#define SCB_2D_INDEX(i,iSide,iel,iVar,N,nEl) i+(N+1)*(iSide + 4*(iel + nEl*iVar))

#define VE_2D_INDEX(i,j,iel,iVar,idir,N,nEl,nVar) i + (N+1)*(j + (N+1)*(iel + nEl*(iVar + nVar*idir)))
#define VEB_2D_INDEX(i,iSide,iel,iVar,idir,N,nEl,nVar) i + (N+1)*(iSide + 4*(iel + nEl*(iVar + nVar*idir)))


#define P2VE_2D_INDEX(dir,n,i,j,iel,iVar,N,nVar) dir-1 + 2*(n + (N+1)*(i + (N+1)*(j + (N+1)*(iel + nEl*iVar))))
#define P2PVE_2D_INDEX(row,col,n,i,j,iel,iVar,N,nVar) row-1 + 2*(col-1 + 2*(n + (N+1)*(i + (N+1)*(j + (N+1)*(iel + nEl*iVar)))))
#define TE_2D_INDEX(row,col,i,j,iel,iVar,N,nVar) row-1 + 2*(col-1 + 2*(i + (N+1)*(j + (N+1)*(iel + nEl*iVar))))

#define TEB_2D_INDEX(row,col,i,iVar,iSide,iel,N,nVar) row-1 + 2*(col-1 + 2*(i + (N+1)*(iVar + nVar*(iSide-1 + 4*iel))))


#define SC_3D_INDEX(i,j,k,iel,iVar,N,nVar) i+(N+1)*(j + (N+1)*(k + (N+1)*(iel + nEl*iVar)))
#define VE_3D_INDEX(dir,i,j,k,iel,iVar,N,nVar) dir-1 + 3*(i + (N+1)*(j + (N+1)*(k + (N+1)*(iel + nEl*iVar))))
#define TE_3D_INDEX(row,col,i,j,k,iel,iVar,N,nVar) row-1 + 3*(col-1 + 3*(i + (N+1)*(j + (N+1)*(k + (N+1)*(iel + nEl*iVar)))))
#define SCB_3D_INDEX(i,j,iVar,iSide,iel,N,nVar) i+(N+1)*(j + (N+1)*(iVar + nVar*(iSide-1 + 6*iel)))
#define VEB_3D_INDEX(dir,i,j,iVar,iSide,iel,N,nVar) dir-1 + 3*(i + (N+1)*(j + (N+1)*(iVar + nVar*(iSide-1 + 6*iel))))
#define TEB_3D_INDEX(row,col,i,j,iVar,iSide,iel,N,nVar) row-1 + 3*(col-1 + 3*(i + (N+1)*(j + (N+1)*(iVar + nVar*(iSide-1 + 6*iel)))))

// Boundary condition flags //
//
//  Conditions on the solution
#define SELF_BC_PRESCRIBED 100
#define SELF_BC_RADIATION 101
#define SELF_BC_NONORMALFLOW 102

// Conditions on the solution gradients
#define SELF_BC_PRESCRIBED_STRESS 200
#define SELF_BC_NOSTRESS 201
