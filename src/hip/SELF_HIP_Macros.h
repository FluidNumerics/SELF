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

#define INDEX(i,j,N) i+j*(N+1)
#define INDEX3(i,j,k,Ni,Nj) i+Ni*(j+Nj*k)
#define SC_1D_INDEX(i,iVar,iel,N,nVar) i+(N+1)*(iVar + nVar*iel)
#define SCB_1D_INDEX(iVar,iSide,iel,N,nVar) iVar + nVar*(iSide + 2*iel)

#define SC_2D_INDEX(i,j,iVar,iel,N,nVar) i+(N+1)*(j + (N+1)*(iVar + nVar*iel)) 
#define VE_2D_INDEX(dir,i,j,iVar,iel,N,nVar) dir-1 + 2*(i + (N+1)*(j + (N+1)*(iVar + nVar*iel)))
#define P2VE_2D_INDEX(dir,n,i,j,iVar,iel,N,nVar) dir-1 + 2*(n + (N+1)*(i + (N+1)*(j + (N+1)*(iVar + nVar*iel))))
#define P2PVE_2D_INDEX(row,col,n,i,j,iVar,iel,N,nVar) row-1 + 2*(col-1 + 2*(n + (N+1)*(i + (N+1)*(j + (N+1)*(iVar + nVar*iel)))))
#define TE_2D_INDEX(row,col,i,j,iVar,iel,N,nVar) row-1 + 2*(col-1 + 2*(i + (N+1)*(j + (N+1)*(iVar + nVar*iel))))
#define SCB_2D_INDEX(i,iVar,iSide,iel,N,nVar) i+(N+1)*(iVar + nVar*(iSide-1 + 4*iel))
#define VEB_2D_INDEX(dir,i,iVar,iSide,iel,N,nVar) dir-1 + 2*(i + (N+1)*(iVar + nVar*(iSide-1 + 4*iel)))
#define TEB_2D_INDEX(row,col,i,iVar,iSide,iel,N,nVar) row-1 + 2*(col-1 + 2*(i + (N+1)*(iVar + nVar*(iSide-1 + 4*iel))))


#define SC_3D_INDEX(i,j,k,iVar,iel,N,nVar) i+(N+1)*(j + (N+1)*(k + (N+1)*(iVar + nVar*iel)))
#define VE_3D_INDEX(dir,i,j,k,iVar,iel,N,nVar) dir-1 + 3*(i + (N+1)*(j + (N+1)*(k + (N+1)*(iVar + nVar*iel))))
#define TE_3D_INDEX(row,col,i,j,k,iVar,iel,N,nVar) row-1 + 3*(col-1 + 3*(i + (N+1)*(j + (N+1)*(k + (N+1)*(iVar + nVar*iel)))))
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
