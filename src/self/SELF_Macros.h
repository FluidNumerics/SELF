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
#define SC_1D_INDEX(i,iVar,iel,N,nVar) i+(N+1)*(iVar + nVar*iel)
#define SC_2D_INDEX(i,j,iVar,iel,N,nVar) i+(N+1)*(j + (N+1)*(iVar + nVar*iel)) 
#define VE_2D_INDEX(dir,i,j,iVar,iel,N,nVar) dir + 2*(i + (N+1)*(j + (N+1)*(iVar + nVar*iel)))
#define TE_2D_INDEX(row,col,i,j,iVar,iel,N,nVar) row + 2*(col + 2*(i + (N+1)*(j + (N+1)*(iVar + nVar*iel))))
#define SC_3D_INDEX(i,j,k,iVar,iel,N,nVar) i+(N+1)*(j + (N+1)*(k + (N+1)*(iVar + nVar*iel)))
#define VE_3D_INDEX(dir,i,j,k,iVar,iel,N,nVar) dir + 3*(i + (N+1)*(j + (N+1)*(k + (N+1)*(iVar + nVar*iel))))
#define TE_3D_INDEX(row,col,i,j,k,iVar,iel,N,nVar) row + 2*(col + 2*(i + (N+1)*(j + (N+1)*(k + (N+1)*(iVar + nVar*iel)))))

