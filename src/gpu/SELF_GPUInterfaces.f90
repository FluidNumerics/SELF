! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2024 Fluid Numerics LLC
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
!    the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_GPUInterfaces

  use iso_c_binding
  use SELF_GPU

  ! Data
  interface
    subroutine Average_gpu(favg,f1,f2,ndof) bind(c,name="Average_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: favg,f1,f2
      integer(c_int),value :: ndof
    endsubroutine Average_gpu
  endinterface

  interface
    subroutine BoundaryInterp_2D_gpu(bMatrix_dev,f_dev,bf_dev,N,nVar,nEl) &
      bind(c,name="BoundaryInterp_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: bMatrix_dev,f_dev,bf_dev
      integer(c_int),value :: N,nVar,nEl
    endsubroutine BoundaryInterp_2D_gpu
  endinterface

  interface
    subroutine BoundaryInterp_3D_gpu(bMatrix_dev,f_dev,bf_dev,N,nVar,nEl) &
      bind(c,name="BoundaryInterp_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: bMatrix_dev,f_dev,bf_dev
      integer(c_int),value :: N,nVar,nEl
    endsubroutine BoundaryInterp_3D_gpu
  endinterface

  interface
    subroutine Divergence_2D_gpu(f,df,dmat,N,nVar,nEl) &
      bind(c,name="Divergence_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,df,dmat
      integer(c_int),value :: N,nVar,nEl
    endsubroutine Divergence_2D_gpu
  endinterface

  ! MappedData

  ! Model
  interface
    subroutine UpdateSolution_gpu(solution,dsdt,dt,ndof) bind(c,name="UpdateSolution_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution,dsdt
      real(c_prec),value :: dt
      integer(c_int),value :: ndof
    endsubroutine UpdateSolution_gpu
  endinterface

  interface
    subroutine UpdateGRK_gpu(grk,solution,dsdt,rk_a,rk_g,dt,ndof) bind(c,name="UpdateGRK_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: grk,solution,dsdt
      real(c_prec),value :: rk_a,rk_g,dt
      integer(c_int),value :: ndof
    endsubroutine UpdateGRK_gpu
  endinterface

  interface
    subroutine CalculateDSDt_gpu(fluxDivergence,source,dsdt,ndof) bind(c,name="CalculateDSDt_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: fluxDivergence,source,dsdt
      integer(c_int),value :: ndof
    endsubroutine CalculateDSDt_gpu
  endinterface

  interface
    ! Fused tendency (source - fluxDivergence) + low-storage RK stage update,
    ! replacing CalculateDSDt_gpu + UpdateGRK_gpu.
    subroutine UpdateGRK_CalculateDSDt_gpu(grk,solution,fluxDivergence,source,rk_a,rk_g,dt,ndof) &
      bind(c,name="UpdateGRK_CalculateDSDt_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: grk,solution,fluxDivergence,source
      real(c_prec),value :: rk_a,rk_g,dt
      integer(c_int),value :: ndof
    endsubroutine UpdateGRK_CalculateDSDt_gpu
  endinterface

  interface
    ! Fused tendency + Euler update, replacing CalculateDSDt_gpu + UpdateSolution_gpu.
    subroutine UpdateSolution_CalculateDSDt_gpu(solution,fluxDivergence,source,dt,ndof) &
      bind(c,name="UpdateSolution_CalculateDSDt_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution,fluxDivergence,source
      real(c_prec),value :: dt
      integer(c_int),value :: ndof
    endsubroutine UpdateSolution_CalculateDSDt_gpu
  endinterface

  interface
    subroutine AccumulateField_gpu(a,b,ndof) bind(c,name="AccumulateField_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: a,b
      integer(c_int),value :: ndof
    endsubroutine AccumulateField_gpu
  endinterface

  interface
    subroutine GradientNormal_1D_gpu(fbn,fbavg,ndof) bind(c,name="GradientNormal_1d_gpu")
      use iso_c_binding
      type(c_ptr),value :: fbn,fbavg
      integer(c_int),value :: ndof
    endsubroutine GradientNormal_1D_gpu
  endinterface

  interface
    subroutine SideExchange_2D_gpu(extboundary,boundary,sideinfo,elemToRank,rankid,offset,n,nvar,nel) &
      bind(c,name="SideExchange_2D_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,sideinfo,elemToRank
      integer(c_int),value :: rankId,offset,N,nVar,nEl
    endsubroutine SideExchange_2D_gpu
  endinterface

  interface
    subroutine DG_BoundaryContribution_2D_gpu(bmatrix,qweights,bf,df,N,nvar,nel) &
      bind(c,name="DG_BoundaryContribution_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: bmatrix,qweights,bf,df
      integer(c_int),value :: N,nvar,nel
    endsubroutine DG_BoundaryContribution_2D_gpu
  endinterface

  interface
    subroutine JacobianWeight_2D_gpu(scalar,J,N,nVar,nEl) &
      bind(c,name="JacobianWeight_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: scalar,J
      integer(c_int),value :: N,nVar,nEl
    endsubroutine JacobianWeight_2D_gpu
  endinterface

  interface
    subroutine SideExchange_3D_gpu(extboundary,boundary,sideinfo,elemToRank,rankid,offset,n,nvar,nel) &
      bind(c,name="SideExchange_3D_gpu")
      use iso_c_binding
      type(c_ptr),value :: extboundary,boundary,sideinfo,elemToRank
      integer(c_int),value :: rankId,offset,N,nVar,nEl
    endsubroutine SideExchange_3D_gpu
  endinterface

  interface
    subroutine HaloPack_2D_gpu(boundary,sendbuf,halosides,n,nvar,nel,nhalo) &
      bind(c,name="HaloPack_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: boundary,sendbuf,halosides
      integer(c_int),value :: n,nvar,nel,nhalo
    endsubroutine HaloPack_2D_gpu
  endinterface

  interface
    subroutine HaloUnpack_2D_gpu(recvbuf,extboundary,halosides,n,nvar,nel,nhalo) &
      bind(c,name="HaloUnpack_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: recvbuf,extboundary,halosides
      integer(c_int),value :: n,nvar,nel,nhalo
    endsubroutine HaloUnpack_2D_gpu
  endinterface

  interface
    subroutine HaloPack_3D_gpu(boundary,sendbuf,halosides,n,nvar,nel,nhalo) &
      bind(c,name="HaloPack_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: boundary,sendbuf,halosides
      integer(c_int),value :: n,nvar,nel,nhalo
    endsubroutine HaloPack_3D_gpu
  endinterface

  interface
    subroutine HaloUnpack_3D_gpu(recvbuf,extboundary,halosides,n,nvar,nel,nhalo) &
      bind(c,name="HaloUnpack_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: recvbuf,extboundary,halosides
      integer(c_int),value :: n,nvar,nel,nhalo
    endsubroutine HaloUnpack_3D_gpu
  endinterface

  interface
    subroutine DG_BoundaryContribution_3D_gpu(bmatrix,qweights,bf,df,N,nvar,nel) &
      bind(c,name="DG_BoundaryContribution_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: bmatrix,qweights,bf,df
      integer(c_int),value :: N,nvar,nel
    endsubroutine DG_BoundaryContribution_3D_gpu
  endinterface

  interface
    ! Boundary contribution with the Jacobian weight (/jacobian) folded into the
    ! same pass -- replaces DG_BoundaryContribution_3D_gpu + JacobianWeight_3D_gpu.
    subroutine DG_BoundaryContribution_JacobianWeight_3D_gpu(bmatrix,qweights,bf,df,jacobian,N,nvar,nel) &
      bind(c,name="DG_BoundaryContribution_JacobianWeight_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: bmatrix,qweights,bf,df,jacobian
      integer(c_int),value :: N,nvar,nel
    endsubroutine DG_BoundaryContribution_JacobianWeight_3D_gpu
  endinterface

  interface
    subroutine JacobianWeight_3D_gpu(scalar,J,N,nVar,nEl) &
      bind(c,name="JacobianWeight_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: scalar,J
      integer(c_int),value :: N,nVar,nEl
    endsubroutine JacobianWeight_3D_gpu
  endinterface

  ! TwoPointVector

  interface
    subroutine TwoPointVectorDivergence_2D_gpu(f,df,dmat,N,nVar,nEl) &
      bind(c,name="TwoPointVectorDivergence_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,df,dmat
      integer(c_int),value :: N,nVar,nEl
    endsubroutine TwoPointVectorDivergence_2D_gpu
  endinterface

  interface
    subroutine TwoPointVectorDivergence_3D_gpu(f,df,dmat,N,nVar,nEl) &
      bind(c,name="TwoPointVectorDivergence_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,df,dmat
      integer(c_int),value :: N,nVar,nEl
    endsubroutine TwoPointVectorDivergence_3D_gpu
  endinterface

  interface
    subroutine MappedTwoPointVectorDivergence_2D_gpu(f,df,dmat,dsdx,jacobian,N,nVar,nEl) &
      bind(c,name="MappedTwoPointVectorDivergence_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,df,dmat,dsdx,jacobian
      integer(c_int),value :: N,nVar,nEl
    endsubroutine MappedTwoPointVectorDivergence_2D_gpu
  endinterface

  interface
    subroutine MappedTwoPointVectorDivergence_3D_gpu(f,df,dmat,dsdx,jacobian,N,nVar,nEl) &
      bind(c,name="MappedTwoPointVectorDivergence_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,df,dmat,dsdx,jacobian
      integer(c_int),value :: N,nVar,nEl
    endsubroutine MappedTwoPointVectorDivergence_3D_gpu
  endinterface

  interface
    subroutine ECDGSurfaceContribution_2D_gpu(fbn,jacobian,bmatrix,qweights,df,N,nvar,nel) &
      bind(c,name="ECDGSurfaceContribution_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: fbn,jacobian,bmatrix,qweights,df
      integer(c_int),value :: N,nvar,nel
    endsubroutine ECDGSurfaceContribution_2D_gpu
  endinterface

  interface
    subroutine ECDGSurfaceContribution_3D_gpu(fbn,jacobian,bmatrix,qweights,df,N,nvar,nel) &
      bind(c,name="ECDGSurfaceContribution_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: fbn,jacobian,bmatrix,qweights,df
      integer(c_int),value :: N,nvar,nel
    endsubroutine ECDGSurfaceContribution_3D_gpu
  endinterface

  ! MatrixMultiply : hand-written tensor-product contraction kernels that
  ! replace the cuBLAS/hipBLAS matrix operators (see SELF_MatrixMultiply.cpp).

  interface
    subroutine MatrixOp_1D_gpu(A,f,Af,opArows,opAcols,ncol) &
      bind(c,name="MatrixOp_1D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: A,f,Af
      integer(c_int),value :: opArows,opAcols,ncol
    endsubroutine MatrixOp_1D_gpu
  endinterface

  interface
    subroutine GridInterp_2D_gpu(A,f,fInterp,N,M,nVar,nEl) &
      bind(c,name="GridInterp_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: A,f,fInterp
      integer(c_int),value :: N,M,nVar,nEl
    endsubroutine GridInterp_2D_gpu
  endinterface

  interface
    subroutine GridInterp_3D_gpu(A,f,fInterp,N,M,nVar,nEl) &
      bind(c,name="GridInterp_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: A,f,fInterp
      integer(c_int),value :: N,M,nVar,nEl
    endsubroutine GridInterp_3D_gpu
  endinterface

  interface
    subroutine ScalarGradient_2D_gpu(A,f,df,N,nVar,nEl) &
      bind(c,name="ScalarGradient_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: A,f,df
      integer(c_int),value :: N,nVar,nEl
    endsubroutine ScalarGradient_2D_gpu
  endinterface

  interface
    subroutine ScalarGradient_3D_gpu(A,f,df,N,nVar,nEl) &
      bind(c,name="ScalarGradient_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: A,f,df
      integer(c_int),value :: N,nVar,nEl
    endsubroutine ScalarGradient_3D_gpu
  endinterface

  interface
    subroutine VectorDivergence_2D_gpu(A,f,df,N,nVar,nEl) &
      bind(c,name="VectorDivergence_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: A,f,df
      integer(c_int),value :: N,nVar,nEl
    endsubroutine VectorDivergence_2D_gpu
  endinterface

  interface
    subroutine VectorDivergence_3D_gpu(A,f,df,N,nVar,nEl) &
      bind(c,name="VectorDivergence_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: A,f,df
      integer(c_int),value :: N,nVar,nEl
    endsubroutine VectorDivergence_3D_gpu
  endinterface

  interface
    ! Fused contravariant-projection + interior divergence for 3-D mapped
    ! vectors. Handles the high-N LDS-overflow case internally (falls back to the
    ! two-kernel projection + divergence path), so callers use it unconditionally.
    ! When jacobian is non-null the /J weight is folded into the epilogue (used by
    ! the strong-form path); pass c_null_ptr for the DG path (which applies /J in
    ! DG_BoundaryContribution_JacobianWeight_3D_gpu after the boundary terms).
    subroutine MappedContravariantDivergence_3D_gpu(dsdx,A,f,df,jacobian,N,nVar,nEl) &
      bind(c,name="MappedContravariantDivergence_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: dsdx,A,f,df,jacobian
      integer(c_int),value :: N,nVar,nEl
    endsubroutine MappedContravariantDivergence_3D_gpu
  endinterface

endmodule SELF_GPUInterfaces
