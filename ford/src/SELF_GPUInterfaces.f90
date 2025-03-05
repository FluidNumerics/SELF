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

  interface
    subroutine Divergence_3D_gpu(f,df,dmat,N,nVar,nEl) &
      bind(c,name="Divergence_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: f,df,dmat
      integer(c_int),value :: N,nVar,nEl
    endsubroutine Divergence_3D_gpu
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
    subroutine ApplyFlip_2D_gpu(extBoundary,sideInfo,elemToRank,rankId,offset,N,nVar,nEl) &
      bind(c,name="ApplyFlip_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: extBoundary,sideInfo,elemToRank
      integer(c_int),value :: rankId,offset,N,nVar,nEl
    endsubroutine ApplyFlip_2D_gpu
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
    subroutine ApplyFlip_3D_gpu(extBoundary,sideInfo,elemToRank,rankId,offset,N,nVar,nEl) &
      bind(c,name="ApplyFlip_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: extBoundary,sideInfo,elemToRank
      integer(c_int),value :: rankId,offset,N,nVar,nEl
    endsubroutine ApplyFlip_3D_gpu
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
    subroutine JacobianWeight_3D_gpu(scalar,J,N,nVar,nEl) &
      bind(c,name="JacobianWeight_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: scalar,J
      integer(c_int),value :: N,nVar,nEl
    endsubroutine JacobianWeight_3D_gpu
  endinterface

endmodule SELF_GPUInterfaces
