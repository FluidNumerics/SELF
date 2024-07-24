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
    subroutine Average_gpu(favg,f1,f2,ndof)bind(c,name="Average_gpu")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: favg,f1,f2
        integer(c_int),value :: ndof
    end subroutine Average_gpu
  end interface

  interface
    subroutine BoundaryInterp_2D_gpu(bMatrix_dev,f_dev,bf_dev,N,nVar,nEl) &
        bind(c,name="BoundaryInterp_2D_gpu")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: bMatrix_dev,f_dev,bf_dev
        integer(c_int),value :: N,nVar,nEl
    end subroutine BoundaryInterp_2D_gpu
  end interface
  

  ! MappedData

  ! Model
  interface
    subroutine UpdateSolution_gpu(solution,dsdt,dt,ndof) bind(c,name="UpdateSolution_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution, dsdt
      real(c_prec),value :: dt
      integer(c_int),value :: ndof
    endsubroutine UpdateSolution_gpu
  endinterface

  interface
    subroutine UpdateGAB2_gpu(prevsol,solution,m,ndof) bind(c,name="UpdateGAB2_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution,prevsol
      integer(c_int),value :: m,ndof
    endsubroutine UpdateGAB2_gpu
  endinterface

  interface
    subroutine UpdateGAB3_gpu(prevsol,solution,m,ndof) bind(c,name="UpdateGAB3_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution,prevsol
      integer(c_int),value :: m,ndof
    endsubroutine UpdateGAB3_gpu
  endinterface

  interface
    subroutine UpdateGAB4_gpu(prevsol,solution,m,ndof) bind(c,name="UpdateGAB4_gpu")
      use iso_c_binding
      use SELF_Constants
      type(c_ptr),value :: solution,prevsol
      integer(c_int),value :: m,ndof
    endsubroutine UpdateGAB4_gpu
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
  
end module SELF_GPUInterfaces