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

program LinearEuler3D_Pickup_Legacy_Sigma
!! Regression test for reading a 3D linear Euler pickup file written before the
!! relaxation rate sigma became a solution variable.
!!
!! Read_DGModel3D_t restores the solution one dataset per variable, keyed on the
!! variable's metadata name. A pickup file written by the six-variable model
!! (u, v, w, P, c, rho0) has no "sigma" dataset, so the reader has to tolerate
!! its absence rather than read through an invalid dataset id. Read_DGModel3D_t
!! is a separate routine from its 2-D counterpart, so it needs its own coverage.
!!
!! A legacy file is synthesised the only way that is faithful to one: write a
!! pickup with the current model, then delete /controlgrid/solution/sigma from
!! it. Reading that file back must restore the six variables it does hold and
!! leave sigma at the value it was initialized with.
!!
!! Assertions:
!!  (a) ReadModel completes on a file with no sigma dataset;
!!  (b) the six variables the file holds are restored bit for bit;
!!  (c) sigma is left at its pre-read value (zero from Init), which is the
!!      undamped configuration - a legacy run restarts with the same dynamics
!!      it was written with.

  use self_data
  use self_lineareuler3d
  use self_mesh_3d
  use HDF5

  implicit none

  integer,parameter :: controlDegree = 3
  integer,parameter :: targetDegree = 6
  integer,parameter :: nElemPerDim = 2
  real(prec),parameter :: dx = 0.1_prec
  real(prec),parameter :: rho0 = 1.2_prec
  real(prec),parameter :: c0 = 2.5_prec
  real(prec),parameter :: sigma_written = 7.0_prec
  character(*),parameter :: pickupFile = "lineareuler3d-legacy-sigma.pickup.h5"

  type(LinearEuler3D) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh3D),target :: mesh
  type(SEMHex),target :: geometry
  integer :: bcids(1:6)
  integer :: i,j,k,iel,ivar
  integer(HID_T) :: fileId
  integer :: hdferr
  real(prec) :: maxdiff,maxsigma
  real(prec),allocatable :: expected(:,:,:,:,:)

  bcids(1:6) = SELF_BC_RADIATION

  call mesh%StructuredMesh(nElemPerDim,nElemPerDim,nElemPerDim, &
                           1,1,1,dx,dx,dx,bcids)

  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  call modelobj%Init(mesh,geometry)
  modelobj%prescribed_bcs_enabled = .false.
  modelobj%tecplot_enabled = .false.

  ! A state whose every variable is distinct per node, so a dataset restored
  ! into the wrong variable slot would be caught.
  do iel = 1,mesh%nElem
    do k = 1,modelobj%solution%N+1
      do j = 1,modelobj%solution%N+1
        do i = 1,modelobj%solution%N+1
          modelobj%solution%interior(i,j,k,iel,1) = 0.01_prec*real(i+j+k+iel,prec) ! u
          modelobj%solution%interior(i,j,k,iel,2) = 0.02_prec*real(i+2*j+k+iel,prec) ! v
          modelobj%solution%interior(i,j,k,iel,3) = 0.03_prec*real(i+j+2*k+iel,prec) ! w
          modelobj%solution%interior(i,j,k,iel,4) = 0.04_prec*real(i+j+k+2*iel,prec) ! P
          modelobj%solution%interior(i,j,k,iel,5) = c0 ! c
          modelobj%solution%interior(i,j,k,iel,6) = rho0 ! rho0
          modelobj%solution%interior(i,j,k,iel,7) = sigma_written ! sigma
        enddo
      enddo
    enddo
  enddo
  call modelobj%solution%UpdateDevice()

  allocate(expected,mold=modelobj%solution%interior)
  expected = modelobj%solution%interior

  call modelobj%WriteModel(pickupFile)

  ! Turn the file into one the six-variable model would have written.
  call h5open_f(hdferr)
  call h5fopen_f(pickupFile,H5F_ACC_RDWR_F,fileId,hdferr)
  if(hdferr /= 0) then
    print*,"Error: could not reopen the pickup file to remove the sigma dataset"
    stop 1
  endif
  call h5ldelete_f(fileId,"/controlgrid/solution/sigma",hdferr)
  if(hdferr /= 0) then
    print*,"Error: could not delete /controlgrid/solution/sigma"
    stop 1
  endif
  call h5fclose_f(fileId,hdferr)
  call h5close_f(hdferr)

  ! Clear the solution so a variable that fails to be restored shows up as a
  ! mismatch rather than as a leftover of what was written.
  modelobj%solution%interior = 0.0_prec
  call modelobj%solution%UpdateDevice()

  call modelobj%ReadModel(pickupFile)
  call modelobj%solution%UpdateHost()

  ! (b) the six variables present in the file round-trip exactly
  maxdiff = 0.0_prec
  do ivar = 1,6
    maxdiff = max(maxdiff, &
                  maxval(abs(modelobj%solution%interior(:,:,:,:,ivar)- &
                             expected(:,:,:,:,ivar))))
  enddo
  print*,"max |restored - written| over variables 1:6 :",maxdiff

  ! (c) sigma is untouched by the read
  maxsigma = maxval(abs(modelobj%solution%interior(:,:,:,:,7)))
  print*,"max |sigma| after reading a file without it   :",maxsigma

  if(maxdiff /= 0.0_prec) then
    print*,"Error: variables present in the legacy pickup file were not restored exactly ",maxdiff
    stop 1
  endif
  if(maxsigma /= 0.0_prec) then
    print*,"Error: sigma was modified by reading a pickup file that has no sigma dataset ",maxsigma
    stop 1
  endif

  deallocate(expected)
  call modelobj%Free()
  call geometry%Free()
  call mesh%Free()
  call interp%Free()

endprogram LinearEuler3D_Pickup_Legacy_Sigma
