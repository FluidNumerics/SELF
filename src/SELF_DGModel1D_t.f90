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

module SELF_DGModel1D_t

  use SELF_SupportRoutines
  use SELF_Metadata
  use SELF_Mesh_1D
  use SELF_MappedScalar_1D
  use SELF_HDF5
  use HDF5
  use FEQParse
  use SELF_Model

  implicit none

#include "SELF_Macros.h"

  type,extends(Model) :: DGModel1D_t
    type(MappedScalar1D) :: solution
    type(MappedScalar1D) :: solutionGradient
    type(MappedScalar1D) :: flux
    type(MappedScalar1D) :: source
    type(MappedScalar1D) :: fluxDivergence
    type(MappedScalar1D) :: dSdt
    type(MappedScalar1D) :: workSol
    type(Mesh1D),pointer :: mesh
    type(Geometry1D),pointer :: geometry

  contains

    procedure :: Init => Init_DGModel1D_t
    procedure :: SetMetadata => SetMetadata_DGModel1D_t
    procedure :: Free => Free_DGModel1D_t

    procedure :: CalculateEntropy => CalculateEntropy_DGModel1D_t
    procedure :: BoundaryFlux => BoundaryFlux_DGModel1D_t
    procedure :: FluxMethod => fluxmethod_DGModel1D_t
    procedure :: SourceMethod => sourcemethod_DGModel1D_t
    procedure :: SetBoundaryCondition => setboundarycondition_DGModel1D_t
    procedure :: SetGradientBoundaryCondition => setgradientboundarycondition_DGModel1D_t

    procedure :: UpdateSolution => UpdateSolution_DGModel1D_t

    procedure :: UpdateGRK2 => UpdateGRK2_DGModel1D_t
    procedure :: UpdateGRK3 => UpdateGRK3_DGModel1D_t
    procedure :: UpdateGRK4 => UpdateGRK4_DGModel1D_t

    procedure :: CalculateSolutionGradient => CalculateSolutionGradient_DGModel1D_t
    procedure :: CalculateTendency => CalculateTendency_DGModel1D_t

    generic :: SetSolution => SetSolutionFromChar_DGModel1D_t, &
      SetSolutionFromEqn_DGModel1D_t
    procedure,private :: SetSolutionFromChar_DGModel1D_t
    procedure,private :: SetSolutionFromEqn_DGModel1D_t

    procedure :: ReadModel => Read_DGModel1D_t
    procedure :: WriteModel => Write_DGModel1D_t
    procedure :: WriteTecplot => WriteTecplot_DGModel1D_t

  endtype DGModel1D_t

contains

  subroutine Init_DGModel1D_t(this,mesh,geometry)
    implicit none
    class(DGModel1D_t),intent(out) :: this
    type(Mesh1D),intent(in),target :: mesh
    type(Geometry1D),intent(in),target :: geometry

    this%mesh => mesh
    this%geometry => geometry
    call this%SetNumberOfVariables()

    call this%solution%Init(geometry%x%interp,this%nvar,this%mesh%nElem)
    call this%workSol%Init(geometry%x%interp,this%nvar,this%mesh%nElem)
    call this%dSdt%Init(geometry%x%interp,this%nvar,this%mesh%nElem)
    call this%solutionGradient%Init(geometry%x%interp,this%nvar,this%mesh%nElem)
    call this%flux%Init(geometry%x%interp,this%nvar,this%mesh%nElem)
    call this%source%Init(geometry%x%interp,this%nvar,this%mesh%nElem)
    call this%fluxDivergence%Init(geometry%x%interp,this%nvar,this%mesh%nElem)

    call this%solution%AssociateGeometry(geometry)
    call this%solutionGradient%AssociateGeometry(geometry)
    call this%flux%AssociateGeometry(geometry)
    call this%fluxDivergence%AssociateGeometry(geometry)

    call this%SetMetadata()

  endsubroutine Init_DGModel1D_t

  subroutine SetMetadata_DGModel1D_t(this)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    ! Local
    integer :: ivar
    character(LEN=3) :: ivarChar
    character(LEN=25) :: varname

    do ivar = 1,this%nvar
      write(ivarChar,'(I3.3)') ivar
      varname = "solution"//trim(ivarChar)
      call this%solution%SetName(ivar,varname)
      call this%solution%SetUnits(ivar,"[null]")
    enddo

  endsubroutine SetMetadata_DGModel1D_t

  subroutine Free_DGModel1D_t(this)
    implicit none
    class(DGModel1D_t),intent(inout) :: this

    call this%solution%DissociateGeometry()
    call this%solutionGradient%DissociateGeometry()
    call this%flux%DissociateGeometry()
    call this%fluxDivergence%DissociateGeometry()

    call this%solution%Free()
    call this%workSol%Free()
    call this%dSdt%Free()
    call this%solutionGradient%Free()
    call this%flux%Free()
    call this%source%Free()
    call this%fluxDivergence%Free()

  endsubroutine Free_DGModel1D_t

  subroutine SetSolutionFromEqn_DGModel1D_t(this,eqn)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    type(EquationParser),intent(in) :: eqn(1:this%solution%nVar)
    ! Local
    integer :: iVar

    ! Copy the equation parser
    do iVar = 1,this%solution%nVar
      call this%solution%SetEquation(ivar,eqn(iVar)%equation)
    enddo

    call this%solution%SetInteriorFromEquation(this%t)
    call this%solution%BoundaryInterp()

  endsubroutine SetSolutionFromEqn_DGModel1D_t

  subroutine SetSolutionFromChar_DGModel1D_t(this,eqnChar)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    character(LEN=SELF_EQUATION_LENGTH),intent(in) :: eqnChar(1:this%solution%nVar)
    ! Local
    integer :: iVar

    do iVar = 1,this%solution%nVar
      print*,iVar,eqnChar(iVar)
      call this%solution%SetEquation(ivar,eqnChar(iVar))
    enddo

    call this%solution%SetInteriorFromEquation(this%t)
    call this%solution%BoundaryInterp()

  endsubroutine SetSolutionFromChar_DGModel1D_t

  subroutine UpdateSolution_DGModel1D_t(this,dt)
    !! Computes a solution update as , where dt is either provided through the interface
    !! or taken as the Model's stored time step size (model % dt)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    real(prec),optional,intent(in) :: dt
    ! Local
    real(prec) :: dtLoc
    integer :: i,iEl,iVar

    if(present(dt)) then
      dtLoc = dt
    else
      dtLoc = this%dt
    endif

    do concurrent(i=1:this%solution%N+1,iel=1:this%mesh%nElem,ivar=1:this%solution%nVar)

      this%solution%interior(i,iEl,iVar) = &
        this%solution%interior(i,iEl,iVar)+ &
        dtLoc*this%dSdt%interior(i,iEl,iVar)

    enddo

  endsubroutine UpdateSolution_DGModel1D_t

  subroutine UpdateGRK2_DGModel1D_t(this,m)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,iEl,iVar

    do concurrent(i=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem,ivar=1:this%solution%nVar)

      this%workSol%interior(i,iEl,iVar) = rk2_a(m)* &
                                          this%workSol%interior(i,iEl,iVar)+ &
                                          this%dSdt%interior(i,iEl,iVar)

      this%solution%interior(i,iEl,iVar) = &
        this%solution%interior(i,iEl,iVar)+ &
        rk2_g(m)*this%dt*this%workSol%interior(i,iEl,iVar)

    enddo

  endsubroutine UpdateGRK2_DGModel1D_t

  subroutine UpdateGRK3_DGModel1D_t(this,m)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,iEl,iVar

    do concurrent(i=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem,ivar=1:this%solution%nVar)

      this%workSol%interior(i,iEl,iVar) = rk3_a(m)* &
                                          this%workSol%interior(i,iEl,iVar)+ &
                                          this%dSdt%interior(i,iEl,iVar)

      this%solution%interior(i,iEl,iVar) = &
        this%solution%interior(i,iEl,iVar)+ &
        rk3_g(m)*this%dt*this%workSol%interior(i,iEl,iVar)

    enddo

  endsubroutine UpdateGRK3_DGModel1D_t

  subroutine UpdateGRK4_DGModel1D_t(this,m)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,iEl,iVar

    do concurrent(i=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem,ivar=1:this%solution%nVar)

      this%workSol%interior(i,iEl,iVar) = rk4_a(m)* &
                                          this%workSol%interior(i,iEl,iVar)+ &
                                          this%dSdt%interior(i,iEl,iVar)

      this%solution%interior(i,iEl,iVar) = &
        this%solution%interior(i,iEl,iVar)+ &
        rk4_g(m)*this%dt*this%workSol%interior(i,iEl,iVar)

    enddo

  endsubroutine UpdateGRK4_DGModel1D_t

  subroutine CalculateSolutionGradient_DGModel1D_t(this)
    implicit none
    class(DGModel1D_t),intent(inout) :: this

    call this%solution%AverageSides()

    this%solution%boundarynormal(1,:,:) = -this%solution%avgBoundary(1,:,:) ! Account for left facing normal
    this%solution%boundarynormal(2,:,:) = this%solution%avgBoundary(2,:,:) ! Account for right facing normal

    call this%solution%MappedDGDerivative(this%solutionGradient%interior)

    ! interpolate the solutiongradient to the element boundaries
    call this%solutionGradient%BoundaryInterp()

    ! perform the side exchange to populate the
    ! solutionGradient % extBoundary attribute
    call this%solutionGradient%SideExchange(this%mesh)

  endsubroutine CalculateSolutionGradient_DGModel1D_t

  subroutine CalculateEntropy_DGModel1D_t(this)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    ! Local
    integer :: iel,i,ivar
    real(prec) :: e,s(1:this%solution%nvar),J

    e = 0.0_prec
    do iel = 1,this%geometry%nelem
      do i = 1,this%solution%interp%N+1
        J = this%geometry%dxds%interior(i,iel,1)
        s(1:this%solution%nvar) = this%solution%interior(i,iel,1:this%solution%nvar)
        e = e+this%entropy_func(s)*J
      enddo
    enddo

    this%entropy = e

  endsubroutine CalculateEntropy_DGModel1D_t

  subroutine setboundarycondition_DGModel1D_t(this)
    ! Here, we use the pre-tendency method to calculate the
    ! derivative of the solution using a bassi-rebay method
    ! We then do a boundary interpolation and side exchange
    ! on the gradient field
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    ! local
    integer :: ivar
    integer :: N,nelem
    real(prec) :: x

    nelem = this%geometry%nelem ! number of elements in the mesh
    N = this%solution%interp%N ! polynomial degree

    ! left-most boundary
    if(this%mesh%bcid(1) == SELF_BC_PRESCRIBED) then

      x = this%geometry%x%boundary(1,1,1)
      this%solution%extBoundary(1,1,1:this%nvar) = &
        this%hbc1d_Prescribed(x,this%t)

    elseif(this%mesh%bcid(1) == SELF_BC_RADIATION) then

      this%solution%extBoundary(1,1,1:this%nvar) = &
        this%hbc1d_Radiation(this%solution%boundary(1,1,1:this%nvar),-1.0_prec)

    elseif(this%mesh%bcid(1) == SELF_BC_NONORMALFLOW) then

      this%solution%extBoundary(1,1,1:this%nvar) = &
        this%hbc1d_NoNormalFlow(this%solution%boundary(1,1,1:this%nvar),-1.0_prec)

    else ! Periodic

      this%solution%extBoundary(1,1,1:this%nvar) = this%solution%boundary(2,nelem,1:this%nvar)

    endif

    ! right-most boundary
    if(this%mesh%bcid(1) == SELF_BC_PRESCRIBED) then

      x = this%geometry%x%boundary(2,nelem,1)
      this%solution%extBoundary(2,nelem,1:this%nvar) = &
        this%hbc1d_Prescribed(x,this%t)

    elseif(this%mesh%bcid(1) == SELF_BC_RADIATION) then

      this%solution%extBoundary(2,nelem,1:this%nvar) = &
        this%hbc1d_Radiation(this%solution%boundary(2,nelem,1:this%nvar),-1.0_prec)

    elseif(this%mesh%bcid(1) == SELF_BC_NONORMALFLOW) then

      this%solution%extBoundary(2,nelem,1:this%nvar) = &
        this%hbc1d_NoNormalFlow(this%solution%boundary(2,nelem,1:this%nvar),-1.0_prec)

    else ! Periodic

      this%solution%extBoundary(2,nelem,1:this%nvar) = this%solution%boundary(1,1,1:this%nvar)

    endif

  endsubroutine setboundarycondition_DGModel1D_t

  subroutine setgradientboundarycondition_DGModel1D_t(this)
    ! Here, we set the boundary conditions for the
    ! solution and the solution gradient at the left
    ! and right most boundaries.
    !
    ! Here, we use periodic boundary conditions
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    ! local
    real(prec) :: x
    integer :: nelem

    nelem = this%geometry%nelem ! number of elements in the mesh

    ! left-most boundary
    if(this%mesh%bcid(1) == SELF_BC_PRESCRIBED) then

      x = this%geometry%x%boundary(1,1,1)
      this%solutionGradient%extBoundary(1,1,1:this%nvar) = &
        this%pbc1d_Prescribed(x,this%t)

    elseif(this%mesh%bcid(1) == SELF_BC_RADIATION) then

      this%solutionGradient%extBoundary(1,1,1:this%nvar) = &
        this%pbc1d_Radiation(this%solutionGradient%boundary(1,1,1:this%nvar),-1.0_prec)

    elseif(this%mesh%bcid(1) == SELF_BC_NONORMALFLOW) then

      this%solutionGradient%extBoundary(1,1,1:this%nvar) = &
        this%pbc1d_NoNormalFlow(this%solutionGradient%boundary(1,1,1:this%nvar),-1.0_prec)

    else ! Periodic

      this%solutionGradient%extBoundary(1,1,1:this%nvar) = this%solutionGradient%boundary(2,nelem,1:this%nvar)

    endif

    ! right-most boundary
    if(this%mesh%bcid(1) == SELF_BC_PRESCRIBED) then

      x = this%geometry%x%boundary(2,nelem,1)
      this%solutionGradient%extBoundary(2,nelem,1:this%nvar) = &
        this%pbc1d_Prescribed(x,this%t)

    elseif(this%mesh%bcid(1) == SELF_BC_RADIATION) then

      this%solutionGradient%extBoundary(2,nelem,1:this%nvar) = &
        this%pbc1d_Radiation(this%solutionGradient%boundary(2,nelem,1:this%nvar),-1.0_prec)

    elseif(this%mesh%bcid(1) == SELF_BC_NONORMALFLOW) then

      this%solutionGradient%extBoundary(2,nelem,1:this%nvar) = &
        this%pbc1d_NoNormalFlow(this%solutionGradient%boundary(2,nelem,1:this%nvar),-1.0_prec)

    else ! Periodic

      this%solutionGradient%extBoundary(2,nelem,1:this%nvar) = this%solutionGradient%boundary(1,1,1:this%nvar)

    endif

  endsubroutine setgradientboundarycondition_DGModel1D_t

  subroutine BoundaryFlux_DGModel1D_t(this)
    ! this method uses an linear upwind solver for the
    ! advective flux and the bassi-rebay method for the
    ! diffusive fluxes
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    real(prec) :: fin(1:this%solution%nvar)
    real(prec) :: fout(1:this%solution%nvar)
    real(prec) :: dfdx(1:this%solution%nvar),nhat

    do concurrent(iside=1:2,iel=1:this%mesh%nElem)

      ! set the normal velocity
      if(iside == 1) then
        nhat = -1.0_prec
      else
        nhat = 1.0_prec
      endif

      fin = this%solution%boundary(iside,iel,1:this%solution%nvar) ! interior solution
      fout = this%solution%extboundary(iside,iel,1:this%solution%nvar) ! exterior solution
      dfdx = this%solutionGradient%avgboundary(iside,iel,1:this%solution%nvar) ! average solution gradient (with direction taken into account)
      this%flux%boundarynormal(iside,iel,1:this%solution%nvar) = &
        this%riemannflux1d(fin,fout,dfdx,nhat)

    enddo

  endsubroutine BoundaryFlux_DGModel1D_t

  subroutine fluxmethod_DGModel1D_t(this)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: i
    real(prec) :: f(1:this%solution%nvar),dfdx(1:this%solution%nvar)

    do concurrent(i=1:this%solution%N+1,iel=1:this%mesh%nElem)

      f = this%solution%interior(i,iel,1:this%solution%nvar)
      dfdx = this%solutionGradient%interior(i,iel,1:this%solution%nvar)

      this%flux%interior(i,iel,1:this%solution%nvar) = &
        this%flux1d(f,dfdx)

    enddo

  endsubroutine fluxmethod_DGModel1D_t

  subroutine sourcemethod_DGModel1D_t(this)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: i
    real(prec) :: f(1:this%solution%nvar),dfdx(1:this%solution%nvar)

    do concurrent(i=1:this%solution%N+1,iel=1:this%mesh%nElem)

      f = this%solution%interior(i,iel,1:this%solution%nvar)
      dfdx = this%solutionGradient%interior(i,iel,1:this%solution%nvar)

      this%source%interior(i,iel,1:this%solution%nvar) = &
        this%source1d(f,dfdx)

    enddo

  endsubroutine sourcemethod_DGModel1D_t

  subroutine CalculateTendency_DGModel1D_t(this)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    ! Local
    integer :: i,iEl,iVar

    call this%solution%BoundaryInterp()
    call this%solution%SideExchange(this%mesh)

    call this%PreTendency() ! User-supplied
    call this%SetBoundaryCondition() ! User-supplied

    if(this%gradient_enabled) then
      call this%CalculateSolutionGradient()
      call this%SetGradientBoundaryCondition() ! User-supplied
      call this%solutionGradient%AverageSides()
    endif

    call this%SourceMethod() ! User supplied
    call this%BoundaryFlux() ! User supplied
    call this%FluxMethod() ! User supplied

    call this%flux%MappedDGDerivative(this%fluxDivergence%interior)
    do concurrent(i=1:this%solution%N+1, &
                  iel=1:this%mesh%nElem,ivar=1:this%solution%nVar)

      this%dSdt%interior(i,iEl,iVar) = &
        this%source%interior(i,iEl,iVar)- &
        this%fluxDivergence%interior(i,iEl,iVar)

    enddo

  endsubroutine CalculateTendency_DGModel1D_t

  subroutine Write_DGModel1D_t(this,fileName)
#undef __FUNC__
#define __FUNC__ "Write_DGModel1D_t"
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    character(*),optional,intent(in) :: fileName
    ! Local
    integer(HID_T) :: fileId
    type(Scalar1D) :: solution
    type(Scalar1D) :: x
    type(Lagrange),target :: interp
    character(LEN=self_FileNameLength) :: pickupFile
    character(13) :: timeStampString

    write(timeStampString,'(I13.13)') this%ioIterate
    if(present(filename)) then
      pickupFile = trim(filename)
    else
      pickupFile = 'solution.'//timeStampString//'.h5'
    endif

    INFO("Writing pickup file : "//trim(pickupFile))
    call this%solution%UpdateHost()

    call Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId)

    ! Write the interpolant to the file
    call this%solution%interp%WriteHDF5(fileId)

    ! In this section, we write the solution and geometry on the control (quadrature) grid
    ! which can be used for model pickup runs or post-processing

    ! Write the model state to file
    call CreateGroup_HDF5(fileId,'/controlgrid')
    call this%solution%WriteHDF5(fileId,'/controlgrid/solution')

    ! Write the geometry to file
    call CreateGroup_HDF5(fileId,'/controlgrid/geometry')
    call this%geometry%x%WriteHDF5(fileId,'/controlgrid/geometry/x')
    ! -- END : writing solution on control grid -- !

    ! Interpolate the solution to a grid for plotting results
    ! Create an interpolant for the uniform grid
    call interp%Init(this%solution%interp%M, &
                     this%solution%interp%targetNodeType, &
                     this%solution%interp%N, &
                     this%solution%interp%controlNodeType)

    call solution%Init(interp, &
                       this%solution%nVar,this%solution%nElem)

    call x%Init(interp,1,this%solution%nElem)

    ! Map the mesh positions to the target grid
    call this%geometry%x%GridInterp(x%interior)

    ! Map the solution to the target grid
    call this%solution%GridInterp(solution%interior)

    ! Write the model state to file
    call CreateGroup_HDF5(fileId,'/targetgrid')
    call solution%WriteHDF5(fileId,'/targetgrid/solution')

    ! Write the geometry to file
    call CreateGroup_HDF5(fileId,'/targetgrid/geometry')
    call x%WriteHDF5(fileId,'/targetgrid/geometry/x')

    call Close_HDF5(fileId)

    call x%Free()
    call solution%Free()
    call interp%Free()

  endsubroutine Write_DGModel1D_t

  subroutine Read_DGModel1D_t(this,fileName)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    character(*),intent(in) :: fileName
    ! Local
    integer(HID_T) :: fileId
    integer(HID_T) :: solOffset(1:3)
    integer :: firstElem
    integer :: N

    call Open_HDF5(fileName,H5F_ACC_RDWR_F,fileId)
    call ReadArray_HDF5(fileId,'/controlgrid/solution/interior',this%solution%interior)
    call Close_HDF5(fileId)
    call this%solution%UpdateDevice()

  endsubroutine Read_DGModel1D_t

  subroutine WriteTecplot_DGModel1D_t(this,filename)
    implicit none
    class(DGModel1D_t),intent(inout) :: this
    character(*),intent(in),optional :: filename
    ! Local
    character(8) :: zoneID
    integer :: fUnit
    integer :: iEl,i,iVar
    character(LEN=self_FileNameLength) :: tecFile
    character(LEN=self_TecplotHeaderLength) :: tecHeader
    character(LEN=self_FormatLength) :: fmat
    character(13) :: timeStampString
    character(5) :: rankString
    type(Scalar1D) :: solution
    type(Scalar1D) :: x
    type(Lagrange),target :: interp

    if(present(filename)) then
      tecFile = filename
    else
      ! Create a 0-padded integer for the output iterate
      write(timeStampString,'(I13.13)') this%ioIterate
      ! Increment the ioIterate
      this%ioIterate = this%ioIterate+1
      tecFile = 'solution.'//timeStampString//'.curve'

    endif
    call this%solution%UpdateHost()
    ! Create an interpolant for the uniform grid
    call interp%Init(this%solution%interp%M, &
                     this%solution%interp%targetNodeType, &
                     this%solution%interp%N, &
                     this%solution%interp%controlNodeType)

    call solution%Init(interp, &
                       this%solution%nVar,this%solution%nElem)

    call x%Init(interp,1,this%solution%nElem)

    ! Map the mesh positions to the target grid
    call this%geometry%x%GridInterp(x%Interior)

    ! Map the solution to the target grid
    call this%solution%GridInterp(solution%interior)

    fmat = '(2(ES16.7E3,1x))'
    ! Let's write some tecplot!!
    open(UNIT=NEWUNIT(fUnit), &
         FILE=trim(tecFile), &
         FORM='formatted', &
         STATUS='replace')

    do iVar = 1,this%solution%nVar
      write(tecHeader,'(E15.6)') this%t
      tecHeader = "#TIME "//trim(tecHeader)
      write(fUnit,*) trim(tecHeader)

      tecHeader = "#"//trim(this%solution%meta(iVar)%name)//" vs position"
      write(fUnit,*) trim(tecHeader)
      do iEl = 1,this%solution%nElem

        do i = 1,this%solution%interp%M+1

          write(fUnit,fmat) x%interior(i,iEl,1), &
            solution%interior(i,iEl,iVar)

        enddo

      enddo
    enddo

    close(UNIT=fUnit)

    call x%Free()
    call solution%Free()
    call interp%Free()

  endsubroutine WriteTecplot_DGModel1D_t

endmodule SELF_DGModel1D_t
