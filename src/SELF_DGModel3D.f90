!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
module SELF_DGModel3D

  use SELF_SupportRoutines
  use SELF_Metadata
  use SELF_Mesh
  use SELF_MappedData_3D
  use SELF_HDF5
  use HDF5
  use FEQParse
  use SELF_Model

  implicit none

#include "SELF_Macros.h"

  type,extends(Model) :: DGModel3D
    type(MappedScalar3D)   :: solution
    type(MappedVector3D)   :: solutionGradient
    type(MappedVector3D)   :: flux
    type(MappedScalar3D)   :: source
    type(MappedScalar3D)   :: fluxDivergence
    type(MappedScalar3D)   :: dSdt
    type(MappedScalar3D)   :: workSol
    type(MappedScalar3D)   :: prevSol
    type(Mesh3D),pointer   :: mesh
    type(SEMHex),pointer  :: geometry

  contains

    procedure :: Init => Init_DGModel3D
    procedure :: Free => Free_DGModel3D

    procedure :: UpdateSolution => UpdateSolution_DGModel3D

    procedure :: ResizePrevSol => ResizePrevSol_DGModel3D

    procedure :: UpdateGAB2 => UpdateGAB2_DGModel3D
    procedure :: UpdateGAB3 => UpdateGAB3_DGModel3D
    procedure :: UpdateGAB4 => UpdateGAB4_DGModel3D

    procedure :: UpdateGRK2 => UpdateGRK2_DGModel3D
    procedure :: UpdateGRK3 => UpdateGRK3_DGModel3D
    procedure :: UpdateGRK4 => UpdateGRK4_DGModel3D

    procedure :: CalculateTendency => CalculateTendency_DGModel3D

    generic :: SetSolution => SetSolutionFromChar_DGModel3D, &
      SetSolutionFromEqn_DGModel3D
    procedure,private :: SetSolutionFromChar_DGModel3D
    procedure,private :: SetSolutionFromEqn_DGModel3D

    procedure :: ReadModel => Read_DGModel3D
    procedure :: WriteModel => Write_DGModel3D
    procedure :: WriteTecplot => WriteTecplot_DGModel3D

  endtype DGModel3D

contains

  subroutine Init_DGModel3D(this,nvar,mesh,geometry,decomp)
    implicit none
    class(DGModel3D),intent(out) :: this
    integer,intent(in) :: nvar
    type(Mesh3D),intent(in),target :: mesh
    type(SEMHex),intent(in),target :: geometry
    type(MPILayer),intent(in),target :: decomp
    ! Local
    integer :: ivar
    character(LEN=3) :: ivarChar
    character(LEN=25) :: varname

    this%decomp => decomp
    this%mesh => mesh
    this%geometry => geometry

    call this%solution%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%workSol%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%prevSol%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%dSdt%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%solutionGradient%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%flux%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%source%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%fluxDivergence%Init(geometry%x%interp,nVar,this%mesh%nElem)

    ! set default metadata
    do ivar = 1,nvar
      write(ivarChar,'(I3.3)') ivar
      varname = "solution"//trim(ivarChar)
      call this%solution%SetName(ivar,varname)
      call this%solution%SetUnits(ivar,"[null]")
    enddo

  endsubroutine Init_DGModel3D

  subroutine Free_DGModel3D(this)
    implicit none
    class(DGModel3D),intent(inout) :: this

    call this%solution%Free()
    call this%workSol%Free()
    call this%prevSol%Free()
    call this%dSdt%Free()
    call this%solutionGradient%Free()
    call this%flux%Free()
    call this%source%Free()
    call this%fluxDivergence%Free()

  endsubroutine Free_DGModel3D

  subroutine ResizePrevSol_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: nVar

    ! Free space, if necessary
    call this%prevSol%Free()

    ! Reallocate with increased variable dimension for
    ! storing "m" copies of solution data
    nVar = this%solution%nVar
    call this%prevSol%Init(this%geometry%x%interp,m*nVar,this%mesh%nElem)

  endsubroutine ResizePrevSol_DGModel3D

  subroutine SetSolutionFromEqn_DGModel3D(this,eqn)
    implicit none
    class(DGModel3D),intent(inout) :: this
    type(EquationParser),intent(in) :: eqn(1:this%solution%nVar)
    ! Local
    integer :: iVar

    ! Copy the equation parser
    do iVar = 1,this%solution%nVar
      call this%solution%SetEquation(ivar,eqn(iVar)%equation)
    enddo

    call this%solution%SetInteriorFromEquation(this%geometry,this%t)

    call this%solution%BoundaryInterp()

  endsubroutine SetSolutionFromEqn_DGModel3D

  subroutine SetSolutionFromChar_DGModel3D(this,eqnChar)
    implicit none
    class(DGModel3D),intent(inout) :: this
    character(*),intent(in) :: eqnChar(1:this%solution%nVar)
    ! Local
    integer :: iVar

    do iVar = 1,this%solution%nVar
      call this%solution%SetEquation(ivar,trim(eqnChar(iVar)))
    enddo

    call this%solution%SetInteriorFromEquation(this%geometry,this%t)

    call this%solution%BoundaryInterp()

  endsubroutine SetSolutionFromChar_DGModel3D

  subroutine UpdateSolution_DGModel3D(this,dt)
    !! Computes a solution update as , where dt is either provided through the interface
    !! or taken as the Model's stored time step size (model % dt)
    implicit none
    class(DGModel3D),intent(inout) :: this
    real(prec),optional,intent(in) :: dt
    ! Local
    real(prec) :: dtLoc
    integer :: i,j,k,iVar,iEl

    if(present(dt)) then
      dtLoc = dt
    else
      dtLoc = this%dt
    endif

    !$omp target map(to:this % dsdt % interior) map(tofrom:this % solution)
    !$omp teams distribute parallel do collapse(4) num_threads(256)
    do iEl = 1,this%solution%nElem
      do iVar = 1,this%solution%nVar
        do k = 1,this%solution%interp%N+1
          do j = 1,this%solution%interp%N+1
            do i = 1,this%solution%interp%N+1

              this%solution%interior(i,j,k,iEl,iVar) = &
                this%solution%interior(i,j,k,iEl,iVar)+ &
                dtLoc*this%dSdt%interior(i,j,k,iEl,iVar)

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine UpdateSolution_DGModel3D

  subroutine UpdateGAB2_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,j,k,nVar,iVar,iEl

    ! ab2_weight
    if(m == 0) then ! Initialization step - store the solution in the prevSol

      !$omp target map(tofrom: this % solution % interior) map(from:this % prevSol % interior)
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do k = 1,this%solution%interp%N+1
            do j = 1,this%solution%interp%N+1
              do i = 1,this%solution%interp%N+1

                this%prevSol%interior(i,j,k,iEl,iVar) = this%solution%interior(i,j,k,iEl,iVar)

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end target

    elseif(m == 1) then ! Reset solution

      !$omp target map(from: this % solution % interior) map(to:this % prevSol % interior)
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do k = 1,this%solution%interp%N+1
            do j = 1,this%solution%interp%N+1
              do i = 1,this%solution%interp%N+1

                this%solution%interior(i,j,k,iEl,iVar) = this%prevSol%interior(i,j,k,iEl,iVar)

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end target

    else ! Main looping section - nVar the previous solution, store the new solution, and
      ! create an interpolated solution to use for tendency calculation

      nVar = this%solution%nVar
      !$omp target map(tofrom: this % solution % interior, this % prevSol % interior)
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do k = 1,this%solution%interp%N+1
            do j = 1,this%solution%interp%N+1
              do i = 1,this%solution%interp%N+1

                ! Bump the last solution
                this%prevSol%interior(i,j,k,iEl,nVar+iVar) = this%prevSol%interior(i,j,k,iEl,iVar)

                ! Store the new solution
                this%prevSol%interior(i,j,k,iEl,iVar) = this%solution%interior(i,j,k,iEl,iVar)

                this%solution%interior(i,j,k,iEl,iVar) = &
                  1.5_prec*this%prevSol%interior(i,j,k,iEl,iVar)- &
                  0.5_prec*this%prevSol%interior(i,j,k,iEl,nVar+iVar)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end target

    endif

  endsubroutine UpdateGAB2_DGModel3D

  subroutine UpdateGAB3_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,j,k,nVar,iVar,iEl

    if(m == 0) then ! Initialization step - store the solution in the prevSol at nvar+ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      nVar = this%solution%nVar
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do k = 1,this%solution%interp%N+1
            do j = 1,this%solution%interp%N+1
              do i = 1,this%solution%interp%N+1

                this%prevSol%interior(i,j,k,iEl,nVar+iVar) = this%solution%interior(i,j,k,iEl,iVar)

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end target

    elseif(m == 1) then ! Initialization step - store the solution in the prevSol at ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do k = 1,this%solution%interp%N+1
            do j = 1,this%solution%interp%N+1
              do i = 1,this%solution%interp%N+1

                this%prevSol%interior(i,j,k,iEl,iVar) = this%solution%interior(i,j,k,iEl,iVar)

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end target

    elseif(m == 2) then ! Copy the solution back from the most recent prevsol

      !$omp target map(from: this % solution % interior) map(to: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do k = 1,this%solution%interp%N+1
            do j = 1,this%solution%interp%N+1
              do i = 1,this%solution%interp%N+1

                this%solution%interior(i,j,k,iEl,iVar) = this%prevSol%interior(i,j,k,iEl,iVar)

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end target

    else ! Main looping section - nVar the previous solution, store the new solution, and
      ! create an interpolated solution to use for tendency calculation

      nVar = this%solution%nVar
      !$omp target map(tofrom: this % solution % interior, this % prevSol % interior)
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do k = 1,this%solution%interp%N+1
            do j = 1,this%solution%interp%N+1
              do i = 1,this%solution%interp%N+1

                ! Bump the last two stored solutions
                this%prevSol%interior(i,j,k,iEl,2*nVar+iVar) = this%prevSol%interior(i,j,k,iEl,nVar+iVar)
                this%prevSol%interior(i,j,k,iEl,nVar+iVar) = this%prevSol%interior(i,j,k,iEl,iVar)

                ! Store the new solution
                this%prevSol%interior(i,j,k,iEl,iVar) = this%solution%interior(i,j,k,iEl,iVar)

                this%solution%interior(i,j,k,iEl,iVar) = &
                  (23.0_prec*this%prevSol%interior(i,j,k,iEl,iVar)- &
                   16.0_prec*this%prevSol%interior(i,j,k,iEl,nVar+iVar)+ &
                   5.0_prec*this%prevSol%interior(i,j,k,iEl,2*nVar+iVar))/12.0_prec

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end target

    endif

  endsubroutine UpdateGAB3_DGModel3D

  subroutine UpdateGAB4_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,j,k,nVar,iVar,iEl

    if(m == 0) then ! Initialization step - store the solution in the prevSol at nvar+ivar

      nVar = this%solution%nVar
      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do k = 1,this%solution%interp%N+1
            do j = 1,this%solution%interp%N+1
              do i = 1,this%solution%interp%N+1

                this%prevSol%interior(i,j,k,iEl,2*nVar+iVar) = this%solution%interior(i,j,k,iEl,iVar)

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end target

    elseif(m == 1) then ! Initialization step - store the solution in the prevSol at ivar

      nVar = this%solution%nVar
      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do k = 1,this%solution%interp%N+1
            do j = 1,this%solution%interp%N+1
              do i = 1,this%solution%interp%N+1

                this%prevSol%interior(i,j,k,iEl,nVar+iVar) = this%solution%interior(i,j,k,iEl,iVar)

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end target

    elseif(m == 2) then ! Initialization step - store the solution in the prevSol at ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do k = 1,this%solution%interp%N+1
            do j = 1,this%solution%interp%N+1
              do i = 1,this%solution%interp%N+1

                this%prevSol%interior(i,j,k,iEl,iVar) = this%solution%interior(i,j,k,iEl,iVar)

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end target

    elseif(m == 3) then ! Copy the solution back from the most recent prevsol

      !$omp target map(from: this % solution % interior) map(to: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do k = 1,this%solution%interp%N+1
            do j = 1,this%solution%interp%N+1
              do i = 1,this%solution%interp%N+1

                this%solution%interior(i,j,k,iEl,iVar) = this%prevSol%interior(i,j,k,iEl,iVar)

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end target

    else ! Main looping section - nVar the previous solution, store the new solution, and
      ! create an interpolated solution to use for tendency calculation

      nVar = this%solution%nVar
      !$omp target map(tofrom: this % solution % interior, this % prevSol % interior)
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do k = 1,this%solution%interp%N+1
            do j = 1,this%solution%interp%N+1
              do i = 1,this%solution%interp%N+1

                ! Bump the last two stored solutions
                this%prevSol%interior(i,j,k,iEl,3*nVar+iVar) = this%prevSol%interior(i,j,k,iEl,2*nVar+iVar)
                this%prevSol%interior(i,j,k,iEl,2*nVar+iVar) = this%prevSol%interior(i,j,k,iEl,nVar+iVar)
                this%prevSol%interior(i,j,k,iEl,nVar+iVar) = this%prevSol%interior(i,j,k,iEl,iVar)

                ! Store the new solution
                this%prevSol%interior(i,j,k,iEl,iVar) = this%solution%interior(i,j,k,iEl,iVar)

                this%solution%interior(i,j,k,iEl,iVar) = &
                  (55.0_prec*this%prevSol%interior(i,j,k,iEl,iVar)- &
                   59.0_prec*this%prevSol%interior(i,j,k,iEl,nVar+iVar)+ &
                   37.0_prec*this%prevSol%interior(i,j,k,iEl,2*nVar+iVar)- &
                   9.0_prec*this%prevSol%interior(i,j,k,iEl,3*nVar+iVar))/24.0_prec

              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end target

    endif

  endsubroutine UpdateGAB4_DGModel3D

  subroutine UpdateGRK2_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,j,k,iVar,iEl

    !$omp target map(tofrom: this % solution % interior, this % workSol % interior) map(to:this % dsdt % interior)
    !$omp teams distribute parallel do collapse(5) num_threads(256)
    do iEl = 1,this%solution%nElem
      do iVar = 1,this%solution%nVar
        do k = 1,this%solution%interp%N+1
          do j = 1,this%solution%interp%N+1
            do i = 1,this%solution%interp%N+1

              this%workSol%interior(i,j,k,iEl,iVar) = rk2_a(m)* &
                                                      this%workSol%interior(i,j,k,iEl,iVar)+ &
                                                      this%dSdt%interior(i,j,k,iEl,iVar)

              this%solution%interior(i,j,k,iEl,iVar) = &
                this%solution%interior(i,j,k,iEl,iVar)+ &
                rk2_g(m)*this%dt*this%workSol%interior(i,j,k,iEl,iVar)

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine UpdateGRK2_DGModel3D

  subroutine UpdateGRK3_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,j,k,iVar,iEl

    !$omp target map(tofrom: this % solution % interior, this % workSol % interior) map(to:this % dsdt % interior)
    !$omp teams distribute parallel do collapse(5) num_threads(256)
    do iEl = 1,this%solution%nElem
      do iVar = 1,this%solution%nVar
        do k = 1,this%solution%interp%N+1
          do j = 1,this%solution%interp%N+1
            do i = 1,this%solution%interp%N+1

              this%workSol%interior(i,j,k,iEl,iVar) = rk3_a(m)* &
                                                      this%workSol%interior(i,j,k,iEl,iVar)+ &
                                                      this%dSdt%interior(i,j,k,iEl,iVar)

              this%solution%interior(i,j,k,iEl,iVar) = &
                this%solution%interior(i,j,k,iEl,iVar)+ &
                rk3_g(m)*this%dt*this%workSol%interior(i,j,k,iEl,iVar)

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine UpdateGRK3_DGModel3D

  subroutine UpdateGRK4_DGModel3D(this,m)
    implicit none
    class(DGModel3D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,j,k,iVar,iEl

    !$omp target map(tofrom: this % solution % interior, this % workSol % interior) map(to:this % dsdt % interior)
    !$omp teams distribute parallel do collapse(5) num_threads(256)
    do iEl = 1,this%solution%nElem
      do iVar = 1,this%solution%nVar
        do k = 1,this%solution%interp%N+1
          do j = 1,this%solution%interp%N+1
            do i = 1,this%solution%interp%N+1

              this%workSol%interior(i,j,k,iEl,iVar) = rk4_a(m)* &
                                                      this%workSol%interior(i,j,k,iEl,iVar)+ &
                                                      this%dSdt%interior(i,j,k,iEl,iVar)

              this%solution%interior(i,j,k,iEl,iVar) = &
                this%solution%interior(i,j,k,iEl,iVar)+ &
                rk4_g(m)*this%dt*this%workSol%interior(i,j,k,iEl,iVar)

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine UpdateGRK4_DGModel3D

  subroutine CalculateTendency_DGModel3D(this)
    implicit none
    class(DGModel3D),intent(inout) :: this
    ! Local
    integer :: i,j,k,iVar,iEl

    call this%PreTendency()
    call this%solution%BoundaryInterp()
    call this%solution%SideExchange(this%mesh,this%decomp)
    call this%SetBoundaryCondition()
    call this%SourceMethod()
    call this%RiemannSolver()
    call this%FluxMethod()
    call this%flux%DGDivergence(this%geometry,this%fluxDivergence)

    !$omp target map(to: this % source, this % fluxDivergence) map(from:this % dSdt)
    !$omp teams distribute parallel do collapse(5) num_threads(256)
    do iEl = 1,this%solution%nElem
      do iVar = 1,this%solution%nVar
        do k = 1,this%solution%interp%N+1
          do j = 1,this%solution%interp%N+1
            do i = 1,this%solution%interp%N+1

              this%dSdt%interior(i,j,k,iEl,iVar) = &
                this%source%interior(i,j,k,iEl,iVar)- &
                this%fluxDivergence%interior(i,j,k,iEl,iVar)

            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine CalculateTendency_DGModel3D

  subroutine Write_DGModel3D(this,fileName)
#undef __FUNC__
#define __FUNC__ "Write_DGModel3D"
    implicit none
    class(DGModel3D),intent(inout) :: this
    character(*),optional,intent(in) :: fileName
    ! Local
    integer(HID_T) :: fileId
    type(Scalar3D) :: solution
    type(Vector3D) :: x
    type(Lagrange),target :: interp
    character(LEN=self_FileNameLength) :: pickupFile
    character(13) :: timeStampString

    if(present(filename)) then
      pickupFile = filename
    else
      write(timeStampString,'(I13.13)') this%ioIterate
      pickupFile = 'solution.'//timeStampString//'.h5'
    endif

    INFO("Writing pickup file : "//trim(pickupFile))

    if(this%decomp%mpiEnabled) then

      call Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId,this%decomp%mpiComm)

      ! Write the interpolant to the file
      INFO("Writing interpolant data to file")
      call this%solution%interp%WriteHDF5(fileId)

      ! In this section, we write the solution and geometry on the control (quadrature) grid
      ! which can be used for model pickup runs or post-processing
      ! Write the model state to file
      INFO("Writing control grid solution to file")
      call CreateGroup_HDF5(fileId,'/controlgrid')
      call this%solution%WriteHDF5(fileId,'/controlgrid/solution', &
                                   this%decomp%offsetElem(this%decomp%rankId),this%decomp%nElem)

      ! Write the geometry to file
      INFO("Writing control grid geometry to file")
      call CreateGroup_HDF5(fileId,'/controlgrid/geometry')
      call this%geometry%x%WriteHDF5(fileId,'/controlgrid/geometry/x', &
                                     this%decomp%offsetElem(this%decomp%rankId),this%decomp%nElem)

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
      call this%geometry%x%GridInterp(x)

      ! Map the solution to the target grid
      call this%solution%GridInterp(solution)

      ! Write the model state to file
      call CreateGroup_HDF5(fileId,'/targetgrid')
      call solution%WriteHDF5(fileId,'/targetgrid/solution', &
                              this%decomp%offsetElem(this%decomp%rankId),this%decomp%nElem)

      ! Write the geometry to file
      call CreateGroup_HDF5(fileId,'/targetgrid/mesh')
      call x%WriteHDF5(fileId,'/targetgrid/mesh/coords', &
                       this%decomp%offsetElem(this%decomp%rankId),this%decomp%nElem)

      call Close_HDF5(fileId)

    else

      call Open_HDF5(pickupFile,H5F_ACC_TRUNC_F,fileId)

      ! Write the interpolant to the file
      INFO("Writing interpolant data to file")
      call this%solution%interp%WriteHDF5(fileId)

      ! In this section, we write the solution and geometry on the control (quadrature) grid
      ! which can be used for model pickup runs or post-processing

      ! Write the model state to file
      INFO("Writing control grid solution to file")
      call CreateGroup_HDF5(fileId,'/controlgrid')
      call this%solution%WriteHDF5(fileId,'/controlgrid/solution')

      ! Write the geometry to file
      INFO("Writing control grid  geometry to file")
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
      call this%geometry%x%GridInterp(x)

      ! Map the solution to the target grid
      call this%solution%GridInterp(solution)

      ! Write the model state to file
      INFO("Writing target grid solution to file")
      call CreateGroup_HDF5(fileId,'/targetgrid')
      call solution%WriteHDF5(fileId,'/targetgrid/solution')

      ! Write the geometry to file
      INFO("Writing target grid geometry to file")
      call CreateGroup_HDF5(fileId,'/targetgrid/geometry')
      call x%WriteHDF5(fileId,'/targetgrid/geometry/x')

      call Close_HDF5(fileId)

    endif

    call x%Free()
    call solution%Free()
    call interp%Free()

  endsubroutine Write_DGModel3D

  subroutine Read_DGModel3D(this,fileName)
    implicit none
    class(DGModel3D),intent(inout) :: this
    character(*),intent(in) :: fileName
    ! Local
    integer(HID_T) :: fileId
    integer(HID_T) :: solOffset(1:5)
    integer :: firstElem
    integer :: N

    if(this%decomp%mpiEnabled) then
      call Open_HDF5(fileName,H5F_ACC_RDWR_F,fileId, &
                     this%decomp%mpiComm)
    else
      call Open_HDF5(fileName,H5F_ACC_RDWR_F,fileId)
    endif

    ! CALL ReadAttribute_HDF5(fileId,'N',N)

    ! IF (this % solution % interp % N /= N) THEN
    !   STOP 'Error : Solution polynomial degree does not match input file'
    ! END IF

    if(this%decomp%mpiEnabled) then
      firstElem = this%decomp%offsetElem(this%decomp%rankId)+1
      solOffset(1:5) = (/0,0,0,1,firstElem/)
      call ReadArray_HDF5(fileId,'/controlgrid/solution/interior', &
                          this%solution%interior,solOffset)
    else
      call ReadArray_HDF5(fileId,'/controlgrid/solution/interior',this%solution%interior)
    endif

    call Close_HDF5(fileId)

  endsubroutine Read_DGModel3D

  subroutine WriteTecplot_DGModel3D(this,filename)
    implicit none
    class(DGModel3D),intent(inout) :: this
    character(*),intent(in),optional :: filename
    ! Local
    character(8) :: zoneID
    integer :: fUnit
    integer :: iEl,i,j,k,iVar
    character(LEN=self_FileNameLength) :: tecFile
    character(LEN=self_TecplotHeaderLength) :: tecHeader
    character(LEN=self_FormatLength) :: fmat
    character(13) :: timeStampString
    character(5) :: rankString
    type(Scalar3D) :: solution
    type(Vector3D) :: solutionGradient
    type(Vector3D) :: x
    type(Lagrange),target :: interp

    if(present(filename)) then
      tecFile = filename
    else
      write(timeStampString,'(I13.13)') this%ioIterate

      if(this%decomp%mpiEnabled) then
        write(rankString,'(I5.5)') this%decomp%rankId
        tecFile = 'solution.'//rankString//'.'//timeStampString//'.tec'
      else
        tecFile = 'solution.'//timeStampString//'.tec'
      endif

    endif

    ! Create an interpolant for the uniform grid
    call interp%Init(this%solution%interp%M, &
                     this%solution%interp%targetNodeType, &
                     this%solution%interp%N, &
                     this%solution%interp%controlNodeType)

    call solution%Init(interp, &
                       this%solution%nVar,this%solution%nElem)

    call solutionGradient%Init(interp, &
                               this%solution%nVar,this%solution%nElem)

    call x%Init(interp,1,this%solution%nElem)

    ! Map the mesh positions to the target grid
    call this%geometry%x%GridInterp(x)

    ! Map the solution to the target grid
    call this%solution%GridInterp(solution)

    ! Map the solution to the target grid
    call this%solutionGradient%GridInterp(solutionGradient)

    open(UNIT=NEWUNIT(fUnit), &
         FILE=trim(tecFile), &
         FORM='formatted', &
         STATUS='replace')

    tecHeader = 'VARIABLES = "X", "Y", "Z"'
    do iVar = 1,this%solution%nVar
      tecHeader = trim(tecHeader)//', "'//trim(this%solution%meta(iVar)%name)//'"'
    enddo

    do iVar = 1,this%solution%nVar
      tecHeader = trim(tecHeader)//', "d/dx('//trim(this%solution%meta(iVar)%name)//')"'
    enddo

    do iVar = 1,this%solution%nVar
      tecHeader = trim(tecHeader)//', "d/dy('//trim(this%solution%meta(iVar)%name)//')"'
    enddo

    write(fUnit,*) trim(tecHeader)

    ! Create format statement
    write(fmat,*) 3*this%solution%nvar+3
    fmat = '('//trim(fmat)//'(ES16.7E3,1x))'

    do iEl = 1,this%solution%nElem

      ! TO DO :: Get the global element ID
      write(zoneID,'(I8.8)') iEl
      write(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',this%solution%interp%M+1, &
        ', J=',this%solution%interp%M+1

      do k = 1,this%solution%interp%M+1
        do j = 1,this%solution%interp%M+1
          do i = 1,this%solution%interp%M+1

            write(fUnit,fmat) x%interior(1,i,j,k,iEl,1), &
              x%interior(2,i,j,k,iEl,1), &
              x%interior(3,i,j,k,iEl,1), &
              solution%interior(i,j,k,iEl,1:this%solution%nvar), &
              solutionGradient%interior(1,i,j,k,iEl,1:this%solution%nvar), &
              solutionGradient%interior(2,i,j,k,iEl,1:this%solution%nvar), &
              solutionGradient%interior(3,i,j,k,iEl,1:this%solution%nvar)

          enddo
        enddo
      enddo

    enddo

    close(UNIT=fUnit)

    call x%Free()
    call solution%Free()
    call interp%Free()

  endsubroutine WriteTecplot_DGModel3D

endmodule SELF_DGModel3D
