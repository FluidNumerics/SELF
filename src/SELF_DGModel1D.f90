!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
module SELF_DGModel1D

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

  type,extends(Model) :: DGModel1D
    type(MappedScalar1D) :: solution
    type(MappedScalar1D) :: solutionGradient
    type(MappedScalar1D) :: velocity
    type(MappedScalar1D) :: flux
    type(MappedScalar1D) :: source
    type(MappedScalar1D) :: fluxDivergence
    type(MappedScalar1D) :: dSdt
    type(MappedScalar1D) :: workSol
    type(MappedScalar1D) :: prevSol
    type(Mesh1D),pointer :: mesh
    type(Geometry1D),pointer :: geometry

  contains

    procedure :: Init => Init_DGModel1D
    procedure :: Free => Free_DGModel1D

    procedure :: UpdateSolution => UpdateSolution_DGModel1D

    procedure :: ResizePrevSol => ResizePrevSol_DGModel1D

    procedure :: UpdateGAB2 => UpdateGAB2_DGModel1D
    procedure :: UpdateGAB3 => UpdateGAB3_DGModel1D
    procedure :: UpdateGAB4 => UpdateGAB4_DGModel1D

    procedure :: UpdateGRK2 => UpdateGRK2_DGModel1D
    procedure :: UpdateGRK3 => UpdateGRK3_DGModel1D
    procedure :: UpdateGRK4 => UpdateGRK4_DGModel1D
    procedure :: CalculateTendency => CalculateTendency_DGModel1D

    generic :: SetSolution => SetSolutionFromChar_DGModel1D, &
      SetSolutionFromEqn_DGModel1D
    procedure,private :: SetSolutionFromChar_DGModel1D
    procedure,private :: SetSolutionFromEqn_DGModel1D

    procedure :: ReadModel => Read_DGModel1D
    procedure :: WriteModel => Write_DGModel1D
    procedure :: WriteTecplot => WriteTecplot_DGModel1D

  endtype DGModel1D

contains

  subroutine Init_DGModel1D(this,nvar,mesh,geometry,decomp)
    implicit none
    class(DGModel1D),intent(out) :: this
    integer,intent(in) :: nvar
    type(Mesh1D),intent(in),target :: mesh
    type(Geometry1D),intent(in),target :: geometry
    type(MPILayer),intent(in),target :: decomp

    this%decomp => decomp
    this%mesh => mesh
    this%geometry => geometry

    call this%solution%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%workSol%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%prevSol%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%velocity%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%dSdt%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%solutionGradient%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%flux%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%source%Init(geometry%x%interp,nVar,this%mesh%nElem)
    call this%fluxDivergence%Init(geometry%x%interp,nVar,this%mesh%nElem)

  endsubroutine Init_DGModel1D

  subroutine Free_DGModel1D(this)
    implicit none
    class(DGModel1D),intent(inout) :: this

    call this%solution%Free()
    call this%workSol%Free()
    call this%prevSol%Free()
    call this%velocity%Free()
    call this%dSdt%Free()
    call this%solutionGradient%Free()
    call this%flux%Free()
    call this%source%Free()
    call this%fluxDivergence%Free()

  endsubroutine Free_DGModel1D

  subroutine ResizePrevSol_DGModel1D(this,m)
    implicit none
    class(DGModel1D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: nVar

    ! Free space, if necessary
    call this%prevSol%Free()

    ! Reallocate with increased variable dimension for
    ! storing "m" copies of solution data
    nVar = this%solution%nVar
    call this%prevSol%Init(this%geometry%x%interp,m*nVar,this%mesh%nElem)

  endsubroutine ResizePrevSol_DGModel1D

  subroutine SetSolutionFromEqn_DGModel1D(this,eqn)
    implicit none
    class(DGModel1D),intent(inout) :: this
    type(EquationParser),intent(in) :: eqn(1:this%solution%nVar)
    ! Local
    integer :: iVar

    ! Copy the equation parser
    do iVar = 1,this%solution%nVar
      call this%solution%SetEquation(ivar,eqn(iVar)%equation)
    enddo

    call this%solution%SetInteriorFromEquation(this%geometry,this%t)
    call this%solution%BoundaryInterp()

  endsubroutine SetSolutionFromEqn_DGModel1D

  subroutine SetSolutionFromChar_DGModel1D(this,eqnChar)
    implicit none
    class(DGModel1D),intent(inout) :: this
    character(LEN=SELF_EQUATION_LENGTH),intent(in) :: eqnChar(1:this%solution%nVar)
    ! Local
    integer :: iVar

    do iVar = 1,this%solution%nVar
      print*,iVar,eqnChar(iVar)
      call this%solution%SetEquation(ivar,eqnChar(iVar))
    enddo

    call this%solution%SetInteriorFromEquation(this%geometry,this%t)
    call this%solution%BoundaryInterp()

  endsubroutine SetSolutionFromChar_DGModel1D

  subroutine UpdateSolution_DGModel1D(this,dt)
    !! Computes a solution update as , where dt is either provided through the interface
    !! or taken as the Model's stored time step size (model % dt)
    implicit none
    class(DGModel1D),intent(inout) :: this
    real(prec),optional,intent(in) :: dt
    ! Local
    real(prec) :: dtLoc
    integer :: i,iEl,iVar

    if(present(dt)) then
      dtLoc = dt
    else
      dtLoc = this%dt
    endif

    !$omp target map(to:this % dsdt % interior) map(tofrom:this % solution)
    !$omp teams distribute parallel do collapse(3) num_threads(256)
    do iEl = 1,this%solution%nElem
      do iVar = 1,this%solution%nVar
        do i = 1,this%solution%interp%N+1

          this%solution%interior(i,iEl,iVar) = &
            this%solution%interior(i,iEl,iVar)+ &
            dtLoc*this%dSdt%interior(i,iEl,iVar)

        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine UpdateSolution_DGModel1D

  subroutine UpdateGAB2_DGModel1D(this,m)
    implicit none
    class(DGModel1D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,nVar,iEl,iVar

    ! ab2_weight
    if(m == 0) then ! Initialization step - store the solution in the prevSol

      !$omp target map(tofrom: this % solution % interior) map(from:this % prevSol % interior)
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do i = 1,this%solution%interp%N+1

            this%prevSol%interior(i,iEl,iVar) = this%solution%interior(i,iEl,iVar)

          enddo
        enddo
      enddo
      !$omp end target

    elseif(m == 1) then ! Copy the solution back from prevsol

      !$omp target map(from: this % solution % interior) map(to:this % prevSol % interior)
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do i = 1,this%solution%interp%N+1

            this%solution%interior(i,iEl,iVar) = this%prevSol%interior(i,iEl,iVar)

          enddo
        enddo
      enddo
      !$omp end target

    else ! Main looping section - nVar the previous solution, store the new solution, and
      ! create an interpolated solution to use for tendency calculation

      nVar = this%solution%nVar
      !$omp target map(tofrom: this % solution % interior, this % prevSol % interior)
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do i = 1,this%solution%interp%N+1

            ! Bump the last solution
            this%prevSol%interior(i,iEl,nVar+iVar) = this%prevSol%interior(i,iEl,iVar)

            ! Store the new solution
            this%prevSol%interior(i,iEl,iVar) = this%solution%interior(i,iEl,iVar)

            this%solution%interior(i,iEl,iVar) = &
              1.5_prec*this%prevSol%interior(i,iEl,iVar)- &
              0.5_prec*this%prevSol%interior(i,iEl,nVar+iVar)
          enddo
        enddo
      enddo
      !$omp end target
    endif

  endsubroutine UpdateGAB2_DGModel1D

  subroutine UpdateGAB3_DGModel1D(this,m)
    implicit none
    class(DGModel1D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,nVar,iEl,iVar

    if(m == 0) then ! Initialization step - store the solution in the prevSol at nvar+ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      nVar = this%solution%nVar
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do i = 1,this%solution%interp%N+1

            this%prevSol%interior(i,iEl,nVar+iVar) = this%solution%interior(i,iEl,iVar)

          enddo
        enddo
      enddo

    elseif(m == 1) then ! Initialization step - store the solution in the prevSol at ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do i = 1,this%solution%interp%N+1

            this%prevSol%interior(i,iEl,iVar) = this%solution%interior(i,iEl,iVar)

          enddo
        enddo
      enddo
      !$omp end target

    elseif(m == 2) then ! Copy the solution back from the most recent prevsol

      !$omp target map(from: this % solution % interior) map(to: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do i = 1,this%solution%interp%N+1

            this%solution%interior(i,iEl,iVar) = this%prevSol%interior(i,iEl,iVar)

          enddo
        enddo
      enddo
      !$omp end target

    else ! Main looping section - nVar the previous solution, store the new solution, and
      ! create an interpolated solution to use for tendency calculation

      !$omp target map(tofrom: this % solution % interior, this % prevSol % interior)
      nVar = this%solution%nVar
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do i = 1,this%solution%interp%N+1

            ! Bump the last two stored solutions
            nVar = this%solution%nVar
            this%prevSol%interior(i,iEl,2*nVar+iVar) = this%prevSol%interior(i,iEl,nVar+iVar)
            this%prevSol%interior(i,iEl,nVar+iVar) = this%prevSol%interior(i,iEl,iVar)

            ! Store the new solution
            this%prevSol%interior(i,iEl,iVar) = this%solution%interior(i,iEl,iVar)

            this%solution%interior(i,iEl,iVar) = &
              (23.0_prec*this%prevSol%interior(i,iEl,iVar)- &
               16.0_prec*this%prevSol%interior(i,iEl,nVar+iVar)+ &
               5.0_prec*this%prevSol%interior(i,iEl,2*nVar+iVar))/12.0_prec

          enddo
        enddo
      enddo
      !$omp end target

    endif

  endsubroutine UpdateGAB3_DGModel1D

  subroutine UpdateGAB4_DGModel1D(this,m)
    implicit none
    class(DGModel1D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,nVar,iEl,iVar

    if(m == 0) then ! Initialization step - store the solution in the prevSol at nvar+ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      nVar = this%solution%nVar
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do i = 1,this%solution%interp%N+1

            this%prevSol%interior(i,iEl,2*nVar+iVar) = this%solution%interior(i,iEl,iVar)

          enddo
        enddo
      enddo
      !$omp end target

    elseif(m == 1) then ! Initialization step - store the solution in the prevSol at ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      nVar = this%solution%nVar
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do i = 1,this%solution%interp%N+1

            this%prevSol%interior(i,iEl,nVar+iVar) = this%solution%interior(i,iEl,iVar)

          enddo
        enddo
      enddo
      !$omp end target

    elseif(m == 2) then ! Initialization step - store the solution in the prevSol at ivar

      !$omp target map(to: this % solution % interior) map(from: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do i = 1,this%solution%interp%N+1

            this%prevSol%interior(i,iEl,iVar) = this%solution%interior(i,iEl,iVar)

          enddo
        enddo
      enddo
      !$omp end target

    elseif(m == 3) then ! Copy the solution back from the most recent prevsol

      !$omp target map(from: this % solution % interior) map(to: this % prevSol % interior)
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do i = 1,this%solution%interp%N+1

            this%solution%interior(i,iEl,iVar) = this%prevSol%interior(i,iEl,iVar)

          enddo
        enddo
      enddo

    else ! Main looping section - nVar the previous solution, store the new solution, and
      ! create an interpolated solution to use for tendency calculation

      !$omp target map(tofrom: this % solution % interior, this % prevSol % interior)
      nVar = this%solution%nVar
      !$omp teams distribute parallel do collapse(3) num_threads(256)
      do iEl = 1,this%solution%nElem
        do iVar = 1,this%solution%nVar
          do i = 1,this%solution%interp%N+1

            ! Bump the last two stored solutions
            this%prevSol%interior(i,iEl,3*nVar+iVar) = this%prevSol%interior(i,iEl,2*nVar+iVar)
            this%prevSol%interior(i,iEl,2*nVar+iVar) = this%prevSol%interior(i,iEl,nVar+iVar)
            this%prevSol%interior(i,iEl,nVar+iVar) = this%prevSol%interior(i,iEl,iVar)

            ! Store the new solution
            this%prevSol%interior(i,iEl,iVar) = this%solution%interior(i,iEl,iVar)

            this%solution%interior(i,iEl,iVar) = &
              (55.0_prec*this%prevSol%interior(i,iEl,iVar)- &
               59.0_prec*this%prevSol%interior(i,iEl,nVar+iVar)+ &
               37.0_prec*this%prevSol%interior(i,iEl,2*nVar+iVar)- &
               9.0_prec*this%prevSol%interior(i,iEl,3*nVar+iVar))/24.0_prec

          enddo
        enddo
      enddo
      !$omp end target

    endif

  endsubroutine UpdateGAB4_DGModel1D

  subroutine UpdateGRK2_DGModel1D(this,m)
    implicit none
    class(DGModel1D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,iEl,iVar

    !$omp target map(tofrom: this % solution % interior, this % workSol % interior) map(to:this % dsdt % interior)
    !$omp teams distribute parallel do collapse(3) num_threads(256)
    do iEl = 1,this%solution%nElem
      do iVar = 1,this%solution%nVar
        do i = 1,this%solution%interp%N+1

          this%workSol%interior(i,iEl,iVar) = rk2_a(m)* &
                                              this%workSol%interior(i,iEl,iVar)+ &
                                              this%dSdt%interior(i,iEl,iVar)

          this%solution%interior(i,iEl,iVar) = &
            this%solution%interior(i,iEl,iVar)+ &
            rk2_g(m)*this%dt*this%workSol%interior(i,iEl,iVar)

        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine UpdateGRK2_DGModel1D

  subroutine UpdateGRK3_DGModel1D(this,m)
    implicit none
    class(DGModel1D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,iEl,iVar

    !$omp target map(tofrom: this % solution % interior, this % workSol % interior) map(to:this % dsdt % interior)
    !$omp teams distribute parallel do collapse(3) num_threads(256)
    do iEl = 1,this%solution%nElem
      do iVar = 1,this%solution%nVar
        do i = 1,this%solution%interp%N+1

          this%workSol%interior(i,iEl,iVar) = rk3_a(m)* &
                                              this%workSol%interior(i,iEl,iVar)+ &
                                              this%dSdt%interior(i,iEl,iVar)

          this%solution%interior(i,iEl,iVar) = &
            this%solution%interior(i,iEl,iVar)+ &
            rk3_g(m)*this%dt*this%workSol%interior(i,iEl,iVar)

        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine UpdateGRK3_DGModel1D

  subroutine UpdateGRK4_DGModel1D(this,m)
    implicit none
    class(DGModel1D),intent(inout) :: this
    integer,intent(in) :: m
    ! Local
    integer :: i,iEl,iVar

    !$omp target map(tofrom: this % solution % interior, this % workSol % interior) map(to:this % dsdt % interior)
    !$omp teams distribute parallel do collapse(3) num_threads(256)
    do iEl = 1,this%solution%nElem
      do iVar = 1,this%solution%nVar
        do i = 1,this%solution%interp%N+1

          this%workSol%interior(i,iEl,iVar) = rk4_a(m)* &
                                              this%workSol%interior(i,iEl,iVar)+ &
                                              this%dSdt%interior(i,iEl,iVar)

          this%solution%interior(i,iEl,iVar) = &
            this%solution%interior(i,iEl,iVar)+ &
            rk4_g(m)*this%dt*this%workSol%interior(i,iEl,iVar)

        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine UpdateGRK4_DGModel1D

  subroutine CalculateTendency_DGModel1D(this)
    implicit none
    class(DGModel1D),intent(inout) :: this
    ! Local
    integer :: i,iEl,iVar

    call this%solution%BoundaryInterp()
    call this%solution%SideExchange(this%mesh,this%decomp)
    call this%PreTendency()
    call this%SetBoundaryCondition()
    call this%SourceMethod()
    call this%RiemannSolver()
    call this%FluxMethod()
    call this%flux%DGDerivative(this%geometry,this%fluxDivergence%interior)

    !$omp target map(to: this % source, this % fluxDivergence) map(from:this % dSdt)
    !$omp teams distribute parallel do collapse(3) num_threads(256)
    do iEl = 1,this%solution%nElem
      do iVar = 1,this%solution%nVar
        do i = 1,this%solution%interp%N+1

          this%dSdt%interior(i,iEl,iVar) = &
            this%source%interior(i,iEl,iVar)- &
            this%fluxDivergence%interior(i,iEl,iVar)

        enddo
      enddo
    enddo
    !$omp end target

  endsubroutine CalculateTendency_DGModel1D

  subroutine Write_DGModel1D(this,fileName)
#undef __FUNC__
#define __FUNC__ "Write_DGModel1D"
    implicit none
    class(DGModel1D),intent(inout) :: this
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
      pickupFile = trim(filename)//timeStampString//'.h5'
    else
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

  endsubroutine Write_DGModel1D

  subroutine Read_DGModel1D(this,fileName)
    implicit none
    class(DGModel1D),intent(inout) :: this
    character(*),intent(in) :: fileName
    ! Local
    integer(HID_T) :: fileId
    integer(HID_T) :: solOffset(1:3)
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
      solOffset(1:3) = (/0,1,firstElem/)
      call ReadArray_HDF5(fileId,'/controlgrid/solution/interior', &
                          this%solution%interior,solOffset)
    else
      call ReadArray_HDF5(fileId,'/controlgrid/solution/interior',this%solution%interior)
    endif

    call Close_HDF5(fileId)

  endsubroutine Read_DGModel1D

  subroutine WriteTecplot_DGModel1D(this,filename)
    implicit none
    class(DGModel1D),intent(inout) :: this
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

        !WRITE (zoneID,'(I8.8)') iEl
        !WRITE (fUnit,*) 'ZONE T="el'//TRIM(zoneID)//'", I=',this % solution % interp % M + 1

        do i = 0,this%solution%interp%M

          write(fUnit,fmat) x%interior(i,1,iEl), &
            solution%interior(i,iEl,iVar)

        enddo

      enddo
    enddo

    close(UNIT=fUnit)

    call x%Free()
    call solution%Free()
    call interp%Free()

  endsubroutine WriteTecplot_DGModel1D

endmodule SELF_DGModel1D
