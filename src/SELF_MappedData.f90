! SELF_MappedData.F90
!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_MappedData

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Data
  use SELF_Mesh
  use SELF_Geometry
  use SELF_HDF5
  use HDF5

  use FEQParse

  use iso_c_binding

  implicit none

#include "SELF_Macros.h"

  type,extends(Scalar1D),public :: MappedScalar1D

  contains
    procedure,public :: SideExchange => SideExchange_MappedScalar1D
    procedure,public :: BassiRebaySides => BassiRebaySides_MappedScalar1D

    generic,public :: Derivative => Derivative_MappedScalar1D
    procedure,private :: Derivative_MappedScalar1D
    generic,public :: DGDerivative => DGDerivative_MappedScalar1D
    procedure,private :: DGDerivative_MappedScalar1D
    generic,public :: BRDerivative => BRDerivative_MappedScalar1D
    procedure,private :: BRDerivative_MappedScalar1D

    procedure,public :: JacobianWeight => JacobianWeight_MappedScalar1D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar1D

  end type MappedScalar1D

  type,extends(Scalar2D),public :: MappedScalar2D

    type(Tensor2D) :: JaScalar ! contravariant weighted scalar
  contains

    procedure,public :: Init => Init_MappedScalar2D
    procedure,public :: Free => Free_MappedScalar2D
    procedure,public :: SideExchange => SideExchange_MappedScalar2D
    procedure,public :: BassiRebaySides => BassiRebaySides_MappedScalar2D

    procedure,public :: ContravariantWeightInterior => ContravariantWeightInterior_MappedScalar2D
    procedure,public :: ContravariantWeightAvgBoundary => ContravariantWeightAvgBoundary_MappedScalar2D

    generic,public :: Gradient => Gradient_MappedScalar2D
    procedure,private :: Gradient_MappedScalar2D

    generic,public :: BRGradient => BRGradient_MappedScalar2D
    procedure,private :: BRGradient_MappedScalar2D

    procedure,public :: JacobianWeight => JacobianWeight_MappedScalar2D

    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedScalar2D
    procedure,private :: ApplyFlip => ApplyFlip_MappedScalar2D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar2D

    !PROCEDURE,PUBLIC :: Integral => Integral_MappedScalar2D

  end type MappedScalar2D

  type,extends(Scalar3D),public :: MappedScalar3D

    type(Tensor3D) :: JaScalar ! contravariant weighted scalar
  contains

    procedure,public :: Init => Init_MappedScalar3D
    procedure,public :: Free => Free_MappedScalar3D
    procedure,public :: SideExchange => SideExchange_MappedScalar3D
    procedure,public :: BassiRebaySides => BassiRebaySides_MappedScalar3D

    procedure,public :: ContravariantWeightInterior => ContravariantWeightInterior_MappedScalar3D
    procedure,public :: ContravariantWeightAvgBoundary => ContravariantWeightAvgBoundary_MappedScalar3D

    generic,public :: Gradient => Gradient_MappedScalar3D
    procedure,private :: Gradient_MappedScalar3D

    generic,public :: BRGradient => BRGradient_MappedScalar3D
    procedure,private :: BRGradient_MappedScalar3D

    procedure,public :: JacobianWeight => JacobianWeight_MappedScalar3D

    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedScalar3D
    procedure,private :: ApplyFlip => ApplyFlip_MappedScalar3D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedScalar3D

  end type MappedScalar3D

  type,extends(Vector2D),public :: MappedVector2D

  contains

    procedure,public :: SideExchange => SideExchange_MappedVector2D
    procedure,public :: BassiRebaySides => BassiRebaySides_MappedVector2D

    procedure,public :: ContravariantProjection => ContravariantProjection_MappedVector2D

    generic,public :: Divergence => Divergence_MappedVector2D
    procedure,private :: Divergence_MappedVector2D

    generic,public :: DGDivergence => DGDivergence_MappedVector2D
    procedure,private :: DGDivergence_MappedVector2D

    procedure,public :: JacobianWeight => JacobianWeight_MappedVector2D

    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedVector2D
    procedure,private :: ApplyFlip => ApplyFlip_MappedVector2D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedVector2D

  end type MappedVector2D

  type,extends(Vector3D),public :: MappedVector3D

  contains

    procedure,public :: SideExchange => SideExchange_MappedVector3D
    procedure,public :: BassiRebaySides => BassiRebaySides_MappedVector3D

    procedure,public :: ContravariantProjection => ContravariantProjection_MappedVector3D

    generic,public :: Divergence => Divergence_MappedVector3D
    procedure,private :: Divergence_MappedVector3D

    generic,public :: DGDivergence => DGDivergence_MappedVector3D
    procedure,private :: DGDivergence_MappedVector3D
    
    procedure,public :: JacobianWeight => JacobianWeight_MappedVector3D

    procedure,private :: MPIExchangeAsync => MPIExchangeAsync_MappedVector3D
    procedure,private :: ApplyFlip => ApplyFlip_MappedVector3D

    procedure,public :: SetInteriorFromEquation => SetInteriorFromEquation_MappedVector3D

  end type MappedVector3D

contains

! ---------------------- Scalars ---------------------- !

  subroutine SetInteriorFromEquation_MappedScalar1D(this,geometry,time)
    !!  Sets the this % interior attribute using the eqn attribute,
    !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar1D),intent(inout) :: this
    type(Geometry1D),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,iEl,iVar

    do ivar = 1,this % nvar
      this % interior(:,:,ivar) = this % eqn(ivar) % evaluate(geometry % x % interior)
    end do

  end subroutine SetInteriorFromEquation_MappedScalar1D

  subroutine SideExchange_MappedScalar1D(this,mesh,decomp)
    implicit none
    class(MappedScalar1D),intent(inout) :: this
    type(Mesh1D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,ivar
    integer :: neighborRank
    integer :: rankId,offset

    !$omp target map(to:decomp % offsetElem, this % boundary) map(tofrom:this % extBoundary)
    !$omp teams distribute parallel do collapse(2) num_threads(256)
    do ivar = 1,this % nvar
      do e1 = 1,mesh % nElem

        if (e1 == 1) then

          s1 = 2
          e2 = e1 + 1
          s2 = 1
          this % extBoundary(s1,e1,ivar) = this % boundary(s2,e2,ivar)

        elseif (e1 == mesh % nElem) then

          s1 = 1
          e2 = e1 - 1
          s2 = 2
          this % extBoundary(s1,e1,ivar) = this % boundary(s2,e2,ivar)

        else

          s1 = 1
          e2 = e1 - 1
          s2 = 2
          this % extBoundary(s1,e1,ivar) = this % boundary(s2,e2,ivar)

          s1 = 2
          e2 = e1 + 1
          s2 = 1
          this % extBoundary(s1,e1,ivar) = this % boundary(s2,e2,ivar)

        end if

      end do
    end do
    !$omp end target

  end subroutine SideExchange_MappedScalar1D

  subroutine BassiRebaySides_MappedScalar1D(this)
    implicit none
    class(MappedScalar1D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i

    !$omp target map(to:this % boundary, this % extBoundary) map(from:this % avgBoundary)
    !$omp teams distribute parallel do collapse(2)
    do iel = 1,this % nElem
      do ivar = 1,this % nVar

        ! Left side - we account for the -\hat{x} normal
        this % avgBoundary(1,iel,ivar) = -0.5_prec*( &
                                           this % boundary(1,iel,ivar) + &
                                           this % extBoundary(1,iel,ivar))

        ! Right side - we account for the +\hat{x} normal
        this % avgBoundary(2,iel,ivar) = 0.5_prec*( &
                                           this % boundary(2,iel,ivar) + &
                                           this % extBoundary(2,iel,ivar))
      end do
    end do
    !$omp end do

  end subroutine BassiRebaySides_MappedScalar1D

  subroutine Derivative_MappedScalar1D(this,geometry,dF)
    implicit none
    class(MappedScalar1D),intent(in) :: this
    type(Geometry1D),intent(in) :: geometry
    type(MappedScalar1D),intent(inout) :: dF

    call this % interp % Derivative_1D(this % interior, &
                                           df % interior, &
                                           this % nVar, &
                                           this % nElem)
    call df % JacobianWeight(geometry)

  end subroutine Derivative_MappedScalar1D

  subroutine DGDerivative_MappedScalar1D(this,geometry,dF)
    implicit none
    class(MappedScalar1D),intent(in) :: this
    type(Geometry1D),intent(in) :: geometry
    type(MappedScalar1D),intent(inout) :: dF


    call this % interp % DGDerivative_1D(this % interior, &
                                             this % boundary, &
                                             df % interior, &
                                             this % nVar, &
                                             this % nElem)

    call df % JacobianWeight(geometry)

  end subroutine DGDerivative_MappedScalar1D

  subroutine BRDerivative_MappedScalar1D(this,geometry,dF)
    implicit none
    class(MappedScalar1D),intent(in) :: this
    type(Geometry1D),intent(in) :: geometry
    type(MappedScalar1D),intent(inout) :: dF


    call this % interp % DGDerivative_1D(this % interior, &
                                             this % avgboundary, &
                                             df % interior, &
                                             this % nVar, &
                                             this % nElem)
    call df % JacobianWeight(geometry)

  end subroutine BRDerivative_MappedScalar1D

  subroutine JacobianWeight_MappedScalar1D(this,geometry)
#undef __FUNC__
#define __FUNC__ "JacobianWeight_MappedScalar1D"
    ! Applies the inverse jacobian
    implicit none
    class(MappedScalar1D),intent(inout) :: this
    type(Geometry1D),intent(in) :: geometry
    ! Local
    integer :: iEl,iVar,i

    !$omp target map(to:geometry % dxds % interior) map(tofrom:this % interior)
    !$omp teams distribute parallel do collapse(3)
      do iEl = 1,this % nElem
        do iVar = 1,this % nVar
          do i = 1,this % interp % N + 1
            this % interior(i,iEl,iVar) = this % interior(i,iEl,iVar)/ &
                                            geometry % dxds % interior(i,iEl,1)
          end do
        end do
      end do
      !$omp end target

  end subroutine JacobianWeight_MappedScalar1D

  subroutine Init_MappedScalar2D(this,interp,nVar,nElem)
    implicit none
    class(MappedScalar2D),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this % interp => interp
    this % nVar = nVar
    this % nElem = nElem

    allocate( this % interior(1:interp % N + 1,interp % N + 1,nelem,nvar),&
    this % interpWork(1:interp % M + 1,1:interp % N + 1,1:nelem,1:nvar),&
    this % boundary(1:interp % N + 1,1:4,1:nelem,1:nvar),&
    this % extBoundary(1:interp % N + 1,1:4,1:nelem,1:nvar),&
    this % avgBoundary(1:interp % N + 1,1:4,1:nelem,1:nvar),&
    this % jumpBoundary(1:interp % N + 1,1:4,1:nelem,1:nvar) )

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % interpWork)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)
    !$omp target enter data map(alloc: this % avgBoundary)
    !$omp target enter data map(alloc: this % jumpBoundary)

    allocate (this % meta(1:nVar))
    allocate (this % eqn(1:nVar))

    call this % JaScalar % Init(interp,nVar,nElem)

  end subroutine Init_MappedScalar2D

  subroutine Free_MappedScalar2D(this)
    implicit none
    class(MappedScalar2D),intent(inout) :: this

    this % nVar = 0
    this % nElem = 0
    this % interp => null()
    deallocate(this % interior)
    deallocate(this % interpWork)
    deallocate(this % boundary)
    deallocate(this % extBoundary)
    deallocate(this % avgBoundary)
    deallocate(this % jumpBoundary)
    deallocate (this % meta)
    deallocate (this % eqn)

    !$omp target exit data map(delete: this % interior)
    !$omp target exit data map(delete: this % interpWork)
    !$omp target exit data map(delete: this % boundary)
    !$omp target exit data map(delete: this % extBoundary)
    !$omp target exit data map(delete: this % avgBoundary)
    !$omp target exit data map(delete: this % jumpBoundary)
    call this % JaScalar % Free()

  end subroutine Free_MappedScalar2D

  subroutine SetInteriorFromEquation_MappedScalar2D(this,geometry,time)
  !!  Sets the this % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,j,iEl,iVar
    real(prec) :: x
    real(prec) :: y

    do iVar = 1,this % nVar
      do iEl = 1,this % nElem
        do j = 1,this % interp % N + 1
          do i = 1,this % interp % N + 1

            ! Get the mesh positions
            x = geometry % x % interior(i,j,iEl,1,1)
            y = geometry % x % interior(i,j,iEl,1,2)

            this % interior(i,j,iEl,iVar) = &
              this % eqn(iVar) % Evaluate((/x,y,0.0_prec,time/))

          end do
        end do
      end do
    end do

  end subroutine SetInteriorFromEquation_MappedScalar2D

  subroutine MPIExchangeAsync_MappedScalar2D(this,decomp,mesh,resetCount)
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh2D),intent(in) :: mesh
    logical,intent(in) :: resetCount
    ! Local
    integer :: e1,s1,e2,s2,ivar
    integer :: globalSideId,r2
    integer :: iError
    integer :: msgCount

    if (decomp % mpiEnabled) then
      if (resetCount) then
        msgCount = 0
      else
        msgCount = decomp % msgCount
      end if

      do ivar = 1,this % nvar
        do e1 = 1,this % nElem
          do s1 = 1,4

            e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
            if (e2 > 0) then
              r2 = decomp % elemToRank(e2) ! Neighbor Rank

              if (r2 /= decomp % rankId) then

                ! to do : create unique tag for each side and each variable
                ! tag = globalsideid + nglobalsides*ivar
                s2 = mesh % sideInfo(4,s1,e1)/10
                globalSideId = abs(mesh % sideInfo(2,s1,e1))

                msgCount = msgCount + 1
                call MPI_IRECV(this % extBoundary(:,s1,e1,ivar), &
                               (this % interp % N + 1), &
                               decomp % mpiPrec, &
                               r2,globalSideId, &
                               decomp % mpiComm, &
                               decomp % requests(msgCount),iError)

                msgCount = msgCount + 1
                call MPI_ISEND(this % boundary(:,s1,e1,ivar), &
                               (this % interp % N + 1), &
                               decomp % mpiPrec, &
                               r2,globalSideId, &
                               decomp % mpiComm, &
                               decomp % requests(msgCount),iError)
              end if
            end if

          end do
        end do
      end do

      decomp % msgCount = msgCount
    end if

  end subroutine MPIExchangeAsync_MappedScalar2D

  subroutine ApplyFlip_MappedScalar2D(this,decomp,mesh)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh2D),intent(in) :: mesh
    ! Local
    integer :: e1,s1,e2,s2
    integer :: i,i2
    integer :: r2,flip,ivar
    integer :: globalSideId
    integer :: bcid
    real(prec) :: extBuff(1:this % interp % N + 1)

    if (decomp % mpiEnabled) then
      !$omp target map(to:mesh % sideInfo, decomp % elemToRank) map(tofrom:this % extBoundary)
      !$omp teams distribute parallel do collapse(3) 
      do ivar = 1,this % nvar
        do e1 = 1,this % nElem
          do s1 = 1,4

            e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
            s2 = mesh % sideInfo(4,s1,e1)/10
            bcid = mesh % sideInfo(5,s1,e1)
            if (s2 > 0 .or. bcid == 0) then ! Interior Element
              r2 = decomp % elemToRank(e2) ! Neighbor Rank

              if (r2 /= decomp % rankId) then

                flip = mesh % sideInfo(4,s1,e1) - s2*10
                globalSideId = mesh % sideInfo(2,s1,e1)

                ! Need to update extBoundary with flip applied
                if (flip == 1) then

                  do i = 1,this % interp % N + 1
                    i2 = this % interp % N + 2 - i
                    extBuff(i) = this % extBoundary(i2,s1,e1,ivar)
                  end do
                  do i = 1,this % interp % N + 1
                    this % extBoundary(i,s1,e1,ivar) = extBuff(i)
                  end do

                end if
              end if

            end if

          end do
        end do
      end do
      !$omp end target

    end if

  end subroutine ApplyFlip_MappedScalar2D

  subroutine SideExchange_MappedScalar2D(this,mesh,decomp)
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(Mesh2D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,ivar
    integer :: neighborRank
    integer :: rankId,offset

    rankId = decomp % rankId
    offset = decomp % offsetElem(rankId + 1)

    call this % MPIExchangeAsync(decomp,mesh,resetCount=.true.)

    !$omp target map(to: mesh % sideInfo, decomp % elemToRank) map(from: this % boundary) map(tofrom: this % extBoundary)
    !$omp teams distribute parallel do collapse(3)
    do ivar = 1,this % nvar
      do e1 = 1,mesh % nElem
        do s1 = 1,4
          e2Global = mesh % sideInfo(3,s1,e1)
          e2 = e2Global - offset
          s2 = mesh % sideInfo(4,s1,e1)/10
          flip = mesh % sideInfo(4,s1,e1) - s2*10
          bcid = mesh % sideInfo(5,s1,e1)

          if (s2 > 0 .or. bcid == 0) then
            neighborRank = decomp % elemToRank(e2Global)

            if (neighborRank == decomp % rankId) then

              if (flip == 0) then

                do i1 = 1,this % interp % N + 1
                  this % extBoundary(i1,s1,e1,ivar) = &
                    this % boundary(i1,s2,e2,ivar)
                end do

              elseif (flip == 1) then

                do i1 = 1,this % interp % N + 1
                  i2 = this % interp % N + 2 - i1
                  this % extBoundary(i1,s1,e1,ivar) = &
                    this % boundary(i2,s2,e2,ivar)
                end do

              end if

            end if

          end if

        end do
      end do
    end do
    !$omp end target

    call decomp % FinalizeMPIExchangeAsync()

    ! Apply side flips for data exchanged with MPI
    call this % ApplyFlip(decomp,mesh)

  end subroutine SideExchange_MappedScalar2D

  subroutine ContravariantWeightInterior_MappedScalar2D(this,geometry)
    !! Computes the scalar multiplied by the contravariant basis vectors
    !!
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    ! Local
    integer    :: i,j,ii,iEl,iVar,row,col

    !$omp target map(to: geometry % dsdx % interior, this % interior) map(from: this % JaScaalar % interior)
    !$omp teams distribute parallel do collapse(6) num_threads(256)
    do col = 1,2
      do row = 1,2
        do iVar = 1,this % nVar
          do iEl = 1,this % nElem
            do j = 1,this % interp % N + 1
              do i = 1,this % interp % N + 1

                this % JaScalar % interior(i,j,iel,ivar,row,col) = geometry % dsdx % interior(i,j,iel,1,row,col)* &
                                                                      this % interior(i,j,iel,ivar)
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end target

  end subroutine ContravariantWeightInterior_MappedScalar2D

  subroutine ContravariantWeightAvgBoundary_MappedScalar2D(this,geometry)
    !! Computes the scalar multiplied by the contravariant basis vectors
    !!
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    ! Local
    integer    :: i,j,ii,iEl,iVar,row,col

    !$omp target map(to:geometry % dsdx % boundary, this % avgBoundary) map(from: this % JaScalar % boundary)
    !$omp teams distribute parallel do collapse(6) num_threads(256)
    do col = 1,2
      do row = 1,2
        do iVar = 1,this % nVar
          do iEl = 1,this % nElem
            do j = 1,4
              do i = 1,this % interp % N + 1

                this % JaScalar % boundary(i,j,iel,ivar,row,col) = geometry % dsdx % boundary(i,j,iel,1,row,col)* &
                                                                      this % avgBoundary(i,j,iel,ivar)
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end target

  end subroutine ContravariantWeightAvgBoundary_MappedScalar2D

  subroutine BassiRebaySides_MappedScalar2D(this)
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i

    !$omp target map(to:this % boundary, this % extBoundary) map(from:this % avgBoundary)
    !$omp teams distribute parallel do collapse(4) num_threads(256)
    do ivar = 1,this % nVar
      do iel = 1,this % nElem
        do iside = 1,4
          do i = 1,this % interp % N + 1
            this % avgBoundary(i,iside,iel,ivar) = 0.5_prec*( &
                                                      this % boundary(i,iside,iel,ivar) + &
                                                      this % extBoundary(i,iside,iel,ivar))
          end do
        end do
      end do
    end do
    !$omp end target

  end subroutine BassiRebaySides_MappedScalar2D

  subroutine JacobianWeight_MappedScalar2D(this,geometry)

    ! Applies the inverse jacobian
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    ! Local
    integer :: iEl,iVar,i,j

    !$omp target map(to: geometry % J % interior) map(tofrom: this % interior)
    !$omp teams distribute parallel do collapse(4) num_threads(256)
    do iVar = 1,this % nVar
      do iEl = 1,this % nElem
        do j = 1,this % interp % N + 1
          do i = 1,this % interp % N + 1
            this % interior(i,j,iEl,iVar) = this % interior(i,j,iEl,iVar)/ &
                                              geometry % J % interior(i,j,iEl,1)
          end do
        end do
      end do
    end do
    !$omp end target

  end subroutine JacobianWeight_MappedScalar2D

  subroutine Gradient_MappedScalar2D(this,geometry,df)
  !! Calculates the gradient of a function using the strong form of the gradient
  !! in mapped coordinates.
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    type(MappedVector2D),intent(inout) :: df

    call this % ContravariantWeightInterior(geometry)
    call this % JaScalar % Divergence(df)
    call df % JacobianWeight(geometry)

  end subroutine Gradient_MappedScalar2D

  subroutine BRGradient_MappedScalar2D(this,geometry,df)
    !! Calculates the gradient of a function using the weak form of the gradient
    !! and the average boundary state.
    !! This method will compute the average boundary state from the
    !! and  attributes of
    implicit none
    class(MappedScalar2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    type(MappedVector2D),intent(inout) :: df

    call this % BassiRebaySides()
    call this % ContravariantWeightInterior(geometry)
    call this % ContravariantWeightAvgBoundary(geometry)
    call this % JaScalar % DGDivergence(df)
    call df % JacobianWeight(geometry)

  end subroutine BRGradient_MappedScalar2D

  subroutine Init_MappedScalar3D(this,interp,nVar,nElem)
    implicit none
    class(MappedScalar3D),intent(out) :: this
    type(Lagrange),intent(in),target :: interp
    integer,intent(in) :: nVar
    integer,intent(in) :: nElem

    this % interp => interp
    this % nVar = nVar
    this % nElem = nElem

    allocate( this % interior(1:interp % N + 1,1:interp % N + 1,1:interp % N + 1,1:nelem,1:nvar),&
              this % interpWork1(1:interp % M + 1,1:interp % N + 1,1:interp % N + 1,1:nelem,1:nvar),&
              this % interpWork2(1:interp % M + 1,1:interp % M + 1,1:interp % N + 1,1:nelem,1:nvar),&
              this % boundary(1:interp % N + 1,1:interp % N + 1,1:6,1:nelem,1:nvar),&
              this % extBoundary(1:interp % N + 1,1:interp % N + 1,1:6,1:nelem,1:nvar),&
              this % avgBoundary(1:interp % N + 1,1:interp % N + 1,1:6,1:nelem,1:nvar),&
              this % jumpBoundary(1:interp % N + 1,1:interp % N + 1,1:6,1:nelem,1:nvar))

    !$omp target enter data map(alloc: this % interior)
    !$omp target enter data map(alloc: this % interpWork1)
    !$omp target enter data map(alloc: this % interpWork2)
    !$omp target enter data map(alloc: this % boundary)
    !$omp target enter data map(alloc: this % extBoundary)
    !$omp target enter data map(alloc: this % avgBoundary)
    !$omp target enter data map(alloc: this % jumpBoundary)

    allocate (this % meta(1:nVar))
    allocate (this % eqn(1:nVar))

    call this % JaScalar % Init(interp,nVar,nElem)

  end subroutine Init_MappedScalar3D

  subroutine Free_MappedScalar3D(this)
    implicit none
    class(MappedScalar3D),intent(inout) :: this

    this % nVar = 0
    this % nElem = 0
    this % interp => null()
    deallocate(this % interior)
    deallocate(this % interpWork1)
    deallocate(this % interpWork2)
    deallocate(this % boundary)
    deallocate(this % extBoundary)
    deallocate(this % avgBoundary)
    deallocate(this % jumpBoundary)
    deallocate (this % meta)
    deallocate (this % eqn)

    !$omp target exit data map(delete: this % interior)
    !$omp target exit data map(delete: this % interpWork1)
    !$omp target exit data map(delete: this % interpWork2)
    !$omp target exit data map(delete: this % boundary)
    !$omp target exit data map(delete: this % extBoundary)
    !$omp target exit data map(delete: this % avgBoundary)
    !$omp target exit data map(delete: this % jumpBoundary)

    call this % JaScalar % Free()

  end subroutine Free_MappedScalar3D

  subroutine SetInteriorFromEquation_MappedScalar3D(this,geometry,time)
  !!  Sets the this % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,j,k,iEl,iVar
    real(prec) :: x
    real(prec) :: y
    real(prec) :: z

    do iVar = 1,this % nVar
      do iEl = 1,this % nElem
        do k = 1,this % interp % N + 1
          do j = 1,this % interp % N + 1
            do i = 1,this % interp % N + 1

              ! Get the mesh positions
              x = geometry % x % interior(i,j,k,iEl,1,1)
              y = geometry % x % interior(i,j,k,iEl,1,2)
              z = geometry % x % interior(i,j,k,iEl,1,3)

              this % interior(i,j,k,iEl,iVar) = &
                this % eqn(iVar) % Evaluate((/x,y,z,time/))

            end do
          end do
        end do
      end do
    end do

  end subroutine SetInteriorFromEquation_MappedScalar3D

  subroutine MPIExchangeAsync_MappedScalar3D(this,decomp,mesh,resetCount)
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh3D),intent(in) :: mesh
    logical,intent(in) :: resetCount
    ! Local
    integer :: e1,s1,e2,s2,ivar
    integer :: globalSideId,r2
    integer :: iError
    integer :: msgCount

    if (decomp % mpiEnabled) then
      if (resetCount) then
        msgCount = 0
      else
        msgCount = decomp % msgCount
      end if

      do ivar = 1,this % nvar
        do e1 = 1,this % nElem
          do s1 = 1,6

            e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
            if (e2 > 0) then
              r2 = decomp % elemToRank(e2) ! Neighbor Rank

              if (r2 /= decomp % rankId) then

                ! to do : create unique tag for each side and each variable
                ! tag = globalsideid + nglobalsides*ivar
                s2 = mesh % sideInfo(4,s1,e1)/10
                globalSideId = abs(mesh % sideInfo(2,s1,e1))

                msgCount = msgCount + 1
                call MPI_IRECV(this % extBoundary(:,:,s1,e1,ivar), &
                               (this % interp % N + 1)*(this % interp % N + 1), &
                               decomp % mpiPrec, &
                               r2,globalSideId, &
                               decomp % mpiComm, &
                               decomp % requests(msgCount),iError)

                msgCount = msgCount + 1
                call MPI_ISEND(this % boundary(:,:,s1,e1,ivar), &
                               (this % interp % N + 1)*(this % interp % N + 1), &
                               decomp % mpiPrec, &
                               r2,globalSideId, &
                               decomp % mpiComm, &
                               decomp % requests(msgCount),iError)
              end if
            end if

          end do
        end do
      end do

      decomp % msgCount = msgCount
    end if

  end subroutine MPIExchangeAsync_MappedScalar3D

  subroutine ApplyFlip_MappedScalar3D(this,decomp,mesh)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh3D),intent(in) :: mesh
    ! Local
    integer :: e1,s1,e2,s2
    integer :: i,i2,j,j2
    integer :: r2,flip,ivar
    integer :: globalSideId
    integer :: bcid
    real(prec) :: extBuff(1:this % interp % N + 1,1:this % interp % N + 1)

    if (decomp % mpiEnabled) then
      !$omp target map(to:mesh % sideInfo, decomp % elemToRank) map(tofrom:this % extBoundary)
      !$omp teams distribute parallel do collapse(3) 
      do ivar = 1,this % nvar
        do e1 = 1,this % nElem
          do s1 = 1,6

            e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
            s2 = mesh % sideInfo(4,s1,e1)/10
            bcid = mesh % sideInfo(5,s1,e1)
            if (s2 > 0 .or. bcid == 0) then ! Interior Element
              r2 = decomp % elemToRank(e2) ! Neighbor Rank

              if (r2 /= decomp % rankId) then

                flip = mesh % sideInfo(4,s1,e1) - s2*10
                globalSideId = mesh % sideInfo(2,s1,e1)

                ! Need to update extBoundary with flip applied
                if (flip == 1) then

                  DO j = 1,this % interp % N+1
                    DO i = 1,this % interp % N+1
                      i2 = j
                      j2 = this % interp % N+2-i
                      extBuff(i,j) = this % extBoundary(i2,j2,s1,e1,ivar)
                    END DO
                  END DO


                else if (flip == 2) then

                  DO j = 1,this % interp % N+1
                    DO i = 1,this % interp % N+1
                      i2 = this % interp % N + 2 - i
                      j2 = this % interp % N + 2 - j
                      extBuff(i,j) = this % extBoundary(i2,j2,s1,e1,ivar)
                    END DO
                  END DO

                else if (flip == 3) then

                  DO j = 1,this % interp % N+1
                    DO i = 1,this % interp % N+1
                      i2 = this % interp % N + 2 - j
                      j2 = i
                      extBuff(i,j) = this % extBoundary(i2,j2,s1,e1,ivar)
                    END DO
                  END DO

                else if (flip == 4) then

                  DO j = 1,this % interp % N+1
                    DO i = 1,this % interp % N+1
                      i2 = j
                      j2 = i
                      extBuff(i,j) = this % extBoundary(i2,j2,s1,e1,ivar)
                    END DO
                  END DO

                end if

                DO j = 1,this % interp % N + 1
                  DO i = 1,this % interp % N + 1
                    this % extBoundary(i,j,s1,e1,ivar) = extBuff(i,j)
                  END DO
                END DO

              end if

            end if

          end do
        end do
      end do
      !$omp end target

    end if

  end subroutine ApplyFlip_MappedScalar3D

  subroutine SideExchange_MappedScalar3D(this,mesh,decomp)
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(Mesh3D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,j1,j2,ivar
    integer :: neighborRank
    integer :: rankId,offset

    rankId = decomp % rankId
    offset = decomp % offsetElem(rankId + 1)

    call this % MPIExchangeAsync(decomp,mesh,resetCount=.true.)

    !$omp target map(to: mesh % sideInfo, decomp % elemToRank) map(from: this % boundary) map(tofrom: this % extBoundary)
    !$omp teams distribute parallel do collapse(3)
    do ivar = 1,this % nvar
      do e1 = 1,mesh % nElem
        do s1 = 1,6
          e2Global = mesh % sideInfo(3,s1,e1)
          e2 = e2Global - offset
          s2 = mesh % sideInfo(4,s1,e1)/10
          flip = mesh % sideInfo(4,s1,e1) - s2*10
          bcid = mesh % sideInfo(5,s1,e1)

          ! If either s2 or e2 are equal to zero, then this is an exterior boundary and
          ! the extBoundary attribute is assigned by a boundary condition
          if (s2 > 0 .or. bcid == 0) then

            neighborRank = decomp % elemToRank(e2Global)

            if (neighborRank == decomp % rankId) then

              if (flip == 0)then

                DO j1 = 1,this % interp % N+1
                  DO i1 = 1,this % interp % N+1
                    this % extBoundary(i1,j1,s1,e1,ivar) = &
                      this % boundary(i1,j1,s2,e2,ivar)
                  END DO
                END DO

              else if (flip == 1)then

                DO j1 = 1,this % interp % N+1
                  DO i1 = 1,this % interp % N+1

                    i2 = j1
                    j2 = this % interp % N + 2 - i1
                    this % extBoundary(i1,j1,s1,e1,ivar) = &
                      this % boundary(i2,j2,s2,e2,ivar)

                  END DO
                END DO

              else if (flip == 2)then

                DO j1 = 1,this % interp % N+1
                  DO i1 = 1,this % interp % N+1
                    i2 = this % interp % N + 2 - i1
                    j2 = this % interp % N + 2 - j1
                    this % extBoundary(i1,j1,s1,e1,ivar) = &
                      this % boundary(i2,j2,s2,e2,ivar)
                  END DO
                END DO

              else if (flip == 3)then

                DO j1 = 1,this % interp % N+1
                  DO i1 = 1,this % interp % N+1
                    i2 = this % interp % N + 2 - j1
                    j2 = i1
                    this % extBoundary(i1,j1,s1,e1,ivar) = &
                      this % boundary(i2,j2,s2,e2,ivar)
                  END DO
                END DO

              else if (flip == 4)then

                DO j1 = 1,this % interp % N+1
                  DO i1 = 1,this % interp % N+1
                    i2 = j1
                    j2 = i1
                    this % extBoundary(i1,j1,s1,e1,ivar) = &
                      this % boundary(i2,j2,s2,e2,ivar)
                  END DO
                END DO

              end if

            end if

          end if

        end do
      end do
    end do
    !$omp end target

    call decomp % FinalizeMPIExchangeAsync()

    ! Apply side flips for data exchanged with MPI
    call this % ApplyFlip(decomp,mesh)

  end subroutine SideExchange_MappedScalar3D

  subroutine ContravariantWeightInterior_MappedScalar3D(this,geometry)
    !! Computes the scalar multiplied by the contravariant basis vectors
    !!
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    ! Local
    integer    :: i,j,k,iEl,iVar,row,col


    ! Interior
    !$omp target map(to: geometry % dsdx % interior, this % interior) map(from: this % JaScaalar % interior)
    !$omp teams distribute parallel do collapse(7) num_threads(256)
    do col = 1,3
      do row = 1,3
        do iVar = 1,this % nVar
          do iEl = 1,this % nElem
            do k = 1,this % interp % N + 1
              do j = 1,this % interp % N + 1
                do i = 1,this % interp % N + 1

                  this % JaScalar % interior(i,j,k,iel,ivar,row,col) = geometry % dsdx % interior(i,j,k,iel,1,row,col)* &
                                                                      this % interior(i,j,k,iel,ivar)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end target

  end subroutine ContravariantWeightInterior_MappedScalar3D

  subroutine ContravariantWeightAvgBoundary_MappedScalar3D(this,geometry)
    !! Computes the scalar multiplied by the contravariant basis vectors
    !!
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    ! Local
    integer    :: i,j,k,iEl,iVar,row,col

    ! Interior
    !$omp target map(to:geometry % dsdx % boundary, this % avgBoundary) map(from: this % JaScalar % boundary)
    !$omp teams distribute parallel do collapse(7) num_threads(256)
    do col = 1,3
      do row = 1,3
        do iVar = 1,this % nVar
          do iEl = 1,this % nElem
            do k = 1,6
              do j = 1,this % interp % N + 1
                do i = 1,this % interp % N + 1

                  this % JaScalar % boundary(i,j,k,iel,ivar,row,col) = geometry % dsdx % boundary(i,j,k,iel,1,row,col)* &
                                                                          this % avgBoundary(i,j,k,iel,ivar)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end target
    
  end subroutine ContravariantWeightAvgBoundary_MappedScalar3D

  subroutine BassiRebaySides_MappedScalar3D(this)
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i, j

    !$omp target map(to:this % boundary, this % extBoundary) map(from:this % avgBoundary)
    !$omp teams distribute parallel do collapse(5) num_threads(256)
    do ivar = 1,this % nVar
      do iel = 1,this % nElem
        do iside = 1,6
          do j = 1,this % interp % N + 1
            do i = 1,this % interp % N + 1
              this % avgBoundary(i,j,iside,iel,ivar) = 0.5_prec*( &
                                                      this % boundary(i,j,iside,iel,ivar) + &
                                                      this % extBoundary(i,j,iside,iel,ivar))
            end do
          end do
        end do
      end do
    end do
    !$omp end target

  end subroutine BassiRebaySides_MappedScalar3D

  subroutine JacobianWeight_MappedScalar3D(this,geometry)

    ! Applies the inverse jacobian
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    ! Local
    integer :: iEl,iVar,i,j,k

    !$omp target map(to: geometry % J % interior) map(tofrom: this % interior)
    !$omp teams distribute parallel do collapse(5) num_threads(256)
    do iEl = 1,this % nElem
      do iVar = 1,this % nVar
        do k = 1,this % interp % N + 1
          do j = 1,this % interp % N + 1
            do i = 1,this % interp % N + 1
              this % interior(i,j,k,iEl,iVar) = this % interior(i,j,k,iEl,iVar)/ &
                                                  geometry % J % interior(i,j,k,iEl,1)
            end do
          end do
        end do
      end do
    end do
    !$omp end target

  end subroutine JacobianWeight_MappedScalar3D

  subroutine Gradient_MappedScalar3D(this,geometry,df)
  !! Calculates the gradient of a function using the strong form of the gradient
  !! in mapped coordinates.
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    type(MappedVector3D),intent(inout) :: df

    call this % ContravariantWeightInterior(geometry)
    call this % JaScalar % Divergence(df)
    call df % JacobianWeight(geometry)

  end subroutine Gradient_MappedScalar3D

  subroutine BRGradient_MappedScalar3D(this,geometry,df)
    !! Calculates the gradient of a function using the weak form of the gradient
    !! and the average boundary state.
    !! This method will compute the average boundary state from the
    !! and  attributes of
    implicit none
    class(MappedScalar3D),intent(inout) :: this
    type(SEMHex),intent(in) :: geometry
    type(MappedVector3D),intent(inout) :: df

    call this % BassiRebaySides()
    call this % ContravariantWeightInterior(geometry)
    call this % ContravariantWeightAvgBoundary(geometry)
    call this % JaScalar % DGDivergence(df)
    call df % JacobianWeight(geometry)

  end subroutine BRGradient_MappedScalar3D

  subroutine SetInteriorFromEquation_MappedVector2D(this,geometry,time)
  !!  Sets the this % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,j,iEl,iVar
    real(prec) :: x
    real(prec) :: y

    do iVar = 1,this % nVar
      do iEl = 1,this % nElem
        do j = 1,this % interp % N + 1
          do i = 1,this % interp % N + 1

            ! Get the mesh positions
            x = geometry % x % interior(i,j,iEl,1,1)
            y = geometry % x % interior(i,j,iEl,1,2)

            this % interior(i,j,iEl,iVar,1) = &
              this % eqn(1 + 2*(iVar - 1)) % Evaluate((/x,y,0.0_prec,time/))

            this % interior(i,j,iEl,iVar,2) = &
              this % eqn(2 + 2*(iVar - 1)) % Evaluate((/x,y,0.0_prec,time/))

          end do
        end do
      end do
    end do

  end subroutine SetInteriorFromEquation_MappedVector2D

  subroutine MPIExchangeAsync_MappedVector2D(this,decomp,mesh,resetCount)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh2D),intent(in) :: mesh
    logical,intent(in) :: resetCount
    ! Local
    integer :: e1,s1,e2,s2,ivar,idir
    integer :: globalSideId,r2
    integer :: iError
    integer :: msgCount

    if (decomp % mpiEnabled) then
      if (resetCount) then
        msgCount = 0
      else
        msgCount = decomp % msgCount
      end if

      do idir = 1,2
        do ivar = 1,this % nvar
          do e1 = 1,this % nElem
            do s1 = 1,4

              e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
              if (e2 > 0) then
                r2 = decomp % elemToRank(e2) ! Neighbor Rank

                if (r2 /= decomp % rankId) then

                  ! to do : create unique tag for each side and each variable
                  ! tag = globalsideid + nglobalsides*(ivar + nvar*idir)
                  s2 = mesh % sideInfo(4,s1,e1)/10
                  globalSideId = abs(mesh % sideInfo(2,s1,e1))

                  msgCount = msgCount + 1
                  call MPI_IRECV(this % extBoundary(:,s1,e1,ivar,idir), &
                                 (this % interp % N + 1), &
                                 decomp % mpiPrec, &
                                 r2,globalSideId, &
                                 decomp % mpiComm, &
                                 decomp % requests(msgCount),iError)

                  msgCount = msgCount + 1
                  call MPI_ISEND(this % boundary(:,s1,e1,ivar,idir), &
                                 (this % interp % N + 1), &
                                 decomp % mpiPrec, &
                                 r2,globalSideId, &
                                 decomp % mpiComm, &
                                 decomp % requests(msgCount),iError)
                end if
              end if

            end do
          end do
        end do
      end do

      decomp % msgCount = msgCount
    end if

  end subroutine MPIExchangeAsync_MappedVector2D

  subroutine ApplyFlip_MappedVector2D(this,decomp,mesh)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(MPILayer),intent(inout) :: decomp
    type(Mesh2D),intent(in) :: mesh
    ! Local
    integer :: e1,s1,e2,s2
    integer :: i,i2
    integer :: r2,flip,ivar,idir
    integer :: globalSideId
    integer :: bcid
    real(prec) :: extBuff(1:this % interp % N + 1)

    if (decomp % mpiEnabled) then
      !$omp target map(to:mesh % sideInfo, decomp % elemToRank) map(tofrom:this % extBoundary)
      !$omp teams distribute parallel do collapse(4) 
      do idir = 1,2
        do ivar = 1,this % nvar
          do e1 = 1,this % nElem
            do s1 = 1,4

              e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
              s2 = mesh % sideInfo(4,s1,e1)/10
              bcid = mesh % sideInfo(5,s1,e1)
              if (s2 > 0 .or. bcid == 0) then ! Interior Element
                r2 = decomp % elemToRank(e2) ! Neighbor Rank

                if (r2 /= decomp % rankId) then

                  flip = mesh % sideInfo(4,s1,e1) - s2*10
                  globalSideId = mesh % sideInfo(2,s1,e1)

                  ! Need to update extBoundary with flip applied
                  if (flip == 1) then

                    do i = 1,this % interp % N + 1
                      i2 = this % interp % N + 2 - i
                      extBuff(i) = this % extBoundary(i2,s1,e1,ivar,idir)
                    end do
                    do i = 1,this % interp % N + 1
                      this % extBoundary(i,s1,e1,ivar,idir) = extBuff(i)
                    end do

                  end if
                end if

              end if

            end do
          end do
        end do
      end do
      !$omp end target

    end if

  end subroutine ApplyFlip_MappedVector2D

  subroutine SideExchange_MappedVector2D(this,mesh,decomp)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(Mesh2D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,ivar,idir
    integer :: neighborRank
    integer :: rankId,offset

    rankId = decomp % rankId
    offset = decomp % offsetElem(rankId + 1)

    call this % MPIExchangeAsync(decomp,mesh,resetCount=.true.)
    !$omp target map(to: mesh % sideInfo, decomp % elemToRank) map(from: this % boundary) map(tofrom: this % extBoundary)
    !$omp teams distribute parallel do collapse(4)
    do idir = 1,2
      do ivar = 1,this % nvar
        do e1 = 1,mesh % nElem
          do s1 = 1,4
            e2Global = mesh % sideInfo(3,s1,e1)
            e2 = e2Global - offset
            s2 = mesh % sideInfo(4,s1,e1)/10
            flip = mesh % sideInfo(4,s1,e1) - s2*10
            bcid = mesh % sideInfo(5,s1,e1)

            if (s2 > 0 .or. bcid == 0) then

              neighborRank = decomp % elemToRank(e2Global)

              if (neighborRank == decomp % rankId) then

                if (flip == 0) then

                  do i1 = 1,this % interp % N + 1
                    this % extBoundary(i1,s1,e1,ivar,idir) = &
                      this % boundary(i1,s2,e2,ivar,idir)
                  end do

                elseif (flip == 1) then

                  do i1 = 1,this % interp % N + 1
                    i2 = this % interp % N + 2 - i1
                    this % extBoundary(i1,s1,e1,ivar,idir) = &
                      this % boundary(i2,s2,e2,ivar,idir)
                  end do

                end if

              end if

            end if

          end do
        end do
      end do
    end do
    !$omp end target

    call decomp % FinalizeMPIExchangeAsync()

    ! Apply side flips for data exchanged with MPI
    call this % ApplyFlip(decomp,mesh)

  end subroutine SideExchange_MappedVector2D

  subroutine BassiRebaySides_MappedVector2D(this)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i
    integer :: idir

    !$omp target map(to:this % boundary, this % extBoundary) map(from:this % avgBoundary)
    !$omp teams distribute parallel do collapse(5) num_threads(256)
    do idir = 1,2
      do ivar = 1,this % nVar
        do iel = 1,this % nElem
          do iside = 1,4
            do i = 1,this % interp % N + 1
              this % avgBoundary(i,iside,iel,ivar,idir) = 0.5_prec*( &
                                                          this % boundary(i,iside,iel,ivar,idir) + &
                                                          this % extBoundary(i,iside,iel,ivar,idir))
            end do
          end do
        end do
      end do
    end do
    !$omp end target

  end subroutine BassiRebaySides_MappedVector2D

  subroutine JacobianWeight_MappedVector2D(this,geometry)

    ! Applies the inverse jacobian
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    ! Local
    integer :: iEl,iVar,i,j,idir

    !$omp target map(to: geometry % J % interior) map(tofrom: this % interior)
    !$omp teams distribute parallel do collapse(5) num_threads(256)
    do idir = 1,2
      do iVar = 1,this % nVar
        do iEl = 1,this % nElem
          do j = 1,this % interp % N + 1
            do i = 1,this % interp % N + 1
              this % interior(i,j,iEl,iVar,idir) = this % interior(i,j,iEl,iVar,idir)/ &
                                                geometry % J % interior(i,j,iEl,1)
              
            end do
          end do
        end do
      end do
    end do
    !$omp end target

  end subroutine JacobianWeight_MappedVector2D

  subroutine ContravariantProjection_MappedVector2D(this,geometry)
    ! Takes a vector that has physical space coordinate directions (x,y,z) and projects the vector
    ! into the the contravariant basis vector directions. Keep in mind that the contravariant basis
    ! vectors are really the Jacobian weighted contravariant basis vectors
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    ! Local
    integer :: i,j,iEl,iVar
    real(prec) :: Fx,Fy

    ! Assume that tensor(j,i) is vector i, component j
    ! => dot product is done along first dimension
    ! to project onto computational space
    !$omp target map(to:geometry % dsdx % interior) map(tofrom:this % interior)
    !$omp teams distribute parallel do collapse(4) num_threads(256)
    do iel = 1,this % nElem
      do ivar = 1,this % nVar
        do j = 1,this % interp % N + 1
          do i = 1,this % interp % N + 1

            Fx = this % interior(i,j,iEl,iVar,1)
            Fy = this % interior(i,j,iEl,iVar,2)

            this % interior(i,j,iEl,iVar,1) = &
              geometry % dsdx % interior(i,j,iEl,1,1,1)*Fx + &
              geometry % dsdx % interior(i,j,iEl,1,2,1)*Fy

            this % interior(i,j,iEl,iVar,2) = &
              geometry % dsdx % interior(i,j,iEl,1,1,2)*Fx + &
              geometry % dsdx % interior(i,j,iEl,1,2,2)*Fy

          end do
        end do
      end do
    end do
    !$omp end target

  end subroutine ContravariantProjection_MappedVector2D

  subroutine Divergence_MappedVector2D(this,geometry,divVector)
    ! Strong Form Operator
    !    !
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    type(MappedScalar2D),intent(inout) :: divVector

    ! Convert from physical to computational space
    call this % ContravariantProjection(geometry)

    ! Compute the divergence
    call this % interp % VectorDivergence_2D(this % interior, &
                                              divVector % interior, &
                                              this % nvar, &
                                              this % nelem)
    ! Divide by the jacobian
    call divVector % JacobianWeight(geometry)

  end subroutine Divergence_MappedVector2D

  subroutine DGDivergence_MappedVector2D(this,geometry,divVector)
    !! Computes the divergence of a 2-D vector using the weak form
    !! On input, the  attribute of the vector
    !! is assigned and the  attribute is set to the physical
    !! directions of the vector. This method will project the vector
    !! onto the contravariant basis vectors.
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    type(MappedScalar2D),intent(inout) :: divVector

      ! Convert from physical to computational space
      call this % ContravariantProjection(geometry)

      ! Compute the divergence
      call this % interp % VectorDGDivergence_2D(this % interior, &
                                                 this % boundaryNormal, &
                                                 divVector % interior, &
                                                 this % nvar, &
                                                 this % nelem)
      ! Divide by the jacobian
      call divVector % JacobianWeight(geometry)

  end subroutine DGDivergence_MappedVector2D

  subroutine SetInteriorFromEquation_MappedVector3D(this,geometry,time)
    !!  Sets the this % interior attribute using the eqn attribute,
    !!  geometry (for physical positions), and provided simulation time.
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(SEMHex),intent(in) :: geometry
      real(prec),intent(in) :: time
      ! Local
      integer :: i,j,k,iEl,iVar
      real(prec) :: x
      real(prec) :: y
      real(prec) :: z
  
      do iVar = 1,this % nVar
        do iEl = 1,this % nElem
          do k = 1,this % interp % N + 1
            do j = 1,this % interp % N + 1
              do i = 1,this % interp % N + 1
    
                ! Get the mesh positions
                x = geometry % x % interior(i,j,k,iEl,1,1)
                y = geometry % x % interior(i,j,k,iEl,1,2)
                z = geometry % x % interior(i,j,k,iEl,1,3)
    
                this % interior(i,j,k,iEl,iVar,1) = &
                  this % eqn(1 + 3*(iVar - 1)) % Evaluate((/x,y,z,time/))
    
                this % interior(i,j,k,iEl,iVar,2) = &
                  this % eqn(2 + 3*(iVar - 1)) % Evaluate((/x,y,z,time/))

                this % interior(i,j,k,iEl,iVar,3) = &
                  this % eqn(3 + 3*(iVar - 1)) % Evaluate((/x,y,z,time/))
    
              end do
            end do
          end do
        end do
      end do
  
    end subroutine SetInteriorFromEquation_MappedVector3D

    subroutine MPIExchangeAsync_MappedVector3D(this,decomp,mesh,resetCount)
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(MPILayer),intent(inout) :: decomp
      type(Mesh3D),intent(in) :: mesh
      logical,intent(in) :: resetCount
      ! Local
      integer :: e1,s1,e2,s2,ivar,idir
      integer :: globalSideId,r2
      integer :: iError
      integer :: msgCount
  
      if (decomp % mpiEnabled) then
        if (resetCount) then
          msgCount = 0
        else
          msgCount = decomp % msgCount
        end if
  
        do idir = 1,3
          do ivar = 1,this % nvar
            do e1 = 1,this % nElem
              do s1 = 1,6
    
                e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
                if (e2 > 0) then
                  r2 = decomp % elemToRank(e2) ! Neighbor Rank
    
                  if (r2 /= decomp % rankId) then
    
                    ! to do : create unique tag for each side and each variable
                    ! tag = globalsideid + nglobalsides*ivar
                    s2 = mesh % sideInfo(4,s1,e1)/10
                    globalSideId = abs(mesh % sideInfo(2,s1,e1))
    
                    msgCount = msgCount + 1
                    call MPI_IRECV(this % extBoundary(:,:,s1,e1,ivar,idir), &
                                  (this % interp % N + 1)*(this % interp % N + 1), &
                                  decomp % mpiPrec, &
                                  r2,globalSideId, &
                                  decomp % mpiComm, &
                                  decomp % requests(msgCount),iError)
    
                    msgCount = msgCount + 1
                    call MPI_ISEND(this % boundary(:,:,s1,e1,ivar,idir), &
                                  (this % interp % N + 1)*(this % interp % N + 1), &
                                  decomp % mpiPrec, &
                                  r2,globalSideId, &
                                  decomp % mpiComm, &
                                  decomp % requests(msgCount),iError)
                  end if
                end if
    
              end do
            end do
          end do
        end do
  
        decomp % msgCount = msgCount
      end if
  
    end subroutine MPIExchangeAsync_MappedVector3D
  
    subroutine ApplyFlip_MappedVector3D(this,decomp,mesh)
      ! Apply side flips to sides where MPI exchanges took place.
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(MPILayer),intent(inout) :: decomp
      type(Mesh3D),intent(in) :: mesh
      ! Local
      integer :: e1,s1,e2,s2
      integer :: i,j,i2,j2
      integer :: r2,flip,ivar
      integer :: globalSideId
      integer :: bcid, idir
      real(prec) :: extBuff(1:this % interp % N + 1,1:this % interp % N + 1)
  
      if (decomp % mpiEnabled) then

        !$omp target map(to:mesh % sideInfo, decomp % elemToRank) map(tofrom:this % extBoundary)
        !$omp teams distribute parallel do collapse(4) 
        do idir = 1,3
          do ivar = 1,this % nvar
            do e1 = 1,this % nElem
              do s1 = 1,6
  
                e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
                s2 = mesh % sideInfo(4,s1,e1)/10
                bcid = mesh % sideInfo(5,s1,e1)
                if (s2 > 0 .or. bcid == 0) then ! Interior Element
                  r2 = decomp % elemToRank(e2) ! Neighbor Rank
  
                  if (r2 /= decomp % rankId) then
  
                    flip = mesh % sideInfo(4,s1,e1) - s2*10
                    globalSideId = mesh % sideInfo(2,s1,e1)
  
                    ! Need to update extBoundary with flip applied
                    if (flip == 1) then
  
                      DO j = 1,this % interp % N+1
                        DO i = 1,this % interp % N+1
                          i2 = j
                          j2 = this % interp % N+2-i
                          extBuff(i,j) = this % extBoundary(i2,j2,s1,e1,ivar,idir)
                        END DO
                      END DO
  
  
                    else if (flip == 2) then
  
                      DO j = 1,this % interp % N+1
                        DO i = 1,this % interp % N+1
                          i2 = this % interp % N + 2 - i
                          j2 = this % interp % N + 2 - j
                          extBuff(i,j) = this % extBoundary(i2,j2,s1,e1,ivar,idir)
                        END DO
                      END DO
  
                    else if (flip == 3) then
  
                      DO j = 1,this % interp % N+1
                        DO i = 1,this % interp % N+1
                          i2 = this % interp % N + 2 - j
                          j2 = i
                          extBuff(i,j) = this % extBoundary(i2,j2,s1,e1,ivar,idir)
                        END DO
                      END DO
  
                    else if (flip == 4) then
  
                      DO j = 1,this % interp % N+1
                        DO i = 1,this % interp % N+1
                          i2 = j
                          j2 = i
                          extBuff(i,j) = this % extBoundary(i2,j2,s1,e1,ivar,idir)
                        END DO
                      END DO
  
                    end if
  
                    DO j = 1,this % interp % N + 1
                      DO i = 1,this % interp % N + 1
                        this % extBoundary(i,j,s1,e1,ivar,idir) = extBuff(i,j)
                      END DO
                    END DO
  
                  end if
  
                end if
  
              end do
            end do
          end do
        end do
        !$omp end target

      end if
  
    end subroutine ApplyFlip_MappedVector3D
  
    subroutine SideExchange_MappedVector3D(this,mesh,decomp)
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(Mesh3D),intent(in) :: mesh
      type(MPILayer),intent(inout) :: decomp
      ! Local
      integer :: e1,e2,s1,s2,e2Global
      integer :: flip,bcid
      integer :: i1,i2,j1,j2,ivar
      integer :: neighborRank
      integer :: rankId,offset
      integer :: idir
  
      rankId = decomp % rankId
      offset = decomp % offsetElem(rankId + 1)
  
      call this % MPIExchangeAsync(decomp,mesh,resetCount=.true.)

      !$omp target map(to: mesh % sideInfo, decomp % elemToRank) map(from: this % boundary) map(tofrom: this % extBoundary)
      !$omp teams distribute parallel do collapse(4)
      do idir = 1, 3
        do ivar = 1,this % nvar
          do e1 = 1,mesh % nElem
            do s1 = 1,6
              e2Global = mesh % sideInfo(3,s1,e1)
              e2 = e2Global - offset
              s2 = mesh % sideInfo(4,s1,e1)/10
              flip = mesh % sideInfo(4,s1,e1) - s2*10
              bcid = mesh % sideInfo(5,s1,e1)
  
              if (s2 > 0 .or. bcid == 0) then
  
                neighborRank = decomp % elemToRank(e2Global)
  
                if (neighborRank == decomp % rankId) then
  
                  if (flip == 0)then
  
                    DO j1 = 1,this % interp % N+1
                      DO i1 = 1,this % interp % N+1
                        this % extBoundary(i1,j1,s1,e1,ivar,idir) = &
                          this % boundary(i1,j1,s2,e2,ivar,idir)
                      END DO
                    END DO
  
                  else if (flip == 1)then
  
                    DO j1 = 1,this % interp % N+1
                      DO i1 = 1,this % interp % N+1
  
                        i2 = j1
                        j2 = this % interp % N + 2 - i1
                        this % extBoundary(i1,j1,s1,e1,ivar,idir) = &
                          this % boundary(i2,j2,s2,e2,ivar,idir)
  
                      END DO
                    END DO
  
                  else if (flip == 2)then
  
                    DO j1 = 1,this % interp % N+1
                      DO i1 = 1,this % interp % N+1
                        i2 = this % interp % N + 2 - i1
                        j2 = this % interp % N + 2 - j1
                        this % extBoundary(i1,j1,s1,e1,ivar,idir) = &
                          this % boundary(i2,j2,s2,e2,ivar,idir)
                      END DO
                    END DO
  
                  else if (flip == 3)then
  
                    DO j1 = 1,this % interp % N+1
                      DO i1 = 1,this % interp % N+1
                        i2 = this % interp % N + 2 - j1
                        j2 = i1
                        this % extBoundary(i1,j1,s1,e1,ivar,idir) = &
                          this % boundary(i2,j2,s2,e2,ivar,idir)
                      END DO
                    END DO
  
                  else if (flip == 4)then
  
                    DO j1 = 1,this % interp % N+1
                      DO i1 = 1,this % interp % N+1
                        i2 = j1
                        j2 = i1
                        this % extBoundary(i1,j1,s1,e1,ivar,idir) = &
                          this % boundary(i2,j2,s2,e2,ivar,idir)
                      END DO
                    END DO
  
                  end if
  
                end if
  
              end if
  
            end do
          end do
        end do
      end do
      !$omp end target

      call decomp % FinalizeMPIExchangeAsync()
  
      ! Apply side flips for data exchanged with MPI
      call this % ApplyFlip(decomp,mesh)
  
    end subroutine SideExchange_MappedVector3D

    subroutine BassiRebaySides_MappedVector3D(this)
      implicit none
      class(MappedVector3D),intent(inout) :: this
      ! Local
      integer :: iel
      integer :: iside
      integer :: ivar
      integer :: i, j
      integer :: idir
  
      !$omp target map(to:this % boundary, this % extBoundary) map(from:this % avgBoundary)
      !$omp teams distribute parallel do collapse(6) num_threads(256)
      do idir = 1,3
        do ivar = 1,this % nVar
          do iel = 1,this % nElem
            do iside = 1,6
              do j = 1,this % interp % N + 1
                do i = 1,this % interp % N + 1
                  this % avgBoundary(i,j,iside,iel,ivar,idir) = 0.5_prec*( &
                                                              this % boundary(i,j,iside,iel,ivar,idir) + &
                                                              this % extBoundary(i,j,iside,iel,ivar,idir))
                end do
              end do
            end do
          end do
        end do
      end do
      !$omp end target
    
    end subroutine BassiRebaySides_MappedVector3D

    subroutine JacobianWeight_MappedVector3D(this,geometry)

      ! Applies the inverse jacobian
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(SEMHex),intent(in) :: geometry
      ! Local
      integer :: iEl,iVar,i,j,k,idir
  
      !$omp target map(to: geometry % J % interior) map(tofrom: this % interior)
      !$omp teams distribute parallel do collapse(6) num_threads(256)
      do idir = 1,3
        do iVar = 1,this % nVar
          do iEl = 1,this % nElem
            do k = 1,this % interp % N + 1
              do j = 1,this % interp % N + 1
                do i = 1,this % interp % N + 1
                  this % interior(i,j,k,iEl,iVar,idir) = this % interior(i,j,k,iEl,iVar,idir)/ &
                                                    geometry % J % interior(i,j,k,iEl,1)
                end do
              end do
            end do
          end do
        end do
      end do
      !$omp end target

    end subroutine JacobianWeight_MappedVector3D

    subroutine ContravariantProjection_MappedVector3D(this,geometry)
      ! Takes a vector that has physical space coordinate directions (x,y,z) and projects the vector
      ! into the the contravariant basis vector directions. Keep in mind that the contravariant basis
      ! vectors are really the Jacobian weighted contravariant basis vectors
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(SEMHex),intent(in) :: geometry
      ! Local
      integer :: i,j,k,iEl,iVar
      real(prec) :: Fx,Fy,Fz
  
      ! Assume that tensor(j,i) is vector i, component j
      ! => dot product is done along first dimension
      ! to project onto computational space
      !$omp target map(to:geometry % dsdx % interior) map(tofrom:this % interior)
      !$omp teams distribute parallel do collapse(5) num_threads(256)
      do iel = 1,this % nElem
        do ivar = 1,this % nVar
          do k = 1,this % interp % N + 1
            do j = 1,this % interp % N + 1
              do i = 1,this % interp % N + 1
  
                Fx = this % interior(i,j,k,iEl,iVar,1)
                Fy = this % interior(i,j,k,iEl,iVar,2)
                Fz = this % interior(i,j,k,iEl,iVar,3)
  
                this % interior(i,j,k,iEl,iVar,1) = &
                  geometry % dsdx % interior(i,j,k,iEl,1,1,1)*Fx + &
                  geometry % dsdx % interior(i,j,k,iEl,1,2,1)*Fy + &
                  geometry % dsdx % interior(i,j,k,iEl,1,3,1)*Fz
  
                this % interior(i,j,k,iEl,iVar,2) = &
                  geometry % dsdx % interior(i,j,k,iEl,1,1,2)*Fx + &
                  geometry % dsdx % interior(i,j,k,iEl,1,2,2)*Fy + &
                  geometry % dsdx % interior(i,j,k,iEl,1,3,2)*Fz

                this % interior(i,j,k,iEl,iVar,3) = &
                  geometry % dsdx % interior(i,j,k,iEl,1,1,3)*Fx + &
                  geometry % dsdx % interior(i,j,k,iEl,1,2,3)*Fy + &
                  geometry % dsdx % interior(i,j,k,iEl,1,3,3)*Fz
  
              end do
            end do
          end do
        end do
      end do
      !$omp end target
    
    end subroutine ContravariantProjection_MappedVector3D
  
    subroutine Divergence_MappedVector3D(this,geometry,divVector)
      ! Strong Form Operator
      !    !
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(SEMHex),intent(in) :: geometry
      type(MappedScalar3D),intent(inout) :: divVector
  
      ! Convert from physical to computational space
      call this % ContravariantProjection(geometry)

      ! Compute the divergence
      call this % interp % VectorDivergence_3D(this % interior, &
                                                divVector % interior, &
                                                this % nvar, &
                                                this % nelem)
      ! Divide by the jacobian
      call divVector % JacobianWeight(geometry)
  
    end subroutine Divergence_MappedVector3D
  
    subroutine DGDivergence_MappedVector3D(this,geometry,divVector)
      !! Computes the divergence of a 3-D vector using the weak form
      !! On input, the  attribute of the vector
      !! is assigned and the  attribute is set to the physical
      !! directions of the vector. This method will project the vector
      !! onto the contravariant basis vectors.
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(SEMHex),intent(in) :: geometry
      type(MappedScalar3D),intent(inout) :: divVector
    
      ! Convert from physical to computational space
      call this % ContravariantProjection(geometry)

      ! Compute the divergence
      call this % interp % VectorDGDivergence_3D(this % interior, &
                                                  this % boundaryNormal, &
                                                  divVector % interior, &
                                                  this % nvar, &
                                                  this % nelem)
      ! Divide by the jacobian
      call divVector % JacobianWeight(geometry)
  
    end subroutine DGDivergence_MappedVector3D

end module SELF_MappedData
