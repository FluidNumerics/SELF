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

    !PROCEDURE,PUBLIC :: Integral => Integral_MappedScalar2D

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

  interface
    subroutine JacobianWeight_1D_gpu(scalar,dxds,N,nVar,nEl) &
      bind(c,name="JacobianWeight_1D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: scalar,dxds
      integer(c_int),value :: N,nVar,nEl
    end subroutine JacobianWeight_1D_gpu
  end interface

  interface
    subroutine JacobianWeight_2D_gpu(scalar,jacobian,N,nVar,nEl) &
      bind(c,name="JacobianWeight_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: scalar,jacobian
      integer(c_int),value :: N,nVar,nEl
    end subroutine JacobianWeight_2D_gpu
  end interface

  interface
    subroutine JacobianWeight_3D_gpu(scalar,jacobian,N,nVar,nEl) &
      bind(c,name="JacobianWeight_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: scalar,jacobian
      integer(c_int),value :: N,nVar,nEl
    end subroutine JacobianWeight_3D_gpu
  end interface

  interface
    subroutine ContravariantWeight_MappedScalar2D_gpu(scalar,dsdx,tensor,isize,jsize,nvar,nel) &
      bind(c,name="ContravariantWeight_MappedScalar2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: scalar,dsdx,tensor
      integer(c_int),value :: isize,jsize,nvar,nel
    end subroutine ContravariantWeight_MappedScalar2D_gpu
  end interface

  interface
    subroutine ContravariantWeight_MappedScalar3D_gpu(scalar,dsdx,tensor,isize,jsize,ksize,nvar,nel) &
      bind(c,name="ContravariantWeight_MappedScalar3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: scalar,dsdx,tensor
      integer(c_int),value :: isize,jsize,ksize,nvar,nel
    end subroutine ContravariantWeight_MappedScalar3D_gpu
  end interface

  interface
    subroutine ContravariantProjection_MappedVector2D_gpu_wrapper(vector,dsdx,N,nvar,nel) &
      bind(c,name="ContravariantProjection_MappedVector2D_gpu_wrapper")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: vector,dsdx
      integer(c_int),value :: N,nvar,nel
    end subroutine ContravariantProjection_MappedVector2D_gpu_wrapper
  end interface

  interface
    subroutine ContravariantProjection_MappedVector3D_gpu_wrapper(vector,dsdx,N,nvar,nel) &
      bind(c,name="ContravariantProjection_MappedVector3D_gpu_wrapper")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: vector,dsdx
      integer(c_int),value :: N,nvar,nel
    end subroutine ContravariantProjection_MappedVector3D_gpu_wrapper
  end interface

  interface
    subroutine SideExchange_2D_gpu(extBoundary,boundary,sideInfo,elemToRank,rankId,offset,N,nVar,nEl) &
      bind(c,name="SideExchange_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: extBoundary,boundary,sideInfo,elemToRank
      integer(c_int),value :: rankId,offset,N,nVar,nEl
    end subroutine SideExchange_2D_gpu
  end interface

  interface
    subroutine SideExchange_3D_gpu(extBoundary,boundary,sideInfo,elemToRank,rankId,offset,N,nVar,nEl) &
      bind(c,name="SideExchange_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: extBoundary,boundary,sideInfo,elemToRank
      integer(c_int),value :: rankId,offset,N,nVar,nEl
    end subroutine SideExchange_3D_gpu
  end interface

  interface
    subroutine BassiRebaySides_2D_gpu(avgBoundary,boundary,extBoundary,N,nVar,nEl) &
      bind(c,name="BassiRebaySides_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: extBoundary,boundary,avgBoundary
      integer(c_int),value :: N,nVar,nEl
    end subroutine BassiRebaySides_2D_gpu
  end interface

  interface
    subroutine BassiRebaySides_3D_gpu(avgBoundary,boundary,extBoundary,N,nVar,nEl) &
      bind(c,name="BassiRebaySides_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: extBoundary,boundary,avgBoundary
      integer(c_int),value :: N,nVar,nEl
    end subroutine BassiRebaySides_3D_gpu
  end interface

  interface
    subroutine ApplyFlip_2D_gpu(extBoundary,selfSideInfo,elemToRank,rankId,N,nVar,nEl) &
      bind(c,name="ApplyFlip_2D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: selfSideInfo,elemToRank,extBoundary
      integer(c_int),value :: rankId,N,nVar,nEl
    end subroutine ApplyFlip_2D_gpu
  end interface

  interface
    subroutine ApplyFlip_3D_gpu(extBoundary,selfSideInfo,elemToRank,rankId,N,nVar,nEl) &
      bind(c,name="ApplyFlip_3D_gpu")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: selfSideInfo,elemToRank,extBoundary
      integer(c_int),value :: rankId,N,nVar,nEl
    end subroutine ApplyFlip_3D_gpu
  end interface

contains

! ---------------------- Scalars ---------------------- !

  subroutine SetInteriorFromEquation_MappedScalar1D(scalar,geometry,time)
    !!  Sets the scalar % interior attribute using the eqn attribute,
    !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar1D),intent(inout) :: scalar
    type(Geometry1D),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,iEl,iVar

    do ivar = 1,scalar % nvar
      scalar % interior(:,:,ivar) = scalar % eqn(ivar) % evaluate(geometry % x % interior)
    end do

  end subroutine SetInteriorFromEquation_MappedScalar1D

  subroutine SideExchange_MappedScalar1D(scalar,mesh,decomp)
    implicit none
    class(MappedScalar1D),intent(inout) :: scalar
    type(Mesh1D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,ivar
    integer :: neighborRank
    integer :: rankId,offset

    rankId = decomp % rankId
    offset = decomp % offsetElem(rankId + 1)

    do e1 = 1,mesh % nElem

      if (e1 == 1) then

        s1 = 2
        e2 = e1 + 1
        s2 = 1
        !neighborRank = decomp % elemToRank(e2Global)
        do ivar = 1,scalar % nvar
          scalar % extBoundary(s1,e1,ivar) = scalar % boundary(s2,e2,ivar)
        end do

      elseif (e1 == mesh % nElem) then

        s1 = 1
        e2 = e1 - 1
        s2 = 2
        !neighborRank = decomp % elemToRank(e2Global)
        do ivar = 1,scalar % nvar
          scalar % extBoundary(s1,e1,ivar) = scalar % boundary(s2,e2,ivar)
        end do

      else

        s1 = 1
        e2 = e1 - 1
        s2 = 2
        !neighborRank = decomp % elemToRank(e2Global)
        do ivar = 1,scalar % nvar
          scalar % extBoundary(s1,e1,ivar) = scalar % boundary(s2,e2,ivar)
        end do

        s1 = 2
        e2 = e1 + 1
        s2 = 1
        !neighborRank = decomp % elemToRank(e2Global)
        do ivar = 1,scalar % nvar
          scalar % extBoundary(s1,e1,ivar) = scalar % boundary(s2,e2,ivar)
        end do

      end if

    end do

  end subroutine SideExchange_MappedScalar1D

  subroutine BassiRebaySides_MappedScalar1D(scalar)
    implicit none
    class(MappedScalar1D),intent(inout) :: scalar
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i

    do iel = 1,scalar % nElem
      do ivar = 1,scalar % nVar

        ! Left side - we account for the -\hat{x} normal
        scalar % avgBoundary(1,iel,ivar) = -0.5_prec*( &
                                           scalar % boundary(1,iel,ivar) + &
                                           scalar % extBoundary(1,iel,ivar))

        ! Right side - we account for the +\hat{x} normal
        scalar % avgBoundary(2,iel,ivar) = 0.5_prec*( &
                                           scalar % boundary(2,iel,ivar) + &
                                           scalar % extBoundary(2,iel,ivar))
      end do
    end do

  end subroutine BassiRebaySides_MappedScalar1D

  subroutine Derivative_MappedScalar1D(scalar,geometry,dF,handle)
    implicit none
    class(MappedScalar1D),intent(in) :: scalar
    type(Geometry1D),intent(in) :: geometry
    type(MappedScalar1D),intent(inout) :: dF
    type(c_ptr),intent(inout),optional :: handle

    if (present(handle)) then
      call scalar % interp % Derivative_1D(scalar % interior, &
                                           df % interior, &
                                           scalar % nVar, &
                                           scalar % nElem, &
                                           handle)
      call df % JacobianWeight(geometry,handle)

    else

      call scalar % interp % Derivative_1D(scalar % interior, &
                                           df % interior, &
                                           scalar % nVar, &
                                           scalar % nElem)
      call df % JacobianWeight(geometry)

    end if

  end subroutine Derivative_MappedScalar1D

  subroutine DGDerivative_MappedScalar1D(scalar,geometry,dF,handle)
    implicit none
    class(MappedScalar1D),intent(in) :: scalar
    type(Geometry1D),intent(in) :: geometry
    type(MappedScalar1D),intent(inout) :: dF
    type(c_ptr),intent(inout),optional :: handle

    if (present(handle)) then

      call scalar % interp % DGDerivative_1D(scalar % interior, &
                                             scalar % boundary, &
                                             df % interior, &
                                             scalar % nVar, &
                                             scalar % nElem, &
                                             handle)

      call df % JacobianWeight(geometry,handle)

    else

      call scalar % interp % DGDerivative_1D(scalar % interior, &
                                             scalar % boundary, &
                                             df % interior, &
                                             scalar % nVar, &
                                             scalar % nElem)
      call df % JacobianWeight(geometry)

    end if

  end subroutine DGDerivative_MappedScalar1D

  subroutine BRDerivative_MappedScalar1D(scalar,geometry,dF,handle)
    implicit none
    class(MappedScalar1D),intent(in) :: scalar
    type(Geometry1D),intent(in) :: geometry
    type(MappedScalar1D),intent(inout) :: dF
    type(c_ptr),intent(inout),optional :: handle

    if (present(handle)) then

      call scalar % interp % DGDerivative_1D(scalar % interior, &
                                             scalar % avgboundary, &
                                             df % interior, &
                                             scalar % nVar, &
                                             scalar % nElem, &
                                             handle)

      call df % JacobianWeight(geometry,handle)

    else

      call scalar % interp % DGDerivative_1D(scalar % interior, &
                                             scalar % avgboundary, &
                                             df % interior, &
                                             scalar % nVar, &
                                             scalar % nElem)
      call df % JacobianWeight(geometry)

    end if

  end subroutine BRDerivative_MappedScalar1D

  subroutine JacobianWeight_MappedScalar1D(scalar,geometry,handle)
#undef __FUNC__
#define __FUNC__ "JacobianWeight_MappedScalar1D"
    ! Applies the inverse jacobian
    implicit none
    class(MappedScalar1D),intent(inout) :: scalar
    type(Geometry1D),intent(in) :: geometry
    type(c_ptr),intent(inout),optional :: handle
    ! Local
    integer :: iEl,iVar,i

    if (present(handle)) then

      call JacobianWeight_1D_gpu(c_loc(scalar % interior), &
                                 c_loc(geometry % dxds % interior), &
                                 scalar % interp % N, &
                                 scalar % nVar, &
                                 scalar % nElem)
    else

      do iEl = 1,scalar % nElem
        do iVar = 1,scalar % nVar
          do i = 1,scalar % interp % N + 1
            scalar % interior(i,iEl,iVar) = scalar % interior(i,iEl,iVar)/ &
                                            geometry % dxds % interior(i,iEl,1)
          end do
        end do
      end do

    end if

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

    call hipcheck(hipMallocManaged(this % interior,interp % N + 1,interp % N + 1,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % interpWork,interp % M + 1,interp % N + 1,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % boundary,interp % N + 1,4,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % extBoundary,interp % N + 1,4,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % avgBoundary,interp % N + 1,4,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % jumpBoundary,interp % N + 1,4,nelem,nvar,hipMemAttachGlobal))

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
    call hipcheck(hipFree(this % interior))
    call hipcheck(hipFree(this % interpWork))
    call hipcheck(hipFree(this % boundary))
    call hipcheck(hipFree(this % extBoundary))
    call hipcheck(hipFree(this % avgBoundary))
    call hipcheck(hipFree(this % jumpBoundary))
    deallocate (this % meta)
    deallocate (this % eqn)
    call this % JaScalar % Free()

  end subroutine Free_MappedScalar2D

  subroutine SetInteriorFromEquation_MappedScalar2D(scalar,geometry,time)
  !!  Sets the scalar % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar2D),intent(inout) :: scalar
    type(SEMQuad),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,j,iEl,iVar
    real(prec) :: x
    real(prec) :: y

    do iVar = 1,scalar % nVar
      do iEl = 1,scalar % nElem
        do j = 1,scalar % interp % N + 1
          do i = 1,scalar % interp % N + 1

            ! Get the mesh positions
            x = geometry % x % interior(i,j,iEl,1,1)
            y = geometry % x % interior(i,j,iEl,1,2)

            scalar % interior(i,j,iEl,iVar) = &
              scalar % eqn(iVar) % Evaluate((/x,y,0.0_prec,time/))

          end do
        end do
      end do
    end do

  end subroutine SetInteriorFromEquation_MappedScalar2D

  subroutine MPIExchangeAsync_MappedScalar2D(scalar,mpiHandler,mesh,resetCount)
    implicit none
    class(MappedScalar2D),intent(inout) :: scalar
    type(MPILayer),intent(inout) :: mpiHandler
    type(Mesh2D),intent(in) :: mesh
    logical,intent(in) :: resetCount
    ! Local
    integer :: e1,s1,e2,s2,ivar
    integer :: globalSideId,r2
    integer :: iError
    integer :: msgCount

    if (mpiHandler % mpiEnabled) then
      if (resetCount) then
        msgCount = 0
      else
        msgCount = mpiHandler % msgCount
      end if

      do ivar = 1,scalar % nvar
        do e1 = 1,scalar % nElem
          do s1 = 1,4

            e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
            if (e2 > 0) then
              r2 = mpiHandler % elemToRank(e2) ! Neighbor Rank

              if (r2 /= mpiHandler % rankId) then

                ! to do : create unique tag for each side and each variable
                ! tag = globalsideid + nglobalsides*ivar
                s2 = mesh % sideInfo(4,s1,e1)/10
                globalSideId = abs(mesh % sideInfo(2,s1,e1))

                msgCount = msgCount + 1
                call MPI_IRECV(scalar % extBoundary(:,s1,e1,ivar), &
                               (scalar % interp % N + 1), &
                               mpiHandler % mpiPrec, &
                               r2,globalSideId, &
                               mpiHandler % mpiComm, &
                               mpiHandler % requests(msgCount),iError)

                msgCount = msgCount + 1
                call MPI_ISEND(scalar % boundary(:,s1,e1,ivar), &
                               (scalar % interp % N + 1), &
                               mpiHandler % mpiPrec, &
                               r2,globalSideId, &
                               mpiHandler % mpiComm, &
                               mpiHandler % requests(msgCount),iError)
              end if
            end if

          end do
        end do
      end do

      mpiHandler % msgCount = msgCount
    end if

  end subroutine MPIExchangeAsync_MappedScalar2D

  subroutine ApplyFlip_MappedScalar2D(scalar,mpiHandler,mesh,handle)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedScalar2D),intent(inout) :: scalar
    type(MPILayer),intent(inout) :: mpiHandler
    type(Mesh2D),intent(in) :: mesh
    type(c_ptr),intent(inout),optional :: handle
    ! Local
    integer :: e1,s1,e2,s2
    integer :: i,i2
    integer :: r2,flip,ivar
    integer :: globalSideId
    integer :: bcid
    real(prec) :: extBuff(1:scalar % interp % N + 1)

    if (mpiHandler % mpiEnabled) then
      if (present(handle)) then

        call ApplyFlip_2D_gpu(c_loc(scalar % extBoundary), &
                              c_loc(mesh % sideInfo), &
                              c_loc(mpiHandler % elemToRank), &
                              mpiHandler % rankId, &
                              scalar % interp % N, &
                              scalar % nVar, &
                              scalar % nElem)
      else

        do ivar = 1,scalar % nvar
          do e1 = 1,scalar % nElem
            do s1 = 1,4

              e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
              bcid = mesh % sideInfo(5,s1,e1)
              if (bcid == 0) then ! Interior Element
                r2 = mpiHandler % elemToRank(e2) ! Neighbor Rank

                if (r2 /= mpiHandler % rankId) then

                  s2 = mesh % sideInfo(4,s1,e1)/10
                  flip = mesh % sideInfo(4,s1,e1) - s2*10
                  globalSideId = mesh % sideInfo(2,s1,e1)

                  ! Need to update extBoundary with flip applied
                  if (flip == 1) then

                    do i = 1,scalar % interp % N + 1
                      i2 = scalar % interp % N + 2 - i
                      extBuff(i) = scalar % extBoundary(i2,s1,e1,ivar)
                    end do
                    do i = 1,scalar % interp % N + 1
                      scalar % extBoundary(i,s1,e1,ivar) = extBuff(i)
                    end do

                  end if
                end if

              end if

            end do
          end do
        end do
      end if
    end if

  end subroutine ApplyFlip_MappedScalar2D

  subroutine SideExchange_MappedScalar2D(scalar,mesh,decomp,handle)
    implicit none
    class(MappedScalar2D),intent(inout) :: scalar
    type(Mesh2D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    type(c_ptr),intent(inout),optional :: handle
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,ivar
    integer :: neighborRank
    integer :: rankId,offset

    rankId = decomp % rankId
    offset = decomp % offsetElem(rankId + 1)

    call scalar % MPIExchangeAsync(decomp,mesh,resetCount=.true.)

    if (present(handle)) then

      call SideExchange_2D_gpu(c_loc(scalar % extBoundary), &
                               c_loc(scalar % boundary), &
                               c_loc(mesh % sideInfo), &
                               c_loc(decomp % elemToRank), &
                               decomp % rankId, &
                               offset, &
                               scalar % interp % N, &
                               scalar % nvar, &
                               scalar % nElem)

    else

      do ivar = 1,scalar % nvar
        do e1 = 1,mesh % nElem
          do s1 = 1,4
            e2Global = mesh % sideInfo(3,s1,e1)
            e2 = e2Global - offset
            s2 = mesh % sideInfo(4,s1,e1)/10
            flip = mesh % sideInfo(4,s1,e1) - s2*10
            bcid = mesh % sideInfo(5,s1,e1)

            if (bcid == 0) then

              neighborRank = decomp % elemToRank(e2Global)

              if (neighborRank == decomp % rankId) then

                if (flip == 0) then

                  do i1 = 1,scalar % interp % N + 1
                    scalar % extBoundary(i1,s1,e1,ivar) = &
                      scalar % boundary(i1,s2,e2,ivar)
                  end do

                elseif (flip == 1) then

                  do i1 = 1,scalar % interp % N + 1
                    i2 = scalar % interp % N + 2 - i1
                    scalar % extBoundary(i1,s1,e1,ivar) = &
                      scalar % boundary(i2,s2,e2,ivar)
                  end do

                end if

              end if

            end if

          end do
        end do
      end do

    end if

    call decomp % FinalizeMPIExchangeAsync()

    ! Apply side flips for data exchanged with MPI
    if (present(handle)) then
      call scalar % ApplyFlip(decomp,mesh,handle)
    else
      call scalar % ApplyFlip(decomp,mesh)
    end if

  end subroutine SideExchange_MappedScalar2D

  subroutine ContravariantWeightInterior_MappedScalar2D(scalar,geometry,handle)
    !! Computes the scalar multiplied by the contravariant basis vectors
    !!
    implicit none
    class(MappedScalar2D),intent(inout) :: scalar
    type(SEMQuad),intent(in) :: geometry
    type(c_ptr),intent(in),optional :: handle
    ! Local
    integer    :: i,j,ii,iEl,iVar,row,col

    if (present(handle)) then

      call ContravariantWeight_MappedScalar2D_gpu(c_loc(scalar % interior), &
                                                  c_loc(geometry % dsdx % interior), &
                                                  c_loc(scalar % JaScalar % interior), &
                                                  scalar % interp % N + 1,scalar % interp % N + 1, &
                                                  scalar % nVar,scalar % nElem)
    else

      ! Interior
      do col = 1,2
        do row = 1,2
          do iVar = 1,scalar % nVar
            do iEl = 1,scalar % nElem
              do j = 1,scalar % interp % N + 1
                do i = 1,scalar % interp % N + 1

                  scalar % JaScalar % interior(i,j,iel,ivar,row,col) = geometry % dsdx % interior(i,j,iel,1,row,col)* &
                                                                       scalar % interior(i,j,iel,ivar)
                end do
              end do
            end do
          end do
        end do
      end do

    end if
  end subroutine ContravariantWeightInterior_MappedScalar2D

  subroutine ContravariantWeightAvgBoundary_MappedScalar2D(scalar,geometry,handle)
    !! Computes the scalar multiplied by the contravariant basis vectors
    !!
    implicit none
    class(MappedScalar2D),intent(inout) :: scalar
    type(SEMQuad),intent(in) :: geometry
    type(c_ptr),intent(in),optional :: handle
    ! Local
    integer    :: i,j,ii,iEl,iVar,row,col

    if (present(handle)) then

      call ContravariantWeight_MappedScalar2D_gpu(c_loc(scalar % avgBoundary), &
                                                          c_loc(geometry % dsdx % boundary), &
                                                          c_loc(scalar % JaScalar % boundary), &
                                                          scalar % interp % N + 1,4,scalar % nVar,scalar % nElem)
    else

      ! Interior
      do col = 1,2
        do row = 1,2
          do iVar = 1,scalar % nVar
            do iEl = 1,scalar % nElem
              do j = 1,4
                do i = 1,scalar % interp % N + 1

                  scalar % JaScalar % boundary(i,j,iel,ivar,row,col) = geometry % dsdx % boundary(i,j,iel,1,row,col)* &
                                                                       scalar % avgBoundary(i,j,iel,ivar)
                end do
              end do
            end do
          end do
        end do
      end do

    end if

  end subroutine ContravariantWeightAvgBoundary_MappedScalar2D

  subroutine BassiRebaySides_MappedScalar2D(scalar,handle)
    implicit none
    class(MappedScalar2D),intent(inout) :: scalar
    type(c_ptr),intent(inout),optional :: handle
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i

    if (present(handle)) then

      call BassiRebaySides_2D_gpu(c_loc(scalar % avgBoundary), &
                                  c_loc(scalar % boundary), &
                                  c_loc(scalar % extBoundary), &
                                  scalar % interp % N, &
                                  scalar % nvar, &
                                  scalar % nElem)

    else

      do ivar = 1,scalar % nVar
        do iel = 1,scalar % nElem
          do iside = 1,4
            do i = 1,scalar % interp % N + 1
              scalar % avgBoundary(i,iside,iel,ivar) = 0.5_prec*( &
                                                       scalar % boundary(i,iside,iel,ivar) + &
                                                       scalar % extBoundary(i,iside,iel,ivar))
            end do
          end do
        end do
      end do

    end if

  end subroutine BassiRebaySides_MappedScalar2D

  subroutine JacobianWeight_MappedScalar2D(scalar,geometry,handle)

    ! Applies the inverse jacobian
    implicit none
    class(MappedScalar2D),intent(inout) :: scalar
    type(SEMQuad),intent(in) :: geometry
    type(c_ptr),intent(in),optional :: handle
    ! Local
    integer :: iEl,iVar,i,j

    if (present(handle)) then

      call JacobianWeight_2D_gpu(c_loc(scalar % interior), &
                                 c_loc(geometry % J % interior), &
                                 scalar % interp % N, &
                                 scalar % nVar, &
                                 scalar % nElem)
    else

      do iVar = 1,scalar % nVar
        do iEl = 1,scalar % nElem
          do j = 1,scalar % interp % N + 1
            do i = 1,scalar % interp % N + 1
              scalar % interior(i,j,iEl,iVar) = scalar % interior(i,j,iEl,iVar)/ &
                                                geometry % J % interior(i,j,iEl,1)
            end do
          end do
        end do
      end do

    end if

  end subroutine JacobianWeight_MappedScalar2D

  subroutine Gradient_MappedScalar2D(scalar,geometry,df,handle)
  !! Calculates the gradient of a function using the strong form of the gradient
  !! in mapped coordinates.
    implicit none
    class(MappedScalar2D),intent(inout) :: scalar
    type(SEMQuad),intent(in) :: geometry
    type(MappedVector2D),intent(inout) :: df
    type(c_ptr),intent(inout),optional :: handle

    if (present(handle)) then
      call scalar % ContravariantWeightInterior(geometry,handle)
      call scalar % JaScalar % Divergence(df,handle)
      call df % JacobianWeight(geometry,handle)
    else
      call scalar % ContravariantWeightInterior(geometry)
      call scalar % JaScalar % Divergence(df)
      call df % JacobianWeight(geometry)
    end if

  end subroutine Gradient_MappedScalar2D

  subroutine BRGradient_MappedScalar2D(scalar,geometry,df,handle)
    !! Calculates the gradient of a function using the weak form of the gradient
    !! and the average boundary state.
    !! This method will compute the average boundary state from the
    !! and  attributes of
    implicit none
    class(MappedScalar2D),intent(inout) :: scalar
    type(SEMQuad),intent(in) :: geometry
    type(MappedVector2D),intent(inout) :: df
    type(c_ptr),intent(inout),optional :: handle

    if (present(handle)) then
      call scalar % BassiRebaySides(handle)
      call scalar % ContravariantWeightInterior(geometry,handle)
      call scalar % ContravariantWeightAvgBoundary(geometry,handle)
      call scalar % JaScalar % DGDivergence(df,handle)
      call df % JacobianWeight(geometry,handle)
    else
      call scalar % BassiRebaySides()
      call scalar % ContravariantWeightInterior(geometry)
      call scalar % ContravariantWeightAvgBoundary(geometry)
      call scalar % JaScalar % DGDivergence(df)
      call df % JacobianWeight(geometry)
    end if

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

    call hipcheck(hipMallocManaged(this % interior,interp % N + 1,interp % N + 1,interp % N + 1,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % interpWork1,interp % M + 1,interp % N + 1,interp % N + 1,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % interpWork2,interp % M + 1,interp % M + 1,interp % N + 1,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % boundary,interp % N + 1,interp % N + 1,6,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % extBoundary,interp % N + 1,interp % N + 1,6,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % avgBoundary,interp % N + 1,interp % N + 1,6,nelem,nvar,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % jumpBoundary,interp % N + 1,interp % N + 1,6,nelem,nvar,hipMemAttachGlobal))

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
    call hipcheck(hipFree(this % interior))
    call hipcheck(hipFree(this % interpWork1))
    call hipcheck(hipFree(this % interpWork2))
    call hipcheck(hipFree(this % boundary))
    call hipcheck(hipFree(this % extBoundary))
    call hipcheck(hipFree(this % avgBoundary))
    call hipcheck(hipFree(this % jumpBoundary))
    deallocate (this % meta)
    deallocate (this % eqn)
    call this % JaScalar % Free()

  end subroutine Free_MappedScalar3D

  subroutine SetInteriorFromEquation_MappedScalar3D(scalar,geometry,time)
  !!  Sets the scalar % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedScalar3D),intent(inout) :: scalar
    type(SEMHex),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,j,k,iEl,iVar
    real(prec) :: x
    real(prec) :: y
    real(prec) :: z

    do iVar = 1,scalar % nVar
      do iEl = 1,scalar % nElem
        do k = 1,scalar % interp % N + 1
          do j = 1,scalar % interp % N + 1
            do i = 1,scalar % interp % N + 1

              ! Get the mesh positions
              x = geometry % x % interior(i,j,k,iEl,1,1)
              y = geometry % x % interior(i,j,k,iEl,1,2)
              z = geometry % x % interior(i,j,k,iEl,1,3)

              scalar % interior(i,j,k,iEl,iVar) = &
                scalar % eqn(iVar) % Evaluate((/x,y,z,time/))

            end do
          end do
        end do
      end do
    end do

  end subroutine SetInteriorFromEquation_MappedScalar3D

  subroutine MPIExchangeAsync_MappedScalar3D(scalar,mpiHandler,mesh,resetCount)
    implicit none
    class(MappedScalar3D),intent(inout) :: scalar
    type(MPILayer),intent(inout) :: mpiHandler
    type(Mesh3D),intent(in) :: mesh
    logical,intent(in) :: resetCount
    ! Local
    integer :: e1,s1,e2,s2,ivar
    integer :: globalSideId,r2
    integer :: iError
    integer :: msgCount

    if (mpiHandler % mpiEnabled) then
      if (resetCount) then
        msgCount = 0
      else
        msgCount = mpiHandler % msgCount
      end if

      do ivar = 1,scalar % nvar
        do e1 = 1,scalar % nElem
          do s1 = 1,6

            e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
            if (e2 > 0) then
              r2 = mpiHandler % elemToRank(e2) ! Neighbor Rank

              if (r2 /= mpiHandler % rankId) then

                ! to do : create unique tag for each side and each variable
                ! tag = globalsideid + nglobalsides*ivar
                s2 = mesh % sideInfo(4,s1,e1)/10
                globalSideId = abs(mesh % sideInfo(2,s1,e1))

                msgCount = msgCount + 1
                call MPI_IRECV(scalar % extBoundary(:,:,s1,e1,ivar), &
                               (scalar % interp % N + 1)*(scalar % interp % N + 1), &
                               mpiHandler % mpiPrec, &
                               r2,globalSideId, &
                               mpiHandler % mpiComm, &
                               mpiHandler % requests(msgCount),iError)

                msgCount = msgCount + 1
                call MPI_ISEND(scalar % boundary(:,:,s1,e1,ivar), &
                               (scalar % interp % N + 1)*(scalar % interp % N + 1), &
                               mpiHandler % mpiPrec, &
                               r2,globalSideId, &
                               mpiHandler % mpiComm, &
                               mpiHandler % requests(msgCount),iError)
              end if
            end if

          end do
        end do
      end do

      mpiHandler % msgCount = msgCount
    end if

  end subroutine MPIExchangeAsync_MappedScalar3D

  subroutine ApplyFlip_MappedScalar3D(scalar,mpiHandler,mesh,handle)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedScalar3D),intent(inout) :: scalar
    type(MPILayer),intent(inout) :: mpiHandler
    type(Mesh3D),intent(in) :: mesh
    type(c_ptr),intent(inout),optional :: handle
    ! Local
    integer :: e1,s1,e2,s2
    integer :: i,i2,j,j2
    integer :: r2,flip,ivar
    integer :: globalSideId
    integer :: bcid
    real(prec) :: extBuff(1:scalar % interp % N + 1,1:scalar % interp % N + 1)

    if (mpiHandler % mpiEnabled) then
      if (present(handle)) then

        call ApplyFlip_3D_gpu(c_loc(scalar % extBoundary), &
                              c_loc(mesh % sideInfo), &
                              c_loc(mpiHandler % elemToRank), &
                              mpiHandler % rankId, &
                              scalar % interp % N, &
                              scalar % nVar, &
                              scalar % nElem)
      else

        do ivar = 1,scalar % nvar
          do e1 = 1,scalar % nElem
            do s1 = 1,6

              e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
              bcid = mesh % sideInfo(5,s1,e1)
              if (bcid == 0) then ! Interior Element
                r2 = mpiHandler % elemToRank(e2) ! Neighbor Rank

                if (r2 /= mpiHandler % rankId) then

                  s2 = mesh % sideInfo(4,s1,e1)/10
                  flip = mesh % sideInfo(4,s1,e1) - s2*10
                  globalSideId = mesh % sideInfo(2,s1,e1)

                  ! Need to update extBoundary with flip applied
                  if (flip == 1) then

                    DO j = 1,scalar % interp % N+1
                      DO i = 1,scalar % interp % N+1
                        i2 = j
                        j2 = scalar % interp % N+2-i
                        extBuff(i,j) = scalar % extBoundary(i2,j2,s1,e1,ivar)
                      END DO
                    END DO


                  else if (flip == 2) then

                    DO j = 1,scalar % interp % N+1
                      DO i = 1,scalar % interp % N+1
                        i2 = scalar % interp % N + 2 - i
                        j2 = scalar % interp % N + 2 - j
                        extBuff(i,j) = scalar % extBoundary(i2,j2,s1,e1,ivar)
                      END DO
                    END DO

                  else if (flip == 3) then

                    DO j = 1,scalar % interp % N+1
                      DO i = 1,scalar % interp % N+1
                        i2 = scalar % interp % N + 2 - j
                        j2 = i
                        extBuff(i,j) = scalar % extBoundary(i2,j2,s1,e1,ivar)
                      END DO
                    END DO

                  else if (flip == 4) then

                    DO j = 1,scalar % interp % N+1
                      DO i = 1,scalar % interp % N+1
                        i2 = j
                        j2 = i
                        extBuff(i,j) = scalar % extBoundary(i2,j2,s1,e1,ivar)
                      END DO
                    END DO

                  end if

                  DO j = 1,scalar % interp % N + 1
                    DO i = 1,scalar % interp % N + 1
                      scalar % extBoundary(i,j,s1,e1,ivar) = extBuff(i,j)
                    END DO
                  END DO

                end if

              end if

            end do
          end do
        end do
      end if
    end if

  end subroutine ApplyFlip_MappedScalar3D

  subroutine SideExchange_MappedScalar3D(scalar,mesh,decomp,handle)
    implicit none
    class(MappedScalar3D),intent(inout) :: scalar
    type(Mesh3D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    type(c_ptr),intent(inout),optional :: handle
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,j1,j2,ivar
    integer :: neighborRank
    integer :: rankId,offset

    rankId = decomp % rankId
    offset = decomp % offsetElem(rankId + 1)

    call scalar % MPIExchangeAsync(decomp,mesh,resetCount=.true.)

    if (present(handle)) then

      call SideExchange_3D_gpu(c_loc(scalar % extBoundary), &
                               c_loc(scalar % boundary), &
                               c_loc(mesh % sideInfo), &
                               c_loc(decomp % elemToRank), &
                               decomp % rankId, &
                               offset, &
                               scalar % interp % N, &
                               scalar % nvar, &
                               scalar % nElem)

    else

      do ivar = 1,scalar % nvar
        do e1 = 1,mesh % nElem
          do s1 = 1,6
            e2Global = mesh % sideInfo(3,s1,e1)
            e2 = e2Global - offset
            s2 = mesh % sideInfo(4,s1,e1)/10
            flip = mesh % sideInfo(4,s1,e1) - s2*10
            bcid = mesh % sideInfo(5,s1,e1)

            if (bcid == 0) then

              neighborRank = decomp % elemToRank(e2Global)

              if (neighborRank == decomp % rankId) then

                if (flip == 0)then

                  DO j1 = 1,scalar % interp % N+1
                    DO i1 = 1,scalar % interp % N+1
                      scalar % extBoundary(i1,j1,s1,e1,ivar) = &
                        scalar % boundary(i1,j1,s2,e2,ivar)
                    END DO
                  END DO

                else if (flip == 1)then

                  DO j1 = 1,scalar % interp % N+1
                    DO i1 = 1,scalar % interp % N+1

                      i2 = j1
                      j2 = scalar % interp % N + 2 - i1
                      scalar % extBoundary(i1,j1,s1,e1,ivar) = &
                        scalar % boundary(i2,j2,s2,e2,ivar)

                    END DO
                  END DO

                else if (flip == 2)then

                  DO j1 = 1,scalar % interp % N+1
                    DO i1 = 1,scalar % interp % N+1
                      i2 = scalar % interp % N + 2 - i1
                      j2 = scalar % interp % N + 2 - j1
                      scalar % extBoundary(i1,j1,s1,e1,ivar) = &
                        scalar % boundary(i2,j2,s2,e2,ivar)
                    END DO
                  END DO

                else if (flip == 3)then

                  DO j1 = 1,scalar % interp % N+1
                    DO i1 = 1,scalar % interp % N+1
                      i2 = scalar % interp % N + 2 - j1
                      j2 = i1
                      scalar % extBoundary(i1,j1,s1,e1,ivar) = &
                        scalar % boundary(i2,j2,s2,e2,ivar)
                    END DO
                  END DO

                else if (flip == 4)then

                  DO j1 = 1,scalar % interp % N+1
                    DO i1 = 1,scalar % interp % N+1
                      i2 = j1
                      j2 = i1
                      scalar % extBoundary(i1,j1,s1,e1,ivar) = &
                        scalar % boundary(i2,j2,s2,e2,ivar)
                    END DO
                  END DO

                end if

              end if

            end if

          end do
        end do
      end do

    end if

    call decomp % FinalizeMPIExchangeAsync()

    ! Apply side flips for data exchanged with MPI
    if (present(handle)) then
      call scalar % ApplyFlip(decomp,mesh,handle)
    else
      call scalar % ApplyFlip(decomp,mesh)
    end if

  end subroutine SideExchange_MappedScalar3D

  subroutine ContravariantWeightInterior_MappedScalar3D(scalar,geometry,handle)
    !! Computes the scalar multiplied by the contravariant basis vectors
    !!
    implicit none
    class(MappedScalar3D),intent(inout) :: scalar
    type(SEMHex),intent(in) :: geometry
    type(c_ptr),intent(in),optional :: handle
    ! Local
    integer    :: i,j,k,iEl,iVar,row,col

    if (present(handle)) then

      call ContravariantWeight_MappedScalar3D_gpu(c_loc(scalar % interior), &
                                                          c_loc(geometry % dsdx % interior), &
                                                          c_loc(scalar % JaScalar % interior), &
                                           scalar % interp % N + 1,scalar % interp % N + 1,scalar % interp % N + 1,scalar % nVar,scalar % nElem)
    else

      ! Interior
      do col = 1,3
        do row = 1,3
          do iVar = 1,scalar % nVar
            do iEl = 1,scalar % nElem
              do k = 1,scalar % interp % N + 1
                do j = 1,scalar % interp % N + 1
                  do i = 1,scalar % interp % N + 1

                    scalar % JaScalar % interior(i,j,k,iel,ivar,row,col) = geometry % dsdx % interior(i,j,k,iel,1,row,col)* &
                                                                        scalar % interior(i,j,k,iel,ivar)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

    end if


  end subroutine ContravariantWeightInterior_MappedScalar3D

  subroutine ContravariantWeightAvgBoundary_MappedScalar3D(scalar,geometry,handle)
    !! Computes the scalar multiplied by the contravariant basis vectors
    !!
    implicit none
    class(MappedScalar3D),intent(inout) :: scalar
    type(SEMHex),intent(in) :: geometry
    type(c_ptr),intent(in),optional :: handle
    ! Local
    integer    :: i,j,k,iEl,iVar,row,col

    if (present(handle)) then

      call ContravariantWeight_MappedScalar3D_gpu(c_loc(scalar % avgBoundary), &
                                                          c_loc(geometry % dsdx % boundary), &
                                                          c_loc(scalar % JaScalar % boundary), &
                                                          scalar % interp % N + 1,scalar % interp % N + 1,6,scalar % nVar,scalar % nElem)
    else

      ! Interior
      do col = 1,3
        do row = 1,3
          do iVar = 1,scalar % nVar
            do iEl = 1,scalar % nElem
              do k = 1,6
                do j = 1,scalar % interp % N + 1
                  do i = 1,scalar % interp % N + 1

                    scalar % JaScalar % boundary(i,j,k,iel,ivar,row,col) = geometry % dsdx % boundary(i,j,k,iel,1,row,col)* &
                                                                           scalar % avgBoundary(i,j,k,iel,ivar)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

    end if

  end subroutine ContravariantWeightAvgBoundary_MappedScalar3D

  subroutine BassiRebaySides_MappedScalar3D(scalar,handle)
    implicit none
    class(MappedScalar3D),intent(inout) :: scalar
    type(c_ptr),intent(inout),optional :: handle
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i, j

    if (present(handle)) then

      call BassiRebaySides_3D_gpu(c_loc(scalar % avgBoundary), &
                                  c_loc(scalar % boundary), &
                                  c_loc(scalar % extBoundary), &
                                  scalar % interp % N, &
                                  scalar % nvar, &
                                  scalar % nElem)

    else

      do ivar = 1,scalar % nVar
        do iel = 1,scalar % nElem
          do iside = 1,6
            do j = 1,scalar % interp % N + 1
              do i = 1,scalar % interp % N + 1
                scalar % avgBoundary(i,j,iside,iel,ivar) = 0.5_prec*( &
                                                        scalar % boundary(i,j,iside,iel,ivar) + &
                                                        scalar % extBoundary(i,j,iside,iel,ivar))
              end do
            end do
          end do
        end do
      end do

    end if

  end subroutine BassiRebaySides_MappedScalar3D

  subroutine JacobianWeight_MappedScalar3D(scalar,geometry,handle)

    ! Applies the inverse jacobian
    implicit none
    class(MappedScalar3D),intent(inout) :: scalar
    type(SEMHex),intent(in) :: geometry
    type(c_ptr),intent(in),optional :: handle
    ! Local
    integer :: iEl,iVar,i,j,k

    if (present(handle)) then

      call JacobianWeight_3D_gpu(c_loc(scalar % interior), &
                                 c_loc(geometry % J % interior), &
                                 scalar % interp % N, &
                                 scalar % nVar, &
                                 scalar % nElem)
    else

      do iEl = 1,scalar % nElem
        do iVar = 1,scalar % nVar
          do k = 1,scalar % interp % N + 1
            do j = 1,scalar % interp % N + 1
              do i = 1,scalar % interp % N + 1
                scalar % interior(i,j,k,iEl,iVar) = scalar % interior(i,j,k,iEl,iVar)/ &
                                                    geometry % J % interior(i,j,k,iEl,1)
              end do
            end do
          end do
        end do
      end do

    end if

  end subroutine JacobianWeight_MappedScalar3D

  subroutine Gradient_MappedScalar3D(scalar,geometry,df,handle)
  !! Calculates the gradient of a function using the strong form of the gradient
  !! in mapped coordinates.
    implicit none
    class(MappedScalar3D),intent(inout) :: scalar
    type(SEMHex),intent(in) :: geometry
    type(MappedVector3D),intent(inout) :: df
    type(c_ptr),intent(inout),optional :: handle

    if (present(handle)) then
      call scalar % ContravariantWeightInterior(geometry,handle)
      call scalar % JaScalar % Divergence(df,handle)
      call df % JacobianWeight(geometry,handle)
    else
      call scalar % ContravariantWeightInterior(geometry)
      call scalar % JaScalar % Divergence(df)
      call df % JacobianWeight(geometry)
    end if

  end subroutine Gradient_MappedScalar3D

  subroutine BRGradient_MappedScalar3D(scalar,geometry,df,handle)
    !! Calculates the gradient of a function using the weak form of the gradient
    !! and the average boundary state.
    !! This method will compute the average boundary state from the
    !! and  attributes of
    implicit none
    class(MappedScalar3D),intent(inout) :: scalar
    type(SEMHex),intent(in) :: geometry
    type(MappedVector3D),intent(inout) :: df
    type(c_ptr),intent(inout),optional :: handle

    if (present(handle)) then
      call scalar % BassiRebaySides(handle)
      call scalar % ContravariantWeightInterior(geometry,handle)
      call scalar % ContravariantWeightAvgBoundary(geometry,handle)
      call scalar % JaScalar % DGDivergence(df,handle)
      call df % JacobianWeight(geometry,handle)
    else
      call scalar % BassiRebaySides()
      call scalar % ContravariantWeightInterior(geometry)
      call scalar % ContravariantWeightAvgBoundary(geometry)
      call scalar % JaScalar % DGDivergence(df)
      call df % JacobianWeight(geometry)
    end if

  end subroutine BRGradient_MappedScalar3D

  subroutine SetInteriorFromEquation_MappedVector2D(vector,geometry,time)
  !!  Sets the scalar % interior attribute using the eqn attribute,
  !!  geometry (for physical positions), and provided simulation time.
    implicit none
    class(MappedVector2D),intent(inout) :: vector
    type(SEMQuad),intent(in) :: geometry
    real(prec),intent(in) :: time
    ! Local
    integer :: i,j,iEl,iVar
    real(prec) :: x
    real(prec) :: y

    do iVar = 1,vector % nVar
      do iEl = 1,vector % nElem
        do j = 1,vector % interp % N + 1
          do i = 1,vector % interp % N + 1

            ! Get the mesh positions
            x = geometry % x % interior(i,j,iEl,1,1)
            y = geometry % x % interior(i,j,iEl,1,2)

            vector % interior(i,j,iEl,iVar,1) = &
              vector % eqn(1 + 2*(iVar - 1)) % Evaluate((/x,y,0.0_prec,time/))

            vector % interior(i,j,iEl,iVar,2) = &
              vector % eqn(2 + 2*(iVar - 1)) % Evaluate((/x,y,0.0_prec,time/))

          end do
        end do
      end do
    end do

  end subroutine SetInteriorFromEquation_MappedVector2D

  subroutine MPIExchangeAsync_MappedVector2D(this,mpiHandler,mesh,resetCount)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(MPILayer),intent(inout) :: mpiHandler
    type(Mesh2D),intent(in) :: mesh
    logical,intent(in) :: resetCount
    ! Local
    integer :: e1,s1,e2,s2,ivar,idir
    integer :: globalSideId,r2
    integer :: iError
    integer :: msgCount

    if (mpiHandler % mpiEnabled) then
      if (resetCount) then
        msgCount = 0
      else
        msgCount = mpiHandler % msgCount
      end if

      do idir = 1,2
        do ivar = 1,this % nvar
          do e1 = 1,this % nElem
            do s1 = 1,4

              e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
              if (e2 > 0) then
                r2 = mpiHandler % elemToRank(e2) ! Neighbor Rank

                if (r2 /= mpiHandler % rankId) then

                  ! to do : create unique tag for each side and each variable
                  ! tag = globalsideid + nglobalsides*(ivar + nvar*idir)
                  s2 = mesh % sideInfo(4,s1,e1)/10
                  globalSideId = abs(mesh % sideInfo(2,s1,e1))

                  msgCount = msgCount + 1
                  call MPI_IRECV(this % extBoundary(:,s1,e1,ivar,idir), &
                                 (this % interp % N + 1), &
                                 mpiHandler % mpiPrec, &
                                 r2,globalSideId, &
                                 mpiHandler % mpiComm, &
                                 mpiHandler % requests(msgCount),iError)

                  msgCount = msgCount + 1
                  call MPI_ISEND(this % boundary(:,s1,e1,ivar,idir), &
                                 (this % interp % N + 1), &
                                 mpiHandler % mpiPrec, &
                                 r2,globalSideId, &
                                 mpiHandler % mpiComm, &
                                 mpiHandler % requests(msgCount),iError)
                end if
              end if

            end do
          end do
        end do
      end do

      mpiHandler % msgCount = msgCount
    end if

  end subroutine MPIExchangeAsync_MappedVector2D

  subroutine ApplyFlip_MappedVector2D(this,mpiHandler,mesh,handle)
    ! Apply side flips to sides where MPI exchanges took place.
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(MPILayer),intent(inout) :: mpiHandler
    type(Mesh2D),intent(in) :: mesh
    type(c_ptr),intent(inout),optional :: handle
    ! Local
    integer :: e1,s1,e2,s2
    integer :: i,i2
    integer :: r2,flip,ivar,idir
    integer :: globalSideId
    integer :: bcid
    real(prec) :: extBuff(1:this % interp % N + 1)

    if (mpiHandler % mpiEnabled) then
      if (present(handle)) then

        call ApplyFlip_2D_gpu(c_loc(this % extBoundary), &
                              c_loc(mesh % sideInfo), &
                              c_loc(mpiHandler % elemToRank), &
                              mpiHandler % rankId, &
                              this % interp % N, &
                              2*this % nVar, &
                              this % nElem)
      else

        do idir = 1,2
          do ivar = 1,this % nvar
            do e1 = 1,this % nElem
              do s1 = 1,4

                e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
                bcid = mesh % sideInfo(5,s1,e1)
                if (bcid == 0) then ! Interior Element
                  r2 = mpiHandler % elemToRank(e2) ! Neighbor Rank

                  if (r2 /= mpiHandler % rankId) then

                    s2 = mesh % sideInfo(4,s1,e1)/10
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
      end if
    end if

  end subroutine ApplyFlip_MappedVector2D

  subroutine SideExchange_MappedVector2D(this,mesh,decomp,handle)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(Mesh2D),intent(in) :: mesh
    type(MPILayer),intent(inout) :: decomp
    type(c_ptr),intent(inout),optional :: handle
    ! Local
    integer :: e1,e2,s1,s2,e2Global
    integer :: flip,bcid
    integer :: i1,i2,ivar,idir
    integer :: neighborRank
    integer :: rankId,offset

    rankId = decomp % rankId
    offset = decomp % offsetElem(rankId + 1)

    call this % MPIExchangeAsync(decomp,mesh,resetCount=.true.)

    if (present(handle)) then

      call SideExchange_2D_gpu(c_loc(this % extBoundary), &
                               c_loc(this % boundary), &
                               c_loc(mesh % sideInfo), &
                               c_loc(decomp % elemToRank), &
                               decomp % rankId, &
                               offset, &
                               this % interp % N, &
                               2*this % nvar, &
                               this % nElem)

    else

      do idir = 1,2
        do ivar = 1,this % nvar
          do e1 = 1,mesh % nElem
            do s1 = 1,4
              e2Global = mesh % sideInfo(3,s1,e1)
              e2 = e2Global - offset
              s2 = mesh % sideInfo(4,s1,e1)/10
              flip = mesh % sideInfo(4,s1,e1) - s2*10
              bcid = mesh % sideInfo(5,s1,e1)

              if (bcid == 0) then

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

    end if

    call decomp % FinalizeMPIExchangeAsync()

    ! Apply side flips for data exchanged with MPI
    if (present(handle)) then
      call this % ApplyFlip(decomp,mesh,handle)
    else
      call this % ApplyFlip(decomp,mesh)
    end if

  end subroutine SideExchange_MappedVector2D

  subroutine BassiRebaySides_MappedVector2D(this,handle)
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(c_ptr),intent(inout),optional :: handle
    ! Local
    integer :: iel
    integer :: iside
    integer :: ivar
    integer :: i
    integer :: idir

    if (present(handle)) then

      call BassiRebaySides_2D_gpu(c_loc(this % avgBoundary), &
                                  c_loc(this % boundary), &
                                  c_loc(this % extBoundary), &
                                  this % interp % N, &
                                  2*this % nvar, &
                                  this % nElem)

    else

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

    end if

  end subroutine BassiRebaySides_MappedVector2D

  subroutine JacobianWeight_MappedVector2D(this,geometry,handle)

    ! Applies the inverse jacobian
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    type(c_ptr),intent(in),optional :: handle
    ! Local
    integer :: iEl,iVar,i,j,idir

    if (present(handle)) then

      call JacobianWeight_2D_gpu(c_loc(this % interior), &
                                 c_loc(geometry % J % interior), &
                                 this % interp % N, &
                                 2*this % nVar, &
                                 this % nElem)
    else

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

    end if

  end subroutine JacobianWeight_MappedVector2D

  subroutine ContravariantProjection_MappedVector2D(vector,geometry,handle)
    ! Takes a vector that has physical space coordinate directions (x,y,z) and projects the vector
    ! into the the contravariant basis vector directions. Keep in mind that the contravariant basis
    ! vectors are really the Jacobian weighted contravariant basis vectors
    implicit none
    class(MappedVector2D),intent(inout) :: vector
    type(SEMQuad),intent(in) :: geometry
    type(c_ptr),intent(in),optional :: handle
    ! Local
    integer :: i,j,iEl,iVar
    real(prec) :: Fx,Fy

    if (present(handle)) then

      call ContravariantProjection_MappedVector2D_gpu_wrapper(c_loc(vector % interior), &
                                                              c_loc(geometry % dsdx % interior), &
                                                              vector % interp % N, &
                                                              vector % nVar, &
                                                              vector % nElem)

    else
      ! Assume that tensor(j,i) is vector i, component j
      ! => dot product is done along first dimension
      ! to project onto computational space
      do iel = 1,vector % nElem
        do ivar = 1,vector % nVar
          do j = 1,vector % interp % N + 1
            do i = 1,vector % interp % N + 1

              Fx = vector % interior(i,j,iEl,iVar,1)
              Fy = vector % interior(i,j,iEl,iVar,2)

              vector % interior(i,j,iEl,iVar,1) = &
                geometry % dsdx % interior(i,j,iEl,1,1,1)*Fx + &
                geometry % dsdx % interior(i,j,iEl,1,2,1)*Fy

              vector % interior(i,j,iEl,iVar,2) = &
                geometry % dsdx % interior(i,j,iEl,1,1,2)*Fx + &
                geometry % dsdx % interior(i,j,iEl,1,2,2)*Fy

            end do
          end do
        end do
      end do

    end if

  end subroutine ContravariantProjection_MappedVector2D

  subroutine Divergence_MappedVector2D(this,geometry,divVector,handle)
    ! Strong Form Operator
    !    !
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    type(MappedScalar2D),intent(inout) :: divVector
    type(c_ptr),intent(inout),optional :: handle

    if (present(handle)) then

      ! Convert from physical to computational space
      call this % ContravariantProjection(geometry,handle)

      ! Compute the divergence
      call this % interp % VectorDivergence_2D(this % interior, &
                                               divVector % interior, &
                                               this % nvar, &
                                               this % nelem, &
                                               handle)

      ! Divide by the jacobian
      call divVector % JacobianWeight(geometry,handle)

    else

      ! Convert from physical to computational space
      call this % ContravariantProjection(geometry)

      ! Compute the divergence
      call this % interp % VectorDivergence_2D(this % interior, &
                                               divVector % interior, &
                                               this % nvar, &
                                               this % nelem)
      ! Divide by the jacobian
      call divVector % JacobianWeight(geometry)

    end if

  end subroutine Divergence_MappedVector2D

  subroutine DGDivergence_MappedVector2D(this,geometry,divVector,handle)
    !! Computes the divergence of a 2-D vector using the weak form
    !! On input, the  attribute of the vector
    !! is assigned and the  attribute is set to the physical
    !! directions of the vector. This method will project the vector
    !! onto the contravariant basis vectors.
    implicit none
    class(MappedVector2D),intent(inout) :: this
    type(SEMQuad),intent(in) :: geometry
    type(MappedScalar2D),intent(inout) :: divVector
    type(c_ptr),intent(inout),optional :: handle

    if (present(handle)) then

      ! Convert from physical to computational space
      call this % ContravariantProjection(geometry,handle)

      ! Compute the divergence
      call this % interp % VectorDGDivergence_2D(this % interior, &
                                                 this % boundaryNormal, &
                                                 divVector % interior, &
                                                 this % nvar, &
                                                 this % nelem, &
                                                 handle)

      ! Divide by the jacobian
      call divVector % JacobianWeight(geometry,handle)

    else

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

    end if

  end subroutine DGDivergence_MappedVector2D

  subroutine SetInteriorFromEquation_MappedVector3D(vector,geometry,time)
    !!  Sets the scalar % interior attribute using the eqn attribute,
    !!  geometry (for physical positions), and provided simulation time.
      implicit none
      class(MappedVector3D),intent(inout) :: vector
      type(SEMHex),intent(in) :: geometry
      real(prec),intent(in) :: time
      ! Local
      integer :: i,j,k,iEl,iVar
      real(prec) :: x
      real(prec) :: y
      real(prec) :: z
  
      do iVar = 1,vector % nVar
        do iEl = 1,vector % nElem
          do k = 1,vector % interp % N + 1
            do j = 1,vector % interp % N + 1
              do i = 1,vector % interp % N + 1
    
                ! Get the mesh positions
                x = geometry % x % interior(i,j,k,iEl,1,1)
                y = geometry % x % interior(i,j,k,iEl,1,2)
                z = geometry % x % interior(i,j,k,iEl,1,3)
    
                vector % interior(i,j,k,iEl,iVar,1) = &
                  vector % eqn(1 + 3*(iVar - 1)) % Evaluate((/x,y,z,time/))
    
                vector % interior(i,j,k,iEl,iVar,2) = &
                  vector % eqn(2 + 3*(iVar - 1)) % Evaluate((/x,y,z,time/))

                vector % interior(i,j,k,iEl,iVar,3) = &
                  vector % eqn(3 + 3*(iVar - 1)) % Evaluate((/x,y,z,time/))
    
              end do
            end do
          end do
        end do
      end do
  
    end subroutine SetInteriorFromEquation_MappedVector3D

    subroutine MPIExchangeAsync_MappedVector3D(this,mpiHandler,mesh,resetCount)
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(MPILayer),intent(inout) :: mpiHandler
      type(Mesh3D),intent(in) :: mesh
      logical,intent(in) :: resetCount
      ! Local
      integer :: e1,s1,e2,s2,ivar,idir
      integer :: globalSideId,r2
      integer :: iError
      integer :: msgCount
  
      if (mpiHandler % mpiEnabled) then
        if (resetCount) then
          msgCount = 0
        else
          msgCount = mpiHandler % msgCount
        end if
  
        do idir = 1,3
          do ivar = 1,this % nvar
            do e1 = 1,this % nElem
              do s1 = 1,6
    
                e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
                if (e2 > 0) then
                  r2 = mpiHandler % elemToRank(e2) ! Neighbor Rank
    
                  if (r2 /= mpiHandler % rankId) then
    
                    ! to do : create unique tag for each side and each variable
                    ! tag = globalsideid + nglobalsides*ivar
                    s2 = mesh % sideInfo(4,s1,e1)/10
                    globalSideId = abs(mesh % sideInfo(2,s1,e1))
    
                    msgCount = msgCount + 1
                    call MPI_IRECV(this % extBoundary(:,:,s1,e1,ivar,idir), &
                                  (this % interp % N + 1)*(this % interp % N + 1), &
                                  mpiHandler % mpiPrec, &
                                  r2,globalSideId, &
                                  mpiHandler % mpiComm, &
                                  mpiHandler % requests(msgCount),iError)
    
                    msgCount = msgCount + 1
                    call MPI_ISEND(this % boundary(:,:,s1,e1,ivar,idir), &
                                  (this % interp % N + 1)*(this % interp % N + 1), &
                                  mpiHandler % mpiPrec, &
                                  r2,globalSideId, &
                                  mpiHandler % mpiComm, &
                                  mpiHandler % requests(msgCount),iError)
                  end if
                end if
    
              end do
            end do
          end do
        end do
  
        mpiHandler % msgCount = msgCount
      end if
  
    end subroutine MPIExchangeAsync_MappedVector3D
  
    subroutine ApplyFlip_MappedVector3D(this,mpiHandler,mesh,handle)
      ! Apply side flips to sides where MPI exchanges took place.
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(MPILayer),intent(inout) :: mpiHandler
      type(Mesh3D),intent(in) :: mesh
      type(c_ptr),intent(inout),optional :: handle
      ! Local
      integer :: e1,s1,e2,s2
      integer :: i,j,i2,j2
      integer :: r2,flip,ivar
      integer :: globalSideId
      integer :: bcid, idir
      real(prec) :: extBuff(1:this % interp % N + 1,1:this % interp % N + 1)
  
      if (mpiHandler % mpiEnabled) then
        if (present(handle)) then
  
          call ApplyFlip_3D_gpu(c_loc(this % extBoundary), &
                                c_loc(mesh % sideInfo), &
                                c_loc(mpiHandler % elemToRank), &
                                mpiHandler % rankId, &
                                this % interp % N, &
                                3*this % nVar, &
                                this % nElem)
        else
  
          do idir = 1,3
            do ivar = 1,this % nvar
              do e1 = 1,this % nElem
                do s1 = 1,6
    
                  e2 = mesh % sideInfo(3,s1,e1) ! Neighbor Element
                  bcid = mesh % sideInfo(5,s1,e1)
                  if (bcid == 0) then ! Interior Element
                    r2 = mpiHandler % elemToRank(e2) ! Neighbor Rank
    
                    if (r2 /= mpiHandler % rankId) then
    
                      s2 = mesh % sideInfo(4,s1,e1)/10
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
        end if
      end if
  
    end subroutine ApplyFlip_MappedVector3D
  
    subroutine SideExchange_MappedVector3D(this,mesh,decomp,handle)
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(Mesh3D),intent(in) :: mesh
      type(MPILayer),intent(inout) :: decomp
      type(c_ptr),intent(inout),optional :: handle
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
  
      if (present(handle)) then
  
        call SideExchange_3D_gpu(c_loc(this % extBoundary), &
                                 c_loc(this % boundary), &
                                 c_loc(mesh % sideInfo), &
                                 c_loc(decomp % elemToRank), &
                                 decomp % rankId, &
                                 offset, &
                                 this % interp % N, &
                                 3*this % nvar, &
                                 this % nElem)
  
      else

        do idir = 1, 3
          do ivar = 1,this % nvar
            do e1 = 1,mesh % nElem
              do s1 = 1,6
                e2Global = mesh % sideInfo(3,s1,e1)
                e2 = e2Global - offset
                s2 = mesh % sideInfo(4,s1,e1)/10
                flip = mesh % sideInfo(4,s1,e1) - s2*10
                bcid = mesh % sideInfo(5,s1,e1)
    
                if (bcid == 0) then
    
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
  
      end if
  
      call decomp % FinalizeMPIExchangeAsync()
  
      ! Apply side flips for data exchanged with MPI
      if (present(handle)) then
        call this % ApplyFlip(decomp,mesh,handle)
      else
        call this % ApplyFlip(decomp,mesh)
      end if
  
    end subroutine SideExchange_MappedVector3D

    subroutine BassiRebaySides_MappedVector3D(this,handle)
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(c_ptr),intent(inout),optional :: handle
      ! Local
      integer :: iel
      integer :: iside
      integer :: ivar
      integer :: i, j
      integer :: idir
  
      if (present(handle)) then
  
        call BassiRebaySides_2D_gpu(c_loc(this % avgBoundary), &
                                    c_loc(this % boundary), &
                                    c_loc(this % extBoundary), &
                                    this % interp % N, &
                                    3*this % nvar, &
                                    this % nElem)
  
      else
  
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
  
      end if
  
    end subroutine BassiRebaySides_MappedVector3D

    subroutine JacobianWeight_MappedVector3D(this,geometry,handle)

      ! Applies the inverse jacobian
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(SEMHex),intent(in) :: geometry
      type(c_ptr),intent(in),optional :: handle
      ! Local
      integer :: iEl,iVar,i,j,k,idir
  
      if (present(handle)) then
  
        call JacobianWeight_3D_gpu(c_loc(this % interior), &
                                   c_loc(geometry % J % interior), &
                                   this % interp % N, &
                                   3*this % nVar, &
                                   this % nElem)
      else
  
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
  
      end if
  
    end subroutine JacobianWeight_MappedVector3D

    subroutine ContravariantProjection_MappedVector3D(vector,geometry,handle)
      ! Takes a vector that has physical space coordinate directions (x,y,z) and projects the vector
      ! into the the contravariant basis vector directions. Keep in mind that the contravariant basis
      ! vectors are really the Jacobian weighted contravariant basis vectors
      implicit none
      class(MappedVector3D),intent(inout) :: vector
      type(SEMHex),intent(in) :: geometry
      type(c_ptr),intent(in),optional :: handle
      ! Local
      integer :: i,j,k,iEl,iVar
      real(prec) :: Fx,Fy,Fz
  
      if (present(handle)) then
  
        call ContravariantProjection_MappedVector3D_gpu_wrapper(c_loc(vector % interior), &
                                                                c_loc(geometry % dsdx % interior), &
                                                                vector % interp % N, &
                                                                vector % nVar, &
                                                                vector % nElem)
  
      else
        ! Assume that tensor(j,i) is vector i, component j
        ! => dot product is done along first dimension
        ! to project onto computational space
        do iel = 1,vector % nElem
          do ivar = 1,vector % nVar
            do k = 1,vector % interp % N + 1
              do j = 1,vector % interp % N + 1
                do i = 1,vector % interp % N + 1
    
                  Fx = vector % interior(i,j,k,iEl,iVar,1)
                  Fy = vector % interior(i,j,k,iEl,iVar,2)
                  Fz = vector % interior(i,j,k,iEl,iVar,3)
    
                  vector % interior(i,j,k,iEl,iVar,1) = &
                    geometry % dsdx % interior(i,j,k,iEl,1,1,1)*Fx + &
                    geometry % dsdx % interior(i,j,k,iEl,1,2,1)*Fy + &
                    geometry % dsdx % interior(i,j,k,iEl,1,3,1)*Fz
    
                  vector % interior(i,j,k,iEl,iVar,2) = &
                    geometry % dsdx % interior(i,j,k,iEl,1,1,2)*Fx + &
                    geometry % dsdx % interior(i,j,k,iEl,1,2,2)*Fy + &
                    geometry % dsdx % interior(i,j,k,iEl,1,3,2)*Fz

                  vector % interior(i,j,k,iEl,iVar,3) = &
                    geometry % dsdx % interior(i,j,k,iEl,1,1,3)*Fx + &
                    geometry % dsdx % interior(i,j,k,iEl,1,2,3)*Fy + &
                    geometry % dsdx % interior(i,j,k,iEl,1,3,3)*Fz
    
                end do
              end do
            end do
          end do
        end do
  
      end if
  
    end subroutine ContravariantProjection_MappedVector3D
  
    subroutine Divergence_MappedVector3D(this,geometry,divVector,handle)
      ! Strong Form Operator
      !    !
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(SEMHex),intent(in) :: geometry
      type(MappedScalar3D),intent(inout) :: divVector
      type(c_ptr),intent(inout),optional :: handle
  
      if (present(handle)) then
  
        ! Convert from physical to computational space
        call this % ContravariantProjection(geometry,handle)
  
        ! Compute the divergence
        call this % interp % VectorDivergence_3D(this % interior, &
                                                 divVector % interior, &
                                                 this % nvar, &
                                                 this % nelem, &
                                                 handle)
  
        ! Divide by the jacobian
        call divVector % JacobianWeight(geometry,handle)
  
      else
  
        ! Convert from physical to computational space
        call this % ContravariantProjection(geometry)
  
        ! Compute the divergence
        call this % interp % VectorDivergence_3D(this % interior, &
                                                 divVector % interior, &
                                                 this % nvar, &
                                                 this % nelem)
        ! Divide by the jacobian
        call divVector % JacobianWeight(geometry)
  
      end if
  
    end subroutine Divergence_MappedVector3D
  
    subroutine DGDivergence_MappedVector3D(this,geometry,divVector,handle)
      !! Computes the divergence of a 3-D vector using the weak form
      !! On input, the  attribute of the vector
      !! is assigned and the  attribute is set to the physical
      !! directions of the vector. This method will project the vector
      !! onto the contravariant basis vectors.
      implicit none
      class(MappedVector3D),intent(inout) :: this
      type(SEMHex),intent(in) :: geometry
      type(MappedScalar3D),intent(inout) :: divVector
      type(c_ptr),intent(inout),optional :: handle
  
      if (present(handle)) then
  
        ! Convert from physical to computational space
        call this % ContravariantProjection(geometry,handle)
  
        ! Compute the divergence
        call this % interp % VectorDGDivergence_3D(this % interior, &
                                                   this % boundaryNormal, &
                                                   divVector % interior, &
                                                   this % nvar, &
                                                   this % nelem, &
                                                   handle)
  
        ! Divide by the jacobian
        call divVector % JacobianWeight(geometry,handle)
  
      else
  
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
  
      end if
  
    end subroutine DGDivergence_MappedVector3D

!   FUNCTION Integral_MappedScalar2D(this, geometry, decomp, gpuAccel) RESULT( fRes )
!     !! Calculates the area integral the scalar over all of the geometry.
!     !! Global reduction is done across all MPI ranks when the domain
!     !! decomposition indicates MPI is enabled.
!     IMPLICIT NONE
!     CLASS(MappedScalar2D) :: this
!     TYPE(SEMQuad) :: geometry
!     TYPE(MPILayer) :: decomp
!     LOGICAL :: gpuAccel
!     REAL(prec) :: fRes
!     ! Local
!     INTEGER :: i, j, iEl
!     REAL(prec) :: wi, wj, fint, Jacobian, f

!       IF( gpuAccel ) THEN
!         CALL this % interior % UpdateHost()
!       ENDIF

!       fint = 0.0_prec

!       DO iEl = 1, geometry % x % nElem
!         DO j = 0, geometry % x % interp % N
!           DO i = 0, geometry % x % interp % N

!             ! Coordinate mapping Jacobian
!             Jacobian = geometry % J % interior(i,j,iEl,1)

!             ! Quadrature weights
!             wi = geometry % x % interp % qWeights(i)
!             wj = geometry % x % interp % qWeights(j)

!             f = this % interior(i,j,4,iEl)

!             fint = fint + f*wi*wj*Jacobian

!           ENDDO
!         ENDDO
!       ENDDO

!       CALL decomp % GlobalReduce( fint, fRes )

!   END FUNCTION Integral_MappedScalar2D


end module SELF_MappedData
