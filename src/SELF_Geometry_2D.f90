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

module SELF_Geometry_2D

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Data
  use SELF_Scalar_2D
  use SELF_Vector_2D
  use SELF_Tensor_2D
  use SELF_SupportRoutines
  use SELF_Mesh_2D

  implicit none

  type,public :: SEMQuad
    type(Vector2D) :: x ! Physical positions
    type(Tensor2D) :: dxds ! Covariant basis vectors
    type(Tensor2D) :: dsdx ! Contavariant basis vectors
    type(Vector2D) :: nHat ! Normal Vectors pointing across coordinate lines
    type(Scalar2D) :: nScale ! Boundary scale
    type(Scalar2D) :: J ! Jacobian of the transformation
    !! Default-initialized so a freshly allocated object can be distinguished from an
    !! initialized one (the controller uses nElem == 0 to decide Init versus Resize).
    integer :: nElem = 0

    !! Cached scratch for GenerateFromMesh (AMR Stage 6c). Both were previously constructed and
    !! destroyed on every call, which the adaptive loop makes once per epoch:
    !!
    !!   meshToModel - the nGeo -> N interpolant. It depends only on
    !!     (mesh%nGeo, mesh%quadrature, interp%N, interp%controlNodeType), all of which are
    !!     invariant across adaptation, yet building it costs a quadrature plus eight matrices
    !!     and, on GPU builds, eight device allocations and uploads.
    !!   xMesh - staging for the mesh node coordinates. Only its element count changes between
    !!     epochs, so it is resized rather than rebuilt.
    !!
    !! meshToModel is a POINTER so that xMesh may hold a valid interp pointer into it: a
    !! derived-type component cannot carry TARGET, and pointing at a component of an object that
    !! is not itself a target is not conforming, whereas an allocated pointer is always a valid
    !! target. Same reasoning as the storage pools in SELF_DataPool.
    type(Lagrange),pointer :: meshToModel => null()
    type(Vector2D) :: xMesh
    logical :: scratchReady = .false.
    integer :: scratchNGeo = -1
  contains

    procedure,public :: Init => Init_SEMQuad
    procedure,public :: Resize => Resize_SEMQuad
    procedure,public :: Free => Free_SEMQuad
    procedure,public :: GenerateFromMesh => GenerateFromMesh_SEMQuad
    procedure,private :: EnsureScratch => EnsureScratch_SEMQuad
    procedure,public :: CalculateMetricTerms => CalculateMetricTerms_SEMQuad
    procedure,private :: CalculateContravariantBasis => CalculateContravariantBasis_SEMQuad
    procedure,public :: WriteTecplot => WriteTecplot_SEMQuad

  endtype SEMQuad

contains

  subroutine Init_SEMQuad(myGeom,interp,nElem)
    implicit none
    class(SEMQuad),intent(out) :: myGeom
    type(Lagrange),pointer,intent(in) :: interp
    integer,intent(in) :: nElem

    myGeom%nElem = nElem

    call myGeom%x%Init(interp=interp, &
                       nVar=1, &
                       nElem=nElem)

    call myGeom%x%meta(1)%SetName("x")

    call myGeom%dxds%Init(interp=interp, &
                          nVar=1, &
                          nElem=nElem)

    call myGeom%dsdx%Init(interp=interp, &
                          nVar=1, &
                          nElem=nElem)

    call myGeom%nHat%Init(interp=interp, &
                          nVar=1, &
                          nElem=nElem)

    call myGeom%nScale%Init(interp=interp, &
                            nVar=1, &
                            nElem=nElem)

    call myGeom%J%Init(interp=interp, &
                       nVar=1, &
                       nElem=nElem)

  endsubroutine Init_SEMQuad

  subroutine Resize_SEMQuad(myGeom,interp,nElem)
    !! Rebind a live geometry to a new element count, reusing storage where it fits (AMR Stage
    !! 6c). This replaces the Free + Init cycle the adaptive loop performed on a freshly allocated
    !! SEMQuad every epoch, which threw away exactly the amortization Stage 6b introduced: each
    !! member Free released its pools and device buffers, and each Init reallocated, zeroed,
    !! rebuilt metadata and equation parsers, and uploaded the zeros.
    !!
    !! Contents are undefined afterwards; GenerateFromMesh (or the incremental reuse path) fills
    !! them. The cached nGeo -> N scratch is preserved, which is what makes caching it worthwhile.
    implicit none
    class(SEMQuad),intent(inout) :: myGeom
    type(Lagrange),pointer,intent(in) :: interp
    integer,intent(in) :: nElem

    myGeom%nElem = nElem

    call myGeom%x%Resize(interp,1,nElem)
    call myGeom%dxds%Resize(interp,1,nElem)
    call myGeom%dsdx%Resize(interp,1,nElem)
    call myGeom%nHat%Resize(interp,1,nElem)
    call myGeom%nScale%Resize(interp,1,nElem)
    call myGeom%J%Resize(interp,1,nElem)

  endsubroutine Resize_SEMQuad

  subroutine Free_SEMQuad(myGeom)
    implicit none
    class(SEMQuad),intent(inout) :: myGeom

    call myGeom%x%Free()
    call myGeom%dxds%Free()
    call myGeom%dsdx%Free()
    call myGeom%nHat%Free()
    call myGeom%nScale%Free()
    call myGeom%J%Free()

    ! Cached GenerateFromMesh scratch (Stage 6c).
    if(myGeom%scratchReady) then
      call myGeom%xMesh%Free()
      call myGeom%meshToModel%Free()
      deallocate(myGeom%meshToModel)
      myGeom%meshToModel => null()
      myGeom%scratchReady = .false.
      myGeom%scratchNGeo = -1
    endif

  endsubroutine Free_SEMQuad

  subroutine GenerateFromMesh_SEMQuad(myGeom,mesh)
    implicit none
    class(SEMQuad),intent(inout) :: myGeom
    type(Mesh2D),intent(in) :: mesh
    ! Local
    integer :: iel
    integer :: i,j

    call myGeom%EnsureScratch(mesh%nGeo,mesh%quadrature,mesh%nElem)

    ! Set the element internal mesh locations
    do iel = 1,mesh%nElem
      do j = 1,mesh%nGeo+1
        do i = 1,mesh%nGeo+1
          myGeom%xMesh%interior(i,j,iel,1,1:2) = mesh%nodeCoords(1:2,i,j,iel)
        enddo
      enddo
    enddo

    call myGeom%xMesh%GridInterp(myGeom%x%interior)
    call myGeom%x%UpdateDevice()
    call myGeom%x%BoundaryInterp() ! Boundary interp will run on GPU if enabled, hence why we close in update host/device
    call myGeom%x%UpdateHost()
    call myGeom%CalculateMetricTerms()

  endsubroutine GenerateFromMesh_SEMQuad

  subroutine EnsureScratch_SEMQuad(myGeom,nGeo,quadrature,nElem)
    !! Prepare the cached GenerateFromMesh scratch for a mesh with nGeo/quadrature and nElem
    !! elements (AMR Stage 6c). The nGeo -> N interpolant is built once and reused; the node
    !! coordinate staging is resized, so an adapting run stops rebuilding either one per epoch.
    !!
    !! The interpolant is rebuilt only if nGeo or the quadrature actually changes, which does not
    !! happen across adaptation but is handled so that reusing one SEMQuad against a different
    !! mesh family stays correct.
    implicit none
    class(SEMQuad),intent(inout) :: myGeom
    integer,intent(in) :: nGeo
    integer,intent(in) :: quadrature
    integer,intent(in) :: nElem

    if(myGeom%scratchReady) then
      if(myGeom%scratchNGeo /= nGeo .or. myGeom%meshToModel%controlNodeType /= quadrature) then
        call myGeom%xMesh%Free()
        call myGeom%meshToModel%Free()
        deallocate(myGeom%meshToModel)
        myGeom%meshToModel => null()
        myGeom%scratchReady = .false.
      endif
    endif

    if(.not. myGeom%scratchReady) then
      allocate(myGeom%meshToModel)
      call myGeom%meshToModel%Init(nGeo, &
                                   quadrature, &
                                   myGeom%x%interp%N, &
                                   myGeom%x%interp%controlNodeType)
      call myGeom%xMesh%Init(myGeom%meshToModel,1,nElem)
      myGeom%scratchReady = .true.
      myGeom%scratchNGeo = nGeo
    elseif(myGeom%xMesh%nElem /= nElem) then
      call myGeom%xMesh%Resize(myGeom%meshToModel,1,nElem)
    endif

  endsubroutine EnsureScratch_SEMQuad

  subroutine CalculateContravariantBasis_SEMQuad(myGeom)
    implicit none
    class(SEMQuad),intent(inout) :: myGeom
    ! Local
    integer :: iEl,i,j,k
    real(prec) :: fac
    real(prec) :: mag

    ! Now calculate the contravariant basis vectors
    ! In this convention, dsdx(j,i) is contravariant vector i, component j
    ! To project onto contravariant vector i, dot vector along the first dimension
    do iEl = 1,myGeom%nElem
      do j = 1,myGeom%dxds%interp%N+1
        do i = 1,myGeom%dxds%interp%N+1

          myGeom%dsdx%interior(i,j,iel,1,1,1) = myGeom%dxds%interior(i,j,iel,1,2,2)
          myGeom%dsdx%interior(i,j,iel,1,2,1) = -myGeom%dxds%interior(i,j,iel,1,1,2)
          myGeom%dsdx%interior(i,j,iel,1,1,2) = -myGeom%dxds%interior(i,j,iel,1,2,1)
          myGeom%dsdx%interior(i,j,iel,1,2,2) = myGeom%dxds%interior(i,j,iel,1,1,1)

        enddo
      enddo
    enddo

    ! Interpolate the contravariant tensor to the boundaries
    call myGeom%dsdx%BoundaryInterp() ! Tensor boundary interp is not offloaded

    ! Now, modify the sign of dsdx so that
    ! myGeom % dsdx % boundary is equal to the outward pointing normal vector
    do iEl = 1,myGeom%nElem
      do k = 1,4
        do i = 1,myGeom%J%interp%N+1
          if(k == selfSide2D_East .or. k == selfSide2D_North) then
            fac = sign(1.0_prec,myGeom%J%boundary(i,k,iEl,1))
          else
            fac = -sign(1.0_prec,myGeom%J%boundary(i,k,iEl,1))
          endif

          if(k == 1) then ! South

            mag = sqrt(myGeom%dsdx%boundary(i,k,iEl,1,1,2)**2+ &
                       myGeom%dsdx%boundary(i,k,iEl,1,2,2)**2)

            myGeom%nScale%boundary(i,k,iEl,1) = mag

            myGeom%nHat%boundary(i,k,iEl,1,1:2) = &
              fac*myGeom%dsdx%boundary(i,k,iEl,1,1:2,2)/mag

          elseif(k == 2) then ! East

            mag = sqrt(myGeom%dsdx%boundary(i,k,iEl,1,1,1)**2+ &
                       myGeom%dsdx%boundary(i,k,iEl,1,2,1)**2)

            myGeom%nScale%boundary(i,k,iEl,1) = mag

            myGeom%nHat%boundary(i,k,iEl,1,1:2) = &
              fac*myGeom%dsdx%boundary(i,k,iEl,1,1:2,1)/mag

          elseif(k == 3) then ! North

            mag = sqrt(myGeom%dsdx%boundary(i,k,iEl,1,1,2)**2+ &
                       myGeom%dsdx%boundary(i,k,iEl,1,2,2)**2)

            myGeom%nScale%boundary(i,k,iEl,1) = mag

            myGeom%nHat%boundary(i,k,iEl,1,1:2) = &
              fac*myGeom%dsdx%boundary(i,k,iEl,1,1:2,2)/mag

          elseif(k == 4) then ! West

            mag = sqrt(myGeom%dsdx%boundary(i,k,iEl,1,1,1)**2+ &
                       myGeom%dsdx%boundary(i,k,iEl,1,2,1)**2)

            myGeom%nScale%boundary(i,k,iEl,1) = mag

            myGeom%nHat%boundary(i,k,iEl,1,1:2) = &
              fac*myGeom%dsdx%boundary(i,k,iEl,1,1:2,1)/mag

          endif

          ! Set the directionality for dsdx on the boundaries
          myGeom%dsdx%boundary(i,k,iEl,1,1:2,1:2) = &
            myGeom%dsdx%boundary(i,k,iEl,1,1:2,1:2)*fac

        enddo
      enddo
    enddo

    call myGeom%dsdx%UpdateDevice()
    call myGeom%nHat%UpdateDevice()
    call myGeom%nScale%UpdateDevice()

  endsubroutine CalculateContravariantBasis_SEMQuad

  subroutine CalculateMetricTerms_SEMQuad(myGeom)
    implicit none
    class(SEMQuad),intent(inout) :: myGeom

    call myGeom%x%Gradient(myGeom%dxds%interior)
    ! No boundary interpolation of dxds, and no device upload of it (AMR Stage 6c). In 2-D dxds
    ! is scratch consumed only inside this module: J comes from its interior via Determinant, and
    ! dsdx from its interior via the adjugate in CalculateContravariantBasis, which then fills
    ! dsdx%boundary with dsdx's own BoundaryInterp. A whole-tree search for %dxds%boundary and
    ! for any dxds device pointer finds exactly one hit, in SELF_Geometry_1D (the SEMLine type,
    ! which is distinct), and none respectively. The removed host tensor boundary interpolation
    ! had no GPU override and ran over the largest arrays in SEMQuad.
    !
    ! Note the contrast with x, just above: x%boundary has 37 consumers across SELF_Points,
    ! ESAtmo2D, examples and tests, so its boundary interpolation and the device round trip that
    ! wraps it are load-bearing and stay.

    call myGeom%dxds%Determinant(myGeom%J%interior)

    call myGeom%J%UpdateDevice()
    call myGeom%J%BoundaryInterp()
    call myGeom%J%UpdateHost()

    call myGeom%CalculateContravariantBasis()

  endsubroutine CalculateMetricTerms_SEMQuad

  subroutine WriteTecplot_SEMQuad(this,filename)
    implicit none
    class(SEMQuad),intent(inout) :: this
    character(*),intent(in) :: filename
    ! Local
    character(8) :: zoneID
    integer :: fUnit
    integer :: iEl,i,j
    character(LEN=self_TecplotHeaderLength) :: tecHeader
    character(LEN=self_FormatLength) :: fmat

    open(UNIT=NEWUNIT(fUnit), &
         FILE=trim(filename), &
         FORM='formatted', &
         STATUS='replace')

    tecHeader = 'VARIABLES = "X", "Y", "eID"'

    write(fUnit,*) trim(tecHeader)

    ! Create format statement
    write(fmat,*) 3
    fmat = '('//trim(fmat)//'(ES16.7E3,1x))'

    do iEl = 1,this%x%nElem

      ! TO DO :: Get the global element ID
      write(zoneID,'(I8.8)') iEl
      write(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',this%x%interp%N+1, &
        ', J=',this%x%interp%N+1

      do j = 1,this%x%interp%N+1
        do i = 1,this%x%interp%N+1

          write(fUnit,fmat) this%x%interior(i,j,iEl,1,1), &
            this%x%interior(i,j,iEl,1,2),real(iEl,prec)

        enddo
      enddo

    enddo

    close(UNIT=fUnit)

  endsubroutine WriteTecplot_SEMQuad

endmodule SELF_Geometry_2D
