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

module SELF_Geometry_3D

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Data
  use SELF_Scalar_3D
  use SELF_Vector_3D
  use SELF_Tensor_3D
  use SELF_SupportRoutines
  use SELF_Mesh_3D

  implicit none

  type,public :: SEMHex
    type(Vector3D) :: x ! Physical positions
    type(Tensor3D) :: dxds ! Covariant basis vectors
    type(Tensor3D) :: dsdx ! Contavariant basis vectors
    type(Vector3D) :: nHat ! Normal Vectors pointing across coordinate lines
    type(Scalar3D) :: nScale ! Boundary scale
    type(Scalar3D) :: J ! Jacobian of the transformation
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
    type(Vector3D) :: xMesh
    logical :: scratchReady = .false.
    integer :: scratchNGeo = -1
  contains

    procedure,public :: Init => Init_SEMHex
    procedure,public :: Resize => Resize_SEMHex
    procedure,public :: Free => Free_SEMHex
    procedure,public :: GenerateFromMesh => GenerateFromMesh_SEMHex
    procedure,public :: GenerateFromNodeCoords => GenerateFromNodeCoords_SEMHex
    procedure,public :: CopyElements => CopyElements_SEMHex
    procedure,public :: UploadGeometry => UploadGeometry_SEMHex
    procedure,private :: EnsureScratch => EnsureScratch_SEMHex
    procedure,public :: CalculateMetricTerms => CalculateMetricTerms_SEMHex
    procedure,private :: CalculateContravariantBasis => CalculateContravariantBasis_SEMHex
    procedure,public :: WriteTecplot => WriteTecplot_SEMHex

  endtype SEMHex

contains

  subroutine Init_SEMHex(myGeom,interp,nElem)
    implicit none
    class(SEMHex),intent(out) :: myGeom
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

  endsubroutine Init_SEMHex

  subroutine Resize_SEMHex(myGeom,interp,nElem)
    !! Rebind a live geometry to a new element count, reusing storage where it fits (AMR Stage
    !! 6c). This replaces the Free + Init cycle the adaptive loop performed on a freshly allocated
    !! SEMHex every epoch, which threw away exactly the amortization Stage 6b introduced: each
    !! member Free released its pools and device buffers, and each Init reallocated, zeroed,
    !! rebuilt metadata and equation parsers, and uploaded the zeros.
    !!
    !! Contents are undefined afterwards; GenerateFromMesh (or the incremental reuse path) fills
    !! them. The cached nGeo -> N scratch is preserved, which is what makes caching it worthwhile.
    implicit none
    class(SEMHex),intent(inout) :: myGeom
    type(Lagrange),pointer,intent(in) :: interp
    integer,intent(in) :: nElem

    myGeom%nElem = nElem

    call myGeom%x%Resize(interp,1,nElem)
    call myGeom%dxds%Resize(interp,1,nElem)
    call myGeom%dsdx%Resize(interp,1,nElem)
    call myGeom%nHat%Resize(interp,1,nElem)
    call myGeom%nScale%Resize(interp,1,nElem)
    call myGeom%J%Resize(interp,1,nElem)

  endsubroutine Resize_SEMHex

  subroutine CopyElements_SEMHex(myGeom,src,srcIdx,dstIdx,n)
    !! Copy whole-element geometry blocks from src into myGeom: element srcIdx(k) of src becomes
    !! element dstIdx(k) of myGeom, for k = 1..n (AMR Stage 6c).
    !!
    !! This is exact, not an interpolation, and it is what lets an adaptation epoch skip
    !! regenerating the elements it did not change. Every geometry quantity for an element depends
    !! only on that element's own mesh node coordinates - GenerateFromMesh, CalculateMetricTerms
    !! and CalculateContravariantBasis contain no neighbour coupling, no side pairing and no
    !! reduction, and the ±sign convention for normals is element-local - so moving an element's
    !! block between two geometries preserves it exactly.
    !!
    !! Element is dimension 4 of every array, so each element's data is contiguous within a given
    !! set of trailing indices; the whole-slice assignments below are the natural expression of
    !! that and let the compiler emit block copies.
    implicit none
    class(SEMHex),intent(inout) :: myGeom
    type(SEMHex),intent(in) :: src
    integer,intent(in) :: srcIdx(:)
    integer,intent(in) :: dstIdx(:)
    integer,intent(in) :: n
    ! Local
    integer :: k,s,d

    do k = 1,n
      s = srcIdx(k)
      d = dstIdx(k)

      myGeom%x%interior(:,:,:,d,:,:) = src%x%interior(:,:,:,s,:,:)
      myGeom%x%boundary(:,:,:,d,:,:) = src%x%boundary(:,:,:,s,:,:)

      myGeom%dxds%interior(:,:,:,d,:,:,:) = src%dxds%interior(:,:,:,s,:,:,:)
      myGeom%dxds%boundary(:,:,:,d,:,:,:) = src%dxds%boundary(:,:,:,s,:,:,:)

      myGeom%dsdx%interior(:,:,:,d,:,:,:) = src%dsdx%interior(:,:,:,s,:,:,:)
      myGeom%dsdx%boundary(:,:,:,d,:,:,:) = src%dsdx%boundary(:,:,:,s,:,:,:)

      myGeom%nHat%interior(:,:,:,d,:,:) = src%nHat%interior(:,:,:,s,:,:)
      myGeom%nHat%boundary(:,:,:,d,:,:) = src%nHat%boundary(:,:,:,s,:,:)

      myGeom%nScale%interior(:,:,:,d,:) = src%nScale%interior(:,:,:,s,:)
      myGeom%nScale%boundary(:,:,:,d,:) = src%nScale%boundary(:,:,:,s,:)

      myGeom%J%interior(:,:,:,d,:) = src%J%interior(:,:,:,s,:)
      myGeom%J%boundary(:,:,:,d,:) = src%J%boundary(:,:,:,s,:)
    enddo

  endsubroutine CopyElements_SEMHex

  subroutine UploadGeometry_SEMHex(myGeom)
    !! Push the geometry the solver kernels read to the device. Mirrors the uploads that
    !! GenerateFromMesh's own path performs, for use when geometry was assembled by element copy
    !! rather than generated (AMR Stage 6c).
    !!
    !! dxds is included, unlike the 2-D UploadGeometry: the 3-D CalculateMetricTerms uploads it
    !! after its boundary interpolation, so an assembled geometry mirrors that exactly.
    implicit none
    class(SEMHex),intent(inout) :: myGeom

    call myGeom%x%UpdateDevice()
    call myGeom%dxds%UpdateDevice()
    call myGeom%dsdx%UpdateDevice()
    call myGeom%nHat%UpdateDevice()
    call myGeom%nScale%UpdateDevice()
    call myGeom%J%UpdateDevice()

  endsubroutine UploadGeometry_SEMHex

  subroutine Free_SEMHex(myGeom)
    implicit none
    class(SEMHex),intent(inout) :: myGeom

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

  endsubroutine Free_SEMHex

  subroutine GenerateFromMesh_SEMHex(myGeom,mesh)
    implicit none
    class(SEMHex),intent(inout) :: myGeom
    type(Mesh3D),intent(in) :: mesh

    call myGeom%GenerateFromNodeCoords(mesh%nodeCoords,mesh%nGeo,mesh%quadrature,mesh%nElem)
    call myGeom%x%UpdateDevice()
    call myGeom%x%BoundaryInterp() ! Boundary interp will run on GPU if enabled, hence why we close in update host/device
    call myGeom%x%UpdateHost()
    call myGeom%CalculateMetricTerms()

  endsubroutine GenerateFromMesh_SEMHex

  subroutine GenerateFromNodeCoords_SEMHex(myGeom,nodeCoords,nGeo,quadrature,nElem)
    !! Generate geometry for nElem elements directly from their mesh node coordinates (AMR Stage
    !! 6c). GenerateFromMesh is a thin wrapper over this.
    !!
    !! Taking the coordinates rather than a Mesh3D is what allows the adaptive loop to generate a
    !! COMPACTED set of elements - just the ones an epoch actually changed - without teaching the
    !! shared data classes about element subsets. The generation loops still run over every element
    !! they are given; the saving comes from being given fewer.
    implicit none
    class(SEMHex),intent(inout) :: myGeom
    real(prec),intent(in) :: nodeCoords(1:3,1:nGeo+1,1:nGeo+1,1:nGeo+1,1:nElem)
    integer,intent(in) :: nGeo
    integer,intent(in) :: quadrature
    integer,intent(in) :: nElem
    ! Local
    integer :: iel,i,j,k

    if(nElem <= 0) return

    call myGeom%EnsureScratch(nGeo,quadrature,nElem)

    ! Set the element internal mesh locations
    do iel = 1,nElem
      do k = 1,nGeo+1
        do j = 1,nGeo+1
          do i = 1,nGeo+1
            myGeom%xMesh%interior(i,j,k,iel,1,1:3) = nodeCoords(1:3,i,j,k,iel)
          enddo
        enddo
      enddo
    enddo

    call myGeom%xMesh%GridInterp(myGeom%x%interior)

  endsubroutine GenerateFromNodeCoords_SEMHex

  subroutine EnsureScratch_SEMHex(myGeom,nGeo,quadrature,nElem)
    !! Prepare the cached GenerateFromMesh scratch for a mesh with nGeo/quadrature and nElem
    !! elements (AMR Stage 6c). The nGeo -> N interpolant is built once and reused; the node
    !! coordinate staging is resized, so an adapting run stops rebuilding either one per epoch.
    !!
    !! The interpolant is rebuilt only if nGeo or the quadrature actually changes, which does not
    !! happen across adaptation but is handled so that reusing one SEMHex against a different
    !! mesh family stays correct.
    implicit none
    class(SEMHex),intent(inout) :: myGeom
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

  endsubroutine EnsureScratch_SEMHex

  subroutine CalculateContravariantBasis_SEMHex(myGeom)
    implicit none
    class(SEMHex),intent(inout) :: myGeom
    ! Local
    integer :: iEl,i,j,k
    real(prec) :: fac
    real(prec) :: mag
    type(Vector3D) :: xlgradxm,xmgradxl
    type(Vector3D) :: curl_xlgradxm,curl_xmgradxl

    ! Here we use the curl invariant form from Kopriva (2006)
    ! to calculate the contravariant basis vectors
    call xlgradxm%Init(myGeom%x%interp,1,myGeom%x%nElem)
    call xmgradxl%Init(myGeom%x%interp,1,myGeom%x%nElem)

    call curl_xlgradxm%Init(myGeom%x%interp,1,myGeom%x%nElem)
    call curl_xmgradxl%Init(myGeom%x%interp,1,myGeom%x%nElem)

    ! Ja^{1:3}_1 (n=1, m=2, l=3) First component of the contravariant basis vectors
    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1
            xlgradxm%interior(i,j,k,iel,1,1:3) = myGeom%x%interior(i,j,k,iel,1,3)*myGeom%dxds%interior(i,j,k,iel,1,2,1:3) ! x(...,l)*dxds(...,m,1:3) ; l=3,m=2
            xmgradxl%interior(i,j,k,iel,1,1:3) = myGeom%x%interior(i,j,k,iel,1,2)*myGeom%dxds%interior(i,j,k,iel,1,3,1:3) ! x(...,m)*dxds(...,l,1:3) ; l=3,m=2
          enddo
        enddo
      enddo
    enddo

    call xlgradxm%Curl(curl_xlgradxm%interior)
    call xmgradxl%Curl(curl_xmgradxl%interior)

    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1
            ! In our convention, dsdx(i,j) is contravariant vector j, component i
            ! dsdx(...,n,i) = Ja^{i}_{n} = contravariant vector i, component n;
            ! Here, i = 1:3, and n=1
            myGeom%dsdx%interior(i,j,k,iel,1,1,1:3) = 0.5_prec*( &
                                                      curl_xmgradxl%interior(i,j,k,iel,1,1:3)- &
                                                      curl_xlgradxm%interior(i,j,k,iel,1,1:3))
          enddo
        enddo
      enddo
    enddo

    ! Ja^{1:3}_2 (n=2, m=3, l=1) Second component of the contravariant basis vectors
    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1
            xlgradxm%interior(i,j,k,iel,1,1:3) = myGeom%x%interior(i,j,k,iel,1,1)*myGeom%dxds%interior(i,j,k,iel,1,3,1:3) ! x(...,l)*dxds(...,m,1:3) ; l=1,m=3
            xmgradxl%interior(i,j,k,iel,1,1:3) = myGeom%x%interior(i,j,k,iel,1,3)*myGeom%dxds%interior(i,j,k,iel,1,1,1:3) ! x(...,m)*dxds(...,l,1:3) ; l=1,m=3
          enddo
        enddo
      enddo
    enddo

    call xlgradxm%Curl(curl_xlgradxm%interior)
    call xmgradxl%Curl(curl_xmgradxl%interior)

    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1
            ! In our convention, dsdx(i,j) is contravariant vector j, component i
            ! dsdx(...,n,i) = Ja^{i}_{n} = contravariant vector i, component n;
            ! Here, i = 1:3, and n=2
            myGeom%dsdx%interior(i,j,k,iel,1,2,1:3) = 0.5_prec*( &
                                                      curl_xmgradxl%interior(i,j,k,iel,1,1:3)- &
                                                      curl_xlgradxm%interior(i,j,k,iel,1,1:3))
          enddo
        enddo
      enddo
    enddo

    ! Ja^{1:3}_3 (n=3, m=1, l=2) Third component of the contravariant basis vectors
    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1
            xlgradxm%interior(i,j,k,iel,1,1:3) = myGeom%x%interior(i,j,k,iel,1,2)*myGeom%dxds%interior(i,j,k,iel,1,1,1:3) ! x(...,l)*dxds(...,m,1:3) ; l=2,m=1
            xmgradxl%interior(i,j,k,iel,1,1:3) = myGeom%x%interior(i,j,k,iel,1,1)*myGeom%dxds%interior(i,j,k,iel,1,2,1:3) ! x(...,m)*dxds(...,l,1:3) ; l=2,m=1
          enddo
        enddo
      enddo
    enddo

    call xlgradxm%Curl(curl_xlgradxm%interior)
    call xmgradxl%Curl(curl_xmgradxl%interior)

    do iEl = 1,myGeom%nElem
      do k = 1,myGeom%dxds%interp%N+1
        do j = 1,myGeom%dxds%interp%N+1
          do i = 1,myGeom%dxds%interp%N+1
            ! In our convention, dsdx(i,j) is contravariant vector j, component i
            ! dsdx(...,n,i) = Ja^{i}_{n} = contravariant vector i, component n;
            ! Here, i = 1:3, and n=3
            myGeom%dsdx%interior(i,j,k,iel,1,3,1:3) = 0.5_prec*( &
                                                      curl_xmgradxl%interior(i,j,k,iel,1,1:3)- &
                                                      curl_xlgradxm%interior(i,j,k,iel,1,1:3))
          enddo
        enddo
      enddo
    enddo

    call xlgradxm%Free()
    call xmgradxl%Free()
    call curl_xlgradxm%Free()
    call curl_xmgradxl%Free()

    ! Interpolate the contravariant tensor to the boundaries
    call myGeom%dsdx%BoundaryInterp() ! Tensor boundary interp is not offloaded

    ! Now, calculate nHat (outward pointing normal)
    do iEl = 1,myGeom%nElem
      do k = 1,6
        do j = 1,myGeom%J%interp%N+1
          do i = 1,myGeom%J%interp%N+1
            if(k == selfSide3D_Top .or. k == selfSide3D_East .or. k == selfSide3D_North) then
              fac = sign(1.0_prec,myGeom%J%boundary(i,j,k,iEl,1))
            else
              fac = -sign(1.0_prec,myGeom%J%boundary(i,j,k,iEl,1))
            endif

            if(k == 1) then ! Bottom

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,3)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,3)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,3)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3)/mag
              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3) = &
                myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3)*fac

            elseif(k == 2) then ! South

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,2)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,2)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,2)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2)/mag

              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2) = &
                myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2)*fac

            elseif(k == 3) then ! East

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,1)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,1)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,1)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1)/mag
              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1) = &
                myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1)*fac

            elseif(k == 4) then ! North

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,2)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,2)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,2)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2)/mag

              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2) = &
                myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,2)*fac

            elseif(k == 5) then ! West

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,1)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,1)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,1)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1)/mag

              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1) = &
                myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,1)*fac

            elseif(k == 6) then ! Top

              mag = sqrt(myGeom%dsdx%boundary(i,j,k,iEl,1,1,3)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,2,3)**2+ &
                         myGeom%dsdx%boundary(i,j,k,iEl,1,3,3)**2)

              myGeom%nScale%boundary(i,j,k,iEl,1) = mag

              myGeom%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3)/mag

              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3) = &
                myGeom%dsdx%boundary(i,j,k,iEl,1,1:3,3)*fac

            endif

          enddo
        enddo
      enddo
    enddo

    call myGeom%dsdx%UpdateDevice()
    call myGeom%nHat%UpdateDevice()
    call myGeom%nScale%UpdateDevice()

  endsubroutine CalculateContravariantBasis_SEMHex

  subroutine CalculateMetricTerms_SEMHex(myGeom)
    implicit none
    class(SEMHex),intent(inout) :: myGeom

    call myGeom%x%Gradient(myGeom%dxds%interior)
    call myGeom%dxds%BoundaryInterp() ! Tensor boundary interp is not offloaded to GPU
    call myGeom%dxds%UpdateDevice()

    call myGeom%dxds%Determinant(myGeom%J%interior)

    call myGeom%J%UpdateDevice()
    call myGeom%J%BoundaryInterp()
    call myGeom%J%UpdateHost()

    call myGeom%CalculateContravariantBasis()

  endsubroutine CalculateMetricTerms_SEMHex

  subroutine WriteTecplot_SEMHex(this,filename)
    implicit none
    class(SEMHex),intent(inout) :: this
    character(*),intent(in) :: filename
    ! Local
    character(8) :: zoneID
    integer :: fUnit
    integer :: iEl,i,j,k
    character(LEN=self_TecplotHeaderLength) :: tecHeader
    character(LEN=self_FormatLength) :: fmat

    open(UNIT=NEWUNIT(fUnit), &
         FILE=trim(filename), &
         FORM='formatted', &
         STATUS='replace')

    tecHeader = 'VARIABLES = "X", "Y", "Z", "eID"'

    write(fUnit,*) trim(tecHeader)

    ! Create format statement
    write(fmat,*) 4
    fmat = '('//trim(fmat)//'(ES16.7E3,1x))'

    do iEl = 1,this%x%nElem

      ! TO DO :: Get the global element ID
      write(zoneID,'(I8.8)') iEl
      write(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',this%x%interp%N+1, &
        ', J=',this%x%interp%N+1,', K=',this%x%interp%N+1

      do k = 1,this%x%interp%N+1
        do j = 1,this%x%interp%N+1
          do i = 1,this%x%interp%N+1

            write(fUnit,fmat) this%x%interior(i,j,k,iEl,1,1), &
              this%x%interior(i,j,k,iEl,1,2),this%x%interior(i,j,k,iel,1,3),real(iEl,prec)

          enddo
        enddo
      enddo

    enddo

    close(UNIT=fUnit)

  endsubroutine WriteTecplot_SEMHex

endmodule SELF_Geometry_3D
