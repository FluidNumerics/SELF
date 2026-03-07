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
  use SELF_Geometry

  implicit none

  type,extends(SEMGeometry),public :: SEMHex
    type(Vector3D) :: x ! Physical positions
    type(Tensor3D) :: dxds ! Covariant basis vectors
    type(Tensor3D) :: dsdx ! Contavariant basis vectors
    type(Vector3D) :: nHat ! Normal Vectors pointing across coordinate lines
    type(Scalar3D) :: nScale ! Boundary scale
    type(Scalar3D) :: J ! Jacobian of the transformation

  contains

    procedure,public :: Init => Init_SEMHex
    procedure,public :: Free => Free_SEMHex
    procedure,public :: GenerateFromMesh => GenerateFromMesh_SEMHex
    procedure,public :: CalculateMetricTerms => CalculateMetricTerms_SEMHex
    procedure,private :: CalculateContravariantBasis => CalculateContravariantBasis_SEMHex
    procedure,public :: WriteTecplot => WriteTecplot_SEMHex

  endtype SEMHex

contains

  subroutine Init_SEMHex(this,interp,nElem)
    implicit none
    class(SEMHex),intent(out) :: this
    type(Lagrange),pointer,intent(in) :: interp
    integer,intent(in) :: nElem

    this%nElem = nElem

    call this%x%Init(interp=interp, &
                     nVar=1, &
                     nElem=nElem)

    call this%x%meta(1)%SetName("x")

    call this%dxds%Init(interp=interp, &
                        nVar=1, &
                        nElem=nElem)

    call this%dsdx%Init(interp=interp, &
                        nVar=1, &
                        nElem=nElem)

    call this%nHat%Init(interp=interp, &
                        nVar=1, &
                        nElem=nElem)

    call this%nScale%Init(interp=interp, &
                          nVar=1, &
                          nElem=nElem)

    call this%J%Init(interp=interp, &
                     nVar=1, &
                     nElem=nElem)

  endsubroutine Init_SEMHex

  subroutine Free_SEMHex(this)
    implicit none
    class(SEMHex),intent(inout) :: this

    call this%x%Free()
    call this%dxds%Free()
    call this%dsdx%Free()
    call this%nHat%Free()
    call this%nScale%Free()
    call this%J%Free()

  endsubroutine Free_SEMHex

  subroutine GenerateFromMesh_SEMHex(this,mesh)
    implicit none
    class(SEMHex),intent(inout) :: this
    type(Mesh3D),intent(in) :: mesh
    ! Local
    integer :: iel
    integer :: i,j,k
    type(Lagrange),target :: meshToModel
    type(Vector3D) :: xMesh

    call meshToModel%Init(mesh%nGeo,mesh%quadrature, &
                          this%x%interp%N, &
                          this%x%interp%controlNodeType)

    call xMesh%Init(meshToModel, &
                    1,mesh%nElem)

    ! Set the element internal mesh locations
    do iel = 1,mesh%nElem
      do k = 1,mesh%nGeo+1
        do j = 1,mesh%nGeo+1
          do i = 1,mesh%nGeo+1
            xMesh%interior(i,j,k,iel,1,1:3) = mesh%nodeCoords(1:3,i,j,k,iel)
          enddo
        enddo
      enddo
    enddo

    call xMesh%GridInterp(this%x%interior)
    call this%x%UpdateDevice()
    call this%x%BoundaryInterp()
    call this%x%UpdateHost()
    call this%CalculateMetricTerms()

    call xMesh%Free()
    call meshToModel%Free()

  endsubroutine GenerateFromMesh_SEMHex

  subroutine CalculateContravariantBasis_SEMHex(this)
    implicit none
    class(SEMHex),intent(inout) :: this
    ! Local
    integer :: iEl,i,j,k
    real(prec) :: fac
    real(prec) :: mag
    type(Vector3D) :: xlgradxm,xmgradxl
    type(Vector3D) :: curl_xlgradxm,curl_xmgradxl

    ! Here we use the curl invariant form from Kopriva (2006)
    ! to calculate the contravariant basis vectors
    call xlgradxm%Init(this%x%interp,1,this%x%nElem)
    call xmgradxl%Init(this%x%interp,1,this%x%nElem)

    call curl_xlgradxm%Init(this%x%interp,1,this%x%nElem)
    call curl_xmgradxl%Init(this%x%interp,1,this%x%nElem)

    ! Ja^{1:3}_1 (n=1, m=2, l=3) First component of the contravariant basis vectors
    do iEl = 1,this%nElem
      do k = 1,this%dxds%interp%N+1
        do j = 1,this%dxds%interp%N+1
          do i = 1,this%dxds%interp%N+1
            xlgradxm%interior(i,j,k,iel,1,1:3) = this%x%interior(i,j,k,iel,1,3)*this%dxds%interior(i,j,k,iel,1,2,1:3) ! x(...,l)*dxds(...,m,1:3) ; l=3,m=2
            xmgradxl%interior(i,j,k,iel,1,1:3) = this%x%interior(i,j,k,iel,1,2)*this%dxds%interior(i,j,k,iel,1,3,1:3) ! x(...,m)*dxds(...,l,1:3) ; l=3,m=2
          enddo
        enddo
      enddo
    enddo

    call xlgradxm%Curl(curl_xlgradxm%interior)
    call xmgradxl%Curl(curl_xmgradxl%interior)

    do iEl = 1,this%nElem
      do k = 1,this%dxds%interp%N+1
        do j = 1,this%dxds%interp%N+1
          do i = 1,this%dxds%interp%N+1
            ! In our convention, dsdx(i,j) is contravariant vector j, component i
            ! dsdx(...,n,i) = Ja^{i}_{n} = contravariant vector i, component n;
            ! Here, i = 1:3, and n=1
            this%dsdx%interior(i,j,k,iel,1,1,1:3) = 0.5_prec*( &
                                                    curl_xmgradxl%interior(i,j,k,iel,1,1:3)- &
                                                    curl_xlgradxm%interior(i,j,k,iel,1,1:3))
          enddo
        enddo
      enddo
    enddo

    ! Ja^{1:3}_2 (n=2, m=3, l=1) Second component of the contravariant basis vectors
    do iEl = 1,this%nElem
      do k = 1,this%dxds%interp%N+1
        do j = 1,this%dxds%interp%N+1
          do i = 1,this%dxds%interp%N+1
            xlgradxm%interior(i,j,k,iel,1,1:3) = this%x%interior(i,j,k,iel,1,1)*this%dxds%interior(i,j,k,iel,1,3,1:3) ! x(...,l)*dxds(...,m,1:3) ; l=1,m=3
            xmgradxl%interior(i,j,k,iel,1,1:3) = this%x%interior(i,j,k,iel,1,3)*this%dxds%interior(i,j,k,iel,1,1,1:3) ! x(...,m)*dxds(...,l,1:3) ; l=1,m=3
          enddo
        enddo
      enddo
    enddo

    call xlgradxm%Curl(curl_xlgradxm%interior)
    call xmgradxl%Curl(curl_xmgradxl%interior)

    do iEl = 1,this%nElem
      do k = 1,this%dxds%interp%N+1
        do j = 1,this%dxds%interp%N+1
          do i = 1,this%dxds%interp%N+1
            ! In our convention, dsdx(i,j) is contravariant vector j, component i
            ! dsdx(...,n,i) = Ja^{i}_{n} = contravariant vector i, component n;
            ! Here, i = 1:3, and n=2
            this%dsdx%interior(i,j,k,iel,1,2,1:3) = 0.5_prec*( &
                                                    curl_xmgradxl%interior(i,j,k,iel,1,1:3)- &
                                                    curl_xlgradxm%interior(i,j,k,iel,1,1:3))
          enddo
        enddo
      enddo
    enddo

    ! Ja^{1:3}_3 (n=3, m=1, l=2) Third component of the contravariant basis vectors
    do iEl = 1,this%nElem
      do k = 1,this%dxds%interp%N+1
        do j = 1,this%dxds%interp%N+1
          do i = 1,this%dxds%interp%N+1
            xlgradxm%interior(i,j,k,iel,1,1:3) = this%x%interior(i,j,k,iel,1,2)*this%dxds%interior(i,j,k,iel,1,1,1:3) ! x(...,l)*dxds(...,m,1:3) ; l=2,m=1
            xmgradxl%interior(i,j,k,iel,1,1:3) = this%x%interior(i,j,k,iel,1,1)*this%dxds%interior(i,j,k,iel,1,2,1:3) ! x(...,m)*dxds(...,l,1:3) ; l=2,m=1
          enddo
        enddo
      enddo
    enddo

    call xlgradxm%Curl(curl_xlgradxm%interior)
    call xmgradxl%Curl(curl_xmgradxl%interior)

    do iEl = 1,this%nElem
      do k = 1,this%dxds%interp%N+1
        do j = 1,this%dxds%interp%N+1
          do i = 1,this%dxds%interp%N+1
            ! In our convention, dsdx(i,j) is contravariant vector j, component i
            ! dsdx(...,n,i) = Ja^{i}_{n} = contravariant vector i, component n;
            ! Here, i = 1:3, and n=3
            this%dsdx%interior(i,j,k,iel,1,3,1:3) = 0.5_prec*( &
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
    call this%dsdx%BoundaryInterp() ! Tensor boundary interp is not offloaded

    ! Now, calculate nHat (outward pointing normal)
    do iEl = 1,this%nElem
      do k = 1,6
        do j = 1,this%J%interp%N+1
          do i = 1,this%J%interp%N+1
            if(k == selfSide3D_Top .or. k == selfSide3D_East .or. k == selfSide3D_North) then
              fac = sign(1.0_prec,this%J%boundary(i,j,k,iEl,1))
            else
              fac = -sign(1.0_prec,this%J%boundary(i,j,k,iEl,1))
            endif

            if(k == 1) then ! Bottom

              mag = sqrt(this%dsdx%boundary(i,j,k,iEl,1,1,3)**2+ &
                         this%dsdx%boundary(i,j,k,iEl,1,2,3)**2+ &
                         this%dsdx%boundary(i,j,k,iEl,1,3,3)**2)

              this%nScale%boundary(i,j,k,iEl,1) = mag

              this%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*this%dsdx%boundary(i,j,k,iEl,1,1:3,3)/mag
              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              this%dsdx%boundary(i,j,k,iEl,1,1:3,3) = &
                this%dsdx%boundary(i,j,k,iEl,1,1:3,3)*fac

            elseif(k == 2) then ! South

              mag = sqrt(this%dsdx%boundary(i,j,k,iEl,1,1,2)**2+ &
                         this%dsdx%boundary(i,j,k,iEl,1,2,2)**2+ &
                         this%dsdx%boundary(i,j,k,iEl,1,3,2)**2)

              this%nScale%boundary(i,j,k,iEl,1) = mag

              this%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*this%dsdx%boundary(i,j,k,iEl,1,1:3,2)/mag

              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              this%dsdx%boundary(i,j,k,iEl,1,1:3,2) = &
                this%dsdx%boundary(i,j,k,iEl,1,1:3,2)*fac

            elseif(k == 3) then ! East

              mag = sqrt(this%dsdx%boundary(i,j,k,iEl,1,1,1)**2+ &
                         this%dsdx%boundary(i,j,k,iEl,1,2,1)**2+ &
                         this%dsdx%boundary(i,j,k,iEl,1,3,1)**2)

              this%nScale%boundary(i,j,k,iEl,1) = mag

              this%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*this%dsdx%boundary(i,j,k,iEl,1,1:3,1)/mag
              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              this%dsdx%boundary(i,j,k,iEl,1,1:3,1) = &
                this%dsdx%boundary(i,j,k,iEl,1,1:3,1)*fac

            elseif(k == 4) then ! North

              mag = sqrt(this%dsdx%boundary(i,j,k,iEl,1,1,2)**2+ &
                         this%dsdx%boundary(i,j,k,iEl,1,2,2)**2+ &
                         this%dsdx%boundary(i,j,k,iEl,1,3,2)**2)

              this%nScale%boundary(i,j,k,iEl,1) = mag

              this%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*this%dsdx%boundary(i,j,k,iEl,1,1:3,2)/mag

              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              this%dsdx%boundary(i,j,k,iEl,1,1:3,2) = &
                this%dsdx%boundary(i,j,k,iEl,1,1:3,2)*fac

            elseif(k == 5) then ! West

              mag = sqrt(this%dsdx%boundary(i,j,k,iEl,1,1,1)**2+ &
                         this%dsdx%boundary(i,j,k,iEl,1,2,1)**2+ &
                         this%dsdx%boundary(i,j,k,iEl,1,3,1)**2)

              this%nScale%boundary(i,j,k,iEl,1) = mag

              this%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*this%dsdx%boundary(i,j,k,iEl,1,1:3,1)/mag

              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              this%dsdx%boundary(i,j,k,iEl,1,1:3,1) = &
                this%dsdx%boundary(i,j,k,iEl,1,1:3,1)*fac

            elseif(k == 6) then ! Top

              mag = sqrt(this%dsdx%boundary(i,j,k,iEl,1,1,3)**2+ &
                         this%dsdx%boundary(i,j,k,iEl,1,2,3)**2+ &
                         this%dsdx%boundary(i,j,k,iEl,1,3,3)**2)

              this%nScale%boundary(i,j,k,iEl,1) = mag

              this%nHat%boundary(i,j,k,iEl,1,1:3) = &
                fac*this%dsdx%boundary(i,j,k,iEl,1,1:3,3)/mag

              ! Set the directionality for dsdx on the boundaries
              ! This is primarily used for DG gradient calculations,
              ! which do not use nHat for the boundary terms.
              this%dsdx%boundary(i,j,k,iEl,1,1:3,3) = &
                this%dsdx%boundary(i,j,k,iEl,1,1:3,3)*fac

            endif

          enddo
        enddo
      enddo
    enddo

    call this%dsdx%UpdateDevice()
    call this%nHat%UpdateDevice()
    call this%nScale%UpdateDevice()

  endsubroutine CalculateContravariantBasis_SEMHex

  subroutine CalculateMetricTerms_SEMHex(this)
    implicit none
    class(SEMHex),intent(inout) :: this

    call this%x%Gradient(this%dxds%interior)
    call this%dxds%BoundaryInterp() ! Tensor boundary interp is not offloaded to GPU
    call this%dxds%UpdateDevice()

    call this%dxds%Determinant(this%J%interior)

    call this%J%UpdateDevice()
    call this%J%BoundaryInterp()
    call this%J%UpdateHost()

    call this%CalculateContravariantBasis()

  endsubroutine CalculateMetricTerms_SEMHex

  subroutine WriteTecplot_SEMHex(this,filename)
    implicit none
    class(SEMHex),intent(inout) :: this
    character(*),intent(in) :: filename
    ! Local
    character(8) :: zoneID
    integer :: fUnit
    integer :: iEl,i,j,k,iVar
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
