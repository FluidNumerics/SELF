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
  use SELF_Geometry

  implicit none

  type,extends(SEMGeometry),public :: SEMQuad
    type(Vector2D) :: x ! Physical positions
    type(Tensor2D) :: dxds ! Covariant basis vectors
    type(Tensor2D) :: dsdx ! Contavariant basis vectors
    type(Vector2D) :: nHat ! Normal Vectors pointing across coordinate lines
    type(Scalar2D) :: nScale ! Boundary scale
    type(Scalar2D) :: J ! Jacobian of the transformation
  contains

    procedure,public :: Init => Init_SEMQuad
    procedure,public :: Free => Free_SEMQuad
    procedure,public :: GenerateFromMesh => GenerateFromMesh_SEMQuad
    procedure,public :: CalculateMetricTerms => CalculateMetricTerms_SEMQuad
    procedure,private :: CalculateContravariantBasis => CalculateContravariantBasis_SEMQuad
    procedure,public :: WriteTecplot => WriteTecplot_SEMQuad

  endtype SEMQuad

contains

  subroutine Init_SEMQuad(this,interp,nElem)
    implicit none
    class(SEMQuad),intent(out) :: this
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

  endsubroutine Init_SEMQuad

  subroutine Free_SEMQuad(this)
    implicit none
    class(SEMQuad),intent(inout) :: this

    call this%x%Free()
    call this%dxds%Free()
    call this%dsdx%Free()
    call this%nHat%Free()
    call this%nScale%Free()
    call this%J%Free()

  endsubroutine Free_SEMQuad

  subroutine GenerateFromMesh_SEMQuad(this,mesh)
    implicit none
    class(SEMQuad),intent(inout) :: this
    type(Mesh2D),intent(in) :: mesh
    ! Local
    integer :: iel
    integer :: i,j
    type(Lagrange),target :: meshToModel
    type(Vector2D) :: xMesh

    call meshToModel%Init(mesh%nGeo, &
                          mesh%quadrature, &
                          this%x%interp%N, &
                          this%x%interp%controlNodeType)

    call xMesh%Init(meshToModel,1,mesh%nElem)

    ! Set the element internal mesh locations
    do iel = 1,mesh%nElem
      do j = 1,mesh%nGeo+1
        do i = 1,mesh%nGeo+1
          xMesh%interior(i,j,iel,1,1:2) = mesh%nodeCoords(1:2,i,j,iel)
        enddo
      enddo
    enddo

    call xMesh%GridInterp(this%x%interior)
    call this%x%UpdateDevice()
    call this%x%BoundaryInterp() ! Boundary interp will run on GPU if enabled, hence why we close in update host/device
    call this%x%UpdateHost()
    call this%CalculateMetricTerms()

    call xMesh%Free()
    call meshToModel%Free()

  endsubroutine GenerateFromMesh_SEMQuad

  subroutine CalculateContravariantBasis_SEMQuad(this)
    implicit none
    class(SEMQuad),intent(inout) :: this
    ! Local
    integer :: iEl,i,j,k
    real(prec) :: fac
    real(prec) :: mag

    ! Now calculate the contravariant basis vectors
    ! In this convention, dsdx(j,i) is contravariant vector i, component j
    ! To project onto contravariant vector i, dot vector along the first dimension
    do iEl = 1,this%nElem
      do j = 1,this%dxds%interp%N+1
        do i = 1,this%dxds%interp%N+1

          this%dsdx%interior(i,j,iel,1,1,1) = this%dxds%interior(i,j,iel,1,2,2)
          this%dsdx%interior(i,j,iel,1,2,1) = -this%dxds%interior(i,j,iel,1,1,2)
          this%dsdx%interior(i,j,iel,1,1,2) = -this%dxds%interior(i,j,iel,1,2,1)
          this%dsdx%interior(i,j,iel,1,2,2) = this%dxds%interior(i,j,iel,1,1,1)

        enddo
      enddo
    enddo

    ! Interpolate the contravariant tensor to the boundaries
    call this%dsdx%BoundaryInterp() ! Tensor boundary interp is not offloaded

    ! Now, modify the sign of dsdx so that
    ! this % dsdx % boundary is equal to the outward pointing normal vector
    do iEl = 1,this%nElem
      do k = 1,4
        do i = 1,this%J%interp%N+1
          if(k == selfSide2D_East .or. k == selfSide2D_North) then
            fac = sign(1.0_prec,this%J%boundary(i,k,iEl,1))
          else
            fac = -sign(1.0_prec,this%J%boundary(i,k,iEl,1))
          endif

          if(k == 1) then ! South

            mag = sqrt(this%dsdx%boundary(i,k,iEl,1,1,2)**2+ &
                       this%dsdx%boundary(i,k,iEl,1,2,2)**2)

            this%nScale%boundary(i,k,iEl,1) = mag

            this%nHat%boundary(i,k,iEl,1,1:2) = &
              fac*this%dsdx%boundary(i,k,iEl,1,1:2,2)/mag

          elseif(k == 2) then ! East

            mag = sqrt(this%dsdx%boundary(i,k,iEl,1,1,1)**2+ &
                       this%dsdx%boundary(i,k,iEl,1,2,1)**2)

            this%nScale%boundary(i,k,iEl,1) = mag

            this%nHat%boundary(i,k,iEl,1,1:2) = &
              fac*this%dsdx%boundary(i,k,iEl,1,1:2,1)/mag

          elseif(k == 3) then ! North

            mag = sqrt(this%dsdx%boundary(i,k,iEl,1,1,2)**2+ &
                       this%dsdx%boundary(i,k,iEl,1,2,2)**2)

            this%nScale%boundary(i,k,iEl,1) = mag

            this%nHat%boundary(i,k,iEl,1,1:2) = &
              fac*this%dsdx%boundary(i,k,iEl,1,1:2,2)/mag

          elseif(k == 4) then ! West

            mag = sqrt(this%dsdx%boundary(i,k,iEl,1,1,1)**2+ &
                       this%dsdx%boundary(i,k,iEl,1,2,1)**2)

            this%nScale%boundary(i,k,iEl,1) = mag

            this%nHat%boundary(i,k,iEl,1,1:2) = &
              fac*this%dsdx%boundary(i,k,iEl,1,1:2,1)/mag

          endif

          ! Set the directionality for dsdx on the boundaries
          this%dsdx%boundary(i,k,iEl,1,1:2,1:2) = &
            this%dsdx%boundary(i,k,iEl,1,1:2,1:2)*fac

        enddo
      enddo
    enddo

    call this%dsdx%UpdateDevice()
    call this%nHat%UpdateDevice()
    call this%nScale%UpdateDevice()

  endsubroutine CalculateContravariantBasis_SEMQuad

  subroutine CalculateMetricTerms_SEMQuad(this)
    implicit none
    class(SEMQuad),intent(inout) :: this

    call this%x%Gradient(this%dxds%interior)
    call this%dxds%BoundaryInterp() ! Tensor boundary interp is not offloaded to GPU
    call this%dxds%UpdateDevice()

    call this%dxds%Determinant(this%J%interior)

    call this%J%UpdateDevice()
    call this%J%BoundaryInterp()
    call this%J%UpdateHost()

    call this%CalculateContravariantBasis()

  endsubroutine CalculateMetricTerms_SEMQuad

  subroutine WriteTecplot_SEMQuad(this,filename)
    implicit none
    class(SEMQuad),intent(inout) :: this
    character(*),intent(in) :: filename
    ! Local
    character(8) :: zoneID
    integer :: fUnit
    integer :: iEl,i,j,iVar
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
