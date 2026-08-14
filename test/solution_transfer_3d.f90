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

program test

  implicit none
  integer :: exit_code

  exit_code = solution_transfer_3d()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains

  integer function solution_transfer_3d() result(r)
    !! Validates AMR Stage 3 solution transfer against the Stage 2 geometry, end to end:
    !!
    !!   1. Prolongation exactness: for a trilinear field (degree <= N in each direction, so its
    !!      degree-N interpolant is the field itself), the prolonged children reproduce the
    !!      parent polynomial pointwise at the fine mesh's (mapped) nodes to roundoff.
    !!   2. Reversibility: prolonging a coarse-element solution to its eight children and then
    !!      restricting it back reproduces the coarse solution to roundoff (sum_k P_k R_k = I).
    !!   3. Reference-cell conservation: sum_ijk w_i w_j w_k u_parent(i,j,k) equals
    !!      (1/8) sum_c sum_ijk w_i w_j w_k u_child_c(i,j,k) to roundoff (each child covers an
    !!      eighth of the parent reference cube).
    !!   4. Physical conservation: the cell-integrated quantity int u dV (with the geometry
    !!      Jacobian) is identical on the coarse mesh and on the uniformly refined mesh whose
    !!      children carry the prolonged solution. Because Stage-2 refinement is exact and the
    !!      coarse solution is its own degree-N interpolant, the two integrals agree to roundoff.
    !!
    !! The child-element ordering used here (fine element 8*(p-1)+c is child c of coarse element
    !! p) is the one produced by UniformRefineMesh.
    use SELF_Constants
    use SELF_Lagrange
    use SELF_Mesh_3D
    use SELF_Geometry_3D
    use SELF_MeshRefinement_3D
    use SELF_SolutionTransfer_3D

    implicit none

    integer,parameter :: controlDegree = 4
    integer,parameter :: nvar = 1
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 1.0e-10_prec
#else
    real(prec),parameter :: tolerance = 1.0e-3_prec
#endif
    type(Lagrange),target :: interp
    type(Mesh3D),target :: coarseMesh,fineMesh
    type(SEMHex),target :: coarseGeom,fineGeom
    integer :: bcids(1:6)
    integer :: p,c,i,j,k,e,nc,nf
    real(prec) :: x,y,z,expected
    real(prec),allocatable :: uCoarse(:,:,:,:,:),uFine(:,:,:,:,:),uBack(:,:,:,:,:)
    real(prec),allocatable :: uChildren(:,:,:,:,:),uParent(:,:,:,:)
    real(prec) :: intCoarse,intFine,idErr,exErr,refErr,sParent,sChildren

    r = 0
    bcids(1:6) = [SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED, &
                  SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED]

    call coarseMesh%StructuredMesh(2,2,2,1,1,1,0.5_prec,0.5_prec,0.5_prec,bcids)
    call UniformRefineMesh(coarseMesh,fineMesh)

    call interp%Init(N=controlDegree,controlNodeType=GAUSS,M=controlDegree,targetNodeType=UNIFORM)
    call coarseGeom%Init(interp,coarseMesh%nElem)
    call coarseGeom%GenerateFromMesh(coarseMesh)
    call fineGeom%Init(interp,fineMesh%nElem)
    call fineGeom%GenerateFromMesh(fineMesh)

    nc = coarseMesh%nElem
    nf = fineMesh%nElem
    allocate(uCoarse(1:controlDegree+1,1:controlDegree+1,1:controlDegree+1,1:nc,1:nvar))
    allocate(uBack(1:controlDegree+1,1:controlDegree+1,1:controlDegree+1,1:nc,1:nvar))
    allocate(uFine(1:controlDegree+1,1:controlDegree+1,1:controlDegree+1,1:nf,1:nvar))
    allocate(uChildren(1:controlDegree+1,1:controlDegree+1,1:controlDegree+1,1:nvar,1:8))
    allocate(uParent(1:controlDegree+1,1:controlDegree+1,1:controlDegree+1,1:nvar))

    ! ---- 1. Prolongation exactness on a trilinear field ----
    ! A trilinear field is its own degree-N interpolant on every (affine) element, so the
    ! prolonged children must reproduce it pointwise at the fine mesh's mapped nodes.
    do p = 1,nc
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            x = coarseGeom%x%interior(i,j,k,p,1,1)
            y = coarseGeom%x%interior(i,j,k,p,1,2)
            z = coarseGeom%x%interior(i,j,k,p,1,3)
            uCoarse(i,j,k,p,1) = 2.0_prec+3.0_prec*x-5.0_prec*y+4.0_prec*z+ &
                                 0.25_prec*x*y-0.5_prec*y*z+x*z+0.75_prec*x*y*z
          enddo
        enddo
      enddo
    enddo
    do p = 1,nc
      call ProlongToChildren(interp,nvar,uCoarse(:,:,:,p,:),uChildren)
      do c = 1,8
        uFine(:,:,:,8*(p-1)+c,1) = uChildren(:,:,:,1,c)
      enddo
    enddo
    exErr = 0.0_prec
    do e = 1,nf
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            x = fineGeom%x%interior(i,j,k,e,1,1)
            y = fineGeom%x%interior(i,j,k,e,1,2)
            z = fineGeom%x%interior(i,j,k,e,1,3)
            expected = 2.0_prec+3.0_prec*x-5.0_prec*y+4.0_prec*z+ &
                       0.25_prec*x*y-0.5_prec*y*z+x*z+0.75_prec*x*y*z
            exErr = max(exErr,abs(uFine(i,j,k,e,1)-expected))
          enddo
        enddo
      enddo
    enddo
    print*,"trilinear-field prolongation max error :",exErr
    if(exErr > tolerance) then
      print*,"FAIL: prolongation is not exact for a trilinear field"
      r = 1
    endif

    ! Coarse solution: a smooth field sampled at the coarse nodes (its degree-N interpolant).
    do p = 1,nc
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            x = coarseGeom%x%interior(i,j,k,p,1,1)
            y = coarseGeom%x%interior(i,j,k,p,1,2)
            z = coarseGeom%x%interior(i,j,k,p,1,3)
            uCoarse(i,j,k,p,1) = sin(1.7_prec*x)*cos(1.3_prec*y)*cos(1.1_prec*z)+ &
                                 0.5_prec*x*y*z
          enddo
        enddo
      enddo
    enddo

    ! Prolong each coarse element onto its eight fine children (fine elem 8*(p-1)+c).
    do p = 1,nc
      call ProlongToChildren(interp,nvar,uCoarse(:,:,:,p,:),uChildren)
      do c = 1,8
        uFine(:,:,:,8*(p-1)+c,1) = uChildren(:,:,:,1,c)
      enddo
    enddo

    ! ---- 2. Restrict back and check reversibility ----
    idErr = 0.0_prec
    do p = 1,nc
      do c = 1,8
        uChildren(:,:,:,1,c) = uFine(:,:,:,8*(p-1)+c,1)
      enddo
      call RestrictFromChildren(interp,nvar,uChildren,uParent)
      uBack(:,:,:,p,:) = uParent(:,:,:,:)
      idErr = max(idErr,maxval(abs(uBack(:,:,:,p,1)-uCoarse(:,:,:,p,1))))
    enddo
    print*,"prolong-restrict reversibility max error :",idErr
    if(idErr > tolerance) then
      print*,"FAIL: prolong-restrict not reversible"
      r = 1
    endif

    ! ---- 3. Reference-cell conservation per parent ----
    refErr = 0.0_prec
    do p = 1,nc
      sParent = 0.0_prec
      sChildren = 0.0_prec
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            sParent = sParent+uCoarse(i,j,k,p,1)* &
                      interp%qWeights(i)*interp%qWeights(j)*interp%qWeights(k)
            do c = 1,8
              sChildren = sChildren+uFine(i,j,k,8*(p-1)+c,1)* &
                          interp%qWeights(i)*interp%qWeights(j)*interp%qWeights(k)
            enddo
          enddo
        enddo
      enddo
      refErr = max(refErr,abs(sParent-0.125_prec*sChildren))
    enddo
    print*,"reference-cell conservation max error :",refErr
    if(refErr > tolerance) then
      print*,"FAIL: reference-cell integrals are not conserved"
      r = 1
    endif

    ! ---- 4. Physical conservation: int u dV on coarse vs fine ----
    intCoarse = 0.0_prec
    do e = 1,nc
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            intCoarse = intCoarse+uCoarse(i,j,k,e,1)*coarseGeom%J%interior(i,j,k,e,1)* &
                        interp%qWeights(i)*interp%qWeights(j)*interp%qWeights(k)
          enddo
        enddo
      enddo
    enddo
    intFine = 0.0_prec
    do e = 1,nf
      do k = 1,controlDegree+1
        do j = 1,controlDegree+1
          do i = 1,controlDegree+1
            intFine = intFine+uFine(i,j,k,e,1)*fineGeom%J%interior(i,j,k,e,1)* &
                      interp%qWeights(i)*interp%qWeights(j)*interp%qWeights(k)
          enddo
        enddo
      enddo
    enddo
    print*,"integral coarse, fine :",intCoarse,intFine
    if(abs(intCoarse-intFine) > tolerance*max(abs(intCoarse),1.0_prec)) then
      print*,"FAIL: prolongation is not conservative"
      r = 1
    endif

    if(r == 0) print*,"SOLUTION TRANSFER CHECKS PASSED"

    deallocate(uCoarse,uBack,uFine,uChildren,uParent)
    call coarseGeom%Free()
    call fineGeom%Free()
    call interp%Free()
    call fineMesh%Free()
    call coarseMesh%Free()

  endfunction solution_transfer_3d

endprogram test
