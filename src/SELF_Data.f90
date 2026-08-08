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

module SELF_Data

  use SELF_Constants
  use SELF_Lagrange
  use SELF_Metadata
  use FEQParse
  use SELF_HDF5

  use HDF5
  use iso_c_binding

  implicit none

  type,public :: SELF_DataObj
  !! The SELF_DataObj class is a base class for all data objects in SELF.
  !! A data object in SELF is a multidimensional array of data, represented
  !! on both host and device, that is associated with an interpolant, metadata,
  !! and (optionally) an equation string.
  !! Type extensions of the SELF_DataObj include scalars, vectors, and tensors
  !! in 1-D, 2-D, and 3-D using the storage patterns that are expected for
  !! derivative and interpolation operations defined in SELF_Lagrange.f90
  !! Additionally, each extended type has the necessary attributes to store
  !! information on element interiors and element boundaries, both of which
  !! are commonly used for spectral element solvers.

    integer :: nVar
    integer :: nElem
    integer :: N
    integer :: M
    type(Lagrange),pointer :: interp
    type(Metadata),allocatable :: meta(:)
    type(EquationParser),allocatable :: eqn(:)

    !! Identity element list (/1,2,...,nElem/), the default argument for the element-subset
    !! form of the tensor-product operators. Operators that accept an optional elems(:)
    !! argument loop over a resolved index list; when the caller passes nothing they loop
    !! over this list instead, which reproduces the plain 1:nElem traversal exactly (same
    !! trip count, same order, same arithmetic) and so is bitwise identical to the
    !! subset-free code. Kept as a member rather than rebuilt per call so that operators
    !! declaring intent(in) on the object can still use it and no allocation happens in a
    !! hot path. Built by MapArrays (i.e. by both Init and Resize).
    integer,pointer,contiguous :: allElem(:) => null()

    !! Identity mortar list (/1,2,...,nMortars/), the same idea for the operators whose
    !! subset is a list of 2:1 nonconforming interfaces rather than of elements. Unlike
    !! allElem this cannot be built by MapArrays, because the mortar count belongs to the
    !! mesh and not to the field; it is built on first use next to the mortar staging
    !! buffer, and released here so that a pointer component does not outlive the object.
    integer,pointer,contiguous :: allMortar(:) => null()

  contains

    ! Procedures for setting metadata for
    procedure,public :: SetName => SetName_DataObj
    procedure,public :: SetDescription => SetDescription_DataObj
    procedure,public :: SetUnits => SetUnits_DataObj
    generic,public :: SetEquation => SetEquation_DataObj
    procedure,private :: SetEquation_DataObj

  endtype SELF_DataObj

  integer,parameter :: selfStrongForm = 0
  integer,parameter :: selfWeakDGForm = 1
  integer,parameter :: selfWeakCGForm = 2
  integer,parameter :: selfWeakBRForm = 3

contains

! -- Subset helpers -- !
!
! The subset-aware form of the tensor-product operators replaces a `do ... = 1, n` loop over
! every element (or every mortar) with a loop over an explicit index list. Local time
! stepping uses this to evaluate a tendency on one refinement level at a time; every
! pre-existing caller passes no list at all and gets the full traversal, unchanged.

  subroutine EnsureIndexList(allIdx,n)
    !! Make allIdx the identity list (/1,...,n/), reallocating only when the length is wrong.
    !! Called from each data type's MapArrays for the element list, and lazily for the mortar
    !! list (whose length is a property of the mesh, not of the field), so that no allocation
    !! happens in a hot path.
    implicit none
    integer,pointer,contiguous,intent(inout) :: allIdx(:)
    integer,intent(in) :: n
    ! Local
    integer :: i

    if(associated(allIdx)) then
      if(size(allIdx) /= n) then
        deallocate(allIdx)
        allIdx => null()
      endif
    endif

    if(.not. associated(allIdx)) then
      allocate(allIdx(1:n))
      do i = 1,n
        allIdx(i) = i
      enddo
    endif

  endsubroutine EnsureIndexList

  subroutine ResolveIndexList(allIdx,n,subset,idx,ns)
    !! Resolve the optional subset argument of a subset-aware operator into the index list
    !! idx(1:ns) that its loops traverse.
    !!
    !! When subset is present and associated the operator touches only those entries. When it
    !! is absent - or present but disassociated - the identity list is used, so `idx(i) == i`
    !! and `ns == n`: the loop is the original 1:n loop, with the same trip count, order and
    !! arithmetic, hence bitwise unchanged. Treating a null pointer as "no subset" lets a
    !! caller thread an optional subset straight through without branching on it.
    !!
    !! Indices are RANK-LOCAL for elements and GLOBAL (replicated) for mortars, matching the
    !! arrays they index. No bounds check is made; these are inner-loop paths.
    !!
    !! subset is a POINTER dummy rather than a TARGET dummy on purpose: associating idx with
    !! the target of a pointer is well defined for as long as that target lives, whereas
    !! associating it with a plain TARGET dummy would leave idx undefined the moment this
    !! routine returns. Callers therefore hold their lists as pointers (see the local time
    !! stepping schedule), and callers that want everything simply omit the argument.
    implicit none
    integer,pointer,contiguous,intent(in) :: allIdx(:)
    integer,intent(in) :: n
    integer,pointer,contiguous,intent(in),optional :: subset(:)
    integer,pointer,contiguous,intent(out) :: idx(:)
    integer,intent(out) :: ns

    if(present(subset)) then
      if(associated(subset)) then
        idx => subset
        ns = size(subset)
        return
      endif
    endif

    idx => allIdx
    ns = n

  endsubroutine ResolveIndexList

  subroutine RejectSubset(subset,file,line)
    !! Guard for backends that have no subset kernel. The GPU operators flatten the element
    !! dimension into a single degree-of-freedom count, so restricting them needs an index
    !! argument threaded through the kernels; until that exists, silently ignoring the subset
    !! would produce a wrong answer, so it is an error. Passing no subset is always accepted,
    !! which is every pre-existing call site.
    implicit none
    integer,pointer,contiguous,intent(in),optional :: subset(:)
    character(*),intent(in) :: file
    integer,intent(in) :: line

    if(present(subset)) then
      if(associated(subset)) then
        print*,file,':',line, &
          ' : Error : subset operators (local time stepping) are not implemented '// &
          'for this backend. Build for CPU to use local time stepping.'
        stop 1
      endif
    endif

  endsubroutine RejectSubset

! -- DataObj -- !

  subroutine SetName_DataObj(this,ivar,name)
    !! Set the name of the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: name

    call this%meta(ivar)%SetName(name)

  endsubroutine SetName_DataObj

  subroutine SetDescription_DataObj(this,ivar,description)
    !! Set the description of the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: description

    call this%meta(ivar)%SetDescription(description)

  endsubroutine SetDescription_DataObj

  subroutine SetUnits_DataObj(this,ivar,units)
    !! Set the units of the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: units

    call this%meta(ivar)%SetUnits(units)

  endsubroutine SetUnits_DataObj

  subroutine SetEquation_DataObj(this,ivar,eqnChar)
    !! Sets the equation parser for the `ivar-th` variable
    implicit none
    class(SELF_DataObj),intent(inout) :: this
    integer,intent(in) :: ivar
    character(*),intent(in) :: eqnChar

    this%eqn(ivar) = EquationParser(trim(eqnChar), &
                                    (/'x','y','z','t'/))

  endsubroutine SetEquation_DataObj

endmodule SELF_Data
