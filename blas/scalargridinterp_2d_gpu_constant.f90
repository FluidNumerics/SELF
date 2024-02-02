program blas_program
  implicit NONE
  INTEGER::exit_code
  
  exit_code=scalargridinterp_2d_gpu_constant()
  print*,exit_code
  exit_code=blas_scalargridinterp_2d_gpu_constant()
  print*,exit_code
  
      CONTAINS

integer function scalargridinterp_2d_gpu_constant() result(r)
  use SELF_Constants
  use SELF_Memory
  use SELF_Lagrange
  use SELF_Data

  implicit none

  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  integer,parameter :: nvar = 1
  integer,parameter :: nelem = 1000
#ifdef DOUBLE_PRECISION
  real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
  real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
  type(Scalar2D) :: f
  type(Scalar2D) :: fTarget
  type(Lagrange),target :: interp
  type(Lagrange),target :: interpTarget

  ! Create an interpolant
  call interp % Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

  call interpTarget % Init(N=targetDegree, &
                           controlNodeType=UNIFORM, &
                           M=targetDegree, &
                           targetNodeType=UNIFORM)

  ! Initialize scalars
  call f % Init(interp,nvar,nelem)
  call fTarget % Init(interpTarget,nvar,nelem)

  ! Set the source scalar (on the control grid) to a non-zero constant
  f % interior % hostdata = 1.0_prec

  call f % interior % updatedevice()

  ! Interpolate with gpuAccel = .true.
  call f % GridInterp(fTarget,.true.)

  call fTarget % interior % updatehost()

  ! Calculate diff from exact
  fTarget % interior % hostdata = abs(fTarget % interior % hostdata - 1.0_prec)

  if (maxval(fTarget % interior % hostdata) <= tolerance) then
    r = 0
  else
    r = 1
  end if

  call f % free()
  call fTarget % free()
  call interp % free()
  call interpTarget % free()

end function scalargridinterp_2d_gpu_constant

integer function blas_scalargridinterp_2d_gpu_constant() result(r)
  use SELF_Constants
  use SELF_Memory
  use SELF_Lagrange
  use SELF_Data
  use hipfort_hipblas
  use self_hip
  use self_hip_enums

  implicit none

  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  integer,parameter :: nvar = 1
  integer,parameter :: nelem = 1000
#ifdef DOUBLE_PRECISION
  real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
  real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
  type(Scalar2D) :: f
  type(Scalar2D) :: fTarget
  TYPE(hfReal_r4) :: fIntermediate
  real(prec), pointer :: fIntermediate_(:,:,:,:)
  real(prec), pointer :: fTarget_(:,:,:,:)
  type(Lagrange),target :: interp
  type(Lagrange),target :: interpTarget
  type(c_ptr) :: handle
  integer(c_int) :: m
  integer(c_int) :: n
  integer(c_int) :: k
  real(c_prec) :: alpha
  integer(c_int) :: lda
  integer(c_int) :: ldb
  integer(c_int) :: ldc
  real(c_prec) :: beta
  ! for gemvstridedbatch
  integer :: i
  integer(c_int64_t) :: strideA
  integer(c_int) :: incx
  integer(c_int64_t) :: stridex
  integer(c_int) :: incy
  integer(c_int64_t) :: stridey
  integer(c_int) :: batchCount
  integer(kind(HIPBLAS_STATUS_SUCCESS)) :: status


  call hipblasCheck(hipblasCreate(handle))

  ! Create an interpolant
  call interp % Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

  call interpTarget % Init(N=targetDegree, &
                           controlNodeType=UNIFORM, &
                           M=targetDegree, &
                           targetNodeType=UNIFORM)

  ! Initialize scalars
  call f % Init(interp,nvar,nelem)
  call fTarget % Init(interpTarget,nvar,nelem)
  call fIntermediate % Alloc(loBound=(/0,0,1,1/), &
                              upBound=(/interp % M,interp % N,nvar,nelem/))
  call c_f_pointer(fIntermediate % deviceData, fIntermediate_, [targetDegree+1,controlDegree+1,nvar,nelem])
  call c_f_pointer(fTarget % interior % deviceData, fTarget_, [targetDegree+1,targetDegree+1,nvar,nelem])
  ! print*, lbound(fIntermediate_)
  ! print*, ubound(fIntermediate_)
  ! print*, lbound(fTarget_)
  ! print*, ubound(fTarget_)
  ! Set the source scalar (on the control grid) to a non-zero constant
  f % interior % hostdata = 1.0_prec

  call f % interior % updatedevice()

  ! Interpolate with gpuAccel = .true.
  !call f % GridInterp(fTarget,.true.)
! Interpolate with gpuAccel = .true.
  m = targetDegree+1 ! number of rows of A^T
  n = nvar*nelem*(controlDegree+1) ! number of columns of B
  k = controlDegree+1! number of columns of A^T
  alpha = 1.0_c_prec
  lda = k ! leading dimension of A (interoplation matrix) (ControlDegree+1)
  ldb = k ! leading dimension of B (f) (controlDegree+1)
  ldc = m ! leading dimension of C (fTarget) (targetDegree+1)
  beta = 0.0_c_prec

! do i = 1,1000
#ifdef DOUBLE_PRECISION
    ! First pass interpolates in the first quadrature dimension
    status = hipblasDgemm(handle,HIPBLAS_OP_T, HIPBLAS_OP_N, &
      m, n, k, alpha, interp % iMatrix % deviceData, lda, &
      f % interior % deviceData, ldb, beta, &
      fIntermediate % deviceData, ldc)
#else
    status = hipblasSgemm(handle,HIPBLAS_OP_T, HIPBLAS_OP_N, &
      m, n, k, alpha, interp % iMatrix % deviceData, lda, &
      f % interior % deviceData, ldb, beta, &
      fIntermediate % deviceData, ldc)
#endif

m = controlDegree+1 ! number of rows of A
n = targetDegree+1 ! number of columns of A
alpha = 1.0_c_prec
lda = m ! leading dimension of A
strideA = 0 ! stride for the batches of A (no stride)
incx = targetDegree+1 !
stridex = (controlDegree+1)*(targetDegree+1)
beta = 0.0_c_prec
incy = targetDegree+1
stridey = (targetDegree+1)*(targetDegree+1)
batchCount = nvar*nelem
do i = 0, targetDegree
#ifdef DOUBLE_PRECISION
    status = hipblasDgemvStridedBatched(handle,HIPBLAS_OP_T, m, &
      n,alpha,interp % iMatrix % deviceData, lda, strideA, &
      c_loc(fIntermediate_(1+i,1,1,1)), incx, stridex, beta, &
      c_loc(fTarget_(1+i,1,1,1)), incy, stridey, batchCount)
#else
    status = hipblasSgemvStridedBatched(handle,HIPBLAS_OP_T, m, &
      n,alpha,interp % iMatrix % deviceData, lda, strideA, &
      c_loc(fIntermediate_(1+i,1,1,1)), incx, stridex, beta, &
      c_loc(fTarget_(1+i,1,1,1)), incy, stridey, batchCount)
#endif
enddo

! enddo

  call fTarget % interior % updatehost()
  ! print*, fTarget % interior % hostdata(:,:,1,1)
  ! Calculate diff from exact
  fTarget % interior % hostdata = abs(fTarget % interior % hostdata - 1.0_prec)

  if (maxval(fTarget % interior % hostdata) <= tolerance) then
    r = 0
  else
    r = 1
  end if

  call f % free()
  call fTarget % free()
  call interp % free()
  call interpTarget % free()
  call fIntermediate % free()

end function blas_scalargridinterp_2d_gpu_constant

end program blas_program
