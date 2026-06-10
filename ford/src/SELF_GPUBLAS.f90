module SELF_GPUBLAS

  use SELF_GPU_enums
  use SELF_Constants
  use iso_c_binding
  use iso_fortran_env

  implicit none

  interface hipblasCreate
#ifdef HAVE_CUDA
    function hipblasCreate_(handle) bind(c,name="cublasCreate_v2")
#else
      function hipblasCreate_(handle) bind(c,name="hipblasCreate")
#endif
        use iso_c_binding
        use SELF_GPU_enums
        implicit none
        integer(kind(HIPBLAS_STATUS_SUCCESS)) :: hipblasCreate_
        type(c_ptr) :: handle
      endfunction
      endinterface

      interface hipblasDestroy
#ifdef HAVE_CUDA
        function hipblasDestroy_(handle) bind(c,name="cublasDestroy_v2")
#else
          function hipblasDestroy_(handle) bind(c,name="hipblasDestroy")
#endif
            use iso_c_binding
            use SELF_GPU_enums
            implicit none
            integer(kind(HIPBLAS_STATUS_SUCCESS)) :: hipblasDestroy_
            type(c_ptr),value :: handle
          endfunction
          endinterface

          interface hipblasSgemm
#ifdef HAVE_CUDA
            function hipblasSgemm_(handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) bind(c,name="cublasSgemm_v2")
#else
              function hipblasSgemm_(handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) bind(c,name="hipblasSgemm")
#endif
                use iso_c_binding
                use SELF_GPU_enums
                implicit none
                integer(kind(HIPBLAS_STATUS_SUCCESS)) :: hipblasSgemm_
                type(c_ptr),value :: handle
                integer(kind(HIPBLAS_OP_N)),value :: transa
                integer(kind(HIPBLAS_OP_N)),value :: transb
                integer(c_int),value :: m
                integer(c_int),value :: n
                integer(c_int),value :: k
                real(c_float) :: alpha
                type(c_ptr),value :: A
                integer(c_int),value :: lda
                type(c_ptr),value :: B
                integer(c_int),value :: ldb
                real(c_float) :: beta
                type(c_ptr),value :: C
                integer(c_int),value :: ldc
              endfunction
              endinterface

              interface hipblasDgemm
#ifdef HAVE_CUDA
                function hipblasDgemm_(handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) bind(c,name="cublasDgemm_v2")
#else
                  function hipblasDgemm_(handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) bind(c,name="hipblasDgemm")
#endif
                    use iso_c_binding
                    use SELF_GPU_enums
                    implicit none
                    integer(kind(HIPBLAS_STATUS_SUCCESS)) :: hipblasDgemm_
                    type(c_ptr),value :: handle
                    integer(kind(HIPBLAS_OP_N)),value :: transa
                    integer(kind(HIPBLAS_OP_N)),value :: transb
                    integer(c_int),value :: m
                    integer(c_int),value :: n
                    integer(c_int),value :: k
                    real(c_double) :: alpha
                    type(c_ptr),value :: A
                    integer(c_int),value :: lda
                    type(c_ptr),value :: B
                    integer(c_int),value :: ldb
                    real(c_double) :: beta
                    type(c_ptr),value :: C
                    integer(c_int),value :: ldc
                  endfunction
                  endinterface

                  interface hipblasSgemvStridedBatched
#ifdef HAVE_CUDA
                    function hipblasSgemvStridedBatched_(handle,trans,m,n,alpha,A,lda,strideA,x, &
                                                         incx,stridex,beta,y,incy,stridey,batchCount) &
                      bind(c,name="cublasSgemvStridedBatched")
#else
                      function hipblasSgemvStridedBatched_(handle,trans,m,n,alpha,A,lda,strideA,x, &
                                                           incx,stridex,beta,y,incy,stridey,batchCount) &
                        bind(c,name="hipblasSgemvStridedBatched")
#endif
                        use iso_c_binding
                        use SELF_GPU_enums
                        implicit none
                        integer(kind(HIPBLAS_STATUS_SUCCESS)) :: hipblasSgemvStridedBatched_
                        type(c_ptr),value :: handle
                        integer(kind(HIPBLAS_OP_N)),value :: trans
                        integer(c_int),value :: m
                        integer(c_int),value :: n
                        real(c_float) :: alpha
                        type(c_ptr),value :: A
                        integer(c_int),value :: lda
                        integer(c_int64_t),value :: strideA
                        type(c_ptr),value :: x
                        integer(c_int),value :: incx
                        integer(c_int64_t),value :: stridex
                        real(c_float) :: beta
                        type(c_ptr),value :: y
                        integer(c_int),value :: incy
                        integer(c_int64_t),value :: stridey
                        integer(c_int),value :: batchCount
                      endfunction

                      endinterface

                      interface hipblasDgemvStridedBatched
#ifdef HAVE_CUDA
                        function hipblasDgemvStridedBatched_(handle,trans,m,n,alpha,A,lda,strideA,x, &
                                                             incx,stridex,beta,y,incy,stridey,batchCount) &
                          bind(c,name="cublasDgemvStridedBatched")
#else
                          function hipblasDgemvStridedBatched_(handle,trans,m,n,alpha,A,lda,strideA,x, &
                                                               incx,stridex,beta,y,incy,stridey,batchCount) &
                            bind(c,name="hipblasDgemvStridedBatched")
#endif
                            use iso_c_binding
                            use SELF_GPU_enums
                            implicit none
                            integer(kind(HIPBLAS_STATUS_SUCCESS)) :: hipblasDgemvStridedBatched_
                            type(c_ptr),value :: handle
                            integer(kind(HIPBLAS_OP_N)),value :: trans
                            integer(c_int),value :: m
                            integer(c_int),value :: n
                            real(c_double) :: alpha
                            type(c_ptr),value :: A
                            integer(c_int),value :: lda
                            integer(c_int64_t),value :: strideA
                            type(c_ptr),value :: x
                            integer(c_int),value :: incx
                            integer(c_int64_t),value :: stridex
                            real(c_double) :: beta
                            type(c_ptr),value :: y
                            integer(c_int),value :: incy
                            integer(c_int64_t),value :: stridey
                            integer(c_int),value :: batchCount
                          endfunction
                          endinterface

                          contains

#ifdef DOUBLE_PRECISION
#define hipblasgemm hipblasDgemm
#define hipblasgemvStridedBatched hipblasDgemvStridedBatched
#else
#define hipblasgemm hipblasSgemm
#define hipblasgemvStridedBatched hipblasSgemvStridedBatched
#endif

                          subroutine hipblasCheck(hipblasError_t)
                            use SELF_GPU_enums
                            implicit none
                            integer(kind(HIPBLAS_STATUS_SUCCESS)) :: hipblasError_t

                            if(hipblasError_t /= HIPBLAS_STATUS_SUCCESS) then
                              write(*,*) "GPUBLAS ERROR: Error code = ",hipblasError_t
                              call exit(hipblasError_t)
                            endif
                          endsubroutine hipblasCheck

                          subroutine self_blas_matrixop_1d(A,f,Af,opArows,opAcols,bcols,handle)
                            type(c_ptr),intent(in) :: A
                            real(prec),pointer,intent(in) :: f(:,:,:)
                            real(prec),pointer,intent(inout) :: Af(:,:,:)
                            integer,intent(in) :: opArows,opAcols,bcols
                            type(c_ptr),intent(inout) :: handle
                            ! Local
                            integer(c_int) :: m
                            integer(c_int) :: n
                            integer(c_int) :: k
                            real(c_prec) :: alpha
                            integer(c_int) :: lda
                            integer(c_int) :: ldb
                            integer(c_int) :: ldc
                            real(c_prec) :: beta

                            m = opArows ! number of rows of A^T
                            n = bcols ! number of columns of B
                            k = opAcols ! number of columns of A^T
                            alpha = 1.0_c_prec
                            lda = k ! leading dimension of A (matrix)
                            ldb = k ! leading dimension of B (f)
                            ldc = m ! leading dimension of C (Af)
                            beta = 0.0_c_prec

                            call hipblasCheck(hipblasgemm(handle, &
                                                          HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                                          m,n,k,alpha, &
                                                          A,lda, &
                                                          c_loc(f),ldb, &
                                                          beta, &
                                                          c_loc(Af),ldc))

                          endsubroutine self_blas_matrixop_1d

                          subroutine self_blas_matrixop_dim1_2d(A,f,Af,controldegree,targetdegree,nvars,nelems,handle)
                            type(c_ptr),intent(in) :: A
                            real(prec),pointer,intent(in) :: f(:,:,:,:)
                            real(prec),pointer,intent(inout) :: Af(:,:,:,:)
                            integer,intent(in) :: controldegree,targetdegree,nvars,nelems
                            type(c_ptr),intent(inout) :: handle
                            ! Local
                            integer(c_int) :: m
                            integer(c_int) :: n
                            integer(c_int) :: k
                            real(c_prec) :: alpha
                            integer(c_int) :: lda
                            integer(c_int) :: ldb
                            integer(c_int) :: ldc
                            real(c_prec) :: beta

                            m = targetdegree+1 ! number of rows of A^T
                            n = nvars*nelems*(controldegree+1) ! number of columns of B
                            k = controldegree+1 ! number of columns of A^T
                            alpha = 1.0_c_prec
                            lda = k ! leading dimension of A (interpolation/derivative matrix)
                            ldb = k ! leading dimension of B (f)
                            ldc = m ! leading dimension of C (fTarget)
                            beta = 0.0_c_prec

                            ! First pass interpolates in the first quadrature dimension
                            call hipblasCheck(hipblasgemm(handle, &
                                                          HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                                          m,n,k,alpha, &
                                                          A,lda, &
                                                          c_loc(f),ldb,beta, &
                                                          c_loc(Af),ldc))

                          endsubroutine self_blas_matrixop_dim1_2d

                          subroutine self_blas_matrixop_dim2_2d(A,f,Af,beta,controldegree,targetdegree,nvars,nelems,handle)
                            type(c_ptr),intent(in) :: A
                            real(prec),pointer,intent(in) :: f(:,:,:,:)
                            real(prec),pointer,intent(inout) :: Af(:,:,:,:)
                            real(c_prec),intent(in) :: beta
                            integer,intent(in) :: controldegree,targetdegree,nvars,nelems
                            type(c_ptr),intent(inout) :: handle
                            ! Local
                            integer(c_int) :: m
                            integer(c_int) :: n
                            real(c_prec) :: alpha
                            integer(c_int) :: lda

                            integer :: i
                            integer(c_int64_t) :: strideA
                            integer(c_int) :: incx
                            integer(c_int64_t) :: stridex
                            integer(c_int) :: incy
                            integer(c_int64_t) :: stridey
                            integer(c_int) :: batchCount

                            m = controldegree+1 ! number of rows of A
                            n = targetdegree+1 ! number of columns of A
                            alpha = 1.0_c_prec
                            lda = m ! leading dimension of A
                            strideA = 0 ! stride for the batches of A (no stride)
                            incx = targetdegree+1 !
                            stridex = (controldegree+1)*(targetdegree+1)
                            incy = targetdegree+1
                            stridey = (targetdegree+1)*(targetdegree+1)
                            batchCount = nvars*nelems
                            do i = 0,targetdegree
                              call hipblasCheck(hipblasgemvStridedBatched(handle, &
                                                                          HIPBLAS_OP_T, &
                                                                          m,n,alpha, &
                                                                          A,lda,strideA, &
                                                                          c_loc(f(1+i,1,1,1)),incx,stridex,beta, &
                                                                          c_loc(Af(1+i,1,1,1)),incy,stridey,batchCount))
                            enddo

                          endsubroutine self_blas_matrixop_dim2_2d

                          subroutine self_blas_matrixop_dim1_3d(A,f,Af,controldegree,targetdegree,nvars,nelems,handle)
                            type(c_ptr),intent(in) :: A
                            real(prec),pointer,intent(in) :: f(:,:,:,:,:)
                            real(prec),pointer,intent(inout) :: Af(:,:,:,:,:)
                            integer,intent(in) :: controldegree,targetdegree,nvars,nelems
                            type(c_ptr),intent(inout) :: handle
                            ! Local
                            integer(c_int) :: m
                            integer(c_int) :: n
                            integer(c_int) :: k
                            real(c_prec) :: alpha
                            integer(c_int) :: lda
                            integer(c_int) :: ldb
                            integer(c_int) :: ldc
                            real(c_prec) :: beta

                            m = targetdegree+1 ! number of rows of A^T
                            n = nvars*nelems*(controldegree+1)*(controldegree+1) ! number of columns of B
                            k = controldegree+1 ! number of columns of A^T
                            alpha = 1.0_c_prec
                            lda = k ! leading dimension of A (interoplation matrix)
                            ldb = k ! leading dimension of B (f)
                            ldc = m ! leading dimension of C (fTarget)
                            beta = 0.0_c_prec

                            ! First pass interpolates in the first quadrature dimension
                            call hipblasCheck(hipblasgemm(handle, &
                                                          HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                                          m,n,k,alpha, &
                                                          A,lda, &
                                                          c_loc(f),ldb,beta, &
                                                          c_loc(Af),ldc))

                          endsubroutine self_blas_matrixop_dim1_3d

                          subroutine self_blas_matrixop_dim2_3d(A,f,Af,beta,controldegree,targetdegree,nvars,nelems,handle)
                            type(c_ptr),intent(in) :: A
                            real(prec),pointer,intent(in) :: f(:,:,:,:,:)
                            real(prec),pointer,intent(inout) :: Af(:,:,:,:,:)
                            real(c_prec),intent(in) :: beta
                            integer,intent(in) :: controldegree,targetdegree,nvars,nelems
                            type(c_ptr),intent(inout) :: handle
                            ! Local
                            integer(c_int) :: m
                            integer(c_int) :: n
                            real(c_prec) :: alpha
                            integer(c_int) :: lda

                            integer :: i
                            integer(c_int64_t) :: strideA
                            integer(c_int) :: incx
                            integer(c_int64_t) :: stridex
                            integer(c_int) :: incy
                            integer(c_int64_t) :: stridey
                            integer(c_int) :: batchCount

                            m = controldegree+1 ! number of rows of A
                            n = targetdegree+1 ! number of columns of A
                            alpha = 1.0_c_prec
                            lda = m ! leading dimension of A
                            strideA = 0 ! stride for the batches of A (no stride)
                            incx = targetdegree+1 !
                            stridex = (controldegree+1)*(targetdegree+1)
                            !beta = 0.0_c_prec
                            incy = targetdegree+1
                            stridey = (targetdegree+1)*(targetdegree+1)
                            batchCount = (controldegree+1)*nvars*nelems
                            do i = 0,targetdegree
                              call hipblasCheck(hipblasgemvStridedBatched(handle, &
                                                                          HIPBLAS_OP_T, &
                                                                          m,n,alpha, &
                                                                          A,lda,strideA, &
                                                                          c_loc(f(1+i,1,1,1,1)),incx,stridex,beta, &
                                                                          c_loc(Af(1+i,1,1,1,1)),incy,stridey,batchCount))
                            enddo

                          endsubroutine self_blas_matrixop_dim2_3d

                          subroutine self_blas_matrixop_dim3_3d(A,f,Af,beta,controldegree,targetdegree,nvars,nelems,handle)
                            type(c_ptr),intent(in) :: A
                            real(prec),pointer,intent(in) :: f(:,:,:,:,:)
                            real(prec),pointer,intent(inout) :: Af(:,:,:,:,:)
                            real(c_prec),intent(in) :: beta
                            integer,intent(in) :: controldegree,targetdegree,nvars,nelems
                            type(c_ptr),intent(inout) :: handle
                            ! Local
                            integer(c_int) :: m
                            integer(c_int) :: n
                            real(c_prec) :: alpha
                            integer(c_int) :: lda

                            integer :: i,j
                            integer(c_int64_t) :: strideA
                            integer(c_int) :: incx
                            integer(c_int64_t) :: stridex
                            integer(c_int) :: incy
                            integer(c_int64_t) :: stridey
                            integer(c_int) :: batchCount

                            m = controldegree+1 ! number of rows of A
                            n = targetdegree+1 ! number of columns of A
                            alpha = 1.0_c_prec
                            lda = m ! leading dimension of A
                            strideA = 0 ! stride for the batches of A (no stride)
                            incx = (targetdegree+1)*(targetdegree+1) !
                            stridex = (controldegree+1)*(targetdegree+1)*(targetdegree+1)
                            incy = (targetdegree+1)*(targetdegree+1)
                            stridey = (targetdegree+1)*(targetdegree+1)*(targetdegree+1)
                            batchCount = nvars*nelems
                            do j = 0,targetdegree
                              do i = 0,targetdegree
                                call hipblasCheck(hipblasgemvStridedBatched(handle, &
                                                                            HIPBLAS_OP_T, &
                                                                            m,n,alpha, &
                                                                            A,lda,strideA, &
                                                                            c_loc(f(1+i,1+j,1,1,1)),incx,stridex,beta, &
                                                                            c_loc(Af(1+i,1+j,1,1,1)),incy,stridey,batchCount))
                              enddo
                            enddo

                          endsubroutine self_blas_matrixop_dim3_3d
                          endmodule SELF_GPUBLAS
