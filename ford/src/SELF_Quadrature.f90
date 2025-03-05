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
!
! Contains routines from D.A. Kopriva, 2009, "Implementing Spectral Methods for Partial
! Differential Equations: Algorithms for Scientists and Engineers", Springer.
!
! Routines are defined for computing Legendre and Chebyshev Gauss and Gauss-Lobatto
! quadrature nodes and weights.

module SELF_Quadrature

  use iso_fortran_env
  use SELF_Constants

  implicit none

  public  :: ChebyshevQuadrature,LegendreQuadrature
  private :: ChebyshevGauss,ChebyshevGaussLobatto, &
             LegendreGauss,LegendreGaussLobatto, &
             LegendreQandL

contains

! =============================================================================================== !
! LegendreQuadrature
!   Returns the specified Legendre quadrature nodes and integration weights.
!
!   Given a polynomial degree, and quadrature type (Gauss or Gauss Lobatto), this subroutine manages
!   the calls to underlying private routines to generate the desired Legendre quadrature.
!
!   Usage :
!
!     INTEGER    :: N, quadType
!     REAL(real64) :: nodes(0:N), weights(0:N)
!
!       CALL LegendreQuadrature( N, quadType, nodes, weights )
!
!   Parameters :
!
!     N (in)
!       Degree of the quadrature
!
!     quadType (in)
!       Flag specifying the quadrature type. Can be set to GAUSS or GAUSS_LOBATTO
!
!     nodes(0:N) (out)
!       Array of quadrature nodes
!
!     weights(0:N) (out)
!       Array of quadrature weights
!
! =============================================================================================== !

  subroutine LegendreQuadrature(N,nodes,weights,QuadType)
    implicit none
    integer,intent(in)     :: N
    real(prec),intent(out) :: nodes(0:N)
    real(prec),intent(out) :: weights(0:N)
    integer,intent(in)     :: QuadType
    ! Local
    real(real64) :: nodesLocal(0:N)
    real(real64) :: weightsLocal(0:N)
    integer :: i

    if(QuadType == GAUSS_LOBATTO) then

      call LegendreGaussLobatto(N,nodesLocal,weightsLocal)

    elseif(QuadType == GAUSS) then

      call LegendreGauss(N,nodesLocal,weightsLocal)

    endif

    do i = 0,N
      nodes(i) = real(nodesLocal(i),prec)
      weights(i) = real(weightsLocal(i),prec)
    enddo

  endsubroutine LegendreQuadrature

! =============================================================================================== !
! ChebyshevQuadrature
!
!   Returns the specified Chebyshev quadrature nodes and integration weights.
!
!   Given a polynomial degree, and quadrature type (Gauss or Gauss Lobatto), this subroutine manages
!   the calls to underlying private routines to generate the desired Chebyshev quadrature.
!
!   Usage :
!
!     INTEGER    :: N, quadType
!     REAL(real64) :: nodes(0:N), weights(0:N)
!
!       CALL ChebyshevQuadrature( N, quadType, nodes, weights )
!
!   Input/Output :
!
!     N (in)
!       Degree of the quadrature
!
!     quadType (in)
!       Flag specifying the quadrature type. Can be set to GAUSS or GAUSS_LOBATTO
!
!     nodes(0:N) (out)
!       Array of quadrature nodes
!
!     weights(0:N) (out)
!       Array of quadrature weights
!
! ================================================================================================ !

  subroutine ChebyshevQuadrature(N,nodes,weights,quadType)
    implicit none
    integer,intent(in)     :: N
    real(prec),intent(out) :: nodes(0:N)
    real(prec),intent(out) :: weights(0:N)
    integer,intent(in)     :: QuadType
    ! Local
    real(real64) :: nodesLocal(0:N)
    real(real64) :: weightsLocal(0:N)
    integer :: i

    if(QuadType == CHEBYSHEV_GAUSS_LOBATTO) then

      call ChebyshevGaussLobatto(N,nodesLocal,weightsLocal)

    elseif(QuadType == CHEBYSHEV_GAUSS) then

      call ChebyshevGauss(N,nodesLocal,weightsLocal)

    endif

    do i = 0,N
      nodes(i) = real(nodesLocal(i),prec)
      weights(i) = real(weightsLocal(i),prec)
    enddo

  endsubroutine ChebyshevQuadrature

! =============================================================================================== !
! S/R ChebyshevGauss
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 67
!   Algorithm 26
! =============================================================================================== !

  subroutine ChebyshevGauss(N,nodes,weights)
    implicit none
    integer       :: N
    real(real64)    :: nodes(0:N)
    real(real64)    :: weights(0:N)
    ! Local
    integer    :: j

    do j = 0,N

      weights(j) = pi/(real(N,real64)+1.0_real64)
      nodes(j) = -cos(pi*(2.0_real64*real(j,real64)+1.0_real64)/(2.0_real64*real(N,real64)+2.0_real64))

    enddo

  endsubroutine ChebyshevGauss

! =============================================================================================== !
! S/R ChebyshevGaussLobatto
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 68
!   Algorithm 27
! =============================================================================================== !

  subroutine ChebyshevGaussLobatto(N,nodes,weights)
    implicit none
    integer       :: N
    real(real64)    :: nodes(0:N)
    real(real64)    :: weights(0:N)
    ! LOCAL
    integer    :: j

    do j = 0,N

      weights(j) = pi/real(N,real64)
      nodes(j) = -cos(pi*real(j,real64)/real(N,real64))

    enddo

    weights(0) = weights(0)*0.5_real64
    weights(N) = weights(N)*0.5_real64

  endsubroutine ChebyshevGaussLobatto

! =============================================================================================== !
! S/R LegendreGauss
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 64
!   Algorithm 23
! =============================================================================================== !

  subroutine LegendreGauss(N,nodes,weights)
    implicit none
    integer    :: N
    real(real64) :: nodes(0:N)
    real(real64) :: weights(0:N)
    ! Local
    real(real64) :: nodes_local(0:N)
    real(real64) :: weights_local(0:N)
    real(real64) :: lN1,dlN1
    real(real64) :: delta
    integer  :: j,kIt

    if(N == 0) then

      nodes_local(0) = 0.0_real64
      weights_local(0) = 2.0_real64

    elseif(N == 1) then

      nodes_local(0) = -sqrt(1.0_real64/3.0_real64)
      weights_local(0) = 1.0_real64
      nodes_local(1) = -nodes(0)
      weights_local(1) = weights(0)

    else

      do j = 0,((N+1)/2)

        nodes_local(j) = -cos((2.0_real64*real(j,real64)+1.0_real64)*pi/(2.0_real64*real(N,real64)+1.0_real64))

        do kIt = 1,newtonMax

          call LegendrePolynomial(N+1,nodes_local(j),lN1,dlN1)
          delta = -lN1/dlN1
          nodes_local(j) = nodes_local(j)+delta
          if(abs(delta) <= TOL*nodes_local(j)) exit

        enddo

        call LegendrePolynomial(N+1,nodes_local(j),lN1,dlN1)
        weights_local(j) = 2.0_real64/((1.0_real64-nodes_local(j)*nodes_local(j))*dlN1*dlN1)
        weights_local(N-j) = weights_local(j)
        nodes_local(N-j) = -nodes_local(j)

      enddo

    endif

    if(mod(real(N,real64),2.0_real64) == 0.0_real64) then

      call LegendrePolynomial(N+1,0.0_real64,lN1,dlN1)
      nodes_local(N/2) = 0.0_real64
      weights_local(N/2) = 2.0/(dlN1*dlN1)

    endif

    do j = 0,N
      nodes(j) = real(nodes_local(j),real64)
      weights(j) = real(weights_local(j),real64)
    enddo

  endsubroutine LegendreGauss

  ! =============================================================================================== !
  ! S/R LegendreGaussLobatto
  !   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 66
  !   Algorithm 25
  ! =============================================================================================== !

  subroutine LegendreGaussLobatto(N,nodes,weights)
    implicit none
    integer    :: N
    real(real64) :: nodes(0:N)
    real(real64) :: weights(0:N)
    ! Local
    real(real64) :: nodes_local(0:N)
    real(real64) :: weights_local(0:N)
    real(real64) :: delta,q,qprime,lN
    integer  :: j,kIt

    if(N == 1) then

      nodes_local(0) = -1.0_real64
      weights_local(0) = 1.0_real64
      nodes_local(1) = 1.0_real64
      weights_local(1) = 1.0_real64

    else

      nodes_local(0) = -1.0_real64
      weights_local(0) = 2.0_real64/(real(N,real64)*(real(N,real64)+1.0_real64))
      nodes_local(N) = 1.0_real64
      weights_local(N) = weights_local(0)

      do j = 1,((N+1)/2-1)

        nodes_local(j) = -cos((real(j,real64)+0.25_real64)*pi/real(N,real64)- &
                              3.0_real64/(8.0_real64*real(N,real64)*pi*(real(j,real64)+0.25_real64)))

        do kIt = 1,newtonMax

          call LegendreQandL(N,nodes_local(j),q,qprime,lN)

          delta = -q/qprime
          nodes_local(j) = nodes_local(j)+delta
          if(abs(delta) <= TOL*nodes_local(j)) exit

        enddo

        call LegendreQandL(N,nodes_local(j),q,qprime,lN)

        weights_local(j) = 2.0_real64/(real(N,real64)*(real(N,real64)+1.0_real64)*lN*lN)
        weights_local(N-j) = weights_local(j)
        nodes_local(N-j) = -nodes_local(j)

      enddo

    endif

    if(mod(real(N,real64),2.0_real64) == 0.0_real64) then

      call LegendreQandL(N,0.0_real64,q,qprime,lN)

      nodes_local(N/2) = 0.0_real64
      weights_local(N/2) = 2.0_real64/(real(N,real64)*(real(N,real64)+1.0_real64)*lN*lN)

    endif

    do j = 0,N
      nodes(j) = real(nodes_local(j),real64)
      weights(j) = real(weights_local(j),real64)
    enddo

  endsubroutine LegendreGaussLobatto

! =============================================================================================== !
! S/R LegendrePolynomial
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 63
!   Algorithm 22
! =============================================================================================== !

  subroutine LegendrePolynomial(N,x,lAtX,dLdxAtX)
    implicit none
    integer     :: N
    real(real64)    :: x
    real(real64)    :: lAtX,dLdxAtX
    ! Local
    real(real64) :: lNm1,lNm2,dlNm1,dlNm2
    integer  :: i

    if(N == 0) then

      lAtX = 1.0_real64
      dLdxAtX = 0.0_real64

    elseif(N == 1) then

      lAtX = x
      dLdxAtX = 1.0_real64

    else

      lnM2 = 1.0_real64
      lnM1 = x
      dlnM2 = 0.0_real64
      dlnM1 = 1.0_real64

      do i = 2,N

        lAtX = ((2.0_real64*real(i,real64)-1.0_real64)*x*lnM1- &
                (real(i,real64)-1.0_real64)*lnM2)/(real(i,real64))

        dldxAtX = dlnM2+(2.0_real64*real(i,real64)-1.0_real64)*lnM1
        lnM2 = lnM1
        lnM1 = lAtX
        dlnM2 = dlnM1
        dlnM1 = dldxAtX

      enddo

    endif

  endsubroutine LegendrePolynomial

! =============================================================================================== !
! S/R LegendreQandL
!   From Kopriva (2009) "Implementing Spectral Methods for Partial Differential Equations", p. 65
!   Algorithm 24
! =============================================================================================== !

  subroutine LegendreQandL(N,x,q,qprime,lN)
    implicit none
    integer    :: N
    real(real64) :: x
    real(real64) :: lN,q,qprime
    ! Local
    real(real64) :: lNm1,lNm2,dlNm1,dlNm2,dlN,lN1,dlN1
    integer    :: i

    lNm2 = 1.0_real64
    lNm1 = x
    dlNm2 = 0.0_real64
    dlNm1 = 1.0_real64

    do i = 2,N

      lN = (2.0_real64*i-1.0_real64)/(real(i,real64))*x*lNm1-(real(i,real64)-1.0_real64)/(real(i,real64))*lNm2
      dlN = dlNm2+(2.0_real64*real(i,real64)-1.0_real64)*lNm1
      lNm2 = lNm1
      lNm1 = lN
      dlNm2 = dlNm1
      dlNm1 = dlN

    enddo

    i = N+1
    lN1 = (2.0_real64*i-1.0_real64)/(real(i,real64))*x*lN-(real(i,real64)-1.0_real64)/(real(i,real64))*lNm2
    dlN1 = dlNm2+(2.0_real64*real(i,real64)-1.0_real64)*lNm1
    q = lN1-lNm2
    qprime = dlN1-dlNm2

  endsubroutine LegendreQandL

endmodule SELF_Quadrature
