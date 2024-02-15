!   MODULE mytime
!   !----------------------------------------------------------------------------
!   USE ISO_C_BINDING
!   !
!   IMPLICIT NONE
!   !
!   SAVE
!   !
!   INTERFACE
!      FUNCTION f_wall ( ) BIND(C,name="cclock") RESULT(t)
!        USE ISO_C_BINDING
!        REAL(kind=c_double) :: t
!      END FUNCTION f_wall
!      FUNCTION f_tcpu ( ) BIND(C,name="scnds") RESULT(t)
!        USE ISO_C_BINDING
!        REAL(kind=c_double) :: t
!      END FUNCTION f_tcpu
!   END INTERFACE
!   !
! END MODULE mytime
  
  !--------------------------------------------------------------------------
  program main
  !--------------------------------------------------------------------------
  !
  USE ISO_C_BINDING
!   USE mytime
  !
  implicit none
  !
  INTERFACE
     FUNCTION f_wall ( ) BIND(C,name="cclock") RESULT(t)
       USE ISO_C_BINDING
       REAL(kind=c_double) :: t
     END FUNCTION f_wall
     FUNCTION f_tcpu ( ) BIND(C,name="scnds") RESULT(t)
       USE ISO_C_BINDING
       REAL(kind=c_double) :: t
     END FUNCTION f_tcpu
  END INTERFACE
  !
  INTEGER, PARAMETER :: DP=16
  !
  INTEGER, PARAMETER :: N=16
  !
  INTEGER, PARAMETER :: Nmat=1000000

  COMPLEX(DP) :: cone=1.0, czero=0.0
  !
  !! number of bands (possibly of the optimal subspace)
  ! 
  INTEGER :: nbnd
  REAL(kind=DP):: eig (N)
  !! interpolated hamiltonian eigenvalues for this kpoint 
  ! 
  !! Hamiltonian in Bloch basis, fine mesh
  COMPLEX(kind=DP) :: cuf(N, N)
  !! Rotation matrix, fine mesh
  !
  ! variables for lapack ZHPEVX 
  !
  INTEGER :: neig, info, ifail( N ), iwork( 5*N )
  REAL(kind=DP) :: w( N )
  REAL(kind=DP) :: rwork( 7*N )
  COMPLEX(kind=DP) :: champ( N*(N+1)/2 )
  COMPLEX(kind=DP) :: cwork( 2*N )
  COMPLEX(kind=DP) :: cz( N, N)
  !
  ! work variables 
  !
  INTEGER :: ibnd, jbnd, iloop
  COMPLEX(kind=DP) :: chf(N, N) 
  real(DP) :: x,y
  !
  real(DP) :: t_start, t_end
!   real(DP), external :: f_wall
  !
  nbnd = N

  call random_seed() 
  DO jbnd = 1, nbnd
    DO ibnd = 1, nbnd
        call random_number( x )
        call random_number( y )
        chf(ibnd,jbnd)=cmplx(x,y)
    End Do
  End Do

  DO jbnd = 1, nbnd
   DO ibnd = 1, jbnd
      champ (ibnd + (jbnd - 1) * jbnd/2 ) = &
      ( chf ( ibnd, jbnd) + conjg ( chf ( jbnd, ibnd) ) ) * 0.5d0
   ENDDO
  ENDDO
  !
  t_start = f_wall()
!   t_start = cclock()
  !
  DO iloop = 1, nmat
    CALL zhpevx ('V', 'A', 'U', nbnd, champ , 0.0, 0.0, &
               0, 0,-1.0, neig, w, cz, nbnd, cwork, &
               rwork, iwork, ifail, info)
  ENDDO
  !
  ! rotation matrix and Ham eigenvalues 
  ! [in Ry, mind when comparing with wannier code]
  ! 
  !
  t_end = f_wall()
!   t_end = cclock()
  !
  write(*,'(" diag ", I12, "    matrices in ", ES20.10, " seconds")') nmat, t_end - t_start
  write(*,'(" average ", ES20.10, " seconds for one." )')(t_end - t_start)/nmat
  !
end program

