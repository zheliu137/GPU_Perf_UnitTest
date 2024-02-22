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
  ! INTEGER, PARAMETER :: N=32
  INTEGER :: N=32
  !
  INTEGER, PARAMETER :: Nmat=100000

  ! COMPLEX(DP) :: cone=1.0, czero=0.0
  !
  ! variables for lapack ZHPEVX 
  !
  INTEGER :: neig, info
  INTEGER, ALLOCATABLE :: ifail( : ), iwork( : )
  REAL(kind=DP), ALLOCATABLE :: w( : )
  REAL(kind=DP), ALLOCATABLE :: rwork( : )
  COMPLEX(kind=DP), ALLOCATABLE :: champ( : )
  COMPLEX(kind=DP), ALLOCATABLE :: cwork( : )
  COMPLEX(kind=DP), ALLOCATABLE :: cz( :, :)
  COMPLEX(kind=DP), ALLOCATABLE :: chf(:, :) 
  !
  ! work variables 
  !
  INTEGER :: ibnd, jbnd, iloop
  INTEGER :: nbnd
  real(DP) :: x, y
  !
  real(DP) :: t_start, t_end, t
  !   real(DP), external :: f_wall
  !
  CHARACTER(len=12) :: N_read_c
  INTEGER :: num_args
  !   real(DP), external :: f_wall
  !

  num_args = command_argument_count()
  if(num_args==1)then
    CALL get_command_argument(1, N_read_c)
    READ(N_read_c,"(BN,I12)") N
  endif

  nbnd = N

  ALLOCATE(ifail( N ), iwork( 5*N ), &
           w( N ), rwork( 7*N ), &
           champ( N*(N+1)/2 ), cwork( 2*N ), &
           cz( N, N), chf(N, N))

  call random_seed() 
  DO jbnd = 1, nbnd
    DO ibnd = 1, nbnd
        call random_number( x )
        call random_number( y )
        chf(ibnd,jbnd)=cmplx(x, y, kind=DP)
    End Do
  End Do

  DO jbnd = 1, nbnd
   DO ibnd = 1, jbnd
      champ (ibnd + (jbnd - 1) * jbnd/2 ) = &
      ( chf ( ibnd, jbnd) + conjg ( chf ( jbnd, ibnd) ) ) * 0.5d0
   ENDDO
  ENDDO
  !
  !   t_start = cclock()
  !
  t=0.0D0
  DO iloop = 1, nmat
    DO jbnd = 1, nbnd
      DO ibnd = 1, nbnd
          call random_number( x )
          call random_number( y )
          chf(ibnd,jbnd)=cmplx(x, y, kind=DP)
      End Do
    End Do
    DO jbnd = 1, nbnd
      DO ibnd = 1, jbnd
        champ (ibnd + (jbnd - 1) * jbnd/2 ) = &
        ( chf ( ibnd, jbnd) + conjg ( chf ( jbnd, ibnd) ) ) * 0.5d0
      ENDDO
    ENDDO
    t_start = f_wall()
    CALL zhpevx ('V', 'A', 'U', nbnd, champ , 0.0, 0.0, &
               0, 0,-1.0, neig, w, cz, nbnd, cwork, &
               rwork, iwork, ifail, info)
    t_end = f_wall()
    t = t + t_end - t_start
  ENDDO
  !
  ! rotation matrix and Ham eigenvalues 
  ! [in Ry, mind when comparing with wannier code]
  ! 
  !
!   t_end = cclock()
  !
  write(*,'(" diag ", I10, I4, " by ", I4, "    matrices in ", ES12.4, "s. Average ", ES20.10, "s for each.")') &
               nmat, N, N, t, t/nmat
  !
end program

