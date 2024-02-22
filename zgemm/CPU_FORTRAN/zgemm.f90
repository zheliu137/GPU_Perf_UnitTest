MODULE mytime
  !----------------------------------------------------------------------------
  USE ISO_C_BINDING
  !
  IMPLICIT NONE
  !
  SAVE
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
END MODULE mytime
  
!--------------------------------------------------------------------------
program main
!--------------------------------------------------------------------------
!
  ! USE ISO_C_BINDING
  USE mytime
  !
  implicit none
  !
  ! INTERFACE
  !    FUNCTION f_wall ( ) BIND(C,name="cclock") RESULT(t)
  !      USE ISO_C_BINDING
  !      REAL(kind=c_double) :: t
  !    END FUNCTION f_wall
  !    FUNCTION f_tcpu ( ) BIND(C,name="scnds") RESULT(t)
  !      USE ISO_C_BINDING
  !      REAL(kind=c_double) :: t
  !    END FUNCTION f_tcpu
  ! END INTERFACE
  !
  INTEGER, PARAMETER :: DP=16
  !
  ! INTEGER, PARAMETER :: N=16
  INTEGER :: N=16
  !
  INTEGER, PARAMETER :: Nmat=1000000

  COMPLEX(DP) :: cone=1.0, czero=0.0
  !
  !! number of bands (possibly of the optimal subspace)
  ! 
  INTEGER :: nbnd
  !
  ! work variables 
  !
  INTEGER :: ibnd, jbnd, iloop
  COMPLEX(kind=DP), allocatable :: chf(:,:), chf2(:,:), chf3(:,:)
  real(DP) :: x,y
  !
  real(DP) :: t_start, t_end

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

  allocate(chf(nbnd,nbnd),  &
           chf2(nbnd,nbnd), &
           chf3(nbnd,nbnd))

  call random_seed() 
  DO jbnd = 1, nbnd
    DO ibnd = 1, nbnd
        call random_number( x )
        call random_number( y )
        chf(ibnd,jbnd)=cmplx(x,y,kind=DP)
    End Do
  End Do

  DO jbnd = 1, nbnd
    DO ibnd = 1, nbnd
        call random_number( x )
        call random_number( y )
        chf2(ibnd,jbnd)=cmplx(x,y,kind=DP)
    End Do
  End Do
  !
  t_start = f_wall()
!   t_start = cclock()
  !
  DO iloop = 1, nmat

    CALL zgemm('N', 'N',  N, N, N, cone, chf, N, chf2, N, czero, chf3, N );

    ! CALL zhpevx ('V', 'A', 'U', nbnd, champ , 0.0, 0.0, &
    !            0, 0,-1.0, neig, w, cz, nbnd, cwork, &
    !            rwork, iwork, ifail, info)
  ENDDO
  !
  ! rotation matrix and Ham eigenvalues 
  ! [in Ry, mind when comparing with wannier code]
  ! 
  !
  t_end = f_wall()
!   t_end = cclock()
  !
  ! write(*,'(I4, " by ", I4, "   matrices multiplication ", I12, " times finished in ", ES20.10, " seconds.")')  &
  !       N, N, nmat, t_end - t_start
  ! write(*,'(" average ", ES20.10, " seconds for each." )')(t_end - t_start)/nmat
  write(*,'(I4, " by ", I4, "   matmul ", I10, " times in ", ES12.4, " seconds.", " average ", ES20.10, "s for each.")')  &
        N, N, nmat, t_end - t_start, (t_end - t_start)/nmat
  !
end program

