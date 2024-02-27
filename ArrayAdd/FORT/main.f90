program ArrayAdd
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
    integer,parameter :: DP = selected_real_kind(14,200)
    integer :: arraylen = 1000000
    integer :: nloop = 100000
    complex(DP),allocatable :: A(:), B(:), C(:)
    integer :: i, iloop
    real(DP) :: x, y
    real(DP) :: t_start, t_end, t
    !
    CHARACTER(len=12) :: N_read_c
    INTEGER :: num_args
    !   real(DP), external :: f_wall
    !
    num_args = command_argument_count()
    if(num_args==1)then
      CALL get_command_argument(1, N_read_c)
      READ(N_read_c,"(BN,I12)") arraylen
    endif

    ALLOCATE( A(arraylen), &
              B(arraylen), &
              C(arraylen)  )

    DO iloop = 1, nloop
        DO i = 1, arraylen
            call random_number( x )
            call random_number( y )
            B(i) = cmplx(x,y,kind=DP)
            call random_number( x )
            call random_number( y )
            C(i) = cmplx(x,y,kind=DP)
        End Do
    End Do

    t_start = f_wall()
    do iloop = 1, nloop
        do i = 1, arraylen
            A(i)=B(i)+C(i)
        end do
    end do
    t_end = f_wall()
    t = t_end - t_start

    ! write(*,'(" arrayadd ", I10,  " arrays", I12, " times in ", ES15.5, "s. Average ", ES20.10, "s for each array.")') &
    !     arraylen, nloop, t, t/nloop
    write(*,'(" arrayadd ", I10,  " arrays", I12, " times in ", ES15.5, "s. Average ", ES20.10, "s for each sum.")') &
        arraylen, nloop, t, t/nloop/arraylen



end program