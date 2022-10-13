program hello
    implicit NONE
    REAL(8)     :: U(3,3)
    REAL(8)     :: r0(3,12)
    REAL(8)     :: t1(3,12)
    REAL(8)     :: t2(3,12)
    REAl(8)     :: a
    REAL(8), PARAMETER :: RT2 = 1.41421356237310d0    ! SQRT(2)
    REAL(8), PARAMETER :: RT3 = 1.73205080756888d0    ! SQRT(3)
    REAL(8), PARAMETER :: RT6 = 2.44948974278318d0    ! SQRT(6)

    a=4.174999999
    U = RESHAPE ( (/ -1d0/RT2,0d0,1d0/RT2, 1d0/RT6,-2d0/RT6,1d0/RT6, -1d0/RT3,-1d0/RT3,-1d0/RT3 /), (/3,3/) )
    t1=RESHAPE ( (/0,1,1,0,1,-1,1,0,1,-1,0,1,1,1,0,1,-1,0,0, &
    -1,-1,0,-1,1,-1,0,-1,1,0,-1,-1,-1,0,-1,1,0/),(/3,12/) )
    t2=t1/2.0d0*a
    r0 = MATMUL(TRANSPOSE(U),t2)

    print *, t1(1,:)
    print *, ''
    print *, t2(1,:)
    print *, ''
    print *, r0(1,:)
end program hello