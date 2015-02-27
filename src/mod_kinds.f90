module kinds_module

  implicit none

! Define some useful constants / parameters
  integer,           parameter :: rkind = 8, &!kind(1.d0), &
                                  ikind = 4, &
                                  ckind = 4

  real (kind=rkind), parameter :: pi = 4d0 * atan(1.0d0)

contains

!-------------------------------------------------------------------

  subroutine test_double( x )
    real (kind=rkind), intent(in) :: x

!-- F2PY VARIABLE BINDINGS
    ! f2py intent(in) :: x

    print *, "x is", x

  end subroutine test_double

  subroutine test_integer( x )
    integer (kind=ikind), intent(in) :: x

!-- F2PY VARIABLE BINDINGS
    ! f2py intent(in) :: x

    print *, "x is", x

  end subroutine test_integer

  subroutine test_complex( x )
    complex (kind=ckind), intent(in) :: x

!-- F2PY VARIABLE BINDINGS
    ! f2py intent(in) :: x

    print *, "x is", x
    print *, "  real(x) =", real(x)
    print *, "  imag(x) =", aimag(x)

  end subroutine test_complex

end module kinds_module