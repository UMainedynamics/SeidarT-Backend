module density_averaging
  implicit none
contains

  pure function face_density(a, b, code) result(r)
    !-- pick one of ’harmonic’, ’geometric’, ’arithmetic’, or ’none’
    real(kind=8), intent(in) :: a, b
    integer, intent(in) :: code
    real(kind=8) :: r

    select case (code)
    case (1) !'harmonic'
      r = 2d0 * a*b / (a + b)
    case (2) !'geometric'
      r = sqrt(a*b)
    case (3) !'arithmetic'
      r = 0.5d0*(a + b)
    case (4) !'none'
      r = a    ! sentinel—means “use deltarho fallback”
    end select
  end function face_density


end module density_averaging