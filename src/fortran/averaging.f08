module averaging
  implicit none
contains

  pure function scalar_mean(a, b, code) result(r)
    !-- pick one of ’harmonic’, ’geometric’, ’arithmetic’, or ’none’
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in)  :: a, b
    integer,      intent(in)  :: code
    real(real64)              :: r

    select case (code)
    case (1) !'harmonic'
      r = 2d0 * a*b / (a + b)
    case (2) !'geometric'
      r = sqrt(a*b)
    case (3) !'arithmetic'
      r = 0.5d0*(a + b)
    case (4) !'none'
      r = a    ! sentinel—means “use deltarho fallback”
    case default 
      error stop "face_density in density_averaging.f08: invalid code"
    end select
  end function scalar_mean

  subroutine array_averaging2d(array, code, direction)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none 
    real(real64), intent(inout)   :: array(:,:)
    integer, intent(in)           :: code 
    character(len=*), intent(in)  :: direction 
    integer                       :: ni, nj, i, j
    
    ni = size(array, 1) 
    nj = size(array, 2) 
    
    select case (trim(direction))
    case('+i')
      do j = 1,nj 
        do i = 2,ni
          array(i,j) = scalar_mean(array(i,j), array(i-1,j), code)
        end do
      end do
      
    case('-i')
      do j = 1, nj
        do i = 1, ni-1
          array(i,j) = scalar_mean(array(i,j), array(i+1,j), code)
        end do
      end do
    
    case('+j')
      do j = 2, nj
        do i = 1, ni
          array(i,j) = scalar_mean(array(i,j), array(i,j-1), code)
      end do
    end do

    case('-j')
      do j = 1, nj-1
        do i = 1, ni
          array(i,j) = scalar_mean(array(i,j), array(i,j+1), code)
        end do
      end do
    
    case('cc')   ! 4-cell average: (i,j), (i-1,j), (i,j-1), (i-1,j-1)
      do j = 2, nj
        do i = 2, ni
          array(i,j) = 0.25_real64 * ( &
                array(i,j)   + array(i-1,j) + &
                array(i,j-1) + array(i-1,j-1) )
        end do
      end do
      
    case default
      error stop "array_averaging2d: invalid direction"
    end select
    
  end subroutine array_averaging2d

  subroutine array_averaging3d(array, code, direction)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    real(real64), intent(inout)   :: array(:,:,:)
    integer,      intent(in)      :: code
    character(len=*), intent(in)  :: direction
    integer                       :: ni, nj, nk, i, j, k

    ni = size(array,1)
    nj = size(array,2)
    nk = size(array,3)

    select case (trim(direction))
    case('+i')
      do k = 1, nk
        do j = 1, nj
          do i = 2, ni
            array(i,j,k) = scalar_mean(array(i,j,k), array(i-1,j,k), code)
          end do
        end do
      end do

    case('-i')
      do k = 1, nk
        do j = 1, nj
          do i = 1, ni-1
            array(i,j,k) = scalar_mean(array(i,j,k), array(i+1,j,k), code)
          end do
        end do
      end do

    case('+j')
      do k = 1, nk
        do j = 2, nj
          do i = 1, ni
            array(i,j,k) = scalar_mean(array(i,j,k), array(i,j-1,k), code)
          end do
        end do
      end do

    case('-j')
      do k = 1, nk
        do j = 1, nj-1
          do i = 1, ni
            array(i,j,k) = scalar_mean(array(i,j,k), array(i,j+1,k), code)
          end do
        end do
      end do

    case('+k')
      do k = 2, nk
        do j = 1, nj
          do i = 1, ni
            array(i,j,k) = scalar_mean(array(i,j,k), array(i,j,k-1), code)
          end do
        end do 
      end do

    case('-k')
      do k = 1, nk-1 
        do j = 1, nj
          do i = 1, ni
            array(i,j,k) = scalar_mean(array(i,j,k), array(i,j,k+1), code)
          end do
        end do 
      end do
    
    case('cc')
      do k = 2, nk
        do j = 2, nj
          do i = 2, ni
            array(i,j,k) = 0.125_real64 * ( &
                array(i,  j,  k  ) + array(i-1,j,  k  ) + &
                array(i,  j-1,k  ) + array(i-1,j-1,k  ) + &
                array(i,  j,  k-1) + array(i-1,j,  k-1) + &
                array(i,  j-1,k-1) + array(i-1,j-1,k-1) )
          end do
        end do
      end do
      
    case default
      error stop "array_averaging3d: invalid direction"
    end select
  end subroutine array_averaging3d

end module averaging
