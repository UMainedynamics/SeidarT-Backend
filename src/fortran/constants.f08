module constants

    use iso_fortran_env, only: real64
    
    implicit none
    
    real(real64), parameter :: pi = 3.141592653589793238462643d0
    real(real64), parameter :: clight = 2.9979458d+8
    real(real64), parameter :: mu0 = 4.0d0 * pi * 1.0d-7
    real(real64), parameter :: eps0 = 8.85418782d-12
    real(real64), parameter :: mu = 1.d0
    real(real64), parameter :: STABILITY_THRESHOLD = 1.d+25

end module constants
