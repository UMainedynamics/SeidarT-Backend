module plane_wave_source

use iso_fortran_env, only: real64

implicit none 
private 

public :: select_injection_faces 
public :: rotmat 
public :: direction_polarization_vector 
public :: boundary_em_field 

real(real64), parameter :: PI = 3.1415926535897932384626433832795_real64

contains 

pure function cross_product(a, b) result(c)
    real(real64), intent(in) :: a(3), b(3)
    real(real64)             :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
end function cross_product

pure subroutine select_injection_faces(p, active)
    implicit none
    real(real64), intent(in)  :: p(3)
    real(real64), parameter :: eps = 1.0_real64 * 1.0e-6_real64
    logical,  intent(out) :: active(6)

    ! IDs:
    ! 1=XMIN, 2=XMAX, 3=YMIN, 4=YMAX, 5=ZMIN, 6=ZMAX

    active = .false.

    if (p(1) >  eps) active(1) = .true.
    if (p(1) < -eps) active(2) = .true.

    if (p(2) >  eps) active(3) = .true.
    if (p(2) < -eps) active(4) = .true.
  
    if (p(3) >  eps) active(5) = .true.
    if (p(3) < -eps) active(6) = .true.
end subroutine select_injection_faces

pure function rotmat(angle_xz, angle_xy, angle_yz) result(R)
    real(real64), intent(in) :: angle_xz, angle_xy, angle_yz
    real(real64) :: Rx(3,3), Ry(3,3), Rz(3,3)
    real(real64) :: cx, cy, cz, sx, sy, sz
    real(real64) :: R(3,3)
    real(real64), parameter :: deg2rad = pi/180

    cx = cos(deg2rad * angle_yz)
    cy = cos(deg2rad * angle_xz)
    cz = cos(deg2rad * angle_xy)
    sx = sin(deg2rad * angle_yz)
    sy = sin(deg2rad * angle_xz)
    sz = sin(deg2rad * angle_xy)

    Rx = reshape([ & 
        1.0_real64, 0.0_real64,  0.0_real64, &
        0.0_real64, cx,     -sx, &
        0.0_real64, sx,      cx  &
    ], [3,3])
    Ry = reshape([ & 
        cy,     0.0_real64, sy, &
        0.0_real64, 1.0_real64, 0.0_real64, &
        -sy,    0.0_real64, cy &
    ], [3,3])
    Rz = reshape([ & 
        cz,     -sz,    0.0_real64, &
        sz,     cz,     0.0_real64, &
        0.0_real64, 0.0_real64, 1.0_real64 &
    ], [3,3])
    R = matmul(Ry, matmul(Rz, Rx) )

end function rotmat

pure subroutine direction_polarization_vector(angle_xz, angle_xy, angle_yz, p, ehat)
    real(real64), intent(in) :: angle_xz, angle_xy, angle_yz
    real(real64), intent(out) :: p(3), ehat(3) 
    real(real64) :: nrm
    real(real64) :: R(3,3) 
    real(real64) :: k0(3), e0(3), tmp(3) 

    ! Reference propagation and polarization
    k0 = [0.0_real64, 0.0_real64, 1.0_real64]
    e0 = [1.0_real64, 0.0_real64, 0.0_real64]

    R = rotmat(angle_xz, angle_xy, angle_yz) 
    p = matmul(R, k0)
    nrm = norm2(p) 
    if (nrm > 0.0_real64) p = p / nrm
    
    tmp = matmul(R, e0) 

    ! Remove any component parallel to p so E is transverse
    tmp = tmp - dot_product(tmp, p) * p
    nrm = norm2(tmp)
    if (nrm > 0.0_real64) then 
        ehat = tmp / nrm 
    else
        if (abs(p(1)) < 0.9_real64) then
            tmp = [1.0_real64, 0.0_real64, 0.0_real64]
        else
            tmp = [0.0_real64, 1.0_real64, 0.0_real64]
        end if

      tmp = tmp - dot_product(tmp, p) * p
      nrm = norm2(tmp) 
      if (nrm > 0.0_real64) then 
        ehat = tmp / nrm 
      else
        ehat = [0.0_real64, 0.0_real64, 0.0_real64]
      end if 
    end if

end subroutine direction_polarization_vector 

pure real(real64) function wavelet(t, f0, type) result(w) 
    real(real64), intent(in) :: t, f0 
    character(len=*), intent(in) :: type 
    real(real64) :: t0, a, g

    t0 = 1.5_real64 / f0 
    a = PI * f0 * (t - t0)
    g = exp(-a * a) 

    select case (type) 
        case ('gaus0')
            w = g 
        case ('gaus1')
            w = -2.0_real64 * a * g 
        case ('gaus2')
            w = (1.0_real64 - 2.0_real64 * a * a) * g 
    end select
    
end function wavelet

subroutine boundary_em_field(r, r0, t, f0, &
                                    angle_xz, angle_xy, angle_yz, &
                                    v, eta, amp, type, E, H)
    real(real64), intent(in) :: t 
    real(real64), intent(in) :: r(3), r0(3)
    real(real64), intent(in) :: angle_xz, angle_xy, angle_yz 
    real(real64), intent(in) :: v ! background wave speed 
    real(real64), intent(in) :: amp, f0, eta ! source amplitude, center frequency, impedence (ohms)
    character(len=*), intent(in) :: type ! Wavelet type
    real(real64), intent(out) :: E(3), H(3)
    

    real(real64) :: p(3), ehat(3), hhat(3)  
    real(real64) :: tau, s, hnrm

    call direction_polarization_vector(angle_xz, angle_xy, angle_yz, p, ehat)
    
    ! tau = (sin(angle_xz) / v ) * ( (r(1) - r0(1) ) * cos(angle_xy) * (r(2) - r0(2) ) ) + cos(angle_xz) * (r(3) - r0(3) ) / v
    tau = t - dot_product(p, r-r0) / v
    ! print*, t - tau
    s = amp * wavelet(tau, f0, type) 
    
    E = s*ehat 
    hhat = cross_product(p, ehat)
    hnrm = norm2(hhat) 
    if (hnrm > 0.0_real64) then 
        hhat = hhat / hnrm 
    end if 

    H = (s / eta) * hhat

end subroutine boundary_em_field

end module plane_wave_source