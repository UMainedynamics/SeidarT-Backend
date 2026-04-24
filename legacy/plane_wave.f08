program source

implicit none 
contains 

    ! ==========================================================================
    subroutine order_modes(V, e)
    !! Sorts the three wave modes so that:
    !!   V(3) = P  (largest velocity)
    !!   V(1),V(2) = S1,S2  (smaller, stable ordering)
    !! Also enforces a consistent polarization handedness.
    use iso_fortran_env, only: real64
    implicit none
    real(real64), intent(inout) :: V(3)
    real(real64), intent(inout) :: e(3,3)   ! columns are eigenvectors

    integer :: idx(3)
    real(real64) :: Vtmp(3), etmp(3,3)
    integer :: i

    !--- sort indices ascending by velocity
    idx = [(i, i=1,3)]
    call sort3(V, idx)   ! small helper (below)
    Vtmp = V(idx)
    etmp = e(:,idx)

    !--- enforce P = largest velocity = 3rd column
    V = (/ Vtmp(1), Vtmp(2), Vtmp(3) /)
    e = etmp

    !--- ensure right-handed coordinate system:
    if (dot_product(cross_product(e(:,1), e(:,2)), e(:,3)) < 0.0_real64) then
        e(:,1:2) = -e(:,1:2)
    end if
    

    end subroutine order_modes

    ! --------------------------------------------------------------------------
    subroutine sort3(arr, idx)
        real(real64), intent(in)  :: arr(3)
        integer,       intent(inout) :: idx(3)
        real(real64) :: tmp
        integer :: i,j,k,imin
        do i=1,2
            imin=i
            do j=i+1,3
                if (arr(idx(j)) < arr(idx(imin))) imin=j
            end do
            if (imin /= i) then
                k=idx(i); idx(i)=idx(imin); idx(imin)=k
            end if
        end do
    end subroutine sort3

    ! --------------------------------------------------------------------------
    pure function cross_product(a,b) result(c)
        real(real64), intent(in) :: a(3), b(3)
        real(real64) :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product

    ! --------------------------------------------------------------------------
    subroutine christoffel_modes(gamma, rho, n, V, e, q, use_precomputed, lam_in, e_in)
        ! Calculate the different modes and velocities for a given direction of
        ! propagation via the Christoffel matrix 
        
        implicit none 
        
        real(real64), intent(in) :: rho, n(3), C(6,6)
        logical, intent(in), optional :: use_precomputed
        real(real64), intent(out) :: e(3,3), q(6,3), V(3,3)
        real(real64), intent(in), optional :: lam_in(3), e_in(3,3)
        real(real64) :: lam(3), Z(3,3), gamma(3,3), nk(6) 
        
        ! Build the christoffel matris 
        nk = ( n(1)**2, n(2)**2, n(3)**2, sqrt(2) * n(2)*n(3), sqrt(2) * n(1)*n(3), sqrt(2) * n(1)*n(2) )
        
        if (present(use_precomputed) .and. use_precomputed .and. present(lam_in) .and. present(e_in)) then 
            lam = lam_in 
            e = e_in 
        else 
            call dsyevr_wrapper(gamma, lam, Z)
        end if 
        
        V = sqrt(max(lam, 0.0_real64) / rho )
        e = Z / spread( sqrt(sum(Z**2, dim=1)), 1,3) ! Normalize the eigenvectors/columns
        
        call order_modes(V,e) 
        call build_q(C, n, V, e, q)
        
    end subroutine christoffel_modes

    ! ==========================================================================
    elemental pure function ricker_s(t, fc, t0_opt) result(s)
        ! Zero phase Ricker wavelet 
        ! s(t) = (1-2 a^2) * exp(-a^2), a = pi*fc*(t-t0)
        real(real64), intent(in) :: t
        real(real64), intent(in) :: fc
        real(real64), intent(in), optional :: t0_opt        
        real(real64) :: s
        real(real64), parameter :: pi = acos(-1.0_real64) 
        real(real64) :: t0, a, e 
        
        
        if (present(t0_opt)) then 
            t0 = t0_opt
        else
            t0 = 1.0_real64/fc 
        endif 
        
        a = pi * fc * (t-t0) 
        e = exp(-a*a) 
        s = (1.0_real64 - 2.0_real64*a*a) * e
        
    end function ricker_s 
    
    ! ==========================================================================
    elemental pure function ricker_sdot(t, fc, t0_opt) result(sd)
        ! Time derivative of Ricker wavelet:
        ! ṡ(t) = 2 pi^2 fc^2 (t - t0) [ 2 pi^2 fc^2 (t - t0)^2 - 3 ] exp(-a^2)
        real(real64), intent(in) :: t
        real(real64), intent(in) :: fc
        real(real64), intent(in), optional :: t0_opt        
        real(real64) :: sd
        real(real64), parameter :: pi = acos(-1.0_real64) 
        real(real64) :: t0, dt, k, e 
        
        
        if (present(t0_opt)) then 
            t0 = t0_opt
        else
            t0 = 1.0_real64/fc 
        endif 
        
        dt = (t-t0) 
        k = pi * pi * fc * fc 
        e = exp(-k*dt*dt) 
        sd = 2.0_real64*k*dt * (2.0_real64*k*dt*dt - 3.0_real64) * e
        
    end function ricker_sdot
    
    
    ! ==========================================================================
    elemental pure function gaussian_s(t, fc, t0_opt, bw_opt) result(s)
        ! Baseband Gaussian: s(t) = exp(-((t - t0)/bw)^2)
        ! If bw (bandwidth parameter) not given, derive from fc (choose ONE 
        ! mapping below).
        real(real64), intent(in) :: t, fc
        real(real64), intent(in), optional :: t0_opt, bw_opt 
        real(real64) :: s 
        real(real64) :: t0, bw, a 
        real(real64), parameter :: pi = acos(-1.0_real64) 
        
        if (present(t0_opt)) then  
            t0 = t0_opt 
        else
            t0 = 1.d0/fc
        endif 
            
        if (present(bw_opt)) then
            bw = bw_opt
        else
            bw = 0.2_real64/ ( pi * fc) 
        endif 
        
        if (bw <= 0.0_real64) bw = tiny(1.0_real64) 
    
        a = (t-t0) / bw
        s = exp(-a*a)
        
    end function gaussian_s
    
    ! ==========================================================================
    elemental pure function gaussian_sdot(t, fc, t0_opt, bw_opt) result(sd)
        ! Derivative: ṡ(t) = -(2/bw^2) (t - t0) exp(-(t - t0)^2 / bw^2)
        real(real64), intent(in) :: t, fc
        real(real64), intent(in), optional :: t0_opt, bw_opt 
        real(real64) :: sd 
        real(real64) :: t0, bw, dt
        real(real64), parameter :: pi = acos(-1.0_real64) 
        
        if (present(bw_opt)) then
            bw = bw_opt
        else
            bw = 0.2_real64/(pi * fc) 
        endif 
        
        if (present(t0_opt)) then  
            t0 = t0_opt 
        else
            t0 = 1.d0/fc
        endif 

        dt = (t-t0)
        invbw2 = 1.0_real64 / (bw * bw)
        sd = -2.0_real64*dt*invbw2 * exp( -(dt*dt)*invbw2 )
        
    end function gaussian_sdot

    ! ==========================================================================
    ! function klauder 
        
    !     s = exp(-(t-t0)**2 / (tau**2) ) * sin(2*pi * fc *(t-t0) )
    ! end function klauder

end program source