module pseudospectral_stencils
    use iso_fortran_env, only: real64
    implicit none
    public :: explicit_constitutive_kernel, implicit_drag_kernel
    
    contains
    
    ! --------------------------------------------------------------------------
    subroutine explicit_constitutive_kernel3_rk4(nx, ny, nz, Q, &
                                            dvx_dx, dvx_dy, dvx_dz, &
                                            dvy_dx, dvy_dy, dvy_dz, &
                                            dvz_dx, dvz_dy, dvz_dz, &
                                            dsxx_dx, dsyy_dy, dszz_dz, &
                                            dsyz_dy, dsyz_dz, &
                                            dsxz_dx, dsxz_dz, &
                                            dsxy_dx, dsxy_dy, &
                                            dp_dx,   dp_dy,   dp_dz,   &
                                            C, gamma_visco, density_s, density_f, &
                                            phi, lwc, bulk_mod_fluid, fluid_pressure_dot, T)
    
        integer, intent(in) :: nx, ny, nz 
        real(real64), intent(in)  :: Q(21, nx, ny, nz)
        real(real64), intent(in)  :: dvx_dx(nx,ny,nz), dvx_dy(nx,ny,nz), dvx_dz(nx,ny,nz)
        real(real64), intent(in)  :: dvy_dx(nx,ny,nz), dvy_dy(nx,ny,nz), dvy_dz(nx,ny,nz)
        real(real64), intent(in)  :: dvz_dx(nx,ny,nz), dvz_dy(nx,ny,nz), dvz_dz(nx,ny,nz)
        real(real64), intent(in)  :: dsxx_dx(nx,ny,nz), dsyy_dy(nx,ny,nz), dszz_dz(nx,ny,nz)
        real(real64), intent(in)  :: dsyz_dy(nx,ny,nz), dsyz_dz(nx,ny,nz)
        real(real64), intent(in)  :: dsxz_dx(nx,ny,nz), dsxz_dz(nx,ny,nz)
        real(real64), intent(in)  :: dsxy_dx(nx,ny,nz), dsxy_dy(nx,ny,nz)
        real(real64), intent(in)  :: dp_dx(nx,ny,nz),   dp_dy(nx,ny,nz),   dp_dz(nx,ny,nz)
        real(real64), intent(in)  :: C(21, nx, ny, nz)
        real(real64), intent(in)  :: gamma_visco(6, nx, ny, nz)
        real(real64), intent(in)  :: density_s(nx, ny, nz), density_f(nx, ny, nz)
        real(real64), intent(in)  :: phi(nx, ny, nz), lwc(nx, ny, nz)
        real(real64), intent(in)  :: bulk_mod_fluid
        real(real64), intent(out) :: fluid_pressure_dot(nx, ny, nz)
        real(real64), intent(out) :: T(21, nx, ny, nz)
        
        ! local variables 
        integer :: i, j, k
        real(real64) :: exx, eyy, ezz, eyz, exz, exy 
        real(real64) :: gyz, gxz, gxy, eff_rho_f
        
        !$omp target teams distribute parallel do collapse(3) 
        do k = 1,nz 
            do j = 1,ny
                do i = 1,nx
                    exx = dvx_dx(i,j,k)
                    eyy = dvy_dy(i,j,k)
                    ezz = dvz_dz(i,j,k)
                    eyz = 0.5_real64 * (dvy_dz(i,j,k) + dvz_dy(i,j,k) )
                    exz = 0.5_real64 * (dvx_dz(i,j,k) + dvz_dx(i,j,k) )
                    exy = 0.5_real64 * (dvx_dy(i,j,k) + dvy_dx(i,j,k) )
                    gyz = 2.0_real64 * eyz
                    gxz = 2.0_real64 * exz
                    gxy = 2.0_real64 * exy
                    
                    ! Memory Variable 
                    T(16,i,j,k) = gamma_visco(1, i,j,k) * exx
                    T(17,i,j,k) = gamma_visco(2, i,j,k) * eyy
                    T(18,i,j,k) = gamma_visco(3, i,j,k) * ezz
                    T(19,i,j,k) = gamma_visco(4, i,j,k) * eyz
                    T(20,i,j,k) = gamma_visco(5, i,j,k) * exz
                    T(21,i,j,k) = gamma_visco(6, i,j,k) * exy
                    
                    ! Stress Rates 
                    !sxx
                    T(7, i,j,k) =  (  C(1,i,j,k) * exx +  C(2, i,j,k) * eyy +  &
                                      C(3,i,j,k) * ezz +  C(4, i,j,k) * gyz +  &
                                      C(5,i,j,k) * gxz +  C(6, i,j,k) * gxy ) - Q(16,i,j,k)
                    T(8, i,j,k) =  (  C(2,i,j,k) * exx +  C(7, i,j,k) * eyy +  &
                                      C(8,i,j,k) * ezz +  C(9, i,j,k) * gyz + &
                                     C(10,i,j,k) * gxz + C(11, i,j,k) * gxy ) - Q(17,i,j,k)
                    T(9, i,j,k) =  (  C(3,i,j,k) * exx +  C(8, i,j,k) * eyy + &
                                     C(12,i,j,k) * ezz + C(13, i,j,k) * gyz + &
                                     C(14,i,j,k) * gxz + C(15, i,j,k) * gxy ) - Q(18,i,j,k)
                    T(10, i,j,k) = (  C(4,i,j,k) * exx +  C(9, i,j,k) * eyy + &
                                     C(13,i,j,k) * ezz + C(16, i,j,k) * gyz + &
                                     C(17,i,j,k) * gxz + C(18, i,j,k) * gxy ) - Q(19,i,j,k)
                    T(11, i,j,k) = (  C(5,i,j,k) * exx + C(10, i,j,k) * eyy + &
                                     C(14,i,j,k) * ezz + C(17, i,j,k) * gyz + &
                                     C(19,i,j,k) * gxz + C(20, i,j,k) * gxy ) - Q(20,i,j,k)
                    T(12, i,j,k) = (  C(6,i,j,k) * exx + C(11, i,j,k) * eyy + &
                                     C(15,i,j,k) * ezz + C(18, i,j,k) * gyz + &
                                     C(20,i,j,k) * gxz + C(21, i,j,k) * gxy ) - Q(21,i,j,k)
                    
                    T(13:15,i,j,k) = 0.0_real64
                    
                    ! Solid Acceleration
                    T(1,i,j,k) = ( dsxx_dx(i,j,k) + dsxy_dy(i,j,k) + dsxz_dz(i,j,k) ) / density_s(i,j,k)
                    T(2,i,j,k) = ( dsxy_dx(i,j,k) + dsyy_dy(i,j,k) + dsyz_dz(i,j,k) ) / density_s(i,j,k)
                    T(3,i,j,k) = ( dsxz_dx(i,j,k) + dsyz_dy(i,j,k) + dszz_dz(i,j,k) ) / density_s(i,j,k)
                    
                    ! Fluid Acceleration
                    eff_rho_f = max(density_f(i,j,k) * lwc(i,j,k), 1.0e-12_real64)
                    
                    T(4,i,j,k) = - dp_dx(i,j,k) / eff_rho_f
                    T(5,i,j,k) = - dp_dy(i,j,k) / eff_rho_f
                    T(6,i,j,k) = - dp_dz(i,j,k) / eff_rho_f
                    
                    ! Fluid Pressure Rate; Volumetric Divergence 
                    fluid_pressure_dot(i,j,k) = - ( bulk_mod_fluid / max(phi(i,j,k), 1.0e-5_real64) ) * &
                                                (exx + eyy + ezz)
                end do
            end do
        end do
        
    end subroutine explicit_constitutive_kernel3_rk4
        
    ! --------------------------------------------------------------------------
    subroutine implicit_drag_kernel3(nx, ny, nz, Q_star, Q_out, &
                                    drag_tensor, tau_relax, &
                                    density_s, density_f, lwc, dt_gamma)
        
        integer, intent(in) :: nx, ny, nz
        real(real64), intent(in) :: Q_star(21, nx, ny, nz)
        real(real64), intent(out) :: Q_out(21, nx, ny, nz)
        real(real64), intent(in) :: drag_tensor(6, nx, ny, nz), tau_relax(6, nx, ny, nz)
        real(real64), intent(in) :: density_s(nx, ny, nz), density_f(nx, ny, nz), lwc(nx, ny, nz)
        real(real64), intent(in) :: dt_gamma
        
        ! local variables
        integer :: i,j,k,l
        real(real64) :: eff_rho_f, beta, alpha, bxx, byy, bzz, byz, bxz, bxy
        real(real64) :: m11, m12, m13, m22, m23, m33, det, inv_det
        real(real64) :: wx_s, wy_s, wz_s, wx, wy, wz, bw_x, bw_y, bw_z
        
        ! ------------
        alpha = dt_gamma
        
        Q_out(7:15,:,:,:) = Q_star(7:15,:,:,:)
    
        !$omp target teams distribute parallel do collapse(3)
        do k = 1,nz 
            do j = 1,ny 
                do i = 1,nx 
                    ! Vicsoelastic Memory Relaxation 
                    do l = 1,6
                        Q_out(15+l, i,j,k) = Q_star(15+l, i,j,k) / &
                        (1.0_real64 + alpha / max(tau_relax(l,i,j,k), 1.0e-12_real64))
                    end do
                    
                    ! Anisotroopic fluid-solid velocity drag 
                    eff_rho_f = max(density_f(i,j,k) * lwc(i,j,k), 1.0e-12_real64) 
                    
                    bxx = drag_tensor(1, i, j, k)
                    byy = drag_tensor(2, i, j, k)
                    bzz = drag_tensor(3, i, j, k)
                    byz = drag_tensor(4, i, j, k)
                    bxz = drag_tensor(5, i, j, k)
                    bxy = drag_tensor(6, i, j, k)
                    
                    if ( (bxx + byy + bzz) <= 1.0e-16_real64 ) then 
                        Q_out(1:6, i,j,k) = Q_star(1:6, i,j,k)
                    else
                        beta = ( density_s(i,j,k) + eff_rho_f ) / ( density_s(i,j,k) * eff_rho_f )
                        
                        ! Construct symmetric 3x3 matrix M = I_3 + alpha * beta * b 
                        m11 = 1.0_real64 + alpha * beta * bxx
                        m22 = 1.0_real64 + alpha * beta * byy
                        m33 = 1.0_real64 + alpha * beta * bzz
                        m12 = alpha * beta * bxy
                        m13 = alpha * beta * bxz
                        m23 = alpha * beta * byz
                        
                        ! Compute determinant and inverse of M
                        det =   m11 * (m22 * m33 - m23 * m23) - &
                                m12 * (m12 * m33 - m13 * m23) + &
                                m13 * (m12 * m23 - m13 * m22)
                        
                        inv_det = 1.0_real64 / det
                        
                        ! Relative velocity w* = v_f* - v_s*
                        wx_s = Q_star(4, i,j,k) - Q_star(1, i,j,k)
                        wy_s = Q_star(5, i,j,k) - Q_star(2, i,j,k)
                        wz_s = Q_star(6, i,j,k) - Q_star(3, i,j,k)
                        
                        ! Solve M * w = w*
                        wx = inv_det * (    (m22*m33 - m23*m23) * wx_s + &
                                            (m13*m23 - m12*m33) * wy_s + &
                                            (m12*m23 - m13*m22) * wz_s )
                        
                        wy = inv_det * (    (m13*m23 - m12*m33) * wx_s + &
                                            (m11*m33 - m13*m13) * wy_s + &
                                            (m12*m13 - m11*m23) * wz_s )
                                        
                        wz = inv_det * (    (m12*m23 - m13*m22) * wx_s + &
                                            (m12*m13 - m11*m23) * wy_s + &
                                            (m11*m22 - m12*m12) * wz_s )
                        
                        ! Compute b * w vector 
                        bw_x = bxx * wx + bxy * wy + bxz * wz
                        bw_y = bxy * wx + byy * wy + byz * wz
                        bw_z = bxz * wx + byz * wy + bzz * wz 
                        
                        ! Update solid velocities
                        Q_out(1, i,j,k) = Q_star(1, i,j,k) + (alpha / density_s(i,j,k) ) * bw_x
                        Q_out(2, i,j,k) = Q_star(2, i,j,k) + (alpha / density_s(i,j,k) ) * bw_y
                        Q_out(3, i,j,k) = Q_star(3, i,j,k) + (alpha / density_s(i,j,k) ) * bw_z
                        
                        ! Update fluid velocities
                        Q_out(4, i,j,k) = Q_star(4, i,j,k) - (alpha / eff_rho_f ) * bw_x
                        Q_out(5, i,j,k) = Q_star(5, i,j,k) - (alpha / eff_rho_f ) * bw_y
                        Q_out(6, i,j,k) = Q_star(6, i,j,k) - (alpha / eff_rho_f ) * bw_z
                    end if
                end do
            end do
        end do
        
    end subroutine implicit_drag_kernel3
    
    ! --------------------------------------------------------------------------
    subroutine explicit_em_kernel3_rk4( &
        nx, ny, nz, &
        Ex, Ey, Ez, Hx, Hy, Hz, &
        dEz_dy, dEy_dz, dEx_dz, dEz_dx, dEy_dx, dEx_dy, &
        dHz_dy, dHy_dz, dHz_dx, dHx_dz, dHy_dx, dHx_dy, &
        aEx, bEx, cEx, aEy, bEy, cEy, aEz, bEz, cEz, &
        sig11, sig12, sig13, sig22, sig23, sig33, &
        mu_inv_x, mu_inv_y, mu_inv_z, &
        dEx_dt, dEy_dt, dEz_dt, dHx_dt, dHy_dt, dHz_dt )

        use, intrinsic :: iso_fortran_env, only : real64
        implicit none

        integer, intent(in) :: nx, ny, nz

        real(real64), intent(in) :: Ex(nx,ny,nz), Ey(nx,ny,nz), Ez(nx,ny,nz)
        real(real64), intent(in) :: Hx(nx,ny,nz), Hy(nx,ny,nz), Hz(nx,ny,nz)

        real(real64), intent(in) :: dEz_dy(nx,ny,nz), dEy_dz(nx,ny,nz)
        real(real64), intent(in) :: dEx_dz(nx,ny,nz), dEz_dx(nx,ny,nz)
        real(real64), intent(in) :: dEy_dx(nx,ny,nz), dEx_dy(nx,ny,nz)

        real(real64), intent(in) :: dHz_dy(nx,ny,nz), dHy_dz(nx,ny,nz)
        real(real64), intent(in) :: dHz_dx(nx,ny,nz), dHx_dz(nx,ny,nz)
        real(real64), intent(in) :: dHy_dx(nx,ny,nz), dHx_dy(nx,ny,nz)

        real(real64), intent(in) :: aEx(nx,ny,nz), bEx(nx,ny,nz), cEx(nx,ny,nz)
        real(real64), intent(in) :: aEy(nx,ny,nz), bEy(nx,ny,nz), cEy(nx,ny,nz)
        real(real64), intent(in) :: aEz(nx,ny,nz), bEz(nx,ny,nz), cEz(nx,ny,nz)

        real(real64), intent(in) :: sig11(nx,ny,nz), sig12(nx,ny,nz)
        real(real64), intent(in) :: sig13(nx,ny,nz), sig22(nx,ny,nz)
        real(real64), intent(in) :: sig23(nx,ny,nz), sig33(nx,ny,nz)

        real(real64), intent(in) :: mu_inv_x(nx,ny,nz)
        real(real64), intent(in) :: mu_inv_y(nx,ny,nz)
        real(real64), intent(in) :: mu_inv_z(nx,ny,nz)

        real(real64), intent(out) :: dEx_dt(nx,ny,nz), dEy_dt(nx,ny,nz), dEz_dt(nx,ny,nz)
        real(real64), intent(out) :: dHx_dt(nx,ny,nz), dHy_dt(nx,ny,nz), dHz_dt(nx,ny,nz)

        integer :: i, j, k
        real(real64) :: rhs_x, rhs_y, rhs_z
        real(real64) :: curl_e_x, curl_e_y, curl_e_z

#ifdef SEIDART_OPENMP_GPU
    !$omp target teams loop collapse(3) &
    !$omp& private(rhs_x, rhs_y, rhs_z, curl_e_x, curl_e_y, curl_e_z)
#else
    !$omp parallel do collapse(3) schedule(static) &
    !$omp& private(i, j, k, rhs_x, rhs_y, rhs_z, curl_e_x, curl_e_y, curl_e_z)
#endif
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    rhs_x = dHz_dy(i,j,k) - dHy_dz(i,j,k) - &
                        (sig11(i,j,k) * Ex(i,j,k) + &
                        sig12(i,j,k) * Ey(i,j,k) + &
                        sig13(i,j,k) * Ez(i,j,k))

                    rhs_y = dHx_dz(i,j,k) - dHz_dx(i,j,k) - &
                        (sig12(i,j,k) * Ex(i,j,k) + &
                        sig22(i,j,k) * Ey(i,j,k) + &
                        sig23(i,j,k) * Ez(i,j,k))

                    rhs_z = dHy_dx(i,j,k) - dHx_dy(i,j,k) - &
                        (sig13(i,j,k) * Ex(i,j,k) + &
                        sig23(i,j,k) * Ey(i,j,k) + &
                        sig33(i,j,k) * Ez(i,j,k))

                    dEx_dt(i,j,k) = aEx(i,j,k) * rhs_x + &
                                    bEx(i,j,k) * rhs_y + &
                                    cEx(i,j,k) * rhs_z

                    dEy_dt(i,j,k) = aEy(i,j,k) * rhs_x + &
                                    bEy(i,j,k) * rhs_y + &
                                    cEy(i,j,k) * rhs_z

                    dEz_dt(i,j,k) = aEz(i,j,k) * rhs_x + &
                                    bEz(i,j,k) * rhs_y + &
                                    cEz(i,j,k) * rhs_z

                    curl_e_x = dEz_dy(i,j,k) - dEy_dz(i,j,k)
                    curl_e_y = dEx_dz(i,j,k) - dEz_dx(i,j,k)
                    curl_e_z = dEy_dx(i,j,k) - dEx_dy(i,j,k)

                    dHx_dt(i,j,k) = -mu_inv_x(i,j,k) * curl_e_x
                    dHy_dt(i,j,k) = -mu_inv_y(i,j,k) * curl_e_y
                    dHz_dt(i,j,k) = -mu_inv_z(i,j,k) * curl_e_z

                end do
            end do
        end do
#ifndef SEIDART_OPENMP_GPU
        !$omp end parallel do
#endif

    end subroutine explicit_em_kernel3_rk4
    
    ! --------------------------------------------------------------------------
    subroutine implicit_em_conduction( &
        nx, ny, nz, &
        Ex_star, Ey_star, Ez_star, &
        Hx_star, Hy_star, Hz_star, &
        Ex_out, Ey_out, Ez_out, &
        Hx_out, Hy_out, Hz_out, &
        aEx, bEx, cEx, &
        aEy, bEy, cEy, &
        aEz, bEz, cEz, &
        sig11, sig12, sig13, sig22, sig23, sig33, &
        dt_gamma )

        use, intrinsic :: iso_fortran_env, only : real64
        implicit none

        integer, intent(in) :: nx, ny, nz
        real(real64), intent(in) :: dt_gamma

        ! IMEX predictor state U_star.
        real(real64), intent(in) :: Ex_star(nx,ny,nz)
        real(real64), intent(in) :: Ey_star(nx,ny,nz)
        real(real64), intent(in) :: Ez_star(nx,ny,nz)

        real(real64), intent(in) :: Hx_star(nx,ny,nz)
        real(real64), intent(in) :: Hy_star(nx,ny,nz)
        real(real64), intent(in) :: Hz_star(nx,ny,nz)

        ! Implicit stage solution U_stage.
        real(real64), intent(out) :: Ex_out(nx,ny,nz)
        real(real64), intent(out) :: Ey_out(nx,ny,nz)
        real(real64), intent(out) :: Ez_out(nx,ny,nz)

        real(real64), intent(out) :: Hx_out(nx,ny,nz)
        real(real64), intent(out) :: Hy_out(nx,ny,nz)
        real(real64), intent(out) :: Hz_out(nx,ny,nz)

        ! Entries of epsilon^{-1}.
        real(real64), intent(in) :: aEx(nx,ny,nz), bEx(nx,ny,nz), cEx(nx,ny,nz)
        real(real64), intent(in) :: aEy(nx,ny,nz), bEy(nx,ny,nz), cEy(nx,ny,nz)
        real(real64), intent(in) :: aEz(nx,ny,nz), bEz(nx,ny,nz), cEz(nx,ny,nz)

        ! Symmetric electrical-conductivity tensor sigma.
        real(real64), intent(in) :: sig11(nx,ny,nz), sig12(nx,ny,nz)
        real(real64), intent(in) :: sig13(nx,ny,nz), sig22(nx,ny,nz)
        real(real64), intent(in) :: sig23(nx,ny,nz), sig33(nx,ny,nz)

        integer :: i, j, k

        real(real64) :: p11, p12, p13, p21, p22, p23, p31, p32, p33
        real(real64) :: m11, m12, m13, m21, m22, m23, m31, m32, m33
        real(real64) :: det, inv_det
        real(real64) :: x11, x12, x13, x21, x22, x23, x31, x32, x33

#ifdef SEIDART_OPENMP_GPU
    !$omp target teams loop collapse(3) &
    !$omp& private(p11, p12, p13, p21, p22, p23, p31, p32, p33, &
    !$omp&         m11, m12, m13, m21, m22, m23, m31, m32, m33, &
    !$omp&         x11, x12, x13, x21, x22, x23, x31, x32, x33, det, inv_det)
#else
    !$omp parallel do collapse(3) schedule(static) &
    !$omp& private(i, j, k, p11, p12, p13, p21, p22, p23, p31, p32, p33, &
    !$omp&         m11, m12, m13, m21, m22, m23, m31, m32, m33, &
    !$omp&         x11, x12, x13, x21, x22, x23, x31, x32, x33, det, inv_det)
#endif
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    ! ------------------------------------------------------------
                    ! P = epsilon^{-1} * sigma
                    ! ------------------------------------------------------------
                    p11 = aEx(i,j,k) * sig11(i,j,k) + &
                        bEx(i,j,k) * sig12(i,j,k) + &
                        cEx(i,j,k) * sig13(i,j,k)

                    p12 = aEx(i,j,k) * sig12(i,j,k) + &
                        bEx(i,j,k) * sig22(i,j,k) + &
                        cEx(i,j,k) * sig23(i,j,k)

                    p13 = aEx(i,j,k) * sig13(i,j,k) + &
                        bEx(i,j,k) * sig23(i,j,k) + &
                        cEx(i,j,k) * sig33(i,j,k)

                    p21 = aEy(i,j,k) * sig11(i,j,k) + &
                        bEy(i,j,k) * sig12(i,j,k) + &
                        cEy(i,j,k) * sig13(i,j,k)

                    p22 = aEy(i,j,k) * sig12(i,j,k) + &
                        bEy(i,j,k) * sig22(i,j,k) + &
                        cEy(i,j,k) * sig23(i,j,k)

                    p23 = aEy(i,j,k) * sig13(i,j,k) + &
                        bEy(i,j,k) * sig23(i,j,k) + &
                        cEy(i,j,k) * sig33(i,j,k)

                    p31 = aEz(i,j,k) * sig11(i,j,k) + &
                        bEz(i,j,k) * sig12(i,j,k) + &
                        cEz(i,j,k) * sig13(i,j,k)

                    p32 = aEz(i,j,k) * sig12(i,j,k) + &
                        bEz(i,j,k) * sig22(i,j,k) + &
                        cEz(i,j,k) * sig23(i,j,k)

                    p33 = aEz(i,j,k) * sig13(i,j,k) + &
                        bEz(i,j,k) * sig23(i,j,k) + &
                        cEz(i,j,k) * sig33(i,j,k)

                    ! ------------------------------------------------------------
                    ! M = I + dt_gamma * P
                    ! ------------------------------------------------------------
                    m11 = 1.0_real64 + dt_gamma * p11
                    m12 =               dt_gamma * p12
                    m13 =               dt_gamma * p13

                    m21 =               dt_gamma * p21
                    m22 = 1.0_real64 + dt_gamma * p22
                    m23 =               dt_gamma * p23

                    m31 =               dt_gamma * p31
                    m32 =               dt_gamma * p32
                    m33 = 1.0_real64 + dt_gamma * p33

                    ! ------------------------------------------------------------
                    ! det(M)
                    ! ------------------------------------------------------------
                    det = m11 * (m22*m33 - m23*m32) - &
                        m12 * (m21*m33 - m23*m31) + &
                        m13 * (m21*m32 - m22*m31)

                    ! A physically valid positive epsilon and positive-semidefinite
                    ! conductivity should make M nonsingular for dt_gamma >= 0.
                    inv_det = 1.0_real64 / det

                    ! ------------------------------------------------------------
                    ! adjugate(M), evaluated explicitly for a local 3x3 solve.
                    ! X = M^{-1} * E_star
                    ! ------------------------------------------------------------
                    x11 =  (m22*m33 - m23*m32) * inv_det
                    x12 = -(m12*m33 - m13*m32) * inv_det
                    x13 =  (m12*m23 - m13*m22) * inv_det

                    x21 = -(m21*m33 - m23*m31) * inv_det
                    x22 =  (m11*m33 - m13*m31) * inv_det
                    x23 = -(m11*m23 - m13*m21) * inv_det

                    x31 =  (m21*m32 - m22*m31) * inv_det
                    x32 = -(m11*m32 - m12*m31) * inv_det
                    x33 =  (m11*m22 - m12*m21) * inv_det

                    Ex_out(i,j,k) = x11 * Ex_star(i,j,k) + &
                                    x12 * Ey_star(i,j,k) + &
                                    x13 * Ez_star(i,j,k)

                    Ey_out(i,j,k) = x21 * Ex_star(i,j,k) + &
                                    x22 * Ey_star(i,j,k) + &
                                    x23 * Ez_star(i,j,k)

                    Ez_out(i,j,k) = x31 * Ex_star(i,j,k) + &
                                    x32 * Ey_star(i,j,k) + &
                                    x33 * Ez_star(i,j,k)

                    ! There is no Ohmic conduction term in Faraday's law.
                    Hx_out(i,j,k) = Hx_star(i,j,k)
                    Hy_out(i,j,k) = Hy_star(i,j,k)
                    Hz_out(i,j,k) = Hz_star(i,j,k)

                end do
            end do
        end do
#ifndef SEIDART_OPENMP_GPU
    !$omp end parallel do
#endif

end subroutine implicit_em_conduction

    ! --------------------------------------------------------------------------
    subroutine electrokinetic_current_kernel( &
        nx, ny, nz, vsx, vsy, vsz, vfx, vfy, vfz, &
        lek, jekx, jeky, jekz )

        use, intrinsic :: iso_fortran_env, only : real64
        implicit none

        integer, intent(in) :: nx, ny, nz

        real(real64), intent(in) :: vsx(nx,ny,nz), vsy(nx,ny,nz), vsz(nx,ny,nz)
        real(real64), intent(in) :: vfx(nx,ny,nz), vfy(nx,ny,nz), vfz(nx,ny,nz)

        real(real64), intent(in) :: lek(6,nx,ny,nz)

        real(real64), intent(out) :: jekx(nx,ny,nz)
        real(real64), intent(out) :: jeky(nx,ny,nz)
        real(real64), intent(out) :: jekz(nx,ny,nz)

        integer :: i, j, k
        real(real64) :: wx, wy, wz
        
#ifdef SEIDART_OPENMP_GPU
    !$omp target teams loop collapse(3) private(wx, wy, wz)
#else
    !$omp parallel do collapse(3) schedule(static) private(i, j, k, wx, wy, wz)
#endif
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    wx = vfx(i,j,k) - vsx(i,j,k)
                    wy = vfy(i,j,k) - vsy(i,j,k)
                    wz = vfz(i,j,k) - vsz(i,j,k)

                    jekx(i,j,k) = lek(1,i,j,k) * wx + &
                                lek(6,i,j,k) * wy + &
                                lek(5,i,j,k) * wz

                    jeky(i,j,k) = lek(6,i,j,k) * wx + &
                                lek(2,i,j,k) * wy + &
                                lek(4,i,j,k) * wz

                    jekz(i,j,k) = lek(5,i,j,k) * wx + &
                                lek(4,i,j,k) * wy + &
                                lek(3,i,j,k) * wz

                end do
            end do
        end do
#ifndef SEIDART_OPENMP_GPU
    !$omp end parallel do
#endif

end subroutine electrokinetic_current_kernel
    
end module pseudospectral_stencils