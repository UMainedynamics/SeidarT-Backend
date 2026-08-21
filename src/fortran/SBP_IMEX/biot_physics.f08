module biot_kernels
    use iso_fortran_env, only: real64
    implicit none
    public :: explicit_constitutive_kernel, implicit_drag_kernel
    
    contains
    
    ! --------------------------------------------------------------------------
    subroutine explicit_constitutive_kernel(nx, ny, nz, Q, &
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
        
    end subroutine explicit_constitutive_kernel
        
    ! --------------------------------------------------------------------------
    subroutine implicit_drag_kernel(nx, ny, nz, Q_star, Q_out, &
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
                        wx_star = Q_star(4, i,j,k) - Q_star(1, i,j,k)
                        wy_star = Q_star(5, i,j,k) - Q_star(2, i,j,k)
                        wz_star = Q_star(6, i,j,k) - Q_star(3, i,j,k)
                        
                        ! Solve M * w = w*
                        wx = inv_det * (    (m22*m33 - m23*m23) * wx_star + &
                                            (m13*m23 - m12*m33) * wy_star + &
                                            (m12*m23 - m13*m22) * wz_star )
                        
                        wy = inv_det * (    (m13*m23 - m12*m33) * wx_star + &
                                            (m11*m33 - m13*m13) * wy_star + &
                                            (m12*m13 - m11*m23) * wz_star )
                                        
                        wz = inv_det * (    (m12*m23 - m13*m22) * wx_star + &
                                            (m12*m13 - m11*m23) * wy_star + &
                                            (m11*m22 - m12*m12) * wz_star )
                        
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
        
    end subroutine implicit_drag_kernel
    
end module biot_kernels