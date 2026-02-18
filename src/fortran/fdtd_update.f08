module fdtd_update

    ! --------------------------------------------------------------------------
    ! Viscoelastic updates
    subroutine update_sigma2(
            sigma, c1, c2, c3, gamma, 
            dvx_dx, dvz_dz, dvx_dz, dvz_dx, dt
        ) 

        use iso_fortran_env, only: real64
        implicit none 

        real(real64), intent(inout) :: sigma 
        real(real64), intent(in) :: c1, c2, c3
        real(real64), intent(in) :: gamma
        real(real64), intent(in) :: dvx_dx, dvz_dz, dvx_dz, dvz_dx
        real(real64), intent(in) :: dt

        sigma = ( sigma + &
            (   c1 * dvx_dx + c2 * dvz_dz + c3 * (dvz_dx + dvx_dz)  ) * dt ) / & 
                    (1.0_real64 + gamma * dt )

    end subroutine update_sigma2

    subroutine update_velocity2(
            v, dsigx1, dsigx2, rho, dt 
        )

        use iso_fortran_env, only: real64
        implicit none 

        real(real64), intent(inout) :: v 
        real(real64), intent(in) :: dsig1, dsig2, 
        real(real64), intent(in) :: rho 
        real(real64), intent(in) :: dt 
        
        v = v + dt * (dsig1 + dsig2) / rho

    end subroutine update_velocity2

    ! --------------------------------------------------------------------------
    ! Poroelastic updates
    subroutine update_sigma2_poro(
        sigma, p, 
        c1, c2, c3,
        gamma,
        alpha_comp, M,
        dvx_dx, dvz_dz, dvx_dz, dvz_dx,
        div_vf, dt
    )
    real(real64), intent(inout) :: sigma      ! stress component
    real(real64), intent(inout) :: p          ! pore pressure
    
    real(real64), intent(in)    :: c1, c2, c3
    real(real64), intent(in)    :: gamma
    real(real64), intent(in)    :: alpha_comp
    real(real64), intent(in)    :: M

    real(real64), intent(in)    :: dvx_dx, dvz_dz
    real(real64), intent(in)    :: dvx_dz, dvz_dx
    real(real64), intent(in)    :: div_vf       ! ∇·v_f
    real(real64), intent(in)    :: dt
    
    real(real64) :: strain_rate 
    real(real64) :: dpdt 
    

    strain_rate = dvx_dx + dvz_dz 
    
    dpdt = M * (-div_vf - alpha_comp * strain_rate)

    end subroutine update_sigma2_poro

    !
    subroutine update_velocity2_poro 
    
    end subroutine update_velocity2_poro

end module fdtd_update