module cpmlfdtd 
    
    use seidartio
    use seidart_types
    
    implicit none 
    
    
    contains
    
    subroutine update_elastic2(i, j, dt, c1, c2, c3, gamma, sigmaxx, vx, vz)
        implicit none
        integer, intent(in) :: i, j
        real(8), intent(in) :: dt, c11, c13, c15, gamma_x
        real(8), intent(inout) :: sigmaxx, vx, vz

        sigmaxx = ( sigmaxx + 
                    (c1 * value_dvx_dx + 
                    c2 * value_dvz_dz + 
                    c3 * (value_dvz_dx + value_dvx_dz)) * dt ) /  
                    (1 + gamma * dt )
    end subroutine

    subroutine update_acoustic2(i, j, dt, rho, pressure)
        implicit none
        integer, intent(in) :: i, j
        real(8), intent(in) :: dt, rho
        real(8), intent(inout) :: pressure

        pressure = pressure + dt * (dp_dx + dp_dz) / rho
    end subroutine

    ! =========================================================================
    subroutine model2(domain, source, material_geometry, SINGLE_OUTPUT)
        use constants
        use omp_lib ! Include the OpenMP library 
        
        ! Input arguments
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT
        integer, dimension(nx, nz) :: material_geometry
        
        ! Local variables
        real(real64) :: deltarho, velocnorm, value_dvx_dx, value_dvx_dz, &
            value_dvz_dx, value_dvz_dz, value_dsigmaxx_dx, value_dsigmazz_dz, &
            value_dsigmaxz_dx, value_dsigmaxz_dz

        ! 1D arrays for damping profiles
        real(real64), allocatable :: c11(:,:), c13(:,:), c15(:,:), c33(:,:), c35(:,:), c55(:,:), rho(:,:)
        real(real64), allocatable :: K_x(:), alpha_x(:), a_x(:), b_x(:), K_x_half(:), alpha_x_half(:), a_x_half(:), b_x_half(:)
        real(real64), allocatable :: K_z(:), alpha_z(:), a_z(:), b_z(:), K_z_half(:), alpha_z_half(:), a_z_half(:), b_z_half(:)
        real(real64), allocatable :: gamma_x(:,:), gamma_z(:,:), gamma_xz(:,:)        
        real(real64), allocatable :: srcx(:), srcz(:) ! The vector time series of the source
        
        ! Model variables
        real(real64), allocatable :: vx(:,:),vz(:,:),sigmaxx(:,:),sigmazz(:,:),sigmaxz(:,:)
        real(real64), allocatable :: pressure(:,:)
        real(real64), allocatable :: memory_dvx_dx(:,:), memory_dvx_dz(:,:), &
                                     memory_dvz_dx(:,:), memory_dvz_dz(:,:), &
                                     memory_dsigmaxx_dx(:,:), memory_dsigmazz_dz(:,:), &
                                     memory_dsigmaxz_dx(:,:), memory_dsigmaxz_dz(:,:)
        
        integer :: nx, nz 
        real(real64) :: dx, dz, dt 
        
        integer :: i, j, it, isource, jsource
        logical :: SINGLE

        ! The default data output is single precision unless SINGLE_OUTPUT is 
        ! set to .FALSE.
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        end if
        
        nx = domain%nx 
        nz = domain%nz 
        dt = source%dt 
        dx = domain%dx 
        dz = domain%dz 
        
        ! Allocate the arrays based on runtime values of nx and nz
        allocate(c11(nx, nz), c13(nx, nz), c15(nx, nz), &
                 c33(nx, nz), c35(nx, nz), c55(nx, nz), rho(nx, nz))
        allocate(K_x(nx), alpha_x(nx), a_x(nx), b_x(nx), &
                 K_x_half(nx), alpha_x_half(nx), a_x_half(nx), b_x_half(nx))
        allocate(K_z(nz), alpha_z(nz), a_z(nz), b_z(nz), &
                K_z_half(nz), alpha_z_half(nz), a_z_half(nz), b_z_half(nz))
        allocate(gamma_x(nx, nz), gamma_z(nx, nz), gamma_xz(nx, nz))
        allocate(srcx(source%time_steps), srcz(source%time_steps))
        
        ! Allocate more
        allocate(memory_dvx_dx(nx, nz), memory_dvx_dz(nx, nz))
        allocate(memory_dvz_dx(nx, nz), memory_dvz_dz(nx, nz))
        allocate(memory_dsigmaxx_dx(nx, nz), memory_dsigmazz_dz(nx, nz))
        allocate(memory_dsigmaxz_dx(nx, nz), memory_dsigmaxz_dz(nx, nz))
        allocate(vx(nx, nz), vz(nx, nz))
        allocate(pressure(nx,nz))
        allocate(sigmaxx(nx, nz), sigmazz(nx, nz), sigmaxz(nx, nz))
        
        ! -------------------- Load Stiffness Coefficients --------------------
    
        call material_rw2('c11.dat', c11, .TRUE.)
        call material_rw2('c13.dat', c13, .TRUE.)
        call material_rw2('c15.dat', c15, .TRUE.)
        call material_rw2('c33.dat', c33, .TRUE.)
        call material_rw2('c35.dat', c35, .TRUE.)
        call material_rw2('c55.dat', c55, .TRUE.)
        call material_rw2('rho.dat', rho, .TRUE.)
                
        ! ------------------- Load Attenuation Coefficients --------------------
        call material_rw2('gamma_x.dat', gamma_x, .TRUE.)
        call material_rw2('gamma_z.dat', gamma_z, .TRUE.)
        call material_rw2('gamma_xz.dat', gamma_xz, .TRUE.)
        
        ! ------------------------ Assign some constants -----------------------
        
        isource = source%xind + domain%cpml
        jsource = source%zind + domain%cpml
    
        ! ================================ LOAD SOURCE =========================
    
        call loadsource('seismicsourcex.dat', source%time_steps, srcx)
        call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        
        ! -----------------------------------------------------------------------------
        !--- define profile of absorption in PML region

        ! Initialize PML 
        K_x(:) = 1.d0
        K_x_half(:) = 1.d0
        alpha_x(:) = 0.d0
        alpha_x_half(:) = 0.d0
        a_x(:) = 0.d0
        a_x_half(:) = 0.d0
        b_x(:) = 0.d0 
        b_x_half(:) = 0.d0

        K_z(:) = 1.d0
        K_z_half(:) = 1.d0 
        alpha_z(:) = 0.d0
        alpha_z_half(:) = 0.d0
        a_z(:) = 0.d0
        a_z_half(:) = 0.d0
        b_z(:) = 0.d0
        b_z_half(:) = 0.d0

        ! ------------------------------ Load the boundary ----------------------------        
        call loadcpml('kappax_cpml.dat', K_x)
        call loadcpml('alphax_cpml.dat', alpha_x)
        call loadcpml('acoefx_cpml.dat', a_x)
        call loadcpml('bcoefx_cpml.dat', b_x)
        
        call loadcpml('kappaz_cpml.dat', K_z)
        call loadcpml('alphaz_cpml.dat', alpha_z)
        call loadcpml('acoefz_cpml.dat', a_z)
        call loadcpml('bcoefz_cpml.dat', b_z)
        
        call loadcpml('kappax_half_cpml.dat', K_x_half)
        call loadcpml('alphax_half_cpml.dat', alpha_x_half)
        call loadcpml('acoefx_half_cpml.dat', a_x_half)
        call loadcpml('bcoefx_half_cpml.dat', b_x_half)
        
        call loadcpml('kappaz_half_cpml.dat', K_z_half)
        call loadcpml('alphaz_half_cpml.dat', alpha_z_half)
        call loadcpml('acoefz_half_cpml.dat', a_z_half)
        call loadcpml('bcoefz_half_cpml.dat', b_z_half)

        ! =============================================================================
        vx(:,:) = 0.d0
        vz(:,:) = 0.d0
        
        ! Load initial condition
        call material_rw2('initialconditionVx.dat', vx, .TRUE.)
        call material_rw2('initialconditionVz.dat', vz, .TRUE.)
        
        
        sigmaxx(:,:) = 0.d0
        sigmazz(:,:) = 0.d0
        sigmaxz(:,:) = 0.d0

        ! PML
        memory_dvx_dx(:,:) = 0.d0
        memory_dvx_dz(:,:) = 0.d0
        memory_dvz_dx(:,:) = 0.d0
        memory_dvz_dz(:,:) = 0.d0
        memory_dsigmaxx_dx(:,:) = 0.d0
        memory_dsigmazz_dz(:,:) = 0.d0
        memory_dsigmaxz_dx(:,:) = 0.d0
        memory_dsigmaxz_dz(:,:) = 0.d0

        ! Load the geometry masked array 
        call material_rw2('geometry_mask.dat', material_geometry, .TRUE.)
        
        !---
        !---  beginning of time loop
        !---
        
        do it = 1,source%time_steps
            !$omp parallel private(i, j, deltarho, value_dvx_dx, value_dvx_dz, value_dvz_dx, value_dvz_dz, value_dsigmaxx_dx, value_dsigmazz_dz, value_dsigmaxz_dx, value_dsigmaxz_dz)
            ! ------------------------------------------------------------
            !  compute stress sigma and update memory variables for C-PML
            ! ------------------------------------------------------------
            
            !$omp do
            do j = 2,nz
                do i = 1,nx-1
                    
                    if ( material_geometry(i,j) == 0) then
                        value_dvx_dx = (vx(i+1,j) - vx(i,j)) / dx
                        value_dvz_dz = (vz(i,j) - vz(i,j-1)) / dz
                        value_dvz_dx = (vz(i+1,j) - vz(i,j)) / dx
                        value_dvx_dz = (vx(i,j) - vx(i,j-1)) / dz
                        
                        memory_dvx_dx(i,j) = b_x_half(j) * memory_dvx_dx(i,j) + &
                                                a_x_half(i) * value_dvx_dx
                        memory_dvz_dz(i,j) = b_z(j) * memory_dvz_dz(i,j) + &
                                                a_z(j) * value_dvz_dz
                        memory_dvx_dz(i,j) = b_z_half(j) * memory_dvx_dz(i,j) + &
                                                a_z_half(j) * value_dvx_dz 
                        memory_dvz_dx(i,j) = b_x(i) * memory_dvz_dx(i,j) + &
                                                a_x(i) * value_dvz_dx
                        
                        value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j)
                        value_dvz_dz = value_dvz_dz / K_z(j) + memory_dvz_dz(i,j)
                        value_dvz_dx = value_dvz_dx / K_x(i) + memory_dvz_dx(i,j)
                        value_dvx_dz = value_dvx_dz / K_z_half(j) + memory_dvx_dz(i,j)
                        
                        sigmaxx(i,j) = ( sigmaxx(i,j) + &
                            (   c11(i,j) * value_dvx_dx + &
                                c13(i,j) * value_dvz_dz + &
                                c15(i,j) * (value_dvz_dx + value_dvx_dz) ) * dt ) / & 
                                    (1 + gamma_x(i,j) * dt )
                        sigmazz(i,j) = ( sigmazz(i,j) + &
                            (   c13(i,j) * value_dvx_dx + &
                                c33(i,j) * value_dvz_dz + &
                                c35(i,j) * (value_dvz_dx + value_dvx_dz) ) * dt) / & 
                                    (1 + gamma_z(i,j) * dt )
                    else
                        dp_dx = (vx(i+1,j) - vx(i,j) ) / dx 
                        dp_dz = (vz(i,j) - vz(i,j-1)) / dz 
                        pressure(i,j) = pressure(i,j) - rho(i,j) * c11(i,j) * dt * (dp_dx + dp_dz)
                    
                    end if
                enddo
            enddo
            !$omp end do
            
            !$omp do
            do j = 1,nz-1
                do i = 2,nx
                    
                    if (material_geometry(i,j) == 0) then
                        value_dvx_dx = (vx(i,j) - vx(i-1,j)) / dx
                        value_dvz_dz = (vz(i,j+1) - vz(i,j)) / dz
                        value_dvz_dx = (vz(i,j) - vz(i-1,j)) / dx
                        value_dvx_dz = (vx(i,j+1) - vx(i,j)) / dz
                        
                        memory_dvx_dx(i,j) = b_x_half(i) * memory_dvx_dx(i,j) + &
                                                a_x_half(i) * value_dvx_dx
                        memory_dvz_dz(i,j) = b_z(j) * memory_dvz_dz(i,j) + &
                                                a_z(j) * value_dvz_dz
                        memory_dvx_dz(i,j) = b_z_half(j) * memory_dvx_dz(i,j) + &
                                                a_z_half(j) * value_dvx_dz 
                        memory_dvz_dx(i,j) = b_x(i) * memory_dvz_dx(i,j) + &
                                                a_x(i) * value_dvz_dx
                        
                        value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j)
                        value_dvz_dz = value_dvz_dz / K_z(j) + memory_dvz_dz(i,j)
                        value_dvz_dx = value_dvz_dx / K_x(i) + memory_dvz_dx(i,j)
                        value_dvx_dz = value_dvx_dz / K_z_half(j) + memory_dvx_dz(i,j)
                        
                        sigmaxz(i,j) = ( sigmaxz(i,j) + &
                            (   c15(i,j)  * value_dvx_dx + & 
                                c35(i,j)  * value_dvz_dz + &
                                c55(i,j) * (value_dvz_dx + value_dvx_dz) ) * dt ) / &
                                    (1 + gamma_xz(i,j) * dt )
                    end if
                enddo
            enddo
            !$omp end do
            
            ! Run pressure update
            !$omp do
            do j = 2,nz-1
                do i = 2,nx-1
                    if (material_geometry(i,j) == 1) then
                        ! Check neighbors
                        if (material_geometry(i+1,j) == 0) then
                            pressure(i,j) = -sigmaxx(i+1,j)
                            vx(i,j) = vx(i+1,j)
                        else if (material_geometry(i-1,j) == 0) then
                            pressure(i,j) = -sigmaxx(i-1,j)
                            vx(i,j) = vx(i-1,j)
                        end if
                        if (material_geometry(i,j+1) == 0) then
                            pressure(i,j) = -sigmazz(i,j+1)
                            vz(i,j) = vz(i,j+1)
                        else if (material_geometry(i,j-1) == 0) then
                            pressure(i,j) = -sigmazz(i,j-1)
                            vz(i,j) = vz(i,j-1)
                        end if
                    end if
                end do
            end do
            !$omp end do

            ! --------------------------------------------------------
            !  compute velocity and update memory variables for C-PML
            ! --------------------------------------------------------
            
            !$omp do 
            do j = 2,nz
                do i = 2,nx
                    
                    deltarho = ( 2*rho(i,j) + rho(i,j-1) + rho(i-1,j) )/4
                    
                    value_dsigmaxx_dx = (sigmaxx(i,j) - sigmaxx(i-1,j)) / dx
                    value_dsigmaxz_dz = (sigmaxz(i,j) - sigmaxz(i,j-1)) / dz
                    
                    memory_dsigmaxx_dx(i,j) = b_x(i) * memory_dsigmaxx_dx(i,j) + &
                                a_x(i) * value_dsigmaxx_dx
                    memory_dsigmaxz_dz(i,j) = b_z(j) * memory_dsigmaxz_dz(i,j) + &
                                a_z(j) * value_dsigmaxz_dz
                    
                    value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + &
                                memory_dsigmaxx_dx(i,j)
                    value_dsigmaxz_dz = value_dsigmaxz_dz / K_z(j) + &
                                memory_dsigmaxz_dz(i,j)
                    vx(i,j) = vx(i,j) + dt * (value_dsigmaxx_dx + value_dsigmaxz_dz) / deltarho
                    
                enddo
            enddo
            !$omp end do
            
            !$omp do
            do j = 1,nz-1
                do i = 1,nx-1
        
                    deltarho = ( 2*rho(i,j) + rho(i+1,j) + rho(i,j+1) )/4
                    value_dsigmaxz_dx = (sigmaxz(i+1,j) - sigmaxz(i,j)) / dx
                    value_dsigmazz_dz = (sigmazz(i,j+1) - sigmazz(i,j)) / dz
            
                    memory_dsigmaxz_dx(i,j) = b_x_half(i) * memory_dsigmaxz_dx(i,j) + &
                                a_x_half(i) * value_dsigmaxz_dx
                    memory_dsigmazz_dz(i,j) = b_z_half(j) * memory_dsigmazz_dz(i,j) + &
                                a_z_half(j) * value_dsigmazz_dz
            
                    value_dsigmaxz_dx = value_dsigmaxz_dx / K_x_half(i) + memory_dsigmaxz_dx(i,j)
                    value_dsigmazz_dz = value_dsigmazz_dz / K_z_half(j) + memory_dsigmazz_dz(i,j)

                    vz(i,j) = vz(i,j) + dt * (value_dsigmaxz_dx + value_dsigmazz_dz) / deltarho
                enddo
            enddo
            !$omp end do
            
            !$omp end parallel
            
            ! Add the source term
            ! vx(isource,jsource) = vx(isource,jsource) + srcx(it) * dt / rho(isource,jsource)
            ! vz(isource,jsource) = vz(isource,jsource) + srcz(it) * dt / rho(isource,jsource)
            vx(isource,jsource) = vx(isource,jsource) + srcx(it) / rho(isource,jsource)
            vz(isource,jsource) = vz(isource,jsource) + srcz(it) / rho(isource,jsource)
        
            ! Dirichlet conditions (rigid boundaries) on the edges or at the 
            ! bottom of the PML layers
            vx(1,:) = 0.d0
            vx(nx,:) = 0.d0
        
            vx(:,1) = 0.d0
            vx(:,nz) = 0.d0
        
            vz(1,:) = 0.d0
            vz(nx,:) = 0.d0
        
            vz(:,1) = 0.d0
            vz(:,nz) = 0.d0
        
            ! print maximum of norm of velocity
            velocnorm = maxval(sqrt(vx**2 + vz**2))
            if (velocnorm > stability_threshold) stop 'code became unstable and blew up'
        
            call write_image2(vx, nx, nz, source, it, 'Vx', SINGLE)
            call write_image2(vz, nx, nz, source, it, 'Vz', SINGLE)

        enddo   ! end of time loop
        
        deallocate(c11, c13, c15, c33, c35, c55, rho)
        deallocate(K_x, alpha_x, a_x, b_x, K_x_half, alpha_x_half, a_x_half, b_x_half)
        deallocate(K_z, alpha_z, a_z, b_z, K_z_half, alpha_z_half, a_z_half, b_z_half)
        deallocate(gamma_x, gamma_z, gamma_xz, srcx, srcz)
        
        ! Allocate more
        deallocate(memory_dvx_dx, memory_dvx_dz)
        deallocate(memory_dvz_dx, memory_dvz_dz)
        deallocate(memory_dsigmaxx_dx, memory_dsigmazz_dz)
        deallocate(memory_dsigmaxz_dx, memory_dsigmaxz_dz)
        deallocate(vx, vz, sigmaxx, sigmazz, sigmaxz)
        deallocate()
    end subroutine model2 
    
    ! =========================================================================
    subroutine model25
    
    end subroutine model25 

    
end module seismoacoustic