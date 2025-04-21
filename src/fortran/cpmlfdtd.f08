module cpmlfdtd 
    
    use seidartio
    use seidart_types
    
    implicit none 
    
    
    contains
    
    ! =========================================================================    
    ! Computations are done in double precision and written to binary as single
    ! precision unless specified by the optional logical, SINGLE_OUTPUT.
    subroutine seismic2(domain, source, SINGLE_OUTPUT)
        
        use constants
        use omp_lib ! Include the OpenMP library 
        
        ! Input arguments
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

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
        ! Load initial condition
        call material_rw2('initialconditionVx.dat', vx, .TRUE.)
        call material_rw2('initialconditionVz.dat', vz, .TRUE.)
        
        vx(:,:) = 0.d0
        vz(:,:) = 0.d0
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
                    
                    value_dvx_dx = (vx(i+1,j) - vx(i,j)) / dx
                    value_dvz_dz = (vz(i,j) - vz(i,j-1)) / dz
                    value_dvz_dx = (vz(i+1,j) - vz(i,j)) / dx
                    value_dvx_dz = (vx(i,j) - vx(i,j-1)) / dz
                    
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
                enddo
            enddo
            !$omp end do
            
            !$omp do
            do j = 1,nz-1
                do i = 2,nx
                    
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
                enddo
            enddo
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
        
    end subroutine seismic2

    ! =========================================================================
    subroutine seismic25(domain, source, SINGLE_OUTPUT)
        !--------------------------------------------------------------------------------------
        use constants
        
        implicit none

        ! Input arguments
        type(Domain_Type), intent(in) :: domain 
        type(Source_Type), intent(in) :: source 
        logical, intent(in), optional :: SINGLE_OUTPUT
        
        ! Local variables
        real(real64) :: deltarho, velocnorm

        
        real(real64), allocatable :: c11(:,:), c12(:,:), c13(:,:), c14(:,:), c15(:,:), c16(:,:), &
                                    c22(:,:), c23(:,:), c24(:,:), c25(:,:), c26(:,:), &
                                    c33(:,:), c34(:,:), c35(:,:), c36(:,:), &
                                    c44(:,:), c45(:,:), c46(:,:), &
                                    c55(:,:), c56(:,:), &
                                    c66(:,:), &
                                    rho(:,:)
        ! real(real64) :: DT
        integer :: i, j, k, it, isource, jsource, ksource

        ! Values of the velocity and stress differentials
        real(real64) :: dvx_dx, dvx_dy, dvx_dz, &
                        dvy_dx, dvy_dy, dvy_dz, &
                        dvz_dx, dvz_dy, dvz_dz, &
                        dsigmaxx_dx, dsigmayy_dy, dsigmazz_dz, &
                        dsigmaxy_dx, dsigmaxy_dy, &
                        dsigmaxz_dx, dsigmaxz_dz, &
                        dsigmayz_dy, dsigmayz_dz

        ! 1D arrays for the damping profiles in each direction
        real(real64), allocatable :: K_x(:), alpha_x(:), a_x(:), b_x(:), &
                        K_x_half(:), alpha_x_half(:), a_x_half(:), b_x_half(:)
        real(real64), allocatable :: K_y(:), alpha_y(:), a_y(:), b_y(:), &
                        K_y_half(:), alpha_y_half(:), a_y_half(:), b_y_half(:)
        real(real64), allocatable :: K_z(:), alpha_z(:), a_z(:), b_z(:), &
                        K_z_half(:), alpha_z_half(:), a_z_half(:), b_z_half(:)
        
        ! Arrays for the PML damping factors
        real(real64), allocatable :: gamma_x(:,:), gamma_y(:,:), gamma_z(:,:)
        real(real64), allocatable :: gamma_xz(:,:), gamma_xy(:,:), gamma_yz(:,:)

        ! Source arrays
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:)
        
        real(real64), allocatable :: vx(:,:,:), vy(:,:,:), vz(:,:,:), &
                        sigmaxx(:,:,:), sigmaxy(:,:,:), sigmaxz(:,:,:), &
                        sigmayy(:,:,:), sigmayz(:,:,:), sigmazz(:,:,:)
        real(real64), allocatable :: &
                memory_dvx_dx(:,:,:), memory_dvx_dy(:,:,:), memory_dvx_dz(:,:,:), &
                memory_dvy_dx(:,:,:), memory_dvy_dy(:,:,:), memory_dvy_dz(:,:,:), &
                memory_dvz_dx(:,:,:), memory_dvz_dy(:,:,:), memory_dvz_dz(:,:,:), &
                memory_dsigmaxx_dx(:,:,:), memory_dsigmayy_dy(:,:,:), memory_dsigmazz_dz(:,:,:), &
                memory_dsigmaxy_dx(:,:,:), memory_dsigmaxy_dy(:,:,:), memory_dsigmaxz_dx(:,:,:), &
                memory_dsigmaxz_dz(:,:,:), memory_dsigmayz_dy(:,:,:), memory_dsigmayz_dz(:,:,:)
        
        integer :: nx, ny, nz
        real(real64) :: dx, dy, dz, dt    
        
        ! Boolean flag to save as double precision or single precision
        logical :: SINGLE
    
        ! The default data output is single precision unless SINGLE_OUTPUT is 
        ! set to .FALSE.
        if (present(SINGLE_OUTPUT)) then 
            SINGLE = SINGLE_OUTPUT 
        else
            SINGLE = .TRUE.
        endif
        
         
        
        nx = domain%nx
        ny = domain%ny 
        nz = domain%nz 
        dt = source%dt 
        dx = domain%dx
        dy = domain%dy 
        dz = domain%dz
        
        ! ----------------------- Allocate Arrays ----------------------------
        allocate(c11(nx, nz), c12(nx, nz), c13(nx,nz), c14(nx,nz), c15(nx,nz), c16(nx,nz), &
                              c22(nx, nz), c23(nx,nz), c24(nx,nz), c25(nx,nz), c26(nx,nz), &
                                           c33(nx,nz), c34(nx,nz), c35(nx,nz), c36(nx,nz), &
                                                       c44(nx,nz), c45(nx,nz), c46(nx,nz), &
                                                                   c55(nx,nz), c56(nx,nz), &
                                                                               c66(nx,nz) )
        allocate(rho(nx,nz))
                
        allocate(K_x(nx), alpha_x(nx), a_x(nx), b_x(nx), &
                K_x_half(nx), alpha_x_half(nx), a_x_half(nx), b_x_half(nx))
        allocate(K_y(ny), alpha_y(ny), a_y(ny), b_y(ny), &
                K_y_half(ny), alpha_y_half(ny), a_y_half(ny), b_y_half(ny))
        allocate(K_z(nz), alpha_z(nz), a_z(nz), b_z(nz), &
                K_z_half(nz), alpha_z_half(nz), a_z_half(nz), b_z_half(nz))
        allocate(gamma_x(nx, nz), gamma_y(nx, nz), gamma_z(nx, nz))
        allocate(gamma_xy(nx, nz), gamma_yz(nx, nz), gamma_xz(nx, nz))
        
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps))
                
        ! Allocate more
        allocate(memory_dvx_dx(nx,ny,nz), memory_dvx_dy(nx,ny,nz), memory_dvx_dz(nx,ny,nz) )
        allocate(memory_dvy_dx(nx,ny,nz), memory_dvy_dy(nx,ny,nz), memory_dvy_dz(nx,ny,nz) )
        allocate(memory_dvz_dx(nx,ny,nz), memory_dvz_dy(nx,ny,nz), memory_dvz_dz(nx,ny,nz) )
        allocate(memory_dsigmaxx_dx(nx,ny,nz), memory_dsigmayy_dy(nx,ny,nz), memory_dsigmazz_dz(nx,ny,nz) )
        allocate(memory_dsigmaxy_dx(nx,ny,nz), memory_dsigmaxy_dy(nx,ny,nz), memory_dsigmaxz_dx(nx,ny,nz) )
        allocate(memory_dsigmaxz_dz(nx,ny,nz), memory_dsigmayz_dy(nx,ny,nz), memory_dsigmayz_dz(nx,ny,nz) )
        
        
        allocate(vx(nx, ny, nz), vy(nx, ny, nz), vz(nx, ny, nz))
        allocate(sigmaxx(nx, ny, nz), sigmaxy(nx, ny, nz), sigmaxz(nx, ny, nz))
        allocate(sigmayy(nx, ny, nz), sigmayz(nx, ny, nz), sigmazz(nx, ny, nz))
                
        
        ! ------------------------ Load Stiffness Coefficients ------------------------
            call material_rw2('c11.dat', c11, .TRUE.)
            call material_rw2('c12.dat', c12, .TRUE.)
            call material_rw2('c13.dat', c13, .TRUE.)
            call material_rw2('c14.dat', c14, .TRUE.)
            call material_rw2('c15.dat', c15, .TRUE.)
            call material_rw2('c16.dat', c16, .TRUE.)
            call material_rw2('c22.dat', c22, .TRUE.)
            call material_rw2('c23.dat', c23, .TRUE.)
            call material_rw2('c24.dat', c24, .TRUE.)
            call material_rw2('c25.dat', c25, .TRUE.)
            call material_rw2('c26.dat', c26, .TRUE.)
            call material_rw2('c33.dat', c33, .TRUE.)
            call material_rw2('c34.dat', c34, .TRUE.)
            call material_rw2('c35.dat', c35, .TRUE.)
            call material_rw2('c36.dat', c36, .TRUE.)
            call material_rw2('c44.dat', c44, .TRUE.)
            call material_rw2('c45.dat', c45, .TRUE.)
            call material_rw2('c46.dat', c46, .TRUE.)
            call material_rw2('c55.dat', c55, .TRUE.)
            call material_rw2('c56.dat', c56, .TRUE.)
            call material_rw2('c66.dat', c66, .TRUE.)
            call material_rw2('rho.dat', rho, .TRUE.)
        
        ! ------------------- Load Attenuation Coefficients --------------------
            call material_rw2('gamma_x.dat', gamma_x, .TRUE.)
            call material_rw2('gamma_z.dat', gamma_z, .TRUE.)
            call material_rw2('gamma_y.dat', gamma_y, .TRUE.)
            call material_rw2('gamma_xz.dat', gamma_xz, .TRUE.)
            call material_rw2('gamma_yz.dat', gamma_yz, .TRUE.)
            call material_rw2('gamma_xy.dat', gamma_xy, .TRUE.)
        ! ------------------------ Assign some constants -----------------------
        isource = source%xind + domain%cpml
        jsource = source%yind + domain%cpml
        ksource = source%zind + domain%cpml

        ! ================================ LOAD SOURCE ================================

        call loadsource('seismicsourcex.dat', source%time_steps, srcx)
        call loadsource('seismicsourcey.dat', source%time_steps, srcy)
        call loadsource('seismicsourcez.dat', source%time_steps, srcz)

        ! ==================================== PML ====================================
        ! Initialize PML 
        K_x(:) = 1.d0
        K_x_half(:) = 1.d0
        alpha_x(:) = 0.d0
        alpha_x_half(:) = 0.d0
        a_x(:) = 0.d0
        a_x_half(:) = 0.d0

        K_y(:) = 1.d0
        K_y_half(:) = 1.d0
        alpha_y(:) = 0.d0
        alpha_y_half(:) = 0.d0
        a_y(:) = 0.d0
        a_y_half(:) = 0.d0

        K_z(:) = 1.d0
        K_z_half(:) = 1.d0 
        alpha_z(:) = 0.d0
        alpha_z_half(:) = 0.d0
        a_z(:) = 0.d0
        a_z_half(:) = 0.d0

        ! ------------------------- Boundary Conditions -------------------------
        call loadcpml('kappax_cpml.dat', K_x)
        call loadcpml('alphax_cpml.dat', alpha_x)
        call loadcpml('acoefx_cpml.dat', a_x)
        call loadcpml('bcoefx_cpml.dat', b_x)

        call loadcpml('kappay_cpml.dat', K_y)
        call loadcpml('alphay_cpml.dat', alpha_y)
        call loadcpml('acoefy_cpml.dat', a_y)
        call loadcpml('bcoefy_cpml.dat', b_y)

        call loadcpml('kappaz_cpml.dat', K_z)
        call loadcpml('alphaz_cpml.dat', alpha_z)
        call loadcpml('acoefz_cpml.dat', a_z)
        call loadcpml('bcoefz_cpml.dat', b_z)

        call loadcpml('kappax_half_cpml.dat', K_x_half)
        call loadcpml('alphax_half_cpml.dat', alpha_x_half)
        call loadcpml('acoefx_half_cpml.dat', a_x_half)
        call loadcpml('bcoefx_half_cpml.dat', b_x_half)

        call loadcpml('kappay_half_cpml.dat', K_y_half)
        call loadcpml('alphay_half_cpml.dat', alpha_y_half)
        call loadcpml('acoefy_half_cpml.dat', a_y_half)
        call loadcpml('bcoefy_half_cpml.dat', b_y_half)

        call loadcpml('kappaz_half_cpml.dat', K_z_half)
        call loadcpml('alphaz_half_cpml.dat', alpha_z_half)
        call loadcpml('acoefz_half_cpml.dat', a_z_half)
        call loadcpml('bcoefz_half_cpml.dat', b_z_half)

        ! Load initial condition
        call material_rw3('initialconditionVx.dat', vx, .TRUE.)
        call material_rw3('initialconditionVy.dat', vy, .TRUE.)
        call material_rw3('initialconditionVz.dat', vz, .TRUE.)

        ! Initialize the stress values
        sigmaxx(:,:,:) = 0.d0
        sigmayy(:,:,:) = 0.d0
        sigmazz(:,:,:) = 0.d0
        sigmaxy(:,:,:) = 0.d0
        sigmaxz(:,:,:) = 0.d0
        sigmayz(:,:,:) = 0.d0

        ! PML
        memory_dvx_dx(:,:,:) = 0.d0
        memory_dvx_dy(:,:,:) = 0.d0
        memory_dvx_dz(:,:,:) = 0.d0

        memory_dvy_dx(:,:,:) = 0.d0
        memory_dvy_dy(:,:,:) = 0.d0
        memory_dvy_dz(:,:,:) = 0.d0

        memory_dvz_dx(:,:,:) = 0.d0
        memory_dvz_dy(:,:,:) = 0.d0 
        memory_dvz_dz(:,:,:) = 0.d0

        memory_dsigmaxx_dx(:,:,:) = 0.d0
        memory_dsigmayy_dy(:,:,:) = 0.d0
        memory_dsigmazz_dz(:,:,:) = 0.d0
        memory_dsigmaxy_dx(:,:,:) = 0.d0
        memory_dsigmaxy_dy(:,:,:) = 0.d0
        memory_dsigmaxz_dx(:,:,:) = 0.d0
        memory_dsigmaxz_dz(:,:,:) = 0.d0
        memory_dsigmayz_dy(:,:,:) = 0.d0
        memory_dsigmayz_dz(:,:,:) = 0.d0

        ! Do it 
        
        ! =============================== Forward Model ===============================
        do it = 1,source%time_steps
            !------------------------------------------------------------
            ! compute stress sigma and update memory variables for C-PML
            !------------------------------------------------------------
            ! Update in the x direction
            do k = 2,nz
                do j = 2,ny
                    do i = 1,nx-1

                        dvx_dx = (vx(i+1,j,k) - vx(i,j,k) ) / dx
                        dvy_dx = (vy(i+1,j,k) - vy(i,j,k) ) / dx
                        dvz_dx = (vz(i+1,j,k) - vz(i,j,k) ) / dx 
                        dvy_dy = (vy(i,j,k) - vy(i,j-1,k) ) / dy
                        dvx_dy = (vx(i,j,k) - vx(i,j-1,k) ) / dy
                        dvz_dy = (vz(i,j,k) - vz(i,j-1,k) ) / dy
                        dvz_dz = (vz(i,j,k) - vz(i,j,k-1) ) / dz
                        dvx_dz = (vx(i,j,k) - vx(i,j,k-1) ) / dz
                        dvy_dz = (vy(i,j,k) - vy(i,j,k-1) ) / dz

                        memory_dvx_dx(i,j,k) = b_x_half(i) * memory_dvx_dx(i,j,k) + a_x_half(i) * dvx_dx
                        memory_dvy_dx(i,j,k) = b_x_half(i) * memory_dvy_dx(i,j,k) + a_x_half(i) * dvy_dx
                        memory_dvz_dx(i,j,k) = b_x_half(i) * memory_dvz_dx(i,j,k) + a_x_half(i) * dvz_dx
                        memory_dvy_dy(i,j,k) = b_y(j) * memory_dvy_dy(i,j,k) + a_y(j) * dvy_dy
                        memory_dvx_dy(i,j,k) = b_y(j) * memory_dvx_dy(i,j,k) + a_y(j) * dvx_dy
                        memory_dvz_dy(i,j,k) = b_y(j) * memory_dvz_dy(i,j,k) + a_y(j) * dvz_dy
                        memory_dvz_dz(i,j,k) = b_z(k) * memory_dvz_dz(i,j,k) + a_z(k) * dvz_dz
                        memory_dvx_dz(i,j,k) = b_z(k) * memory_dvx_dz(i,j,k) + a_z(k) * dvx_dz
                        memory_dvy_dz(i,j,k) = b_z(k) * memory_dvy_dz(i,j,k) + a_z(k) * dvy_dz

                        dvx_dx = dvx_dx / K_x_half(i) + memory_dvx_dx(i,j,k)
                        dvy_dx = dvy_dx / K_x_half(i) + memory_dvy_dx(i,j,k)
                        dvz_dx = dvz_dx / K_x_half(i) + memory_dvz_dx(i,j,k)
                        dvy_dy = dvy_dy / K_y(j) + memory_dvy_dy(i,j,k)
                        dvx_dy = dvx_dy / K_y(j) + memory_dvx_dy(i,j,k)
                        dvz_dy = dvz_dy / K_y(j) + memory_dvz_dy(i,j,k)
                        dvz_dz = dvz_dz / K_z(k) + memory_dvz_dz(i,j,k)
                        dvx_dz = dvx_dz / K_z(k) + memory_dvx_dz(i,j,k)
                        dvy_dz = dvy_dz / K_z(k) + memory_dvy_dz(i,j,k)

                        sigmaxx(i,j,k) = ( sigmaxx(i,j,k) + &
                        (   c11(i,k) * dvx_dx + c12(i,k) * dvy_dy + c13(i,k) * dvz_dz + &
                            c14(i,k) * (dvy_dz + dvz_dy) + c15(i,k) * (dvx_dz + dvz_dx) + &
                            c16(i,k) * (dvx_dy + dvz_dy) ) * dt ) / &
                                (1 + gamma_x(i,j) * dt )

                        ! Full 3D will need a gradient in the y-direction
                        sigmayy(i,j,k) = ( sigmayy(i,j,k) + &
                        (   c12(i,k) * dvx_dx + c22(i,k) * dvy_dy + c23(i,k) * dvz_dz + &
                            c24(i,k) * (dvy_dz + dvz_dy) + c25(i,k) * (dvx_dz + dvz_dx) + &
                            c26(i,k) * (dvy_dx + dvx_dy) ) * dt ) / &
                                (1 + gamma_y(i,j) * dt )

                        sigmazz(i,j,k) = ( sigmazz(i,j,k) + &
                        (   c13(i,k) * dvx_dx + c23(i,k) * dvy_dy + c33(i,k) * dvz_dz + &
                            c34(i,k) * (dvy_dz + dvz_dy) + c35(i,k) * (dvx_dz + dvz_dx) + &
                            c36(i,k) * (dvy_dx + dvx_dy) ) * dt ) / &
                                (1 + gamma_z(i,j) * dt )

                    enddo
                enddo
            enddo

            ! Update sigmaxy, x-direction is full nodes
            do k = 2,nz
                do j = 1,ny-1
                    do i = 2,nx

                        dvx_dx = (vx(i,j,k) - vx(i-1,j,k)) / dx
                        dvy_dx = (vy(i,j,k) - vy(i-1,j,k)) / dx
                        dvz_dx = (vz(i,j,k) - vz(i-1,j,k)) / dx
                        dvy_dy = (vy(i,j+1,k) - vy(i,j,k)) / dy
                        dvx_dy = (vx(i,j+1,k) - vx(i,j,k)) / dy
                        dvz_dy = (vz(i,j+1,k) - vz(i,j,k)) / dy
                        dvz_dz = (vz(i,j,k) - vz(i,j,k-1)) / dz
                        dvx_dz = (vx(i,j,k) - vx(i,j,k-1)) / dz
                        dvy_dz = (vy(i,j,k) - vy(i,j,k-1)) / dz

                        memory_dvx_dx(i,j,k) = b_x(i) * memory_dvx_dx(i,j,k) + a_x(i) * dvx_dx
                        memory_dvy_dx(i,j,k) = b_x(i) * memory_dvy_dx(i,j,k) + a_x(i) * dvy_dx
                        memory_dvz_dx(i,j,k) = b_x(i) * memory_dvz_dx(i,j,k) + a_x(i) * dvz_dx
                        memory_dvy_dy(i,j,k) = b_y_half(j) * memory_dvy_dy(i,j,k) + a_y_half(j) * dvy_dy
                        memory_dvx_dy(i,j,k) = b_y_half(j) * memory_dvx_dy(i,j,k) + a_y_half(j) * dvx_dy
                        memory_dvz_dy(i,j,k) = b_y_half(j) * memory_dvz_dy(i,j,k) + a_y_half(j) * dvz_dy
                        memory_dvz_dz(i,j,k) = b_z_half(k) * memory_dvz_dz(i,j,k) + a_z_half(k) * dvz_dz
                        memory_dvx_dz(i,j,k) = b_z_half(k) * memory_dvx_dz(i,j,k) + a_z_half(k) * dvx_dz
                        memory_dvy_dz(i,j,k) = b_z_half(k) * memory_dvy_dz(i,j,k) + a_z_half(k) * dvy_dz

                        dvx_dx = dvx_dx / K_x(i) + memory_dvx_dx(i,j,k)
                        dvy_dx = dvy_dx / K_x(i) + memory_dvy_dx(i,j,k)
                        dvz_dx = dvz_dx / K_x(i) + memory_dvz_dx(i,j,k)
                        dvy_dy = dvy_dy / K_y_half(j) + memory_dvy_dy(i,j,k)
                        dvx_dy = dvx_dy / K_y_half(j) + memory_dvx_dy(i,j,k)
                        dvz_dy = dvz_dy / K_y_half(j) + memory_dvz_dy(i,j,k)
                        dvz_dz = dvz_dz / K_z_half(k) + memory_dvz_dz(i,j,k)
                        dvy_dz = dvy_dz / K_z_half(k) + memory_dvy_dz(i,j,k)

                        sigmaxy(i,j,k) = ( sigmaxy(i,j,k) + &
                        (   c16(i,k) * dvx_dx + c26(i,k) * dvy_dy + c36(i,k) * dvz_dz + &
                            c46(i,k) * (dvz_dy + dvy_dz) + c56(i,k) * (dvz_dx + dvx_dz) + &
                            c66(i,k) * (dvy_dx + dvx_dy) ) * dt ) / &
                                (1 + gamma_xy(i,j) * dt )

                    enddo
                enddo
            enddo

            ! Update sigmaxz, z-direction is full nodes
            do k = 1,nz-1
                do j = 2,ny
                    do i = 2,nx

                        dvx_dx = (vx(i,j,k) - vx(i-1,j,k)) / dx
                        dvy_dx = (vy(i,j,k) - vy(i-1,j,k)) / dx
                        dvz_dx = (vz(i,j,k) - vz(i-1,j,k)) / dx
                        dvy_dy = (vy(i,j,k) - vy(i,j-1,k)) / dy
                        dvz_dy = (vz(i,j,k) - vz(i,j-1,k)) / dy
                        dvx_dy = (vx(i,j,k) - vx(i,j-1,k)) / dy
                        dvz_dz = (vz(i,j,k+1) - vz(i,j,k)) / dz
                        dvx_dz = (vx(i,j,k+1) - vx(i,j,k)) / dz
                        dvy_dz = (vy(i,j,k+1) - vy(i,j,k)) / dz

                        memory_dvx_dx(i,j,k) = b_x(i) * memory_dvx_dx(i,j,k) + a_x(i) * dvx_dx
                        memory_dvy_dx(i,j,k) = b_x(i) * memory_dvy_dx(i,j,k) + a_x(i) * dvy_dx
                        memory_dvz_dx(i,j,k) = b_x(i) * memory_dvz_dx(i,j,k) + a_x(i) * dvz_dx
                        memory_dvy_dy(i,j,k) = b_y(j) * memory_dvy_dy(i,j,k) + a_y(j) * dvy_dy
                        memory_dvx_dy(i,j,k) = b_y(j) * memory_dvx_dy(i,j,k) + a_y(j) * dvx_dy
                        memory_dvz_dy(i,j,k) = b_y(j) * memory_dvz_dy(i,j,k) + a_y(j) * dvz_dy
                        memory_dvz_dz(i,j,k) = b_z_half(k) * memory_dvz_dz(i,j,k) + a_z_half(k) * dvz_dz
                        memory_dvx_dz(i,j,k) = b_z_half(k) * memory_dvx_dz(i,j,k) + a_z_half(k) * dvx_dz
                        memory_dvy_dz(i,j,k) = b_z_half(k) * memory_dvy_dz(i,j,k) + a_z_half(k) * dvy_dz

                        dvx_dx = dvx_dx / K_x(i) + memory_dvx_dx(i,j,k)
                        dvy_dx = dvy_dx / K_x(i) + memory_dvy_dx(i,j,k)
                        dvz_dx = dvz_dx / K_x(i) + memory_dvz_dx(i,j,k) 
                        dvy_dy = dvy_dy / K_y(j) + memory_dvy_dy(i,j,k)
                        dvx_dy = dvx_dy / K_y(j) + memory_dvx_dy(i,j,k)
                        dvz_dy = dvz_dy / K_y(j) + memory_dvz_dy(i,j,k)
                        dvz_dz = dvz_dz / K_z_half(k) + memory_dvz_dz(i,j,k)
                        dvx_dz = dvx_dz / K_z_half(k) + memory_dvx_dz(i,j,k)
                        dvy_dz = dvy_dz / K_z_half(k) + memory_dvy_dz(i,j,k)

                        sigmaxz(i,j,k) = ( sigmaxz(i,j,k) + &
                            (   c15(i,k) * dvx_dx + c25(i,k) * dvy_dy + c35(i,k) * dvz_dz + &
                                c45(i,k) * ( dvx_dz + dvz_dx) + c55(i,k) * ( dvx_dz + dvz_dx) + &
                                c56(i,k) * ( dvx_dy + dvy_dx) ) * dt  ) / &
                                    (1 + gamma_xz(i,j) * dt )

                    enddo
                enddo

                !   ! update sigmayz, y-direction is full nodes
                do j = 1,ny-1
                    do i = 1,nx-1

                        dvx_dx = (vx(i+1,j,k) - vx(i,j,k)) / dx
                        dvy_dx = (vy(i+1,j,k) - vy(i,j,k)) / dx
                        dvz_dx = (vz(i+1,j,k) - vz(i,j,k)) / dx
                        dvy_dy = (vy(i,j+1,k) - vy(i,j,k)) / dy
                        dvx_dy = (vx(i,j+1,k) - vx(i,j,k)) / dy
                        dvz_dy = (vz(i,j+1,k) - vz(i,j,k)) / dy
                        dvz_dz = (vz(i,j,k+1) - vz(i,j,k)) / dz
                        dvx_dz = (vx(i,j,k+1) - vx(i,j,k)) / dz 
                        dvy_dz = (vy(i,j,k+1) - vy(i,j,k)) / dz

                        memory_dvx_dx(i,j,k) = b_x_half(i) * memory_dvx_dx(i,j,k) + a_x_half(i) * dvx_dx
                        memory_dvy_dx(i,j,k) = b_x_half(i) * memory_dvy_dx(i,j,k) + a_x_half(i) * dvy_dx
                        memory_dvz_dx(i,j,k) = b_x_half(i) * memory_dvz_dx(i,j,k) + a_x_half(i) * dvz_dx
                        memory_dvy_dy(i,j,k) = b_y_half(j) * memory_dvy_dy(i,j,k) + a_y_half(j) * dvy_dy
                        memory_dvx_dy(i,j,k) = b_y_half(j) * memory_dvx_dy(i,j,k) + a_y_half(j) * dvx_dy
                        memory_dvz_dy(i,j,k) = b_y_half(j) * memory_dvz_dy(i,j,k) + a_y_half(j) * dvz_dy
                        memory_dvz_dz(i,j,k) = b_z_half(k) * memory_dvz_dz(i,j,k) + a_z_half(k) * dvz_dz
                        memory_dvx_dz(i,j,k) = b_z_half(k) * memory_dvx_dz(i,j,k) + a_z_half(k) * dvx_dz
                        memory_dvy_dz(i,j,k) = b_z_half(k) * memory_dvy_dz(i,j,k) + a_z_half(k) * dvy_dz

                        dvx_dx = dvx_dx / K_x_half(i) + memory_dvx_dx(i,j,k)
                        dvy_dx = dvy_dx / K_x_half(i) + memory_dvy_dx(i,j,k)
                        dvz_dx = dvz_dx / K_x_half(i) + memory_dvz_dx(i,j,k)
                        dvy_dy = dvy_dy / K_y_half(j) + memory_dvy_dy(i,j,k)
                        dvx_dy = dvx_dy / K_y_half(j) + memory_dvx_dy(i,j,k)
                        dvz_dy = dvz_dy / K_y_half(j) + memory_dvz_dy(i,j,k)
                        dvz_dz = dvz_dz / K_z_half(k) + memory_dvz_dz(i,j,k)
                        dvx_dz = dvx_dz / K_z_half(k) + memory_dvx_dz(i,j,k)
                        dvy_dz = dvy_dz / K_z_half(k) + memory_dvy_dz(i,j,k)

                        sigmayz(i,j,k) = ( sigmayz(i,j,k)  + &
                            (   c14(i,k) * dvx_dx + c24(i,k) * dvy_dy + c34(i,k) * dvz_dz + &
                                c44(i,k) * ( dvy_dz + dvz_dy) + c45(i,k) * ( dvx_dz + dvz_dx) + &
                                c46(i,k) * ( dvy_dx + dvx_dy) ) * dt  ) / & 
                                    (1 + gamma_yz(i,j) * dt )
                    enddo
                enddo
            enddo

            !--------------------------------------------------------
            ! compute velocity and update memory variables for C-PML
            !--------------------------------------------------------
            do k = 2,nz
                do j = 2,ny
                    do i = 2,nx
                        ! ds1/dx, ds6/dy, ds5,dz
                        deltarho = (4 * rho(i,k) + rho(i-1,k) + rho(i,k-1) )/6

                        dsigmaxx_dx = (sigmaxx(i,j,k) - sigmaxx(i-1,j,k) ) / dx
                        dsigmaxy_dy = (sigmaxy(i,j,k) - sigmaxy(i,j-1,k) ) / dy
                        dsigmaxz_dz = (sigmaxz(i,j,k) - sigmaxz(i,j,k-1) ) / dz

                        memory_dsigmaxx_dx(i,j,k) = b_x(i) * &
                            memory_dsigmaxx_dx(i,j,k) + a_x(i) * dsigmaxx_dx
                        memory_dsigmaxy_dy(i,j,k) = b_y(j) * &
                            memory_dsigmaxy_dy(i,j,k) + a_y(j) * dsigmaxy_dy
                        memory_dsigmaxz_dz(i,j,k) = b_z(k) * &
                            memory_dsigmaxz_dz(i,j,k) + a_z(k) * dsigmaxz_dz

                        dsigmaxx_dx = dsigmaxx_dx / K_x(i) + memory_dsigmaxx_dx(i,j,k)
                        dsigmaxy_dy = dsigmaxy_dy / K_y(j) + memory_dsigmaxy_dy(i,j,k)
                        dsigmaxz_dz = dsigmaxz_dz / K_z(k) + memory_dsigmaxz_dz(i,j,k) 

                        vx(i,j,k) = vx(i,j,k) + &
                            (dsigmaxx_dx + dsigmaxy_dy + dsigmaxz_dz) * &
                            dt / deltarho !rho(i,k)
                    enddo
                enddo

                do j = 1,ny-1
                    do i = 1,nx-1
                        ! ds6/dx, ds2/dy, ds4/dz
                        deltarho = (4*rho(i,k) + rho(i+1,k) + rho(i,k-1) )/6

                        dsigmaxy_dx = ( sigmaxy(i+1,j,k) - sigmaxy(i,j,k) ) / dx
                        dsigmayy_dy = ( sigmayy(i,j+1,k) - sigmayy(i,j,k) ) / dy
                        dsigmayz_dz = ( sigmayz(i,j,k) - sigmayz(i,j,k-1) ) / dz

                        memory_dsigmaxy_dx(i,j,k) = b_x_half(i) * memory_dsigmaxy_dx(i,j,k) + a_x_half(i) * dsigmaxy_dx
                        memory_dsigmayy_dy(i,j,k) = b_y_half(j) * memory_dsigmayy_dy(i,j,k) + a_y_half(j) * dsigmayy_dy
                        memory_dsigmayz_dz(i,j,k) = b_z(k) * memory_dsigmayz_dz(i,j,k) + a_z(k) * dsigmayz_dz

                        dsigmaxy_dx = dsigmaxy_dx / K_x_half(i) + memory_dsigmaxy_dx(i,j,k)
                        dsigmayy_dy = dsigmayy_dy / K_y_half(j) + memory_dsigmayy_dy(i,j,k)
                        dsigmayz_dz = dsigmayz_dz / K_z(k) + memory_dsigmayz_dz(i,j,k)

                        vy(i,j,k) = vy(i,j,k) + &
                            (dsigmaxy_dx + dsigmayy_dy + dsigmayz_dz) * &
                            dt / deltarho !rho(i,k)
                    enddo
                enddo
            enddo

            do k = 1,nz-1
                do j = 2,ny
                    do i = 1,nx-1
                        ! ds5/dx, ds4/dy, ds3/dz
                        deltarho = ( rho(i+1,k) + rho(i,k+1) + 4*rho(i,k) )/6

                        dsigmaxz_dx = ( sigmaxz(i+1,j,k) - sigmaxz(i,j,k) ) / dx
                        dsigmayz_dy = ( sigmayz(i,j,k) - sigmayz(i,j-1,k) ) / dy
                        dsigmazz_dz = ( sigmazz(i,j,k+1) - sigmazz(i,j,k) ) / dz

                        memory_dsigmaxz_dx(i,j,k) = b_x_half(i) * memory_dsigmaxz_dx(i,j,k) + a_x_half(i) * dsigmaxz_dx
                        memory_dsigmayz_dy(i,j,k) = b_y(j) * memory_dsigmayz_dy(i,j,k) + a_y(j) * dsigmayz_dy
                        memory_dsigmazz_dz(i,j,k) = b_z_half(k) * memory_dsigmazz_dz(i,j,k) + a_z_half(k) * dsigmazz_dz

                        dsigmaxz_dx = dsigmaxz_dx / K_x_half(i) + memory_dsigmaxz_dx(i,j,k)
                        dsigmayz_dy = dsigmayz_dy / K_y(j) + memory_dsigmayz_dy(i,j,k)
                        dsigmazz_dz = dsigmazz_dz / K_z_half(k) + memory_dsigmazz_dz(i,j,k)

                        vz(i,j,k) = vz(i,j,k) + &
                            (dsigmaxz_dx + dsigmayz_dy + dsigmazz_dz) * &
                            dt / deltarho !rho(i,k)

                    enddo
                enddo
            enddo

            vx(isource,jsource,ksource) = vx(isource,jsource,ksource) + &
                    srcx(it) * dt / rho(isource,ksource)
            vy(isource,jsource,ksource) = vy(isource,jsource,ksource) + &
                    srcy(it) * dt / rho(isource,ksource)
            vz(isource,jsource,ksource) = vz(isource,jsource,ksource) + &
                    srcz(it) * dt / rho(isource,ksource)

            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            vx(1,:,:) = 0.d0
            vy(1,:,:) = 0.d0
            vz(1,:,:) = 0.d0

            vx(:,1,:) = 0.d0
            vy(:,1,:) = 0.d0
            vz(:,1,:) = 0.d0

            vx(:,:,1) = 0.d0
            vy(:,:,1) = 0.d0
            vz(:,:,1) = 0.d0

            vx(nx,:,:) = 0.d0
            vy(nx,:,:) = 0.d0
            vz(nx,:,:) = 0.d0

            vx(:,ny,:) = 0.d0
            vy(:,ny,:) = 0.d0
            vz(:,ny,:) = 0.d0

            vx(:,:,nz) = 0.d0
            vy(:,:,nz) = 0.d0
            vz(:,:,nz) = 0.d0

            ! check norm of velocity to make sure the solution isn't diverging
            velocnorm = maxval( sqrt(vx**2 + vy**2 + vz**2) )
            ! print *,'Time step # ',it,' out of ',time_step
            ! print *,'Time: ',(it-1)*DT,' seconds'
            ! print *,'Max vals for vx, vy, vz: ', maxval(vx), maxval(vy), maxval(vz)

            if (velocnorm > stability_threshold) stop 'code became unstable and blew up'

            ! Write the velocity values to an unformatted binary file
            call write_image3(vx, nx, ny, nz, source, it, 'Vx', SINGLE)
            call write_image3(vy, nx, ny, nz, source, it, 'Vy', SINGLE)
            call write_image3(vz, nx, ny, nz, source, it, 'Vz', SINGLE)
            ! Now write the stress Values
            ! call write_image3(sigmaxx, nx, ny, nz, it, 'S1')
            ! call write_image3(sigmayy, nx, ny, nz, it, 'S2')
            ! call write_image3(sigmazz, nx, ny, nz, it, 'S3')
            ! call write_image3(sigmaxy, nx, ny, nz, it, 'S6')
            ! call write_image3(sigmayz, nx, ny, nz, it, 'S4')
            ! call write_image3(sigmaxz, nx, ny, nz, it, 'S5')

        enddo   ! end of time loop
    end subroutine seismic25

    ! =========================================================================
    ! FORWARD AND BACKWARD DIFFERENCE SCHEME
    subroutine electromag2(domain, source, SINGLE_OUTPUT)
        !--------------------------------------------------------------------------------------
        ! Electromagnetic wave propagation in a 2D grid with Convolutional-PML (C-PML)
        ! absorbing conditions for an anisotropic medium.
        !
        ! This subroutine solves electromagnetic wave propagation using a finite-difference
        ! time-domain (FDTD) method in a 2D grid, with PML absorbing conditions.
        !--------------------------------------------------------------------------------------
        
        use constants
        
        implicit none

        ! Input arguments
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT
         
        ! Local variabless
        real(real64), allocatable :: epsilonx(:,:), epsilonz(:,:), &
                                            sigmax(:,:), sigmaz(:,:)

        ! real(real64) :: DT
        real(real64), allocatable :: srcx(:), srcz(:) ! The vector time series of the source
        integer :: isource, jsource, i, j, it

        ! Coefficients for the finite difference scheme
        real(real64), allocatable :: caEx(:,:), cbEx(:,:)
        real(real64), allocatable :: caEz(:,:), cbEz(:,:)
        real(real64) :: daHy, dbHy
        real(real64) :: value_dEx_dz, value_dEz_dx, value_dHy_dz, value_dHy_dx

        ! 1D arrays for the damping profiles
        real(real64), allocatable :: K_x(:), alpha_x(:), a_x(:), b_x(:), K_x_half(:), alpha_x_half(:), a_x_half(:), b_x_half(:)
        real(real64), allocatable :: K_z(:), alpha_z(:), a_z(:), b_z(:), K_z_half(:), alpha_z_half(:), a_z_half(:), b_z_half(:)
        
        real(real64), allocatable :: Ex(:,:), Ez(:,:), Hy(:,:) 
        real(real64), allocatable :: memory_dEx_dz(:,:), memory_dEz_dx(:,:), &
                                        memory_dHy_dx(:,:), memory_dHy_dz(:,:)
        
        real(real64), allocatable :: eps11(:,:), eps13(:,:), eps33(:,:), &
                                            sig11(:,:), sig13(:,:), sig33(:,:)
        ! Velocity normalization factor
        real(real64) :: velocnorm
        
        integer :: nx, nz
        real(real64) :: dx, dz, dt   
        
        ! Boolean flag to save as double precision or single precision
        logical :: SINGLE

        ! Check if SINGLE_OUTPUT is provided, default to single precision if not
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        endif
        
        ! ----------------------------------------------------------------------
        nx = domain%nx 
        nz = domain%nz 
        dx = domain%dx 
        dz = domain%dz 
        dt = source%dt
        
        allocate(eps11(nx, nz), eps13(nx, nz),  &
                    eps33(nx, nz))
        allocate(sig11(nx, nz), sig13(nx, nz),  &
                    sig33(nx, nz))
        allocate(K_x(nx), alpha_x(nx), a_x(nx), b_x(nx), &
                K_x_half(nx), alpha_x_half(nx), a_x_half(nx), b_x_half(nx))
        allocate(K_z(nz), alpha_z(nz), a_z(nz), b_z(nz), &
                K_z_half(nz), alpha_z_half(nz), a_z_half(nz), b_z_half(nz))
        allocate(srcx(source%time_steps), srcz(source%time_steps))
        
        ! Allocate more
        allocate(epsilonx(nx, nz), epsilonz(nx, nz))
        allocate(sigmax(nx, nz), sigmaz(nx, nz))
        allocate(caEx(nx, nz), cbEx(nx, nz), caEz(nx, nz), cbEz(nx, nz))
        allocate(memory_dEz_dx(nx, nz), memory_dEx_dz(nx, nz))
        allocate(memory_dHy_dx(nx, nz), memory_dHy_dz(nx, nz))
        allocate(Ex(nx, nz), Ez(nx, nz), Hy(nx, nz))
            
        ! ======================================================================
        ! ----------------------- Load Permittivity Coefficients ----------------------
        call material_rw2('e11.dat', eps11, .TRUE.)
        call material_rw2('e13.dat', eps13, .TRUE.)
        call material_rw2('e33.dat', eps33, .TRUE.) ! We will change y to z soon
        call material_rw2('s11.dat', sig11, .TRUE.)
        call material_rw2('s13.dat', sig13, .TRUE.)
        call material_rw2('s33.dat', sig33, .TRUE.)

        ! ------------------------ Assign some constants -----------------------
        ! Assign the source location indices
        isource = source%xind + domain%cpml
        jsource = source%zind + domain%cpml

        ! ================================ LOAD SOURCE ================================
        call loadsource('electromagneticsourcex.dat', source%time_steps, srcx)
        call loadsource('electromagneticsourcez.dat', source%time_steps, srcz)

        ! ----------------------------------------------------------------------
        ! Initialize CPML damping variables
        K_x(:) = 1.0d0
        K_x_half(:) = 1.0d0
        alpha_x(:) = 0.0d0
        alpha_x_half(:) = 0.0d0
        a_x(:) = 0.0d0
        a_x_half(:) = 0.0d0
        b_x(:) = 0.0d0 
        b_x_half(:) = 0.0d0 

        K_z(:) = 1.0d0
        K_z_half(:) = 1.0d0
        alpha_z(:) = 0.0d0
        alpha_z_half(:) = 0.0d0
        a_z(:) = 0.0d0
        a_z_half(:) = 0.0d0

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

        ! ----------------------------------------------------------------------
        ! Load initial conditions and initialize variables
        call material_rw2('initialconditionEx.dat', Ex, .TRUE.)
        ! call material_rw2('initialconditionHy.dat', Hy, .TRUE.)
        call material_rw2('initialconditionEz.dat', Ez, .TRUE.)
    
        Hy(:,:) = 0.d0
        ! PML
        memory_dEx_dz(:,:) = 0.0d0
        memory_dEz_dx(:,:) = 0.0d0
        memory_dHy_dx(:,:) = 0.0d0
        memory_dHy_dz(:,:) = 0.0d0
        
        ! ----------------------------------------------------------------------
        ! Compute the coefficients of the FD scheme. First scale the relative 
        ! permittivity and permeabilities to get the absolute values 
        epsilonx(:,:) = (eps11 + eps13)*eps0
        epsilonz(:,:) = (eps13 + eps33)*eps0
        sigmax(:,:) = sig11 + sig13 
        sigmaz(:,:) = sig13 + sig33 

        ! We need to change sigma to dsigma, same for epsilon
        caEx(:,:) = ( 1.0d0 - sigmax * dt / ( 2.0d0 * epsilonx ) ) / &
                    ( 1.0d0 + sigmax * dt / (2.0d0 * epsilonx ) )
        cbEx(:,:) = (dt / epsilonx ) / ( 1.0d0 + sigmax * dt / ( 2.0d0 * epsilonx ) )

        caEz(:,:) = ( 1.0d0 - sigmaz * dt / ( 2.0d0 * epsilonz ) ) / &
                    ( 1.0d0 + sigmaz * dt / (2.0d0 * epsilonz ) )
        cbEz(:,:) = (dt / epsilonz ) / ( 1.0d0 + sigmaz * dt / ( 2.0d0 * epsilonz ) )
        daHy = dt/(4.0d0*mu0*mu)
        dbHy = dt/mu0 !dt/(mu*mu*dx*(1+daHy) ) 
        daHy = 1.0d0 ! (1-daHy)/(1+daHy) ! 
        
        !---
        !---  beginning of time loop
        !---
        do it = 1,source%time_steps
            !$omp parallel private(i, j, value_dEx_dz, value_dEz_dx, value_dHy_dz, value_dHy_dx) &
            !$omp& shared(Ex, Ez, Hy, memory_dEx_dz, memory_dEz_dx, memory_dHy_dz, memory_dHy_dx, b_z, b_x, a_z, a_x, K_z, K_x, caEz, cbEz, caEx, cbEx, daHy, dbHy) & 
            !$omp& reduction(max:velocnorm)
            
            !--------------------------------------------------------
            ! compute magnetic field and update memory variables for C-PML
            !--------------------------------------------------------
            !$omp do
            do j = 1,nz-1  
                do i = 1,nx-1
                
                    ! Values needed for the magnetic field updates
                    value_dEx_dz = ( Ex(i,j+1) - Ex(i,j) )/dz
                    memory_dEx_dz(i,j) = b_z(j) * memory_dEx_dz(i,j) + a_z(j) * value_dEx_dz
                    value_dEx_dz = value_dEx_dz/ K_z(j) + memory_dEx_dz(i,j)

                    ! The rest of the equation needed for agnetic field updates
                    value_dEz_dx = ( Ez(i+1,j) - Ez(i,j) )/dx
                    memory_dEz_dx(i,j) = b_x(i) * memory_dEz_dx(i,j) + a_x(i) * value_dEz_dx
                    value_dEz_dx = value_dEz_dx/ K_x(i) + memory_dEz_dx(i,j)

                    ! Now update the Magnetic field
                    Hy(i,j) = daHy*Hy(i,j) + dbHy*( value_dEz_dx + value_dEx_dz )
                    velocnorm = max(velocnorm, sqrt(Ex(i, j)**2 + Ez(i, j)**2))

                enddo  
            enddo
            !$omp end do 
            
            ! Electric field and update memory variables for C-PML
            ! Compute the differences in the y-direction
            !$omp do
            do j = 2,nz
                do i = 1,nx
                    ! Update the Ex field
                    value_dHy_dz = ( Hy(i,j) - Hy(i,j-1) )/dz ! this is nz-1 length vector
                    memory_dHy_dz(i,j) = b_z(j) * memory_dHy_dz(i,j) + a_z(j) * value_dHy_dz
                    value_dHy_dz = value_dHy_dz/K_z(j) + memory_dHy_dz(i,j)
                    
                    Ex(i,j) = caEx(i,j) * Ex(i,j) + cbEx(i,j) * value_dHy_dz
                enddo
            enddo
            !$omp end do 
            
            !$omp do
            do j = 1,nz
                do i = 2,nx
                    ! Update the Ez field
                    value_dHy_dx = ( Hy(i,j) - Hy(i-1,j) )/dx
                    memory_dHy_dx(i,j) = b_x_half(i) * memory_dHy_dx(i,j) + a_x_half(i) * value_dHy_dx
                    value_dHy_dx = value_dHy_dx/K_x_half(i) + memory_dHy_dx(i,j)
                    
                    Ez(i,j) = caEz(i,j) * Ez(i,j) + cbEz(i,j) * value_dHy_dx 
                enddo
            enddo
            !$omp end do 
            
            !$omp end parallel
            !----------------------------------------------------------------------------
            Ex(isource,jsource) = Ex(isource,jsource) + &
                            srcx(it) * dt / eps11(isource,jsource)
            Ez(isource,jsource) = Ez(isource,jsource) + &
                            srcz(it) * dt / eps33(isource,jsource) 
            
            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            Ex(1,:) = 0.d0
            Ex(nx,:) = 0.d0
            Ex(:,1) = 0.d0
            Ex(:,nz) = 0.d0

            Ez(1,:) = 0.d0
            Ez(nx,:) = 0.d0
            Ez(:,1) = 0.d0
            Ez(:,nz) = 0.d0

            Hy(1,:) = 0.d0
            Hy(nx,:) = 0.d0
            Hy(:,1) = 0.d0
            Hy(:,nz) = 0.d0

            ! print maximum of norm of velocity
            ! velocnorm = maxval(sqrt(Ex**2 + Ez**2))
            if (velocnorm > stability_threshold) stop 'code became unstable and blew up'

            call write_image2(Ex, nx, nz, source, it, 'Ex', SINGLE)
            call write_image2(Ez, nx, nz, source, it, 'Ez', SINGLE)
        enddo
        
        
        deallocate(eps11, eps13,  eps33, sig11, sig13,  sig33, srcx, srcz)
        deallocate(K_x, alpha_x, a_x, b_x, K_x_half, alpha_x_half, a_x_half, b_x_half)
        deallocate(K_z, alpha_z, a_z, b_z, K_z_half, alpha_z_half, a_z_half, b_z_half)
        deallocate(epsilonx, epsilonz, sigmax, sigmaz, Ex, Ez, Hy)
        deallocate(memory_dEz_dx, memory_dEx_dz, memory_dHy_dx, memory_dHy_dz)
        
    end subroutine electromag2


    ! =========================================================================
    subroutine electromag25(domain, source, SINGLE_OUTPUT)
        !--------------------------------------------------------------------------------------
        ! Electromagnetic wave propagation in a 3D grid with Convolutional-PML (C-PML)
        ! absorbing conditions for an anisotropic medium.
        !
        ! This subroutine solves electromagnetic wave propagation using a finite-difference
        ! time-domain (FDTD) method in a 3D grid, with PML absorbing conditions.
        !
        !--------------------------------------------------------------------------------------
        
        use constants
        
        implicit none

        ! Input arguments
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT
        
        ! Local variables
        real(real64), allocatable :: epsilonx(:,:), epsilony(:,:), epsilonz(:,:), &
                                        sigmax(:,:), sigmay(:,:), sigmaz(:,:)

        ! real(real64) :: DT
        real(real64) :: velocnorm
        integer :: isource, jsource, ksource, i, j, k, it

        ! Coefficients for the finite difference scheme
        real(real64), allocatable :: caEx(:,:), cbEx(:,:), caEy(:,:), cbEy(:,:), caEz(:,:), cbEz(:,:)
        real(real64) :: daHx, dbHx, daHy, dbHy, daHz, dbHz

        real(real64) :: dEx_dy, dEy_dx, dEy_dz, dEz_dy, dEz_dx, dEx_dz, &
                        dHx_dy, dHx_dz, dHy_dx, dHy_dz, dHz_dy, dHz_dx

        ! Source arrays
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:)

        ! 1D arrays for the damping profiles in each direction
        real(real64), allocatable :: K_x(:), alpha_x(:), a_x(:), b_x(:), & 
                        K_x_half(:), alpha_x_half(:), a_x_half(:), b_x_half(:)
        real(real64), allocatable :: K_y(:), alpha_y(:), a_y(:), b_y(:), & 
                        K_y_half(:), alpha_y_half(:), a_y_half(:), b_y_half(:)
        real(real64), allocatable :: K_z(:), alpha_z(:), a_z(:), b_z(:), & 
                        K_z_half(:), alpha_z_half(:), a_z_half(:), b_z_half(:)
        
        real(real64), allocatable :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:), &
                                    Hx(:,:,:), Hy(:,:,:), Hz(:,:,:) 
        real(real64), allocatable :: memory_dEx_dy(:,:,:), memory_dEy_dx(:,:,:), &
                                    memory_dEx_dz(:,:,:), memory_dEz_dx(:,:,:), &
                                    memory_dEz_dy(:,:,:), memory_dEy_dz(:,:,:)
        real(real64), allocatable :: memory_dHz_dx(:,:,:), memory_dHx_dz(:,:,:), &
                                    memory_dHz_dy(:,:,:), memory_dHy_dz(:,:,:), &
                                    memory_dHx_dy(:,:,:), memory_dHy_dx(:,:,:)
        real(real64), allocatable :: eps11(:,:), eps12(:,:), eps13(:,:), &
                                    eps22(:,:), eps23(:,:), eps33(:,:)
        real(real64), allocatable :: sig11(:,:), sig12(:,:), sig13(:,:), &
                                    sig22(:,:), sig23(:,:), sig33(:,:)
        
        
        integer :: nx, ny, nz
        real(real64) :: dx, dy, dz, dt   
            
        ! Boolean flag to save as double precision or single precision
        logical :: SINGLE

        ! Check if SINGLE_OUTPUT is provided, default to single precision if not
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        endif
        
        
        nx = domain%nx
        ny = domain%ny
        nz = domain%nz
        dx = domain%dx
        dy = domain%dy
        dz = domain%dz
        dt = source%dt
        allocate(eps11(nx,nz), eps12(nx,nz), eps13(nx,nz), &
                    eps22(nx,nz), eps23(nx,nz), eps33(nx,nz))
        allocate(sig11(nx,nz), sig12(nx,nz), sig13(nx,nz), &
                    sig22(nx,nz), sig23(nx,nz), sig33(nx,nz))
        allocate(K_x(nx), alpha_x(nx), a_x(nx), b_x(nx), K_x_half(nx), &
                    alpha_x_half(nx), a_x_half(nx), b_x_half(nx))
        allocate(K_y(ny), alpha_y(ny), a_y(ny), b_y(ny), K_y_half(ny), &
                    alpha_y_half(ny), a_y_half(ny), b_y_half(ny))
        allocate(K_z(nz), alpha_z(nz), a_z(nz), b_z(nz), K_z_half(nz), &
                    alpha_z_half(nz), a_z_half(nz), b_z_half(nz))
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps))
        
        ! Allocate more
        allocate(epsilonx(nx,nz), epsilony(nx,nz), epsilonz(nx,nz))
        allocate(sigmax(nx,nz), sigmay(nx,nz), sigmaz(nx, nz))
        allocate(caEx(nx,nz), cbEx(nx,nz), caEy(nx,nz), cbEy(nx,nz), &
                    caEz(nx,nz), cbEz(nx,nz))
        allocate( memory_dEx_dy(nx,ny,nz), memory_dEy_dx(nx,ny,nz), &
                    memory_dEx_dz(nx,ny,nz), memory_dEz_dx(nx,ny,nz), &
                    memory_dEz_dy(nx,ny,nz), memory_dEy_dz(nx,ny,nz) )
        allocate(memory_dHz_dx(nx,ny,nz), memory_dHx_dz(nx,ny,nz), & 
                    memory_dHz_dy(nx,ny,nz), memory_dHy_dz(nx,ny,nz), &
                    memory_dHx_dy(nx,ny,nz), memory_dHy_dx(nx,ny,nz) )
        
        allocate(Ex(nx, ny, nz), Ey(nx, ny, nz), Ez(nx, ny, nz))
        allocate(Hx(nx, ny, nz), Hy(nx, ny, nz), Hz(nx, ny, nz))
        
        ! ------------------------ Load Permittivity Coefficients ------------------------
        ! Load Epsilon
        call material_rw2('e11.dat', eps11, .TRUE.)
        call material_rw2('e12.dat', eps12, .TRUE.)
        call material_rw2('e13.dat', eps13, .TRUE.)
        call material_rw2('e22.dat', eps22, .TRUE.)
        call material_rw2('e23.dat', eps23, .TRUE.)
        call material_rw2('e33.dat', eps33, .TRUE.)
        ! Load Sigma
        call material_rw2('s11.dat', sig11, .TRUE.)
        call material_rw2('s12.dat', sig12, .TRUE.)
        call material_rw2('s13.dat', sig13, .TRUE.)
        call material_rw2('s22.dat', sig22, .TRUE.)
        call material_rw2('s23.dat', sig23, .TRUE.)
        call material_rw2('s33.dat', sig33, .TRUE.)

        ! ------------------------ Assign some constants -----------------------
        ! Assign the source location indices
        isource = source%xind + domain%cpml
        jsource = source%yind + domain%cpml
        ksource = source%zind + domain%cpml

        ! Define the 
        ! DT = minval( (/dx, dy, dz/) )/ ( 2.0d0 * Clight/ sqrt( minval( (/ eps11, eps22, eps33 /) ) ) )

        ! Compute the coefficients of the FD scheme. First scale the relative 
        ! permittivity and permeabilities to get the absolute values 
        epsilonx(:,:) = (eps11 + eps12 + eps13)*eps0 
        epsilony(:,:) = (eps12 + eps22 + eps23)*eps0
        epsilonz(:,:) = (eps13 + eps23 + eps33)*eps0
        sigmax(:,:) = sig11 + sig12 + sig13
        sigmay(:,:) = sig12 + sig22 + sig23
        sigmaz(:,:) = sig13 + sig23 + sig33

        caEx(:,:) = ( 1.0d0 - sigmax * dt / (2.0d0 * epsilonx ) ) / &
                    ( 1.0d0 + sigmax * dt / (2.0d0 * epsilonx ) )
        cbEx(:,:) = (dt / epsilonx ) / ( 1.0d0 + sigmax * dt / ( 2.0d0 * epsilonx ) )

        caEy(:,:) = ( 1.0d0 - sigmay * dt / (2.0d0 * epsilony ) ) / &
                    ( 1.0d0 + sigmay * dt / (2.0d0 * epsilony ) )
        cbEy(:,:) = (dt / epsilony ) / ( 1.0d0 + sigmay * dt / ( 2.0d0 * epsilony ) )

        caEz(:,:) = ( 1.0d0 - sigmaz * dt / ( 2.0d0 * epsilonz ) ) / &
                    ( 1.0d0 + sigmaz * dt / (2.0d0 * epsilonz ) )
        cbEz(:,:) = (dt / epsilonz ) / ( 1.0d0 + sigmaz * dt / (2.0d0 * epsilonz ) )

        daHx = dt/(4.0d0*mu0*mu)
        dbHx = dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
        daHx = 1.0d0 ! (1-daHz)/(1+daHz) ! 

        daHy = dt/(4.0d0*mu0*mu)
        dbHy = dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
        daHy = 1.0d0 ! (1-daHz)/(1+daHz) ! 

        daHz = dt/(4.0d0*mu0*mu)
        dbHz = dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
        daHz = 1.0d0 ! (1-daHz)/(1+daHz) ! 


        ! ----------------------------------------------------------------------
        !---
        !--- program starts here
        !---
        
        ! ================================ LOAD SOURCE ================================

        call loadsource('electromagneticsourcex.dat', source%time_steps, srcx)
        call loadsource('electromagneticsourcey.dat', source%time_steps, srcy)
        call loadsource('electromagneticsourcez.dat', source%time_steps, srcz)

        ! =============================================================================

        !--- define profile of absorption in PML region

        ! Initialize CPML damping variables
        K_x(:) = 1.0d0
        K_x_half(:) = 1.0d0
        alpha_x(:) = 0.0d0
        alpha_x_half(:) = 0.0d0
        a_x(:) = 0.0d0
        a_x_half(:) = 0.0d0
        b_x(:) = 0.0d0 
        b_x_half(:) = 0.0d0 

        K_y(:) = 1.0d0
        K_y_half(:) = 1.0d0
        alpha_y(:) = 0.0d0
        alpha_y_half(:) = 0.0d0
        a_y(:) = 0.0d0
        a_y_half(:) = 0.0d0
        b_y(:) = 0.d0
        K_z(:) = 1.0d0
        K_z_half(:) = 1.0d0
        alpha_z(:) = 0.0d0
        alpha_z_half(:) = 0.0d0
        a_z(:) = 0.0d0
        a_z_half(:) = 0.0d0

        ! ------------------------------ Load the boundary ----------------------------
        call loadcpml('kappax_cpml.dat', K_x)
        call loadcpml('alphax_cpml.dat', alpha_x)
        call loadcpml('acoefx_cpml.dat', a_x)
        call loadcpml('bcoefx_cpml.dat', b_x)

        call loadcpml('kappay_cpml.dat', K_y)
        call loadcpml('alphay_cpml.dat', alpha_y)
        call loadcpml('acoefy_cpml.dat', a_y)
        call loadcpml('bcoefy_cpml.dat', b_y)

        call loadcpml('kappaz_cpml.dat', K_z)
        call loadcpml('alphaz_cpml.dat', alpha_z)
        call loadcpml('acoefz_cpml.dat', a_z)
        call loadcpml('bcoefz_cpml.dat', b_z)

        call loadcpml('kappax_half_cpml.dat', K_x_half)
        call loadcpml('alphax_half_cpml.dat', alpha_x_half)
        call loadcpml('acoefx_half_cpml.dat', a_x_half)
        call loadcpml('bcoefx_half_cpml.dat', b_x_half)

        call loadcpml('kappay_half_cpml.dat', K_y_half)
        call loadcpml('alphay_half_cpml.dat', alpha_y_half)
        call loadcpml('acoefy_half_cpml.dat', a_y_half)
        call loadcpml('bcoefy_half_cpml.dat', b_y_half)

        call loadcpml('kappaz_half_cpml.dat', K_z_half)
        call loadcpml('alphaz_half_cpml.dat', alpha_z_half)
        call loadcpml('acoefz_half_cpml.dat', a_z_half)
        call loadcpml('bcoefz_half_cpml.dat', b_z_half)

        ! do i = 1,nz
        !   print *, K_z(i), alpha_z(i), a_z(i), b_z(i)
        ! enddo

        ! -----------------------------------------------------------------------------
        ! Load initial conditions
        call material_rw3('initialconditionEx.dat', Ex, .TRUE.)
        call material_rw3('initialconditionEy.dat', Ey, .TRUE.)
        call material_rw3('initialconditionEz.dat', Ez, .TRUE.)
        
        ! call material_rw3('initialconditionHx.dat', Hx, .TRUE.)
        ! call material_rw3('initialconditionHy.dat', Hy, .TRUE.)
        ! call material_rw3('initialconditionHz.dat', Hz, .TRUE.)
        Hx(:,:,:) = 0.d0  
        Hy(:,:,:) = 0.d0 
        Hz(:,:,:) = 0.d0 
        
        ! PML
        memory_dEx_dy(:,:,:) = 0.0d0
        memory_dEy_dx(:,:,:) = 0.0d0
        memory_dEx_dz(:,:,:) = 0.0d0
        memory_dEz_dx(:,:,:) = 0.0d0
        memory_dEz_dy(:,:,:) = 0.0d0
        memory_dEy_dz(:,:,:) = 0.0d0

        memory_dHz_dx(:,:,:) = 0.0d0
        memory_dHx_dz(:,:,:) = 0.0d0
        memory_dHz_dy(:,:,:) = 0.0d0
        memory_dHy_dz(:,:,:) = 0.0d0
        memory_dHx_dy(:,:,:) = 0.0d0
        memory_dHy_dx(:,:,:) = 0.0d0

        ! ---
        ! ---  beginning of time loop
        ! ---
        do it = 1,source%time_steps
            !--------------------------------------------------------
            ! compute magnetic field and update memory variables for C-PML
            !--------------------------------------------------------
            ! Update Hx
            do k = 1,nz-1
                do i = 1,nx-1  
                    do j = 1,ny-1
                        ! Values needed for the magnetic field updates
                        dEz_dy = ( Ez(i,j,k) - Ez(i,j+1,k) )/dy
                        memory_dEz_dy(i,j,k) = b_y_half(j) * memory_dEz_dy(i,j,k) + a_y_half(j) * dEz_dy
                        dEz_dy = dEz_dy/ K_y_half(j) + memory_dEz_dy(i,j,k)

                        ! The rest of the equation needed for agnetic field updates
                        dEy_dz = ( Ey(i,j,k+1) - Ey(i,j,k) )/dz
                        memory_dEy_dz(i,j,k) = b_z_half(k) * memory_dEy_dz(i,j,k) + a_z_half(k) * dEy_dz
                        dEy_dz = dEy_dz/ K_z_half(k) + memory_dEy_dz(i,j,k)

                        ! Now update the Magnetic field
                        Hx(i,j,k) = daHx*Hx(i,j,k) + dbHx*( dEy_dz + dEz_dy )
                    enddo
                enddo  
            enddo

                ! Update Hy
            do k = 1,nz-1
                do i = 1,nx-1      
                    do j = 1,ny-1
                    
                        ! Values needed for the magnetic field updates
                        dEx_dz = ( Ex(i,j,k) - Ex(i,j,k+1) )/dz
                        memory_dEx_dz(i,j,k) = b_z(k) * memory_dEx_dz(i,j,k) + &
                            a_z(k) * dEx_dz
                        dEx_dz = dEx_dz/ K_z(k) + memory_dEx_dz(i,j,k)

                        ! The rest of the equation needed for agnetic field updates
                        dEz_dx = ( Ez(i+1,j,k) - Ez(i,j,k) )/dx
                        memory_dEz_dx(i,j,k) = b_x(i) * memory_dEz_dx(i,j,k) + &
                            a_x(i) * dEz_dx
                        dEz_dx = dEz_dx/ K_x(i) + memory_dEz_dx(i,j,k)

                        ! Now update the Magnetic field
                        Hy(i,j,k) = daHy*Hy(i,j,k) + dbHy*( dEz_dx + dEx_dz )

                    enddo
                enddo  
            enddo

                ! Update Hz
            do k = 2,nz-1
                do i = 1,nx-1      
                    do j = 1,ny-1
                        ! Values needed for the magnetic field updates
                        dEx_dy = ( Ex(i,j+1,k) - Ex(i,j,k) )/dy
                        memory_dEx_dy(i,j,k) = b_y(j) * memory_dEx_dy(i,j,k) + & 
                            a_y(j) * dEx_dy
                        dEx_dy = dEx_dy/ K_y(j) + memory_dEx_dy(i,j,k)

                        ! The rest of the equation needed for agnetic field updates
                        dEy_dx = ( Ey(i,j,k) - Ey(i+1,j,k) )/dx
                        memory_dEy_dx(i,j,k) = b_x(i) * memory_dEy_dx(i,j,k) + & 
                            a_x(i) * dEy_dx
                        dEy_dx = dEy_dx/ K_x(i) + memory_dEy_dx(i,j,k)

                        ! Now update the Magnetic field
                        Hz(i,j,k) = daHz*Hz(i,j,k) + dbHz*( dEy_dx + dEx_dy )
                    enddo
                enddo  
            enddo

            !--------------------------------------------------------
            ! compute electric field and update memory variables for C-PML
            !--------------------------------------------------------
            ! Compute the differences in the x-direction
            do k = 2,nz-1
                do i = 1,nx-1
                    do j = 2,ny-1  
                        ! Update the Ex field
                        dHz_dy = ( Hz(i,j,k) - Hz(i,j-1,k) )/dy
                        memory_dHz_dy(i,j,k) = b_y_half(j) * memory_dHz_dy(i,j,k) + & 
                            a_y_half(j) * dHz_dy
                        dHz_dy = dHz_dy/K_y_half(j) + memory_dHz_dy(i,j,k)

                        ! Changed from half to full node positions 
                        dHy_dz = ( Hy(i,j,k-1) - Hy(i,j,k) )/dz
                        memory_dHy_dz(i,j,k) = b_z(k) * memory_dHy_dz(i,j,k) + &
                            a_z(k) * dHy_dz
                        dHy_dz = dHy_dz/K_z(k) + memory_dHy_dz(i,j,k)
                        
                        Ex(i,j,k) = caEx(i,k)*Ex(i,j,k) + & 
                        cbEx(i,k)*(dHz_dy + dHy_dz) 
                    enddo
                enddo

                ! ! Compute the differences in the y-direction
                do i = 2,nx-1 
                    do j = 1,ny-1 
                        ! Update the Ey field
                        dHz_dx = ( Hz(i-1,j,k) - Hz(i,j,k) )/dx ! this is ny-1 length vector
                        memory_dHz_dx(i,j,k) = b_x_half(i) * memory_dHz_dx(i,j,k) + & 
                            a_x_half(i) * dHz_dx
                        dHz_dx = dHz_dx/K_x_half(i) + memory_dHz_dx(i,j,k)

                        dHx_dz = ( Hx(i,j,k) - Hx(i,j,k-1) )/dz ! this is ny-1 length vector
                        memory_dHx_dz(i,j,k) = b_z_half(k) * memory_dHx_dz(i,j,k) + &
                            a_z_half(k) * dHx_dz
                        dHx_dz = dHx_dz/K_z_half(k) + memory_dHx_dz(i,j,k)

                        ! Ey(i,j,k) = ( ( 4*caEy(i,k) + caEy(i-1,k) + caEy(i,k-1) )/6) * Ey(i,j,k) + & 
                        ! ( ( 4*cbEy(i,k) + cbEy(i-1,k) + cbEy(i,k-1) )/6 ) * & 
                        ! (dHz_dx + dHx_dz)
                        Ey(i,j,k) = caEy(i,k) * Ey(i,j,k) + cbEy(i,k) * (dHz_dx + dHx_dz)
                    enddo
                enddo
            enddo 

                ! Compute the differences in the z-direction
            do k = 1,nz-1
                do i = 2,nx-1  
                    do j = 2,ny-1
                        ! Update the Ez field
                        dHx_dy = ( Hx(i,j-1,k) - Hx(i,j,k) )/dy
                        memory_dHx_dy(i,j,k) = b_y_half(j) * memory_dHx_dy(i,j,k) + &
                            a_y_half(j) * dHx_dy
                        dHx_dy = dHx_dy/K_y_half(j) + memory_dHx_dy(i,j,k)

                        dHy_dx = ( Hy(i,j,k) - Hy(i-1,j,k) )/dx
                        memory_dHy_dx(i,j,k) = b_x_half(i) * memory_dHy_dx(i,j,k) + &
                            a_x_half(i) * dHy_dx
                        dHy_dx = dHy_dx/K_x_half(i) + memory_dHy_dx(i,j,k)
                        
                        Ez(i,j,k) = ( ( 4*caEz(i,k) + caEz(i-1,k) + caEz(i,k+1) )/6 ) * &
                            Ez(i,j,k) + ( ( 4*cbEz(i,k) + cbEz(i-1,k) + cbEz(i,k+1) )/6 ) * & 
                        (dHx_dy + dHy_dx)
                    enddo
                enddo
            enddo


            ! add the source (force vector located at a given grid point)
            Ex(isource,jsource,ksource) = Ex(isource,jsource,ksource) + & 
                        srcx(it) * dt / eps11(isource,ksource)
            Ey(isource,jsource,ksource) = Ey(isource,jsource,ksource) + & 
                        srcy(it) * dt / eps22(isource,ksource) 
            Ez(isource,jsource,ksource) = Ez(isource,jsource,ksource) + & 
                        srcz(it) * dt / eps33(isource,ksource)
            
            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            Ex(1,:,:) = 0.0d0
            Ex(:,1,:) = 0.0d0
            Ex(:,:,1) = 0.0d0
            Ex(nx,:,:) = 0.0d0
            Ex(:,ny,:) = 0.0d0
            Ex(:,:,nz) = 0.0d0 

            Ey(1,:,:) = 0.0d0
            Ey(:,1,:) = 0.0d0
            Ey(:,:,1) = 0.0d0
            Ey(nx,:,:) = 0.0d0
            Ey(:,ny,:) = 0.0d0
            Ey(:,:,nz) = 0.0d0
            
            Ez(1,:,:) = 0.0d0
            Ez(:,1,:) = 0.0d0
            Ez(:,:,1) = 0.0d0
            Ez(nx,:,:) = 0.0d0
            Ez(:,ny,:) = 0.0d0
            Ez(:,:,nz) = 0.0d0
            
            Hx(1,:,:) = 0.0d0
            Hx(:,1,:) = 0.0d0
            Hx(:,:,1) = 0.0d0
            Hx(nx,:,:) = 0.0d0
            Hx(:,ny,:) = 0.0d0
            Hx(:,:,nz) = 0.0d0

            Hy(1,:,:) = 0.0d0
            Hy(:,1,:) = 0.0d0
            Hy(:,:,1) = 0.0d0
            Hy(nx,:,:) = 0.0d0
            Hy(:,ny,:) = 0.0d0
            Hy(:,:,nz) = 0.0d0
            
            Hz(1,:,:) = 0.0d0
            Hz(:,1,:) = 0.0d0
            Hz(:,:,1) = 0.0d0
            Hz(nx,:,:) = 0.0d0
            Hz(:,ny,:) = 0.0d0
            Hz(:,:,nz) = 0.0d0

            ! check norm of velocity to make sure the solution isn't diverging
            velocnorm = maxval(sqrt(Ex**2.0d0 + Ey**2.0d0 + Ez**2.0d0) )
            if (velocnorm > stability_threshold) stop 'code became unstable and blew up'
            ! print *,'Max vals for Ex, Ey, Ez: ', maxval(Ex), maxval(Ey), maxval(Ez)

            ! print *, maxval(Ex), maxval(Ey), maxval(Ez)
            call write_image(Ex, domain, source, it, 'Ex', SINGLE)
            call write_image(Ey, domain, source, it, 'Ey', SINGLE)
            call write_image(Ez, domain, source, it, 'Ez', SINGLE)

        enddo   ! end of time loop
    end subroutine electromag25
    
    ! ==========================================================================
    ! =========================================================================
    subroutine seismoacoustic2(domain, source, SINGLE_OUTPUT)
        use constants
        use omp_lib ! Include the OpenMP library 
        
        ! Input arguments
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT
        
        ! Local variables
        real(real64) :: deltarho, velocnorm, value_dvx_dx, value_dvx_dz, &
            value_dvz_dx, value_dvz_dz, value_dsigmaxx_dx, value_dsigmazz_dz, &
            value_dsigmaxz_dx, value_dsigmaxz_dz, dp_dx, dp_dz

        ! 1D arrays for damping profiles
        real(real64), allocatable :: c11(:,:), c13(:,:), c15(:,:), c33(:,:), c35(:,:), c55(:,:), rho(:,:)
        real(real64), allocatable :: K_x(:), alpha_x(:), a_x(:), b_x(:), K_x_half(:), alpha_x_half(:), a_x_half(:), b_x_half(:)
        real(real64), allocatable :: K_z(:), alpha_z(:), a_z(:), b_z(:), K_z_half(:), alpha_z_half(:), a_z_half(:), b_z_half(:)
        real(real64), allocatable :: gamma_x(:,:), gamma_z(:,:), gamma_xz(:,:)        
        real(real64), allocatable :: srcx(:), srcz(:) ! The vector time series of the source
        
        ! Model variables
        real(real64), allocatable :: vx(:,:),vz(:,:),sigmaxx(:,:),sigmazz(:,:),sigmaxz(:,:)
        real(real64), allocatable :: pressure(:,:)
        integer, allocatable :: material_geometry(:,:)
        real(real64), allocatable :: memory_dvx_dx(:,:), memory_dvx_dz(:,:), &
                                     memory_dvz_dx(:,:), memory_dvz_dz(:,:), &
                                     memory_dsigmaxx_dx(:,:), memory_dsigmazz_dz(:,:), &
                                     memory_dsigmaxz_dx(:,:), memory_dsigmaxz_dz(:,:)
        
        real(real64), allocatable :: pmemory_dvx_dx(:,:), pmemory_dvz_dz(:,:)
        
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
        allocate(pressure(nx,nz), material_geometry(nx,nz))
        allocate(pmemory_dvx_dx(nx,nz), pmemory_dvz_dz(nx,nz) )
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
        pressure(:,:) = 0.d0
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
        
        pmemory_dvx_dx(:,:) = 0.d0
        pmemory_dvz_dz(:,:) = 0.d0
        
        ! Load the geometry masked array 
        call read_geometry('geometry_mask.dat', material_geometry)
        
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
                    
                    value_dvx_dx = (vx(i+1,j) - vx(i,j)) / dx
                    value_dvz_dz = (vz(i,j) - vz(i,j-1)) / dz
                    
                    if ( material_geometry(i,j) == 0) then 
                        value_dvz_dx = (vz(i+1,j) - vz(i,j)) / dx
                        value_dvx_dz = (vx(i,j) - vx(i,j-1)) / dz
                        
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
                        pmemory_dvx_dx(i,j) = b_x_half(i) * pmemory_dvx_dx(i,j) + &
                                a_x_half(i) * value_dvx_dx
                        pmemory_dvz_dz(i,j) = b_z(j) * pmemory_dvz_dz(i,j) + &
                                a_z(j) * value_dvz_dz
                        pressure(i,j) = pressure(i,j) - rho(i,j) * c11(i,j) * dt * (value_dvx_dx + value_dvz_dz)
                    
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
            ! do j = 2,nz-1
            !     do i = 2,nx-1
            !         if (material_geometry(i,j) == 1) then
            !             ! Check neighbors
            !             if (material_geometry(i+1,j) == 0) then
            !                 pressure(i,j) = -sigmaxx(i+1,j)
            !                 vx(i,j) = vx(i+1,j)
            !             else if (material_geometry(i-1,j) == 0) then
            !                 pressure(i,j) = -sigmaxx(i-1,j)
            !                 vx(i,j) = vx(i-1,j)
            !             end if
            !             if (material_geometry(i,j+1) == 0) then
            !                 pressure(i,j) = -sigmazz(i,j+1)
            !                 vz(i,j) = vz(i,j+1)
            !             else if (material_geometry(i,j-1) == 0) then
            !                 pressure(i,j) = -sigmazz(i,j-1)
            !                 vz(i,j) = vz(i,j-1)
            !             end if
            !         end if
            !     end do
            ! end do
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
            if (material_geometry(isource, jsource) == 1) then
                pressure(isource, jsource) = pressure(isource, jsource) + sqrt(srcx(it)**2 + srcz(it)**2)
            else
                vx(isource, jsource) = vx(isource, jsource) + srcx(it) / rho(isource, jsource)
                vz(isource, jsource) = vz(isource, jsource) + srcz(it) / rho(isource, jsource)
            end if

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
        deallocate(pressure, material_geometry)
    end subroutine seismoacoustic2 
    
end module cpmlfdtd