module cpmlfdtd 
    
    use seidartio
    use seidart_types
    use density_averaging, only: face_density   
     
    implicit none 
    
    contains
    
    ! =========================================================================    
    ! Computations are done in double precision and written to binary as single
    ! precision unless specified by the optional logical, SINGLE_OUTPUT.
    ! ==========================================================================
    
    subroutine seismic2(domain, source, density_method, SINGLE_OUTPUT)
        ! Standard Staggered Grid Formulation
        use constants
        use omp_lib ! Include the OpenMP library 
        
        implicit none 
        
        ! Input arguments
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        character(len=256), intent(in) :: density_method
        logical, intent(in), optional :: SINGLE_OUTPUT

        ! Local variables
        real(real64) :: velocnorm, value_dvx_dx, value_dvx_dz, &
            value_dvz_dx, value_dvz_dz, value_dsigmaxx_dx, value_dsigmazz_dz, &
            value_dsigmaxz_dx, value_dsigmaxz_dz, &
            rhoxx, rhozx, rhoxz, rhozz !deltarho,  

        ! 1D arrays for damping profiles
        real(real64), allocatable :: c11(:,:), c13(:,:), c15(:,:), c33(:,:), c35(:,:), c55(:,:), rho(:,:)

        real(real64), allocatable :: kappa(:,:), alpha(:,:), acoef(:,:), bcoef(:,:)
        real(real64), allocatable :: kappa_half(:,:), alpha_half(:,:), acoef_half(:,:), bcoef_half(:,:) 
        real(real64), allocatable :: gamma_x(:,:), gamma_z(:,:), gamma_xz(:,:) 
        real(real64), allocatable :: srcx(:), srcz(:)       ! Velocity injection
        real(real64), allocatable :: srcxx(:), srczz(:), srcxz(:) ! Stress injection
        
        ! Model variables
        real(real64), allocatable :: vx(:,:),vz(:,:),sigmaxx(:,:),sigmazz(:,:),sigmaxz(:,:)
        real(real64), allocatable :: memory_dvx_dx1(:,:), memory_dvx_dz1(:,:), &
                                     memory_dvz_dx1(:,:), memory_dvz_dz1(:,:), &
                                     memory_dvx_dx2(:,:), memory_dvx_dz2(:,:), &
                                     memory_dvz_dx2(:,:), memory_dvz_dz2(:,:), &
                                     memory_dsigmaxx_dx(:,:), memory_dsigmazz_dz(:,:), &
                                     memory_dsigmaxz_dx(:,:), memory_dsigmaxz_dz(:,:)
        
        integer :: nx, nz 
        real(real64) :: dx, dz, dt 
        
        integer :: i, j, it, isource, jsource
        logical :: SINGLE
        integer :: density_code 

        ! The default data output is single precision unless SINGLE_OUTPUT is 
        ! set to .FALSE.
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        end if
        
        select case (trim(adjustl(density_method)))
            case ('none'      ); density_code = 0
            case ('harmonic'  ); density_code = 1
            case ('geometric'); density_code = 2
            case ('arithmetic'); density_code = 3
        end select
        
        nx = domain%nx 
        nz = domain%nz 
        dt = source%dt 
        dx = domain%dx 
        dz = domain%dz 
        
        ! Allocate the arrays based on runtime values of nx and nz
        allocate(c11(nx, nz), c13(nx, nz), c15(nx, nz), &
                 c33(nx, nz), c35(nx, nz), c55(nx, nz), rho(nx, nz))
        allocate(kappa(nx,nz), alpha(nx,nz), acoef(nx,nz), bcoef(nx,nz), &
                 kappa_half(nx,nz), alpha_half(nx,nz), acoef_half(nx,nz), bcoef_half(nx,nz) )
        allocate(gamma_x(nx, nz), gamma_z(nx, nz), gamma_xz(nx, nz))
        allocate(srcx(source%time_steps), srcz(source%time_steps))
        allocate(srcxx(source%time_steps), srcxz(source%time_steps),srczz(source%time_steps))
        
        
        ! Allocate more
        allocate(memory_dvx_dx1(nx, nz), memory_dvx_dz1(nx, nz))
        allocate(memory_dvz_dx1(nx, nz), memory_dvz_dz1(nx, nz))
        allocate(memory_dvx_dx2(nx, nz), memory_dvx_dz2(nx, nz))
        allocate(memory_dvz_dx2(nx, nz), memory_dvz_dz2(nx, nz))
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
        srcxx(:) = 0.0_real64
        srcxz(:) = 0.0_real64
        srczz(:) = 0.0_real64
        srcx(:) = 0.0_real64
        srcz(:) = 0.0_real64
        
        if ( source%source_type == 'ac') then 
            call loadsource('seismicsourcex.dat', source%time_steps, srcx)
            call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        else
            call loadsource('seismicsourcexx.dat', source%time_steps, srcxx)
            call loadsource('seismicsourcexz.dat', source%time_steps, srcxz)
            call loadsource('seismicsourcezz.dat', source%time_steps, srczz)
        endif 
        
        ! -----------------------------------------------------------------------------
        !--- define profile of absorption in PML region

        ! Initialize PML 
        kappa(:,:) = 1.0_real64
        kappa_half(:,:) = 1.0_real64
        alpha(:,:) = 0.0_real64
        alpha_half(:,:) = 0.0_real64
        acoef(:,:) = 0.0_real64
        acoef_half(:,:) = 0.0_real64
        bcoef(:,:) = 0.0_real64 
        bcoef_half(:,:) = 0.0_real64

        ! ------------------------------ Load the boundary ----------------------------        
        call material_rw2('kappa_cpml.dat', kappa, .TRUE.)
        call material_rw2('alpha_cpml.dat', alpha, .TRUE.)
        call material_rw2('acoef_cpml.dat', acoef, .TRUE.)
        call material_rw2('bcoef_cpml.dat', bcoef, .TRUE.)
        
        call material_rw2('kappa_half_cpml.dat', kappa_half, .TRUE.)
        call material_rw2('alpha_half_cpml.dat', alpha_half, .TRUE.)
        call material_rw2('acoef_half_cpml.dat', acoef_half, .TRUE.)
        call material_rw2('bcoef_half_cpml.dat', bcoef_half, .TRUE.)
        
        ! =============================================================================
        ! Load initial condition
        call material_rw2('initialconditionVx.dat', vx, .TRUE.)
        call material_rw2('initialconditionVz.dat', vz, .TRUE.)
        
        vx(:,:) = 0.0_real64
        vz(:,:) = 0.0_real64
        sigmaxx(:,:) = 0.0_real64
        sigmazz(:,:) = 0.0_real64
        sigmaxz(:,:) = 0.0_real64

        ! PML
        memory_dvx_dx1(:,:) = 0.0_real64
        memory_dvx_dz1(:,:) = 0.0_real64
        memory_dvz_dx1(:,:) = 0.0_real64
        memory_dvz_dz1(:,:) = 0.0_real64
        memory_dvx_dx2(:,:) = 0.0_real64
        memory_dvx_dz2(:,:) = 0.0_real64
        memory_dvz_dx2(:,:) = 0.0_real64
        memory_dvz_dz2(:,:) = 0.0_real64
        memory_dsigmaxx_dx(:,:) = 0.0_real64
        memory_dsigmazz_dz(:,:) = 0.0_real64
        memory_dsigmaxz_dx(:,:) = 0.0_real64
        memory_dsigmaxz_dz(:,:) = 0.0_real64

        !---
        !---  beginning of time loop
        !---
        
        do it = 1,source%time_steps
            !$omp parallel private(i, j, deltarho, value_dvx_dx, value_dvx_dz, value_dvz_dx, value_dvz_dz, value_dsigmaxx_dx, value_dsigmazz_dz, value_dsigmaxz_dx, value_dsigmaxz_dz)
            ! ------------------------------------------------------------
            !  compute stress sigma and update memory variables for C-PML
            ! ------------------------------------------------------------
            
            !$omp do
            do j = 3,nz-1
                do i = 2,nx-2
                    
                    value_dvx_dx = ( 27.0_real64*(vx(i+1,j) - vx(i,j)) - (vx(i+2,j) - vx(i-1,j))) / (24.0_real64*dx)
                    value_dvz_dz = ( 27.0_real64*(vz(i,j) - vz(i,j-1)) - (vz(i,j+1) - vz(i,j-2))) / (24.0_real64*dz)
                    value_dvz_dx = ( 27.0_real64*(vz(i+1,j) - vz(i,j)) - (vz(i+2,j) - vz(i-1,j))) / (24.0_real64*dx)
                    value_dvx_dz = ( 27.0_real64*(vx(i,j) - vx(i,j-1)) - (vx(i,j+1) - vx(i,j-2))) / (24.0_real64*dz)

                    memory_dvx_dx1(i,j) = bcoef_half(i,j) * memory_dvx_dx1(i,j) + acoef_half(i,j) * value_dvx_dx
                    memory_dvx_dz1(i,j) = bcoef_half(i,j) * memory_dvx_dz1(i,j) + acoef_half(i,j) * value_dvx_dz 
                    memory_dvz_dz1(i,j) = bcoef(i,j) * memory_dvz_dz1(i,j) + acoef(i,j) * value_dvz_dz
                    memory_dvz_dx1(i,j) = bcoef(i,j) * memory_dvz_dx1(i,j) + acoef(i,j) * value_dvz_dx
                    
                    value_dvx_dx = value_dvx_dx / kappa_half(i,j) + memory_dvx_dx1(i,j)
                    value_dvx_dz = value_dvx_dz / kappa_half(i,j) + memory_dvx_dz1(i,j)
                    value_dvz_dz = value_dvz_dz / kappa(i,j) + memory_dvz_dz1(i,j)
                    value_dvz_dx = value_dvz_dx / kappa(i,j) + memory_dvz_dx1(i,j)
                    
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
            do j = 2,nz-2
                do i = 3,nx-1
                    
                    value_dvx_dx = ( 27.0_real64*(vx(i,j) - vx(i-1,j)) - (vx(i+1,j) - vx(i-2,j))) / (24.0_real64*dx)
                    value_dvz_dz = ( 27.0_real64*(vz(i,j+1) - vz(i,j)) - (vz(i,j+2) - vz(i,j-1))) / (24.0_real64*dz)
                    value_dvz_dx = ( 27.0_real64*(vz(i,j) - vz(i-1,j)) - (vz(i+1,j) - vz(i-2,j))) / (24.0_real64*dx)
                    value_dvx_dz = ( 27.0_real64*(vx(i,j+1) - vx(i,j)) - (vx(i,j+2) - vx(i,j-1))) / (24.0_real64*dz)

                    memory_dvx_dx2(i,j) = bcoef_half(i,j) * memory_dvx_dx2(i,j) + acoef_half(i,j) * value_dvx_dx
                    memory_dvx_dz2(i,j) = bcoef_half(i,j) * memory_dvx_dz2(i,j) + acoef_half(i,j) * value_dvx_dz 
                    memory_dvz_dz2(i,j) = bcoef_half(i,j) * memory_dvz_dz2(i,j) + acoef_half(i,j) * value_dvz_dz
                    memory_dvz_dx2(i,j) = bcoef_half(i,j) * memory_dvz_dx2(i,j) + acoef_half(i,j) * value_dvz_dx
                    
                    value_dvx_dx = value_dvx_dx / kappa_half(i,j) + memory_dvx_dx2(i,j)
                    value_dvx_dz = value_dvx_dz / kappa_half(i,j) + memory_dvx_dz2(i,j)
                    value_dvz_dz = value_dvz_dz / kappa_half(i,j) + memory_dvz_dz2(i,j)
                    value_dvz_dx = value_dvz_dx / kappa_half(i,j) + memory_dvz_dx2(i,j)
                    
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
            do j = 3,nz-1
                do i = 3,nx-1
                    
                    rhoxx = face_density(rho(i,j), rho(i-1,j), density_code) 
                    rhozx = face_density(rho(i,j), rho(i,j-1), density_code) 

                    value_dsigmaxx_dx = (27.0_real64*sigmaxx(i,j) - 27.0_real64*sigmaxx(i-1,j) + sigmaxx(i-2,j)) / (24.0_real64*dx)
                    value_dsigmaxz_dz = (27.0_real64*sigmaxz(i,j) - 27.0_real64*sigmaxz(i,j-1) + sigmaxz(i,j-2)) / (24.0_real64*dz)
                    
                    memory_dsigmaxx_dx(i,j) = bcoef(i,j) * memory_dsigmaxx_dx(i,j) + acoef(i,j) * value_dsigmaxx_dx
                    memory_dsigmaxz_dz(i,j) = bcoef(i,j) * memory_dsigmaxz_dz(i,j) + acoef(i,j) * value_dsigmaxz_dz
                    
                    value_dsigmaxx_dx = value_dsigmaxx_dx / kappa(i,j) + memory_dsigmaxx_dx(i,j)
                    value_dsigmaxz_dz = value_dsigmaxz_dz / kappa(i,j) + memory_dsigmaxz_dz(i,j)
                    vx(i,j) = vx(i,j) + dt * (value_dsigmaxx_dx/rhoxx + value_dsigmaxz_dz/rhozx ) !deltarho
                    
                enddo
            enddo
            !$omp end do
            
            !$omp do
            do j = 2,nz-2
                do i = 2,nx-2
                    
                    rhoxz = face_density(rho(i,j), rho(i+1,j), density_code) 
                    rhozz = face_density(rho(i,j), rho(i,j+1), density_code) 

                    value_dsigmaxz_dx = (-27.0_real64*sigmaxz(i,j) + 27.0_real64*sigmaxz(i+1,j) - sigmaxz(i+2,j)) / (24.0_real64*dx)
                    value_dsigmazz_dz = (-27.0_real64*sigmazz(i,j) + 27.0_real64*sigmazz(i,j+1) - sigmazz(i,j+2)) / (24.0_real64*dz)

            
                    memory_dsigmaxz_dx(i,j) = bcoef_half(i,j) * memory_dsigmaxz_dx(i,j) + acoef_half(i,j) * value_dsigmaxz_dx
                    memory_dsigmazz_dz(i,j) = bcoef_half(i,j) * memory_dsigmazz_dz(i,j) + acoef_half(i,j) * value_dsigmazz_dz
            
                    value_dsigmaxz_dx = value_dsigmaxz_dx / kappa_half(i,j) + memory_dsigmaxz_dx(i,j)
                    value_dsigmazz_dz = value_dsigmazz_dz / kappa_half(i,j) + memory_dsigmazz_dz(i,j)

                    vz(i,j) = vz(i,j) + dt * (value_dsigmaxz_dx/rhoxz + value_dsigmazz_dz/rhozz) !deltarho
                enddo
            enddo
            !$omp end do
            
            !$omp end parallel
            
            ! Add the source term. If it is an accelerated weight drop the src_ij terms are zero
            ! If it is any other source, the src_i terms are zero
            sigmaxx(isource,jsource) = sigmaxx(isource, jsource) + srcxx(it) / rho(isource,jsource)
            sigmaxz(isource+1,jsource+1) = sigmaxz(isource+1, jsource+1) + srcxz(it) / rho(isource+1,jsource+1)  
            sigmazz(isource,jsource) = sigmaxz(isource, jsource) + srczz(it) / rho(isource,jsource) 
            vx(isource+1,jsource) = vx(isource+1,jsource) + srcx(it) / rho(isource+1,jsource)
            vz(isource,jsource+1) = vz(isource,jsource+1) + srcz(it) / rho(isource,jsource+1)
        
            ! Dirichlet conditions (rigid boundaries) on the edges or at the 
            ! bottom of the PML layers
            vx(1,:) = 0.0_real64
            vx(nx,:) = 0.0_real64
        
            vx(:,1) = 0.0_real64
            vx(:,nz) = 0.0_real64
        
            vz(1,:) = 0.0_real64
            vz(nx,:) = 0.0_real64
        
            vz(:,1) = 0.0_real64
            vz(:,nz) = 0.0_real64
        
            ! print maximum of norm of velocity
            velocnorm = maxval(sqrt(vx**2 + vz**2))
            if (velocnorm > stability_threshold) stop 'code became unstable and blew up'
        
            call write_image2(vx, nx, nz, source, it, 'Vx', SINGLE)
            call write_image2(vz, nx, nz, source, it, 'Vz', SINGLE)

        enddo   ! end of time loop
        
        deallocate(c11, c13, c15, c33, c35, c55, rho)
        deallocate(kappa, alpha, acoef, bcoef, kappa_half, alpha_half, acoef_half, bcoef_half)
        deallocate(gamma_x, gamma_z, gamma_xz, srcxx, srcxz, srczz)
        deallocate(srcx, srcz)
        deallocate(srcxx, srcxz, srczz)
        deallocate(memory_dvx_dx1, memory_dvx_dz1)
        deallocate(memory_dvz_dx1, memory_dvz_dz1)
        deallocate(memory_dvx_dx2, memory_dvx_dz2)
        deallocate(memory_dvz_dx2, memory_dvz_dz2)
        deallocate(memory_dsigmaxx_dx, memory_dsigmazz_dz)
        deallocate(memory_dsigmaxz_dx, memory_dsigmaxz_dz)
        deallocate(vx, vz, sigmaxx, sigmazz, sigmaxz)
        
    end subroutine seismic2
    
    ! ==========================================================================    
    subroutine seismic2iso(domain, source, density_method, SINGLE_OUTPUT)
        !-----------------------------------------------------------------------
        use constants
        
        implicit none

        ! Input arguments
        type(Domain_Type), intent(in) :: domain 
        type(Source_Type), intent(in) :: source 
        logical, intent(in), optional :: SINGLE_OUTPUT
        character(len=256), intent(in) :: density_method
        ! Local variables
        real(real64) :: rhoxx, rhozx, &
                        rhoxz, rhozz
        
        real(real64), allocatable :: velocnorm(:), stressnorm(:)

        real(real64), allocatable :: c11(:,:), c13(:,:), c33(:,:), c55(:,:), rho(:,:)
        ! real(real64) :: DT
        integer :: i, k, it, isource, ksource

        ! Values of the velocity and stress differentials
        real(real64) :: dvx_dx, dvx_dz, &
                        dvz_dx, dvz_dz, &
                        dsigmaxx_dx, dsigmazz_dz, &
                        dsigmaxz_dx, dsigmaxz_dz

        ! 1D arrays for the damping profiles in each direction
        real(real64), allocatable :: kappa(:,:), alpha(:,:), acoef(:,:), bcoef(:,:), &
                        kappa_half(:,:), alpha_half(:,:), acoef_half(:,:), bcoef_half(:,:)
        
        ! Arrays for the PML damping factors
        real(real64), allocatable :: gamma_x(:,:), gamma_z(:,:), gamma_xz(:,:)

        ! Source arrays
        real(real64), allocatable :: srcx(:), srcz(:)
        real(real64), allocatable :: srcxx(:), srcxz(:), srczz(:)
        
        real(real64), allocatable :: vx(:,:), vz(:,:), &
                        sigmaxx(:,:), sigmaxz(:,:), sigmazz(:,:)
        real(real64), allocatable :: &
                memory_dvx_dx(:,:),  memory_dvx_dz(:,:), &
                memory_dvz_dx(:,:),  memory_dvz_dz(:,:), &
                memory_dsigmaxx_dx(:,:),  memory_dsigmazz_dz(:,:), &
                memory_dsigmaxz_dx(:,:), memory_dsigmaxz_dz(:,:)
        
        integer :: nx, nz
        real(real64) :: dx, dz, dt    
        
        ! Boolean flag to save as double precision or single precision
        logical :: SINGLE
        integer :: density_code
        
        ! The default data output is single precision unless SINGLE_OUTPUT is 
        ! set to .FALSE.
        if (present(SINGLE_OUTPUT)) then 
            SINGLE = SINGLE_OUTPUT 
        else
            SINGLE = .TRUE.
        endif
        
        select case (trim(adjustl(density_method)))
            case ('none'      ); density_code = 0
            case ('harmonic'  ); density_code = 1
            case ('geometric'); density_code = 2
            case ('arithmetic'); density_code = 3
        end select
        
        nx = domain%nx
        nz = domain%nz 
        dt = source%dt 
        dx = domain%dx
        dz = domain%dz
        
        ! ----------------------- Allocate Arrays ----------------------------
        allocate(c11(nx, nz), c13(nx, nz), &
                 c33(nx, nz), c55(nx, nz) )
        allocate(rho(nx,nz))
                
        allocate(kappa(nx, nz), alpha(nx, nz), acoef(nx, nz), bcoef(nx, nz), &
                 kappa_half(nx, nz), alpha_half(nx, nz), acoef_half(nx, nz), &
                 bcoef_half(nx, nz))
        allocate(gamma_x(nx, nz), gamma_z(nx, nz), gamma_xz(nx, nz))
        
        allocate(velocnorm(source%time_steps), stressnorm(source%time_steps))
        
        allocate(srcx(source%time_steps), srcz(source%time_steps))
        allocate(srcxx(source%time_steps), srcxz(source%time_steps), &
                 srczz(source%time_steps) )                
        
        ! Allocate more
        allocate(memory_dvx_dx(nx,nz), memory_dvx_dz(nx,nz) )
        allocate(memory_dvz_dx(nx,nz), memory_dvz_dz(nx,nz) )
        
        allocate(memory_dsigmaxx_dx(nx,nz), memory_dsigmazz_dz(nx,nz) )
        allocate(memory_dsigmaxz_dx(nx,nz), memory_dsigmaxz_dz(nx,nz))
        
        
        allocate(vx(nx, nz), vz(nx, nz))
        allocate(sigmaxx(nx, nz), sigmaxz(nx, nz), sigmazz(nx, nz))
                
        
        ! --------------------- Load Stiffness Coefficients --------------------
            call material_rw2('c11.dat', c11, .TRUE.)
            call material_rw2('c13.dat', c13, .TRUE.)
            call material_rw2('c33.dat', c33, .TRUE.)
            call material_rw2('c55.dat', c55, .TRUE.)
            call material_rw2('rho.dat', rho, .TRUE.)
        
        ! ------------------- Load Attenuation Coefficients --------------------
            call material_rw2('gamma_x.dat', gamma_x, .TRUE.)
            call material_rw2('gamma_z.dat', gamma_z, .TRUE.)
            call material_rw2('gamma_xz.dat', gamma_xz, .TRUE.)

        ! ------------------------ Assign some constants -----------------------
        isource = source%xind + domain%cpml
        ksource = source%zind + domain%cpml

        ! ========================== LOAD SOURCE ===============================
        srcxx(:) = 0.0_real64
        srcxz(:) = 0.0_real64
        srczz(:) = 0.0_real64
        srcx(:) = 0.0_real64
        srcz(:) = 0.0_real64
        
        if ( source%source_type == 'ac') then 
            call loadsource('seismicsourcex.dat', source%time_steps, srcx)
            call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        else
            call loadsource('seismicsourcexx.dat', source%time_steps, srcxx)
            call loadsource('seismicsourcexz.dat', source%time_steps, srcxz)
            call loadsource('seismicsourcezz.dat', source%time_steps, srczz)
        endif 
        
        ! =============================== PML ==================================
        ! Initialize PML 
        kappa(:,:) = 1.0_real64
        kappa_half(:,:) = 1.0_real64
        alpha(:,:) = 0.0_real64
        alpha_half(:,:) = 0.0_real64
        acoef(:,:) = 0.0_real64
        acoef_half(:,:) = 0.0_real64

        ! ------------------------- Boundary Conditions ------------------------
        call material_rw2('kappa_cpml.dat', kappa, .TRUE.)
        call material_rw2('alpha_cpml.dat', alpha, .TRUE.)
        call material_rw2('acoef_cpml.dat', acoef, .TRUE.)
        call material_rw2('bcoef_cpml.dat', bcoef, .TRUE.)

        call material_rw2('kappa_half_cpml.dat', kappa_half, .TRUE.)
        call material_rw2('alpha_half_cpml.dat', alpha_half, .TRUE.)
        call material_rw2('acoef_half_cpml.dat', acoef_half, .TRUE.)
        call material_rw2('bcoef_half_cpml.dat', bcoef_half, .TRUE.)

        ! Load initial condition
        call material_rw2('initialconditionVx.dat', vx, .TRUE.)
        call material_rw2('initialconditionVz.dat', vz, .TRUE.)

        ! Initialize the stress values
        sigmaxx(:,:) = 0.0_real64
        sigmazz(:,:) = 0.0_real64
        sigmaxz(:,:) = 0.0_real64

        ! PML
        memory_dvx_dx(:,:) = 0.0_real64
        memory_dvx_dz(:,:) = 0.0_real64
        memory_dvz_dx(:,:) = 0.0_real64
        memory_dvz_dz(:,:) = 0.0_real64

        memory_dsigmaxx_dx(:,:) = 0.0_real64
        memory_dsigmazz_dz(:,:) = 0.0_real64
        memory_dsigmaxz_dx(:,:) = 0.0_real64
        memory_dsigmaxz_dz(:,:) = 0.0_real64

        
        ! Initialize velocnorm and stressnorm 
        velocnorm(:) = 0.0_real64
        stressnorm(:) = 0.0_real64 
        ! Do it 
        
        ! =========================== Forward Model ============================
        do it = 1,source%time_steps
            !------------------------------------------------------------
            ! compute stress sigma and update memory variables for C-PML
            !------------------------------------------------------------
            ! Update in the x direction
            ! Loop 1
            do k = 3,nz-1
                do i = 2,nx-2
                    ! Forward step on half grid
                    dvx_dx = ( 27.0_real64*(vx(i+1,k) - vx(i,k)) - (vx(i+2,k) - vx(i-1,k))) / (24.0_real64*dx)
                
                    ! Backward step on full grid
                    dvz_dz = ( 27.0_real64*(vz(i,k) - vz(i,k-1)) - (vz(i,k+1) - vz(i,k-2))) / (24.0_real64*dz)
                    
                    memory_dvx_dx(i,k) = bcoef_half(i,k) * memory_dvx_dx(i,k) + acoef_half(i,k) * dvx_dx
                    memory_dvz_dz(i,k) = bcoef(i,k) * memory_dvz_dz(i,k) + acoef(i,k) * dvz_dz
                    
                    dvx_dx = dvx_dx / kappa_half(i,k) + memory_dvx_dx(i,k)
                    dvz_dz = dvz_dz / kappa(i,k) + memory_dvz_dz(i,k)
                    
                    sigmaxx(i,k) = ( sigmaxx(i,k) + &
                    ( c11(i,k) * dvx_dx + c13(i,k) * dvz_dz ) * dt ) / &
                            (1 + gamma_x(i,k) * dt )

                    sigmazz(i,k) = ( sigmazz(i,k) + &
                    ( c13(i,k) * dvx_dx + c33(i,k) * dvz_dz ) * dt ) / &
                            (1 + gamma_z(i,k) * dt )
                enddo
            enddo
            ! Loop 2
            ! Update sigmaxz, z-direction is full nodes
            do k = 2,nz-2
                do i = 3,nx-1
                    ! Backward difference full grid
                    dvz_dx = ( 27.0_real64*(vz(i,k) - vz(i-1,k)) - (vz(i+1,k) - vz(i-2,k))) / (24.0_real64*dx)
                    ! Forward difference half grid
                    dvx_dz = ( 27.0_real64*(vx(i,k+1) - vx(i,k)) - (vx(i,k+2) - vx(i,k-1))) / (24.0_real64*dz)
                    
                    memory_dvz_dx(i,k) = bcoef(i,k) * memory_dvz_dx(i,k) + acoef(i,k) * dvz_dx
                    memory_dvx_dz(i,k) = bcoef_half(i,k) * memory_dvx_dz(i,k) + acoef_half(i,k) * dvx_dz

                    dvz_dx = dvz_dx / kappa(i,k) + memory_dvz_dx(i,k) 
                    dvx_dz = dvx_dz / kappa_half(i,k) + memory_dvx_dz(i,k)
                    
                    sigmaxz(i,k) = ( sigmaxz(i,k) + &
                        (  c55(i,k) * ( dvx_dz + dvz_dx) ) * dt  ) / &
                                (1 + gamma_xz(i,k) * dt )

                enddo
            enddo

            !--------------------------------------------------------
            ! compute velocity and update memory variables for C-PML
            !--------------------------------------------------------
            ! Loop 5
            do k = 3,nz-1
                do i = 3,nx-1
                    ! ds1/dx, ds6/dy, ds5,dz
                    rhoxx = face_density(rho(i,k), rho(i-1,k), density_code)
                    rhozx = face_density(rho(i,k), rho(i,k-1), density_code) 
                    ! Backward difference half grid
                    dsigmaxx_dx = (27.0_real64*sigmaxx(i,k) - 27.0_real64*sigmaxx(i-1,k) + sigmaxx(i-2,k)) / (24.0_real64*dx)
                    dsigmaxz_dz = (27.0_real64*sigmaxz(i,k) - 27.0_real64*sigmaxz(i,k-1) + sigmaxz(i,k-2)) / (24.0_real64*dz)
                    
                    memory_dsigmaxx_dx(i,k) = bcoef_half(i,k) * memory_dsigmaxx_dx(i,k) + acoef_half(i,k) * dsigmaxx_dx
                    memory_dsigmaxz_dz(i,k) = bcoef_half(i,k) * memory_dsigmaxz_dz(i,k) + acoef_half(i,k) * dsigmaxz_dz

                    dsigmaxx_dx = dsigmaxx_dx / kappa_half(i,k) + memory_dsigmaxx_dx(i,k)
                    dsigmaxz_dz = dsigmaxz_dz / kappa_half(i,k) + memory_dsigmaxz_dz(i,k) 

                    vx(i,k) = vx(i,k) + &
                        (dsigmaxx_dx/rhoxx + dsigmaxz_dz/rhozx) * dt 
                enddo
            enddo

            do k = 2,nz-2
                do i = 2,nx-2
                    ! ds5/dx, ds4/dy, ds3/dz
                    rhoxz = face_density(rho(i,k), rho(i+1,k), density_code)
                    rhozz = face_density(rho(i,k), rho(i,k+1), density_code)
                    
                    ! Forward difference half grid
                    dsigmaxz_dx = (-27.0_real64*sigmaxz(i,k) + 27.0_real64*sigmaxz(i+1,k) - sigmaxz(i+2,k)) / (24.0_real64*dx)
                    ! Forward difference half grid
                    dsigmazz_dz = (-27.0_real64*sigmazz(i,k) + 27.0_real64*sigmazz(i,k+1) - sigmazz(i,k+2)) / (24.0_real64*dz)
                    
                    memory_dsigmaxz_dx(i,k) = bcoef_half(i,k) * memory_dsigmaxz_dx(i,k) + acoef_half(i,k) * dsigmaxz_dx
                    memory_dsigmazz_dz(i,k) = bcoef_half(i,k) * memory_dsigmazz_dz(i,k) + acoef_half(i,k) * dsigmazz_dz

                    dsigmaxz_dx = dsigmaxz_dx / kappa_half(i,k) + memory_dsigmaxz_dx(i,k)
                    dsigmazz_dz = dsigmazz_dz / kappa_half(i,k) + memory_dsigmazz_dz(i,k)

                    vz(i,k) = vz(i,k) + &
                        (dsigmaxz_dx/rhoxz + dsigmazz_dz/rhozz) * &
                        dt !/ deltarho !rho(i,k)

                enddo
            enddo
            
            sigmaxx(isource,ksource) = sigmaxx(isource,ksource) + srcxx(it) / rho(isource,ksource)  
            sigmaxz(isource+1,ksource+1) = sigmaxz(isource+1,ksource+1) + srcxz(it) / rho(isource+1,ksource)  
            sigmazz(isource,ksource) = sigmazz(isource,ksource) + srczz(it) / rho(isource,ksource)
            vx(isource+1,ksource) = vx(isource+1,ksource) + &
                    srcx(it) * dt / rho(isource+1,ksource)
            vz(isource,ksource+1) = vz(isource,ksource+1) + &
                    srcz(it) * dt / rho(isource,ksource+1)

            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            vx(1,:) = 0.0_real64
            vz(1,:) = 0.0_real64

            vx(:,1) = 0.0_real64
            vz(:,1) = 0.0_real64

            vx(nx,:) = 0.0_real64
            vz(nx,:) = 0.0_real64

            vx(:,nz) = 0.0_real64
            vz(:,nz) = 0.0_real64

            ! check norm of velocity to make sure the solution isn't diverging
            velocnorm(it) = maxval( sqrt(vx**2 + vz**2) )
            stressnorm(it) = maxval( &
                        sqrt(sigmaxx**2 + sigmazz**2 + 2.0_real64*(sigmaxz**2) ) )
            ! print *,'Time step # ',it,' out of ',time_step
            ! print *,'Time: ',(it-1)*DT,' seconds'
            ! print *,'Max vals for vx, vy, vz: ', maxval(vx), maxval(vy), maxval(vz)

            if (velocnorm(it) > stability_threshold) stop 'code became unstable and blew up'

            ! Write the velocity values to an unformatted binary file
            call write_image2(vx, nx, nz, source, it, 'Vx', SINGLE)
            call write_image2(vz, nx, nz, source, it, 'Vz', SINGLE)
            ! Now write the stress Values
            ! call write_image2(sigmaxx, nx, nz, it, 'S1')
            ! call write_image2(sigmayy, nx, nz, it, 'S2')
            ! call write_image2(sigmazz, nx, nz, it, 'S3')
            ! call write_image2(sigmaxy, nx, nz, it, 'S6')
            ! call write_image2(sigmayz, nx, nz, it, 'S4')
            ! call write_image2(sigmaxz, nx, nz, it, 'S5')

        enddo   ! end of time loop
        
        call write_array('velocity_norm.dat', source%time_steps, velocnorm )
        call write_array('stress_norm.dat', source%time_steps, stressnorm )
        
        deallocate(c11, c13, c33, c55 )
        deallocate(rho)
        deallocate(kappa, alpha, acoef, bcoef, kappa_half, alpha_half, acoef_half, bcoef_half)
        deallocate(gamma_x, gamma_z, gamma_xz)
        deallocate(velocnorm, stressnorm)
        deallocate(srcx, srcz)
        deallocate(srcxx, srczz, srcxz)
        deallocate(memory_dvx_dx, memory_dvz_dx, memory_dvx_dz, memory_dvz_dz )
        deallocate(memory_dsigmaxx_dx, memory_dsigmazz_dz, memory_dsigmaxz_dx, memory_dsigmaxz_dz )
        deallocate(vx, vz)
        deallocate(sigmaxx, sigmaxz,  sigmazz)
        
    end subroutine seismic2iso

    ! ==========================================================================    
    ! subroutine seismic2r(domain, source, density_method, SINGLE_OUTPUT)
    !     ! Rotated Staggered Grid to accommodate lower orders of symmetry
    !     use constants
    !     use omp_lib ! Include the OpenMP library 
        
    !     implicit none 
        
    !     ! Input arguments
    !     type(Domain_Type), intent(in) :: domain
    !     type(Source_Type), intent(in) :: source
    !     character(len=256), intent(in) :: density_method
    !     logical, intent(in), optional :: SINGLE_OUTPUT

    !     ! Local variables
    !     real(real64) :: velocnorm, value_dvx_dx, value_dvx_dz, &
    !         value_dvz_dx, value_dvz_dz, value_dsigmaxx_dx, value_dsigmazz_dz, &
    !         value_dsigmaxz_dx, value_dsigmaxz_dz, &
    !         rhoxx, rhozx, rhoxz, rhozz !deltarho,  

    !     ! 1D arrays for damping profiles
    !     real(real64), allocatable :: c11(:,:), c13(:,:), c15(:,:), c33(:,:), c35(:,:), c55(:,:), rho(:,:)

    !     real(real64), allocatable :: kappa(:,:), alpha(:,:), acoef(:,:), bcoef(:,:)
    !     real(real64), allocatable :: kappa_half(:,:), alpha_half(:,:), acoef_half(:,:), bcoef_half(:,:) 
    !     real(real64), allocatable :: gamma_x(:,:), gamma_z(:,:), gamma_xz(:,:)        
    !     real(real64), allocatable :: srcx(:), srcz(:) ! The vector time series of the source
        
    !     ! Model variables
    !     real(real64), allocatable :: vx(:,:),vz(:,:),sigmaxx(:,:),sigmazz(:,:),sigmaxz(:,:)
    !     real(real64), allocatable :: memory_dvx_dx(:,:), memory_dvx_dz(:,:), &
    !                                  memory_dvz_dx(:,:), memory_dvz_dz(:,:), &
    !                                  memory_dsigmaxx_dx(:,:), memory_dsigmazz_dz(:,:), &
    !                                  memory_dsigmaxz_dx(:,:), memory_dsigmaxz_dz(:,:)
        
    !     integer :: nx, nz 
    !     real(real64) :: dx, dz, dt, dl
    !     real(real64) :: fdx_plus, fdx_minus, fdz_plus, fdz_minus
    !     integer :: i, j, it, isource, jsource
    !     logical :: SINGLE
    !     integer :: density_code 

    !     ! The default data output is single precision unless SINGLE_OUTPUT is 
    !     ! set to .FALSE.
    !     if (present(SINGLE_OUTPUT)) then
    !         SINGLE = SINGLE_OUTPUT
    !     else
    !         SINGLE = .TRUE.
    !     end if
        
    !     select case (trim(adjustl(density_method)))
    !         case ('none'      ); density_code = 0
    !         case ('harmonic'  ); density_code = 1
    !         case ('geometric'); density_code = 2
    !         case ('arithmetic'); density_code = 3
    !     end select
        
    !     nx = domain%nx 
    !     nz = domain%nz 
    !     dt = source%dt 
    !     dx = domain%dx 
    !     dz = domain%dz 
        
    !     ! Allocate the arrays based on runtime values of nx and nz
    !     allocate(c11(nx, nz), c13(nx, nz), c15(nx, nz), &
    !              c33(nx, nz), c35(nx, nz), c55(nx, nz), rho(nx, nz))
    !     allocate(kappa(nx,nz), alpha(nx,nz), acoef(nx,nz), bcoef(nx,nz), &
    !              kappa_half(nx,nz), alpha_half(nx,nz), acoef_half(nx,nz), bcoef_half(nx,nz) )
    !     allocate(gamma_x(nx, nz), gamma_z(nx, nz), gamma_xz(nx, nz))
    !     allocate(srcx(source%time_steps), srcz(source%time_steps))
        
    !     ! Allocate more
    !     allocate(memory_dvx_dx(nx, nz), memory_dvx_dz(nx, nz))
    !     allocate(memory_dvz_dx(nx, nz), memory_dvz_dz(nx, nz))
    !     allocate(memory_dsigmaxx_dx(nx, nz), memory_dsigmazz_dz(nx, nz))
    !     allocate(memory_dsigmaxz_dx(nx, nz), memory_dsigmaxz_dz(nx, nz))
    !     allocate(vx(nx, nz), vz(nx, nz))
    !     allocate(sigmaxx(nx, nz), sigmazz(nx, nz), sigmaxz(nx, nz))
        
    !     ! -------------------- Load Stiffness Coefficients --------------------
    
    !     call material_rw2('c11.dat', c11, .TRUE.)
    !     call material_rw2('c13.dat', c13, .TRUE.)
    !     call material_rw2('c15.dat', c15, .TRUE.)
    !     call material_rw2('c33.dat', c33, .TRUE.)
    !     call material_rw2('c35.dat', c35, .TRUE.)
    !     call material_rw2('c55.dat', c55, .TRUE.)
    !     call material_rw2('rho.dat', rho, .TRUE.)
                
    !     ! ------------------- Load Attenuation Coefficients --------------------
    !     call material_rw2('gamma_x.dat', gamma_x, .TRUE.)
    !     call material_rw2('gamma_z.dat', gamma_z, .TRUE.)
    !     call material_rw2('gamma_xz.dat', gamma_xz, .TRUE.)
        
    !     ! ------------------------ Assign some constants -----------------------
        
    !     isource = source%xind + domain%cpml
    !     jsource = source%zind + domain%cpml
    
    !     ! ================================ LOAD SOURCE =========================
    
    !     call loadsource('seismicsourcex.dat', source%time_steps, srcx)
    !     call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        
    !     ! -----------------------------------------------------------------------------
    !     !--- define profile of absorption in PML region

    !     ! Initialize PML 
    !     kappa(:,:) = 1.0_real64
    !     kappa_half(:,:) = 1.0_real64
    !     alpha(:,:) = 0.0_real64
    !     alpha_half(:,:) = 0.0_real64
    !     acoef(:,:) = 0.0_real64
    !     acoef_half(:,:) = 0.0_real64
    !     bcoef(:,:) = 0.0_real64 
    !     bcoef_half(:,:) = 0.0_real64

    !     ! ------------------------------ Load the boundary ----------------------------        
    !     call material_rw2('kappa_cpml.dat', kappa, .TRUE.)
    !     call material_rw2('alpha_cpml.dat', alpha, .TRUE.)
    !     call material_rw2('acoef_cpml.dat', acoef, .TRUE.)
    !     call material_rw2('bcoef_cpml.dat', bcoef, .TRUE.)
        
    !     call material_rw2('kappa_half_cpml.dat', kappa_half, .TRUE.)
    !     call material_rw2('alpha_half_cpml.dat', alpha_half, .TRUE.)
    !     call material_rw2('acoef_half_cpml.dat', acoef_half, .TRUE.)
    !     call material_rw2('bcoef_half_cpml.dat', bcoef_half, .TRUE.)
        
    !     ! =============================================================================
    !     ! Load initial condition
    !     call material_rw2('initialconditionVx.dat', vx, .TRUE.)
    !     call material_rw2('initialconditionVz.dat', vz, .TRUE.)
        
    !     ! vx(:,:) = 0.0_real64
    !     ! vz(:,:) = 0.0_real64
    !     sigmaxx(:,:) = 0.0_real64
    !     sigmazz(:,:) = 0.0_real64
    !     sigmaxz(:,:) = 0.0_real64

    !     ! PML
    !     memory_dvx_dx(:,:) = 0.0_real64
    !     memory_dvx_dz(:,:) = 0.0_real64
    !     memory_dvz_dx(:,:) = 0.0_real64
    !     memory_dvz_dz(:,:) = 0.0_real64
    !     memory_dsigmaxx_dx(:,:) = 0.0_real64
    !     memory_dsigmazz_dz(:,:) = 0.0_real64
    !     memory_dsigmaxz_dx(:,:) = 0.0_real64
    !     memory_dsigmaxz_dz(:,:) = 0.0_real64

    !     !---
    !     !---  beginning of time loop
    !     !---
        
    !     dl = sqrt(dx**2 + dz**2)
        
    !     do it = 1,source%time_steps
    !         !$omp parallel private(i, j, deltarho, value_dvx_dx, value_dvx_dz, value_dvz_dx, value_dvz_dz, value_dsigmaxx_dx, value_dsigmazz_dz, value_dsigmaxz_dx, value_dsigmaxz_dz)
    !         ! ------------------------------------------------------------
    !         !  compute stress sigma and update memory variables for C-PML
    !         ! ------------------------------------------------------------
            
    !         !$omp do
    !         do j = 3,nz-1
    !             do i = 2,nx-2
                    
    !                 fdx_plus = ( 27.0_real64*(vx(i+1,j+1) - vx(i,j)) - (vx(i+2,j+2) - vx(i-1,j-1)) ) / (24.0_real64 * dl)
    !                 fdx_minus = ( 27.0_real64*(vx(i+1,j-1) - vx(i,j)) - (vx(i+2,j-2) - vx(i-1,j+1)) ) / (24.0_real64 * dl)
    !                 value_dvx_dx = (dx/dl) * 0.50_real64 * (fdx_plus + fdx_minus)
    !                 value_dvx_dz = (dz/dl) * 0.50_real64 * (fdx_plus - fdx_minus)

    !                 fdz_plus = ( 27.0_real64*(vz(i+1,j+1) - vz(i,j)) - (vz(i+2,j+2) - vz(i-1,j-1)) ) / (24.0_real64 * dl)
    !                 fdz_minus = ( 27.0_real64*(vz(i+1,j-1) - vz(i,j)) - (vz(i+2,j-2) - vz(i-1,j+1)) ) / (24.0_real64 * dl)                    
    !                 value_dvz_dx = (dx/dl) * 0.50_real64 * (fdz_plus + fdz_minus)
    !                 value_dvz_dz = (dz/dl) * 0.50_real64 * (fdz_plus - fdz_minus)
                    
    !                 memory_dvx_dx(i,j) = bcoef(i,j) * memory_dvx_dx(i,j) + acoef(i,j) * value_dvx_dx
    !                 memory_dvx_dz(i,j) = bcoef(i,j) * memory_dvx_dz(i,j) + acoef(i,j) * value_dvx_dz 
    !                 memory_dvz_dz(i,j) = bcoef(i,j) * memory_dvz_dz(i,j) + acoef(i,j) * value_dvz_dz
    !                 memory_dvz_dx(i,j) = bcoef(i,j) * memory_dvz_dx(i,j) + acoef(i,j) * value_dvz_dx
                    
    !                 value_dvx_dx = value_dvx_dx / kappa(i,j) + memory_dvx_dx(i,j)
    !                 value_dvx_dz = value_dvx_dz / kappa(i,j) + memory_dvx_dz(i,j)
    !                 value_dvz_dz = value_dvz_dz / kappa(i,j) + memory_dvz_dz(i,j)
    !                 value_dvz_dx = value_dvz_dx / kappa(i,j) + memory_dvz_dx(i,j)
                    
    !                 sigmaxx(i,j) = ( sigmaxx(i,j) + &
    !                     (   c11(i,j) * value_dvx_dx + &
    !                         c13(i,j) * value_dvz_dz + &
    !                         c15(i,j) * (value_dvz_dx + value_dvx_dz) ) * dt ) / & 
    !                            (1 + gamma_x(i,j) * dt )
    !                 sigmazz(i,j) = ( sigmazz(i,j) + &
    !                     (   c13(i,j) * value_dvx_dx + &
    !                         c33(i,j) * value_dvz_dz + &
    !                         c35(i,j) * (value_dvz_dx + value_dvx_dz) ) * dt) / & 
    !                             (1 + gamma_z(i,j) * dt )
    !             enddo
    !         enddo
    !         !$omp end do
            
    !         !$omp do
    !         do j = 2,nz-2
    !             do i = 3,nx-1
                    
    !                 ! dvx/dz: use (x+z) ⇄ (x−z) diagonals perpendicular to normal-stress direction
    !                 fdx_plus  = ( 27.0_real64*(vx(i+1,j-1) - vx(i,j)) - (vx(i+2,j-2) - vx(i-1,j+1)) ) / (24.0_real64*dl)
    !                 fdx_minus = ( 27.0_real64*(vx(i+1,j+1) - vx(i,j)) - (vx(i+2,j+2) - vx(i-1,j-1)) ) / (24.0_real64*dl)
    !                 value_dvx_dx = (dx/dl)*0.50_real64*(fdx_plus + fdx_minus)
    !                 value_dvx_dz = (dz/dl)*0.50_real64*(fdx_plus - fdx_minus)

    !                 ! dvz/dx: mirrored diagonals
    !                 fdz_plus  = ( 27.0_real64*(vz(i+1,j+1) - vz(i,j)) - (vz(i+2,j+2) - vz(i-1,j-1)) ) / (24.0_real64*dl)
    !                 fdz_minus = ( 27.0_real64*(vz(i+1,j-1) - vz(i,j)) - (vz(i+2,j-2) - vz(i-1,j+1)) ) / (24.0_real64*dl)
    !                 value_dvz_dx = (dx/dl)*0.50_real64*(fdz_plus + fdz_minus)
    !                 value_dvz_dz = (dz/dl)*0.50_real64*(fdz_plus - fdz_minus)

                    
    !                 memory_dvx_dx(i,j) = bcoef(i,j) * memory_dvx_dx(i,j) + acoef(i,j) * value_dvx_dx
    !                 memory_dvx_dz(i,j) = bcoef(i,j) * memory_dvx_dz(i,j) + acoef(i,j) * value_dvx_dz 
    !                 memory_dvz_dz(i,j) = bcoef(i,j) * memory_dvz_dz(i,j) + acoef(i,j) * value_dvz_dz
    !                 memory_dvz_dx(i,j) = bcoef(i,j) * memory_dvz_dx(i,j) + acoef(i,j) * value_dvz_dx
                    
    !                 value_dvx_dx = value_dvx_dx / kappa(i,j) + memory_dvx_dx(i,j)
    !                 value_dvx_dz = value_dvx_dz / kappa(i,j) + memory_dvx_dz(i,j)
    !                 value_dvz_dz = value_dvz_dz / kappa(i,j) + memory_dvz_dz(i,j)
    !                 value_dvz_dx = value_dvz_dx / kappa(i,j) + memory_dvz_dx(i,j)
                    
    !                 sigmaxz(i,j) = ( sigmaxz(i,j) + &
    !                     (   c15(i,j)  * value_dvx_dx + & 
    !                         c35(i,j)  * value_dvz_dz + &
    !                         c55(i,j) * (value_dvz_dx + value_dvx_dz) ) * dt ) / &
    !                             (1 + gamma_xz(i,j) * dt )
    !             enddo
    !         enddo
    !         !$omp end do
            
    !         ! --------------------------------------------------------
    !         !  compute velocity and update memory variables for C-PML
    !         ! --------------------------------------------------------
            
    !         !$omp do 
    !         do j = 3, nz-1
    !             do i = 3, nx-1

    !                 rhoxx = face_density(rho(i,j), rho(i-1,j), density_code)
    !                 rhozx = face_density(rho(i,j), rho(i,j-1), density_code)

    !                 ! ==== Rotated finite‑difference derivatives for σxx ====
    !                 fdx_plus  = ( 27.0_real64*(sigmaxx(i+1,j+1) - sigmaxx(i,j)) - (sigmaxx(i+2,j+2) - sigmaxx(i-1,j-1)) ) / (24.0_real64 * dl)
    !                 fdx_minus = ( 27.0_real64*(sigmaxx(i+1,j-1) - sigmaxx(i,j)) - (sigmaxx(i+2,j-2) - sigmaxx(i-1,j+1)) ) / (24.0_real64 * dl)
    !                 value_dsigmaxx_dx = (dx/dl) * 0.50_real64 * (fdx_plus + fdx_minus)
    !                 value_dsigmaxx_dz = (dz/dl) * 0.50_real64 * (fdx_plus - fdx_minus)

    !                 ! ==== Rotated finite‑difference derivatives for σxz ====
    !                 fdz_plus  = ( 27.0_real64*(sigmaxz(i+1,j+1) - sigmaxz(i,j)) - (sigmaxz(i+2,j+2) - sigmaxz(i-1,j-1)) ) / (24.0_real64 * dl)
    !                 fdz_minus = ( 27.0_real64*(sigmaxz(i+1,j-1) - sigmaxz(i,j)) - (sigmaxz(i+2,j-2) - sigmaxz(i-1,j+1)) ) / (24.0_real64 * dl)
    !                 value_dsigmaxz_dx = (dx/dl) * 0.50_real64 * (fdz_plus + fdz_minus)
    !                 value_dsigmaxz_dz = (dz/dl) * 0.50_real64 * (fdz_plus - fdz_minus)

    !                 ! ==== Apply CPML memory updates ====
    !                 memory_dsigmaxx_dx(i,j) = bcoef_half(i,j)*memory_dsigmaxx_dx(i,j) + acoef_half(i,j)*value_dsigmaxx_dx
    !                 memory_dsigmaxz_dz(i,j) = bcoef_half(i,j)*memory_dsigmaxz_dz(i,j) + acoef_half(i,j)*value_dsigmaxz_dz

    !                 value_dsigmaxx_dx = value_dsigmaxx_dx / kappa_half(i,j) + memory_dsigmaxx_dx(i,j)
    !                 value_dsigmaxz_dz = value_dsigmaxz_dz / kappa_half(i,j) + memory_dsigmaxz_dz(i,j)

    !                 ! ==== Velocity update ====
    !                 vx(i,j) = vx(i,j) + dt * ( value_dsigmaxx_dx / rhoxx + value_dsigmaxz_dz / rhozx )
    !             end do
    !         end do
    !         !$omp end do

    !         !$omp do
    !         do j = 2, nz-2
    !             do i = 2, nx-2

    !                 rhoxz = face_density(rho(i,j), rho(i+1,j), density_code)
    !                 rhozz = face_density(rho(i,j), rho(i,j+1), density_code)

    !                 ! ===== Rotated derivative of σxz (along x+z and x−z diagonals) =====
    !                 fdx_plus  = ( 27.0_real64*(sigmaxz(i+1,j+1) - sigmaxz(i,j)) - (sigmaxz(i+2,j+2) - sigmaxz(i-1,j-1)) ) / (24.0_real64 * dl)
    !                 fdx_minus = ( 27.0_real64*(sigmaxz(i+1,j-1) - sigmaxz(i,j)) - (sigmaxz(i+2,j-2) - sigmaxz(i-1,j+1)) ) / (24.0_real64 * dl)
    !                 value_dsigmaxz_dx = (dx/dl) * 0.50_real64 * (fdx_plus + fdx_minus)
    !                 value_dsigmaxz_dz = (dz/dl) * 0.50_real64 * (fdx_plus - fdx_minus)

    !                 ! ===== Rotated derivative of σzz (also along diagonals) =====
    !                 fdz_plus  = ( 27.0_real64*(sigmazz(i+1,j+1) - sigmazz(i,j)) - (sigmazz(i+2,j+2) - sigmazz(i-1,j-1)) ) / (24.0_real64 * dl)
    !                 fdz_minus = ( 27.0_real64*(sigmazz(i+1,j-1) - sigmazz(i,j)) - (sigmazz(i+2,j-2) - sigmazz(i-1,j+1)) ) / (24.0_real64 * dl)
    !                 value_dsigmazz_dx = (dx/dl) * 0.50_real64 * (fdz_plus + fdz_minus)
    !                 value_dsigmazz_dz = (dz/dl) * 0.50_real64 * (fdz_plus - fdz_minus)

    !                 ! ===== Apply CPML memory update and projection =====
    !                 memory_dsigmaxz_dx(i,j) = bcoef_half(i,j)*memory_dsigmaxz_dx(i,j) + acoef_half(i,j)*value_dsigmaxz_dx
    !                 memory_dsigmazz_dz(i,j) = bcoef_half(i,j)*memory_dsigmazz_dz(i,j) + acoef_half(i,j)*value_dsigmazz_dz

    !                 value_dsigmaxz_dx = value_dsigmaxz_dx / kappa_half(i,j) + memory_dsigmaxz_dx(i,j)
    !                 value_dsigmazz_dz = value_dsigmazz_dz / kappa_half(i,j) + memory_dsigmazz_dz(i,j)

    !                 ! ===== Vz velocity update =====
    !                 vz(i,j) = vz(i,j) + dt * (value_dsigmaxz_dx / rhoxz + value_dsigmazz_dz / rhozz)
    !             end do
    !         end do
    !         !$omp end do

    !         !$omp end parallel
            
    !         ! Add the source term
    !         ! vx(isource,jsource) = vx(isource,jsource) + srcx(it) * dt / rho(isource,jsource)
    !         ! vz(isource,jsource) = vz(isource,jsource) + srcz(it) * dt / rho(isource,jsource)
    !         vx(isource,jsource) = vx(isource,jsource) + dt * srcx(it) / rho(isource,jsource)
    !         vz(isource,jsource) = vz(isource,jsource) + dt * srcz(it) / rho(isource,jsource)
        
    !         ! Dirichlet conditions (rigid boundaries) on the edges or at the 
    !         ! bottom of the PML layers
    !         vx(1,:) = 0.0_real64
    !         vx(nx,:) = 0.0_real64
        
    !         vx(:,1) = 0.0_real64
    !         vx(:,nz) = 0.0_real64
        
    !         vz(1,:) = 0.0_real64
    !         vz(nx,:) = 0.0_real64
        
    !         vz(:,1) = 0.0_real64
    !         vz(:,nz) = 0.0_real64
        
    !         ! print maximum of norm of velocity
    !         velocnorm = maxval(sqrt(vx**2 + vz**2))
    !         if (velocnorm > stability_threshold) stop 'code became unstable and blew up'
        
    !         call write_image2(vx, nx, nz, source, it, 'Vx', SINGLE)
    !         call write_image2(vz, nx, nz, source, it, 'Vz', SINGLE)

    !     enddo   ! end of time loop
        
    !     deallocate(c11, c13, c15, c33, c35, c55, rho)
    !     deallocate(kappa, alpha, acoef, bcoef, kappa_half, alpha_half, acoef_half, bcoef_half)
    !     deallocate(gamma_x, gamma_z, gamma_xz, srcx, srcz)
        
    !     ! Allocate more
    !     deallocate(memory_dvx_dx, memory_dvx_dz)
    !     deallocate(memory_dvz_dx, memory_dvz_dz)
    !     deallocate(memory_dsigmaxx_dx, memory_dsigmazz_dz)
    !     deallocate(memory_dsigmaxz_dx, memory_dsigmaxz_dz)
    !     deallocate(vx, vz, sigmaxx, sigmazz, sigmaxz)
        
    ! end subroutine seismic2r
    ! =========================================================================
    subroutine seismic25iso(domain, source, density_method, SINGLE_OUTPUT)
        !--------------------------------------------------------------------------------------
        use constants
        
        implicit none

        ! Input arguments
        type(Domain_Type), intent(in) :: domain 
        type(Source_Type), intent(in) :: source 
        logical, intent(in), optional :: SINGLE_OUTPUT
        character(len=256), intent(in) :: density_method
        ! Local variables
        real(real64) :: rhoxx, rhoyx, rhozx, &
                        rhoxy, rhoyy, rhozy, &
                        rhoxz, rhoyz, rhozz
        
        real(real64), allocatable :: velocnorm(:), stressnorm(:)

        real(real64), allocatable :: c11(:,:), c12(:,:), c13(:,:), &
                                     c22(:,:), c23(:,:), c33(:,:), &
                                     c44(:,:), c55(:,:), c66(:,:), rho(:,:)
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
        real(real64), allocatable :: kappa(:,:), alpha(:,:), acoef(:,:), bcoef(:,:), &
                        kappa_half(:,:), alpha_half(:,:), acoef_half(:,:), bcoef_half(:,:)
        
        ! Arrays for the PML damping factors
        real(real64), allocatable :: gamma_x(:,:), gamma_y(:,:), gamma_z(:,:)
        real(real64), allocatable :: gamma_xz(:,:), gamma_xy(:,:), gamma_yz(:,:)

        ! Source arrays
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:)
        real(real64), allocatable :: srcxx(:), srcxy(:), srcxz(:), srcyy(:), srcyz(:), srczz(:)
        
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
        integer :: density_code
        
        ! The default data output is single precision unless SINGLE_OUTPUT is 
        ! set to .FALSE.
        if (present(SINGLE_OUTPUT)) then 
            SINGLE = SINGLE_OUTPUT 
        else
            SINGLE = .TRUE.
        endif
        
        select case (trim(adjustl(density_method)))
            case ('none'      ); density_code = 0
            case ('harmonic'  ); density_code = 1
            case ('geometric'); density_code = 2
            case ('arithmetic'); density_code = 3
        end select
        
        nx = domain%nx
        ny = domain%ny 
        nz = domain%nz 
        dt = source%dt 
        dx = domain%dx
        dy = domain%dy 
        dz = domain%dz
        
        ! ----------------------- Allocate Arrays ----------------------------
        allocate(c11(nx, nz), c12(nx, nz), c13(nx, nz), &
                 c22(nx, nz), c23(nx, nz), c33(nx, nz), &
                 c44(nx, nz), c55(nx, nz), c66(nx, nz) )
        allocate(rho(nx,nz))
                
        allocate(kappa(nx, nz), alpha(nx, nz), acoef(nx, nz), bcoef(nx, nz), &
                kappa_half(nx, nz), alpha_half(nx, nz), acoef_half(nx, nz), bcoef_half(nx, nz))
        allocate(gamma_x(nx, nz), gamma_y(nx, nz), gamma_z(nx, nz))
        allocate(gamma_xy(nx, nz), gamma_yz(nx, nz), gamma_xz(nx, nz))
        
        allocate(velocnorm(source%time_steps), stressnorm(source%time_steps))
        
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps))
        allocate(srcxx(source%time_steps), srcxy(source%time_steps), srcxz(source%time_steps), &
                 srcyy(source%time_steps), srcyz(source%time_steps), srczz(source%time_steps) )                
        
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
            call material_rw2('c22.dat', c22, .TRUE.)
            call material_rw2('c23.dat', c23, .TRUE.)
            call material_rw2('c33.dat', c33, .TRUE.)
            call material_rw2('c44.dat', c44, .TRUE.)
            call material_rw2('c55.dat', c55, .TRUE.)
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
        srcxx(:) = 0.0_real64
        srcxy(:) = 0.0_real64
        srcxz(:) = 0.0_real64
        srcyy(:) = 0.0_real64
        srcyz(:) = 0.0_real64
        srczz(:) = 0.0_real64
        srcx(:) = 0.0_real64
        srcy(:) = 0.0_real64
        srcz(:) = 0.0_real64
        
        if ( source%source_type == 'ac') then 
            call loadsource('seismicsourcex.dat', source%time_steps, srcx)
            call loadsource('seismicsourcey.dat', source%time_steps, srcy)
            call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        else
            call loadsource('seismicsourcexx.dat', source%time_steps, srcxx)
            call loadsource('seismicsourcexy.dat', source%time_steps, srcxy)
            call loadsource('seismicsourcexz.dat', source%time_steps, srcxz)
            call loadsource('seismicsourceyy.dat', source%time_steps, srcyy)
            call loadsource('seismicsourceyz.dat', source%time_steps, srcyz)
            call loadsource('seismicsourcezz.dat', source%time_steps, srczz)
        endif 
        
        ! ==================================== PML ====================================
        ! Initialize PML 
        kappa(:,:) = 1.0_real64
        kappa_half(:,:) = 1.0_real64
        alpha(:,:) = 0.0_real64
        alpha_half(:,:) = 0.0_real64
        acoef(:,:) = 0.0_real64
        acoef_half(:,:) = 0.0_real64

        ! ------------------------- Boundary Conditions -------------------------
        call material_rw2('kappa_cpml.dat', kappa, .TRUE.)
        call material_rw2('alpha_cpml.dat', alpha, .TRUE.)
        call material_rw2('acoef_cpml.dat', acoef, .TRUE.)
        call material_rw2('bcoef_cpml.dat', bcoef, .TRUE.)

        call material_rw2('kappa_half_cpml.dat', kappa_half, .TRUE.)
        call material_rw2('alpha_half_cpml.dat', alpha_half, .TRUE.)
        call material_rw2('acoef_half_cpml.dat', acoef_half, .TRUE.)
        call material_rw2('bcoef_half_cpml.dat', bcoef_half, .TRUE.)

        ! Load initial condition
        call material_rw3('initialconditionVx.dat', vx, .TRUE.)
        call material_rw3('initialconditionVy.dat', vy, .TRUE.)
        call material_rw3('initialconditionVz.dat', vz, .TRUE.)

        ! Initialize the stress values
        sigmaxx(:,:,:) = 0.0_real64
        sigmayy(:,:,:) = 0.0_real64
        sigmazz(:,:,:) = 0.0_real64
        sigmaxy(:,:,:) = 0.0_real64
        sigmaxz(:,:,:) = 0.0_real64
        sigmayz(:,:,:) = 0.0_real64

        ! PML
        memory_dvx_dx(:,:,:) = 0.0_real64
        memory_dvx_dy(:,:,:) = 0.0_real64
        memory_dvx_dz(:,:,:) = 0.0_real64
        memory_dvy_dx(:,:,:) = 0.0_real64
        memory_dvy_dy(:,:,:) = 0.0_real64
        memory_dvy_dz(:,:,:) = 0.0_real64
        memory_dvz_dx(:,:,:) = 0.0_real64
        memory_dvz_dy(:,:,:) = 0.0_real64 
        memory_dvz_dz(:,:,:) = 0.0_real64

        memory_dsigmaxx_dx(:,:,:) = 0.0_real64
        memory_dsigmayy_dy(:,:,:) = 0.0_real64
        memory_dsigmazz_dz(:,:,:) = 0.0_real64
        memory_dsigmaxy_dx(:,:,:) = 0.0_real64
        memory_dsigmaxy_dy(:,:,:) = 0.0_real64
        memory_dsigmaxz_dx(:,:,:) = 0.0_real64
        memory_dsigmaxz_dz(:,:,:) = 0.0_real64
        memory_dsigmayz_dy(:,:,:) = 0.0_real64
        memory_dsigmayz_dz(:,:,:) = 0.0_real64
        
        ! Initialize velocnorm and stressnorm 
        velocnorm(:) = 0.0_real64
        stressnorm(:) = 0.0_real64 
        ! Do it 
        
        ! =============================== Forward Model ===============================
        do it = 1,source%time_steps
            !------------------------------------------------------------
            ! compute stress sigma and update memory variables for C-PML
            !------------------------------------------------------------
            ! Update in the x direction
            ! Loop 1
            do k = 3,nz-1
                do j = 3,ny-1
                    do i = 2,nx-2
                        ! Forward step on half grid
                        dvx_dx = ( 27.0_real64*(vx(i+1,j,k) - vx(i,j,k)) - (vx(i+2,j,k) - vx(i-1,j,k))) / (24.0_real64*dx)
                        ! Backward step on full grid
                        dvy_dy = ( 27.0_real64*(vy(i,j,k) - vy(i,j-1,k)) - (vy(i,j+1,k) - vy(i,j-2,k))) / (24.0_real64*dy)
                        ! Backward step on full grid
                        dvz_dz = ( 27.0_real64*(vz(i,j,k) - vz(i,j,k-1)) - (vz(i,j,k+1) - vz(i,j,k-2))) / (24.0_real64*dz)
                        
                        memory_dvx_dx(i,j,k) = bcoef_half(i,k) * memory_dvx_dx(i,j,k) + acoef_half(i,k) * dvx_dx
                        memory_dvy_dy(i,j,k) = bcoef(i,k) * memory_dvy_dy(i,j,k) + acoef(i,k) * dvy_dy
                        memory_dvz_dz(i,j,k) = bcoef(i,k) * memory_dvz_dz(i,j,k) + acoef(i,k) * dvz_dz
                        
                        dvx_dx = dvx_dx / kappa_half(i,k) + memory_dvx_dx(i,j,k)
                        dvy_dy = dvy_dy / kappa(i,k) + memory_dvy_dy(i,j,k)
                        dvz_dz = dvz_dz / kappa(i,k) + memory_dvz_dz(i,j,k)
                        
                        sigmaxx(i,j,k) = ( sigmaxx(i,j,k) + &
                        ( c11(i,k) * dvx_dx + c12(i,k) * dvy_dy + c13(i,k) * dvz_dz ) * dt ) / &
                                (1 + gamma_x(i,k) * dt )

                        ! Full 3D will need a gradient in the y-direction
                        sigmayy(i,j,k) = ( sigmayy(i,j,k) + &
                        ( c12(i,k) * dvx_dx + c22(i,k) * dvy_dy + c23(i,k) * dvz_dz ) * dt ) / &
                                (1 + gamma_y(i,k) * dt )

                        sigmazz(i,j,k) = ( sigmazz(i,j,k) + &
                        ( c13(i,k) * dvx_dx + c23(i,k) * dvy_dy + c33(i,k) * dvz_dz ) * dt ) / &
                                (1 + gamma_z(i,k) * dt )

                    enddo
                enddo
            enddo
            ! Loop 2
            ! Update sigmaxy, x-direction is full nodes
            do k = 3,nz-1
                do j = 2,ny-2
                    do i = 3,nx-1
                        ! Backward step full grid
                        dvy_dx = ( 27.0_real64*(vy(i,j,k) - vy(i-1,j,k)) - (vy(i+1,j,k) - vy(i-2,j,k))) / (24.0_real64*dx)
                        ! Forward step half grid
                        dvx_dy = ( 27.0_real64*(vx(i,j+1,k) - vx(i,j,k)) - (vx(i,j+2,k) - vx(i,j-1,k))) / (24.0_real64*dy)
                        
                        memory_dvy_dx(i,j,k) = bcoef(i,k) * memory_dvy_dx(i,j,k) + acoef(i,k) * dvy_dx
                        memory_dvx_dy(i,j,k) = bcoef_half(i,k) * memory_dvx_dy(i,j,k) + acoef_half(i,k) * dvx_dy
                        
                        dvy_dx = dvy_dx / kappa(i,k) + memory_dvy_dx(i,j,k)
                        dvx_dy = dvx_dy / kappa_half(i,k) + memory_dvx_dy(i,j,k)
                        
                        sigmaxy(i,j,k) = ( sigmaxy(i,j,k) + &
                        (  c66(i,k) * (dvy_dx + dvx_dy) ) * dt ) / &
                                (1 + gamma_xy(i,k) * dt )

                    enddo
                enddo
            enddo
            ! Loop 3
            ! Update sigmaxz, z-direction is full nodes
            do k = 2,nz-2
                do j = 3,ny-1
                    do i = 3,nx-1
                        ! Backward difference full grid
                        dvz_dx = ( 27.0_real64*(vz(i,j,k) - vz(i-1,j,k)) - (vz(i+1,j,k) - vz(i-2,j,k))) / (24.0_real64*dx)
                        ! Forward difference half grid
                        dvx_dz = ( 27.0_real64*(vx(i,j,k+1) - vx(i,j,k)) - (vx(i,j,k+2) - vx(i,j,k-1))) / (24.0_real64*dz)
                        
                        memory_dvz_dx(i,j,k) = bcoef(i,k) * memory_dvz_dx(i,j,k) + acoef(i,k) * dvz_dx
                        memory_dvx_dz(i,j,k) = bcoef_half(i,k) * memory_dvx_dz(i,j,k) + acoef_half(i,k) * dvx_dz

                        dvz_dx = dvz_dx / kappa(i,k) + memory_dvz_dx(i,j,k) 
                        dvx_dz = dvx_dz / kappa_half(i,k) + memory_dvx_dz(i,j,k)
                        
                        sigmaxz(i,j,k) = ( sigmaxz(i,j,k) + &
                            (  c55(i,k) * ( dvx_dz + dvz_dx) ) * dt  ) / &
                                    (1 + gamma_xz(i,k) * dt )

                    enddo
                enddo
                ! Loop 4
                !   ! update sigmayz, y-direction is full nodes
                do j = 2,ny-2
                    do i = 2,nx-2
                        ! Forward difference half grid
                        dvz_dy = ( 27.0_real64*(vz(i,j+1,k) - vz(i,j,k)) - (vz(i,j+2,k) - vz(i,j-1,k))) / (24.0_real64*dy)
                        ! Forwar difference half grid
                        dvy_dz = ( 27.0_real64*(vy(i,j,k+1) - vy(i,j,k)) - (vy(i,j,k+2) - vy(i,j,k-1))) / (24.0_real64*dz)
                        
                        memory_dvz_dy(i,j,k) = bcoef_half(i,k) * memory_dvz_dy(i,j,k) + acoef_half(i,k) * dvz_dy
                        memory_dvy_dz(i,j,k) = bcoef_half(i,k) * memory_dvy_dz(i,j,k) + acoef_half(i,k) * dvy_dz

                        dvz_dy = dvz_dy / kappa_half(i,k) + memory_dvz_dy(i,j,k)
                        dvy_dz = dvy_dz / kappa_half(i,k) + memory_dvy_dz(i,j,k)

                        sigmayz(i,j,k) = ( sigmayz(i,j,k)  + &
                            ( c44(i,k) * ( dvy_dz + dvz_dy) ) * dt  ) / & 
                                    (1 + gamma_yz(i,k) * dt )
                    enddo
                enddo
            enddo

            !--------------------------------------------------------
            ! compute velocity and update memory variables for C-PML
            !--------------------------------------------------------
            ! Loop 5
            do k = 3,nz-1
                do j = 3,ny-1
                    do i = 3,nx-1
                        ! ds1/dx, ds6/dy, ds5,dz
                        rhoxx = face_density(rho(i,k), rho(i-1,k), density_code)
                        rhoyx = face_density(rho(i,k), rho(i,k), density_code) 
                        rhozx = face_density(rho(i,k), rho(i,k-1), density_code) 
                        ! Backward difference half grid
                        dsigmaxx_dx = (27.0_real64*sigmaxx(i,j,k) - 27.0_real64*sigmaxx(i-1,j,k) + sigmaxx(i-2,j,k)) / (24.0_real64*dx)
                        dsigmaxy_dy = (27.0_real64*sigmaxy(i,j,k) - 27.0_real64*sigmaxy(i,j-1,k) + sigmaxy(i,j-2,k)) / (24.0_real64*dy)
                        dsigmaxz_dz = (27.0_real64*sigmaxz(i,j,k) - 27.0_real64*sigmaxz(i,j,k-1) + sigmaxz(i,j,k-2)) / (24.0_real64*dz)
                        
                        memory_dsigmaxx_dx(i,j,k) = bcoef_half(i,k) * memory_dsigmaxx_dx(i,j,k) + acoef_half(i,k) * dsigmaxx_dx
                        memory_dsigmaxy_dy(i,j,k) = bcoef_half(i,k) * memory_dsigmaxy_dy(i,j,k) + acoef_half(i,k) * dsigmaxy_dy
                        memory_dsigmaxz_dz(i,j,k) = bcoef_half(i,k) * memory_dsigmaxz_dz(i,j,k) + acoef_half(i,k) * dsigmaxz_dz

                        dsigmaxx_dx = dsigmaxx_dx / kappa_half(i,k) + memory_dsigmaxx_dx(i,j,k)
                        dsigmaxy_dy = dsigmaxy_dy / kappa_half(i,k) + memory_dsigmaxy_dy(i,j,k)
                        dsigmaxz_dz = dsigmaxz_dz / kappa_half(i,k) + memory_dsigmaxz_dz(i,j,k) 

                        vx(i,j,k) = vx(i,j,k) + &
                            (dsigmaxx_dx/rhoxx + dsigmaxy_dy/rhoyx + dsigmaxz_dz/rhozx) * dt 
                    enddo
                enddo

                do j = 2,ny-2
                    do i = 2,nx-2
                        ! ds6/dx, ds2/dy, ds4/dz
                        rhoxy = face_density(rho(i,k), rho(i+1,k), density_code)
                        rhoyy = face_density(rho(i,k), rho(i,k), density_code) 
                        rhozy = face_density(rho(i,k), rho(i,k-1), density_code)
                        ! Forward difference half grid
                        dsigmaxy_dx = (-27.0_real64*sigmaxy(i,j,k) + 27.0_real64*sigmaxy(i+1,j,k) - sigmaxy(i+2,j,k)) / (24.0_real64*dx)
                        dsigmayy_dy = (-27.0_real64*sigmayy(i,j,k) + 27.0_real64*sigmayy(i,j+1,k) - sigmayy(i,j+2,k)) / (24.0_real64*dy)
                        ! Backward difference full grid
                        dsigmayz_dz = (27.0_real64*sigmayz(i,j,k) - 27.0_real64*sigmayz(i,j,k-1) + sigmayz(i,j,k-2)) / (24.0_real64*dz)
                        
                        memory_dsigmaxy_dx(i,j,k) = bcoef_half(i,k) * memory_dsigmaxy_dx(i,j,k) + acoef_half(i,k) * dsigmaxy_dx
                        memory_dsigmayy_dy(i,j,k) = bcoef_half(i,k) * memory_dsigmayy_dy(i,j,k) + acoef_half(i,k) * dsigmayy_dy
                        memory_dsigmayz_dz(i,j,k) = bcoef(i,k) * memory_dsigmayz_dz(i,j,k) + acoef(i,k) * dsigmayz_dz
                        
                        dsigmaxy_dx = dsigmaxy_dx / kappa_half(i,k) + memory_dsigmaxy_dx(i,j,k)
                        dsigmayy_dy = dsigmayy_dy / kappa_half(i,k) + memory_dsigmayy_dy(i,j,k)
                        dsigmayz_dz = dsigmayz_dz / kappa(i,k) + memory_dsigmayz_dz(i,j,k)

                        vy(i,j,k) = vy(i,j,k) + &
                            (dsigmaxy_dx/rhoxy + dsigmayy_dy/rhoyy + dsigmayz_dz/rhozy) * dt
                    enddo
                enddo
            enddo

            do k = 2,nz-2
                do j = 3,ny-1
                    do i = 2,nx-2
                        ! ds5/dx, ds4/dy, ds3/dz
                        rhoxz = face_density(rho(i,k), rho(i+1,k), density_code)
                        rhoyz = face_density(rho(i,k), rho(i,k), density_code) 
                        rhozz = face_density(rho(i,k), rho(i,k+1), density_code)
                        
                        ! Forward difference half grid
                        dsigmaxz_dx = (-27.0_real64*sigmaxz(i,j,k) + 27.0_real64*sigmaxz(i+1,j,k) - sigmaxz(i+2,j,k)) / (24.0_real64*dx)
                        ! Backward difference full grid
                        dsigmayz_dy = (27.0_real64*sigmayz(i,j,k) - 27.0_real64*sigmayz(i,j-1,k) + sigmayz(i,j-2,k)) / (24.0_real64*dy)
                        ! Forward difference half grid
                        dsigmazz_dz = (-27.0_real64*sigmazz(i,j,k) + 27.0_real64*sigmazz(i,j,k+1) - sigmazz(i,j,k+2)) / (24.0_real64*dz)
                        
                        memory_dsigmaxz_dx(i,j,k) = bcoef_half(i,k) * memory_dsigmaxz_dx(i,j,k) + acoef_half(i,k) * dsigmaxz_dx
                        memory_dsigmayz_dy(i,j,k) = bcoef(i,k) * memory_dsigmayz_dy(i,j,k) + acoef(i,k) * dsigmayz_dy
                        memory_dsigmazz_dz(i,j,k) = bcoef_half(i,k) * memory_dsigmazz_dz(i,j,k) + acoef_half(i,k) * dsigmazz_dz

                        dsigmaxz_dx = dsigmaxz_dx / kappa_half(i,k) + memory_dsigmaxz_dx(i,j,k)
                        dsigmayz_dy = dsigmayz_dy / kappa(i,k) + memory_dsigmayz_dy(i,j,k)
                        dsigmazz_dz = dsigmazz_dz / kappa_half(i,k) + memory_dsigmazz_dz(i,j,k)

                        vz(i,j,k) = vz(i,j,k) + &
                            (dsigmaxz_dx/rhoxz + dsigmayz_dy/rhoyz + dsigmazz_dz/rhozz) * &
                            dt !/ deltarho !rho(i,k)

                    enddo
                enddo
            enddo
            
            sigmaxx(isource,jsource,ksource) = sigmaxx(isource,jsource,ksource) + srcxx(it) / rho(isource,ksource)  
            sigmaxy(isource+1,jsource+1,ksource) = sigmaxy(isource+1,jsource+1,ksource) + srcxy(it) / rho(isource+1,ksource+1)
            sigmaxz(isource+1,jsource,ksource+1) = sigmaxz(isource+1,jsource,ksource+1) + srcxz(it) / rho(isource+1,ksource)  
            sigmayy(isource,jsource,ksource) = sigmayy(isource,jsource,ksource) + srcyy(it) / rho(isource,ksource)
            sigmayz(isource,jsource+1,ksource+1) = sigmayz(isource,jsource+1,ksource+1) + srcyz(it) / rho(isource,ksource+1)  
            sigmazz(isource,jsource,ksource) = sigmazz(isource,jsource,ksource) + srczz(it) / rho(isource,ksource)
            vx(isource+1,jsource,ksource) = vx(isource+1,jsource,ksource) + &
                    srcx(it) * dt / rho(isource+1,ksource)
            vy(isource,jsource+1,ksource) = vy(isource,jsource+1,ksource) + &
                    srcy(it) * dt / rho(isource,ksource)
            vz(isource,jsource,ksource+1) = vz(isource,jsource,ksource+1) + &
                    srcz(it) * dt / rho(isource,ksource+1)

            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            vx(1,:,:) = 0.0_real64
            vy(1,:,:) = 0.0_real64
            vz(1,:,:) = 0.0_real64

            vx(:,1,:) = 0.0_real64
            vy(:,1,:) = 0.0_real64
            vz(:,1,:) = 0.0_real64

            vx(:,:,1) = 0.0_real64
            vy(:,:,1) = 0.0_real64
            vz(:,:,1) = 0.0_real64

            vx(nx,:,:) = 0.0_real64
            vy(nx,:,:) = 0.0_real64
            vz(nx,:,:) = 0.0_real64

            vx(:,ny,:) = 0.0_real64
            vy(:,ny,:) = 0.0_real64
            vz(:,ny,:) = 0.0_real64

            vx(:,:,nz) = 0.0_real64
            vy(:,:,nz) = 0.0_real64
            vz(:,:,nz) = 0.0_real64

            ! check norm of velocity to make sure the solution isn't diverging
            velocnorm(it) = maxval( sqrt(vx**2 + vy**2 + vz**2) )
            stressnorm(it) = maxval( &
                        sqrt(sigmaxx**2 + sigmayy**2 + sigmazz**2 + &
                        2.0_real64*(sigmaxy**2 + sigmayz**2 + sigmaxz**2) ) )
            ! print *,'Time step # ',it,' out of ',time_step
            ! print *,'Time: ',(it-1)*DT,' seconds'
            ! print *,'Max vals for vx, vy, vz: ', maxval(vx), maxval(vy), maxval(vz)

            if (velocnorm(it) > stability_threshold) stop 'code became unstable and blew up'

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
        
        call write_array('velocity_norm.dat', source%time_steps, velocnorm )
        call write_array('stress_norm.dat', source%time_steps, stressnorm )
        
        deallocate(c11, c12, c13, c22, c23, c33, c44, c55, c66 )
        deallocate(rho)
        deallocate(kappa, alpha, acoef, bcoef, kappa_half, alpha_half, acoef_half, bcoef_half)
        deallocate(gamma_x, gamma_y, gamma_z, gamma_xy, gamma_yz, gamma_xz)
        deallocate(velocnorm, stressnorm)
        deallocate(srcx, srcy, srcz)
        deallocate(srcxx, srcyy, srczz, srcxz, srcxy, srcyz)
        deallocate(memory_dvx_dx, memory_dvx_dy, memory_dvx_dz )
        deallocate(memory_dvy_dx, memory_dvy_dy, memory_dvy_dz )
        deallocate(memory_dvz_dx, memory_dvz_dy, memory_dvz_dz )
        deallocate(memory_dsigmaxx_dx, memory_dsigmayy_dy, memory_dsigmazz_dz )
        deallocate(memory_dsigmaxy_dx, memory_dsigmaxy_dy, memory_dsigmaxz_dx )
        deallocate(memory_dsigmaxz_dz, memory_dsigmayz_dy, memory_dsigmayz_dz )
        deallocate(vx, vy, vz)
        deallocate(sigmaxx, sigmaxy, sigmaxz)
        deallocate(sigmayy, sigmayz, sigmazz)
        
    end subroutine seismic25iso

    ! =========================================================================
    subroutine seismic25(domain, source, density_method, SINGLE_OUTPUT)
        !--------------------------------------------------------------------------------------
        use constants
        
        implicit none

        ! Input arguments
        type(Domain_Type), intent(in) :: domain 
        type(Source_Type), intent(in) :: source 
        logical, intent(in), optional :: SINGLE_OUTPUT
        character(len=256), intent(in) :: density_method
        ! Local variables
        real(real64) :: rhoxx, rhoyx, rhozx, &
                        rhoxy, rhoyy, rhozy, &
                        rhoxz, rhoyz, rhozz
        
        real(real64), allocatable :: velocnorm(:), stressnorm(:)

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
        real(real64), allocatable :: kappa(:,:), alpha(:,:), acoef(:,:), bcoef(:,:), &
                        kappa_half(:,:), alpha_half(:,:), acoef_half(:,:), bcoef_half(:,:)
        
        ! Arrays for the PML damping factors
        real(real64), allocatable :: gamma_x(:,:), gamma_y(:,:), gamma_z(:,:)
        real(real64), allocatable :: gamma_xz(:,:), gamma_xy(:,:), gamma_yz(:,:)

        ! Source arrays
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:)
        real(real64), allocatable :: srcxx(:), srcxy(:), srcxz(:), srcyy(:), srcyz(:), srczz(:)
        
        real(real64), allocatable :: vx(:,:,:), vy(:,:,:), vz(:,:,:), &
                        sigmaxx(:,:,:), sigmaxy(:,:,:), sigmaxz(:,:,:), &
                        sigmayy(:,:,:), sigmayz(:,:,:), sigmazz(:,:,:)
        real(real64), allocatable :: &
                memory_dvx_dx1(:,:,:), memory_dvx_dy1(:,:,:), memory_dvx_dz1(:,:,:), &
                memory_dvy_dx1(:,:,:), memory_dvy_dy1(:,:,:), memory_dvy_dz1(:,:,:), &
                memory_dvz_dx1(:,:,:), memory_dvz_dy1(:,:,:), memory_dvz_dz1(:,:,:), &
                memory_dvx_dx2(:,:,:), memory_dvx_dy2(:,:,:), memory_dvx_dz2(:,:,:), &
                memory_dvy_dx2(:,:,:), memory_dvy_dy2(:,:,:), memory_dvy_dz2(:,:,:), &
                memory_dvz_dx2(:,:,:), memory_dvz_dy2(:,:,:), memory_dvz_dz2(:,:,:), &
                memory_dvx_dx3(:,:,:), memory_dvx_dy3(:,:,:), memory_dvx_dz3(:,:,:), &
                memory_dvy_dx3(:,:,:), memory_dvy_dy3(:,:,:), memory_dvy_dz3(:,:,:), &
                memory_dvz_dx3(:,:,:), memory_dvz_dy3(:,:,:), memory_dvz_dz3(:,:,:), &
                memory_dvx_dx4(:,:,:), memory_dvx_dy4(:,:,:), memory_dvx_dz4(:,:,:), &
                memory_dvy_dx4(:,:,:), memory_dvy_dy4(:,:,:), memory_dvy_dz4(:,:,:), &
                memory_dvz_dx4(:,:,:), memory_dvz_dy4(:,:,:), memory_dvz_dz4(:,:,:), &
                memory_dsigmaxx_dx(:,:,:), memory_dsigmayy_dy(:,:,:), memory_dsigmazz_dz(:,:,:), &
                memory_dsigmaxy_dx(:,:,:), memory_dsigmaxy_dy(:,:,:), memory_dsigmaxz_dx(:,:,:), &
                memory_dsigmaxz_dz(:,:,:), memory_dsigmayz_dy(:,:,:), memory_dsigmayz_dz(:,:,:)
        
        integer :: nx, ny, nz
        real(real64) :: dx, dy, dz, dt    
        
        ! Boolean flag to save as double precision or single precision
        logical :: SINGLE
        integer :: density_code
        
        ! The default data output is single precision unless SINGLE_OUTPUT is 
        ! set to .FALSE.
        if (present(SINGLE_OUTPUT)) then 
            SINGLE = SINGLE_OUTPUT 
        else
            SINGLE = .TRUE.
        endif
        
        select case (trim(adjustl(density_method)))
            case ('none'      ); density_code = 0
            case ('harmonic'  ); density_code = 1
            case ('geometric'); density_code = 2
            case ('arithmetic'); density_code = 3
        end select
        
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
                
        allocate(kappa(nx, nz), alpha(nx, nz), acoef(nx, nz), bcoef(nx, nz), &
                kappa_half(nx, nz), alpha_half(nx, nz), acoef_half(nx, nz), bcoef_half(nx, nz))
        allocate(gamma_x(nx, nz), gamma_y(nx, nz), gamma_z(nx, nz))
        allocate(gamma_xy(nx, nz), gamma_yz(nx, nz), gamma_xz(nx, nz))
        
        allocate(velocnorm(source%time_steps), stressnorm(source%time_steps))
        
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps))
        allocate(srcxx(source%time_steps), srcxy(source%time_steps), srcxz(source%time_steps), &
                 srcyy(source%time_steps), srcyz(source%time_steps), srczz(source%time_steps) )                
        
        ! Allocate more
        allocate(memory_dvx_dx1(nx,ny,nz), memory_dvx_dy1(nx,ny,nz), memory_dvx_dz1(nx,ny,nz) )
        allocate(memory_dvy_dx1(nx,ny,nz), memory_dvy_dy1(nx,ny,nz), memory_dvy_dz1(nx,ny,nz) )
        allocate(memory_dvz_dx1(nx,ny,nz), memory_dvz_dy1(nx,ny,nz), memory_dvz_dz1(nx,ny,nz) )
        
        allocate(memory_dvx_dx2(nx,ny,nz), memory_dvx_dy2(nx,ny,nz), memory_dvx_dz2(nx,ny,nz) )
        allocate(memory_dvy_dx2(nx,ny,nz), memory_dvy_dy2(nx,ny,nz), memory_dvy_dz2(nx,ny,nz) )
        allocate(memory_dvz_dx2(nx,ny,nz), memory_dvz_dy2(nx,ny,nz), memory_dvz_dz2(nx,ny,nz) )
        
        allocate(memory_dvx_dx3(nx,ny,nz), memory_dvx_dy3(nx,ny,nz), memory_dvx_dz3(nx,ny,nz) )
        allocate(memory_dvy_dx3(nx,ny,nz), memory_dvy_dy3(nx,ny,nz), memory_dvy_dz3(nx,ny,nz) )
        allocate(memory_dvz_dx3(nx,ny,nz), memory_dvz_dy3(nx,ny,nz), memory_dvz_dz3(nx,ny,nz) )
        
        allocate(memory_dvx_dx4(nx,ny,nz), memory_dvx_dy4(nx,ny,nz), memory_dvx_dz4(nx,ny,nz) )
        allocate(memory_dvy_dx4(nx,ny,nz), memory_dvy_dy4(nx,ny,nz), memory_dvy_dz4(nx,ny,nz) )
        allocate(memory_dvz_dx4(nx,ny,nz), memory_dvz_dy4(nx,ny,nz), memory_dvz_dz4(nx,ny,nz) )
        
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
        srcxx(:) = 0.0_real64
        srcxy(:) = 0.0_real64
        srcxz(:) = 0.0_real64
        srcyy(:) = 0.0_real64
        srcyz(:) = 0.0_real64
        srczz(:) = 0.0_real64
        srcx(:) = 0.0_real64
        srcy(:) = 0.0_real64
        srcz(:) = 0.0_real64
        
        if ( source%source_type == 'ac') then 
            call loadsource('seismicsourcex.dat', source%time_steps, srcx)
            call loadsource('seismicsourcey.dat', source%time_steps, srcy)
            call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        else
            call loadsource('seismicsourcexx.dat', source%time_steps, srcxx)
            call loadsource('seismicsourcexy.dat', source%time_steps, srcxy)
            call loadsource('seismicsourcexz.dat', source%time_steps, srcxz)
            call loadsource('seismicsourceyy.dat', source%time_steps, srcyy)
            call loadsource('seismicsourceyz.dat', source%time_steps, srcyz)
            call loadsource('seismicsourcezz.dat', source%time_steps, srczz)
        endif 
        
        ! ==================================== PML ====================================
        ! Initialize PML 
        kappa(:,:) = 1.0_real64
        kappa_half(:,:) = 1.0_real64
        alpha(:,:) = 0.0_real64
        alpha_half(:,:) = 0.0_real64
        acoef(:,:) = 0.0_real64
        acoef_half(:,:) = 0.0_real64

        ! ------------------------- Boundary Conditions -------------------------
        call material_rw2('kappa_cpml.dat', kappa, .TRUE.)
        call material_rw2('alpha_cpml.dat', alpha, .TRUE.)
        call material_rw2('acoef_cpml.dat', acoef, .TRUE.)
        call material_rw2('bcoef_cpml.dat', bcoef, .TRUE.)

        call material_rw2('kappa_half_cpml.dat', kappa_half, .TRUE.)
        call material_rw2('alpha_half_cpml.dat', alpha_half, .TRUE.)
        call material_rw2('acoef_half_cpml.dat', acoef_half, .TRUE.)
        call material_rw2('bcoef_half_cpml.dat', bcoef_half, .TRUE.)

        ! Load initial condition
        call material_rw3('initialconditionVx.dat', vx, .TRUE.)
        call material_rw3('initialconditionVy.dat', vy, .TRUE.)
        call material_rw3('initialconditionVz.dat', vz, .TRUE.)

        ! Initialize the stress values
        sigmaxx(:,:,:) = 0.0_real64
        sigmayy(:,:,:) = 0.0_real64
        sigmazz(:,:,:) = 0.0_real64
        sigmaxy(:,:,:) = 0.0_real64
        sigmaxz(:,:,:) = 0.0_real64
        sigmayz(:,:,:) = 0.0_real64

        ! PML
        memory_dvx_dx1(:,:,:) = 0.0_real64
        memory_dvx_dy1(:,:,:) = 0.0_real64
        memory_dvx_dz1(:,:,:) = 0.0_real64

        memory_dvy_dx1(:,:,:) = 0.0_real64
        memory_dvy_dy1(:,:,:) = 0.0_real64
        memory_dvy_dz1(:,:,:) = 0.0_real64

        memory_dvz_dx1(:,:,:) = 0.0_real64
        memory_dvz_dy1(:,:,:) = 0.0_real64 
        memory_dvz_dz1(:,:,:) = 0.0_real64
        
        memory_dvx_dx2(:,:,:) = 0.0_real64
        memory_dvx_dy2(:,:,:) = 0.0_real64
        memory_dvx_dz2(:,:,:) = 0.0_real64

        memory_dvy_dx2(:,:,:) = 0.0_real64
        memory_dvy_dy2(:,:,:) = 0.0_real64
        memory_dvy_dz2(:,:,:) = 0.0_real64

        memory_dvz_dx2(:,:,:) = 0.0_real64
        memory_dvz_dy2(:,:,:) = 0.0_real64 
        memory_dvz_dz2(:,:,:) = 0.0_real64
        
        memory_dvx_dx3(:,:,:) = 0.0_real64
        memory_dvx_dy3(:,:,:) = 0.0_real64
        memory_dvx_dz3(:,:,:) = 0.0_real64

        memory_dvy_dx3(:,:,:) = 0.0_real64
        memory_dvy_dy3(:,:,:) = 0.0_real64
        memory_dvy_dz3(:,:,:) = 0.0_real64

        memory_dvz_dx3(:,:,:) = 0.0_real64
        memory_dvz_dy3(:,:,:) = 0.0_real64 
        memory_dvz_dz3(:,:,:) = 0.0_real64
        
        memory_dvx_dx4(:,:,:) = 0.0_real64
        memory_dvx_dy4(:,:,:) = 0.0_real64
        memory_dvx_dz4(:,:,:) = 0.0_real64

        memory_dvy_dx4(:,:,:) = 0.0_real64
        memory_dvy_dy4(:,:,:) = 0.0_real64
        memory_dvy_dz4(:,:,:) = 0.0_real64

        memory_dvz_dx4(:,:,:) = 0.0_real64
        memory_dvz_dy4(:,:,:) = 0.0_real64 
        memory_dvz_dz4(:,:,:) = 0.0_real64

        memory_dsigmaxx_dx(:,:,:) = 0.0_real64
        memory_dsigmayy_dy(:,:,:) = 0.0_real64
        memory_dsigmazz_dz(:,:,:) = 0.0_real64
        memory_dsigmaxy_dx(:,:,:) = 0.0_real64
        memory_dsigmaxy_dy(:,:,:) = 0.0_real64
        memory_dsigmaxz_dx(:,:,:) = 0.0_real64
        memory_dsigmaxz_dz(:,:,:) = 0.0_real64
        memory_dsigmayz_dy(:,:,:) = 0.0_real64
        memory_dsigmayz_dz(:,:,:) = 0.0_real64
        
        ! Initialize velocnorm and stressnorm 
        velocnorm(:) = 0.0_real64
        stressnorm(:) = 0.0_real64 
        ! Do it 
        
        ! =============================== Forward Model ===============================
        do it = 1,source%time_steps
            !------------------------------------------------------------
            ! compute stress sigma and update memory variables for C-PML
            !------------------------------------------------------------
            ! Update in the x direction
            ! Loop 1
            do k = 3,nz-1
                do j = 3,ny-1
                    do i = 2,nx-2
                        ! Forward step on half grid
                        dvx_dx = ( 27.0_real64*(vx(i+1,j,k) - vx(i,j,k)) - (vx(i+2,j,k) - vx(i-1,j,k))) / (24.0_real64*dx)
                        dvy_dx = ( 27.0_real64*(vy(i+1,j,k) - vy(i,j,k)) - (vy(i+2,j,k) - vy(i-1,j,k))) / (24.0_real64*dx)
                        dvz_dx = ( 27.0_real64*(vz(i+1,j,k) - vz(i,j,k)) - (vz(i+2,j,k) - vz(i-1,j,k))) / (24.0_real64*dx)
                        ! Backward step on full grid
                        dvx_dy = ( 27.0_real64*(vx(i,j,k) - vx(i,j-1,k)) - (vx(i,j+1,k) - vx(i,j-2,k))) / (24.0_real64*dy)
                        dvy_dy = ( 27.0_real64*(vy(i,j,k) - vy(i,j-1,k)) - (vy(i,j+1,k) - vy(i,j-2,k))) / (24.0_real64*dy)
                        dvz_dy = ( 27.0_real64*(vz(i,j,k) - vz(i,j-1,k)) - (vz(i,j+1,k) - vz(i,j-2,k))) / (24.0_real64*dy)
                        ! Backward step on full grid
                        dvx_dz = ( 27.0_real64*(vx(i,j,k) - vx(i,j,k-1)) - (vx(i,j,k+1) - vx(i,j,k-2))) / (24.0_real64*dz)
                        dvy_dz = ( 27.0_real64*(vy(i,j,k) - vy(i,j,k-1)) - (vy(i,j,k+1) - vy(i,j,k-2))) / (24.0_real64*dz)
                        dvz_dz = ( 27.0_real64*(vz(i,j,k) - vz(i,j,k-1)) - (vz(i,j,k+1) - vz(i,j,k-2))) / (24.0_real64*dz)
                        
                        memory_dvx_dx1(i,j,k) = bcoef_half(i,k) * memory_dvx_dx1(i,j,k) + acoef_half(i,k) * dvx_dx
                        memory_dvy_dx1(i,j,k) = bcoef_half(i,k) * memory_dvy_dx1(i,j,k) + acoef_half(i,k) * dvy_dx
                        memory_dvz_dx1(i,j,k) = bcoef_half(i,k) * memory_dvz_dx1(i,j,k) + acoef_half(i,k) * dvz_dx
                        memory_dvy_dy1(i,j,k) = bcoef(i,k) * memory_dvy_dy1(i,j,k) + acoef(i,k) * dvy_dy
                        memory_dvx_dy1(i,j,k) = bcoef(i,k) * memory_dvx_dy1(i,j,k) + acoef(i,k) * dvx_dy
                        memory_dvz_dy1(i,j,k) = bcoef(i,k) * memory_dvz_dy1(i,j,k) + acoef(i,k) * dvz_dy
                        memory_dvz_dz1(i,j,k) = bcoef(i,k) * memory_dvz_dz1(i,j,k) + acoef(i,k) * dvz_dz
                        memory_dvx_dz1(i,j,k) = bcoef(i,k) * memory_dvx_dz1(i,j,k) + acoef(i,k) * dvx_dz
                        memory_dvy_dz1(i,j,k) = bcoef(i,k) * memory_dvy_dz1(i,j,k) + acoef(i,k) * dvy_dz

                        dvx_dx = dvx_dx / kappa_half(i,k) + memory_dvx_dx1(i,j,k)
                        dvy_dx = dvy_dx / kappa_half(i,k) + memory_dvy_dx1(i,j,k)
                        dvz_dx = dvz_dx / kappa_half(i,k) + memory_dvz_dx1(i,j,k)
                        dvy_dy = dvy_dy / kappa(i,k) + memory_dvy_dy1(i,j,k)
                        dvx_dy = dvx_dy / kappa(i,k) + memory_dvx_dy1(i,j,k)
                        dvz_dy = dvz_dy / kappa(i,k) + memory_dvz_dy1(i,j,k)
                        dvz_dz = dvz_dz / kappa(i,k) + memory_dvz_dz1(i,j,k)
                        dvx_dz = dvx_dz / kappa(i,k) + memory_dvx_dz1(i,j,k)
                        dvy_dz = dvy_dz / kappa(i,k) + memory_dvy_dz1(i,j,k)

                        sigmaxx(i,j,k) = ( sigmaxx(i,j,k) + &
                        (   c11(i,k) * dvx_dx + c12(i,k) * dvy_dy + c13(i,k) * dvz_dz + &
                            c14(i,k) * (dvy_dz + dvz_dy) + c15(i,k) * (dvx_dz + dvz_dx) + &
                            c16(i,k) * (dvy_dx + dvx_dy) ) * dt ) / &
                                (1 + gamma_x(i,k) * dt )

                        ! Full 3D will need a gradient in the y-direction
                        sigmayy(i,j,k) = ( sigmayy(i,j,k) + &
                        (   c12(i,k) * dvx_dx + c22(i,k) * dvy_dy + c23(i,k) * dvz_dz + &
                            c24(i,k) * (dvy_dz + dvz_dy) + c25(i,k) * (dvx_dz + dvz_dx) + &
                            c26(i,k) * (dvy_dx + dvx_dy) ) * dt ) / &
                                (1 + gamma_y(i,k) * dt )

                        sigmazz(i,j,k) = ( sigmazz(i,j,k) + &
                        (   c13(i,k) * dvx_dx + c23(i,k) * dvy_dy + c33(i,k) * dvz_dz + &
                            c34(i,k) * (dvy_dz + dvz_dy) + c35(i,k) * (dvx_dz + dvz_dx) + &
                            c36(i,k) * (dvy_dx + dvx_dy) ) * dt ) / &
                                (1 + gamma_z(i,k) * dt )

                    enddo
                enddo
            enddo
            ! Loop 2
            ! Update sigmaxy, x-direction is full nodes
            do k = 3,nz-1
                do j = 2,ny-2
                    do i = 3,nx-1
                        ! Backward step full grid
                        dvx_dx = ( 27.0_real64*(vx(i,j,k) - vx(i-1,j,k)) - (vx(i+1,j,k) - vx(i-2,j,k))) / (24.0_real64*dx)
                        dvy_dx = ( 27.0_real64*(vy(i,j,k) - vy(i-1,j,k)) - (vy(i+1,j,k) - vy(i-2,j,k))) / (24.0_real64*dx)
                        dvz_dx = ( 27.0_real64*(vz(i,j,k) - vz(i-1,j,k)) - (vz(i+1,j,k) - vz(i-2,j,k))) / (24.0_real64*dx)
                        ! Forward step half grid
                        dvx_dy = ( 27.0_real64*(vx(i,j+1,k) - vx(i,j,k)) - (vx(i,j+2,k) - vx(i,j-1,k))) / (24.0_real64*dy)
                        dvy_dy = ( 27.0_real64*(vy(i,j+1,k) - vy(i,j,k)) - (vy(i,j+2,k) - vy(i,j-1,k))) / (24.0_real64*dy)
                        dvz_dy = ( 27.0_real64*(vz(i,j+1,k) - vz(i,j,k)) - (vz(i,j+2,k) - vz(i,j-1,k))) / (24.0_real64*dy)
                        ! Backward step half grid
                        dvx_dz = ( 27.0_real64*(vx(i,j,k) - vx(i,j,k-1)) - (vx(i,j,k+1) - vx(i,j,k-2))) / (24.0_real64*dz)
                        dvy_dz = ( 27.0_real64*(vy(i,j,k) - vy(i,j,k-1)) - (vy(i,j,k+1) - vy(i,j,k-2))) / (24.0_real64*dz)
                        dvz_dz = ( 27.0_real64*(vz(i,j,k) - vz(i,j,k-1)) - (vz(i,j,k+1) - vz(i,j,k-2))) / (24.0_real64*dz)
                        
                        memory_dvx_dx2(i,j,k) = bcoef(i,k) * memory_dvx_dx2(i,j,k) + acoef(i,k) * dvx_dx
                        memory_dvy_dx2(i,j,k) = bcoef(i,k) * memory_dvy_dx2(i,j,k) + acoef(i,k) * dvy_dx
                        memory_dvz_dx2(i,j,k) = bcoef(i,k) * memory_dvz_dx2(i,j,k) + acoef(i,k) * dvz_dx
                        memory_dvy_dy2(i,j,k) = bcoef_half(i,k) * memory_dvy_dy2(i,j,k) + acoef_half(i,k) * dvy_dy
                        memory_dvx_dy2(i,j,k) = bcoef_half(i,k) * memory_dvx_dy2(i,j,k) + acoef_half(i,k) * dvx_dy
                        memory_dvz_dy2(i,j,k) = bcoef_half(i,k) * memory_dvz_dy2(i,j,k) + acoef_half(i,k) * dvz_dy
                        memory_dvz_dz2(i,j,k) = bcoef_half(i,k) * memory_dvz_dz2(i,j,k) + acoef_half(i,k) * dvz_dz
                        memory_dvx_dz2(i,j,k) = bcoef_half(i,k) * memory_dvx_dz2(i,j,k) + acoef_half(i,k) * dvx_dz
                        memory_dvy_dz2(i,j,k) = bcoef_half(i,k) * memory_dvy_dz2(i,j,k) + acoef_half(i,k) * dvy_dz

                        dvx_dx = dvx_dx / kappa(i,k) + memory_dvx_dx2(i,j,k)
                        dvy_dx = dvy_dx / kappa(i,k) + memory_dvy_dx2(i,j,k)
                        dvz_dx = dvz_dx / kappa(i,k) + memory_dvz_dx2(i,j,k)
                        dvy_dy = dvy_dy / kappa_half(i,k) + memory_dvy_dy2(i,j,k)
                        dvx_dy = dvx_dy / kappa_half(i,k) + memory_dvx_dy2(i,j,k)
                        dvz_dy = dvz_dy / kappa_half(i,k) + memory_dvz_dy2(i,j,k)
                        dvz_dz = dvz_dz / kappa_half(i,k) + memory_dvz_dz2(i,j,k)
                        dvy_dz = dvy_dz / kappa_half(i,k) + memory_dvy_dz2(i,j,k)

                        sigmaxy(i,j,k) = ( sigmaxy(i,j,k) + &
                        (   c16(i,k) * dvx_dx + c26(i,k) * dvy_dy + c36(i,k) * dvz_dz + &
                            c46(i,k) * (dvy_dz + dvz_dy) + c56(i,k) * (dvx_dz + dvz_dx) + &
                            c66(i,k) * (dvy_dx + dvx_dy) ) * dt ) / &
                                (1 + gamma_xy(i,k) * dt )

                    enddo
                enddo
            enddo
            ! Loop 3
            ! Update sigmaxz, z-direction is full nodes
            do k = 2,nz-2
                do j = 3,ny-1
                    do i = 3,nx-1
                        ! Backward difference full grid
                        dvx_dx = ( 27.0_real64*(vx(i,j,k) - vx(i-1,j,k)) - (vx(i+1,j,k) - vx(i-2,j,k))) / (24.0_real64*dx)
                        dvy_dx = ( 27.0_real64*(vy(i,j,k) - vy(i-1,j,k)) - (vy(i+1,j,k) - vy(i-2,j,k))) / (24.0_real64*dx)
                        dvz_dx = ( 27.0_real64*(vz(i,j,k) - vz(i-1,j,k)) - (vz(i+1,j,k) - vz(i-2,j,k))) / (24.0_real64*dx)
                        ! Backward difference full grid
                        dvx_dy = ( 27.0_real64*(vx(i,j,k) - vx(i,j-1,k)) - (vx(i,j+1,k) - vx(i,j-2,k))) / (24.0_real64*dy)
                        dvy_dy = ( 27.0_real64*(vy(i,j,k) - vy(i,j-1,k)) - (vy(i,j+1,k) - vy(i,j-2,k))) / (24.0_real64*dy)
                        dvz_dy = ( 27.0_real64*(vz(i,j,k) - vz(i,j-1,k)) - (vz(i,j+1,k) - vz(i,j-2,k))) / (24.0_real64*dy)
                        ! Forward difference half grid
                        dvx_dz = ( 27.0_real64*(vx(i,j,k+1) - vx(i,j,k)) - (vx(i,j,k+2) - vx(i,j,k-1))) / (24.0_real64*dz)
                        dvy_dz = ( 27.0_real64*(vy(i,j,k+1) - vy(i,j,k)) - (vy(i,j,k+2) - vy(i,j,k-1))) / (24.0_real64*dz)
                        dvz_dz = ( 27.0_real64*(vz(i,j,k+1) - vz(i,j,k)) - (vz(i,j,k+2) - vz(i,j,k-1))) / (24.0_real64*dz)
                        
                        memory_dvx_dx3(i,j,k) = bcoef(i,k) * memory_dvx_dx3(i,j,k) + acoef(i,k) * dvx_dx
                        memory_dvy_dx3(i,j,k) = bcoef(i,k) * memory_dvy_dx3(i,j,k) + acoef(i,k) * dvy_dx
                        memory_dvz_dx3(i,j,k) = bcoef(i,k) * memory_dvz_dx3(i,j,k) + acoef(i,k) * dvz_dx
                        memory_dvy_dy3(i,j,k) = bcoef(i,k) * memory_dvy_dy3(i,j,k) + acoef(i,k) * dvy_dy
                        memory_dvx_dy3(i,j,k) = bcoef(i,k) * memory_dvx_dy3(i,j,k) + acoef(i,k) * dvx_dy
                        memory_dvz_dy3(i,j,k) = bcoef(i,k) * memory_dvz_dy3(i,j,k) + acoef(i,k) * dvz_dy
                        memory_dvz_dz3(i,j,k) = bcoef_half(i,k) * memory_dvz_dz3(i,j,k) + acoef_half(i,k) * dvz_dz
                        memory_dvx_dz3(i,j,k) = bcoef_half(i,k) * memory_dvx_dz3(i,j,k) + acoef_half(i,k) * dvx_dz
                        memory_dvy_dz3(i,j,k) = bcoef_half(i,k) * memory_dvy_dz3(i,j,k) + acoef_half(i,k) * dvy_dz

                        dvx_dx = dvx_dx / kappa(i,k) + memory_dvx_dx3(i,j,k)
                        dvy_dx = dvy_dx / kappa(i,k) + memory_dvy_dx3(i,j,k)
                        dvz_dx = dvz_dx / kappa(i,k) + memory_dvz_dx3(i,j,k) 
                        dvy_dy = dvy_dy / kappa(i,k) + memory_dvy_dy3(i,j,k)
                        dvx_dy = dvx_dy / kappa(i,k) + memory_dvx_dy3(i,j,k)
                        dvz_dy = dvz_dy / kappa(i,k) + memory_dvz_dy3(i,j,k)
                        dvz_dz = dvz_dz / kappa_half(i,k) + memory_dvz_dz3(i,j,k)
                        dvx_dz = dvx_dz / kappa_half(i,k) + memory_dvx_dz3(i,j,k)
                        dvy_dz = dvy_dz / kappa_half(i,k) + memory_dvy_dz3(i,j,k)

                        sigmaxz(i,j,k) = ( sigmaxz(i,j,k) + &
                            (   c15(i,k) * dvx_dx + c25(i,k) * dvy_dy + c35(i,k) * dvz_dz + &
                                c45(i,k) * ( dvy_dz + dvz_dy) + c55(i,k) * ( dvx_dz + dvz_dx) + &
                                c56(i,k) * ( dvy_dx + dvx_dy) ) * dt  ) / &
                                    (1 + gamma_xz(i,k) * dt )

                    enddo
                enddo
                ! Loop 4
                !   ! update sigmayz, y-direction is full nodes
                do j = 2,ny-2
                    do i = 2,nx-2
                        ! Forward difference half grid
                        dvx_dx = ( 27.0_real64*(vx(i+1,j,k) - vx(i,j,k)) - (vx(i+2,j,k) - vx(i-1,j,k))) / (24.0_real64*dx)
                        dvy_dx = ( 27.0_real64*(vy(i+1,j,k) - vy(i,j,k)) - (vy(i+2,j,k) - vy(i-1,j,k))) / (24.0_real64*dx)
                        dvz_dx = ( 27.0_real64*(vz(i+1,j,k) - vz(i,j,k)) - (vz(i+2,j,k) - vz(i-1,j,k))) / (24.0_real64*dx)
                        ! Forward difference half grid
                        dvx_dy = ( 27.0_real64*(vx(i,j+1,k) - vx(i,j,k)) - (vx(i,j+2,k) - vx(i,j-1,k))) / (24.0_real64*dy)
                        dvy_dy = ( 27.0_real64*(vy(i,j+1,k) - vy(i,j,k)) - (vy(i,j+2,k) - vy(i,j-1,k))) / (24.0_real64*dy)
                        dvz_dy = ( 27.0_real64*(vz(i,j+1,k) - vz(i,j,k)) - (vz(i,j+2,k) - vz(i,j-1,k))) / (24.0_real64*dy)
                        ! Forwar difference half grid
                        dvx_dz = ( 27.0_real64*(vx(i,j,k+1) - vx(i,j,k)) - (vx(i,j,k+2) - vx(i,j,k-1))) / (24.0_real64*dz)
                        dvy_dz = ( 27.0_real64*(vy(i,j,k+1) - vy(i,j,k)) - (vy(i,j,k+2) - vy(i,j,k-1))) / (24.0_real64*dz)
                        dvz_dz = ( 27.0_real64*(vz(i,j,k+1) - vz(i,j,k)) - (vz(i,j,k+2) - vz(i,j,k-1))) / (24.0_real64*dz)
                        
                        memory_dvx_dx4(i,j,k) = bcoef_half(i,k) * memory_dvx_dx4(i,j,k) + acoef_half(i,k) * dvx_dx
                        memory_dvy_dx4(i,j,k) = bcoef_half(i,k) * memory_dvy_dx4(i,j,k) + acoef_half(i,k) * dvy_dx
                        memory_dvz_dx4(i,j,k) = bcoef_half(i,k) * memory_dvz_dx4(i,j,k) + acoef_half(i,k) * dvz_dx
                        memory_dvy_dy4(i,j,k) = bcoef_half(i,k) * memory_dvy_dy4(i,j,k) + acoef_half(i,k) * dvy_dy
                        memory_dvx_dy4(i,j,k) = bcoef_half(i,k) * memory_dvx_dy4(i,j,k) + acoef_half(i,k) * dvx_dy
                        memory_dvz_dy4(i,j,k) = bcoef_half(i,k) * memory_dvz_dy4(i,j,k) + acoef_half(i,k) * dvz_dy
                        memory_dvz_dz4(i,j,k) = bcoef_half(i,k) * memory_dvz_dz4(i,j,k) + acoef_half(i,k) * dvz_dz
                        memory_dvx_dz4(i,j,k) = bcoef_half(i,k) * memory_dvx_dz4(i,j,k) + acoef_half(i,k) * dvx_dz
                        memory_dvy_dz4(i,j,k) = bcoef_half(i,k) * memory_dvy_dz4(i,j,k) + acoef_half(i,k) * dvy_dz

                        dvx_dx = dvx_dx / kappa_half(i,k) + memory_dvx_dx4(i,j,k)
                        dvy_dx = dvy_dx / kappa_half(i,k) + memory_dvy_dx4(i,j,k)
                        dvz_dx = dvz_dx / kappa_half(i,k) + memory_dvz_dx4(i,j,k)
                        dvy_dy = dvy_dy / kappa_half(i,k) + memory_dvy_dy4(i,j,k)
                        dvx_dy = dvx_dy / kappa_half(i,k) + memory_dvx_dy4(i,j,k)
                        dvz_dy = dvz_dy / kappa_half(i,k) + memory_dvz_dy4(i,j,k)
                        dvz_dz = dvz_dz / kappa_half(i,k) + memory_dvz_dz4(i,j,k)
                        dvx_dz = dvx_dz / kappa_half(i,k) + memory_dvx_dz4(i,j,k)
                        dvy_dz = dvy_dz / kappa_half(i,k) + memory_dvy_dz4(i,j,k)

                        sigmayz(i,j,k) = ( sigmayz(i,j,k)  + &
                            (   c14(i,k) * dvx_dx + c24(i,k) * dvy_dy + c34(i,k) * dvz_dz + &
                                c44(i,k) * ( dvy_dz + dvz_dy) + c45(i,k) * ( dvx_dz + dvz_dx) + &
                                c46(i,k) * ( dvy_dx + dvx_dy) ) * dt  ) / & 
                                    (1 + gamma_yz(i,k) * dt )
                    enddo
                enddo
            enddo

            !--------------------------------------------------------
            ! compute velocity and update memory variables for C-PML
            !--------------------------------------------------------
            ! Loop 5
            do k = 3,nz-1
                do j = 3,ny-1
                    do i = 3,nx-1
                        ! ds1/dx, ds6/dy, ds5,dz
                        rhoxx = face_density(rho(i,k), rho(i-1,k), density_code)
                        rhoyx = face_density(rho(i,k), rho(i,k), density_code) 
                        rhozx = face_density(rho(i,k), rho(i,k-1), density_code) 
                        ! Backward difference half grid
                        dsigmaxx_dx = (27.0_real64*sigmaxx(i,j,k) - 27.0_real64*sigmaxx(i-1,j,k) + sigmaxx(i-2,j,k)) / (24.0_real64*dx)
                        dsigmaxy_dy = (27.0_real64*sigmaxy(i,j,k) - 27.0_real64*sigmaxy(i,j-1,k) + sigmaxy(i,j-2,k)) / (24.0_real64*dy)
                        dsigmaxz_dz = (27.0_real64*sigmaxz(i,j,k) - 27.0_real64*sigmaxz(i,j,k-1) + sigmaxz(i,j,k-2)) / (24.0_real64*dz)
                        
                        memory_dsigmaxx_dx(i,j,k) = bcoef_half(i,k) * memory_dsigmaxx_dx(i,j,k) + acoef_half(i,k) * dsigmaxx_dx
                        memory_dsigmaxy_dy(i,j,k) = bcoef_half(i,k) * memory_dsigmaxy_dy(i,j,k) + acoef_half(i,k) * dsigmaxy_dy
                        memory_dsigmaxz_dz(i,j,k) = bcoef_half(i,k) * memory_dsigmaxz_dz(i,j,k) + acoef_half(i,k) * dsigmaxz_dz

                        dsigmaxx_dx = dsigmaxx_dx / kappa_half(i,k) + memory_dsigmaxx_dx(i,j,k)
                        dsigmaxy_dy = dsigmaxy_dy / kappa_half(i,k) + memory_dsigmaxy_dy(i,j,k)
                        dsigmaxz_dz = dsigmaxz_dz / kappa_half(i,k) + memory_dsigmaxz_dz(i,j,k) 

                        vx(i,j,k) = vx(i,j,k) + &
                            (dsigmaxx_dx/rhoxx + dsigmaxy_dy/rhoyx + dsigmaxz_dz/rhozx) * dt 
                    enddo
                enddo

                do j = 2,ny-2
                    do i = 2,nx-2
                        ! ds6/dx, ds2/dy, ds4/dz
                        rhoxy = face_density(rho(i,k), rho(i+1,k), density_code)
                        rhoyy = face_density(rho(i,k), rho(i,k), density_code) 
                        rhozy = face_density(rho(i,k), rho(i,k-1), density_code)
                        ! Forward difference half grid
                        dsigmaxy_dx = (-27.0_real64*sigmaxy(i,j,k) + 27.0_real64*sigmaxy(i+1,j,k) - sigmaxy(i+2,j,k)) / (24.0_real64*dx)
                        dsigmayy_dy = (-27.0_real64*sigmayy(i,j,k) + 27.0_real64*sigmayy(i,j+1,k) - sigmayy(i,j+2,k)) / (24.0_real64*dy)
                        ! Backward difference full grid
                        dsigmayz_dz = (27.0_real64*sigmayz(i,j,k) - 27.0_real64*sigmayz(i,j,k-1) + sigmayz(i,j,k-2)) / (24.0_real64*dz)
                        
                        memory_dsigmaxy_dx(i,j,k) = bcoef_half(i,k) * memory_dsigmaxy_dx(i,j,k) + acoef_half(i,k) * dsigmaxy_dx
                        memory_dsigmayy_dy(i,j,k) = bcoef_half(i,k) * memory_dsigmayy_dy(i,j,k) + acoef_half(i,k) * dsigmayy_dy
                        memory_dsigmayz_dz(i,j,k) = bcoef(i,k) * memory_dsigmayz_dz(i,j,k) + acoef(i,k) * dsigmayz_dz
                        
                        dsigmaxy_dx = dsigmaxy_dx / kappa_half(i,k) + memory_dsigmaxy_dx(i,j,k)
                        dsigmayy_dy = dsigmayy_dy / kappa_half(i,k) + memory_dsigmayy_dy(i,j,k)
                        dsigmayz_dz = dsigmayz_dz / kappa(i,k) + memory_dsigmayz_dz(i,j,k)

                        vy(i,j,k) = vy(i,j,k) + &
                            (dsigmaxy_dx/rhoxy + dsigmayy_dy/rhoyy + dsigmayz_dz/rhozy) * dt
                    enddo
                enddo
            enddo

            do k = 2,nz-2
                do j = 3,ny-1
                    do i = 2,nx-2
                        ! ds5/dx, ds4/dy, ds3/dz
                        rhoxz = face_density(rho(i,k), rho(i+1,k), density_code)
                        rhoyz = face_density(rho(i,k), rho(i,k), density_code) 
                        rhozz = face_density(rho(i,k), rho(i,k+1), density_code)
                        
                        ! Forward difference half grid
                        dsigmaxz_dx = (-27.0_real64*sigmaxz(i,j,k) + 27.0_real64*sigmaxz(i+1,j,k) - sigmaxz(i+2,j,k)) / (24.0_real64*dx)
                        ! Backward difference full grid
                        dsigmayz_dy = (27.0_real64*sigmayz(i,j,k) - 27.0_real64*sigmayz(i,j-1,k) + sigmayz(i,j-2,k)) / (24.0_real64*dy)
                        ! Forward difference half grid
                        dsigmazz_dz = (-27.0_real64*sigmazz(i,j,k) + 27.0_real64*sigmazz(i,j,k+1) - sigmazz(i,j,k+2)) / (24.0_real64*dz)
                        
                        memory_dsigmaxz_dx(i,j,k) = bcoef_half(i,k) * memory_dsigmaxz_dx(i,j,k) + acoef_half(i,k) * dsigmaxz_dx
                        memory_dsigmayz_dy(i,j,k) = bcoef(i,k) * memory_dsigmayz_dy(i,j,k) + acoef(i,k) * dsigmayz_dy
                        memory_dsigmazz_dz(i,j,k) = bcoef_half(i,k) * memory_dsigmazz_dz(i,j,k) + acoef_half(i,k) * dsigmazz_dz

                        dsigmaxz_dx = dsigmaxz_dx / kappa_half(i,k) + memory_dsigmaxz_dx(i,j,k)
                        dsigmayz_dy = dsigmayz_dy / kappa(i,k) + memory_dsigmayz_dy(i,j,k)
                        dsigmazz_dz = dsigmazz_dz / kappa_half(i,k) + memory_dsigmazz_dz(i,j,k)

                        vz(i,j,k) = vz(i,j,k) + &
                            (dsigmaxz_dx/rhoxz + dsigmayz_dy/rhoyz + dsigmazz_dz/rhozz) * &
                            dt !/ deltarho !rho(i,k)

                    enddo
                enddo
            enddo
            
            sigmaxx(isource,jsource,ksource) = sigmaxx(isource,jsource,ksource) + srcxx(it) / rho(isource,ksource)  
            sigmaxy(isource+1,jsource+1,ksource) = sigmaxy(isource+1,jsource+1,ksource) + srcxy(it) / rho(isource+1,ksource+1)
            sigmaxz(isource+1,jsource,ksource+1) = sigmaxz(isource+1,jsource,ksource+1) + srcxz(it) / rho(isource+1,ksource)  
            sigmayy(isource,jsource,ksource) = sigmayy(isource,jsource,ksource) + srcyy(it) / rho(isource,ksource)
            sigmayz(isource,jsource+1,ksource+1) = sigmayz(isource,jsource+1,ksource+1) + srcyz(it) / rho(isource,ksource+1)  
            sigmazz(isource,jsource,ksource) = sigmazz(isource,jsource,ksource) + srczz(it) / rho(isource,ksource)
            vx(isource+1,jsource,ksource) = vx(isource+1,jsource,ksource) + &
                    srcx(it) * dt / rho(isource+1,ksource)
            vy(isource,jsource+1,ksource) = vy(isource,jsource+1,ksource) + &
                    srcy(it) * dt / rho(isource,ksource)
            vz(isource,jsource,ksource+1) = vz(isource,jsource,ksource+1) + &
                    srcz(it) * dt / rho(isource,ksource+1)

            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            vx(1,:,:) = 0.0_real64
            vy(1,:,:) = 0.0_real64
            vz(1,:,:) = 0.0_real64

            vx(:,1,:) = 0.0_real64
            vy(:,1,:) = 0.0_real64
            vz(:,1,:) = 0.0_real64

            vx(:,:,1) = 0.0_real64
            vy(:,:,1) = 0.0_real64
            vz(:,:,1) = 0.0_real64

            vx(nx,:,:) = 0.0_real64
            vy(nx,:,:) = 0.0_real64
            vz(nx,:,:) = 0.0_real64

            vx(:,ny,:) = 0.0_real64
            vy(:,ny,:) = 0.0_real64
            vz(:,ny,:) = 0.0_real64

            vx(:,:,nz) = 0.0_real64
            vy(:,:,nz) = 0.0_real64
            vz(:,:,nz) = 0.0_real64

            ! check norm of velocity to make sure the solution isn't diverging
            velocnorm(it) = maxval( sqrt(vx**2 + vy**2 + vz**2) )
            stressnorm(it) = maxval( &
                        sqrt(sigmaxx**2 + sigmayy**2 + sigmazz**2 + &
                        2.0_real64*(sigmaxy**2 + sigmayz**2 + sigmaxz**2) ) )
            ! print *,'Time step # ',it,' out of ',time_step
            ! print *,'Time: ',(it-1)*DT,' seconds'
            ! print *,'Max vals for vx, vy, vz: ', maxval(vx), maxval(vy), maxval(vz)

            if (velocnorm(it) > stability_threshold) stop 'code became unstable and blew up'

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
        
        call write_array('velocity_norm.dat', source%time_steps, velocnorm )
        call write_array('stress_norm.dat', source%time_steps, stressnorm )
        
        deallocate(c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, &
                    c33, c34, c35, c36, c44, c45, c46, c55, c56, c66 )
        deallocate(rho)
        deallocate(kappa, alpha, acoef, bcoef, kappa_half, alpha_half, acoef_half, bcoef_half)
        deallocate(gamma_x, gamma_y, gamma_z, gamma_xy, gamma_yz, gamma_xz)
        deallocate(velocnorm, stressnorm)
        deallocate(srcx, srcy, srcz)
        deallocate(srcxx, srcyy, srczz, srcxz, srcxy, srcyz)
        deallocate(memory_dvx_dx1, memory_dvx_dy1, memory_dvx_dz1 )
        deallocate(memory_dvy_dx1, memory_dvy_dy1, memory_dvy_dz1 )
        deallocate(memory_dvz_dx1, memory_dvz_dy1, memory_dvz_dz1 )
        deallocate(memory_dvx_dx2, memory_dvx_dy2, memory_dvx_dz2 )
        deallocate(memory_dvy_dx2, memory_dvy_dy2, memory_dvy_dz2 )
        deallocate(memory_dvz_dx2, memory_dvz_dy2, memory_dvz_dz2 )
        deallocate(memory_dvx_dx3, memory_dvx_dy3, memory_dvx_dz3 )
        deallocate(memory_dvy_dx3, memory_dvy_dy3, memory_dvy_dz3 )
        deallocate(memory_dvz_dx3, memory_dvz_dy3, memory_dvz_dz3 )
        deallocate(memory_dvx_dx4, memory_dvx_dy4, memory_dvx_dz4 )
        deallocate(memory_dvy_dx4, memory_dvy_dy4, memory_dvy_dz4 )
        deallocate(memory_dvz_dx4, memory_dvz_dy4, memory_dvz_dz4 )
        deallocate(memory_dsigmaxx_dx, memory_dsigmayy_dy, memory_dsigmazz_dz )
        deallocate(memory_dsigmaxy_dx, memory_dsigmaxy_dy, memory_dsigmaxz_dx )
        deallocate(memory_dsigmaxz_dz, memory_dsigmayz_dy, memory_dsigmayz_dz )
        deallocate(vx, vy, vz)
        deallocate(sigmaxx, sigmaxy, sigmaxz)
        deallocate(sigmayy, sigmayz, sigmazz)
        
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
        real(real64), allocatable :: epsilon11(:,:), epsilon13(:,:), epsilon33(:,:), &
                                            sigmax(:,:), sigmaz(:,:)

        ! real(real64) :: DT
        real(real64), allocatable :: srcx(:), srcz(:) ! The vector time series of the source
        integer :: isource, jsource, i, j, it

        ! Coefficients for the finite difference scheme
        ! real(real64), allocatable :: caEx(:,:), cbEx(:,:)
        ! real(real64), allocatable :: caEz(:,:), cbEz(:,:)
        real(real64) :: daHy, dbHy 
        real(real64), allocatable ::det(:,:)
        real(real64), allocatable :: aEx(:,:), bEx(:,:), dEx(:,:), eEx(:,:), &
                                    aEz(:,:), bEz(:,:), dEz(:,:), eEz(:,:), &
                                    caEx(:,:), cbEx(:,:), caEz(:,:), cbEz(:,:)
        real(real64) :: dEx_dz, dEz_dx, dHy_dz, dHy_dx

        ! 1D arrays for the damping profiles
        real(real64), allocatable :: kappa(:,:), alpha(:,:), acoef(:,:), bcoef(:,:), kappa_half(:,:), alpha_half(:,:), acoef_half(:,:), bcoef_half(:,:)
        
        real(real64), allocatable :: Ex(:,:), Ez(:,:), Hy(:,:) 
        real(real64), allocatable :: Ex_old(:,:), Ez_old(:,:) 
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
        
        ! real(real64) :: m11,m12,m21,m22,n11,n12,n21,n22,detM
        ! real(real64) :: rhs1,rhs2, exn, ezn, exnp1, eznp1

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
        allocate(kappa(nx,nz), alpha(nx,nz), acoef(nx,nz), bcoef(nx,nz), &
                kappa_half(nx,nz), alpha_half(nx,nz), acoef_half(nx,nz), bcoef_half(nx,nz))
        allocate(srcx(source%time_steps), srcz(source%time_steps))
        
        ! Allocate more
        allocate(epsilon11(nx, nz), epsilon13(nx,nz), epsilon33(nx, nz))
        allocate(sigmax(nx, nz), sigmaz(nx, nz))
        allocate(aEx(nx,nz), bEx(nx,nz), dEx(nx,nz), eEx(nx,nz) )
        allocate(aEz(nx,nz), bEz(nx,nz), dEz(nx,nz), eEz(nx,nz) )
        allocate(caEx(nx,nz), cbEx(nx,nz), caEz(nx,nz), cbEz(nx,nz))
        allocate(det(nx,nz))
        allocate(memory_dEz_dx(nx, nz), memory_dEx_dz(nx, nz))
        allocate(memory_dHy_dx(nx, nz), memory_dHy_dz(nx, nz))
        allocate(Ex(nx, nz), Ez(nx, nz), Hy(nx, nz))
        allocate(Ex_old(nx,nz), Ez_old(nx,nz) )
            
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
        kappa(:,:) = 1.0_real64
        kappa_half(:,:) = 1.0_real64
        alpha(:,:) = 0.0_real64
        alpha_half(:,:) = 0.0_real64
        acoef(:,:) = 0.0_real64
        acoef_half(:,:) = 0.0_real64
        bcoef(:,:) = 0.0_real64 
        bcoef_half(:,:) = 0.0_real64 
        
        call material_rw2('kappa_cpml.dat', kappa, .TRUE.)
        call material_rw2('alpha_cpml.dat', alpha, .TRUE.)
        call material_rw2('acoef_cpml.dat', acoef, .TRUE.)
        call material_rw2('bcoef_cpml.dat', bcoef, .TRUE.)

        call material_rw2('kappa_half_cpml.dat', kappa_half, .TRUE.)
        call material_rw2('alpha_half_cpml.dat', alpha_half, .TRUE.)
        call material_rw2('acoef_half_cpml.dat', acoef_half, .TRUE.)
        call material_rw2('bcoef_half_cpml.dat', bcoef_half, .TRUE.)


        ! ----------------------------------------------------------------------
        ! Load initial conditions and initialize variables
        call material_rw2('initialconditionEx.dat', Ex, .TRUE.)
        call material_rw2('initialconditionEz.dat', Ez, .TRUE.)
    
        Hy(:,:) = 0._real64
        ! PML
        memory_dEx_dz(:,:) = 0.0_real64
        memory_dEz_dx(:,:) = 0.0_real64
        memory_dHy_dx(:,:) = 0.0_real64
        memory_dHy_dz(:,:) = 0.0_real64
        
        ! ----------------------------------------------------------------------
        ! Compute the coefficients of the FD scheme. First scale the relative 
        ! permittivity and permeabilities to get the absolute values 
        ! daHy = dt/(4.0_real64*mu0*mu)
        dbHy = dt/mu0 !dt/(mu*mu*dx*(1+daHy) ) 
        daHy = 1.0_real64 ! (1-daHy)/(1+daHy) ! 
        
        
        epsilon11 = eps11 * eps0
        epsilon13 = eps13 * eps0
        epsilon33 = eps33 * eps0
        
        det = (epsilon11*epsilon33 - epsilon13*epsilon13)
        where (det <= 1.0d-18*eps0*eps0)
            det = 1.0d-18*eps0*eps0
        end where
        
        aEx = epsilon33 / det 
        bEx = -epsilon13 / det 
        dEx = (epsilon33 * sig11 - epsilon13 * sig13) / det  
        eEx = (epsilon33 * sig13 - epsilon13 * sig33) / det
        
        aEz = -epsilon13 / det 
        bEz = epsilon11 / det 
        dEz = (epsilon11 * sig13 - epsilon13 * sig11) / det  
        eEz = (epsilon11 * sig33 - epsilon13 * sig13) / det
        
        caEx(:,:) = ( 1.0_real64 - sig11 * dt / ( 2.0d0 * epsilon11 ) ) / &
                    ( 1.0_real64 + sig11 * dt / ( 2.0d0 * epsilon11 ) )
        cbEx(:,:) = (dt / epsilon11 ) / ( 1.0_real64 + sig11 * dt / ( 2.0d0 * epsilon11 ) )

        caEz(:,:) = ( 1.0_real64 - sig33 * dt / ( 2.0d0 * epsilon33 ) ) / &
                    ( 1.0_real64 + sig33 * dt / ( 2.0d0 * epsilon33 ) )
        cbEz(:,:) = (dt / epsilon33 ) / ( 1.0_real64 + sig33 * dt / ( 2.0d0 * epsilon33 ) )
        !---
        !---  beginning of time loop
        !---
        
        Ex_old = Ex 
        Ez_old = Ez
        
        do it = 1,source%time_steps
            velocnorm = 0.0_real64
            !$omp parallel private(i, j, dEx_dz, dEz_dx, dHy_dz, dHy_dx) &
            !$omp& shared(Ex, Ez, Hy, memory_dEx_dz, memory_dEz_dx, memory_dHy_dz, memory_dHy_dx, bcoef, acoef, kappa, aEz, bEz, aEx, bEx, dEz, eEz, dEx, eEx, daHy, dbHy) & 
            !--------------------------------------------------------
            ! compute magnetic field and update memory variables for C-PML
            !--------------------------------------------------------
            !$omp do
            do j = 1,nz-1 
                do i = 1,nx-1
                
                    ! Values needed for the magnetic field updates
                    dEx_dz = ( Ex(i,j+1) - Ex(i,j) )/dz
                    dEz_dx = ( Ez(i+1,j) - Ez(i,j) )/dx
                    
                    memory_dEx_dz(i,j) = bcoef(i,j) * memory_dEx_dz(i,j) + acoef(i,j) * dEx_dz
                    memory_dEz_dx(i,j) = bcoef_half(i,j) * memory_dEz_dx(i,j) + acoef_half(i,j) * dEz_dx
                    
                    dEz_dx = dEz_dx/ kappa_half(i,j) + memory_dEz_dx(i,j)
                    dEx_dz = dEx_dz/ kappa(i,j) + memory_dEx_dz(i,j)
                    
                    ! Now update the Magnetic field
                    Hy(i,j) = daHy*Hy(i,j) + dbHy*( dEz_dx - dEx_dz )
                    ! velocnorm = max(velocnorm, sqrt(Ex(i, j)**2 + Ez(i, j)**2))
                    
                enddo  
            enddo
            !$omp end do 
            
            ! Electric field and update memory variables for C-PML
            ! Compute the differences in the y-direction
            !$omp do
            do j = 2,nz
                do i = 2,nx
                    dHy_dz = ( Hy(i,j) - Hy(i,j-1) )/dz ! this is nz-1 length vector
                    dHy_dx = ( Hy(i,j) - Hy(i-1,j) )/dx
                    
                    memory_dHy_dx(i,j) = bcoef_half(i,j) * memory_dHy_dx(i,j) + acoef_half(i,j) * dHy_dx
                    memory_dHy_dz(i,j) = bcoef(i,j) * memory_dHy_dz(i,j) + acoef(i,j) * dHy_dz
                    
                    dHy_dz = dHy_dz/kappa(i,j) + memory_dHy_dz(i,j)
                    dHy_dx = dHy_dx/kappa_half(i,j) + memory_dHy_dx(i,j)
                    
                    Ez(i,j) = Ez_old(i,j) + ( -aEz(i,j)*dHy_dz + bEz(i,j)*dHy_dx - dEz(i,j)*Ex_old(i,j) - eEz(i,j)*Ez_old(i,j) ) * dt
                    Ex(i,j) = Ex_old(i,j) + ( -aEx(i,j)*dHy_dz + bEx(i,j)*dHy_dx - dEx(i,j)*Ex_old(i,j) - eEx(i,j)*Ez_old(i,j) ) * dt
                
                enddo
            enddo
            ! do j = 2, nz                          ! needs j-1
            !     do i = 1, nx
            !         dHy_dz = (Hy(i,j) - Hy(i,j-1)) / dz
            !         memory_dHy_dz(i,j) = bcoef(i,j) * memory_dHy_dz(i,j) + acoef(i,j) * dHy_dz    ! FULL z
            !         dHy_dz = dHy_dz / kappa(i,j) + memory_dHy_dz(i,j)

            !         ! coupled Ex update (explicit σ-terms as in your code)
            !         Ex(i,j) = Ex_old(i,j) + ( -aEx(i,j)*dHy_dz + bEx(i,j)*0d0 - dEx(i,j)*Ex_old(i,j) - eEx(i,j)*Ez_old(i,j) ) * dt
            !     end do
            ! end do

            ! --- E_z update (uses ∂Hy/∂x → HALF x CPML) ---
            ! do j = 1, nz
            !     do i = 2, nx                        ! needs i-1
            !         dHy_dx = (Hy(i,j) - Hy(i-1,j)) / dx
            !         memory_dHy_dx(i,j) = bcoef_half(i,j) * memory_dHy_dx(i,j) + acoef_half(i,j) * dHy_dx   ! HALF x
            !         dHy_dx = dHy_dx / kappa_half(i,j) + memory_dHy_dx(i,j)

            !         Ez(i,j) = Ez_old(i,j) + ( -aEz(i,j)*0d0 + bEz(i,j)*dHy_dx - dEz(i,j)*Ex_old(i,j) - eEz(i,j)*Ez_old(i,j) ) * dt
            !     end do
            ! end do
            
            ! do j = 2,nz
            !     do i = 1,nx
            !         ! Update the Ex field
            !         dHy_dz = ( Hy(i,j) - Hy(i,j-1) )/dz ! this is nz-1 length vector
            !         memory_dHy_dz(i,j) = bcoef(i,j) * memory_dHy_dz(i,j) + acoef(i,j) * dHy_dz
            !         dHy_dz = dHy_dz/kappa(i,j) + memory_dHy_dz(i,j)

            !         ! Ex(i,j) = (( caEx(i,j) + caEx(i,j-1) )/2) * Ex(i,j) + &
            !         !     (( cbEx(i,j) + cbEx(i,j-1) )/2 ) * value_dHy_dz
            !         Ex(i,j) = caEx(i,j) * Ex(i,j) + cbEx(i,j) * dHy_dz
            !     enddo
            ! enddo

            ! do j = 1,nz
            !     do i = 2,nx
            !         ! Update the Ez field
            !         dHy_dx = ( Hy(i,j) - Hy(i-1,j) )/dx
            !         memory_dHy_dx(i,j) = bcoef_half(i,j) * memory_dHy_dx(i,j) + acoef_half(i,j) * dHy_dx
            !         dHy_dx = dHy_dx/kappa_half(i,j) + memory_dHy_dx(i,j)
                    
            !         ! Ez(i,j) = (( caEz(i,j) + caEz(i-1,j) )/2) * Ez(i,j) + &
            !         !     (( cbEz(i,j) + cbEz(i-1,j) )/2) * value_dHy_dx 
            !         Ez(i,j) = caEz(i,j) * Ez(i,j) + cbEz(i,j) * dHy_dx 
            !     enddo
            ! enddo
            !$omp end do 
            
            !$omp end parallel
            !----------------------------------------------------------------------------
            Ex(isource,jsource) = Ex(isource,jsource) + &
                            srcx(it) * dt / epsilon11(isource,jsource)
            Ez(isource,jsource) = Ez(isource,jsource) + &
                            srcz(it) * dt / epsilon33(isource,jsource) 
            
            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            Ex(1,:) = 0.0_real64
            Ex(nx,:) = 0.0_real64
            Ex(:,1) = 0.0_real64
            Ex(:,nz) = 0.0_real64

            Ez(1,:) = 0.0_real64
            Ez(nx,:) = 0.0_real64
            Ez(:,1) = 0.0_real64
            Ez(:,nz) = 0.0_real64

            Hy(1,:) = 0.0_real64
            Hy(nx,:) = 0.0_real64
            Hy(:,1) = 0.0_real64
            Hy(:,nz) = 0.0_real64
            
            Ex_old = Ex 
            Ez_old = Ez
            ! print maximum of norm of velocity
            velocnorm = maxval(sqrt(Ex**2 + Ez**2))
            if (velocnorm > stability_threshold) stop 'code became unstable and blew up'

            call write_image2(Ex, nx, nz, source, it, 'Ex', SINGLE)
            call write_image2(Ez, nx, nz, source, it, 'Ez', SINGLE)
        enddo
        
        
        deallocate(eps11, eps13,  eps33, sig11, sig13,  sig33, srcx, srcz)
        deallocate(kappa, alpha, acoef, bcoef, kappa_half, alpha_half, acoef_half, bcoef_half)
        deallocate(Ex, Ez, Hy)
        deallocate(memory_dEz_dx, memory_dEx_dz, memory_dHy_dx, memory_dHy_dz)

        deallocate(aEx, bEx, dEx, eEx, aEz, bEz, dEz, eEz, det)
        deallocate(Ex_old, Ez_old )
            
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
        ! real(real64), allocatable :: epsilonx(:,:), epsilony(:,:), epsilonz(:,:), &
                                        ! sigmax(:,:), sigmay(:,:), sigmaz(:,:)

        ! real(real64) :: DT
        real(real64) :: velocnorm
        integer :: isource, jsource, ksource, i, j, k, it

        ! Coefficients for the finite difference scheme
        ! real(real64), allocatable :: caEx(:,:), cbEx(:,:), caEy(:,:), cbEy(:,:), caEz(:,:), cbEz(:,:)
        real(real64), allocatable :: aEx(:,:), bEx(:,:), cEx(:,:), &
                                     aEy(:,:), bEy(:,:), cEy(:,:), &
                                     aEz(:,:), bEz(:,:), cEz(:,:) 
        real(real64), allocatable :: det(:,:) 
        real(real64) :: daHx, dbHx, daHy, dbHy, daHz, dbHz

        real(real64) :: dEx_dy, dEy_dx, dEy_dz, dEz_dy, dEz_dx, dEx_dz, &
                        dHx_dy, dHx_dz, dHy_dx, dHy_dz, dHz_dy, dHz_dx
        
        real(real64) :: rhs_x, rhs_y, rhs_z 

        ! Source arrays
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:)

        ! 1D arrays for the damping profiles in each direction
        real(real64), allocatable :: kappa(:,:), alpha(:,:), acoef(:,:), bcoef(:,:), & 
                        kappa_half(:,:), alpha_half(:,:), acoef_half(:,:), bcoef_half(:,:)

        
        real(real64), allocatable :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:), &
                                    Hx(:,:,:), Hy(:,:,:), Hz(:,:,:) 
        real(real64), allocatable :: Ex_old(:,:,:), Ey_old(:,:,:), Ez_old(:,:,:)
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
        
        allocate(alpha(nx,nz), kappa(nx,nz), acoef(nx,nz), bcoef(nx,nz) )         
        allocate(alpha_half(nx,nz), kappa_half(nx,nz), acoef_half(nx,nz), bcoef_half(nx,nz) )

        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps))
        
        ! Allocate more
        allocate(aEx(nx,nz), bEx(nx,nz), cEx(nx,nz), &
                 aEy(nx,nz), bEy(nx,nz), cEy(nx,nz), &
                 aEz(nx,nz), bEz(nx,nz), cEz(nx,nz))
        allocate(det(nx,nz))
        allocate( memory_dEx_dy(nx,ny,nz), memory_dEy_dx(nx,ny,nz), &
                    memory_dEx_dz(nx,ny,nz), memory_dEz_dx(nx,ny,nz), &
                    memory_dEz_dy(nx,ny,nz), memory_dEy_dz(nx,ny,nz) )
        allocate(memory_dHz_dx(nx,ny,nz), memory_dHx_dz(nx,ny,nz), & 
                    memory_dHz_dy(nx,ny,nz), memory_dHy_dz(nx,ny,nz), &
                    memory_dHx_dy(nx,ny,nz), memory_dHy_dx(nx,ny,nz) )
        
        allocate(Ex(nx, ny, nz), Ey(nx, ny, nz), Ez(nx, ny, nz))
        allocate(Ex_old(nx, ny, nz), Ey_old(nx, ny, nz), Ez_old(nx, ny, nz))
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
        dbHx = dt/mu0 !dt/(mu0*mu ) 
        daHx = 1.0_real64 ! This seems useless, but this saves room for a permeability value that isn't assumed to be unity.

        dbHy = dt/mu0 !dt/(mu0*mu ) 
        daHy = 1.0_real64 ! 

        dbHz = dt/mu0 !dt/(mu0*mu ) 
        daHz = 1.0_real64 


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
        kappa(:,:) = 1.0_real64
        kappa_half(:,:) = 1.0_real64
        alpha(:,:) = 0.0_real64
        alpha_half(:,:) = 0.0_real64
        acoef(:,:) = 0.0_real64
        acoef_half(:,:) = 0.0_real64
        bcoef(:,:) = 0.0_real64 
        bcoef_half(:,:) = 0.0_real64 


        ! ------------------------------ Load the boundary ----------------------------
        call material_rw2('kappa_cpml.dat', kappa, .TRUE.)
        call material_rw2('alpha_cpml.dat', alpha, .TRUE.)
        call material_rw2('acoef_cpml.dat', acoef, .TRUE.)
        call material_rw2('bcoef_cpml.dat', bcoef, .TRUE.)


        call material_rw2('kappa_half_cpml.dat', kappa_half, .TRUE.)
        call material_rw2('alpha_half_cpml.dat', alpha_half, .TRUE.)
        call material_rw2('acoef_half_cpml.dat', acoef_half, .TRUE.)
        call material_rw2('bcoef_half_cpml.dat', bcoef_half, .TRUE.)

        ! -----------------------------------------------------------------------------
        ! Load initial conditions
        call material_rw3('initialconditionEx.dat', Ex, .TRUE.)
        call material_rw3('initialconditionEy.dat', Ey, .TRUE.)
        call material_rw3('initialconditionEz.dat', Ez, .TRUE.)
        
        Hx(:,:,:) = 0.0_real64  
        Hy(:,:,:) = 0.0_real64 
        Hz(:,:,:) = 0.0_real64 
        
        ! PML
        memory_dEx_dy(:,:,:) = 0.0_real64
        memory_dEy_dx(:,:,:) = 0.0_real64
        memory_dEx_dz(:,:,:) = 0.0_real64
        memory_dEz_dx(:,:,:) = 0.0_real64
        memory_dEz_dy(:,:,:) = 0.0_real64
        memory_dEy_dz(:,:,:) = 0.0_real64

        memory_dHz_dx(:,:,:) = 0.0_real64
        memory_dHx_dz(:,:,:) = 0.0_real64
        memory_dHz_dy(:,:,:) = 0.0_real64
        memory_dHy_dz(:,:,:) = 0.0_real64
        memory_dHx_dy(:,:,:) = 0.0_real64
        memory_dHy_dx(:,:,:) = 0.0_real64
        
        ! Scale the permittivity in terms of relative permittivity 
        eps11 = eps11 * eps0
        eps12 = eps12 * eps0 
        eps13 = eps13 * eps0
        eps22 = eps22 * eps0 
        eps23 = eps23 * eps0 
        eps33 = eps33 * eps0
        
        det =   eps11*(eps22*eps33 - eps23*eps23) - &
                eps12*(eps12*eps33 - eps23*eps13) + & 
                eps13*(eps12*eps23 - eps22*eps13)
        
        aEx = (eps22*eps33 - eps23*eps23) / det
        bEx = - (eps12*eps33 - eps23*eps13) / det 
        cEx = (eps12*eps23 - eps22*eps13) / det 
        
        aEy = bEx  
        bEy =   (eps11*eps33 - eps13*eps13) / det
        cEy = - (eps11*eps23 - eps12*eps13) / det 
        
        aEz = cEx
        bEz = cEy
        cEz = (eps11*eps22 - eps12*eps12) / det
        
        ! ---
        ! ---  beginning of time loop
        ! ---
        Ex_old = Ex 
        Ey_old = Ey 
        Ez_old = Ez
        
        do it = 1,source%time_steps
            !--------------------------------------------------------
            ! compute magnetic field and update memory variables for C-PML
            !--------------------------------------------------------
            ! Update Hx
            do k = 1,nz-1
                do i = 1,nx-1  
                    do j = 1,ny-1
                        ! Values needed for the magnetic field updates
                        dEz_dy = ( Ez(i,j+1,k) - Ez(i,j,k) )/dy
                        dEy_dz = ( Ey(i,j,k+1) - Ey(i,j,k) )/dz
            
                        ! The rest of the equation needed for agnetic field updates
                        memory_dEy_dz(i,j,k) = bcoef_half(i,k) * memory_dEy_dz(i,j,k) + acoef_half(i,k) * dEy_dz
                        memory_dEz_dy(i,j,k) = bcoef_half(i,k) * memory_dEz_dy(i,j,k) + acoef_half(i,k) * dEz_dy
                        
                        dEz_dy = dEz_dy/ kappa_half(i,k) + memory_dEz_dy(i,j,k)
                        dEy_dz = dEy_dz/ kappa_half(i,k) + memory_dEy_dz(i,j,k)

                        ! Now update the Magnetic field
                        Hx(i,j,k) = daHx*Hx(i,j,k) + dbHx*( dEy_dz - dEz_dy )
                    enddo
                enddo  
            enddo

                ! Update Hy
            do k = 1,nz-1
                do i = 1,nx-1      
                    do j = 1,ny-1
                    
                        ! Values needed for the magnetic field updates
                        dEx_dz = ( Ex(i,j,k+1) - Ex(i,j,k) )/dz
                        dEz_dx = ( Ez(i+1,j,k) - Ez(i,j,k) )/dx
                        
                        memory_dEx_dz(i,j,k) = bcoef(i,k) * memory_dEx_dz(i,j,k) + acoef(i,k) * dEx_dz
                        memory_dEz_dx(i,j,k) = bcoef(i,k) * memory_dEz_dx(i,j,k) + acoef(i,k) * dEz_dx
                        
                        dEx_dz = dEx_dz/ kappa(i,k) + memory_dEx_dz(i,j,k)
                        dEz_dx = dEz_dx/ kappa(i,k) + memory_dEz_dx(i,j,k)
                        
                        ! Now update the Magnetic field
                        Hy(i,j,k) = daHy*Hy(i,j,k) + dbHy*( dEz_dx - dEx_dz )

                    enddo
                enddo  
            enddo

                ! Update Hz
            do k = 2,nz-1
                do i = 1,nx-1      
                    do j = 1,ny-1
                        ! Values needed for the magnetic field updates
                        dEx_dy = ( Ex(i,j+1,k) - Ex(i,j,k) )/dy
                        dEy_dx = ( Ey(i+1,j,k) - Ey(i,j,k) )/dx
                        
                        memory_dEx_dy(i,j,k) = bcoef(i,k) * memory_dEx_dy(i,j,k) + acoef(i,k) * dEx_dy
                        memory_dEy_dx(i,j,k) = bcoef(i,k) * memory_dEy_dx(i,j,k) + acoef(i,k) * dEy_dx
                        
                        dEx_dy = dEx_dy/ kappa(i,k) + memory_dEx_dy(i,j,k)
                        dEy_dx = dEy_dx/ kappa(i,k) + memory_dEy_dx(i,j,k)

                        ! Now update the Magnetic field
                        Hz(i,j,k) = daHz*Hz(i,j,k) + dbHz*( dEx_dy - dEy_dx )
                    enddo
                enddo  
            enddo

            !--------------------------------------------------------
            ! compute electric field and update memory variables for C-PML
            !--------------------------------------------------------
            do k =2,nz-1 
                do i = 2,nx-1
                    do j = 2,ny-1 
                        dHz_dy = ( Hz(i,j,k) - Hz(i,j-1,k) )/dy
                        dHy_dz = ( Hy(i,j,k) - Hy(i,j,k-1) )/dz
                        dHz_dx = ( Hz(i,j,k) - Hz(i-1,j,k) )/dx
                        dHx_dz = ( Hx(i,j,k) - Hx(i,j,k-1) )/dz
                        dHy_dx = ( Hy(i,j,k) - Hy(i-1,j,k) )/dx
                        dHx_dy = ( Hx(i,j,k) - Hx(i,j-1,k) )/dy
                        
                        memory_dHz_dy(i,j,k) = bcoef_half(i,k) * memory_dHz_dy(i,j,k) + acoef_half(i,k) * dHz_dy
                        memory_dHy_dz(i,j,k) = bcoef(i,k) * memory_dHy_dz(i,j,k) + acoef(i,k) * dHy_dz
                        memory_dHz_dx(i,j,k) = bcoef_half(i,k) * memory_dHz_dx(i,j,k) + acoef_half(i,k) * dHz_dx
                        memory_dHx_dz(i,j,k) = bcoef_half(i,k) * memory_dHx_dz(i,j,k) + acoef_half(i,k) * dHx_dz
                        memory_dHx_dy(i,j,k) = bcoef_half(i,k) * memory_dHx_dy(i,j,k) + acoef_half(i,k) * dHx_dy
                        memory_dHy_dx(i,j,k) = bcoef_half(i,k) * memory_dHy_dx(i,j,k) + acoef_half(i,k) * dHy_dx
                        
                        dHz_dy = dHz_dy/kappa_half(i,k) + memory_dHz_dy(i,j,k)
                        dHy_dz = dHy_dz/kappa(i,k) + memory_dHy_dz(i,j,k)
                        dHz_dx = dHz_dx/kappa_half(i,k) + memory_dHz_dx(i,j,k)
                        dHx_dz = dHx_dz/kappa_half(i,k) + memory_dHx_dz(i,j,k)
                        dHx_dy = dHx_dy/kappa_half(i,k) + memory_dHx_dy(i,j,k)
                        dHy_dx = dHy_dx/kappa_half(i,k) + memory_dHy_dx(i,j,k)
                        
                        rhs_x = dHz_dy - dHy_dz - (sig11(i,k) * Ex_old(i,j,k) + sig12(i,k)*Ey_old(i,j,k) + sig13(i,k) * Ez_old(i,j,k))
                        rhs_y = dHx_dz - dHz_dx - (sig12(i,k) * Ex_old(i,j,k) + sig22(i,k)*Ey_old(i,j,k) + sig23(i,k) * Ez_old(i,j,k))
                        rhs_z = dHy_dx - dHx_dy - (sig13(i,k) * Ex_old(i,j,k) + sig23(i,k)*Ey_old(i,j,k) + sig33(i,k) * Ez_old(i,j,k))
                        
                        Ex(i,j,k) = Ex_old(i,j,k) + ( aEx(i,k) * rhs_x + bEx(i,k) * rhs_y + cEx(i,k) * rhs_z ) * dt
                        Ey(i,j,k) = Ey_old(i,j,k) + ( aEy(i,k) * rhs_x + bEy(i,k) * rhs_y + cEy(i,k) * rhs_z ) * dt
                        Ez(i,j,k) = Ez_old(i,j,k) + ( aEz(i,k) * rhs_x + bEz(i,k) * rhs_y + cEz(i,k) * rhs_z ) * dt
                    end do 
                end do 
            end do 

            ! add the source (force vector located at a given grid point)
            Ex(isource,jsource,ksource) = Ex(isource,jsource,ksource) + & 
                        srcx(it) * dt / eps11(isource,ksource)
            Ey(isource,jsource,ksource) = Ey(isource,jsource,ksource) + & 
                        srcy(it) * dt / eps22(isource,ksource) 
            Ez(isource,jsource,ksource) = Ez(isource,jsource,ksource) + & 
                        srcz(it) * dt / eps33(isource,ksource)
            
            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            Ex(1,:,:) = 0.0_real64
            Ex(:,1,:) = 0.0_real64
            Ex(:,:,1) = 0.0_real64
            Ex(nx,:,:) = 0.0_real64
            Ex(:,ny,:) = 0.0_real64
            Ex(:,:,nz) = 0.0_real64 

            Ey(1,:,:) = 0.0_real64
            Ey(:,1,:) = 0.0_real64
            Ey(:,:,1) = 0.0_real64
            Ey(nx,:,:) = 0.0_real64
            Ey(:,ny,:) = 0.0_real64
            Ey(:,:,nz) = 0.0_real64
            
            Ez(1,:,:) = 0.0_real64
            Ez(:,1,:) = 0.0_real64
            Ez(:,:,1) = 0.0_real64
            Ez(nx,:,:) = 0.0_real64
            Ez(:,ny,:) = 0.0_real64
            Ez(:,:,nz) = 0.0_real64
            
            Hx(1,:,:) = 0.0_real64
            Hx(:,1,:) = 0.0_real64
            Hx(:,:,1) = 0.0_real64
            Hx(nx,:,:) = 0.0_real64
            Hx(:,ny,:) = 0.0_real64
            Hx(:,:,nz) = 0.0_real64

            Hy(1,:,:) = 0.0_real64
            Hy(:,1,:) = 0.0_real64
            Hy(:,:,1) = 0.0_real64
            Hy(nx,:,:) = 0.0_real64
            Hy(:,ny,:) = 0.0_real64
            Hy(:,:,nz) = 0.0_real64
            
            Hz(1,:,:) = 0.0_real64
            Hz(:,1,:) = 0.0_real64
            Hz(:,:,1) = 0.0_real64
            Hz(nx,:,:) = 0.0_real64
            Hz(:,ny,:) = 0.0_real64
            Hz(:,:,nz) = 0.0_real64
            
            Ex_old = Ex 
            Ey_old = Ey 
            Ez_old = Ez 
            
            ! check norm of velocity to make sure the solution isn't diverging
            velocnorm = maxval(sqrt(Ex**2.0d0 + Ey**2.0d0 + Ez**2.0d0) )
            if (velocnorm > stability_threshold) stop 'code became unstable and blew up'
            ! print *,'Max vals for Ex, Ey, Ez: ', maxval(Ex), maxval(Ey), maxval(Ez)

            ! print *, maxval(Ex), maxval(Ey), maxval(Ez)
            call write_image(Ex, domain, source, it, 'Ex', SINGLE)
            call write_image(Ey, domain, source, it, 'Ey', SINGLE)
            call write_image(Ez, domain, source, it, 'Ez', SINGLE)

        enddo   ! end of time loop
    
        deallocate( eps11, eps12, eps13, eps22, eps23, eps33)
        deallocate( sig11, sig12, sig13, sig22, sig23, sig33)
        deallocate( kappa, alpha, acoef, bcoef, kappa_half)
        deallocate( alpha_half, acoef_half, bcoef_half)
        deallocate( srcx, srcy, srcz)
        deallocate( aEx, bEx, cEx, aEy, bEy, cEy, aEz, bEz, cEz, det)
        deallocate( memory_dEx_dy, memory_dEy_dx, &
                    memory_dEx_dz, memory_dEz_dx, &
                    memory_dEz_dy, memory_dEy_dz )
        deallocate( memory_dHz_dx, memory_dHx_dz, & 
                    memory_dHz_dy, memory_dHy_dz, &
                    memory_dHx_dy, memory_dHy_dx )
        
        deallocate(Ex, Ey, Ez, Hx, Hy, Hz, Ex_old, Ey_old, Ez_old)

    end subroutine electromag25
    

    
end module cpmlfdtd