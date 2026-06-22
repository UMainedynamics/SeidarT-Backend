

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
            case ('harmonic'  ); density_code = 1
            case ('geometric'); density_code = 2
            case ('arithmetic'); density_code = 3
            case ('none'      ); density_code = 4
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
                    rhoxx = scalar_mean(rho(i,k), rho(i-1,k), density_code)
                    rhozx = scalar_mean(rho(i,k), rho(i,k-1), density_code) 
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
                    rhoxz = scalar_mean(rho(i,k), rho(i+1,k), density_code)
                    rhozz = scalar_mean(rho(i,k), rho(i,k+1), density_code)
                    
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
                        dt 

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
            case ('harmonic'  ); density_code = 1
            case ('geometric'); density_code = 2
            case ('arithmetic'); density_code = 3
            case ('none'      ); density_code = 4
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
                        rhoxx = scalar_mean(rho(i,k), rho(i-1,k), density_code)
                        rhoyx = scalar_mean(rho(i,k), rho(i,k), density_code) 
                        rhozx = scalar_mean(rho(i,k), rho(i,k-1), density_code) 
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
                        rhoxy = scalar_mean(rho(i,k), rho(i+1,k), density_code)
                        rhoyy = scalar_mean(rho(i,k), rho(i,k), density_code) 
                        rhozy = scalar_mean(rho(i,k), rho(i,k-1), density_code)
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
                        rhoxz = scalar_mean(rho(i,k), rho(i+1,k), density_code)
                        rhoyz = scalar_mean(rho(i,k), rho(i,k), density_code) 
                        rhozz = scalar_mean(rho(i,k), rho(i,k+1), density_code)
                        
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
                            dt 

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
    subroutine seismic3iso(domain, source, density_method, SINGLE_OUTPUT)
        !--------------------------------------------------------------------------------------
        use constants
        use omp_lib
        
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

        real(real64), allocatable :: c11(:,:,:), c12(:,:,:), c13(:,:,:), &
                                     c22(:,:,:), c23(:,:,:), c33(:,:,:), &
                                     c44(:,:,:), c55(:,:,:), c66(:,:,:), rho(:,:,:)
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
        real(real64), allocatable :: kappa(:,:,:), alpha(:,:,:), acoef(:,:,:), bcoef(:,:,:), &
                        kappa_half(:,:,:), alpha_half(:,:,:), acoef_half(:,:,:), bcoef_half(:,:,:)
        
        ! Arrays for the PML damping factors
        real(real64), allocatable :: gamma_x(:,:,:), gamma_y(:,:,:), gamma_z(:,:,:)
        real(real64), allocatable :: gamma_xz(:,:,:), gamma_xy(:,:,:), gamma_yz(:,:,:)

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
            case ('harmonic'  ); density_code = 1
            case ('geometric'); density_code = 2
            case ('arithmetic'); density_code = 3
            case ('none'      ); density_code = 4
        end select
        
        nx = domain%nx
        ny = domain%ny 
        nz = domain%nz 
        dt = source%dt 
        dx = domain%dx
        dy = domain%dy 
        dz = domain%dz
        
        ! ----------------------- Allocate Arrays ----------------------------
        allocate(c11(nx, ny, nz), c12(nx, ny, nz), c13(nx, ny, nz), &
                 c22(nx, ny, nz), c23(nx, ny, nz), c33(nx, ny, nz), &
                 c44(nx, ny, nz), c55(nx, ny, nz), c66(nx, ny, nz) )
        allocate(rho(nx,ny,nz))
                
        allocate(kappa(nx, ny, nz), alpha(nx, ny, nz), acoef(nx, ny, nz), bcoef(nx, ny, nz), &
                kappa_half(nx, ny, nz), alpha_half(nx, ny, nz), acoef_half(nx, ny, nz), bcoef_half(nx, ny, nz))
        allocate(gamma_x(nx, ny, nz), gamma_y(nx, ny, nz), gamma_z(nx, ny, nz))
        allocate(gamma_xy(nx, ny, nz), gamma_yz(nx, ny, nz), gamma_xz(nx, ny, nz))
        
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
            call material_rw3('c11.dat', c11, .TRUE.)
            call material_rw3('c12.dat', c12, .TRUE.)
            call material_rw3('c13.dat', c13, .TRUE.)
            call material_rw3('c22.dat', c22, .TRUE.)
            call material_rw3('c23.dat', c23, .TRUE.)
            call material_rw3('c33.dat', c33, .TRUE.)
            call material_rw3('c44.dat', c44, .TRUE.)
            call material_rw3('c55.dat', c55, .TRUE.)
            call material_rw3('c66.dat', c66, .TRUE.)
            call material_rw3('rho.dat', rho, .TRUE.)
        
        ! ------------------- Load Attenuation Coefficients --------------------
            call material_rw3('gamma_x.dat', gamma_x, .TRUE.)
            call material_rw3('gamma_z.dat', gamma_z, .TRUE.)
            call material_rw3('gamma_y.dat', gamma_y, .TRUE.)
            call material_rw3('gamma_xz.dat', gamma_xz, .TRUE.)
            call material_rw3('gamma_yz.dat', gamma_yz, .TRUE.)
            call material_rw3('gamma_xy.dat', gamma_xy, .TRUE.)
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
        kappa(:,:,:) = 1.0_real64
        kappa_half(:,:,:) = 1.0_real64
        alpha(:,:,:) = 0.0_real64
        alpha_half(:,:,:) = 0.0_real64
        acoef(:,:,:) = 0.0_real64
        acoef_half(:,:,:) = 0.0_real64

        ! ------------------------- Boundary Conditions -------------------------
        call material_rw3('kappa_cpml.dat', kappa, .TRUE.)
        call material_rw3('alpha_cpml.dat', alpha, .TRUE.)
        call material_rw3('acoef_cpml.dat', acoef, .TRUE.)
        call material_rw3('bcoef_cpml.dat', bcoef, .TRUE.)

        call material_rw3('kappa_half_cpml.dat', kappa_half, .TRUE.)
        call material_rw3('alpha_half_cpml.dat', alpha_half, .TRUE.)
        call material_rw3('acoef_half_cpml.dat', acoef_half, .TRUE.)
        call material_rw3('bcoef_half_cpml.dat', bcoef_half, .TRUE.)

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
            ! Loop 1: update sigmaxx, sigmayy, sigmazz
#ifndef SEIDART_OPENMP_GPU
            !!$omp parallel do private(i, j, k, dvx_dx, dvy_dy, dvz_dz) &
            !!$omp& schedule(static)
#endif
            do k = 3,nz-1
                do j = 3,ny-1
                    do i = 2,nx-2
                        ! Forward step on half grid
                        dvx_dx = ( 27.0_real64*(vx(i+1,j,k) - vx(i,j,k)) - (vx(i+2,j,k) - vx(i-1,j,k))) / (24.0_real64*dx)
                        ! Backward step on full grid
                        dvy_dy = ( 27.0_real64*(vy(i,j,k) - vy(i,j-1,k)) - (vy(i,j+1,k) - vy(i,j-2,k))) / (24.0_real64*dy)
                        ! Backward step on full grid
                        dvz_dz = ( 27.0_real64*(vz(i,j,k) - vz(i,j,k-1)) - (vz(i,j,k+1) - vz(i,j,k-2))) / (24.0_real64*dz)
                        
                        memory_dvx_dx(i,j,k) = bcoef_half(i,j,k) * memory_dvx_dx(i,j,k) + acoef_half(i,j,k) * dvx_dx
                        memory_dvy_dy(i,j,k) = bcoef(i,j,k) * memory_dvy_dy(i,j,k) + acoef(i,j,k) * dvy_dy
                        memory_dvz_dz(i,j,k) = bcoef(i,j,k) * memory_dvz_dz(i,j,k) + acoef(i,j,k) * dvz_dz
                        
                        dvx_dx = dvx_dx / kappa_half(i,j,k) + memory_dvx_dx(i,j,k)
                        dvy_dy = dvy_dy / kappa(i,j,k) + memory_dvy_dy(i,j,k)
                        dvz_dz = dvz_dz / kappa(i,j,k) + memory_dvz_dz(i,j,k)
                        
                        sigmaxx(i,j,k) = ( sigmaxx(i,j,k) + &
                        ( c11(i,j,k) * dvx_dx + c12(i,j,k) * dvy_dy + c13(i,j,k) * dvz_dz ) * dt ) / &
                                (1 + gamma_x(i,j,k) * dt )

                        ! Full 3D will need a gradient in the y-direction
                        sigmayy(i,j,k) = ( sigmayy(i,j,k) + &
                        ( c12(i,j,k) * dvx_dx + c22(i,j,k) * dvy_dy + c23(i,j,k) * dvz_dz ) * dt ) / &
                                (1 + gamma_y(i,j,k) * dt )

                        sigmazz(i,j,k) = ( sigmazz(i,j,k) + &
                        ( c13(i,j,k) * dvx_dx + c23(i,j,k) * dvy_dy + c33(i,j,k) * dvz_dz ) * dt ) / &
                                (1 + gamma_z(i,j,k) * dt )

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
                        
                        memory_dvy_dx(i,j,k) = bcoef(i,j,k) * memory_dvy_dx(i,j,k) + acoef(i,j,k) * dvy_dx
                        memory_dvx_dy(i,j,k) = bcoef_half(i,j,k) * memory_dvx_dy(i,j,k) + acoef_half(i,j,k) * dvx_dy
                        
                        dvy_dx = dvy_dx / kappa(i,j,k) + memory_dvy_dx(i,j,k)
                        dvx_dy = dvx_dy / kappa_half(i,j,k) + memory_dvx_dy(i,j,k)
                        
                        sigmaxy(i,j,k) = ( sigmaxy(i,j,k) + &
                        (  c66(i,j,k) * (dvy_dx + dvx_dy) ) * dt ) / &
                                (1 + gamma_xy(i,j,k) * dt )

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
                        
                        memory_dvz_dx(i,j,k) = bcoef(i,j,k) * memory_dvz_dx(i,j,k) + acoef(i,j,k) * dvz_dx
                        memory_dvx_dz(i,j,k) = bcoef_half(i,j,k) * memory_dvx_dz(i,j,k) + acoef_half(i,j,k) * dvx_dz

                        dvz_dx = dvz_dx / kappa(i,j,k) + memory_dvz_dx(i,j,k) 
                        dvx_dz = dvx_dz / kappa_half(i,j,k) + memory_dvx_dz(i,j,k)
                        
                        sigmaxz(i,j,k) = ( sigmaxz(i,j,k) + &
                            (  c55(i,j,k) * ( dvx_dz + dvz_dx) ) * dt  ) / &
                                    (1 + gamma_xz(i,j,k) * dt )

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
                        
                        memory_dvz_dy(i,j,k) = bcoef_half(i,j,k) * memory_dvz_dy(i,j,k) + acoef_half(i,j,k) * dvz_dy
                        memory_dvy_dz(i,j,k) = bcoef_half(i,j,k) * memory_dvy_dz(i,j,k) + acoef_half(i,j,k) * dvy_dz

                        dvz_dy = dvz_dy / kappa_half(i,j,k) + memory_dvz_dy(i,j,k)
                        dvy_dz = dvy_dz / kappa_half(i,j,k) + memory_dvy_dz(i,j,k)

                        sigmayz(i,j,k) = ( sigmayz(i,j,k)  + &
                            ( c44(i,j,k) * ( dvy_dz + dvz_dy) ) * dt  ) / & 
                                    (1 + gamma_yz(i,j,k) * dt )
                    enddo
                enddo
            enddo
#ifndef SEIDART_OPENMP_GPU
            !!$omp end parallel do
#endif

            !--------------------------------------------------------
            ! compute velocity and update memory variables for C-PML
            !--------------------------------------------------------
            ! Loop 5: update vx
#ifndef SEIDART_OPENMP_GPU
            !$omp parallel do collapse(2) private(i, j, k, rhoxx, rhoyx, rhozx, &
            !$omp&    dsigmaxx_dx, dsigmaxy_dy, dsigmaxz_dz) schedule(static)
#endif
            do k = 3,nz-1
                do j = 3,ny-1
                    do i = 3,nx-1
                        ! ds1/dx, ds6/dy, ds5,dz
                        rhoxx = scalar_mean(rho(i,j,k), rho(i-1,j,k), density_code)
                        rhoyx = scalar_mean(rho(i,j,k), rho(i,j,k), density_code) 
                        rhozx = scalar_mean(rho(i,j,k), rho(i,j,k-1), density_code) 
                        ! Backward difference half grid
                        dsigmaxx_dx = (27.0_real64*sigmaxx(i,j,k) - 27.0_real64*sigmaxx(i-1,j,k) + sigmaxx(i-2,j,k)) / (24.0_real64*dx)
                        dsigmaxy_dy = (27.0_real64*sigmaxy(i,j,k) - 27.0_real64*sigmaxy(i,j-1,k) + sigmaxy(i,j-2,k)) / (24.0_real64*dy)
                        dsigmaxz_dz = (27.0_real64*sigmaxz(i,j,k) - 27.0_real64*sigmaxz(i,j,k-1) + sigmaxz(i,j,k-2)) / (24.0_real64*dz)
                        
                        memory_dsigmaxx_dx(i,j,k) = bcoef_half(i,j,k) * memory_dsigmaxx_dx(i,j,k) + acoef_half(i,j,k) * dsigmaxx_dx
                        memory_dsigmaxy_dy(i,j,k) = bcoef_half(i,j,k) * memory_dsigmaxy_dy(i,j,k) + acoef_half(i,j,k) * dsigmaxy_dy
                        memory_dsigmaxz_dz(i,j,k) = bcoef_half(i,j,k) * memory_dsigmaxz_dz(i,j,k) + acoef_half(i,j,k) * dsigmaxz_dz

                        dsigmaxx_dx = dsigmaxx_dx / kappa_half(i,j,k) + memory_dsigmaxx_dx(i,j,k)
                        dsigmaxy_dy = dsigmaxy_dy / kappa_half(i,j,k) + memory_dsigmaxy_dy(i,j,k)
                        dsigmaxz_dz = dsigmaxz_dz / kappa_half(i,j,k) + memory_dsigmaxz_dz(i,j,k) 

                        vx(i,j,k) = vx(i,j,k) + &
                            (dsigmaxx_dx/rhoxx + dsigmaxy_dy/rhoyx + dsigmaxz_dz/rhozx) * dt 
                    enddo
                enddo
            enddo
#ifndef SEIDART_OPENMP_GPU
            !$omp end parallel do
#endif
            ! Loop 6: update vy
#ifndef SEIDART_OPENMP_GPU
            !!$omp parallel do collapse(2) private(i, j, k, rhoxy, rhoyy, rhozy, &
            !!$omp&    dsigmaxy_dx, dsigmayy_dy, dsigmayz_dz) schedule(static)
#endif
            do k = 3,nz-1
                do j = 2,ny-2
                    do i = 2,nx-2
                        ! ds6/dx, ds2/dy, ds4/dz
                        rhoxy = scalar_mean(rho(i,j,k), rho(i+1,j,k), density_code)
                        rhoyy = scalar_mean(rho(i,j,k), rho(i,j,k), density_code) 
                        rhozy = scalar_mean(rho(i,j,k), rho(i,j,k-1), density_code)
                        ! Forward difference half grid
                        dsigmaxy_dx = (-27.0_real64*sigmaxy(i,j,k) + 27.0_real64*sigmaxy(i+1,j,k) - sigmaxy(i+2,j,k)) / (24.0_real64*dx)
                        dsigmayy_dy = (-27.0_real64*sigmayy(i,j,k) + 27.0_real64*sigmayy(i,j+1,k) - sigmayy(i,j+2,k)) / (24.0_real64*dy)
                        ! Backward difference full grid
                        dsigmayz_dz = (27.0_real64*sigmayz(i,j,k) - 27.0_real64*sigmayz(i,j,k-1) + sigmayz(i,j,k-2)) / (24.0_real64*dz)
                        
                        memory_dsigmaxy_dx(i,j,k) = bcoef_half(i,j,k) * memory_dsigmaxy_dx(i,j,k) + acoef_half(i,j,k) * dsigmaxy_dx
                        memory_dsigmayy_dy(i,j,k) = bcoef_half(i,j,k) * memory_dsigmayy_dy(i,j,k) + acoef_half(i,j,k) * dsigmayy_dy
                        memory_dsigmayz_dz(i,j,k) = bcoef(i,j,k) * memory_dsigmayz_dz(i,j,k) + acoef(i,j,k) * dsigmayz_dz
                        
                        dsigmaxy_dx = dsigmaxy_dx / kappa_half(i,j,k) + memory_dsigmaxy_dx(i,j,k)
                        dsigmayy_dy = dsigmayy_dy / kappa_half(i,j,k) + memory_dsigmayy_dy(i,j,k)
                        dsigmayz_dz = dsigmayz_dz / kappa(i,j,k) + memory_dsigmayz_dz(i,j,k)

                        vy(i,j,k) = vy(i,j,k) + &
                            (dsigmaxy_dx/rhoxy + dsigmayy_dy/rhoyy + dsigmayz_dz/rhozy) * dt
                    enddo
                enddo
            enddo

            do k = 2,nz-2
                do j = 3,ny-1
                    do i = 2,nx-2
                        ! ds5/dx, ds4/dy, ds3/dz
                        rhoxz = scalar_mean(rho(i,j,k), rho(i+1,j,k), density_code)
                        rhoyz = scalar_mean(rho(i,j,k), rho(i,j,k), density_code) 
                        rhozz = scalar_mean(rho(i,j,k), rho(i,j,k+1), density_code)
                        
                        ! Forward difference half grid
                        dsigmaxz_dx = (-27.0_real64*sigmaxz(i,j,k) + 27.0_real64*sigmaxz(i+1,j,k) - sigmaxz(i+2,j,k)) / (24.0_real64*dx)
                        ! Backward difference full grid
                        dsigmayz_dy = (27.0_real64*sigmayz(i,j,k) - 27.0_real64*sigmayz(i,j-1,k) + sigmayz(i,j-2,k)) / (24.0_real64*dy)
                        ! Forward difference half grid
                        dsigmazz_dz = (-27.0_real64*sigmazz(i,j,k) + 27.0_real64*sigmazz(i,j,k+1) - sigmazz(i,j,k+2)) / (24.0_real64*dz)
                        
                        memory_dsigmaxz_dx(i,j,k) = bcoef_half(i,j,k) * memory_dsigmaxz_dx(i,j,k) + acoef_half(i,j,k) * dsigmaxz_dx
                        memory_dsigmayz_dy(i,j,k) = bcoef(i,j,k) * memory_dsigmayz_dy(i,j,k) + acoef(i,j,k) * dsigmayz_dy
                        memory_dsigmazz_dz(i,j,k) = bcoef_half(i,j,k) * memory_dsigmazz_dz(i,j,k) + acoef_half(i,j,k) * dsigmazz_dz

                        dsigmaxz_dx = dsigmaxz_dx / kappa_half(i,j,k) + memory_dsigmaxz_dx(i,j,k)
                        dsigmayz_dy = dsigmayz_dy / kappa(i,j,k) + memory_dsigmayz_dy(i,j,k)
                        dsigmazz_dz = dsigmazz_dz / kappa_half(i,j,k) + memory_dsigmazz_dz(i,j,k)

                        vz(i,j,k) = vz(i,j,k) + &
                            (dsigmaxz_dx/rhoxz + dsigmayz_dy/rhoyz + dsigmazz_dz/rhozz) * dt 

                    enddo
                enddo
            enddo
#ifndef SEIDART_OPENMP_GPU
            !!$omp end parallel do
#endif
            
            sigmaxx(isource,jsource,ksource) = sigmaxx(isource,jsource,ksource) + srcxx(it) / rho(isource,jsource,ksource)
            sigmaxy(isource+1,jsource+1,ksource) = sigmaxy(isource+1,jsource+1,ksource) + srcxy(it) / rho(isource+1,jsource,ksource+1)
            sigmaxz(isource+1,jsource,ksource+1) = sigmaxz(isource+1,jsource,ksource+1) + srcxz(it) / rho(isource+1,jsource,ksource)
            sigmayy(isource,jsource,ksource) = sigmayy(isource,jsource,ksource) + srcyy(it) / rho(isource,jsource,ksource)
            sigmayz(isource,jsource+1,ksource+1) = sigmayz(isource,jsource+1,ksource+1) + srcyz(it) / rho(isource,jsource,ksource+1)
            sigmazz(isource,jsource,ksource) = sigmazz(isource,jsource,ksource) + srczz(it) / rho(isource,jsource,ksource)
            vx(isource+1,jsource,ksource) = vx(isource+1,jsource,ksource) + &
                    srcx(it) * dt / rho(isource+1,jsource,ksource)
            vy(isource,jsource+1,ksource) = vy(isource,jsource+1,ksource) + &
                    srcy(it) * dt / rho(isource,jsource,ksource)
            vz(isource,jsource,ksource+1) = vz(isource,jsource,ksource+1) + &
                    srcz(it) * dt / rho(isource,jsource,ksource+1)

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
        
    end subroutine seismic3iso