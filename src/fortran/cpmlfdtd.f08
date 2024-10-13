module cpmlfdtd 
    
    use seidartio
    use seidart_types
    
    implicit none 
    ! --------------------------------------------------------------------------
    ! FD 
    ! subroutine fd2(vx, vz, sigmaxx, sigmazz, )
        
    !     do j = 2,domain%nz
    !         do i = 1,domain%nx-1
    
    !             value_dvx_dx = (vx(i+1,j) - vx(i,j)) / DX
    !             value_dvz_dz = (vz(i,j) - vz(i,j-1)) / DZ
    !             value_dvz_dx = (vz(i+1,j) - vz(i,j)) / DX
    !             value_dvx_dz = (vx(i,j) - vx(i,j-1)) / DZ

    !             memdvx_dx(i,j) = b_x_half(j) * memdvx_dx(i,j) + &
    !                                     a_x_half(i) * value_dvx_dx
    !             memdvz_dz(i,j) = b_z(j) * memdvz_dz(i,j) + &
    !                                     a_z(j) * value_dvz_dz
    !             memdvx_dz(i,j) = b_z_half(j) * memdvx_dz(i,j) + &
    !                                     a_z_half(j) * value_dvx_dz 
    !             memdvz_dx(i,j) = b_x(i) * memdvz_dx(i,j) + &
    !                                     a_x(i) * value_dvz_dx

    !             value_dvx_dx = value_dvx_dx / K_x_half(i) + memdvx_dx(i,j)
    !             value_dvz_dz = value_dvz_dz / K_z(j) + memdvz_dz(i,j)
    !             value_dvz_dx = value_dvz_dx / K_x(i) + memdvz_dx(i,j)
    !             value_dvx_dz = value_dvx_dz / K_z_half(j) + memdvx_dz(i,j)
                
    !             sigmaxx(i,j) = sigmaxx(i,j) + &
    !                 (   c11(i,j) * value_dvx_dx + &
    !                     c13(i,j) * value_dvz_dz + &
    !                     c15(i,j) * (value_dvz_dx + value_dvx_dz) ) * time_params%dt
    !             sigmazz(i,j) = sigmazz(i,j) + &
    !                 (   c13(i,j) * value_dvx_dx + &
    !                     c33(i,j) * value_dvz_dz + &
    !                     c35(i,j) * (value_dvz_dx + value_dvx_dz) ) * time_params%dt
    
    !         enddo
    !     enddo
    
    !     do j = 1,domain%nz-1
    !         do i = 2,domain%nx
    
    !         value_dvx_dx = (vx(i,j) - vx(i-1,j)) / DX
    !         value_dvz_dz = (vz(i,j+1) - vz(i,j)) / DZ
    !         value_dvz_dx = (vz(i,j) - vz(i-1,j)) / DX
    !         value_dvx_dz = (vx(i,j+1) - vx(i,j)) / DZ
            
    !         memdvx_dx(i,j) = b_x_half(i) * memdvx_dx(i,j) + &
    !                                 a_x_half(i) * value_dvx_dx
    !         memdvz_dz(i,j) = b_z(j) * memdvz_dz(i,j) + &
    !                                 a_z(j) * value_dvz_dz
    !         memdvx_dz(i,j) = b_z_half(j) * memdvx_dz(i,j) + &
    !                                 a_z_half(j) * value_dvx_dz 
    !         memdvz_dx(i,j) = b_x(i) * memdvz_dx(i,j) + &
    !                                 a_x(i) * value_dvz_dx
            
    !         value_dvx_dx = value_dvx_dx / K_x_half(i) + memdvx_dx(i,j)
    !         value_dvz_dz = value_dvz_dz / K_z(j) + memdvz_dz(i,j)
    !         value_dvz_dx = value_dvz_dx / K_x(i) + memdvz_dx(i,j)
    !         value_dvx_dz = value_dvx_dz / K_z_half(j) + memdvx_dz(i,j)
    
    !         sigmaxz(i,j) = sigmaxz(i,j) + &
    !             (   c15(i,j)  * value_dvx_dx + & 
    !                 c35(i,j)  * value_dvz_dz + &
    !                 c55(i,j) * (value_dvz_dx + value_dvx_dz) ) * time_params%dt
    
    !         enddo
    !     enddo
    
    !     !--------------------------------------------------------
    !     ! compute velocity and update memory variables for C-PML
    !     !--------------------------------------------------------
    !     do j = 2,domain%nz
    !         do i = 2,domain%nx
    
    !         deltarho = ( rho(i,j) + rho(i,j+1) + rho(i+1,j) + rho(i+1,j+1) )/4

    !         value_dsigmaxx_dx = (sigmaxx(i,j) - sigmaxx(i-1,j)) / DX
    !         value_dsigmaxz_dz = (sigmaxz(i,j) - sigmaxz(i,j-1)) / DZ
    
    !         memdsigmaxx_dx(i,j) = b_x(i) * memdsigmaxx_dx(i,j) + &
    !                     a_x(i) * value_dsigmaxx_dx
    !         memdsigmaxz_dz(i,j) = b_z(j) * memdsigmaxz_dz(i,j) + &
    !                     a_z(j) * value_dsigmaxz_dz
    
    !         value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + &
    !                     memdsigmaxx_dx(i,j)
    !         value_dsigmaxz_dz = value_dsigmaxz_dz / K_z(j) + &
    !                     memdsigmaxz_dz(i,j)
    
    !         vx(i,j) = vx(i,j)*(1 - gammax(i,j) ) + (value_dsigmaxx_dx + value_dsigmaxz_dz) * time_params%dt / rho(i,j)
    
    !         enddo
    !     enddo
    
    ! end subroutine fd2
    
    contains
    
    ! =========================================================================    
    ! Computations are done in double precision and written to binary as single
    ! precision unless specified by the optional logical, SINGLE_OUTPUT.
    subroutine seismic2(domain, time_params, source, seisvar, SINGLE_OUTPUT)
        
        use constants
        
        ! Input arguments
        type(Domain_Type), intent(in) :: domain
        type(Time_Parameters_Type), intent(in) :: time_params
        type(Source_Type), intent(in) :: source
        type(Seismic2_Variables_Type), intent(inout) :: seisvar
        logical, intent(in), optional :: SINGLE_OUTPUT

        ! Local variables
        real(real64) :: deltarho, velocnorm, value_dvx_dx, value_dvx_dz, &
            value_dvz_dx, value_dvz_dz, value_dsigmaxx_dx, value_dsigmazz_dz, &
            value_dsigmaxz_dx, value_dsigmaxz_dz

        ! 1D arrays for damping profiles
        real(real64), allocatable :: c11(:,:), c13(:,:), c15(:,:), c33(:,:), c35(:,:), c55(:,:), rho(:,:)
        real(real64), allocatable :: K_x(:), alpha_x(:), a_x(:), b_x(:), K_x_half(:), alpha_x_half(:), a_x_half(:), b_x_half(:)
        real(real64), allocatable :: K_z(:), alpha_z(:), a_z(:), b_z(:), K_z_half(:), alpha_z_half(:), a_z_half(:), b_z_half(:)

        real(real64), allocatable :: gammax(:,:), gammaz(:,:)        
        real(real64), allocatable :: srcx(:), srcz(:) ! The vector time series of the source

        ! real(real64) :: dt
        integer :: i, j, it, isource, jsource
        logical :: SINGLE

        ! The default data output is single precision unless SINGLE_OUTPUT is 
        ! set to .FALSE.
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        end if

        
        
        ! Allocate the arrays based on runtime values of domain%nx and domain%nz
        allocate(c11(domain%nx, domain%nz), c13(domain%nx, domain%nz), c15(domain%nx, domain%nz), &
                c33(domain%nx, domain%nz), c35(domain%nx, domain%nz), c55(domain%nx, domain%nz), rho(domain%nx, domain%nz))
        allocate(K_x(domain%nx), alpha_x(domain%nx), a_x(domain%nx), b_x(domain%nx), &
                K_x_half(domain%nx), alpha_x_half(domain%nx), a_x_half(domain%nx), b_x_half(domain%nx))
        allocate(K_z(domain%nz), alpha_z(domain%nz), a_z(domain%nz), b_z(domain%nz), &
                K_z_half(domain%nz), alpha_z_half(domain%nz), a_z_half(domain%nz), b_z_half(domain%nz))
        allocate(gammax(domain%nx, domain%nz), gammaz(domain%nx, domain%nz))
        allocate(srcx(time_params%time_steps), srcz(time_params%time_steps))
        
        ! -------------------- Load Stiffness Coefficients --------------------
    
        call material_rw('c11.dat', c11, .TRUE.)
        call material_rw('c13.dat', c13, .TRUE.)
        call material_rw('c15.dat', c15, .TRUE.)
        call material_rw('c33.dat', c33, .TRUE.)
        call material_rw('c35.dat', c35, .TRUE.)
        call material_rw('c55.dat', c55, .TRUE.)
        call material_rw('rho.dat', rho, .TRUE.)
        
        ! ------------------- Load Attenuation Coefficients --------------------
        call material_rw('gammax.dat', gammax, .TRUE.)
        call material_rw('gammaz.dat', gammaz, .TRUE.)
        
        ! ------------------------ Assign some constants -----------------------
    
        isource = source%xind + domain%cpml
        jsource = source%zind + domain%cpml
    
        ! DT = minval( (/dx,dz/) )/ &
            ! (sqrt( 3.d0*( maxval( (/ c11/rho, c33/rho /) ) ) ) ) 

            ! ================================ LOAD SOURCE ================================
    
        call loadsource('seismicsourcex.dat', time_params%time_steps, srcx)
        ! We are using the coordinate names x, Z but the math computes the source in 
        ! the x-z plane
        call loadsource('seismicsourcez.dat', time_params%time_steps, srcz)
    
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
        call material_rw('initialconditionVx.dat', seisvar%vx, .TRUE.)
        call material_rw('initialconditionVz.dat', seisvar%vz, .TRUE.)
        
        seisvar%sigxx(:,:) = 0.d0
        seisvar%sigzz(:,:) = 0.d0
        seisvar%sigxz(:,:) = 0.d0

        ! PML
        seisvar%memdvx_dx(:,:) = 0.d0
        seisvar%memdvx_dz(:,:) = 0.d0
        seisvar%memdvz_dx(:,:) = 0.d0
        seisvar%memdvz_dz(:,:) = 0.d0
        seisvar%memdsigxx_dx(:,:) = 0.d0
        seisvar%memdsigzz_dz(:,:) = 0.d0
        seisvar%memdsigxz_dx(:,:) = 0.d0
        seisvar%memdsigxz_dz(:,:) = 0.d0

        !---
        !---  beginning of time loop
        !---

        do it = 1,time_params%time_steps
            !------------------------------------------------------------
            ! compute stress sigma and update memory variables for C-PML
            !------------------------------------------------------------
            do j = 2,domain%nz
                do i = 1,domain%nx-1
        
                    value_dvx_dx = (seisvar%vx(i+1,j) - seisvar%vx(i,j)) / domain%dx
                    value_dvz_dz = (seisvar%vz(i,j) - seisvar%vz(i,j-1)) / domain%dz
                    value_dvz_dx = (seisvar%vz(i+1,j) - seisvar%vz(i,j)) / domain%dx
                    value_dvx_dz = (seisvar%vx(i,j) - seisvar%vx(i,j-1)) / domain%dz

                    seisvar%memdvx_dx(i,j) = b_x_half(j) * seisvar%memdvx_dx(i,j) + &
                                            a_x_half(i) * value_dvx_dx
                    seisvar%memdvz_dz(i,j) = b_z(j) * seisvar%memdvz_dz(i,j) + &
                                            a_z(j) * value_dvz_dz
                    seisvar%memdvx_dz(i,j) = b_z_half(j) * seisvar%memdvx_dz(i,j) + &
                                            a_z_half(j) * value_dvx_dz 
                    seisvar%memdvz_dx(i,j) = b_x(i) * seisvar%memdvz_dx(i,j) + &
                                            a_x(i) * value_dvz_dx

                    value_dvx_dx = value_dvx_dx / K_x_half(i) + seisvar%memdvx_dx(i,j)
                    value_dvz_dz = value_dvz_dz / K_z(j) + seisvar%memdvz_dz(i,j)
                    value_dvz_dx = value_dvz_dx / K_x(i) + seisvar%memdvz_dx(i,j)
                    value_dvx_dz = value_dvx_dz / K_z_half(j) + seisvar%memdvx_dz(i,j)
                    
                    seisvar%sigxx(i,j) = seisvar%sigxx(i,j) + &
                        (   c11(i,j) * value_dvx_dx + &
                            c13(i,j) * value_dvz_dz + &
                            c15(i,j) * (value_dvz_dx + value_dvx_dz) ) * time_params%dt
                    seisvar%sigzz(i,j) = seisvar%sigzz(i,j) + &
                        (   c13(i,j) * value_dvx_dx + &
                            c33(i,j) * value_dvz_dz + &
                            c35(i,j) * (value_dvz_dx + value_dvx_dz) ) * time_params%dt
        
                enddo
            enddo
        
            do j = 1,domain%nz-1
                do i = 2,domain%nx
        
                    value_dvx_dx = (seisvar%vx(i,j) - seisvar%vx(i-1,j)) / domain%DX
                    value_dvz_dz = (seisvar%vz(i,j+1) - seisvar%vz(i,j)) / domain%DZ
                    value_dvz_dx = (seisvar%vz(i,j) - seisvar%vz(i-1,j)) / domain%DX
                    value_dvx_dz = (seisvar%vx(i,j+1) - seisvar%vx(i,j)) / domain%DZ
                    
                    seisvar%memdvx_dx(i,j) = b_x_half(i) * seisvar%memdvx_dx(i,j) + &
                                            a_x_half(i) * value_dvx_dx
                    seisvar%memdvz_dz(i,j) = b_z(j) * seisvar%memdvz_dz(i,j) + &
                                            a_z(j) * value_dvz_dz
                    seisvar%memdvx_dz(i,j) = b_z_half(j) * seisvar%memdvx_dz(i,j) + &
                                            a_z_half(j) * value_dvx_dz 
                    seisvar%memdvz_dx(i,j) = b_x(i) * seisvar%memdvz_dx(i,j) + &
                                            a_x(i) * value_dvz_dx
                    
                    value_dvx_dx = value_dvx_dx / K_x_half(i) + seisvar%memdvx_dx(i,j)
                    value_dvz_dz = value_dvz_dz / K_z(j) + seisvar%memdvz_dz(i,j)
                    value_dvz_dx = value_dvz_dx / K_x(i) + seisvar%memdvz_dx(i,j)
                    value_dvx_dz = value_dvx_dz / K_z_half(j) + seisvar%memdvx_dz(i,j)
            
                    seisvar%sigxz(i,j) = seisvar%sigxz(i,j) + &
                    (   c15(i,j)  * value_dvx_dx + & 
                        c35(i,j)  * value_dvz_dz + &
                        c55(i,j) * (value_dvz_dx + value_dvx_dz) ) * time_params%dt
        
                enddo
            enddo
        
            !--------------------------------------------------------
            ! compute velocity and update memory variables for C-PML
            !--------------------------------------------------------
            do j = 2,domain%nz
                do i = 2,domain%nx
        
                    deltarho = ( rho(i,j) + rho(i,j+1) + rho(i+1,j) + rho(i+1,j+1) )/4

                    value_dsigmaxx_dx = (seisvar%sigxx(i,j) - seisvar%sigxx(i-1,j)) / domain%DX
                    value_dsigmaxz_dz = (seisvar%sigxz(i,j) - seisvar%sigxz(i,j-1)) / domain%DZ
            
                    seisvar%memdsigxx_dx(i,j) = b_x(i) * seisvar%memdsigxx_dx(i,j) + &
                                a_x(i) * value_dsigmaxx_dx
                    seisvar%memdsigxz_dz(i,j) = b_z(j) * seisvar%memdsigxz_dz(i,j) + &
                                a_z(j) * value_dsigmaxz_dz
            
                    value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + &
                                seisvar%memdsigxx_dx(i,j)
                    value_dsigmaxz_dz = value_dsigmaxz_dz / K_z(j) + &
                                seisvar%memdsigxz_dz(i,j)
            
                    seisvar%vx(i,j) = seisvar%vx(i,j)*(1 - gammax(i,j) ) + &
                        (value_dsigmaxx_dx + value_dsigmaxz_dz) * time_params%dt / rho(i,j)
        
                enddo
            enddo
        
            do j = 1,domain%nz-1
                do i = 1,domain%nx-1
        
                    deltarho = ( 2*rho(i,j) + rho(i+1,j) + rho(i,j+1) )/4
                    value_dsigmaxz_dx = (seisvar%sigxz(i+1,j) - seisvar%sigxz(i,j)) / domain%dx
                    value_dsigmazz_dz = (seisvar%sigzz(i,j+1) - seisvar%sigzz(i,j)) / domain%dz
            
                    seisvar%memdsigxz_dx(i,j) = b_x_half(i) * seisvar%memdsigxz_dx(i,j) + &
                                a_x_half(i) * value_dsigmaxz_dx
                    seisvar%memdsigzz_dz(i,j) = b_z_half(j) * seisvar%memdsigzz_dz(i,j) + &
                                a_z_half(j) * value_dsigmazz_dz
            
                    value_dsigmaxz_dx = value_dsigmaxz_dx / K_x_half(i) + seisvar%memdsigxz_dx(i,j)
                    value_dsigmazz_dz = value_dsigmazz_dz / K_z_half(j) + seisvar%memdsigzz_dz(i,j)
            
                    seisvar%vz(i,j) = seisvar%vz(i,j)*(1 - gammaz(i,j) ) + &
                            (value_dsigmaxz_dx + value_dsigmazz_dz) * time_params%dt / deltarho
        
                enddo
            enddo

            ! Add the source term
            seisvar%vx(isource,jsource) = seisvar%vx(isource,jsource) + srcx(it) * time_params%dt / rho(isource,jsource)
            seisvar%vz(isource,jsource) = seisvar%vz(isource,jsource) + srcz(it) * time_params%dt / rho(isource,jsource)
        
            ! Dirichlet conditions (rigid boundaries) on the edges or at the 
            ! bottom of the PML layers
            seisvar%vx(1,:) = 0.d0
            seisvar%vx(domain%nx,:) = 0.d0
        
            seisvar%vx(:,1) = 0.d0
            seisvar%vx(:,domain%nz) = 0.d0
        
            seisvar%vz(1,:) = 0.d0
            seisvar%vz(domain%nx,:) = 0.d0
        
            seisvar%vz(:,1) = 0.d0
            seisvar%vz(:,domain%nz) = 0.d0
        
            ! print maximum of norm of velocity
            velocnorm = maxval(sqrt(seisvar%vx**2 + seisvar%vz**2))
            if (velocnorm > stability_threshold) stop 'code became unstable and blew up'
        
            call write_image(seisvar%vx, domain, source, it, 'Vx', SINGLE)
            call write_image(seisvar%vz, domain, source, it, 'Vz', SINGLE)

        enddo   ! end of time loop
    end subroutine seismic2

    ! =========================================================================
    subroutine seismic25(domain, time_params, source, seisvar, SINGLE_OUTPUT)
        !--------------------------------------------------------------------------------------
        use constants
        
        implicit none

        ! Input arguments
        type(Domain_Type), intent(in) :: domain 
        type(Time_Parameters_Type), intent(in) :: time_params 
        type(Source_Type), intent(in) :: source 
        type(Seismic3_Variables_Type), intent(inout) :: seisvar
        logical, intent(in), optional :: SINGLE_OUTPUT

        ! Local variables
        real(real64), dimension(domain%nx,domain%nz) :: c11, c12, c13, c14, c15, c16, &
                                            c22, c23, c24, c25, c26, &
                                            c33, c34, c35, c36, &
                                            c44, c45, c46, &
                                            c55, c56, &
                                            c66, &
                                            rho
        real(real64) :: deltarho, velocnorm
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
        real(real64), dimension(domain%nx) :: K_x, alpha_x, a_x, b_x, K_x_half, & 
                                                alpha_x_half, a_x_half, b_x_half
        real(real64), dimension(domain%ny) :: K_y, alpha_y, a_y, b_y, K_y_half, & 
                                                alpha_y_half, a_y_half, b_y_half
        real(real64), dimension(domain%nz) :: K_z, alpha_z, a_z, b_z, K_z_half, & 
                                                alpha_z_half, a_z_half, b_z_half

        ! Arrays for the PML damping factors
        real(real64), dimension(domain%nx,domain%nz) :: gammax, gammay, gammaz

        ! Source arrays
        real(real64), dimension(time_params%time_steps) :: srcx, srcy, srcz

        ! Boolean flag to save as double precision or single precision
        logical :: SINGLE

        ! The default data output is single precision unless SINGLE_OUTPUT is 
        ! set to .FALSE.
        if (present(SINGLE_OUTPUT)) then 
            SINGLE = SINGLE_OUTPUT 
        else
            SINGLE = .TRUE.
        endif
        
        ! ------------------------ Load Stiffness Coefficients ------------------------
        call material_rw('c11.dat', c11, .TRUE.)
        call material_rw('c12.dat', c12, .TRUE.)
        call material_rw('c13.dat', c13, .TRUE.)
        call material_rw('c14.dat', c14, .TRUE.)
        call material_rw('c15.dat', c15, .TRUE.)
        call material_rw('c16.dat', c16, .TRUE.)
        call material_rw('c22.dat', c22, .TRUE.)
        call material_rw('c23.dat', c23, .TRUE.)
        call material_rw('c24.dat', c24, .TRUE.)
        call material_rw('c25.dat', c25, .TRUE.)
        call material_rw('c26.dat', c26, .TRUE.)
        call material_rw('c33.dat', c33, .TRUE.)
        call material_rw('c34.dat', c34, .TRUE.)
        call material_rw('c35.dat', c35, .TRUE.)
        call material_rw('c36.dat', c36, .TRUE.)
        call material_rw('c44.dat', c44, .TRUE.)
        call material_rw('c45.dat', c45, .TRUE.)
        call material_rw('c46.dat', c46, .TRUE.)
        call material_rw('c55.dat', c55, .TRUE.)
        call material_rw('c56.dat', c56, .TRUE.)
        call material_rw('c66.dat', c66, .TRUE.)
        call material_rw('rho.dat', rho, .TRUE.)
        
        ! ------------------- Load Attenuation Coefficients --------------------
        call material_rw('gammax.dat', gammax, .TRUE.)
        call material_rw('gammaz.dat', gammaz, .TRUE.)
        call material_rw('gammay.dat', gammay, .TRUE.)
        
        ! ------------------------ Assign some constants -----------------------
        isource = source%xind + domain%cpml
        jsource = source%yind + domain%cpml
        ksource = source%zind + domain%cpml

        ! To ensure a courant number <= 1.0, we can calculate the time step from
        ! the velocity
        ! DT = 0.7 * minval( (/domain%dx,domain%dy,domain%dz/) )/ &
        ! ( sqrt( 3.d0 * ( maxval( (/ c11/rho, c22/rho, c33/rho /) ) ) ) )

        ! ================================ LOAD SOURCE ================================

        call loadsource('seismicsourcex.dat', time_params%time_steps, srcx)
        call loadsource('seismicsourcey.dat', time_params%time_steps, srcy)
        call loadsource('seismicsourcez.dat', time_params%time_steps, srcz)

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

        ! =============================== Forward Model ===============================
        ! Load initial condition
        call material_rw('initialconditionVx.dat', seisvar%vx, .TRUE.)
        call material_rw('initialconditionVy.dat', seisvar%vy, .TRUE.)
        call material_rw('initialconditionVz.dat', seisvar%vz, .TRUE.)

        ! Initialize the stress values
        seisvar%sigxx(:,:,:) = 0.d0
        seisvar%sigyy(:,:,:) = 0.d0
        seisvar%sigzz(:,:,:) = 0.d0
        seisvar%sigxy(:,:,:) = 0.d0
        seisvar%sigxz(:,:,:) = 0.d0
        seisvar%sigyz(:,:,:) = 0.d0

        ! PML
        seisvar%memdvx_dx(:,:,:) = 0.d0
        seisvar%memdvx_dy(:,:,:) = 0.d0
        seisvar%memdvx_dz(:,:,:) = 0.d0

        seisvar%memdvy_dx(:,:,:) = 0.d0
        seisvar%memdvy_dy(:,:,:) = 0.d0
        seisvar%memdvy_dz(:,:,:) = 0.d0

        seisvar%memdvz_dx(:,:,:) = 0.d0
        seisvar%memdvz_dy(:,:,:) = 0.d0 
        seisvar%memdvz_dz(:,:,:) = 0.d0

        seisvar%memdsigxx_dx(:,:,:) = 0.d0
        seisvar%memdsigyy_dy(:,:,:) = 0.d0
        seisvar%memdsigzz_dz(:,:,:) = 0.d0

        seisvar%memdsigxy_dx(:,:,:) = 0.d0
        seisvar%memdsigxy_dy(:,:,:) = 0.d0
        seisvar%memdsigxz_dx(:,:,:) = 0.d0
        seisvar%memdsigxz_dz(:,:,:) = 0.d0
        seisvar%memdsigyz_dy(:,:,:) = 0.d0
        seisvar%memdsigyz_dz(:,:,:) = 0.d0

        ! Do it 
        do it = 1,time_params%time_steps
            !------------------------------------------------------------
            ! compute stress sigma and update memory variables for C-PML
            !------------------------------------------------------------
            ! Update in the x direction
            do k = 2,domain%nz
                do j = 2,domain%ny
                    do i = 1,domain%nx-1

                        dvx_dx = (seisvar%vx(i+1,j,k) - seisvar%vx(i,j,k) ) / domain%dx
                        dvy_dx = (seisvar%vy(i+1,j,k) - seisvar%vy(i,j,k) ) / domain%dx
                        dvz_dx = (seisvar%vz(i+1,j,k) - seisvar%vz(i,j,k) ) / domain%dx 
                        dvy_dy = (seisvar%vy(i,j,k) - seisvar%vy(i,j-1,k) ) / domain%dy
                        dvx_dy = (seisvar%vx(i,j,k) - seisvar%vx(i,j-1,k) ) / domain%dy
                        dvz_dy = (seisvar%vz(i,j,k) - seisvar%vz(i,j-1,k) ) / domain%dy
                        dvz_dz = (seisvar%vz(i,j,k) - seisvar%vz(i,j,k-1) ) / domain%dz
                        dvx_dz = (seisvar%vx(i,j,k) - seisvar%vx(i,j,k-1) ) / domain%dz
                        dvy_dz = (seisvar%vy(i,j,k) - seisvar%vy(i,j,k-1) ) / domain%dz

                        seisvar%memdvx_dx(i,j,k) = b_x_half(i) * seisvar%memdvx_dx(i,j,k) + a_x_half(i) * dvx_dx
                        seisvar%memdvy_dx(i,j,k) = b_x_half(i) * seisvar%memdvy_dx(i,j,k) + a_x_half(i) * dvy_dx
                        seisvar%memdvz_dx(i,j,k) = b_x_half(i) * seisvar%memdvz_dx(i,j,k) + a_x_half(i) * dvz_dx
                        seisvar%memdvy_dy(i,j,k) = b_y(j) * seisvar%memdvy_dy(i,j,k) + a_y(j) * dvy_dy
                        seisvar%memdvx_dy(i,j,k) = b_y(j) * seisvar%memdvx_dy(i,j,k) + a_y(j) * dvx_dy
                        seisvar%memdvz_dy(i,j,k) = b_y(j) * seisvar%memdvz_dy(i,j,k) + a_y(j) * dvz_dy
                        seisvar%memdvz_dz(i,j,k) = b_z(k) * seisvar%memdvz_dz(i,j,k) + a_z(k) * dvz_dz
                        seisvar%memdvx_dz(i,j,k) = b_z(k) * seisvar%memdvx_dz(i,j,k) + a_z(k) * dvx_dz
                        seisvar%memdvy_dz(i,j,k) = b_z(k) * seisvar%memdvy_dz(i,j,k) + a_z(k) * dvy_dz

                        dvx_dx = dvx_dx / K_x_half(i) + seisvar%memdvx_dx(i,j,k)
                        dvy_dx = dvy_dx / K_x_half(i) + seisvar%memdvy_dx(i,j,k)
                        dvz_dx = dvz_dx / K_x_half(i) + seisvar%memdvz_dx(i,j,k)
                        dvy_dy = dvy_dy / K_y(j) + seisvar%memdvy_dy(i,j,k)
                        dvx_dy = dvx_dy / K_y(j) + seisvar%memdvx_dy(i,j,k)
                        dvz_dy = dvz_dy / K_y(j) + seisvar%memdvz_dy(i,j,k)
                        dvz_dz = dvz_dz / K_z(k) + seisvar%memdvz_dz(i,j,k)
                        dvx_dz = dvx_dz / K_z(k) + seisvar%memdvx_dz(i,j,k)
                        dvy_dz = dvy_dz / K_z(k) + seisvar%memdvy_dz(i,j,k)

                        seisvar%sigxx(i,j,k) = seisvar%sigxx(i,j,k) + &
                        (   c11(i,k) * dvx_dx + c12(i,k) * dvy_dy + c13(i,k) * dvz_dz + &
                            c14(i,k) * (dvy_dz + dvz_dy) + c15(i,k) * (dvx_dz + dvz_dx) + &
                            c16(i,k) * (dvx_dy + dvz_dy) ) * time_params%dt

                        ! Full 3D will need a gradient in the y-direction
                        seisvar%sigyy(i,j,k) = seisvar%sigyy(i,j,k) + &
                        (   c12(i,k) * dvx_dx + c22(i,k) * dvy_dy + c23(i,k) * dvz_dz + &
                            c24(i,k) * (dvy_dz + dvz_dy) + c25(i,k) * (dvx_dz + dvz_dx) + &
                            c26(i,k) * (dvy_dx + dvx_dy) ) * time_params%dt

                        seisvar%sigzz(i,j,k) = seisvar%sigzz(i,j,k) + &
                        (   c13(i,k) * dvx_dx + c23(i,k) * dvy_dy + c33(i,k) * dvz_dz + &
                            c34(i,k) * (dvy_dz + dvz_dy) + c35(i,k) * (dvx_dz + dvz_dx) + &
                            c36(i,k) * (dvy_dx + dvx_dy) ) * time_params%dt

                    enddo
                enddo
            enddo

            ! Update sigmaxy, x-direction is full nodes
            do k = 2,domain%nz
                do j = 1,domain%ny-1
                    do i = 2,domain%nx

                        dvx_dx = (seisvar%vx(i,j,k) - seisvar%vx(i-1,j,k)) / domain%dx
                        dvy_dx = (seisvar%vy(i,j,k) - seisvar%vy(i-1,j,k)) / domain%dx
                        dvz_dx = (seisvar%vz(i,j,k) - seisvar%vz(i-1,j,k)) / domain%dx
                        dvy_dy = (seisvar%vy(i,j+1,k) - seisvar%vy(i,j,k)) / domain%dy
                        dvx_dy = (seisvar%vx(i,j+1,k) - seisvar%vx(i,j,k)) / domain%dy
                        dvz_dy = (seisvar%vz(i,j+1,k) - seisvar%vz(i,j,k)) / domain%dy
                        dvz_dz = (seisvar%vz(i,j,k) - seisvar%vz(i,j,k-1)) / domain%dz
                        dvx_dz = (seisvar%vx(i,j,k) - seisvar%vx(i,j,k-1)) / domain%dz
                        dvy_dz = (seisvar%vy(i,j,k) - seisvar%vy(i,j,k-1)) / domain%dz

                        seisvar%memdvx_dx(i,j,k) = b_x(i) * seisvar%memdvx_dx(i,j,k) + a_x(i) * dvx_dx
                        seisvar%memdvy_dx(i,j,k) = b_x(i) * seisvar%memdvy_dx(i,j,k) + a_x(i) * dvy_dx
                        seisvar%memdvz_dx(i,j,k) = b_x(i) * seisvar%memdvz_dx(i,j,k) + a_x(i) * dvz_dx
                        seisvar%memdvy_dy(i,j,k) = b_y_half(j) * seisvar%memdvy_dy(i,j,k) + a_y_half(j) * dvy_dy
                        seisvar%memdvx_dy(i,j,k) = b_y_half(j) * seisvar%memdvx_dy(i,j,k) + a_y_half(j) * dvx_dy
                        seisvar%memdvz_dy(i,j,k) = b_y_half(j) * seisvar%memdvz_dy(i,j,k) + a_y_half(j) * dvz_dy
                        seisvar%memdvz_dz(i,j,k) = b_z_half(k) * seisvar%memdvz_dz(i,j,k) + a_z_half(k) * dvz_dz
                        seisvar%memdvx_dz(i,j,k) = b_z_half(k) * seisvar%memdvx_dz(i,j,k) + a_z_half(k) * dvx_dz
                        seisvar%memdvy_dz(i,j,k) = b_z_half(k) * seisvar%memdvy_dz(i,j,k) + a_z_half(k) * dvy_dz

                        dvx_dx = dvx_dx / K_x(i) + seisvar%memdvx_dx(i,j,k)
                        dvy_dx = dvy_dx / K_x(i) + seisvar%memdvy_dx(i,j,k)
                        dvy_dy = dvy_dy / K_y_half(j) + seisvar%memdvy_dy(i,j,k)
                        dvx_dy = dvx_dy / K_y_half(j) + seisvar%memdvx_dy(i,j,k)
                        dvz_dy = dvz_dy / K_y_half(j) + seisvar%memdvz_dy(i,j,k)
                        dvz_dz = dvz_dz / K_z_half(k) + seisvar%memdvz_dz(i,j,k)
                        dvy_dz = dvy_dz / K_z_half(k) + seisvar%memdvy_dz(i,j,k)

                        seisvar%sigxy(i,j,k) = seisvar%sigxy(i,j,k) + &
                        (   c16(i,k) * dvx_dx + c26(i,k) * dvy_dy + c36(i,k) * dvz_dz + &
                            c46(i,k) * (dvz_dy + dvy_dz) + c56(i,k) * (dvz_dx + dvx_dz) + &
                            c66(i,k) * (dvy_dx + dvx_dy) ) * time_params%dt

                    enddo
                enddo
            enddo

            ! Update sigmaxz, z-direction is full nodes
            do k = 1,domain%nz-1
                do j = 2,domain%ny
                    do i = 2,domain%nx

                        dvx_dx = (seisvar%vx(i,j,k) - seisvar%vx(i-1,j,k)) / domain%dx
                        dvy_dx = (seisvar%vy(i,j,k) - seisvar%vy(i-1,j,k)) / domain%dx
                        dvz_dx = (seisvar%vz(i,j,k) - seisvar%vz(i-1,j,k)) / domain%dx
                        dvy_dy = (seisvar%vy(i,j,k) - seisvar%vy(i,j-1,k)) / domain%dy
                        dvz_dy = (seisvar%vz(i,j,k) - seisvar%vz(i,j-1,k)) / domain%dy
                        dvx_dy = (seisvar%vx(i,j,k) - seisvar%vx(i,j-1,k)) / domain%dy
                        dvz_dz = (seisvar%vz(i,j,k+1) - seisvar%vz(i,j,k)) / domain%dz
                        dvx_dz = (seisvar%vx(i,j,k+1) - seisvar%vx(i,j,k)) / domain%dz
                        dvy_dz = (seisvar%vy(i,j,k+1) - seisvar%vy(i,j,k)) / domain%dz

                        seisvar%memdvx_dx(i,j,k) = b_x(i) * seisvar%memdvx_dx(i,j,k) + a_x(i) * dvx_dx
                        seisvar%memdvy_dx(i,j,k) = b_x(i) * seisvar%memdvy_dx(i,j,k) + a_x(i) * dvy_dx
                        seisvar%memdvz_dx(i,j,k) = b_x(i) * seisvar%memdvz_dx(i,j,k) + a_x(i) * dvz_dx
                        seisvar%memdvy_dy(i,j,k) = b_y(j) * seisvar%memdvy_dy(i,j,k) + a_y(j) * dvy_dy
                        seisvar%memdvx_dy(i,j,k) = b_y(j) * seisvar%memdvx_dy(i,j,k) + a_y(j) * dvx_dy
                        seisvar%memdvz_dy(i,j,k) = b_y(j) * seisvar%memdvz_dy(i,j,k) + a_y(j) * dvz_dy
                        seisvar%memdvz_dz(i,j,k) = b_z_half(k) * seisvar%memdvz_dz(i,j,k) + a_z_half(k) * dvz_dz
                        seisvar%memdvx_dz(i,j,k) = b_z_half(k) * seisvar%memdvx_dz(i,j,k) + a_z_half(k) * dvx_dz
                        seisvar%memdvy_dz(i,j,k) = b_z_half(k) * seisvar%memdvy_dz(i,j,k) + a_z_half(k) * dvy_dz

                        dvx_dx = dvx_dx / K_x(i) + seisvar%memdvx_dx(i,j,k)
                        dvy_dx = dvy_dx / K_x(i) + seisvar%memdvy_dx(i,j,k)
                        dvz_dx = dvz_dx / K_x(i) + seisvar%memdvz_dx(i,j,k) 
                        dvy_dy = dvy_dy / K_y(j) + seisvar%memdvy_dy(i,j,k)
                        dvx_dy = dvx_dy / K_y(j) + seisvar%memdvx_dy(i,j,k)
                        dvz_dy = dvz_dy / K_y(j) + seisvar%memdvz_dy(i,j,k)
                        dvz_dz = dvz_dz / K_z_half(k) + seisvar%memdvz_dz(i,j,k)
                        dvx_dz = dvx_dz / K_z_half(k) + seisvar%memdvx_dz(i,j,k)
                        dvy_dz = dvy_dz / K_z_half(k) + seisvar%memdvy_dz(i,j,k)

                        seisvar%sigxz(i,j,k) = seisvar%sigxz(i,j,k) + &
                            (   c15(i,k) * dvx_dx + c25(i,k) * dvy_dy + c35(i,k) * dvz_dz + &
                                c45(i,k) * ( dvx_dz + dvz_dx) + c55(i,k) * ( dvx_dz + dvz_dx) + &
                                c56(i,k) * ( dvx_dy + dvy_dx) ) * time_params%dt 

                    enddo
                enddo

                !   ! update sigmayz, y-direction is full nodes
                do j = 1,domain%ny-1
                    do i = 1,domain%nx-1

                        dvx_dx = (seisvar%vx(i+1,j,k) - seisvar%vx(i,j,k)) / domain%DX
                        dvy_dx = (seisvar%vy(i+1,j,k) - seisvar%vy(i,j,k)) / domain%DX
                        dvz_dx = (seisvar%vz(i+1,j,k) - seisvar%vz(i,j,k)) / domain%DX
                        dvy_dy = (seisvar%vy(i,j+1,k) - seisvar%vy(i,j,k)) / domain%DY
                        dvx_dy = (seisvar%vx(i,j+1,k) - seisvar%vx(i,j,k)) / domain%DY
                        dvz_dy = (seisvar%vz(i,j+1,k) - seisvar%vz(i,j,k)) / domain%DY
                        dvz_dz = (seisvar%vz(i,j,k+1) - seisvar%vz(i,j,k)) / domain%DZ
                        dvx_dz = (seisvar%vx(i,j,k+1) - seisvar%vx(i,j,k)) / domain%DZ 
                        dvy_dz = (seisvar%vy(i,j,k+1) - seisvar%vy(i,j,k)) / domain%DZ

                        seisvar%memdvx_dx(i,j,k) = b_x_half(i) * seisvar%memdvx_dx(i,j,k) + a_x_half(i) * dvx_dx
                        seisvar%memdvy_dx(i,j,k) = b_x_half(i) * seisvar%memdvy_dx(i,j,k) + a_x_half(i) * dvy_dx
                        seisvar%memdvz_dx(i,j,k) = b_x_half(i) * seisvar%memdvz_dx(i,j,k) + a_x_half(i) * dvz_dx
                        seisvar%memdvy_dy(i,j,k) = b_y_half(j) * seisvar%memdvy_dy(i,j,k) + a_y_half(j) * dvy_dy
                        seisvar%memdvx_dy(i,j,k) = b_y_half(j) * seisvar%memdvx_dy(i,j,k) + a_y_half(j) * dvx_dy
                        seisvar%memdvz_dy(i,j,k) = b_y_half(j) * seisvar%memdvz_dy(i,j,k) + a_y_half(j) * dvz_dy
                        seisvar%memdvz_dz(i,j,k) = b_z_half(k) * seisvar%memdvz_dz(i,j,k) + a_z_half(k) * dvz_dz
                        seisvar%memdvx_dz(i,j,k) = b_z_half(k) * seisvar%memdvx_dz(i,j,k) + a_z_half(k) * dvx_dz
                        seisvar%memdvy_dz(i,j,k) = b_z_half(k) * seisvar%memdvy_dz(i,j,k) + a_z_half(k) * dvy_dz

                        dvx_dx = dvx_dx / K_x_half(i) + seisvar%memdvx_dx(i,j,k)
                        dvy_dx = dvy_dx / K_x_half(i) + seisvar%memdvy_dx(i,j,k)
                        dvz_dx = dvz_dx / K_x_half(i) + seisvar%memdvz_dx(i,j,k)
                        dvy_dy = dvy_dy / K_y_half(j) + seisvar%memdvy_dy(i,j,k)
                        dvx_dy = dvx_dy / K_y_half(j) + seisvar%memdvx_dy(i,j,k)
                        dvz_dy = dvz_dy / K_y_half(j) + seisvar%memdvz_dy(i,j,k)
                        dvz_dz = dvz_dz / K_z_half(k) + seisvar%memdvz_dz(i,j,k)
                        dvx_dz = dvx_dz / K_z_half(k) + seisvar%memdvx_dz(i,j,k)
                        dvy_dz = dvy_dz / K_z_half(k) + seisvar%memdvy_dz(i,j,k)

                        seisvar%sigyz(i,j,k) = seisvar%sigyz(i,j,k)  + &
                            (   c14(i,k) * dvx_dx + c24(i,k) * dvy_dy + c34(i,k) * dvz_dz + &
                                c44(i,k) * ( dvy_dz + dvz_dy) + c45(i,k) * ( dvx_dz + dvz_dx) + &
                                c46(i,k) * ( dvy_dx + dvx_dy) ) * time_params%dt 
                    enddo
                enddo
            enddo

            !--------------------------------------------------------
            ! compute velocity and update memory variables for C-PML
            !--------------------------------------------------------
            do k = 2,domain%nz
                do j = 2,domain%ny
                    do i = 2,domain%nx
                        ! ds1/dx, ds6/dy, ds5,dz
                        deltarho = (4 * rho(i,k) + rho(i-1,k) + rho(i,k-1) )/6

                        dsigmaxx_dx = (seisvar%sigxx(i,j,k) - seisvar%sigxx(i-1,j,k) ) / domain%dx
                        dsigmaxy_dy = (seisvar%sigxy(i,j,k) - seisvar%sigxy(i,j-1,k) ) / domain%dy
                        dsigmaxz_dz = (seisvar%sigxz(i,j,k) - seisvar%sigxz(i,j,k-1) ) / domain%dz

                        seisvar%memdsigxx_dx(i,j,k) = b_x(i) * &
                            seisvar%memdsigxx_dx(i,j,k) + a_x(i) * dsigmaxx_dx
                        seisvar%memdsigxy_dy(i,j,k) = b_y(j) * &
                            seisvar%memdsigxy_dy(i,j,k) + a_y(j) * dsigmaxy_dy
                        seisvar%memdsigxz_dz(i,j,k) = b_z(k) * &
                            seisvar%memdsigxz_dz(i,j,k) + a_z(k) * dsigmaxz_dz

                        dsigmaxx_dx = dsigmaxx_dx / K_x(i) + seisvar%memdsigxx_dx(i,j,k)
                        dsigmaxy_dy = dsigmaxy_dy / K_y(j) + seisvar%memdsigxy_dy(i,j,k)
                        dsigmaxz_dz = dsigmaxz_dz / K_z(k) + seisvar%memdsigxz_dz(i,j,k) 

                        seisvar%vx(i,j,k) = seisvar%vx(i,j,k) * (1 - gammax(i,j) ) + &
                            (dsigmaxx_dx + dsigmaxy_dy + dsigmaxz_dz) * &
                            time_params%dt / deltarho !rho(i,k)
                    enddo
                enddo

                do j = 1,domain%ny-1
                    do i = 1,domain%nx-1
                        ! ds6/dx, ds2/dy, ds4/dz
                        deltarho = (4*rho(i,k) + rho(i+1,k) + rho(i,k-1) )/6

                        dsigmaxy_dx = ( seisvar%sigxy(i+1,j,k) - seisvar%sigxy(i,j,k) ) / domain%dx
                        dsigmayy_dy = ( seisvar%sigyy(i,j+1,k) - seisvar%sigyy(i,j,k) ) / domain%dy
                        dsigmayz_dz = ( seisvar%sigyz(i,j,k) - seisvar%sigyz(i,j,k-1) ) / domain%dz

                        seisvar%memdsigxy_dx(i,j,k) = b_x_half(i) * seisvar%memdsigxy_dx(i,j,k) + a_x_half(i) * dsigmaxy_dx
                        seisvar%memdsigyy_dy(i,j,k) = b_y_half(j) * seisvar%memdsigyy_dy(i,j,k) + a_y_half(j) * dsigmayy_dy
                        seisvar%memdsigyz_dz(i,j,k) = b_z(k) * seisvar%memdsigyz_dz(i,j,k) + a_z(k) * dsigmayz_dz

                        dsigmaxy_dx = dsigmaxy_dx / K_x_half(i) + seisvar%memdsigxy_dx(i,j,k)
                        dsigmayy_dy = dsigmayy_dy / K_y_half(j) + seisvar%memdsigyy_dy(i,j,k)
                        dsigmayz_dz = dsigmayz_dz / K_z(k) + seisvar%memdsigyz_dz(i,j,k)

                        seisvar%vy(i,j,k) = seisvar%vy(i,j,k) * (1 - gammay(i,j) )+ &
                            (dsigmaxy_dx + dsigmayy_dy + dsigmayz_dz) * &
                            time_params%dt / deltarho !rho(i,k)
                    enddo
                enddo
            enddo

            do k = 1,domain%nz-1
                do j = 2,domain%ny
                    do i = 1,domain%nx-1
                        ! ds5/dx, ds4/dy, ds3/dz
                        deltarho = ( rho(i+1,k) + rho(i,k+1) + 4*rho(i,k) )/6

                        dsigmaxz_dx = ( seisvar%sigxz(i+1,j,k) - seisvar%sigxz(i,j,k) ) / domain%dx
                        dsigmayz_dy = ( seisvar%sigyz(i,j,k) - seisvar%sigyz(i,j-1,k) ) / domain%dy
                        dsigmazz_dz = ( seisvar%sigzz(i,j,k+1) - seisvar%sigzz(i,j,k) ) / domain%dz

                        seisvar%memdsigxz_dx(i,j,k) = b_x_half(i) * seisvar%memdsigxz_dx(i,j,k) + a_x_half(i) * dsigmaxz_dx
                        seisvar%memdsigyz_dy(i,j,k) = b_y(j) * seisvar%memdsigyz_dy(i,j,k) + a_y(j) * dsigmayz_dy
                        seisvar%memdsigzz_dz(i,j,k) = b_z_half(k) * seisvar%memdsigzz_dz(i,j,k) + a_z_half(k) * dsigmazz_dz

                        dsigmaxz_dx = dsigmaxz_dx / K_x_half(i) + seisvar%memdsigxz_dx(i,j,k)
                        dsigmayz_dy = dsigmayz_dy / K_y(j) + seisvar%memdsigyz_dy(i,j,k)
                        dsigmazz_dz = dsigmazz_dz / K_z_half(k) + seisvar%memdsigzz_dz(i,j,k)

                        seisvar%vz(i,j,k) = seisvar%vz(i,j,k) * (1 - gammaz(i,j) )+ &
                            (dsigmaxz_dx + dsigmayz_dy + dsigmazz_dz) * &
                            time_params%dt / deltarho !rho(i,k)

                    enddo
                enddo
            enddo

            seisvar%vx(isource,jsource,ksource) = seisvar%vx(isource,jsource,ksource) + &
                    srcx(it) * time_params%dt / rho(isource,ksource)
            seisvar%vy(isource,jsource,ksource) = seisvar%vy(isource,jsource,ksource) + &
                    srcy(it) * time_params%dt / rho(isource,ksource)
            seisvar%vz(isource,jsource,ksource) = seisvar%vz(isource,jsource,ksource) + &
                    srcz(it) * time_params%dt / rho(isource,ksource)

            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            seisvar%vx(1,:,:) = 0.d0
            seisvar%vy(1,:,:) = 0.d0
            seisvar%vz(1,:,:) = 0.d0

            seisvar%vx(:,1,:) = 0.d0
            seisvar%vy(:,1,:) = 0.d0
            seisvar%vz(:,1,:) = 0.d0

            seisvar%vx(:,:,1) = 0.d0
            seisvar%vy(:,:,1) = 0.d0
            seisvar%vz(:,:,1) = 0.d0

            seisvar%vx(domain%nx,:,:) = 0.d0
            seisvar%vy(domain%nx,:,:) = 0.d0
            seisvar%vz(domain%nx,:,:) = 0.d0

            seisvar%vx(:,domain%ny,:) = 0.d0
            seisvar%vy(:,domain%ny,:) = 0.d0
            seisvar%vz(:,domain%ny,:) = 0.d0

            seisvar%vx(:,:,domain%nz) = 0.d0
            seisvar%vy(:,:,domain%nz) = 0.d0
            seisvar%vz(:,:,domain%nz) = 0.d0

            ! check norm of velocity to make sure the solution isn't diverging
            velocnorm = maxval( sqrt(seisvar%vx**2 + seisvar%vy**2 + seisvar%vz**2) )
            ! print *,'Time step # ',it,' out of ',time_step
            ! print *,'Time: ',(it-1)*DT,' seconds'
            ! print *,'Max vals for vx, vy, vz: ', maxval(vx), maxval(vy), maxval(vz)

            if (velocnorm > stability_threshold) stop 'code became unstable and blew up'

            ! Write the velocity values to an unformatted binary file
            call write_image(seisvar%vx, domain, source, it, 'Vx', SINGLE)
            call write_image(seisvar%vy, domain, source, it, 'Vy', SINGLE)
            call write_image(seisvar%vz, domain, source, it, 'Vz', SINGLE)
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
    subroutine electromag2(domain, time_params, source, emvar, SINGLE_OUTPUT)
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
        type(Time_Parameters_Type), intent(in) :: time_params
        type(Source_Type), intent(in) :: source
        type(Electromagnetic2_Variables_Type), intent(inout) :: emvar
        logical, intent(in), optional :: SINGLE_OUTPUT
         
        ! Local variabless
        real(real64), dimension(domain%nx,domain%nz) :: epsilonx, epsilonz, &
                                            sigmax, sigmaz

        ! real(real64) :: DT
        real(real64), dimension(time_params%time_steps) :: srcx, srcz
        integer :: isource, jsource, i, j, it

        ! Coefficients for the finite difference scheme
        real(real64), dimension(domain%nx,domain%nz) :: caEx, cbEx
        real(real64), dimension(domain%nx,domain%nz) :: caEz, cbEz
        real(real64) :: daHy, dbHy
        real(real64) :: value_dEx_dz, value_dEz_dx, value_dHy_dz, value_dHy_dx

        ! 1D arrays for the damping profiles
        real(real64), dimension(domain%nx) :: K_x, alpha_x, a_x, b_x, &
                                        K_x_half, alpha_x_half, a_x_half, b_x_half
        real(real64), dimension(domain%nz) :: K_z, alpha_z, a_z, b_z, &
                                        K_z_half, alpha_z_half, a_z_half, b_z_half

        ! Velocity normalization factor
        real(real64) :: velocnorm

        ! Boolean flag to save as double precision or single precision
        logical :: SINGLE

        ! Check if SINGLE_OUTPUT is provided, default to single precision if not
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        endif
        
        ! ----------------------------------------------------------------------

        ! ======================================================================
        ! ----------------------- Load Permittivity Coefficients ----------------------
        call material_rw('eps11.dat', emvar%eps11, .TRUE.)
        call material_rw('eps13.dat', emvar%eps13, .TRUE.)
        call material_rw('eps33.dat', emvar%eps33, .TRUE.) ! We will change y to z soon
        call material_rw('sig11.dat', emvar%sig11, .TRUE.)
        call material_rw('sig13.dat', emvar%sig13, .TRUE.)
        call material_rw('sig33.dat', emvar%sig33, .TRUE.)

        ! ------------------------ Assign some constants -----------------------
        ! Assign the source location indices
        isource = source%xind + domain%cpml
        jsource = source%zind + domain%cpml

        ! Define the 
        ! DT = minval( (/dx, dz/) )/ ( 2.d0 * Clight/sqrt( minval( (/ eps11, eps33 /) ) ) ) 

        ! Compute the coefficients of the FD scheme. First scale the relative 
        ! permittivity and permeabilities to get the absolute values 
        epsilonx(:,:) = (emvar%eps11 + emvar%eps13)*eps0
        epsilonz(:,:) = (emvar%eps13 + emvar%eps33)*eps0
        sigmax(:,:) = emvar%sig11 + emvar%sig13 
        sigmaz(:,:) = emvar%sig13 + emvar%sig33 

        ! We need to change sigma to dsigma, same for epsilon
        caEx(:,:) = ( 1.0d0 - sigmax * time_params%dt / &
                    ( 2.0d0 * epsilonx ) ) / &
                    ( 1.0d0 + sigmax * time_params%dt / &
                    (2.0d0 * epsilonx ) )
        cbEx(:,:) = (time_params%dt / epsilonx ) / &
                    ( 1.0d0 + sigmax * time_params%dt / &
                    ( 2.0d0 * epsilonx ) )

        caEz(:,:) = ( 1.0d0 - sigmaz * time_params%dt / &
                    ( 2.0d0 * epsilonz ) ) / &
                    ( 1.0d0 + sigmaz * time_params%dt / &
                    (2.0d0 * epsilonz ) )
        cbEz(:,:) = (time_params%dt / epsilonz ) / &
                    ( 1.0d0 + sigmaz * time_params%dt / &
                    ( 2.0d0 * epsilonz ) )
        daHy = time_params%dt/(4.0d0*mu0*mu)
        dbHy = time_params%dt/mu0 !dt/(mu*mu*dx*(1+daHy) ) 
        daHy = 1.0d0 ! (1-daHy)/(1+daHy) ! 

        ! ================================ LOAD SOURCE ================================
        call loadsource('electromagneticsourcex.dat', time_params%time_steps, srcx)
        call loadsource('electromagneticsourcez.dat', time_params%time_steps, srcz)

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
        
        ! Load initial conditions
        call material_rw('initialconditionEx.dat', emvar%Ex, .TRUE.)
        call material_rw('initialconditionHy.dat', emvar%Hy, .TRUE.)
        call material_rw('initialconditionEz.dat', emvar%Ez, .TRUE.)

        ! PML
        emvar%memdEx_dz(:,:) = 0.0d0
        emvar%memdEz_dx(:,:) = 0.0d0
        emvar%memdHy_dx(:,:) = 0.0d0
        emvar%memdHy_dz(:,:) = 0.0d0

        !---
        !---  beginning of time loop
        !---
        do it = 1,time_params%time_steps
            !--------------------------------------------------------
            ! compute magnetic field and update memory variables for C-PML
            !--------------------------------------------------------
            do j = 1,domain%nz-1  
                do i = 1,domain%nx-1
                
                    ! Values needed for the magnetic field updates
                    value_dEx_dz = ( emvar%Ex(i,j+1) - emvar%Ex(i,j) )/domain%dz
                    emvar%memdEx_dz(i,j) = b_z(j) * emvar%memdEx_dz(i,j) + a_z(j) * value_dEx_dz
                    value_dEx_dz = value_dEx_dz/ K_z(j) + emvar%memdEx_dz(i,j)

                    ! The rest of the equation needed for agnetic field updates
                    value_dEz_dx = ( emvar%Ez(i+1,j) - emvar%Ez(i,j) )/domain%dx
                    emvar%memdEz_dx(i,j) = b_x(i) * emvar%memdEz_dx(i,j) + a_x(i) * value_dEz_dx
                    value_dEz_dx = value_dEz_dx/ K_x(i) + emvar%memdEz_dx(i,j)

                    ! Now update the Magnetic field
                    emvar%Hy(i,j) = daHy*emvar%Hy(i,j) + dbHy*( value_dEz_dx + value_dEx_dz )

                enddo  
            enddo

            !--------------------------------------------------------
            ! compute electric field and update memory variables for C-PML
            !--------------------------------------------------------
            ! Compute the differences in the y-direction
            do j = 2,domain%nz
                do i = 1,domain%nx
                    ! Update the Ex field
                    value_dHy_dz = ( emvar%Hy(i,j) - emvar%Hy(i,j-1) )/domain%dz ! this is nz-1 length vector
                    emvar%memdHy_dz(i,j) = b_z(j) * emvar%memdHy_dz(i,j) + a_z(j) * value_dHy_dz
                    value_dHy_dz = value_dHy_dz/K_z(j) + emvar%memdHy_dz(i,j)

                    ! Ex(i,j) = (( caEx(i,j) + caEx(i,j-1) )/2) * Ex(i,j) + &
                    !     (( cbEx(i,j) + cbEx(i,j-1) )/2 ) * value_dHy_dz
                    emvar%Ex(i,j) = caEx(i,j) * emvar%Ex(i,j) + cbEx(i,j) * value_dHy_dz
                enddo
            enddo

            do j = 1,domain%nz
                do i = 2,domain%nx
                    ! Update the Ez field
                    value_dHy_dx = ( emvar%Hy(i,j) - emvar%Hy(i-1,j) )/domain%dx
                    emvar%memdHy_dx(i,j) = b_x_half(i) * emvar%memdHy_dx(i,j) + a_x_half(i) * value_dHy_dx
                    value_dHy_dx = value_dHy_dx/K_x_half(i) + emvar%memdHy_dx(i,j)
                    
                    ! Ez(i,j) = (( caEz(i,j) + caEz(i-1,j) )/2) * Ez(i,j) + &
                    !     (( cbEz(i,j) + cbEz(i-1,j) )/2) * value_dHy_dx 
                    emvar%Ez(i,j) = caEz(i,j) * emvar%Ez(i,j) + cbEz(i,j) * value_dHy_dx 
                enddo
            enddo

            !----------------------------------------------------------------------------
            emvar%Ex(isource,jsource) = emvar%Ex(isource,jsource) + &
                            srcx(it) * time_params%dt / emvar%eps11(isource,jsource)
            emvar%Ez(isource,jsource) = emvar%Ez(isource,jsource) + &
                            srcz(it) * time_params%dt / emvar%eps33(isource,jsource) 
            
            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            emvar%Ex(1,:) = 0.d0
            emvar%Ex(domain%nx,:) = 0.d0
            emvar%Ex(:,1) = 0.d0
            emvar%Ex(:,domain%nz) = 0.d0

            emvar%Ez(1,:) = 0.d0
            emvar%Ez(domain%nx,:) = 0.d0
            emvar%Ez(:,1) = 0.d0
            emvar%Ez(:,domain%nz) = 0.d0

            emvar%Hy(1,:) = 0.d0
            emvar%Hy(domain%nx,:) = 0.d0
            emvar%Hy(:,1) = 0.d0
            emvar%Hy(:,domain%nz) = 0.d0

            ! print maximum of norm of velocity
            velocnorm = maxval(sqrt(emvar%Ex**2 + emvar%Ez**2))
            if (velocnorm > stability_threshold) stop 'code became unstable and blew up'

            call write_image(emvar%Ex, domain, source, it, 'Ex', SINGLE)
            call write_image(emvar%Ez, domain, source, it, 'Ez', SINGLE)
        enddo
    end subroutine electromag2


    ! =========================================================================
    subroutine electromag25(domain, source, time_params, emvar, SINGLE_OUTPUT)
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
        type(Time_Parameters_Type), intent(in) :: time_params
        type(Source_Type), intent(in) :: source
        type(Electromagnetic3_Variables_Type), intent(inout) :: emvar
        logical, intent(in), optional :: SINGLE_OUTPUT
        
        ! Local variables
        real(real64), dimension(domain%nx,domain%nz) :: epsilonx, epsilony, epsilonz, &
                                            sigmax, sigmay, sigmaz

        ! real(real64) :: DT
        real(real64) :: velocnorm
        integer :: isource, jsource, ksource, i, j, k, it

        ! Coefficients for the finite difference scheme
        real(real64), dimension(domain%nx,domain%nz) :: caEx, cbEx, caEy, cbEy, caEz, cbEz
        real(real64) :: daHx, dbHx, daHy, dbHy, daHz, dbHz

        real(real64) :: dEx_dy, dEy_dx, dEy_dz, dEz_dy, dEz_dx, dEx_dz, &
                        dHx_dy, dHx_dz, dHy_dx, dHy_dz, dHz_dy, dHz_dx


        ! Source arrays
        real(real64), dimension(time_params%time_steps) :: srcx, srcy, srcz

        ! 1D arrays for the damping profiles in each direction
        real(real64), dimension(domain%nx) :: K_x, alpha_x, a_x, b_x, K_x_half, alpha_x_half, a_x_half, b_x_half
        real(real64), dimension(domain%ny) :: K_y, alpha_y, a_y, b_y, K_y_half, alpha_y_half, a_y_half, b_y_half
        real(real64), dimension(domain%nz) :: K_z, alpha_z, a_z, b_z, K_z_half, alpha_z_half, a_z_half, b_z_half

        ! Boolean flag to save as double precision or single precision
        logical :: SINGLE

        ! Check if SINGLE_OUTPUT is provided, default to single precision if not
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        endif
        
        ! ------------------------ Load Permittivity Coefficients ------------------------
        ! Load Epsilon
        call material_rw('eps11.dat', emvar%eps11, .TRUE.)
        call material_rw('eps12.dat', emvar%eps12, .TRUE.)
        call material_rw('eps13.dat', emvar%eps13, .TRUE.)
        call material_rw('eps22.dat', emvar%eps22, .TRUE.)
        call material_rw('eps23.dat', emvar%eps23, .TRUE.)
        call material_rw('eps33.dat', emvar%eps33, .TRUE.)
        ! Load Sigma
        call material_rw('sig11.dat', emvar%sig11, .TRUE.)
        call material_rw('sig12.dat', emvar%sig12, .TRUE.)
        call material_rw('sig13.dat', emvar%sig13, .TRUE.)
        call material_rw('sig22.dat', emvar%sig22, .TRUE.)
        call material_rw('sig23.dat', emvar%sig23, .TRUE.)
        call material_rw('sig33.dat', emvar%sig33, .TRUE.)

        ! ------------------------ Assign some constants -----------------------
        ! Assign the source location indices
        isource = source%xind + domain%cpml
        jsource = source%yind + domain%cpml
        ksource = source%zind + domain%cpml

        ! Define the 
        ! DT = minval( (/dx, dy, dz/) )/ ( 2.0d0 * Clight/ sqrt( minval( (/ eps11, eps22, eps33 /) ) ) )

        ! Compute the coefficients of the FD scheme. First scale the relative 
        ! permittivity and permeabilities to get the absolute values 
        epsilonx(:,:) = (emvar%eps11 + emvar%eps12 + emvar%eps13)*eps0 
        epsilony(:,:) = (emvar%eps12 + emvar%eps22 + emvar%eps23)*eps0
        epsilonz(:,:) = (emvar%eps13 + emvar%eps23 + emvar%eps33)*eps0
        sigmax(:,:) = emvar%sig11 + emvar%sig12 + emvar%sig13
        sigmay(:,:) = emvar%sig12 + emvar%sig22 + emvar%sig23
        sigmaz(:,:) = emvar%sig13 + emvar%sig23 + emvar%sig33

        caEx(:,:) = ( 1.0d0 - sigmax * time_params%dt / &
                    (2.0d0 * epsilonx ) ) / &
                    ( 1.0d0 + sigmax * time_params%dt / &
                    (2.0d0 * epsilonx ) )
        cbEx(:,:) = (time_params%dt / epsilonx ) / &
                    ( 1.0d0 + sigmax * time_params%dt / &
                    ( 2.0d0 * epsilonx ) )

        caEy(:,:) = ( 1.0d0 - sigmay * time_params%dt / (2.0d0 * epsilony ) ) / &
                    ( 1.0d0 + sigmay * time_params%dt / (2.0d0 * epsilony ) )
        cbEy(:,:) = (time_params%dt / epsilony ) / &
                    ( 1.0d0 + sigmay * time_params%dt / ( 2.0d0 * epsilony ) )

        caEz(:,:) = ( 1.0d0 - sigmaz * time_params%dt / &
                    ( 2.0d0 * epsilonz ) ) / &
                    ( 1.0d0 + sigmaz * time_params%dt / (2.0d0 * epsilonz ) )
        cbEz(:,:) = (time_params%dt / epsilonz ) / &
                    ( 1.0d0 + sigmaz * time_params%dt / (2.0d0 * epsilonz ) )

        daHx = time_params%dt/(4.0d0*mu0*mu)
        dbHx = time_params%dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
        daHx = 1.0d0 ! (1-daHz)/(1+daHz) ! 

        daHy = time_params%dt/(4.0d0*mu0*mu)
        dbHy = time_params%dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
        daHy = 1.0d0 ! (1-daHz)/(1+daHz) ! 

        daHz = time_params%dt/(4.0d0*mu0*mu)
        dbHz = time_params%dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
        daHz = 1.0d0 ! (1-daHz)/(1+daHz) ! 


        ! ----------------------------------------------------------------------
        !---
        !--- program starts here
        !---

        ! ================================ LOAD SOURCE ================================

        call loadsource('electromagneticsourcex.dat', time_params%time_steps, srcx)
        call loadsource('electromagneticsourcey.dat', time_params%time_steps, srcy)
        call loadsource('electromagneticsourcez.dat', time_params%time_steps, srcz)

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
        call material_rw('initialconditionEx.dat', emvar%Ex, .TRUE.)
        call material_rw('initialconditionEy.dat', emvar%Ey, .TRUE.)
        call material_rw('initialconditionEz.dat', emvar%Ez, .TRUE.)
        
        call material_rw('initialconditionHx.dat', emvar%Hx, .TRUE.)
        call material_rw('initialconditionHy.dat', emvar%Hy, .TRUE.)
        call material_rw('initialconditionHz.dat', emvar%Hz, .TRUE.)

        ! PML
        emvar%memdEx_dy(:,:,:) = 0.0d0
        emvar%memdEy_dx(:,:,:) = 0.0d0
        emvar%memdEx_dz(:,:,:) = 0.0d0
        emvar%memdEz_dx(:,:,:) = 0.0d0
        emvar%memdEz_dy(:,:,:) = 0.0d0
        emvar%memdEy_dz(:,:,:) = 0.0d0

        emvar%memdHz_dx(:,:,:) = 0.0d0
        emvar%memdHx_dz(:,:,:) = 0.0d0
        emvar%memdHz_dy(:,:,:) = 0.0d0
        emvar%memdHy_dz(:,:,:) = 0.0d0
        emvar%memdHx_dy(:,:,:) = 0.0d0
        emvar%memdHy_dx(:,:,:) = 0.0d0

        ! ---
        ! ---  beginning of time loop
        ! ---
        do it = 1,time_params%time_steps
            !--------------------------------------------------------
            ! compute magnetic field and update memory variables for C-PML
            !--------------------------------------------------------
            ! Update Hx
            do k = 1,domain%nz-1
                do i = 1,domain%nx-1  
                    do j = 1,domain%ny-1
                        ! Values needed for the magnetic field updates
                        dEz_dy = ( emvar%Ez(i,j,k) - emvar%Ez(i,j+1,k) )/domain%dy
                        emvar%memdEz_dy(i,j,k) = b_y_half(j) * emvar%memdEz_dy(i,j,k) + a_y_half(j) * dEz_dy
                        dEz_dy = dEz_dy/ K_y_half(j) + emvar%memdEz_dy(i,j,k)

                        ! The rest of the equation needed for agnetic field updates
                        dEy_dz = ( emvar%Ey(i,j,k+1) - emvar%Ey(i,j,k) )/domain%dz
                        emvar%memdEy_dz(i,j,k) = b_z_half(k) * emvar%memdEy_dz(i,j,k) + a_z_half(k) * dEy_dz
                        dEy_dz = dEy_dz/ K_z_half(k) + emvar%memdEy_dz(i,j,k)

                        ! Now update the Magnetic field
                        emvar%Hx(i,j,k) = daHx*emvar%Hx(i,j,k) + dbHx*( dEy_dz + dEz_dy )
                    enddo
                enddo  
            enddo

                ! Update Hy
            do k = 1,domain%nz-1
                do i = 1,domain%nx-1      
                    do j = 1,domain%ny-1
                    
                        ! Values needed for the magnetic field updates
                        dEx_dz = ( emvar%Ex(i,j,k) - emvar%Ex(i,j,k+1) )/domain%dz
                        emvar%memdEx_dz(i,j,k) = b_z(k) * emvar%memdEx_dz(i,j,k) + &
                            a_z(k) * dEx_dz
                        dEx_dz = dEx_dz/ K_z(k) + emvar%memdEx_dz(i,j,k)

                        ! The rest of the equation needed for agnetic field updates
                        dEz_dx = ( emvar%Ez(i+1,j,k) - emvar%Ez(i,j,k) )/domain%dx
                        emvar%memdEz_dx(i,j,k) = b_x(i) * emvar%memdEz_dx(i,j,k) + &
                            a_x(i) * dEz_dx
                        dEz_dx = dEz_dx/ K_x(i) + emvar%memdEz_dx(i,j,k)

                        ! Now update the Magnetic field
                        emvar%Hy(i,j,k) = daHy*emvar%Hy(i,j,k) + dbHy*( dEz_dx + dEx_dz )

                    enddo
                enddo  
            enddo

                ! Update Hz
            do k = 2,domain%nz-1
                do i = 1,domain%nx-1      
                    do j = 1,domain%ny-1
                        ! Values needed for the magnetic field updates
                        dEx_dy = ( emvar%Ex(i,j+1,k) - emvar%Ex(i,j,k) )/domain%dy
                        emvar%memdEx_dy(i,j,k) = b_y(j) * emvar%memdEx_dy(i,j,k) + & 
                            a_y(j) * dEx_dy
                        dEx_dy = dEx_dy/ K_y(j) + emvar%memdEx_dy(i,j,k)

                        ! The rest of the equation needed for agnetic field updates
                        dEy_dx = ( emvar%Ey(i,j,k) - emvar%Ey(i+1,j,k) )/domain%dx
                        emvar%memdEy_dx(i,j,k) = b_x(i) * emvar%memdEy_dx(i,j,k) + & 
                            a_x(i) * dEy_dx
                        dEy_dx = dEy_dx/ K_x(i) + emvar%memdEy_dx(i,j,k)

                        ! Now update the Magnetic field
                        emvar%Hz(i,j,k) = daHz*emvar%Hz(i,j,k) + dbHz*( dEy_dx + dEx_dy )
                    enddo
                enddo  
            enddo

            !--------------------------------------------------------
            ! compute electric field and update memory variables for C-PML
            !--------------------------------------------------------
            ! Compute the differences in the x-direction
            do k = 2,domain%nz-1
                do i = 1,domain%nx-1
                    do j = 2,domain%ny-1  
                        ! Update the Ex field
                        dHz_dy = ( emvar%Hz(i,j,k) - emvar%Hz(i,j-1,k) )/domain%dy
                        emvar%memdHz_dy(i,j,k) = b_y_half(j) * emvar%memdHz_dy(i,j,k) + & 
                            a_y_half(j) * dHz_dy
                        dHz_dy = dHz_dy/K_y_half(j) + emvar%memdHz_dy(i,j,k)

                        ! Changed from half to full node positions 
                        dHy_dz = ( emvar%Hy(i,j,k-1) - emvar%Hy(i,j,k) )/domain%dz
                        emvar%memdHy_dz(i,j,k) = b_z(k) * emvar%memdHy_dz(i,j,k) + &
                            a_z(k) * dHy_dz
                        dHy_dz = dHy_dz/K_z(k) + emvar%memdHy_dz(i,j,k)
                        
                        emvar%Ex(i,j,k) = caEx(i,k)*emvar%Ex(i,j,k) + & 
                        cbEx(i,k)*(dHz_dy + dHy_dz) 
                    enddo
                enddo

                ! ! Compute the differences in the y-direction
                do i = 2,domain%nx-1 
                    do j = 1,domain%ny-1 
                        ! Update the Ey field
                        dHz_dx = ( emvar%Hz(i-1,j,k) - emvar%Hz(i,j,k) )/domain%dx ! this is ny-1 length vector
                        emvar%memdHz_dx(i,j,k) = b_x_half(i) * emvar%memdHz_dx(i,j,k) + & 
                            a_x_half(i) * dHz_dx
                        dHz_dx = dHz_dx/K_x_half(i) + emvar%memdHz_dx(i,j,k)

                        dHx_dz = ( emvar%Hx(i,j,k) - emvar%Hx(i,j,k-1) )/domain%dz ! this is ny-1 length vector
                        emvar%memdHx_dz(i,j,k) = b_z_half(k) * emvar%memdHx_dz(i,j,k) + &
                            a_z_half(k) * dHx_dz
                        dHx_dz = dHx_dz/K_z_half(k) + emvar%memdHx_dz(i,j,k)

                        ! Ey(i,j,k) = ( ( 4*caEy(i,k) + caEy(i-1,k) + caEy(i,k-1) )/6) * Ey(i,j,k) + & 
                        ! ( ( 4*cbEy(i,k) + cbEy(i-1,k) + cbEy(i,k-1) )/6 ) * & 
                        ! (dHz_dx + dHx_dz)
                        emvar%Ey(i,j,k) = caEy(i,k) * emvar%Ey(i,j,k) + cbEy(i,k) * (dHz_dx + dHx_dz)
                    enddo
                enddo
            enddo 

                ! Compute the differences in the z-direction
            do k = 1,domain%nz-1
                do i = 2,domain%nx-1  
                    do j = 2,domain%ny-1
                        ! Update the Ez field
                        dHx_dy = ( emvar%Hx(i,j-1,k) - emvar%Hx(i,j,k) )/domain%dy
                        emvar%memdHx_dy(i,j,k) = b_y_half(j) * emvar%memdHx_dy(i,j,k) + &
                            a_y_half(j) * dHx_dy
                        dHx_dy = dHx_dy/K_y_half(j) + emvar%memdHx_dy(i,j,k)

                        dHy_dx = ( emvar%Hy(i,j,k) - emvar%Hy(i-1,j,k) )/domain%dx
                        emvar%memdHy_dx(i,j,k) = b_x_half(i) * emvar%memdHy_dx(i,j,k) + &
                            a_x_half(i) * dHy_dx
                        dHy_dx = dHy_dx/K_x_half(i) + emvar%memdHy_dx(i,j,k)
                        
                        emvar%Ez(i,j,k) = ( ( 4*caEz(i,k) + caEz(i-1,k) + caEz(i,k+1) )/6 ) * &
                            emvar%Ez(i,j,k) + ( ( 4*cbEz(i,k) + cbEz(i-1,k) + cbEz(i,k+1) )/6 ) * & 
                        (dHx_dy + dHy_dx)
                    enddo
                enddo
            enddo


            ! add the source (force vector located at a given grid point)
            emvar%Ex(isource,jsource,ksource) = emvar%Ex(isource,jsource,ksource) + & 
                        srcx(it) * time_params%dt / emvar%eps11(isource,ksource)
            emvar%Ey(isource,jsource,ksource) = emvar%Ey(isource,jsource,ksource) + & 
                        srcy(it) * time_params%dt / emvar%eps22(isource,ksource) 
            emvar%Ez(isource,jsource,ksource) = emvar%Ez(isource,jsource,ksource) + & 
                        srcz(it) * time_params%dt / emvar%eps33(isource,ksource)
            
            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            emvar%Ex(1,:,:) = 0.0d0
            emvar%Ex(:,1,:) = 0.0d0
            emvar%Ex(:,:,1) = 0.0d0
            emvar%Ex(domain%nx,:,:) = 0.0d0
            emvar%Ex(:,domain%ny,:) = 0.0d0
            emvar%Ex(:,:,domain%nz) = 0.0d0 

            emvar%Ey(1,:,:) = 0.0d0
            emvar%Ey(:,1,:) = 0.0d0
            emvar%Ey(:,:,1) = 0.0d0
            emvar%Ey(domain%nx,:,:) = 0.0d0
            emvar%Ey(:,domain%ny,:) = 0.0d0
            emvar%Ey(:,:,domain%nz) = 0.0d0
            
            emvar%Ez(1,:,:) = 0.0d0
            emvar%Ez(:,1,:) = 0.0d0
            emvar%Ez(:,:,1) = 0.0d0
            emvar%Ez(domain%nx,:,:) = 0.0d0
            emvar%Ez(:,domain%ny,:) = 0.0d0
            emvar%Ez(:,:,domain%nz) = 0.0d0
            
            emvar%Hx(1,:,:) = 0.0d0
            emvar%Hx(:,1,:) = 0.0d0
            emvar%Hx(:,:,1) = 0.0d0
            emvar%Hx(domain%nx,:,:) = 0.0d0
            emvar%Hx(:,domain%ny,:) = 0.0d0
            emvar%Hx(:,:,domain%nz) = 0.0d0

            emvar%Hy(1,:,:) = 0.0d0
            emvar%Hy(:,1,:) = 0.0d0
            emvar%Hy(:,:,1) = 0.0d0
            emvar%Hy(domain%nx,:,:) = 0.0d0
            emvar%Hy(:,domain%ny,:) = 0.0d0
            emvar%Hy(:,:,domain%nz) = 0.0d0
            
            emvar%Hz(1,:,:) = 0.0d0
            emvar%Hz(:,1,:) = 0.0d0
            emvar%Hz(:,:,1) = 0.0d0
            emvar%Hz(domain%nx,:,:) = 0.0d0
            emvar%Hz(:,domain%ny,:) = 0.0d0
            emvar%Hz(:,:,domain%nz) = 0.0d0

            ! check norm of velocity to make sure the solution isn't diverging
            velocnorm = maxval(sqrt(emvar%Ex**2.0d0 + emvar%Ey**2.0d0 + emvar%Ez**2.0d0) )
            if (velocnorm > stability_threshold) stop 'code became unstable and blew up'
            ! print *,'Max vals for Ex, Ey, Ez: ', maxval(Ex), maxval(Ey), maxval(Ez)

            ! print *, maxval(Ex), maxval(Ey), maxval(Ez)
            call write_image(emvar%Ex, domain, source, it, 'Ex', SINGLE)
            call write_image(emvar%Ey, domain, source, it, 'Ey', SINGLE)
            call write_image(emvar%Ez, domain, source, it, 'Ez', SINGLE)

        enddo   ! end of time loop
    end subroutine electromag25



end module cpmlfdtd