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
        real(real64), allocatable :: gammax(:,:), gammaz(:,:)        
        real(real64), allocatable :: srcx(:), srcz(:) ! The vector time series of the source
        
        ! Model variables
        real(real64), allocatable :: vx(:,:),vz(:,:),sigmaxx(:,:),sigmazz(:,:),sigmaxz(:,:)
        real(real64), allocatable :: &
                memory_dvx_dx(:,:), &
                memory_dvx_dz(:,:), &
                memory_dvz_dx(:,:), &
                memory_dvz_dz(:,:), &
                memory_dsigmaxx_dx(:,:), &
                memory_dsigmazz_dz(:,:), &
                memory_dsigmaxz_dx(:,:), &
                memory_dsigmaxz_dz(:,:)
        
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
        
        ! Allocate the arrays based on runtime values of domain%nx and domain%nz
        allocate(c11(nx, nz), c13(nx, nz), c15(nx, nz), &
                 c33(nx, nz), c35(nx, nz), c55(nx, nz), rho(nx, nz))
        allocate(K_x(nx), alpha_x(nx), a_x(nx), b_x(nx), &
                 K_x_half(nx), alpha_x_half(nx), a_x_half(nx), b_x_half(nx))
        allocate(K_z(nz), alpha_z(nz), a_z(nz), b_z(nz), &
                K_z_half(domain%nz), alpha_z_half(domain%nz), a_z_half(domain%nz), b_z_half(domain%nz))
        allocate(gammax(nx, nz), gammaz(nx, nz))
        allocate(srcx(source%time_steps), srcz(source%time_steps))
        
        ! Allocate more
        allocate(memory_dvx_dx(nx, nz), memory_dvx_dz(nx, nz))
        allocate(memory_dvz_dx(nx, nz), memory_dvz_dz(nx, nz))
        allocate(memory_dsigmaxx_dx(nx, nz), memory_dsigmazz_dz(nx, nz))
        allocate(memory_dsigmaxz_dx(nx, nz), memory_dsigmaxz_dz(nx, nz))
        allocate(vx(nx, nz), vz(nx, nz))
        allocate(sigmaxx(nx, nz), sigmazz(nx, nz), &
                 sigmaxz(nx, nz))
            
        ! -------------------- Load Stiffness Coefficients --------------------
    
        call material_rw2('c11.dat', c11, .TRUE.)
        call material_rw2('c13.dat', c13, .TRUE.)
        call material_rw2('c15.dat', c15, .TRUE.)
        call material_rw2('c33.dat', c33, .TRUE.)
        call material_rw2('c35.dat', c35, .TRUE.)
        call material_rw2('c55.dat', c55, .TRUE.)
        call material_rw2('rho.dat', rho, .TRUE.)
                
        ! ------------------- Load Attenuation Coefficients --------------------
        call material_rw2('gamma_x.dat', gammax, .TRUE.)
        call material_rw2('gamma_z.dat', gammaz, .TRUE.)
        
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
        
            !------------------------------------------------------------
            ! compute stress sigma and update memory variables for C-PML
            !------------------------------------------------------------
            !$omp parallel do private(i, j, deltarho, &
            ! value_dvx_dx, value_dvx_dz, value_dvz_dx, value_dvz_dz)
            do j = 2,nz
                do i = 1,nx-1
        
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
                    
                    sigmaxx(i,j) = sigmaxx(i,j) + &
                        (   c11(i,j) * value_dvx_dx + &
                            c13(i,j) * value_dvz_dz + &
                            c15(i,j) * (value_dvz_dx + value_dvx_dz) ) * dt
                    sigmazz(i,j) = sigmazz(i,j) + &
                        (   c13(i,j) * value_dvx_dx + &
                            c33(i,j) * value_dvz_dz + &
                            c35(i,j) * (value_dvz_dx + value_dvx_dz) ) * dt
        
                enddo
            enddo
            
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
            
                    sigmaxz(i,j) = sigmaxz(i,j) + &
                        (   c15(i,j)  * value_dvx_dx + & 
                            c35(i,j)  * value_dvz_dz + &
                            c55(i,j) * (value_dvz_dx + value_dvx_dz) ) * dt
        
                enddo
            enddo
            !$omp end parallel do
            
            !--------------------------------------------------------
            ! compute velocity and update memory variables for C-PML
            !--------------------------------------------------------
            !$omp parallel do private(i, j, deltarho, value_dsigmaxx_dx, value_dsigmazz_dz, &
            !    value_dsigmaxz_dx, value_dsigmaxz_dz)

            do j = 2,domain%nz
                do i = 2,domain%nx
        
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
            
                    vx(i,j) = vx(i,j) + (value_dsigmaxx_dx + value_dsigmaxz_dz) * &
                            ( dt / (deltarho - dt * gammax(i,j) ))
                enddo
            enddo
            
            do j = 1,domain%nz-1
                do i = 1,domain%nx-1
        
                    deltarho = ( 2*rho(i,j) + rho(i+1,j) + rho(i,j+1) )/4
                    value_dsigmaxz_dx = (sigmaxz(i+1,j) - sigmaxz(i,j)) / domain%dx
                    value_dsigmazz_dz = (sigmazz(i,j+1) - sigmazz(i,j)) / domain%dz
            
                    memory_dsigmaxz_dx(i,j) = b_x_half(i) * memory_dsigmaxz_dx(i,j) + &
                                a_x_half(i) * value_dsigmaxz_dx
                    memory_dsigmazz_dz(i,j) = b_z_half(j) * memory_dsigmazz_dz(i,j) + &
                                a_z_half(j) * value_dsigmazz_dz
            
                    value_dsigmaxz_dx = value_dsigmaxz_dx / K_x_half(i) + memory_dsigmaxz_dx(i,j)
                    value_dsigmazz_dz = value_dsigmazz_dz / K_z_half(j) + memory_dsigmazz_dz(i,j)
            
                    vz(i,j) = vz(i,j) + (value_dsigmaxz_dx + value_dsigmazz_dz) * &
                            ( dt / ( deltarho - dt * gammaz(i,j) )  )
        
                enddo
            enddo
            !$omp end parallel do
            
            ! Add the source term
            vx(isource,jsource) = vx(isource,jsource) + srcx(it) * dt / rho(isource,jsource)
            vz(isource,jsource) = vz(isource,jsource) + srcz(it) * dt / rho(isource,jsource)
        
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
        deallocate(gammax, gammaz, srcx, srcz)
        
        ! Allocate more
        deallocate(memory_dvx_dx, memory_dvx_dz)
        deallocate(memory_dvz_dx, memory_dvz_dz)
        deallocate(memory_dsigmaxx_dx, memory_dsigmazz_dz)
        deallocate(memory_dsigmaxz_dx, memory_dsigmaxz_dz)
        deallocate(vx, vz, sigmaxx, sigmazz, sigmaxz)
        
    end subroutine seismic2
    
    ! subroutine seismic2(nx, nz, dx, dz, npoints_pml, src, nstep, OUTPUT_SINGLE)

    !     ! 2D elastic finite-difference code in velocity and stress formulation
    !     ! with Convolutional-PML (C-PML) absorbing conditions for an 
    !     ! anisotropic medium
    !     !
    !     ! If using this program please give credit to the following: 
    !     !
    !     ! Dimitri Komatitsch, University of Pau, France, April 2007.
    !     ! Anisotropic implementation by Roland Martin and Dimitri Komatitsch, 
    !     ! University of Pau, France, April 2007.
        
    !     ! The second-order staggered-grid formulation of Madariaga (1976) and 
    !     ! Virieux (1986) is used:

    !     ! INPUT
    !     !   im (INTEGER)
    !     !   nx, ny (INTEGER)
    !     !   c11, c12, c22, c66, rho (REAL)
    !     !   dx, dy (REAL)
    !     !   npoints_pml (INTEGER) - the thickness of the pml
    !     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
    !         implicit none
        
    !         integer, parameter :: dp=kind(0.d0)

    !         ! total number of grid points in each direction of the grid
    !         integer, intent(in) :: nx
    !         integer, intent(in) :: nz
        
    !         ! thickness of the PML layer in grid points
    !         integer :: npoints_pml
    !         ! integer, dimension(nx,nz)
    !         real(real64), dimension(nx,nz) :: c11, c13, c15, c33, c35, c55, rho
    !         real(real64) :: deltarho
        
    !         ! total number of time steps
    !         integer :: nstep
        
    !         ! time step in seconds 
    !         real(real64) :: DT
    !         real(real64), intent(in) :: dx, dz 

    !         ! source
    !         integer,dimension(:) :: src
    !         integer :: isource, jsource
        
    !         ! velocity threshold above which we consider that the code became unstable
    !         real(real64), parameter :: STABILITY_THRESHOLD = 1.d+25
        
    !         ! main arrays
    !         real(real64), dimension(nx,nz) :: vx,vz,sigmaxx,sigmazz,sigmaxz
                    
    !         ! arrays for the memory variables
    !         ! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
    !         real(real64), dimension(NX,NZ) :: &
    !             memory_dvx_dx, &
    !             memory_dvx_dz, &
    !             memory_dvz_dx, &
    !             memory_dvz_dz, &
    !             memory_dsigmaxx_dx, &
    !             memory_dsigmazz_dz, &
    !             memory_dsigmaxz_dx, &
    !             memory_dsigmaxz_dz
        
    !         real(real64) :: &
    !             value_dvx_dx, &
    !             value_dvx_dz, &
    !             value_dvz_dx, &
    !             value_dvz_dz, &
    !             value_dsigmaxx_dx, &
    !             value_dsigmazz_dz, &
    !             value_dsigmaxz_dx, &
    !             value_dsigmaxz_dz
        
    !         ! 1D arrays for the damping profiles
    !         real(real64), dimension(nx) :: K_x, alpha_x, a_x, b_x, &
    !                                     K_x_half, alpha_x_half, &
    !                                     a_x_half, b_x_half
    !         real(real64), dimension(nz) :: K_z, alpha_z, a_z, b_z, &
    !                                     K_z_half, alpha_z_half, &
    !                                     a_z_half, b_z_half
        
    !         real(real64), dimension(nx,nz) :: gamma_x, gamma_z
    !         ! for the source
    !         real(real64),dimension(nstep) :: srcx, srcz
        
    !         integer :: i,j,it
        
    !         real(real64) :: velocnorm

    !         ! Boolean flag to save as double precision or single precision 
    !         logical :: SINGLE
    !         logical, intent(in), optional :: OUTPUT_SINGLE 

    !         ! The default data output is single precision unless OUTPUT_SINGLE is 
    !         ! set to .FALSE.
    !         if (present(OUTPUT_SINGLE)) then 
    !             SINGLE = OUTPUT_SINGLE 
    !         else
    !             SINGLE = .TRUE.
    !         endif
        
    !         ! -------------------- Load Stiffness Coefficients --------------------
        
    !         call material_rw('c11.dat', c11, .TRUE.)
    !         call material_rw('c13.dat', c13, .TRUE.)
    !         call material_rw('c15.dat', c15, .TRUE.)
    !         call material_rw('c33.dat', c33, .TRUE.)
    !         call material_rw('c35.dat', c35, .TRUE.)
    !         call material_rw('c55.dat', c55, .TRUE.)
    !         call material_rw('rho.dat', rho, .TRUE.)
            
    !         ! ------------------- Load Attenuation Coefficients --------------------
    !         call material_rw('gamma_x.dat', gamma_x, .TRUE.)
    !         call material_rw('gamma_z.dat', gamma_z, .TRUE.)
            
    !         ! ------------------------ Assign some constants -----------------------
        
    !         isource = src(1)+npoints_pml
    !         jsource = src(2)+npoints_pml
        
    !         DT = minval( (/dx,dz/) )/ &
    !             (sqrt( 3.d0*( maxval( (/ c11/rho, c33/rho /) ) ) ) ) 

    !         ! ================================ LOAD SOURCE ================================
        
    !         call loadsource('seismicsourcex.dat', nstep, srcx)
    !         ! We are using the coordinate names x, Z but the math computes the source in 
    !         ! the x-z plane
    !         call loadsource('seismicsourcez.dat', nstep, srcz)
        
    !         ! -----------------------------------------------------------------------------
    !         !--- define profile of absorption in PML region

    !         ! Initialize PML 
    !         K_x(:) = 1.d0
    !         K_x_half(:) = 1.d0
    !         alpha_x(:) = 0.d0
    !         alpha_x_half(:) = 0.d0
    !         a_x(:) = 0.d0
    !         a_x_half(:) = 0.d0
    !         b_x(:) = 0.d0 
    !         b_x_half(:) = 0.d0

    !         K_z(:) = 1.d0
    !         K_z_half(:) = 1.d0 
    !         alpha_z(:) = 0.d0
    !         alpha_z_half(:) = 0.d0
    !         a_z(:) = 0.d0
    !         a_z_half(:) = 0.d0
    !         b_z(:) = 0.d0
    !         b_z_half(:) = 0.d0

    !     ! ------------------------------ Load the boundary ----------------------------        
    !     call loadcpml('kappax_cpml.dat', K_x)
    !     call loadcpml('alphax_cpml.dat', alpha_x)
    !     call loadcpml('acoefx_cpml.dat', a_x)
    !     call loadcpml('bcoefx_cpml.dat', b_x)
        
    !     call loadcpml('kappaz_cpml.dat', K_z)
    !     call loadcpml('alphaz_cpml.dat', alpha_z)
    !     call loadcpml('acoefz_cpml.dat', a_z)
    !     call loadcpml('bcoefz_cpml.dat', b_z)
        
    !     call loadcpml('kappax_half_cpml.dat', K_x_half)
    !     call loadcpml('alphax_half_cpml.dat', alpha_x_half)
    !     call loadcpml('acoefx_half_cpml.dat', a_x_half)
    !     call loadcpml('bcoefx_half_cpml.dat', b_x_half)
        
    !     call loadcpml('kappaz_half_cpml.dat', K_z_half)
    !     call loadcpml('alphaz_half_cpml.dat', alpha_z_half)
    !     call loadcpml('acoefz_half_cpml.dat', a_z_half)
    !     call loadcpml('bcoefz_half_cpml.dat', b_z_half)

    !     ! =============================================================================


    !     ! initialize arrays
    !     vx(:,:) = 0.d0
    !     vz(:,:) = 0.d0
    !     sigmaxx(:,:) = 0.d0
    !     sigmazz(:,:) = 0.d0
    !     sigmaxz(:,:) = 0.d0

    !     ! PML
    !     memory_dvx_dx(:,:) = 0.d0
    !     memory_dvx_dz(:,:) = 0.d0
    !     memory_dvz_dx(:,:) = 0.d0
    !     memory_dvz_dz(:,:) = 0.d0
    !     memory_dsigmaxx_dx(:,:) = 0.d0
    !     memory_dsigmazz_dz(:,:) = 0.d0
    !     memory_dsigmaxz_dx(:,:) = 0.d0
    !     memory_dsigmaxz_dz(:,:) = 0.d0

    !     !---
    !     !---  beginning of time loop
    !     !---

    !     do it = 1,NSTEP
    !         !------------------------------------------------------------
    !         ! compute stress sigma and update memory variables for C-PML
    !         !------------------------------------------------------------
    !         do j = 2,NZ
    !             do i = 1,NX-1
        
    !                 value_dvx_dx = (vx(i+1,j) - vx(i,j)) / DX
    !                 value_dvz_dz = (vz(i,j) - vz(i,j-1)) / DZ
    !                 value_dvz_dx = (vz(i+1,j) - vz(i,j)) / DX
    !                 value_dvx_dz = (vx(i,j) - vx(i,j-1)) / DZ

    !                 memory_dvx_dx(i,j) = b_x_half(j) * memory_dvx_dx(i,j) + &
    !                                         a_x_half(i) * value_dvx_dx
    !                 memory_dvz_dz(i,j) = b_z(j) * memory_dvz_dz(i,j) + &
    !                                         a_z(j) * value_dvz_dz
    !                 memory_dvx_dz(i,j) = b_z_half(j) * memory_dvx_dz(i,j) + &
    !                                         a_z_half(j) * value_dvx_dz 
    !                 memory_dvz_dx(i,j) = b_x(i) * memory_dvz_dx(i,j) + &
    !                                         a_x(i) * value_dvz_dx

    !                 value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j)
    !                 value_dvz_dz = value_dvz_dz / K_z(j) + memory_dvz_dz(i,j)
    !                 value_dvz_dx = value_dvz_dx / K_x(i) + memory_dvz_dx(i,j)
    !                 value_dvx_dz = value_dvx_dz / K_z_half(j) + memory_dvx_dz(i,j)
                    
    !                 sigmaxx(i,j) = sigmaxx(i,j) + &
    !                     (   c11(i,j) * value_dvx_dx + &
    !                         c13(i,j) * value_dvz_dz + &
    !                         c15(i,j) * (value_dvz_dx + value_dvx_dz) ) * DT
    !                 sigmazz(i,j) = sigmazz(i,j) + &
    !                     (   c13(i,j) * value_dvx_dx + &
    !                         c33(i,j) * value_dvz_dz + &
    !                         c35(i,j) * (value_dvz_dx + value_dvx_dz) ) * DT
        
    !             enddo
    !         enddo
        
    !         do j = 1,NZ-1
    !             do i = 2,NX
        
    !             value_dvx_dx = (vx(i,j) - vx(i-1,j)) / DX
    !             value_dvz_dz = (vz(i,j+1) - vz(i,j)) / DZ
    !             value_dvz_dx = (vz(i,j) - vz(i-1,j)) / DX
    !             value_dvx_dz = (vx(i,j+1) - vx(i,j)) / DZ
                
    !             memory_dvx_dx(i,j) = b_x_half(i) * memory_dvx_dx(i,j) + &
    !                                     a_x_half(i) * value_dvx_dx
    !             memory_dvz_dz(i,j) = b_z(j) * memory_dvz_dz(i,j) + &
    !                                     a_z(j) * value_dvz_dz
    !             memory_dvx_dz(i,j) = b_z_half(j) * memory_dvx_dz(i,j) + &
    !                                     a_z_half(j) * value_dvx_dz 
    !             memory_dvz_dx(i,j) = b_x(i) * memory_dvz_dx(i,j) + &
    !                                     a_x(i) * value_dvz_dx
                
    !             value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j)
    !             value_dvz_dz = value_dvz_dz / K_z(j) + memory_dvz_dz(i,j)
    !             value_dvz_dx = value_dvz_dx / K_x(i) + memory_dvz_dx(i,j)
    !             value_dvx_dz = value_dvx_dz / K_z_half(j) + memory_dvx_dz(i,j)
        
    !             sigmaxz(i,j) = sigmaxz(i,j) + &
    !                 (   c15(i,j)  * value_dvx_dx + & 
    !                     c35(i,j)  * value_dvz_dz + &
    !                     c55(i,j) * (value_dvz_dx + value_dvx_dz) ) * DT
        
    !             enddo
    !         enddo
        
    !         !--------------------------------------------------------
    !         ! compute velocity and update memory variables for C-PML
    !         !--------------------------------------------------------
    !         do j = 2,NZ
    !             do i = 2,NX
        
    !             deltarho = ( rho(i,j) + rho(i,j-1) + rho(i-1,j) + rho(i-1,j-1) )/4

    !             value_dsigmaxx_dx = (sigmaxx(i,j) - sigmaxx(i-1,j)) / DX
    !             value_dsigmaxz_dz = (sigmaxz(i,j) - sigmaxz(i,j-1)) / DZ
        
    !             memory_dsigmaxx_dx(i,j) = b_x(i) * memory_dsigmaxx_dx(i,j) + &
    !                         a_x(i) * value_dsigmaxx_dx
    !             memory_dsigmaxz_dz(i,j) = b_z(j) * memory_dsigmaxz_dz(i,j) + &
    !                         a_z(j) * value_dsigmaxz_dz
        
    !             value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + &
    !                         memory_dsigmaxx_dx(i,j)
    !             value_dsigmaxz_dz = value_dsigmaxz_dz / K_z(j) + &
    !                         memory_dsigmaxz_dz(i,j)
        
    !             vx(i,j) = vx(i,j)*(1 - gamma_x(i,j) ) + (value_dsigmaxx_dx + value_dsigmaxz_dz) * DT / deltarho !rho(i,j)
        
    !             enddo
    !         enddo
        
    !         do j = 1,NZ-1
    !             do i = 1,NX-1
        
    !                 deltarho = ( 2*rho(i,j) + rho(i+1,j) + rho(i,j+1) )/4
    !                 value_dsigmaxz_dx = (sigmaxz(i+1,j) - sigmaxz(i,j)) / DX
    !                 value_dsigmazz_dz = (sigmazz(i,j+1) - sigmazz(i,j)) / DZ
            
    !                 memory_dsigmaxz_dx(i,j) = b_x_half(i) * memory_dsigmaxz_dx(i,j) + &
    !                             a_x_half(i) * value_dsigmaxz_dx
    !                 memory_dsigmazz_dz(i,j) = b_z_half(j) * memory_dsigmazz_dz(i,j) + &
    !                             a_z_half(j) * value_dsigmazz_dz
            
    !                 value_dsigmaxz_dx = value_dsigmaxz_dx / K_x_half(i) + memory_dsigmaxz_dx(i,j)
    !                 value_dsigmazz_dz = value_dsigmazz_dz / K_z_half(j) + memory_dsigmazz_dz(i,j)
            
    !                 vz(i,j) = vz(i,j)*(1 - gamma_z(i,j) ) + (value_dsigmaxz_dx + value_dsigmazz_dz) * DT / deltarho
        
    !             enddo
    !         enddo

    !         ! Add the source term
    !         vx(isource,jsource) = vx(isource,jsource) + srcx(it) * DT / rho(isource,jsource)
    !         vz(isource,jsource) = vz(isource,jsource) + srcz(it) * DT / rho(isource,jsource)
        
    !         ! Dirichlet conditions (rigid boundaries) on the edges or at the 
    !         ! bottom of the PML layers
    !         vx(1,:) = 0.d0
    !         vx(NX,:) = 0.d0
        
    !         vx(:,1) = 0.d0
    !         vx(:,NZ) = 0.d0
        
    !         vz(1,:) = 0.d0
    !         vz(NX,:) = 0.d0
        
    !         vz(:,1) = 0.d0
    !         vz(:,NZ) = 0.d0
        
    !         ! print maximum of norm of velocity
    !         velocnorm = maxval(sqrt(vx**2 + vz**2))
    !         if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'
        
    !         call write_image2(vx, nx, nz, src, it, 'Vx', SINGLE)
    !         call write_image2(vz, nx, nz, src, it, 'Vz', SINGLE)

    !     enddo   ! end of time loop
    ! end subroutine seismic2

    ! =========================================================================
    ! subroutine seismic25(domain, source, seisvar, SINGLE_OUTPUT)
    !     !--------------------------------------------------------------------------------------
    !     use constants
        
    !     implicit none

    !     ! Input arguments
    !     type(Domain_Type), intent(in) :: domain 
    !     type(Source_Type), intent(in) :: source 
    !     type(Seismic3_Variables_Type), intent(inout) :: seisvar
    !     logical, intent(in), optional :: SINGLE_OUTPUT

    !     ! Local variables
    !     real(real64), dimension(domain%nx,domain%nz) :: c11, c12, c13, c14, c15, c16, &
    !                                         c22, c23, c24, c25, c26, &
    !                                         c33, c34, c35, c36, &
    !                                         c44, c45, c46, &
    !                                         c55, c56, &
    !                                         c66, &
    !                                         rho
    !     real(real64) :: deltarho, velocnorm
    !     ! real(real64) :: DT
    !     integer :: i, j, k, it, isource, jsource, ksource

    !     ! Values of the velocity and stress differentials
    !     real(real64) :: dvx_dx, dvx_dy, dvx_dz, &
    !                     dvy_dx, dvy_dy, dvy_dz, &
    !                     dvz_dx, dvz_dy, dvz_dz, &
    !                     dsigmaxx_dx, dsigmayy_dy, dsigmazz_dz, &
    !                     dsigmaxy_dx, dsigmaxy_dy, &
    !                     dsigmaxz_dx, dsigmaxz_dz, &
    !                     dsigmayz_dy, dsigmayz_dz

    !     ! 1D arrays for the damping profiles in each direction
    !     real(real64), dimension(domain%nx) :: K_x, alpha_x, a_x, b_x, K_x_half, & 
    !                                             alpha_x_half, a_x_half, b_x_half
    !     real(real64), dimension(domain%ny) :: K_y, alpha_y, a_y, b_y, K_y_half, & 
    !                                             alpha_y_half, a_y_half, b_y_half
    !     real(real64), dimension(domain%nz) :: K_z, alpha_z, a_z, b_z, K_z_half, & 
    !                                             alpha_z_half, a_z_half, b_z_half

    !     ! Arrays for the PML damping factors
    !     real(real64), dimension(domain%nx,domain%nz) :: gammax, gammay, gammaz

    !     ! Source arrays
    !     real(real64), dimension(source%time_steps) :: srcx, srcy, srcz

    !     ! Boolean flag to save as double precision or single precision
    !     logical :: SINGLE

    !     ! The default data output is single precision unless SINGLE_OUTPUT is 
    !     ! set to .FALSE.
    !     if (present(SINGLE_OUTPUT)) then 
    !         SINGLE = SINGLE_OUTPUT 
    !     else
    !         SINGLE = .TRUE.
    !     endif
        
    !     ! ------------------------ Load Stiffness Coefficients ------------------------
    !     call material_rw2('c11.dat', c11, .TRUE.)
    !     call material_rw2('c12.dat', c12, .TRUE.)
    !     call material_rw2('c13.dat', c13, .TRUE.)
    !     call material_rw2('c14.dat', c14, .TRUE.)
    !     call material_rw2('c15.dat', c15, .TRUE.)
    !     call material_rw2('c16.dat', c16, .TRUE.)
    !     call material_rw2('c22.dat', c22, .TRUE.)
    !     call material_rw2('c23.dat', c23, .TRUE.)
    !     call material_rw2('c24.dat', c24, .TRUE.)
    !     call material_rw2('c25.dat', c25, .TRUE.)
    !     call material_rw2('c26.dat', c26, .TRUE.)
    !     call material_rw2('c33.dat', c33, .TRUE.)
    !     call material_rw2('c34.dat', c34, .TRUE.)
    !     call material_rw2('c35.dat', c35, .TRUE.)
    !     call material_rw2('c36.dat', c36, .TRUE.)
    !     call material_rw2('c44.dat', c44, .TRUE.)
    !     call material_rw2('c45.dat', c45, .TRUE.)
    !     call material_rw2('c46.dat', c46, .TRUE.)
    !     call material_rw2('c55.dat', c55, .TRUE.)
    !     call material_rw2('c56.dat', c56, .TRUE.)
    !     call material_rw2('c66.dat', c66, .TRUE.)
    !     call material_rw2('rho.dat', rho, .TRUE.)
        
    !     ! ------------------- Load Attenuation Coefficients --------------------
    !     call material_rw2('gamma_x.dat', gammax, .TRUE.)
    !     call material_rw2('gamma_z.dat', gammaz, .TRUE.)
    !     call material_rw2('gamma_y.dat', gammay, .TRUE.)
        
    !     ! ------------------------ Assign some constants -----------------------
    !     isource = source%xind + domain%cpml
    !     jsource = source%yind + domain%cpml
    !     ksource = source%zind + domain%cpml

    !     ! ================================ LOAD SOURCE ================================

    !     call loadsource('seismicsourcex.dat', source%time_steps, srcx)
    !     call loadsource('seismicsourcey.dat', source%time_steps, srcy)
    !     call loadsource('seismicsourcez.dat', source%time_steps, srcz)

    !     ! ==================================== PML ====================================
    !     ! Initialize PML 
    !     K_x(:) = 1.d0
    !     K_x_half(:) = 1.d0
    !     alpha_x(:) = 0.d0
    !     alpha_x_half(:) = 0.d0
    !     a_x(:) = 0.d0
    !     a_x_half(:) = 0.d0

    !     K_y(:) = 1.d0
    !     K_y_half(:) = 1.d0
    !     alpha_y(:) = 0.d0
    !     alpha_y_half(:) = 0.d0
    !     a_y(:) = 0.d0
    !     a_y_half(:) = 0.d0

    !     K_z(:) = 1.d0
    !     K_z_half(:) = 1.d0 
    !     alpha_z(:) = 0.d0
    !     alpha_z_half(:) = 0.d0
    !     a_z(:) = 0.d0
    !     a_z_half(:) = 0.d0

    !     ! ------------------------- Boundary Conditions -------------------------
    !     call loadcpml('kappax_cpml.dat', K_x)
    !     call loadcpml('alphax_cpml.dat', alpha_x)
    !     call loadcpml('acoefx_cpml.dat', a_x)
    !     call loadcpml('bcoefx_cpml.dat', b_x)

    !     call loadcpml('kappay_cpml.dat', K_y)
    !     call loadcpml('alphay_cpml.dat', alpha_y)
    !     call loadcpml('acoefy_cpml.dat', a_y)
    !     call loadcpml('bcoefy_cpml.dat', b_y)

    !     call loadcpml('kappaz_cpml.dat', K_z)
    !     call loadcpml('alphaz_cpml.dat', alpha_z)
    !     call loadcpml('acoefz_cpml.dat', a_z)
    !     call loadcpml('bcoefz_cpml.dat', b_z)

    !     call loadcpml('kappax_half_cpml.dat', K_x_half)
    !     call loadcpml('alphax_half_cpml.dat', alpha_x_half)
    !     call loadcpml('acoefx_half_cpml.dat', a_x_half)
    !     call loadcpml('bcoefx_half_cpml.dat', b_x_half)

    !     call loadcpml('kappay_half_cpml.dat', K_y_half)
    !     call loadcpml('alphay_half_cpml.dat', alpha_y_half)
    !     call loadcpml('acoefy_half_cpml.dat', a_y_half)
    !     call loadcpml('bcoefy_half_cpml.dat', b_y_half)

    !     call loadcpml('kappaz_half_cpml.dat', K_z_half)
    !     call loadcpml('alphaz_half_cpml.dat', alpha_z_half)
    !     call loadcpml('acoefz_half_cpml.dat', a_z_half)
    !     call loadcpml('bcoefz_half_cpml.dat', b_z_half)

    !     ! =============================== Forward Model ===============================
    !     ! Load initial condition
    !     call material_rw3('initialconditionVx.dat', seisvar%vx, .TRUE.)
    !     call material_rw3('initialconditionVy.dat', seisvar%vy, .TRUE.)
    !     call material_rw3('initialconditionVz.dat', seisvar%vz, .TRUE.)

    !     ! Initialize the stress values
    !     seisvar%sigxx(:,:,:) = 0.d0
    !     seisvar%sigyy(:,:,:) = 0.d0
    !     seisvar%sigzz(:,:,:) = 0.d0
    !     seisvar%sigxy(:,:,:) = 0.d0
    !     seisvar%sigxz(:,:,:) = 0.d0
    !     seisvar%sigyz(:,:,:) = 0.d0

    !     ! PML
    !     seisvar%memdvx_dx(:,:,:) = 0.d0
    !     seisvar%memdvx_dy(:,:,:) = 0.d0
    !     seisvar%memdvx_dz(:,:,:) = 0.d0

    !     seisvar%memdvy_dx(:,:,:) = 0.d0
    !     seisvar%memdvy_dy(:,:,:) = 0.d0
    !     seisvar%memdvy_dz(:,:,:) = 0.d0

    !     seisvar%memdvz_dx(:,:,:) = 0.d0
    !     seisvar%memdvz_dy(:,:,:) = 0.d0 
    !     seisvar%memdvz_dz(:,:,:) = 0.d0

    !     seisvar%memdsigxx_dx(:,:,:) = 0.d0
    !     seisvar%memdsigyy_dy(:,:,:) = 0.d0
    !     seisvar%memdsigzz_dz(:,:,:) = 0.d0

    !     seisvar%memdsigxy_dx(:,:,:) = 0.d0
    !     seisvar%memdsigxy_dy(:,:,:) = 0.d0
    !     seisvar%memdsigxz_dx(:,:,:) = 0.d0
    !     seisvar%memdsigxz_dz(:,:,:) = 0.d0
    !     seisvar%memdsigyz_dy(:,:,:) = 0.d0
    !     seisvar%memdsigyz_dz(:,:,:) = 0.d0

    !     ! Do it 
    !     do it = 1,source%time_steps
    !         !------------------------------------------------------------
    !         ! compute stress sigma and update memory variables for C-PML
    !         !------------------------------------------------------------
    !         ! Update in the x direction
    !         do k = 2,domain%nz
    !             do j = 2,domain%ny
    !                 do i = 1,domain%nx-1

    !                     dvx_dx = (seisvar%vx(i+1,j,k) - seisvar%vx(i,j,k) ) / domain%dx
    !                     dvy_dx = (seisvar%vy(i+1,j,k) - seisvar%vy(i,j,k) ) / domain%dx
    !                     dvz_dx = (seisvar%vz(i+1,j,k) - seisvar%vz(i,j,k) ) / domain%dx 
    !                     dvy_dy = (seisvar%vy(i,j,k) - seisvar%vy(i,j-1,k) ) / domain%dy
    !                     dvx_dy = (seisvar%vx(i,j,k) - seisvar%vx(i,j-1,k) ) / domain%dy
    !                     dvz_dy = (seisvar%vz(i,j,k) - seisvar%vz(i,j-1,k) ) / domain%dy
    !                     dvz_dz = (seisvar%vz(i,j,k) - seisvar%vz(i,j,k-1) ) / domain%dz
    !                     dvx_dz = (seisvar%vx(i,j,k) - seisvar%vx(i,j,k-1) ) / domain%dz
    !                     dvy_dz = (seisvar%vy(i,j,k) - seisvar%vy(i,j,k-1) ) / domain%dz

    !                     seisvar%memdvx_dx(i,j,k) = b_x_half(i) * seisvar%memdvx_dx(i,j,k) + a_x_half(i) * dvx_dx
    !                     seisvar%memdvy_dx(i,j,k) = b_x_half(i) * seisvar%memdvy_dx(i,j,k) + a_x_half(i) * dvy_dx
    !                     seisvar%memdvz_dx(i,j,k) = b_x_half(i) * seisvar%memdvz_dx(i,j,k) + a_x_half(i) * dvz_dx
    !                     seisvar%memdvy_dy(i,j,k) = b_y(j) * seisvar%memdvy_dy(i,j,k) + a_y(j) * dvy_dy
    !                     seisvar%memdvx_dy(i,j,k) = b_y(j) * seisvar%memdvx_dy(i,j,k) + a_y(j) * dvx_dy
    !                     seisvar%memdvz_dy(i,j,k) = b_y(j) * seisvar%memdvz_dy(i,j,k) + a_y(j) * dvz_dy
    !                     seisvar%memdvz_dz(i,j,k) = b_z(k) * seisvar%memdvz_dz(i,j,k) + a_z(k) * dvz_dz
    !                     seisvar%memdvx_dz(i,j,k) = b_z(k) * seisvar%memdvx_dz(i,j,k) + a_z(k) * dvx_dz
    !                     seisvar%memdvy_dz(i,j,k) = b_z(k) * seisvar%memdvy_dz(i,j,k) + a_z(k) * dvy_dz

    !                     dvx_dx = dvx_dx / K_x_half(i) + seisvar%memdvx_dx(i,j,k)
    !                     dvy_dx = dvy_dx / K_x_half(i) + seisvar%memdvy_dx(i,j,k)
    !                     dvz_dx = dvz_dx / K_x_half(i) + seisvar%memdvz_dx(i,j,k)
    !                     dvy_dy = dvy_dy / K_y(j) + seisvar%memdvy_dy(i,j,k)
    !                     dvx_dy = dvx_dy / K_y(j) + seisvar%memdvx_dy(i,j,k)
    !                     dvz_dy = dvz_dy / K_y(j) + seisvar%memdvz_dy(i,j,k)
    !                     dvz_dz = dvz_dz / K_z(k) + seisvar%memdvz_dz(i,j,k)
    !                     dvx_dz = dvx_dz / K_z(k) + seisvar%memdvx_dz(i,j,k)
    !                     dvy_dz = dvy_dz / K_z(k) + seisvar%memdvy_dz(i,j,k)

    !                     seisvar%sigxx(i,j,k) = seisvar%sigxx(i,j,k) + &
    !                     (   c11(i,k) * dvx_dx + c12(i,k) * dvy_dy + c13(i,k) * dvz_dz + &
    !                         c14(i,k) * (dvy_dz + dvz_dy) + c15(i,k) * (dvx_dz + dvz_dx) + &
    !                         c16(i,k) * (dvx_dy + dvz_dy) ) * source%dt

    !                     ! Full 3D will need a gradient in the y-direction
    !                     seisvar%sigyy(i,j,k) = seisvar%sigyy(i,j,k) + &
    !                     (   c12(i,k) * dvx_dx + c22(i,k) * dvy_dy + c23(i,k) * dvz_dz + &
    !                         c24(i,k) * (dvy_dz + dvz_dy) + c25(i,k) * (dvx_dz + dvz_dx) + &
    !                         c26(i,k) * (dvy_dx + dvx_dy) ) * source%dt

    !                     seisvar%sigzz(i,j,k) = seisvar%sigzz(i,j,k) + &
    !                     (   c13(i,k) * dvx_dx + c23(i,k) * dvy_dy + c33(i,k) * dvz_dz + &
    !                         c34(i,k) * (dvy_dz + dvz_dy) + c35(i,k) * (dvx_dz + dvz_dx) + &
    !                         c36(i,k) * (dvy_dx + dvx_dy) ) * source%dt

    !                 enddo
    !             enddo
    !         enddo

    !         ! Update sigmaxy, x-direction is full nodes
    !         do k = 2,domain%nz
    !             do j = 1,domain%ny-1
    !                 do i = 2,domain%nx

    !                     dvx_dx = (seisvar%vx(i,j,k) - seisvar%vx(i-1,j,k)) / domain%dx
    !                     dvy_dx = (seisvar%vy(i,j,k) - seisvar%vy(i-1,j,k)) / domain%dx
    !                     dvz_dx = (seisvar%vz(i,j,k) - seisvar%vz(i-1,j,k)) / domain%dx
    !                     dvy_dy = (seisvar%vy(i,j+1,k) - seisvar%vy(i,j,k)) / domain%dy
    !                     dvx_dy = (seisvar%vx(i,j+1,k) - seisvar%vx(i,j,k)) / domain%dy
    !                     dvz_dy = (seisvar%vz(i,j+1,k) - seisvar%vz(i,j,k)) / domain%dy
    !                     dvz_dz = (seisvar%vz(i,j,k) - seisvar%vz(i,j,k-1)) / domain%dz
    !                     dvx_dz = (seisvar%vx(i,j,k) - seisvar%vx(i,j,k-1)) / domain%dz
    !                     dvy_dz = (seisvar%vy(i,j,k) - seisvar%vy(i,j,k-1)) / domain%dz

    !                     seisvar%memdvx_dx(i,j,k) = b_x(i) * seisvar%memdvx_dx(i,j,k) + a_x(i) * dvx_dx
    !                     seisvar%memdvy_dx(i,j,k) = b_x(i) * seisvar%memdvy_dx(i,j,k) + a_x(i) * dvy_dx
    !                     seisvar%memdvz_dx(i,j,k) = b_x(i) * seisvar%memdvz_dx(i,j,k) + a_x(i) * dvz_dx
    !                     seisvar%memdvy_dy(i,j,k) = b_y_half(j) * seisvar%memdvy_dy(i,j,k) + a_y_half(j) * dvy_dy
    !                     seisvar%memdvx_dy(i,j,k) = b_y_half(j) * seisvar%memdvx_dy(i,j,k) + a_y_half(j) * dvx_dy
    !                     seisvar%memdvz_dy(i,j,k) = b_y_half(j) * seisvar%memdvz_dy(i,j,k) + a_y_half(j) * dvz_dy
    !                     seisvar%memdvz_dz(i,j,k) = b_z_half(k) * seisvar%memdvz_dz(i,j,k) + a_z_half(k) * dvz_dz
    !                     seisvar%memdvx_dz(i,j,k) = b_z_half(k) * seisvar%memdvx_dz(i,j,k) + a_z_half(k) * dvx_dz
    !                     seisvar%memdvy_dz(i,j,k) = b_z_half(k) * seisvar%memdvy_dz(i,j,k) + a_z_half(k) * dvy_dz

    !                     dvx_dx = dvx_dx / K_x(i) + seisvar%memdvx_dx(i,j,k)
    !                     dvy_dx = dvy_dx / K_x(i) + seisvar%memdvy_dx(i,j,k)
    !                     dvy_dy = dvy_dy / K_y_half(j) + seisvar%memdvy_dy(i,j,k)
    !                     dvx_dy = dvx_dy / K_y_half(j) + seisvar%memdvx_dy(i,j,k)
    !                     dvz_dy = dvz_dy / K_y_half(j) + seisvar%memdvz_dy(i,j,k)
    !                     dvz_dz = dvz_dz / K_z_half(k) + seisvar%memdvz_dz(i,j,k)
    !                     dvy_dz = dvy_dz / K_z_half(k) + seisvar%memdvy_dz(i,j,k)

    !                     seisvar%sigxy(i,j,k) = seisvar%sigxy(i,j,k) + &
    !                     (   c16(i,k) * dvx_dx + c26(i,k) * dvy_dy + c36(i,k) * dvz_dz + &
    !                         c46(i,k) * (dvz_dy + dvy_dz) + c56(i,k) * (dvz_dx + dvx_dz) + &
    !                         c66(i,k) * (dvy_dx + dvx_dy) ) * source%dt

    !                 enddo
    !             enddo
    !         enddo

    !         ! Update sigmaxz, z-direction is full nodes
    !         do k = 1,domain%nz-1
    !             do j = 2,domain%ny
    !                 do i = 2,domain%nx

    !                     dvx_dx = (seisvar%vx(i,j,k) - seisvar%vx(i-1,j,k)) / domain%dx
    !                     dvy_dx = (seisvar%vy(i,j,k) - seisvar%vy(i-1,j,k)) / domain%dx
    !                     dvz_dx = (seisvar%vz(i,j,k) - seisvar%vz(i-1,j,k)) / domain%dx
    !                     dvy_dy = (seisvar%vy(i,j,k) - seisvar%vy(i,j-1,k)) / domain%dy
    !                     dvz_dy = (seisvar%vz(i,j,k) - seisvar%vz(i,j-1,k)) / domain%dy
    !                     dvx_dy = (seisvar%vx(i,j,k) - seisvar%vx(i,j-1,k)) / domain%dy
    !                     dvz_dz = (seisvar%vz(i,j,k+1) - seisvar%vz(i,j,k)) / domain%dz
    !                     dvx_dz = (seisvar%vx(i,j,k+1) - seisvar%vx(i,j,k)) / domain%dz
    !                     dvy_dz = (seisvar%vy(i,j,k+1) - seisvar%vy(i,j,k)) / domain%dz

    !                     seisvar%memdvx_dx(i,j,k) = b_x(i) * seisvar%memdvx_dx(i,j,k) + a_x(i) * dvx_dx
    !                     seisvar%memdvy_dx(i,j,k) = b_x(i) * seisvar%memdvy_dx(i,j,k) + a_x(i) * dvy_dx
    !                     seisvar%memdvz_dx(i,j,k) = b_x(i) * seisvar%memdvz_dx(i,j,k) + a_x(i) * dvz_dx
    !                     seisvar%memdvy_dy(i,j,k) = b_y(j) * seisvar%memdvy_dy(i,j,k) + a_y(j) * dvy_dy
    !                     seisvar%memdvx_dy(i,j,k) = b_y(j) * seisvar%memdvx_dy(i,j,k) + a_y(j) * dvx_dy
    !                     seisvar%memdvz_dy(i,j,k) = b_y(j) * seisvar%memdvz_dy(i,j,k) + a_y(j) * dvz_dy
    !                     seisvar%memdvz_dz(i,j,k) = b_z_half(k) * seisvar%memdvz_dz(i,j,k) + a_z_half(k) * dvz_dz
    !                     seisvar%memdvx_dz(i,j,k) = b_z_half(k) * seisvar%memdvx_dz(i,j,k) + a_z_half(k) * dvx_dz
    !                     seisvar%memdvy_dz(i,j,k) = b_z_half(k) * seisvar%memdvy_dz(i,j,k) + a_z_half(k) * dvy_dz

    !                     dvx_dx = dvx_dx / K_x(i) + seisvar%memdvx_dx(i,j,k)
    !                     dvy_dx = dvy_dx / K_x(i) + seisvar%memdvy_dx(i,j,k)
    !                     dvz_dx = dvz_dx / K_x(i) + seisvar%memdvz_dx(i,j,k) 
    !                     dvy_dy = dvy_dy / K_y(j) + seisvar%memdvy_dy(i,j,k)
    !                     dvx_dy = dvx_dy / K_y(j) + seisvar%memdvx_dy(i,j,k)
    !                     dvz_dy = dvz_dy / K_y(j) + seisvar%memdvz_dy(i,j,k)
    !                     dvz_dz = dvz_dz / K_z_half(k) + seisvar%memdvz_dz(i,j,k)
    !                     dvx_dz = dvx_dz / K_z_half(k) + seisvar%memdvx_dz(i,j,k)
    !                     dvy_dz = dvy_dz / K_z_half(k) + seisvar%memdvy_dz(i,j,k)

    !                     seisvar%sigxz(i,j,k) = seisvar%sigxz(i,j,k) + &
    !                         (   c15(i,k) * dvx_dx + c25(i,k) * dvy_dy + c35(i,k) * dvz_dz + &
    !                             c45(i,k) * ( dvx_dz + dvz_dx) + c55(i,k) * ( dvx_dz + dvz_dx) + &
    !                             c56(i,k) * ( dvx_dy + dvy_dx) ) * source%dt 

    !                 enddo
    !             enddo

    !             !   ! update sigmayz, y-direction is full nodes
    !             do j = 1,domain%ny-1
    !                 do i = 1,domain%nx-1

    !                     dvx_dx = (seisvar%vx(i+1,j,k) - seisvar%vx(i,j,k)) / domain%DX
    !                     dvy_dx = (seisvar%vy(i+1,j,k) - seisvar%vy(i,j,k)) / domain%DX
    !                     dvz_dx = (seisvar%vz(i+1,j,k) - seisvar%vz(i,j,k)) / domain%DX
    !                     dvy_dy = (seisvar%vy(i,j+1,k) - seisvar%vy(i,j,k)) / domain%DY
    !                     dvx_dy = (seisvar%vx(i,j+1,k) - seisvar%vx(i,j,k)) / domain%DY
    !                     dvz_dy = (seisvar%vz(i,j+1,k) - seisvar%vz(i,j,k)) / domain%DY
    !                     dvz_dz = (seisvar%vz(i,j,k+1) - seisvar%vz(i,j,k)) / domain%DZ
    !                     dvx_dz = (seisvar%vx(i,j,k+1) - seisvar%vx(i,j,k)) / domain%DZ 
    !                     dvy_dz = (seisvar%vy(i,j,k+1) - seisvar%vy(i,j,k)) / domain%DZ

    !                     seisvar%memdvx_dx(i,j,k) = b_x_half(i) * seisvar%memdvx_dx(i,j,k) + a_x_half(i) * dvx_dx
    !                     seisvar%memdvy_dx(i,j,k) = b_x_half(i) * seisvar%memdvy_dx(i,j,k) + a_x_half(i) * dvy_dx
    !                     seisvar%memdvz_dx(i,j,k) = b_x_half(i) * seisvar%memdvz_dx(i,j,k) + a_x_half(i) * dvz_dx
    !                     seisvar%memdvy_dy(i,j,k) = b_y_half(j) * seisvar%memdvy_dy(i,j,k) + a_y_half(j) * dvy_dy
    !                     seisvar%memdvx_dy(i,j,k) = b_y_half(j) * seisvar%memdvx_dy(i,j,k) + a_y_half(j) * dvx_dy
    !                     seisvar%memdvz_dy(i,j,k) = b_y_half(j) * seisvar%memdvz_dy(i,j,k) + a_y_half(j) * dvz_dy
    !                     seisvar%memdvz_dz(i,j,k) = b_z_half(k) * seisvar%memdvz_dz(i,j,k) + a_z_half(k) * dvz_dz
    !                     seisvar%memdvx_dz(i,j,k) = b_z_half(k) * seisvar%memdvx_dz(i,j,k) + a_z_half(k) * dvx_dz
    !                     seisvar%memdvy_dz(i,j,k) = b_z_half(k) * seisvar%memdvy_dz(i,j,k) + a_z_half(k) * dvy_dz

    !                     dvx_dx = dvx_dx / K_x_half(i) + seisvar%memdvx_dx(i,j,k)
    !                     dvy_dx = dvy_dx / K_x_half(i) + seisvar%memdvy_dx(i,j,k)
    !                     dvz_dx = dvz_dx / K_x_half(i) + seisvar%memdvz_dx(i,j,k)
    !                     dvy_dy = dvy_dy / K_y_half(j) + seisvar%memdvy_dy(i,j,k)
    !                     dvx_dy = dvx_dy / K_y_half(j) + seisvar%memdvx_dy(i,j,k)
    !                     dvz_dy = dvz_dy / K_y_half(j) + seisvar%memdvz_dy(i,j,k)
    !                     dvz_dz = dvz_dz / K_z_half(k) + seisvar%memdvz_dz(i,j,k)
    !                     dvx_dz = dvx_dz / K_z_half(k) + seisvar%memdvx_dz(i,j,k)
    !                     dvy_dz = dvy_dz / K_z_half(k) + seisvar%memdvy_dz(i,j,k)

    !                     seisvar%sigyz(i,j,k) = seisvar%sigyz(i,j,k)  + &
    !                         (   c14(i,k) * dvx_dx + c24(i,k) * dvy_dy + c34(i,k) * dvz_dz + &
    !                             c44(i,k) * ( dvy_dz + dvz_dy) + c45(i,k) * ( dvx_dz + dvz_dx) + &
    !                             c46(i,k) * ( dvy_dx + dvx_dy) ) * source%dt 
    !                 enddo
    !             enddo
    !         enddo

    !         !--------------------------------------------------------
    !         ! compute velocity and update memory variables for C-PML
    !         !--------------------------------------------------------
    !         do k = 2,domain%nz
    !             do j = 2,domain%ny
    !                 do i = 2,domain%nx
    !                     ! ds1/dx, ds6/dy, ds5,dz
    !                     deltarho = (4 * rho(i,k) + rho(i-1,k) + rho(i,k-1) )/6

    !                     dsigmaxx_dx = (seisvar%sigxx(i,j,k) - seisvar%sigxx(i-1,j,k) ) / domain%dx
    !                     dsigmaxy_dy = (seisvar%sigxy(i,j,k) - seisvar%sigxy(i,j-1,k) ) / domain%dy
    !                     dsigmaxz_dz = (seisvar%sigxz(i,j,k) - seisvar%sigxz(i,j,k-1) ) / domain%dz

    !                     seisvar%memdsigxx_dx(i,j,k) = b_x(i) * &
    !                         seisvar%memdsigxx_dx(i,j,k) + a_x(i) * dsigmaxx_dx
    !                     seisvar%memdsigxy_dy(i,j,k) = b_y(j) * &
    !                         seisvar%memdsigxy_dy(i,j,k) + a_y(j) * dsigmaxy_dy
    !                     seisvar%memdsigxz_dz(i,j,k) = b_z(k) * &
    !                         seisvar%memdsigxz_dz(i,j,k) + a_z(k) * dsigmaxz_dz

    !                     dsigmaxx_dx = dsigmaxx_dx / K_x(i) + seisvar%memdsigxx_dx(i,j,k)
    !                     dsigmaxy_dy = dsigmaxy_dy / K_y(j) + seisvar%memdsigxy_dy(i,j,k)
    !                     dsigmaxz_dz = dsigmaxz_dz / K_z(k) + seisvar%memdsigxz_dz(i,j,k) 

    !                     seisvar%vx(i,j,k) = seisvar%vx(i,j,k) * (1 - gammax(i,j) ) + &
    !                         (dsigmaxx_dx + dsigmaxy_dy + dsigmaxz_dz) * &
    !                         source%dt / deltarho !rho(i,k)
    !                 enddo
    !             enddo

    !             do j = 1,domain%ny-1
    !                 do i = 1,domain%nx-1
    !                     ! ds6/dx, ds2/dy, ds4/dz
    !                     deltarho = (4*rho(i,k) + rho(i+1,k) + rho(i,k-1) )/6

    !                     dsigmaxy_dx = ( seisvar%sigxy(i+1,j,k) - seisvar%sigxy(i,j,k) ) / domain%dx
    !                     dsigmayy_dy = ( seisvar%sigyy(i,j+1,k) - seisvar%sigyy(i,j,k) ) / domain%dy
    !                     dsigmayz_dz = ( seisvar%sigyz(i,j,k) - seisvar%sigyz(i,j,k-1) ) / domain%dz

    !                     seisvar%memdsigxy_dx(i,j,k) = b_x_half(i) * seisvar%memdsigxy_dx(i,j,k) + a_x_half(i) * dsigmaxy_dx
    !                     seisvar%memdsigyy_dy(i,j,k) = b_y_half(j) * seisvar%memdsigyy_dy(i,j,k) + a_y_half(j) * dsigmayy_dy
    !                     seisvar%memdsigyz_dz(i,j,k) = b_z(k) * seisvar%memdsigyz_dz(i,j,k) + a_z(k) * dsigmayz_dz

    !                     dsigmaxy_dx = dsigmaxy_dx / K_x_half(i) + seisvar%memdsigxy_dx(i,j,k)
    !                     dsigmayy_dy = dsigmayy_dy / K_y_half(j) + seisvar%memdsigyy_dy(i,j,k)
    !                     dsigmayz_dz = dsigmayz_dz / K_z(k) + seisvar%memdsigyz_dz(i,j,k)

    !                     seisvar%vy(i,j,k) = seisvar%vy(i,j,k) * (1 - gammay(i,j) )+ &
    !                         (dsigmaxy_dx + dsigmayy_dy + dsigmayz_dz) * &
    !                         source%dt / deltarho !rho(i,k)
    !                 enddo
    !             enddo
    !         enddo

    !         do k = 1,domain%nz-1
    !             do j = 2,domain%ny
    !                 do i = 1,domain%nx-1
    !                     ! ds5/dx, ds4/dy, ds3/dz
    !                     deltarho = ( rho(i+1,k) + rho(i,k+1) + 4*rho(i,k) )/6

    !                     dsigmaxz_dx = ( seisvar%sigxz(i+1,j,k) - seisvar%sigxz(i,j,k) ) / domain%dx
    !                     dsigmayz_dy = ( seisvar%sigyz(i,j,k) - seisvar%sigyz(i,j-1,k) ) / domain%dy
    !                     dsigmazz_dz = ( seisvar%sigzz(i,j,k+1) - seisvar%sigzz(i,j,k) ) / domain%dz

    !                     seisvar%memdsigxz_dx(i,j,k) = b_x_half(i) * seisvar%memdsigxz_dx(i,j,k) + a_x_half(i) * dsigmaxz_dx
    !                     seisvar%memdsigyz_dy(i,j,k) = b_y(j) * seisvar%memdsigyz_dy(i,j,k) + a_y(j) * dsigmayz_dy
    !                     seisvar%memdsigzz_dz(i,j,k) = b_z_half(k) * seisvar%memdsigzz_dz(i,j,k) + a_z_half(k) * dsigmazz_dz

    !                     dsigmaxz_dx = dsigmaxz_dx / K_x_half(i) + seisvar%memdsigxz_dx(i,j,k)
    !                     dsigmayz_dy = dsigmayz_dy / K_y(j) + seisvar%memdsigyz_dy(i,j,k)
    !                     dsigmazz_dz = dsigmazz_dz / K_z_half(k) + seisvar%memdsigzz_dz(i,j,k)

    !                     seisvar%vz(i,j,k) = seisvar%vz(i,j,k) * (1 - gammaz(i,j) )+ &
    !                         (dsigmaxz_dx + dsigmayz_dy + dsigmazz_dz) * &
    !                         source%dt / deltarho !rho(i,k)

    !                 enddo
    !             enddo
    !         enddo

    !         seisvar%vx(isource,jsource,ksource) = seisvar%vx(isource,jsource,ksource) + &
    !                 srcx(it) * source%dt / rho(isource,ksource)
    !         seisvar%vy(isource,jsource,ksource) = seisvar%vy(isource,jsource,ksource) + &
    !                 srcy(it) * source%dt / rho(isource,ksource)
    !         seisvar%vz(isource,jsource,ksource) = seisvar%vz(isource,jsource,ksource) + &
    !                 srcz(it) * source%dt / rho(isource,ksource)

    !         ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
    !         seisvar%vx(1,:,:) = 0.d0
    !         seisvar%vy(1,:,:) = 0.d0
    !         seisvar%vz(1,:,:) = 0.d0

    !         seisvar%vx(:,1,:) = 0.d0
    !         seisvar%vy(:,1,:) = 0.d0
    !         seisvar%vz(:,1,:) = 0.d0

    !         seisvar%vx(:,:,1) = 0.d0
    !         seisvar%vy(:,:,1) = 0.d0
    !         seisvar%vz(:,:,1) = 0.d0

    !         seisvar%vx(domain%nx,:,:) = 0.d0
    !         seisvar%vy(domain%nx,:,:) = 0.d0
    !         seisvar%vz(domain%nx,:,:) = 0.d0

    !         seisvar%vx(:,domain%ny,:) = 0.d0
    !         seisvar%vy(:,domain%ny,:) = 0.d0
    !         seisvar%vz(:,domain%ny,:) = 0.d0

    !         seisvar%vx(:,:,domain%nz) = 0.d0
    !         seisvar%vy(:,:,domain%nz) = 0.d0
    !         seisvar%vz(:,:,domain%nz) = 0.d0

    !         ! check norm of velocity to make sure the solution isn't diverging
    !         velocnorm = maxval( sqrt(seisvar%vx**2 + seisvar%vy**2 + seisvar%vz**2) )
    !         ! print *,'Time step # ',it,' out of ',time_step
    !         ! print *,'Time: ',(it-1)*DT,' seconds'
    !         ! print *,'Max vals for vx, vy, vz: ', maxval(vx), maxval(vy), maxval(vz)

    !         if (velocnorm > stability_threshold) stop 'code became unstable and blew up'

    !         ! Write the velocity values to an unformatted binary file
    !         call write_image3(seisvar%vx, domain%nx, domain%ny, domain%nz, source, it, 'Vx', SINGLE)
    !         call write_image3(seisvar%vy, domain%nx, domain%ny, domain%nz, source, it, 'Vy', SINGLE)
    !         call write_image3(seisvar%vz, domain%nx, domain%ny, domain%nz, source, it, 'Vz', SINGLE)
    !         ! Now write the stress Values
    !         ! call write_image3(sigmaxx, nx, ny, nz, it, 'S1')
    !         ! call write_image3(sigmayy, nx, ny, nz, it, 'S2')
    !         ! call write_image3(sigmazz, nx, ny, nz, it, 'S3')
    !         ! call write_image3(sigmaxy, nx, ny, nz, it, 'S6')
    !         ! call write_image3(sigmayz, nx, ny, nz, it, 'S4')
    !         ! call write_image3(sigmaxz, nx, ny, nz, it, 'S5')

    !     enddo   ! end of time loop
    ! end subroutine seismic25

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
        real(real64), dimension(domain%nx,domain%nz) :: caEx, cbEx
        real(real64), dimension(domain%nx,domain%nz) :: caEz, cbEz
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

        ! Boolean flag to save as double precision or single precision
        logical :: SINGLE

        ! Check if SINGLE_OUTPUT is provided, default to single precision if not
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        endif
        
        ! ----------------------------------------------------------------------
        allocate(eps11(domain%nx, domain%nz), eps13(domain%nx, domain%nz),  &
                    eps33(domain%nx, domain%nz))
        allocate(sig11(domain%nx, domain%nz), sig13(domain%nx, domain%nz),  &
                    sig33(domain%nx, domain%nz))
        allocate(K_x(domain%nx), alpha_x(domain%nx), a_x(domain%nx), b_x(domain%nx), &
                K_x_half(domain%nx), alpha_x_half(domain%nx), a_x_half(domain%nx), b_x_half(domain%nx))
        allocate(K_z(domain%nz), alpha_z(domain%nz), a_z(domain%nz), b_z(domain%nz), &
                K_z_half(domain%nz), alpha_z_half(domain%nz), a_z_half(domain%nz), b_z_half(domain%nz))
        allocate(srcx(source%time_steps), srcz(source%time_steps))
        
        ! Allocate more
        allocate(epsilonx(domain%nx, domain%nz), epsilonz(domain%nx, domain%nz))
        allocate(sigmax(domain%nx, domain%nz), sigmaz(domain%nx, domain%nz))
        allocate(memory_dEz_dx(domain%nx, domain%nz), memory_dEx_dz(domain%nx, domain%nz))
        allocate(memory_dHy_dx(domain%nx, domain%nz), memory_dHy_dz(domain%nx, domain%nz))
        allocate(Ex(domain%nx, domain%nz), Ez(domain%nx, domain%nz))
        allocate(Hy(domain%nx, domain%nz))
            
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
        call material_rw2('initialconditionHy.dat', Hy, .TRUE.)
        call material_rw2('initialconditionEz.dat', Ez, .TRUE.)

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
        caEx(:,:) = ( 1.0d0 - sigmax * source%dt / &
                    ( 2.0d0 * epsilonx ) ) / &
                    ( 1.0d0 + sigmax * source%dt / &
                    (2.0d0 * epsilonx ) )
        cbEx(:,:) = (source%dt / epsilonx ) / &
                    ( 1.0d0 + sigmax * source%dt / &
                    ( 2.0d0 * epsilonx ) )

        caEz(:,:) = ( 1.0d0 - sigmaz * source%dt / &
                    ( 2.0d0 * epsilonz ) ) / &
                    ( 1.0d0 + sigmaz * source%dt / &
                    (2.0d0 * epsilonz ) )
        cbEz(:,:) = (source%dt / epsilonz ) / &
                    ( 1.0d0 + sigmaz * source%dt / &
                    ( 2.0d0 * epsilonz ) )
        daHy = source%dt/(4.0d0*mu0*mu)
        dbHy = source%dt/mu0 !dt/(mu*mu*dx*(1+daHy) ) 
        daHy = 1.0d0 ! (1-daHy)/(1+daHy) ! 
        
        !---
        !---  beginning of time loop
        !---
        do it = 1,source%time_steps
            !--------------------------------------------------------
            ! compute magnetic field and update memory variables for C-PML
            !--------------------------------------------------------
            do j = 1,domain%nz-1  
                do i = 1,domain%nx-1
                
                    ! Values needed for the magnetic field updates
                    value_dEx_dz = ( Ex(i,j+1) - Ex(i,j) )/domain%dz
                    memory_dEx_dz(i,j) = b_z(j) * memory_dEx_dz(i,j) + a_z(j) * value_dEx_dz
                    value_dEx_dz = value_dEx_dz/ K_z(j) + memory_dEx_dz(i,j)

                    ! The rest of the equation needed for agnetic field updates
                    value_dEz_dx = ( Ez(i+1,j) - Ez(i,j) )/domain%dx
                    memory_dEz_dx(i,j) = b_x(i) * memory_dEz_dx(i,j) + a_x(i) * value_dEz_dx
                    value_dEz_dx = value_dEz_dx/ K_x(i) + memory_dEz_dx(i,j)

                    ! Now update the Magnetic field
                    Hy(i,j) = daHy*Hy(i,j) + dbHy*( value_dEz_dx + value_dEx_dz )

                enddo  
            enddo

            !--------------------------------------------------------
            ! compute electric field and update memory variables for C-PML
            !--------------------------------------------------------
            ! Compute the differences in the y-direction
            do j = 2,domain%nz
                do i = 1,domain%nx
                    ! Update the Ex field
                    value_dHy_dz = ( Hy(i,j) - Hy(i,j-1) )/domain%dz ! this is nz-1 length vector
                    memory_dHy_dz(i,j) = b_z(j) * memory_dHy_dz(i,j) + a_z(j) * value_dHy_dz
                    value_dHy_dz = value_dHy_dz/K_z(j) + memory_dHy_dz(i,j)

                    ! Ex(i,j) = (( caEx(i,j) + caEx(i,j-1) )/2) * Ex(i,j) + &
                    !     (( cbEx(i,j) + cbEx(i,j-1) )/2 ) * value_dHy_dz
                    Ex(i,j) = caEx(i,j) * Ex(i,j) + cbEx(i,j) * value_dHy_dz
                enddo
            enddo

            do j = 1,domain%nz
                do i = 2,domain%nx
                    ! Update the Ez field
                    value_dHy_dx = ( Hy(i,j) - Hy(i-1,j) )/domain%dx
                    memory_dHy_dx(i,j) = b_x_half(i) * memory_dHy_dx(i,j) + a_x_half(i) * value_dHy_dx
                    value_dHy_dx = value_dHy_dx/K_x_half(i) + memory_dHy_dx(i,j)
                    
                    ! Ez(i,j) = (( caEz(i,j) + caEz(i-1,j) )/2) * Ez(i,j) + &
                    !     (( cbEz(i,j) + cbEz(i-1,j) )/2) * value_dHy_dx 
                    Ez(i,j) = caEz(i,j) * Ez(i,j) + cbEz(i,j) * value_dHy_dx 
                enddo
            enddo

            !----------------------------------------------------------------------------
            Ex(isource,jsource) = Ex(isource,jsource) + &
                            srcx(it) * source%dt / eps11(isource,jsource)
            Ez(isource,jsource) = Ez(isource,jsource) + &
                            srcz(it) * source%dt / eps33(isource,jsource) 
            
            ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
            Ex(1,:) = 0.d0
            Ex(domain%nx,:) = 0.d0
            Ex(:,1) = 0.d0
            Ex(:,domain%nz) = 0.d0

            Ez(1,:) = 0.d0
            Ez(domain%nx,:) = 0.d0
            Ez(:,1) = 0.d0
            Ez(:,domain%nz) = 0.d0

            Hy(1,:) = 0.d0
            Hy(domain%nx,:) = 0.d0
            Hy(:,1) = 0.d0
            Hy(:,domain%nz) = 0.d0

            ! print maximum of norm of velocity
            velocnorm = maxval(sqrt(Ex**2 + Ez**2))
            if (velocnorm > stability_threshold) stop 'code became unstable and blew up'

            call write_image2(Ex, domain%nx, domain%nz, source, it, 'Ex', SINGLE)
            call write_image2(Ez, domain%nx, domain%nz, source, it, 'Ez', SINGLE)
        enddo
        
        
        deallocate(eps11, eps13,  eps33, sig11, sig13,  sig33, srcx, srcz)
        deallocate(K_x, alpha_x, a_x, b_x, K_x_half, alpha_x_half, a_x_half, b_x_half)
        deallocate(K_z, alpha_z, a_z, b_z, K_z_half, alpha_z_half, a_z_half, b_z_half)
        deallocate(epsilonx, epsilonz, sigmax, sigmaz, Ex, Ez, Hy)
        deallocate(memory_dEz_dx, memory_dEx_dz, memory_dHy_dx, memory_dHy_dz)
        
    end subroutine electromag2


    ! =========================================================================
    ! subroutine electromag25(domain, source, emvar, SINGLE_OUTPUT)
    !     !--------------------------------------------------------------------------------------
    !     ! Electromagnetic wave propagation in a 3D grid with Convolutional-PML (C-PML)
    !     ! absorbing conditions for an anisotropic medium.
    !     !
    !     ! This subroutine solves electromagnetic wave propagation using a finite-difference
    !     ! time-domain (FDTD) method in a 3D grid, with PML absorbing conditions.
    !     !
    !     !--------------------------------------------------------------------------------------
        
    !     use constants
        
    !     implicit none

    !     ! Input arguments
    !     type(Domain_Type), intent(in) :: domain
    !     type(Source_Type), intent(in) :: source
    !     type(Electromagnetic3_Variables_Type), intent(inout) :: emvar
    !     logical, intent(in), optional :: SINGLE_OUTPUT
        
    !     ! Local variables
    !     real(real64), dimension(domain%nx,domain%nz) :: epsilonx, epsilony, epsilonz, &
    !                                         sigmax, sigmay, sigmaz

    !     ! real(real64) :: DT
    !     real(real64) :: velocnorm
    !     integer :: isource, jsource, ksource, i, j, k, it

    !     ! Coefficients for the finite difference scheme
    !     real(real64), dimension(domain%nx,domain%nz) :: caEx, cbEx, caEy, cbEy, caEz, cbEz
    !     real(real64) :: daHx, dbHx, daHy, dbHy, daHz, dbHz

    !     real(real64) :: dEx_dy, dEy_dx, dEy_dz, dEz_dy, dEz_dx, dEx_dz, &
    !                     dHx_dy, dHx_dz, dHy_dx, dHy_dz, dHz_dy, dHz_dx


    !     ! Source arrays
    !     real(real64), dimension(source%time_steps) :: srcx, srcy, srcz

    !     ! 1D arrays for the damping profiles in each direction
    !     real(real64), dimension(domain%nx) :: K_x, alpha_x, a_x, b_x, K_x_half, alpha_x_half, a_x_half, b_x_half
    !     real(real64), dimension(domain%ny) :: K_y, alpha_y, a_y, b_y, K_y_half, alpha_y_half, a_y_half, b_y_half
    !     real(real64), dimension(domain%nz) :: K_z, alpha_z, a_z, b_z, K_z_half, alpha_z_half, a_z_half, b_z_half

    !     ! Boolean flag to save as double precision or single precision
    !     logical :: SINGLE

    !     ! Check if SINGLE_OUTPUT is provided, default to single precision if not
    !     if (present(SINGLE_OUTPUT)) then
    !         SINGLE = SINGLE_OUTPUT
    !     else
    !         SINGLE = .TRUE.
    !     endif
        
    !     ! ------------------------ Load Permittivity Coefficients ------------------------
    !     ! Load Epsilon
    !     call material_rw2('eps11.dat', eps11, .TRUE.)
    !     call material_rw2('eps12.dat', eps12, .TRUE.)
    !     call material_rw2('eps13.dat', eps13, .TRUE.)
    !     call material_rw2('eps22.dat', eps22, .TRUE.)
    !     call material_rw2('eps23.dat', eps23, .TRUE.)
    !     call material_rw2('eps33.dat', eps33, .TRUE.)
    !     ! Load Sigma
    !     call material_rw2('sig11.dat', sig11, .TRUE.)
    !     call material_rw2('sig12.dat', sig12, .TRUE.)
    !     call material_rw2('sig13.dat', sig13, .TRUE.)
    !     call material_rw2('sig22.dat', sig22, .TRUE.)
    !     call material_rw2('sig23.dat', sig23, .TRUE.)
    !     call material_rw2('sig33.dat', sig33, .TRUE.)

    !     ! ------------------------ Assign some constants -----------------------
    !     ! Assign the source location indices
    !     isource = source%xind + domain%cpml
    !     jsource = source%yind + domain%cpml
    !     ksource = source%zind + domain%cpml

    !     ! Define the 
    !     ! DT = minval( (/dx, dy, dz/) )/ ( 2.0d0 * Clight/ sqrt( minval( (/ eps11, eps22, eps33 /) ) ) )

    !     ! Compute the coefficients of the FD scheme. First scale the relative 
    !     ! permittivity and permeabilities to get the absolute values 
    !     epsilonx(:,:) = (eps11 + eps12 + eps13)*eps0 
    !     epsilony(:,:) = (eps12 + eps22 + eps23)*eps0
    !     epsilonz(:,:) = (eps13 + eps23 + eps33)*eps0
    !     sigmax(:,:) = sig11 + sig12 + sig13
    !     sigmay(:,:) = sig12 + sig22 + sig23
    !     sigmaz(:,:) = sig13 + sig23 + sig33

    !     caEx(:,:) = ( 1.0d0 - sigmax * source%dt / &
    !                 (2.0d0 * epsilonx ) ) / &
    !                 ( 1.0d0 + sigmax * source%dt / &
    !                 (2.0d0 * epsilonx ) )
    !     cbEx(:,:) = (source%dt / epsilonx ) / &
    !                 ( 1.0d0 + sigmax * source%dt / &
    !                 ( 2.0d0 * epsilonx ) )

    !     caEy(:,:) = ( 1.0d0 - sigmay * source%dt / (2.0d0 * epsilony ) ) / &
    !                 ( 1.0d0 + sigmay * source%dt / (2.0d0 * epsilony ) )
    !     cbEy(:,:) = (source%dt / epsilony ) / &
    !                 ( 1.0d0 + sigmay * source%dt / ( 2.0d0 * epsilony ) )

    !     caEz(:,:) = ( 1.0d0 - sigmaz * source%dt / &
    !                 ( 2.0d0 * epsilonz ) ) / &
    !                 ( 1.0d0 + sigmaz * source%dt / (2.0d0 * epsilonz ) )
    !     cbEz(:,:) = (source%dt / epsilonz ) / &
    !                 ( 1.0d0 + sigmaz * source%dt / (2.0d0 * epsilonz ) )

    !     daHx = source%dt/(4.0d0*mu0*mu)
    !     dbHx = source%dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
    !     daHx = 1.0d0 ! (1-daHz)/(1+daHz) ! 

    !     daHy = source%dt/(4.0d0*mu0*mu)
    !     dbHy = source%dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
    !     daHy = 1.0d0 ! (1-daHz)/(1+daHz) ! 

    !     daHz = source%dt/(4.0d0*mu0*mu)
    !     dbHz = source%dt/mu0 !dt/(mu*mu*dx*(1+daHz) ) 
    !     daHz = 1.0d0 ! (1-daHz)/(1+daHz) ! 


    !     ! ----------------------------------------------------------------------
    !     !---
    !     !--- program starts here
    !     !---

    !     ! ================================ LOAD SOURCE ================================

    !     call loadsource('electromagneticsourcex.dat', source%time_steps, srcx)
    !     call loadsource('electromagneticsourcey.dat', source%time_steps, srcy)
    !     call loadsource('electromagneticsourcez.dat', source%time_steps, srcz)

    !     ! =============================================================================

    !     !--- define profile of absorption in PML region

    !     ! Initialize CPML damping variables
    !     K_x(:) = 1.0d0
    !     K_x_half(:) = 1.0d0
    !     alpha_x(:) = 0.0d0
    !     alpha_x_half(:) = 0.0d0
    !     a_x(:) = 0.0d0
    !     a_x_half(:) = 0.0d0
    !     b_x(:) = 0.0d0 
    !     b_x_half(:) = 0.0d0 

    !     K_y(:) = 1.0d0
    !     K_y_half(:) = 1.0d0
    !     alpha_y(:) = 0.0d0
    !     alpha_y_half(:) = 0.0d0
    !     a_y(:) = 0.0d0
    !     a_y_half(:) = 0.0d0
    !     b_y(:) = 0.d0
    !     K_z(:) = 1.0d0
    !     K_z_half(:) = 1.0d0
    !     alpha_z(:) = 0.0d0
    !     alpha_z_half(:) = 0.0d0
    !     a_z(:) = 0.0d0
    !     a_z_half(:) = 0.0d0

    !     ! ------------------------------ Load the boundary ----------------------------
    !     call loadcpml('kappax_cpml.dat', K_x)
    !     call loadcpml('alphax_cpml.dat', alpha_x)
    !     call loadcpml('acoefx_cpml.dat', a_x)
    !     call loadcpml('bcoefx_cpml.dat', b_x)

    !     call loadcpml('kappay_cpml.dat', K_y)
    !     call loadcpml('alphay_cpml.dat', alpha_y)
    !     call loadcpml('acoefy_cpml.dat', a_y)
    !     call loadcpml('bcoefy_cpml.dat', b_y)

    !     call loadcpml('kappaz_cpml.dat', K_z)
    !     call loadcpml('alphaz_cpml.dat', alpha_z)
    !     call loadcpml('acoefz_cpml.dat', a_z)
    !     call loadcpml('bcoefz_cpml.dat', b_z)

    !     call loadcpml('kappax_half_cpml.dat', K_x_half)
    !     call loadcpml('alphax_half_cpml.dat', alpha_x_half)
    !     call loadcpml('acoefx_half_cpml.dat', a_x_half)
    !     call loadcpml('bcoefx_half_cpml.dat', b_x_half)

    !     call loadcpml('kappay_half_cpml.dat', K_y_half)
    !     call loadcpml('alphay_half_cpml.dat', alpha_y_half)
    !     call loadcpml('acoefy_half_cpml.dat', a_y_half)
    !     call loadcpml('bcoefy_half_cpml.dat', b_y_half)

    !     call loadcpml('kappaz_half_cpml.dat', K_z_half)
    !     call loadcpml('alphaz_half_cpml.dat', alpha_z_half)
    !     call loadcpml('acoefz_half_cpml.dat', a_z_half)
    !     call loadcpml('bcoefz_half_cpml.dat', b_z_half)

    !     ! do i = 1,nz
    !     !   print *, K_z(i), alpha_z(i), a_z(i), b_z(i)
    !     ! enddo

    !     ! -----------------------------------------------------------------------------
    !     ! Load initial conditions
    !     call material_rw3('initialconditionEx.dat', Ex, .TRUE.)
    !     call material_rw3('initialconditionEy.dat', Ey, .TRUE.)
    !     call material_rw3('initialconditionEz.dat', Ez, .TRUE.)
        
    !     call material_rw3('initialconditionHx.dat', Hx, .TRUE.)
    !     call material_rw3('initialconditionHy.dat', Hy, .TRUE.)
    !     call material_rw3('initialconditionHz.dat', Hz, .TRUE.)

    !     ! PML
    !     memory_dEx_dy(:,:,:) = 0.0d0
    !     memory_dEy_dx(:,:,:) = 0.0d0
    !     memory_dEx_dz(:,:,:) = 0.0d0
    !     memory_dEz_dx(:,:,:) = 0.0d0
    !     memory_dEz_dy(:,:,:) = 0.0d0
    !     memory_dEy_dz(:,:,:) = 0.0d0

    !     memory_dHz_dx(:,:,:) = 0.0d0
    !     memory_dHx_dz(:,:,:) = 0.0d0
    !     memory_dHz_dy(:,:,:) = 0.0d0
    !     memory_dHy_dz(:,:,:) = 0.0d0
    !     memory_dHx_dy(:,:,:) = 0.0d0
    !     memory_dHy_dx(:,:,:) = 0.0d0

    !     ! ---
    !     ! ---  beginning of time loop
    !     ! ---
    !     do it = 1,source%time_steps
    !         !--------------------------------------------------------
    !         ! compute magnetic field and update memory variables for C-PML
    !         !--------------------------------------------------------
    !         ! Update Hx
    !         do k = 1,domain%nz-1
    !             do i = 1,domain%nx-1  
    !                 do j = 1,domain%ny-1
    !                     ! Values needed for the magnetic field updates
    !                     dEz_dy = ( Ez(i,j,k) - Ez(i,j+1,k) )/domain%dy
    !                     memory_dEz_dy(i,j,k) = b_y_half(j) * memory_dEz_dy(i,j,k) + a_y_half(j) * dEz_dy
    !                     dEz_dy = dEz_dy/ K_y_half(j) + memory_dEz_dy(i,j,k)

    !                     ! The rest of the equation needed for agnetic field updates
    !                     dEy_dz = ( Ey(i,j,k+1) - Ey(i,j,k) )/domain%dz
    !                     memory_dEy_dz(i,j,k) = b_z_half(k) * memory_dEy_dz(i,j,k) + a_z_half(k) * dEy_dz
    !                     dEy_dz = dEy_dz/ K_z_half(k) + memory_dEy_dz(i,j,k)

    !                     ! Now update the Magnetic field
    !                     Hx(i,j,k) = daHx*Hx(i,j,k) + dbHx*( dEy_dz + dEz_dy )
    !                 enddo
    !             enddo  
    !         enddo

    !             ! Update Hy
    !         do k = 1,domain%nz-1
    !             do i = 1,domain%nx-1      
    !                 do j = 1,domain%ny-1
                    
    !                     ! Values needed for the magnetic field updates
    !                     dEx_dz = ( Ex(i,j,k) - Ex(i,j,k+1) )/domain%dz
    !                     memory_dEx_dz(i,j,k) = b_z(k) * memory_dEx_dz(i,j,k) + &
    !                         a_z(k) * dEx_dz
    !                     dEx_dz = dEx_dz/ K_z(k) + memory_dEx_dz(i,j,k)

    !                     ! The rest of the equation needed for agnetic field updates
    !                     dEz_dx = ( Ez(i+1,j,k) - Ez(i,j,k) )/domain%dx
    !                     memory_dEz_dx(i,j,k) = b_x(i) * memory_dEz_dx(i,j,k) + &
    !                         a_x(i) * dEz_dx
    !                     dEz_dx = dEz_dx/ K_x(i) + memory_dEz_dx(i,j,k)

    !                     ! Now update the Magnetic field
    !                     Hy(i,j,k) = daHy*Hy(i,j,k) + dbHy*( dEz_dx + dEx_dz )

    !                 enddo
    !             enddo  
    !         enddo

    !             ! Update Hz
    !         do k = 2,domain%nz-1
    !             do i = 1,domain%nx-1      
    !                 do j = 1,domain%ny-1
    !                     ! Values needed for the magnetic field updates
    !                     dEx_dy = ( Ex(i,j+1,k) - Ex(i,j,k) )/domain%dy
    !                     memory_dEx_dy(i,j,k) = b_y(j) * memory_dEx_dy(i,j,k) + & 
    !                         a_y(j) * dEx_dy
    !                     dEx_dy = dEx_dy/ K_y(j) + memory_dEx_dy(i,j,k)

    !                     ! The rest of the equation needed for agnetic field updates
    !                     dEy_dx = ( Ey(i,j,k) - Ey(i+1,j,k) )/domain%dx
    !                     memory_dEy_dx(i,j,k) = b_x(i) * memory_dEy_dx(i,j,k) + & 
    !                         a_x(i) * dEy_dx
    !                     dEy_dx = dEy_dx/ K_x(i) + memory_dEy_dx(i,j,k)

    !                     ! Now update the Magnetic field
    !                     Hz(i,j,k) = daHz*Hz(i,j,k) + dbHz*( dEy_dx + dEx_dy )
    !                 enddo
    !             enddo  
    !         enddo

    !         !--------------------------------------------------------
    !         ! compute electric field and update memory variables for C-PML
    !         !--------------------------------------------------------
    !         ! Compute the differences in the x-direction
    !         do k = 2,domain%nz-1
    !             do i = 1,domain%nx-1
    !                 do j = 2,domain%ny-1  
    !                     ! Update the Ex field
    !                     dHz_dy = ( Hz(i,j,k) - Hz(i,j-1,k) )/domain%dy
    !                     memory_dHz_dy(i,j,k) = b_y_half(j) * memory_dHz_dy(i,j,k) + & 
    !                         a_y_half(j) * dHz_dy
    !                     dHz_dy = dHz_dy/K_y_half(j) + memory_dHz_dy(i,j,k)

    !                     ! Changed from half to full node positions 
    !                     dHy_dz = ( Hy(i,j,k-1) - Hy(i,j,k) )/domain%dz
    !                     memory_dHy_dz(i,j,k) = b_z(k) * memory_dHy_dz(i,j,k) + &
    !                         a_z(k) * dHy_dz
    !                     dHy_dz = dHy_dz/K_z(k) + memory_dHy_dz(i,j,k)
                        
    !                     Ex(i,j,k) = caEx(i,k)*Ex(i,j,k) + & 
    !                     cbEx(i,k)*(dHz_dy + dHy_dz) 
    !                 enddo
    !             enddo

    !             ! ! Compute the differences in the y-direction
    !             do i = 2,domain%nx-1 
    !                 do j = 1,domain%ny-1 
    !                     ! Update the Ey field
    !                     dHz_dx = ( Hz(i-1,j,k) - Hz(i,j,k) )/domain%dx ! this is ny-1 length vector
    !                     memory_dHz_dx(i,j,k) = b_x_half(i) * memory_dHz_dx(i,j,k) + & 
    !                         a_x_half(i) * dHz_dx
    !                     dHz_dx = dHz_dx/K_x_half(i) + memory_dHz_dx(i,j,k)

    !                     dHx_dz = ( Hx(i,j,k) - Hx(i,j,k-1) )/domain%dz ! this is ny-1 length vector
    !                     memory_dHx_dz(i,j,k) = b_z_half(k) * memory_dHx_dz(i,j,k) + &
    !                         a_z_half(k) * dHx_dz
    !                     dHx_dz = dHx_dz/K_z_half(k) + memory_dHx_dz(i,j,k)

    !                     ! Ey(i,j,k) = ( ( 4*caEy(i,k) + caEy(i-1,k) + caEy(i,k-1) )/6) * Ey(i,j,k) + & 
    !                     ! ( ( 4*cbEy(i,k) + cbEy(i-1,k) + cbEy(i,k-1) )/6 ) * & 
    !                     ! (dHz_dx + dHx_dz)
    !                     Ey(i,j,k) = caEy(i,k) * Ey(i,j,k) + cbEy(i,k) * (dHz_dx + dHx_dz)
    !                 enddo
    !             enddo
    !         enddo 

    !             ! Compute the differences in the z-direction
    !         do k = 1,domain%nz-1
    !             do i = 2,domain%nx-1  
    !                 do j = 2,domain%ny-1
    !                     ! Update the Ez field
    !                     dHx_dy = ( Hx(i,j-1,k) - Hx(i,j,k) )/domain%dy
    !                     memory_dHx_dy(i,j,k) = b_y_half(j) * memory_dHx_dy(i,j,k) + &
    !                         a_y_half(j) * dHx_dy
    !                     dHx_dy = dHx_dy/K_y_half(j) + memory_dHx_dy(i,j,k)

    !                     dHy_dx = ( Hy(i,j,k) - Hy(i-1,j,k) )/domain%dx
    !                     memory_dHy_dx(i,j,k) = b_x_half(i) * memory_dHy_dx(i,j,k) + &
    !                         a_x_half(i) * dHy_dx
    !                     dHy_dx = dHy_dx/K_x_half(i) + memory_dHy_dx(i,j,k)
                        
    !                     Ez(i,j,k) = ( ( 4*caEz(i,k) + caEz(i-1,k) + caEz(i,k+1) )/6 ) * &
    !                         Ez(i,j,k) + ( ( 4*cbEz(i,k) + cbEz(i-1,k) + cbEz(i,k+1) )/6 ) * & 
    !                     (dHx_dy + dHy_dx)
    !                 enddo
    !             enddo
    !         enddo


    !         ! add the source (force vector located at a given grid point)
    !         Ex(isource,jsource,ksource) = Ex(isource,jsource,ksource) + & 
    !                     srcx(it) * source%dt / eps11(isource,ksource)
    !         Ey(isource,jsource,ksource) = Ey(isource,jsource,ksource) + & 
    !                     srcy(it) * source%dt / eps22(isource,ksource) 
    !         Ez(isource,jsource,ksource) = Ez(isource,jsource,ksource) + & 
    !                     srcz(it) * source%dt / eps33(isource,ksource)
            
    !         ! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
    !         Ex(1,:,:) = 0.0d0
    !         Ex(:,1,:) = 0.0d0
    !         Ex(:,:,1) = 0.0d0
    !         Ex(domain%nx,:,:) = 0.0d0
    !         Ex(:,domain%ny,:) = 0.0d0
    !         Ex(:,:,domain%nz) = 0.0d0 

    !         Ey(1,:,:) = 0.0d0
    !         Ey(:,1,:) = 0.0d0
    !         Ey(:,:,1) = 0.0d0
    !         Ey(domain%nx,:,:) = 0.0d0
    !         Ey(:,domain%ny,:) = 0.0d0
    !         Ey(:,:,domain%nz) = 0.0d0
            
    !         Ez(1,:,:) = 0.0d0
    !         Ez(:,1,:) = 0.0d0
    !         Ez(:,:,1) = 0.0d0
    !         Ez(domain%nx,:,:) = 0.0d0
    !         Ez(:,domain%ny,:) = 0.0d0
    !         Ez(:,:,domain%nz) = 0.0d0
            
    !         Hx(1,:,:) = 0.0d0
    !         Hx(:,1,:) = 0.0d0
    !         Hx(:,:,1) = 0.0d0
    !         Hx(domain%nx,:,:) = 0.0d0
    !         Hx(:,domain%ny,:) = 0.0d0
    !         Hx(:,:,domain%nz) = 0.0d0

    !         Hy(1,:,:) = 0.0d0
    !         Hy(:,1,:) = 0.0d0
    !         Hy(:,:,1) = 0.0d0
    !         Hy(domain%nx,:,:) = 0.0d0
    !         Hy(:,domain%ny,:) = 0.0d0
    !         Hy(:,:,domain%nz) = 0.0d0
            
    !         Hz(1,:,:) = 0.0d0
    !         Hz(:,1,:) = 0.0d0
    !         Hz(:,:,1) = 0.0d0
    !         Hz(domain%nx,:,:) = 0.0d0
    !         Hz(:,domain%ny,:) = 0.0d0
    !         Hz(:,:,domain%nz) = 0.0d0

    !         ! check norm of velocity to make sure the solution isn't diverging
    !         velocnorm = maxval(sqrt(Ex**2.0d0 + Ey**2.0d0 + Ez**2.0d0) )
    !         if (velocnorm > stability_threshold) stop 'code became unstable and blew up'
    !         ! print *,'Max vals for Ex, Ey, Ez: ', maxval(Ex), maxval(Ey), maxval(Ez)

    !         ! print *, maxval(Ex), maxval(Ey), maxval(Ez)
    !         call write_image(Ex, domain, source, it, 'Ex', SINGLE)
    !         call write_image(Ey, domain, source, it, 'Ey', SINGLE)
    !         call write_image(Ez, domain, source, it, 'Ez', SINGLE)

    !     enddo   ! end of time loop
    ! end subroutine electromag25



end module cpmlfdtd