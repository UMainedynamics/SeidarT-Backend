module cpmlfdtd 

    use iso_fortran_env, only: real64
    implicit non 
    
    ! =========================================================================    
    ! Computations are done in double precision and written to binary as single
    ! precision unless specified by the optional logical, SINGLE_OUTPUT.
    subroutine seismic2(domain, time_params, source, SINGLE_OUTPUT)

        implicit none

        ! Input arguments
        type(Domain_Type), intent(in) :: domain
        type(Time_Parameters_Type), intent(in) :: time_params
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        ! Local variables
        
        real(real64), dimension(domain%nx,domain%nz) :: c11, c13, c15, c33, c35, c55, rho
        real(real64), dimension(domain%nx,domain%nz) :: vx, vz, sigmaxx, sigmazz, sigmaxz
        real(real64), dimension(domain%nx,domain%nz) :: &
            memory_dvx_dx, memory_dvx_dz, memory_dvz_dx, memory_dvz_dz, &
            memory_dsigmaxx_dx, memory_dsigmazz_dz, memory_dsigmaxz_dx, memory_dsigmaxz_dz
        real(real64) :: deltarho, velocnorm, value_dvx_dx, value_dvx_dz, &
            value_dvz_dx, value_dvz_dz, value_dsigmaxx_dx, value_dsigmazz_dz, &
            value_dsigmaxz_dx, value_dsigmaxz_dz

        ! 1D arrays for damping profiles
        real(real64), dimension(domain&nx) :: K_x, alpha_x, a_x, b_x, K_x_half, alpha_x_half, a_x_half, b_x_half
        real(real64), dimension(domain&nz) :: K_z, alpha_z, a_z, b_z, K_z_half, alpha_z_half, a_z_half, b_z_half

        real(real64), dimension(domain%nx,domain%nz) :: gamma_x, gamma_z
        real(real64), dimension(time_params%time_step) :: srcx, srcz

        real(real64) :: dt
        integer :: i, j, it, isource, jsource
        logical :: SINGLE

        ! Stability threshold
        real(c_double), parameter :: STABILITY_THRESHOLD = 1.d+25

        ! The default data output is single precision unless SINGLE_OUTPUT is 
        ! set to .FALSE.
        if (present(SINGLE_OUTPUT)) then 
            SINGLE = SINGLE_OUTPUT 
        else
            SINGLE = .TRUE.
        endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Check declarations
        ! -------------------- Load Stiffness Coefficients --------------------
    
        call material_rw('c11.dat', c11, .TRUE.)
        call material_rw('c13.dat', c13, .TRUE.)
        call material_rw('c15.dat', c15, .TRUE.)
        call material_rw('c33.dat', c33, .TRUE.)
        call material_rw('c35.dat', c35, .TRUE.)
        call material_rw('c55.dat', c55, .TRUE.)
        call material_rw('rho.dat', rho, .TRUE.)
        
        ! ------------------- Load Attenuation Coefficients --------------------
        call material_rw('gamma_x.dat', gamma_x, .TRUE.)
        call material_rw('gamma_z.dat', gamma_z, .TRUE.)
        
        ! ------------------------ Assign some constants -----------------------
    
        isource = src(1)+cpml
        jsource = src(2)+cpml
    
        DT = minval( (/dx,dz/) )/ &
            (sqrt( 3.d0*( maxval( (/ c11/rho, c33/rho /) ) ) ) ) 

            ! ================================ LOAD SOURCE ================================
    
        call loadsource('seismicsourcex.dat', time_step, srcx)
        ! We are using the coordinate names x, Z but the math computes the source in 
        ! the x-z plane
        call loadsource('seismicsourcez.dat', time_step, srcz)
    
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


        ! initialize arrays
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

        do it = 1,time_step
            !------------------------------------------------------------
            ! compute stress sigma and update memory variables for C-PML
            !------------------------------------------------------------
            do j = 2,NZ
                do i = 1,NX-1
        
                    value_dvx_dx = (vx(i+1,j) - vx(i,j)) / DX
                    value_dvz_dz = (vz(i,j) - vz(i,j-1)) / DZ
                    value_dvz_dx = (vz(i+1,j) - vz(i,j)) / DX
                    value_dvx_dz = (vx(i,j) - vx(i,j-1)) / DZ

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
                            c15(i,j) * (value_dvz_dx + value_dvx_dz) ) * DT
                    sigmazz(i,j) = sigmazz(i,j) + &
                        (   c13(i,j) * value_dvx_dx + &
                            c33(i,j) * value_dvz_dz + &
                            c35(i,j) * (value_dvz_dx + value_dvx_dz) ) * DT
        
                enddo
            enddo
        
            do j = 1,NZ-1
                do i = 2,NX
        
                value_dvx_dx = (vx(i,j) - vx(i-1,j)) / DX
                value_dvz_dz = (vz(i,j+1) - vz(i,j)) / DZ
                value_dvz_dx = (vz(i,j) - vz(i-1,j)) / DX
                value_dvx_dz = (vx(i,j+1) - vx(i,j)) / DZ
                
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
                        c55(i,j) * (value_dvz_dx + value_dvx_dz) ) * DT
        
                enddo
            enddo
        
            !--------------------------------------------------------
            ! compute velocity and update memory variables for C-PML
            !--------------------------------------------------------
            do j = 2,NZ
                do i = 2,NX
        
                deltarho = ( rho(i,j) + rho(i,j+1) + rho(i+1,j) + rho(i+1,j+1) )/4

                value_dsigmaxx_dx = (sigmaxx(i,j) - sigmaxx(i-1,j)) / DX
                value_dsigmaxz_dz = (sigmaxz(i,j) - sigmaxz(i,j-1)) / DZ
        
                memory_dsigmaxx_dx(i,j) = b_x(i) * memory_dsigmaxx_dx(i,j) + &
                            a_x(i) * value_dsigmaxx_dx
                memory_dsigmaxz_dz(i,j) = b_z(j) * memory_dsigmaxz_dz(i,j) + &
                            a_z(j) * value_dsigmaxz_dz
        
                value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + &
                            memory_dsigmaxx_dx(i,j)
                value_dsigmaxz_dz = value_dsigmaxz_dz / K_z(j) + &
                            memory_dsigmaxz_dz(i,j)
        
                vx(i,j) = vx(i,j)*(1 - gamma_x(i,j) ) + (value_dsigmaxx_dx + value_dsigmaxz_dz) * DT / rho(i,j)
        
                enddo
            enddo
        
            do j = 1,NZ-1
                do i = 1,NX-1
        
                    deltarho = ( 2*rho(i,j) + rho(i+1,j) + rho(i,j+1) )/4
                    value_dsigmaxz_dx = (sigmaxz(i+1,j) - sigmaxz(i,j)) / DX
                    value_dsigmazz_dz = (sigmazz(i,j+1) - sigmazz(i,j)) / DZ
            
                    memory_dsigmaxz_dx(i,j) = b_x_half(i) * memory_dsigmaxz_dx(i,j) + &
                                a_x_half(i) * value_dsigmaxz_dx
                    memory_dsigmazz_dz(i,j) = b_z_half(j) * memory_dsigmazz_dz(i,j) + &
                                a_z_half(j) * value_dsigmazz_dz
            
                    value_dsigmaxz_dx = value_dsigmaxz_dx / K_x_half(i) + memory_dsigmaxz_dx(i,j)
                    value_dsigmazz_dz = value_dsigmazz_dz / K_z_half(j) + memory_dsigmazz_dz(i,j)
            
                    vz(i,j) = vz(i,j)*(1 - gamma_z(i,j) ) + (value_dsigmaxz_dx + value_dsigmazz_dz) * DT / deltarho
        
                enddo
            enddo

            ! Add the source term
            vx(isource,jsource) = vx(isource,jsource) + srcx(it) * DT / rho(isource,jsource)
            vz(isource,jsource) = vz(isource,jsource) + srcz(it) * DT / rho(isource,jsource)
        
            ! Dirichlet conditions (rigid boundaries) on the edges or at the 
            ! bottom of the PML layers
            vx(1,:) = 0.d0
            vx(NX,:) = 0.d0
        
            vx(:,1) = 0.d0
            vx(:,NZ) = 0.d0
        
            vz(1,:) = 0.d0
            vz(NX,:) = 0.d0
        
            vz(:,1) = 0.d0
            vz(:,NZ) = 0.d0
        
            ! print maximum of norm of velocity
            velocnorm = maxval(sqrt(vx**2 + vz**2))
            if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'
        
            call write_image2(vx, nx, nz, src, it, 'Vx', SINGLE)
            call write_image2(vz, nx, nz, src, it, 'Vz', SINGLE)

        enddo   ! end of time loop
    end subroutine seismic2

    ! =========================================================================
    subroutine seismic25(nx, ny, nz, dx, dy, dz, cpml, src, time_step, SINGLE_OUTPUT)
        !--------------------------------------------------------------------------------------
        ! 3D elastic finite-difference code in velocity and stress formulation
        ! with Convolutional-PML (C-PML) absorbing conditions for an anisotropic medium.
        !
        ! This subroutine solves elastic wave propagation using a staggered-grid finite-difference
        ! scheme in a 3D anisotropic medium, including PML (Perfectly Matched Layer) boundaries
        ! for absorbing outgoing waves at the grid boundaries.
        !
        ! INPUT PARAMETERS:
        !   nx, ny, nz (INTEGER, IN)   : Number of grid points in the x, y, and z directions.
        !   dx, dy, dz (REAL, IN)      : Grid spacing in the x, y, and z directions (in meters).
        !   cpml (INTEGER, IN)  : Thickness of the PML absorbing boundary (in grid points).
        !   src (INTEGER, IN)          : Array containing the source location coordinates.
        !   time_step (INTEGER, IN)        : Total number of time steps in the simulation.
        !   SINGLE_OUTPUT (LOGICAL, IN, OPTIONAL) : Flag to specify if the output should be saved 
        !                                           as single precision (default is double precision).
        !
        ! LOCAL VARIABLES:
        !   c11, c12, ..., rho (REAL)  : Elastic stiffness tensors and density matrix for the medium.
        !   vx, vy, vz (REAL)          : Velocity components in the x, y, and z directions.
        !   sigmaxx, sigmayy, sigmazz (REAL) : Stress components in the x, y, and z directions.
        !   gamma_x, gamma_y, gamma_z (REAL) : PML damping coefficients in x, y, and z directions.
        !   srcx, srcy, srcz (REAL)    : Source time history for x, y, and z components.
        !
        !   SINGLE (LOGICAL)           : Boolean flag to determine if output should be in single or 
        !                                double precision. Default is TRUE (single precision).
        !
        ! OUTPUT:
        !   Wave propagation and PML absorbing boundary computations are performed based on the 
        !   given elastic properties and source function.
        !
        !   The results of the simulation can be written to output files as specified by the user.
        !--------------------------------------------------------------------------------------

        use, intrinsic :: iso_c_binding, only: c_int, c_double
        implicit none

        ! Input arguments
        integer(c_int), intent(in) :: nx, ny, nz, cpml, time_step
        real(c_double), intent(in) :: dx, dy, dz
        integer(c_int), dimension(:), intent(in) :: src
        logical, intent(in), optional :: SINGLE_OUTPUT

        ! Local variables
        real(c_double), dimension(nx,nz) :: c11, c12, c13, c14, c15, c16, &
                                            c22, c23, c24, c25, c26, &
                                            c33, c34, c35, c36, &
                                            c44, c45, c46, &
                                            c55, c56, &
                                            c66, &
                                            rho
        real(c_double) :: deltarho, velocnorm, DT
        integer(c_int) :: i, j, k, it, isource, jsource, ksource

        ! Arrays for velocity and stress components
        real(c_double), dimension(nx,ny,nz) :: vx, vy, vz, &
                                            sigmaxx, sigmayy, sigmazz, &
                                            sigmaxy, sigmaxz, sigmayz

        ! Memory arrays for the PML
        real(c_double), dimension(nx,ny,nz) :: memory_dvx_dx, memory_dvx_dy, memory_dvx_dz, &
                                            memory_dvy_dx, memory_dvy_dy, memory_dvy_dz, &
                                            memory_dvz_dx, memory_dvz_dy, memory_dvz_dz, &
                                            memory_dsigmaxx_dx, memory_dsigmayy_dy, memory_dsigmazz_dz, &
                                            memory_dsigmaxy_dx, memory_dsigmaxy_dy, &
                                            memory_dsigmaxz_dx, memory_dsigmaxz_dz, &
                                            memory_dsigmayz_dy, memory_dsigmayz_dz

        ! Values of the velocity and stress differentials
        real(c_double) :: dvx_dx, dvx_dy, dvx_dz, &
                        dvy_dx, dvy_dy, dvy_dz, &
                        dvz_dx, dvz_dy, dvz_dz, &
                        dsigmaxx_dx, dsigmayy_dy, dsigmazz_dz, &
                        dsigmaxy_dx, dsigmaxy_dy, &
                        dsigmaxz_dx, dsigmaxz_dz, &
                        dsigmayz_dy, dsigmayz_dz

        ! 1D arrays for the damping profiles in each direction
        real(c_double), dimension(nx) :: K_x, alpha_x, a_x, b_x, K_x_half, alpha_x_half, a_x_half, b_x_half
        real(c_double), dimension(ny) :: K_y, alpha_y, a_y, b_y, K_y_half, alpha_y_half, a_y_half, b_y_half
        real(c_double), dimension(nz) :: K_z, alpha_z, a_z, b_z, K_z_half, alpha_z_half, a_z_half, b_z_half

        ! Arrays for the PML damping factors
        real(c_double), dimension(nx,nz) :: gamma_x, gamma_y, gamma_z

        ! Source arrays
        real(c_double), dimension(time_step) :: srcx, srcy, srcz

        ! Velocity threshold for stability
        real(c_double), parameter :: STABILITY_THRESHOLD = 1.d+25

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
        call material_rw('gamma_x.dat', gamma_x, .TRUE.)
        call material_rw('gamma_z.dat', gamma_z, .TRUE.)
        call material_rw('gamma_y.dat', gamma_y, .TRUE.)
        
        ! ------------------------ Assign some constants -----------------------
        isource = src(1)+cpml
        jsource = src(2)+cpml
        ksource = src(3)+cpml

        ! To ensure a courant number <= 1.0, we can calculate the time step from
        ! the velocity
        DT = 0.7 * minval( (/dx,dy,dz/) )/ &
        ( sqrt( 3.d0 * ( maxval( (/ c11/rho, c22/rho, c33/rho /) ) ) ) )

        ! ================================ LOAD SOURCE ================================

        call loadsource('seismicsourcex.dat', time_step, srcx)
        call loadsource('seismicsourcey.dat', time_step, srcy)
        call loadsource('seismicsourcez.dat', time_step, srcz)

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
        ! initialize arrays
        vx(:,:,:) = 0.d0
        vy(:,:,:) = 0.d0
        vz(:,:,:) = 0.d0

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
        do it = 1,time_step
            !------------------------------------------------------------
            ! compute stress sigma and update memory variables for C-PML
            !------------------------------------------------------------
            ! Update in the x direction
            do k = 2,nz
                do j = 2,NY
                    do i = 1,NX-1

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

                        sigmaxx(i,j,k) = sigmaxx(i,j,k) + &
                        (   c11(i,k) * dvx_dx + c12(i,k) * dvy_dy + c13(i,k) * dvz_dz + &
                            c14(i,k) * (dvy_dz + dvz_dy) + c15(i,k) * (dvx_dz + dvz_dx) + &
                            c16(i,k) * (dvx_dy + dvz_dy) ) * dt

                        ! Full 3D will need a gradient in the y-direction
                        sigmayy(i,j,k) = sigmayy(i,j,k) + &
                        (   c12(i,k) * dvx_dx + c22(i,k) * dvy_dy + c23(i,k) * dvz_dz + &
                            c24(i,k) * (dvy_dz + dvz_dy) + c25(i,k) * (dvx_dz + dvz_dx) + &
                            c26(i,k) * (dvy_dx + dvx_dy) ) * dt

                        sigmazz(i,j,k) = sigmazz(i,j,k) + &
                        (   c13(i,k) * dvx_dx + c23(i,k) * dvy_dy + c33(i,k) * dvz_dz + &
                            c34(i,k) * (dvy_dz + dvz_dy) + c35(i,k) * (dvx_dz + dvz_dx) + &
                            c36(i,k) * (dvy_dx + dvx_dy) ) * dt

                    enddo
                enddo
            enddo

            ! Update sigmaxy, x-direction is full nodes
            do k = 2,nz
                do j = 1,ny-1
                    do i = 2,nx

                        dvx_dx = (vx(i,j,k) - vx(i-1,j,k)) / DX
                        dvy_dx = (vy(i,j,k) - vy(i-1,j,k)) / DX
                        dvz_dx = (vz(i,j,k) - vz(i-1,j,k)) / DX

                        dvy_dy = (vy(i,j+1,k) - vy(i,j,k)) / DY
                        dvx_dy = (vx(i,j+1,k) - vx(i,j,k)) / DY
                        dvz_dy = (vz(i,j+1,k) - vz(i,j,k)) / DY

                        dvz_dz = (vz(i,j,k) - vz(i,j,k-1)) / DZ
                        dvx_dz = (vx(i,j,k) - vx(i,j,k-1)) / DZ
                        dvy_dz = (vy(i,j,k) - vy(i,j,k-1)) / DZ

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
                        dvy_dy = dvy_dy / K_y_half(j) + memory_dvy_dy(i,j,k)
                        dvx_dy = dvx_dy / K_y_half(j) + memory_dvx_dy(i,j,k)
                        dvz_dy = dvz_dy / K_y_half(j) + memory_dvz_dy(i,j,k)
                        dvz_dz = dvz_dz / K_z_half(k) + memory_dvz_dz(i,j,k)
                        dvy_dz = dvy_dz / K_z_half(k) + memory_dvy_dz(i,j,k)

                        sigmaxy(i,j,k) = sigmaxy(i,j,k) + &
                        (   c16(i,k) * dvx_dx + c26(i,k) * dvy_dy + c36(i,k) * dvz_dz + &
                            c46(i,k) * (dvz_dy + dvy_dz) + c56(i,k) * (dvz_dx + dvx_dz) + &
                            c66(i,k) * (dvy_dx + dvx_dy) ) * dt

                    enddo
                enddo
            enddo

            ! Update sigmaxz, z-direction is full nodes
            do k = 1,nz-1
                do j = 2,ny
                    do i = 2,nx

                        dvx_dx = (vx(i,j,k) - vx(i-1,j,k)) / DX
                        dvy_dx = (vy(i,j,k) - vy(i-1,j,k)) / DX
                        dvz_dx = (vz(i,j,k) - vz(i-1,j,k)) / DX
                        dvy_dy = (vy(i,j,k) - vy(i,j-1,k)) / DY
                        dvz_dy = (vz(i,j,k) - vz(i,j-1,k)) / DY
                        dvx_dy = (vx(i,j,k) - vx(i,j-1,k)) / DY
                        dvz_dz = (vz(i,j,k+1) - vz(i,j,k)) / DZ
                        dvx_dz = (vx(i,j,k+1) - vx(i,j,k)) /DZ
                        dvy_dz = (vy(i,j,k+1) - vy(i,j,k)) / DZ

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

                        sigmaxz(i,j,k) = sigmaxz(i,j,k) + &
                            (   c15(i,k) * dvx_dx + c25(i,k) * dvy_dy + c35(i,k) * dvz_dz + &
                                c45(i,k) * ( dvx_dz + dvz_dx) + c55(i,k) * ( dvx_dz + dvz_dx) + &
                                c56(i,k) * ( dvx_dy + dvy_dx) ) * dt 

                    enddo
                enddo

                !   ! update sigmayz, y-direction is full nodes
                do j = 1,ny-1
                    do i = 1,nx-1

                        dvx_dx = (vx(i+1,j,k) - vx(i,j,k)) / DX
                        dvy_dx = (vy(i+1,j,k) - vy(i,j,k)) / DX
                        dvz_dx = (vz(i+1,j,k) - vz(i,j,k)) / DX
                        dvy_dy = (vy(i,j+1,k) - vy(i,j,k)) / DY
                        dvx_dy = (vx(i,j+1,k) - vx(i,j,k)) / DY
                        dvz_dy = (vz(i,j+1,k) - vz(i,j,k)) / DY
                        dvz_dz = (vz(i,j,k+1) - vz(i,j,k)) / DZ
                        dvx_dz = (vx(i,j,k+1) - vx(i,j,k)) / DZ 
                        dvy_dz = (vy(i,j,k+1) - vy(i,j,k)) / DZ

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

                        sigmayz(i,j,k) = sigmayz(i,j,k)  + &
                            (   c14(i,k) * dvx_dx + c24(i,k) * dvy_dy + c34(i,k) * dvz_dz + &
                                c44(i,k) * ( dvy_dz + dvz_dy) + c45(i,k) * ( dvx_dz + dvz_dx) + &
                                c46(i,k) * ( dvy_dx + dvx_dy) ) * dt 
                    enddo
                enddo
            enddo

            !--------------------------------------------------------
            ! compute velocity and update memory variables for C-PML
            !--------------------------------------------------------
            do k = 2,nz
                do j = 2,NY
                    do i = 2,NX
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

                        vx(i,j,k) = vx(i,j,k) * (1 - gamma_x(i,j) ) + &
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

                        vy(i,j,k) = vy(i,j,k) * (1 - gamma_y(i,j) )+ &
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

                        vz(i,j,k) = vz(i,j,k) * (1 - gamma_z(i,j) )+ &
                            (dsigmaxz_dx + dsigmayz_dy + dsigmazz_dz) * &
                            dt / deltarho !rho(i,k)

                    enddo
                enddo
            enddo

            vx(isource,jsource,ksource) = vx(isource,jsource,ksource) + srcx(it) * DT / rho(isource,ksource)
            vy(isource,jsource,ksource) = vy(isource,jsource,ksource) + srcy(it) * DT / rho(isource,ksource)
            vz(isource,jsource,ksource) = vz(isource,jsource,ksource) + srcz(it) * DT / rho(isource,ksource)

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

            vx(NX,:,:) = 0.d0
            vy(NX,:,:) = 0.d0
            vz(NX,:,:) = 0.d0

            vx(:,NY,:) = 0.d0
            vy(:,NY,:) = 0.d0
            vz(:,NY,:) = 0.d0

            vx(:,:,NZ) = 0.d0
            vy(:,:,NZ) = 0.d0
            vz(:,:,NZ) = 0.d0

            ! check norm of velocity to make sure the solution isn't diverging
            velocnorm = maxval( sqrt(vx**2 + vy**2 + vz**2) )
            ! print *,'Time step # ',it,' out of ',time_step
            ! print *,'Time: ',(it-1)*DT,' seconds'
            ! print *,'Max vals for vx, vy, vz: ', maxval(vx), maxval(vy), maxval(vz)

            if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

            ! Write the velocity values to an unformatted binary file
            call write_image3(vx, nx, ny, nz, src, it, 'Vx', SINGLE)
            call write_image3(vy, nx, ny, nz, src, it, 'Vy', SINGLE)
            call write_image3(vz, nx, ny, nz, src, it, 'Vz', SINGLE)
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
    subroutine electromag2(nx, nz, dx, dz, cpml, src, time_step, SINGLE_OUTPUT) bind(C, name="electromag2_")
        !--------------------------------------------------------------------------------------
        ! Electromagnetic wave propagation in a 2D grid with Convolutional-PML (C-PML)
        ! absorbing conditions for an anisotropic medium.
        !
        ! This subroutine solves electromagnetic wave propagation using a finite-difference
        ! time-domain (FDTD) method in a 2D grid, with PML absorbing conditions.
        !
        ! INPUT PARAMETERS:
        !   nx, nz (INTEGER, IN)          : Number of grid points in the x and z directions.
        !   dx, dz (REAL, IN)             : Grid spacing in the x and z directions (in meters).
        !   cpml (INTEGER, IN)     : Thickness of the PML absorbing boundary (in grid points).
        !   src (INTEGER ARRAY, IN)       : Source location array.
        !   time_step (INTEGER, IN)           : Total number of time steps in the simulation.
        !   SINGLE_OUTPUT (LOGICAL, IN, OPTIONAL) : Flag to specify if the output should be saved
        !                                           as single precision (default is double precision).
        !
        ! LOCAL VARIABLES:
        !   eps11, sig11, ... (REAL)      : Permittivity and conductivity arrays in the grid.
        !   Ex, Ez, Hy (REAL)             : Electric and magnetic field components.
        !   caEx, cbEx, caEz, cbEz (REAL) : Coefficients for the finite difference scheme.
        !   memory_dEz_dx, ... (REAL)     : Arrays to store memory variables for the PML.
        !   velocnorm (REAL)              : Normalized velocity for stability checks.
        !
        ! OUTPUT:
        !   The results of the simulation are calculated for electric and magnetic fields,
        !   and optionally written to output files or further processed.
        !--------------------------------------------------------------------------------------

        use, intrinsic :: iso_c_binding, only: c_int, c_double
        implicit none

        ! Input arguments
        integer(c_int), intent(in) :: nx, nz, cpml, time_step
        real(c_double), intent(in) :: dx, dz
        integer(c_int), dimension(:), intent(in) :: src
        logical, intent(in), optional :: SINGLE_OUTPUT

        ! Local variables
        real(c_double), dimension(nx,nz) :: eps11, eps13, eps33, &
                                            sig11, sig13, sig33, &
                                            epsilonx, epsilonz, &
                                            sigmax, sigmaz

        real(c_double) :: DT
        real(c_double), dimension(time_step) :: srcx, srcz
        integer(c_int) :: isource, jsource, i, j, it

        ! Constants
        real(c_double), parameter :: PI = 3.141592653589793238462643d0
        real(c_double), parameter :: Clight = 2.9979458d+8
        real(c_double), parameter :: mu0 = 4.0d0 * PI * 1.0d-7
        real(c_double), parameter :: eps0 = 8.85418782d-12
        real(c_double), parameter :: mu = 1.d0
        real(c_double), parameter :: STABILITY_THRESHOLD = 1.d+25

        ! Main arrays for electric and magnetic field components
        real(c_double), dimension(nx,nz) :: Ex, Ez, Hy

        ! Coefficients for the finite difference scheme
        real(c_double), dimension(nx,nz) :: caEx, cbEx
        real(c_double), dimension(nx,nz) :: caEz, cbEz
        real(c_double) :: daHy, dbHy
        real(c_double) :: value_dEx_dz, value_dEz_dx, value_dHy_dz, value_dHy_dx

        ! Arrays for the memory variables in PML
        real(c_double), dimension(nx,nz) :: memory_dEz_dx, memory_dEx_dz
        real(c_double), dimension(nx,nz) :: memory_dHy_dx, memory_dHy_dz

        ! 1D arrays for the damping profiles
        real(c_double), dimension(nx) :: K_x, alpha_x, a_x, b_x, &
                                        K_x_half, alpha_x_half, a_x_half, b_x_half
        real(c_double), dimension(nz) :: K_z, alpha_z, a_z, b_z, &
                                        K_z_half, alpha_z_half, a_z_half, b_z_half

        ! Velocity normalization factor
        real(c_double) :: velocnorm

        ! Boolean flag to save as double precision or single precision
        logical :: SINGLE

        ! Check if SINGLE_OUTPUT is provided, default to single precision if not
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        endif
        ! =============================================================================
        ! ----------------------- Load Permittivity Coefficients ----------------------
        call material_rw('eps11.dat', eps11, .TRUE.)
        call material_rw('eps13.dat', eps13, .TRUE.)
        call material_rw('eps33.dat', eps33, .TRUE.) ! We will change y to z soon
        call material_rw('sig11.dat', sig11, .TRUE.)
        call material_rw('sig13.dat', sig13, .TRUE.)
        call material_rw('sig33.dat', sig33, .TRUE.)

        ! ------------------------ Assign some constants -----------------------
        ! Assign the source location indices
        isource = int(src(1)) + cpml
        jsource = int(src(2)) + cpml

        ! Define the 
        DT = minval( (/dx, dz/) )/ ( 2.d0 * Clight/sqrt( minval( (/ eps11, eps33 /) ) ) ) 

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

        ! ================================ LOAD SOURCE ================================
        call loadsource('electromagneticsourcex.dat', time_step, srcx)
        call loadsource('electromagneticsourcez.dat', time_step, srcz)

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

        ! initialize arrays
        Ex(:,:) = 0.0d0
        Ez(:,:) = 0.0d0
        Hy(:,:) = 0.0d0

        ! PML
        memory_dEx_dz(:,:) = 0.0d0
        memory_dEz_dx(:,:) = 0.0d0
        memory_dHy_dx(:,:) = 0.0d0
        memory_dHy_dz(:,:) = 0.0d0

        !---
        !---  beginning of time loop
        !---
        do it = 1,time_step
            !--------------------------------------------------------
            ! compute magnetic field and update memory variables for C-PML
            !--------------------------------------------------------
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

                enddo  
            enddo

            !--------------------------------------------------------
            ! compute electric field and update memory variables for C-PML
            !--------------------------------------------------------
            ! Compute the differences in the y-direction
            do j = 2,nz
                do i = 1,nx
                    ! Update the Ex field
                    value_dHy_dz = ( Hy(i,j) - Hy(i,j-1) )/dz ! this is nz-1 length vector
                    memory_dHy_dz(i,j) = b_z(j) * memory_dHy_dz(i,j) + a_z(j) * value_dHy_dz
                    value_dHy_dz = value_dHy_dz/K_z(j) + memory_dHy_dz(i,j)

                    ! Ex(i,j) = (( caEx(i,j) + caEx(i,j-1) )/2) * Ex(i,j) + &
                    !     (( cbEx(i,j) + cbEx(i,j-1) )/2 ) * value_dHy_dz
                    Ex(i,j) = caEx(i,j) * Ex(i,j) + cbEx(i,j) * value_dHy_dz
                enddo
            enddo

            do j = 1,nz
                do i = 2,nx
                    ! Update the Ez field
                    value_dHy_dx = ( Hy(i,j) - Hy(i-1,j) )/dx
                    memory_dHy_dx(i,j) = b_x_half(i) * memory_dHy_dx(i,j) + a_x_half(i) * value_dHy_dx
                    value_dHy_dx = value_dHy_dx/K_x_half(i) + memory_dHy_dx(i,j)
                    
                    ! Ez(i,j) = (( caEz(i,j) + caEz(i-1,j) )/2) * Ez(i,j) + &
                    !     (( cbEz(i,j) + cbEz(i-1,j) )/2) * value_dHy_dx 
                    Ez(i,j) = caEz(i,j) * Ez(i,j) + cbEz(i,j) * value_dHy_dx 
                enddo
            enddo

            !----------------------------------------------------------------------------
            Ex(isource,jsource) = Ex(isource,jsource) + srcx(it) * DT / eps11(isource,jsource)
            Ez(isource,jsource) = Ez(isource,jsource) + srcz(it) * DT / eps33(isource,jsource) 
            
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
            velocnorm = maxval(sqrt(Ex**2 + Ez**2))
            if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

            call write_image2(Ex, nx, nz, src, it, 'Ex', SINGLE)
            call write_image2(Ez, nx, nz, src, it, 'Ez', SINGLE)
        enddo
    end subroutine electromag2


    ! =========================================================================
    subroutine electromag25(nx, ny, nz, dx, dy, dz, cpml, src, time_step, SINGLE_OUTPUT) bind(C, name="electromag25_")
        !--------------------------------------------------------------------------------------
        ! Electromagnetic wave propagation in a 3D grid with Convolutional-PML (C-PML)
        ! absorbing conditions for an anisotropic medium.
        !
        ! This subroutine solves electromagnetic wave propagation using a finite-difference
        ! time-domain (FDTD) method in a 3D grid, with PML absorbing conditions.
        !
        ! INPUT PARAMETERS:
        !   nx, ny, nz (INTEGER, IN)      : Number of grid points in x, y, and z directions.
        !   dx, dy, dz (REAL, IN)         : Grid spacing in the x, y, and z directions (in meters).
        !   cpml (INTEGER, IN)     : Thickness of the PML absorbing boundary (in grid points).
        !   src (INTEGER ARRAY, IN)       : Source location array.
        !   time_step (INTEGER, IN)           : Total number of time steps in the simulation.
        !   SINGLE_OUTPUT (LOGICAL, IN, OPTIONAL) : Flag to specify if the output should be saved
        !                                           as single precision (default is double precision).
        !
        ! LOCAL VARIABLES:
        !   eps11, sig11, ... (REAL)      : Permittivity and conductivity arrays in the grid.
        !   Ex, Ey, Ez, Hx, Hy, Hz (REAL) : Electric and magnetic field components.
        !   memory_dEz_dx, ... (REAL)     : Arrays to store memory variables for the PML.
        !   velocnorm (REAL)              : Normalized velocity for stability checks.
        !
        ! OUTPUT:
        !   The results of the simulation are calculated for electric and magnetic fields,
        !   and optionally written to output files or further processed.
        !--------------------------------------------------------------------------------------

        use, intrinsic :: iso_c_binding, only: c_int, c_double
        implicit none

        ! Input arguments
        integer(c_int), intent(in) :: nx, ny, nz, cpml, time_step
        real(c_double), intent(in) :: dx, dy, dz
        integer(c_int), dimension(:), intent(in) :: src
        logical, intent(in), optional :: SINGLE_OUTPUT

        ! Local variables
        real(c_double), dimension(nx,nz) :: eps11, eps22, eps33, &
                                            eps12, eps13, eps23, &
                                            sig11, sig22, sig33, &
                                            sig12, sig13, sig23, &
                                            epsilonx, epsilony, epsilonz, &
                                            sigmax, sigmay, sigmaz

        real(c_double) :: DT
        real(c_double) :: velocnorm
        integer(c_int) :: isource, jsource, ksource, i, j, k, it

        ! Constants
        real(c_double), parameter :: PI = 3.141592653589793238462643d0
        real(c_double), parameter :: Clight = 2.9979458d+8
        real(c_double), parameter :: mu0 = 4.0d0 * PI * 1.0d-7
        real(c_double), parameter :: eps0 = 8.85418782d-12
        real(c_double), parameter :: mu = 1.0d0
        real(c_double), parameter :: STABILITY_THRESHOLD = 1.d+25

        ! Main arrays for electric and magnetic field components
        real(c_double), dimension(nx,ny,nz) :: Ex, Ey, Ez
        real(c_double), dimension(nx,ny,nz) :: Hx, Hy, Hz

        ! Coefficients for the finite difference scheme
        real(c_double), dimension(nx,nz) :: caEx, cbEx, caEy, cbEy, caEz, cbEz
        real(c_double) :: daHx, dbHx, daHy, dbHy, daHz, dbHz

        real(c_double) :: dEx_dy, dEy_dx, dEy_dz, dEz_dy, dEz_dx, dEx_dz, &
                        dHx_dy, dHx_dz, dHy_dx, dHy_dz, dHz_dy, dHz_dx

        ! Arrays for the memory variables in PML
        real(c_double), dimension(nx,ny,nz) :: memory_dEy_dx, memory_dEx_dy, &
                                            memory_dEz_dx, memory_dEx_dz, &
                                            memory_dEy_dz, memory_dEz_dy

        real(c_double), dimension(nx,ny,nz) :: memory_dHz_dx, memory_dHx_dz, &
                                            memory_dHy_dx, memory_dHx_dy, &
                                            memory_dHy_dz, memory_dHz_dy

        ! Source arrays
        real(c_double), dimension(time_step) :: srcx, srcy, srcz

        ! 1D arrays for the damping profiles in each direction
        real(c_double), dimension(nx) :: K_x, alpha_x, a_x, b_x, K_x_half, alpha_x_half, a_x_half, b_x_half
        real(c_double), dimension(ny) :: K_y, alpha_y, a_y, b_y, K_y_half, alpha_y_half, a_y_half, b_y_half
        real(c_double), dimension(nz) :: K_z, alpha_z, a_z, b_z, K_z_half, alpha_z_half, a_z_half, b_z_half

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
        call material_rw('eps11.dat', eps11, .TRUE.)
        call material_rw('eps12.dat', eps12, .TRUE.)
        call material_rw('eps13.dat', eps13, .TRUE.)
        call material_rw('eps22.dat', eps22, .TRUE.)
        call material_rw('eps23.dat', eps23, .TRUE.)
        call material_rw('eps33.dat', eps33, .TRUE.)
        ! Load Sigma
        call material_rw('sig11.dat', sig11, .TRUE.)
        call material_rw('sig12.dat', sig12, .TRUE.)
        call material_rw('sig13.dat', sig13, .TRUE.)
        call material_rw('sig22.dat', sig22, .TRUE.)
        call material_rw('sig23.dat', sig23, .TRUE.)
        call material_rw('sig33.dat', sig33, .TRUE.)

        ! ------------------------ Assign some constants -----------------------
        ! Assign the source location indices
        isource = int(src(1)) + cpml
        jsource = int(src(2)) + cpml
        ksource = int(src(3)) + cpml

        ! Define the 
        DT = minval( (/dx, dy, dz/) )/ ( 2.0d0 * Clight/ sqrt( minval( (/ eps11, eps22, eps33 /) ) ) )

        ! Compute the coefficients of the FD scheme. First scale the relative 
        ! permittivity and permeabilities to get the absolute values 
        epsilonx(:,:) = (eps11 + eps12 + eps13)*eps0 
        epsilony(:,:) = (eps12 + eps22 + eps23)*eps0
        epsilonz(:,:) = (eps13 + eps23 + eps33)*eps0
        sigmaX(:,:) = sig11 + sig12 + sig13
        sigmay(:,:) = sig12 + sig22 + sig23
        sigmaz(:,:) = sig13 + sig23 + sig33

        caEx(:,:) = ( 1.0d0 - sigmaX * dt / ( 2.0d0 * epsilonx ) ) / &
                    ( 1.0d0 + sigmaX * dt / (2.0d0 * epsilonx ) )
        cbEx(:,:) = (dt / epsilonx ) / ( 1.0d0 + sigmax * dt / ( 2.0d0 * epsilonx ) )

        caEy(:,:) = ( 1.0d0 - sigmay * dt / ( 2.0d0 * epsilony ) ) / &
                    ( 1.0d0 + sigmay * dt / (2.0d0 * epsilony ) )
        cbEy(:,:) = (dt / epsilony ) / ( 1.0d0 + sigmay * dt / ( 2.0d0 * epsilony ) )

        caEz(:,:) = ( 1.0d0 - sigmaz * dt / ( 2.0d0 * epsilonz ) ) / &
                    ( 1.0d0 + sigmaz * dt / (2.0d0 * epsilonz ) )
        cbEz(:,:) = (dt / epsilonz ) / ( 1.0d0 + sigmaz * dt / ( 2.0d0 * epsilonz ) )

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

        call loadsource('electromagneticsourcex.dat', time_step, srcx)
        call loadsource('electromagneticsourcey.dat', time_step, srcy)
        call loadsource('electromagneticsourcez.dat', time_step, srcz)

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

        ! initialize arrays
        Ex(:,:,:) = 0.0d0
        Ey(:,:,:) = 0.0d0
        Ez(:,:,:) = 0.0d0

        Hx(:,:,:) = 0.0d0
        Hy(:,:,:) = 0.0d0
        Hz(:,:,:) = 0.0d0


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
        do it = 1,time_step
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
                        
                        Ez(i,j,k) = ( ( 4*caEz(i,k) + caEz(i-1,k) + caEz(i,k+1) )/6 ) * Ez(i,j,k) + & 
                        ( ( 4*cbEz(i,k) + cbEz(i-1,k) + cbEz(i,k+1) )/6 ) * & 
                        (dHx_dy + dHy_dx)
                    enddo
                enddo
            enddo


            ! add the source (force vector located at a given grid point)
            Ex(isource,jsource,ksource) = Ex(isource,jsource,ksource) + srcx(it) * DT / eps11(isource,ksource)
            Ey(isource,jsource,ksource) = Ey(isource,jsource,ksource) + srcy(it) * DT / eps22(isource,ksource) 
            Ez(isource,jsource,ksource) = Ez(isource,jsource,ksource) + srcz(it) * DT / eps33(isource,ksource)
            
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
            if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'
            ! print *,'Max vals for Ex, Ey, Ez: ', maxval(Ex), maxval(Ey), maxval(Ez)

            ! print *, maxval(Ex), maxval(Ey), maxval(Ez)
            call write_image3(Ex, nx, ny, nz, src, it, 'Ex', SINGLE)
            call write_image3(Ey, nx, ny, nz, src, it, 'Ey', SINGLE)
            call write_image3(Ez, nx, ny, nz, src, it, 'Ez', SINGLE)

        enddo   ! end of time loop
    end subroutine electromag25



end module cpmlfdtd