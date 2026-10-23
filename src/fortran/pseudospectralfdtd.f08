module pseudospectralfdtd
    
    use iso_c_binding
    use iso_fortran_env, only: real64
    use spectral_operators
    use absorbing_boundary
    use source_module
    
    use seidartio
    use seidart_types
    use constants 
    use omp_lib
    
    use pseudospectral_stencils
    
    implicit none
    include 'fftw3.f03'
    
    private
    
    public :: load_elastic_coefficients2, load_elastic_coefficients3, &
                seismic2_pseudospectral_rk4, seismic25_pseudospectral_rk4, &
                seismic3_pseudospectral_rk4
                ! electromag2_pseudospectral_rk4, electromag25_pseudospectral_rk4
    
    
    contains 
    
    ! =========================================================================
    subroutine load_elastic_coefficients2(nx, nz, C, gamma, rho)
        integer, intent(in) :: nx, nz
        real(real64), intent(out) :: C(6, nx, nz)
        real(real64), intent(out) :: gamma(3, nx, nz)
        real(real64), intent(out) :: rho(nx, nz)
        
        ! Load stiffness coefficients (Upper-triangular Voigt packing)
        call material_rw2('c11.dat', C(1,:,:), .TRUE.)
        call material_rw2('c13.dat', C(2,:,:), .TRUE.)
        call material_rw2('c15.dat', C(3,:,:), .TRUE.)
        call material_rw2('c33.dat', C(4,:,:), .TRUE.)
        call material_rw2('c35.dat', C(5,:,:), .TRUE.)
        call material_rw2('c55.dat', C(6,:,:), .TRUE.)
        
        ! Load viscoelastic coefficients
        call material_rw2('gamma_x.dat',  gamma(1,:,:), .TRUE.)
        call material_rw2('gamma_z.dat',  gamma(2,:,:), .TRUE.)
        call material_rw2('gamma_xz.dat', gamma(3,:,:), .TRUE.)
        
        call material_rw2('density.dat', rho, .TRUE.)        
    end subroutine load_elastic_coefficients2
    
    ! ------------------------------------------------------------------------
    subroutine load_elastic_coefficients25(nx, nz, C, gamma, rho)
        integer, intent(in) :: nx, nz
        real(real64), intent(out) :: C(21, nx, nz)
        real(real64), intent(out) :: gamma(6, nx, nz)
        real(real64), intent(out) :: rho(nx, nz)
        
        ! Load stiffness coefficients (Upper-triangular Voigt packing)
        call material_rw2('c11.dat', C(1,:,:), .TRUE.)
        call material_rw2('c12.dat', C(2,:,:), .TRUE.)
        call material_rw2('c13.dat', C(3,:,:), .TRUE.)
        call material_rw2('c14.dat', C(4,:,:), .TRUE.)
        call material_rw2('c15.dat', C(5,:,:), .TRUE.)
        call material_rw2('c16.dat', C(6,:,:), .TRUE.)
        call material_rw2('c22.dat', C(7,:,:), .TRUE.)
        call material_rw2('c23.dat', C(8,:,:), .TRUE.)
        call material_rw2('c24.dat', C(9,:,:), .TRUE.)
        call material_rw2('c25.dat', C(10,:,:), .TRUE.)
        call material_rw2('c26.dat', C(11,:,:), .TRUE.)
        call material_rw2('c33.dat', C(12,:,:), .TRUE.)
        call material_rw2('c34.dat', C(13,:,:), .TRUE.)
        call material_rw2('c35.dat', C(14,:,:), .TRUE.)
        call material_rw2('c36.dat', C(15,:,:), .TRUE.)
        call material_rw2('c44.dat', C(16,:,:), .TRUE.)
        call material_rw2('c45.dat', C(17,:,:), .TRUE.)
        call material_rw2('c46.dat', C(18,:,:), .TRUE.)
        call material_rw2('c55.dat', C(19,:,:), .TRUE.)
        call material_rw2('c56.dat', C(20,:,:), .TRUE.)
        call material_rw2('c66.dat', C(21,:,:), .TRUE.)
        
        ! Load viscoelastic coefficients
        call material_rw2('gamma_x.dat',  gamma(1,:,:), .TRUE.)
        call material_rw2('gamma_y.dat',  gamma(2,:,:), .TRUE.)
        call material_rw2('gamma_z.dat',  gamma(3,:,:), .TRUE.)
        call material_rw2('gamma_yz.dat', gamma(4,:,:), .TRUE.)
        call material_rw2('gamma_xz.dat', gamma(5,:,:), .TRUE.)
        call material_rw2('gamma_xy.dat', gamma(6,:,:), .TRUE.)
        
        call material_rw2('density.dat', rho, .TRUE.)        
    end subroutine load_elastic_coefficients25
    
    ! ------------------------------------------------------------------------
    subroutine load_elastic_coefficients3(nx, ny, nz, C, gamma, rho)
        integer, intent(in) :: nx, ny, nz
        real(real64), intent(out) :: C(21, nx, ny, nz)
        real(real64), intent(out) :: gamma(6, nx, ny, nz)
        real(real64), intent(out) :: rho(nx, ny, nz)
        
        ! Load stiffness coefficients (Upper-triangular Voigt packing)
        call material_rw3('c11.dat', C(1,:,:,:), .TRUE.)
        call material_rw3('c12.dat', C(2,:,:,:), .TRUE.)
        call material_rw3('c13.dat', C(3,:,:,:), .TRUE.)
        call material_rw3('c14.dat', C(4,:,:,:), .TRUE.)
        call material_rw3('c15.dat', C(5,:,:,:), .TRUE.)
        call material_rw3('c16.dat', C(6,:,:,:), .TRUE.)
        call material_rw3('c22.dat', C(7,:,:,:), .TRUE.)
        call material_rw3('c23.dat', C(8,:,:,:), .TRUE.)
        call material_rw3('c24.dat', C(9,:,:,:), .TRUE.)
        call material_rw3('c25.dat', C(10,:,:,:), .TRUE.)
        call material_rw3('c26.dat', C(11,:,:,:), .TRUE.)
        call material_rw3('c33.dat', C(12,:,:,:), .TRUE.)
        call material_rw3('c34.dat', C(13,:,:,:), .TRUE.)
        call material_rw3('c35.dat', C(14,:,:,:), .TRUE.)
        call material_rw3('c36.dat', C(15,:,:,:), .TRUE.)
        call material_rw3('c44.dat', C(16,:,:,:), .TRUE.)
        call material_rw3('c45.dat', C(17,:,:,:), .TRUE.)
        call material_rw3('c46.dat', C(18,:,:,:), .TRUE.)
        call material_rw3('c55.dat', C(19,:,:,:), .TRUE.)
        call material_rw3('c56.dat', C(20,:,:,:), .TRUE.)
        call material_rw3('c66.dat', C(21,:,:,:), .TRUE.)
        
        ! Load viscoelastic coefficients
        call material_rw3('gamma_x.dat',  gamma(1,:,:,:), .TRUE.)
        call material_rw3('gamma_y.dat',  gamma(2,:,:,:), .TRUE.)
        call material_rw3('gamma_z.dat',  gamma(3,:,:,:), .TRUE.)
        call material_rw3('gamma_yz.dat', gamma(4,:,:,:), .TRUE.)
        call material_rw3('gamma_xz.dat', gamma(5,:,:,:), .TRUE.)
        call material_rw3('gamma_xy.dat', gamma(6,:,:,:), .TRUE.)
        
        call material_rw3('density.dat', rho, .TRUE.)        
    end subroutine load_elastic_coefficients3
    
    
    ! ----------------------------------------------------------------------
    subroutine seismic2_pseudospectral_rk4(domain, source, SINGLE_OUTPUT)
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT
        
        real(real64), allocatable :: C(:,:,:), rho(:,:), gamma(:,:,:)
        
        type(spectral_grid_t) :: grid
        type(sponge_layer_t) :: sponge

        ! State vector Q(5, nx, nz)
        real(real64), allocatable :: Q(:,:,:)

        ! RK4 stage memory
        real(real64), allocatable :: Q_stage(:,:,:), k1(:,:,:), k2(:,:,:), k3(:,:,:), k4(:,:,:)
        ! real(real64), allocatable :: kx(:), kz(:)
        integer :: nx, nz, it
        real(real64) :: dx, dz, dt
        
        type(Source_Type) :: src
        integer :: n_spec
        real(real64), allocatable :: freq_spec(:)
        complex(real64), allocatable :: F_spec(:)
        
        ! Derivative arrays
        ! dV_dX(i,k,velocity_component,derivative_direction)
        ! Components: 1 = vx, 2 = vz
        ! Directions: 1 = x,  2 = z
        ! dS_dX(i,k,momentum_component,derivative_direction)
        !
        ! (:,:,1,1) = d(sxx)/dx
        ! (:,:,1,2) = d(sxz)/dz
        ! (:,:,2,1) = d(sxz)/dx
        ! (:,:,2,2) = d(szz)/dz
        real(real64), allocatable :: dV_dX(:,:,:,:), dS_dX(:,:,:,:)
        
        logical :: block_output, legacy_output, SINGLE
        integer :: steps_per_block, block_count, ii, jj 
        character(len=256) :: current_block_file
        
        ! -------------------------------------------------------------------------
        nx = domain%nx
        nz = domain%nz
        dx = domain%dx
        dz = domain%dz
        dt = source%dt
        
        ! ----------------------------------------------------------------------
        ! Allocations
        allocate(Q(5, nx, nz), Q_stage(5, nx, nz)) 
        allocate( k1(5, nx, nz), k2(5, nx, nz), k3(5, nx, nz), k4(5, nx, nz))
        allocate(C(6, nx, nz), gamma(3, nx, nz), rho(nx, nz) )
        allocate(dV_dX(nx,nz,2,2))
        allocate(dS_dX(nx,nz,2,2))
        
        call init_spectral_grid(grid, nx, 1, nz, dx, 1.0_real64, dz)
        call init_sponge_layer2(sponge, domain, 3.0_real64)
        
        ! ---------------------------- Source Terms ----------------------------
        ! Load coefficients 
        call load_elastic_coefficients2(nx, nz, C, gamma, rho) 
        if (any(rho <= 0.0_real64)) then
            error stop 'Elastic density must be positive.'
        end if

        call load_spectrum( 'freq_spectrum.dat', 'F_spectrum.dat', &
            freq_spec, F_spec, n_spec)

        call init_source_elastic_xz( src, source, domain, grid, &
            F_spec, freq_spec, n_spec)

        deallocate(freq_spec, F_spec)
        
        ! ----------------------------------------------------------------------
        Q = 0.0_real64
        
        ! Determine output precision flag
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        end if

        ! Setup I/O parameters and initialize block structures
        call setup_io_params_2d(nx, nz, domain%npml, block_output, steps_per_block, legacy_output)
        if (block_output) call init_io_2d(nx, nz, domain%npml, steps_per_block)
        block_count = 0
        
        ! 
        do it = 1, source%time_steps
            
            call rhs_stage_2d(Q, k1, 0.0_real64)
            
            Q_stage = Q + 0.5_real64 * dt * k1
            call rhs_stage_2d(Q_stage, k2, 0.5_real64) 
            
            Q_stage = Q + 0.5_real64 * dt * k2
            call rhs_stage_2d(Q_stage, k3, 0.5_real64)
            
            Q_stage = Q + dt * k3
            call rhs_stage_2d(Q_stage, k4, 1.0_real64)
            
            Q = Q + (dt / 6.0_real64) * (k1 + 2.0_real64 * k2 + 2.0_real64 * k3 + k4)
            
            ! Apply sponge damping
            call elastic_sponge_damping_2d(sponge, domain, Q)
            
            if (mod(it, 100) == 0) then
                print *, "Pseudospectral RK4 time step:", it, "of", source%time_steps
            end if
            

            ! -------- Output --------
            if (block_output) then
                if (current_step_in_block == 0) then
                    block_count = block_count + 1
                    ii = (block_count - 1) * steps_per_block + 1
                    jj = min(block_count * steps_per_block, source%time_steps)
                    write(current_block_file, "(a, i6.6, a, i6.6, a, i0, a, i0, a)") &
                        'S2RK4.', ii, '-', jj, '.', source%xind, '.', source%zind, '.blk.zst'
                endif
                call add_step_to_block_2d(current_block_file, &
                    real(Q(1,:,:), real64), real( Q(2,:,:), real64))
            endif

            if (legacy_output) then
                call write_image2(real(Q(1,:,:), real64), nx, nz, source, it, 'Vx', SINGLE)
                call write_image2(real(Q(2,:,:), real64), nx, nz, source, it, 'Vz', SINGLE)
            endif
            
        end do
        
        if (block_output) call finalize_io(current_block_file)
        
        call free_spectral_grid(grid)
        call free_sponge_layer(sponge)
        
        deallocate(Q, Q_stage, k1, k2, k3, k4)
        deallocate(dV_dX, dS_dX)
        deallocate(C, gamma, rho)
        
        call free_source(src)
        
        contains

        subroutine rhs_stage_2d(Q_in, stage_rate, stage_fraction)

            real(real64), intent(inout) :: Q_in(5, nx, nz)
            real(real64), intent(out) :: stage_rate(5, nx, nz)
            real(real64), intent(in) :: stage_fraction

            real(real64) :: scratch(nx, nz)

            ! Q = [vx, vz, sxx, szz, sxz]

            call grad2d(grid, Q_in(1,:,:), dV_dX(:,:,1,1), dV_dX(:,:,1,2))
            call grad2d(grid, Q_in(2,:,:), dV_dX(:,:,2,1), dV_dX(:,:,2,2))

            call grad2d(grid, Q_in(3,:,:), dS_dX(:,:,1,1), scratch)
            call grad2d(grid, Q_in(4,:,:), scratch, dS_dX(:,:,2,2))
            call grad2d(grid, Q_in(5,:,:), dS_dX(:,:,2,1), dS_dX(:,:,1,2))

            call elastic2_rk4( &
                nx, nz, Q_in, dV_dX, dS_dX, C, gamma, rho, stage_rate)

            call inject_source_explicit_stage_elastic2( src, nx, nz, &
                                        stage_rate, it, stage_fraction, rho)

        end subroutine rhs_stage_2d
        
    end subroutine seismic2_pseudospectral_rk4
    
    ! ----------------------------------------------------------------------
    subroutine seismic25_pseudospectral_rk4(domain, source, SINGLE_OUTPUT)
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT
        
        type(spectral_grid_t) :: grid
        type(sponge_layer_t) :: sponge
        
        real(real64), allocatable :: C(:,:,:), rho(:,:), gamma(:,:,:)
        
        ! 2.5D State vector Q(9, nx, nz) complex wavefield:
        ! Q(1)=Vx, Q(2)=Vy, Q(3)=Vz, Q(4)=Sxx, Q(5)=Syy, Q(6)=Szz, Q(7)=Syz, Q(8)=Sxz, Q(9)=Sxy        
        complex(real64), allocatable :: Q(:,:,:)

        ! RK4 stage memory
        complex(real64), allocatable :: Q_stage(:,:,:), k1(:,:,:), k2(:,:,:), k3(:,:,:), k4(:,:,:)
        
        integer :: nx, nz, it, istat
        real(real64) :: dx, dz, dt
        character(len=32) :: env_val
        
        ! Derivative arrays
        complex(real64), allocatable :: dV_dX(:,:,:,:)
        complex(real64), allocatable :: dS_dX(:,:,:,:)
        
        logical :: block_output, legacy_output, SINGLE
        integer :: steps_per_block, block_count, ii, jj 
        character(len=256) :: current_block_file
        
        type(Source_Type) :: src
        integer :: n_spec
        real(real64), allocatable :: freq_spec(:)
        complex(real64), allocatable :: F_spec(:)
        
        real(real64) :: ky 
        complex(real64) :: imag_unit = (0.0_real64, 1.0_real64)
        
        ! -------------------------------------------------------------------------
        nx = domain%nx
        nz = domain%nz
        dx = domain%dx
        dz = domain%dz
        dt = source%dt
        
        call init_spectral_grid(grid, nx, 1, nz, dx, 1.0_real64, dz)
        call init_sponge_layer2(sponge, domain, 3.0_real64)
        
        ! ----------------------------------------------------------------------
        ! Allocations
        allocate(Q(9, nx, nz), Q_stage(9, nx, nz)) 
        allocate( k1(9, nx, nz), k2(9, nx, nz), k3(9, nx, nz), k4(9, nx, nz))
        allocate(C(21, nx, nz), gamma(6, nx, nz), rho(nx, nz) )
        allocate( dV_dX(nx, nz, 3, 3) )
        allocate( dS_dX(nx, nz, 3, 3) )

         ! ---------------------------- Source Terms ----------------------------
        ! Load coefficients 
        call load_elastic_coefficients25(nx, nz, C, gamma, rho)
        
        if (any(rho <= 0.0_real64)) then
            error stop 'Elastic density must be positive.'
        end if

        call load_spectrum( 'freq_spectrum.dat', 'F_spectrum.dat', &
            freq_spec, F_spec, n_spec)

        call init_source_elastic_xz( src, source, domain, grid, &
            F_spec, freq_spec, n_spec)

        deallocate(freq_spec, F_spec)
        
        ! ----------------------------------------------------------------------        

        Q = (0.0_real64, 0.0_real64)
        
        ! Determine output precision flag
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        end if

        ! Setup I/O parameters and initialize block structures
        call setup_io_params_2d(nx, nz, domain%npml, block_output, steps_per_block, legacy_output)
        if (block_output) call init_io_2d(nx, nz, domain%npml, steps_per_block)
        block_count = 0
        
        ! 
        do it = 1, source%time_steps
            
            call rhs_stage_25d(Q, k1, 0.0_real64)
            
            Q_stage = Q + 0.5_real64 * dt * k1
            call rhs_stage_25d(Q_stage, k2, 0.5_real64) 
            
            Q_stage = Q + 0.5_real64 * dt * k2
            call rhs_stage_25d(Q_stage, k3, 0.5_real64)
            
            Q_stage = Q + dt * k3
            call rhs_stage_25d(Q_stage, k4, 1.0_real64)
            
            Q = Q + (dt / 6.0_real64) * (k1 + 2.0_real64 * k2 + 2.0_real64 * k3 + k4)
            
            ! Apply sponge damping
            call elastic_sponge_damping_2d(sponge, domain, Q)
            
            if (mod(it, 100) == 0) then
                print *, "Pseudospectral RK4 time step:", it, "of", source%time_steps
            end if
            
            ! -------- Output --------
            if (block_output) then
                if (current_step_in_block == 0) then
                    block_count = block_count + 1
                    ii = (block_count - 1) * steps_per_block + 1
                    jj = min(block_count * steps_per_block, source%time_steps)
                    write(current_block_file, "(a, i6.6, a, i6.6, a, i0, a, i0, a)") &
                        'S2RK4.', ii, '-', jj, '.', source%xind, '.', source%zind, '.blk.zst'
                endif
                call add_step_to_block_2d(current_block_file, &
                    real(Q(1,:,:), real64), real( Q(2,:,:), real64))
            endif

            if (legacy_output) then
                call write_image2(real(Q(1,:,:), real64), nx, nz, source, it, 'Vx', SINGLE)
                call write_image2(real(Q(2,:,:), real64), nx, nz, source, it, 'Vz', SINGLE)
            endif
            
        end do
        
        if (block_output) call finalize_io(current_block_file)
        
        call free_spectral_grid(grid)
        call free_sponge_layer(sponge)
        
        deallocate(Q, Q_stage, k1, k2, k3, k4)
        deallocate(dV_dX)
        deallocate(dS_dX)
        deallocate(C, gamma, rho)
        
        call free_source(src)
        
        contains 
        
        subroutine rhs_stage_25d(Q_in, stage_rate, stage_fraction)

            complex(real64), intent(in) :: Q_in(9, nx, nz)
            complex(real64), intent(out) :: stage_rate(9, nx, nz)
            real(real64), intent(in) :: stage_fraction

            complex(real64) :: scratch(nx, nz)
            integer :: component

            ! Velocity gradients.
            do component = 1, 3

                call grad2d_complex(grid, Q_in(component,:,:), &
                                    dV_dX(:,:,component,1), &
                                    dV_dX(:,:,component,3))

                dV_dX(:,:,component,2) = &
                    imag_unit * ky * Q_in(component,:,:)

            end do

            ! x-momentum stress-divergence terms.
            call grad2d_complex(grid, Q_in(4,:,:), &
                                dS_dX(:,:,1,1), scratch)

            dS_dX(:,:,1,2) = imag_unit * ky * Q_in(9,:,:)

            ! xy stress: x derivative contributes to y momentum.
            call grad2d_complex(grid, Q_in(9,:,:), &
                                dS_dX(:,:,2,1), scratch)

            dS_dX(:,:,2,2) = imag_unit * ky * Q_in(5,:,:)

            ! yz stress: z derivative contributes to y momentum.
            call grad2d_complex(grid, Q_in(7,:,:), &
                                scratch, dS_dX(:,:,2,3))

            ! xz stress supplies two divergence terms.
            call grad2d_complex(grid, Q_in(8,:,:), &
                                dS_dX(:,:,3,1), dS_dX(:,:,1,3))

            dS_dX(:,:,3,2) = imag_unit * ky * Q_in(7,:,:)

            call grad2d_complex(grid, Q_in(6,:,:), &
                                scratch, dS_dX(:,:,3,3))

            call elastic25_rk4( &
                nx, nz, Q_in, dV_dX, dS_dX, &
                C, gamma, rho, stage_rate)
            call inject_source_explicit_stage_elastic25( src, nx, nz, &
                stage_rate, it, stage_fraction, rho, ky)
        end subroutine rhs_stage_25d 
        
    end subroutine seismic25_pseudospectral_rk4
    
    ! ----------------------------------------------------------------------
    subroutine seismic3_pseudospectral_rk4(domain, source, SINGLE_OUTPUT)
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        
        logical, intent(in), optional :: SINGLE_OUTPUT
        
        real(real64), allocatable :: C(:,:,:,:), rho(:,:,:), gamma(:,:,:,:)
        
        type(spectral_grid_t) :: grid
        type(sponge_layer_t) :: sponge

        ! State vector Q(5, nx, nz)
        real(real64), allocatable :: Q(:,:,:,:)

        ! RK4 stage memory
        real(real64), allocatable :: Q_stage(:,:,:,:), k1(:,:,:,:), k2(:,:,:,:), k3(:,:,:,:), k4(:,:,:,:)
        integer :: nx, ny, nz, it
        real(real64) :: dx, dy, dz, dt
        
        type(Source_Type) :: src 
        integer :: n_spec 
        real(real64), allocatable :: freq_spec(:)
        complex(real64), allocatable :: F_spec(:)
        
        ! Derivative arrays
        ! dV_dX(i,k,velocity_component,derivative_direction)
        ! Components: 1 = vx, 2 = vy, 3 = vz
        ! Directions: 1 = x,  2 = y, 3 = z
        ! dS_dX(i,k,momentum_component,derivative_direction)
        !
        ! (:,:,1,1) = d(sxx)/dx
        ! (: :,1,2) = d(sxy)/dy
        ! (:,:,1,3) = d(sxz)/dz
        ! (:,:,2,1) = d(sxy)/dx
        ! (:,:,2,2) = d(syy)/dy
        ! (:,:,2,3) = d(syz)/dz
        ! (:,:,3,1) = d(sxz)/dx 
        ! (:,:,3,2) = d(syz)/dy
        ! (:,:,3,3) = d(szz)/dz
        real(real64), allocatable :: dV_dX(:,:,:,:,:), dS_dX(:,:,:,:,:)
        
        logical :: block_output, legacy_output, SINGLE
        integer :: steps_per_block, block_count, ii, jj 
        character(len=256) :: current_block_file
        
        ! -------------------------------------------------------------------------
        nx = domain%nx
        ny = domain%ny
        nz = domain%nz
        dx = domain%dx
        dy = domain%dy
        dz = domain%dz
        dt = source%dt
        
        call init_spectral_grid(grid, nx, ny, nz, dx, dy, dz)
        call init_sponge_layer3(sponge, domain, 3.0_real64)
        
        ! --------------------------- Source Values ----------------------------
        src = source
        call free_source(src) 
        src%force_vec = 0.0_real64 
        src%moment_tensor = 0.0_real64 
        
        call load_spectrum('freq_spectrum.dat', 'F_spectrum.dat', &
                            freq_spec, F_spec, n_spec)
        
        
        src%xind = source%xind + domain%npml
        src%yind = source%yind + domain%npml
        src%zind = source%zind + domain%npml

        src%x = source%x + real(domain%npml,real64)*domain%dx
        src%y = source%y + real(domain%npml,real64)*domain%dy
        src%z = source%z + real(domain%npml,real64)*domain%dz
        
        select case (trim(src%source_type))
            case ('ac')
                call init_source_weight_drop( &
                    src, domain, F_spec, freq_spec, n_spec, grid)

            case ('tnt')
                call init_source_explosive( &
                    src, F_spec, freq_spec, n_spec, grid, domain)

            case ('dc')
                call init_source_double_couple( &
                    src, domain, F_spec, freq_spec, n_spec, grid)

            case ('clvd')
                call init_source_clvd( &
                    src, domain, F_spec, freq_spec, n_spec, grid)

            case ('pw')
                error stop 'Elastic plane-wave injection is not implemented here.'

            case default
                error stop 'Unknown 3D elastic source type.'

        end select

        deallocate(freq_spec, F_spec)
        ! ----------------------------------------------------------------------
        
        ! Allocations
        allocate(Q(9, nx, ny, nz), Q_stage(9, nx, ny, nz)) 
        allocate( k1(9, nx, ny, nz), k2(9, nx, ny, nz), k3(9, nx, ny, nz), k4(9, nx, ny, nz))
        allocate(C(21, nx, ny, nz), gamma(6, nx, ny, nz), rho(nx, ny, nz) )
        allocate(dV_dX(nx, ny,nz,3,3))
        allocate(dS_dX(nx, ny,nz,3,3))
        
        ! Load coefficients 
        call load_elastic_coefficients3(nx, ny, nz, C, gamma, rho) 
        
        if (any(rho <= 0.0_real64)) then
            error stop '3D elastic density must be positive.'
        end if
        
        Q = 0.0_real64
        
        ! Determine output precision flag
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        end if

        ! Setup I/O parameters and initialize block structures
        call setup_io_params_3d(nx, ny, nz, domain%npml, block_output, steps_per_block, legacy_output)
        if (block_output) call init_io_3d(nx, ny, nz, domain%npml, steps_per_block)
        block_count = 0
        
        ! 
        do it = 1, source%time_steps
            
            call rhs_stage_3d(Q, k1, 0.0_real64)
            
            Q_stage = Q + 0.5_real64 * dt * k1
            call rhs_stage_3d(Q_stage, k2, 0.5_real64) 
            
            Q_stage = Q + 0.5_real64 * dt * k2
            call rhs_stage_3d(Q_stage, k3, 0.5_real64)
            
            Q_stage = Q + dt * k3
            call rhs_stage_3d(Q_stage, k4, 1.0_real64)
            
            Q = Q + (dt / 6.0_real64) * (k1 + 2.0_real64 * k2 + 2.0_real64 * k3 + k4)
            
            ! Apply sponge damping
            call elastic_sponge_damping_3d(sponge, domain, Q)
            
            if (mod(it, 100) == 0) then
                print *, "Pseudospectral RK4 time step:", it, "of", source%time_steps
            end if
            

            ! -------- Output --------
            if (block_output) then
                if (current_step_in_block == 0) then
                    block_count = block_count + 1
                    ii = (block_count - 1) * steps_per_block + 1
                    jj = min(block_count * steps_per_block, source%time_steps)
                    write(current_block_file, &
                            "(a,i6.6,a,i6.6,a,i0,a,i0,a,i0,a)") &
                            'S3RK4.', ii, '-', jj, '.', &
                            source%xind, '.', source%yind, '.', source%zind, &
                            '.blk.zst'
                endif
                call add_step_to_block_3d( current_block_file, &
                                Q(1,:,:,:), Q(2,:,:,:), Q(3,:,:,:))
            endif

            if (legacy_output) then
                call write_image3( Q(1,:,:,:), nx, ny, nz, source, it, 'Vx', SINGLE)
                call write_image3( Q(2,:,:,:), nx, ny, nz, source, it, 'Vy', SINGLE)
                call write_image3( Q(3,:,:,:), nx, ny, nz, source, it, 'Vz', SINGLE)
            endif
            
        end do
        
        if (block_output) call finalize_io(current_block_file)
        
        call free_spectral_grid(grid)
        call free_sponge_layer(sponge)
        
        deallocate(Q, Q_stage, k1, k2, k3, k4)
        deallocate(dV_dX, dS_dX)
        deallocate(C, gamma, rho)
        
        contains

        subroutine rhs_stage_3d(Q_in, stage_rate, stage_fraction)

            real(real64), intent(inout) :: Q_in(9,nx,ny,nz)
            real(real64), intent(out) :: stage_rate(9,nx,ny,nz)
            real(real64), intent(in) :: stage_fraction

            integer :: component

            ! Two distinct scratch outputs are needed because grad3d
            ! returns all three derivatives in each call.
            real(real64) :: scratch1(nx,ny,nz)
            real(real64) :: scratch2(nx,ny,nz)

            ! Q = [vx, vy, vz, sxx, syy, szz, syz, sxz, sxy]

            do component = 1, 3
                call grad3d(grid, Q_in(component,:,:,:), &
                            dV_dX(:,:,:,component,1), &
                            dV_dX(:,:,:,component,2), &
                            dV_dX(:,:,:,component,3))
            end do

            ! Normal stresses: retain their corresponding normal derivative.
            call grad3d(grid, Q_in(4,:,:,:), &
                        dS_dX(:,:,:,1,1), scratch1, scratch2)

            call grad3d(grid, Q_in(5,:,:,:), &
                        scratch1, dS_dX(:,:,:,2,2), scratch2)

            call grad3d(grid, Q_in(6,:,:,:), &
                        scratch1, scratch2, dS_dX(:,:,:,3,3))

            ! syz: y derivative -> z momentum; z derivative -> y momentum.
            call grad3d(grid, Q_in(7,:,:,:), &
                        scratch1, dS_dX(:,:,:,3,2), dS_dX(:,:,:,2,3))

            ! sxz: x derivative -> z momentum; z derivative -> x momentum.
            call grad3d(grid, Q_in(8,:,:,:), &
                        dS_dX(:,:,:,3,1), scratch1, dS_dX(:,:,:,1,3))

            ! sxy: x derivative -> y momentum; y derivative -> x momentum.
            call grad3d(grid, Q_in(9,:,:,:), &
                        dS_dX(:,:,:,2,1), dS_dX(:,:,:,1,2), scratch1)

            call elastic3_rk4( &
                nx, ny, nz, Q_in, dV_dX, dS_dX, &
                C, gamma, rho, stage_rate)
            
            call inject_source_explicit_stage_elastic3(src, nx, ny, nz, &
                                stage_rate, it, stage_fraction, rho)
            
        end subroutine rhs_stage_3d
        
    end subroutine seismic3_pseudospectral_rk4
    
end module pseudospectralfdtd