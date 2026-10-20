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
    
    implicit none
    include 'fftw3.f03'
    
    private
    public :: seismic2_pseudospectral_rk4, seismic25_pseudospectral_rk4
    public :: electromag2_pseudospectral_rk4, electromag25_pseudospectral_rk4
    !public :: seismic3_pseudospectral_rk4, electromag3_pseudos
    
    ! =========================================================================
    subroutine load_elastic_coefficients2(nx, nz, C, gamma, rho)
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
    end subroutine load_elastic_coefficients2
    
    ! ------------------------------------------------------------------------
    subroutine load_elastic_coefficients3(nx, ny, nz, C, gamma, rho)
        integer, intent(in) :: nx, ny, nz
        real(real64), intent(out) :: C(21, nx, ny, nz)
        real(real64), intent(out) :: gamma(6, nx, ny, nz)
        real(real64), intent(out) :: rho(nx, ny, nz)
        integer, intent(out) :: dim(3)
        
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
        type(Domain), intent(in) :: domain
        type(Source), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT
        
        real(real64), allocatable :: C(:,:,:), rho(:,:), gamma(:,:,:)
        
        ! State vector Q(5, nx, nz)
        real(real64), allocatable :: Q(:,:,:)

        ! RK4 stage memory
        real(real64), allocatable :: Q_stage(:,:,:), k1(:,:,:), k2(:,:,:), k3(:,:,:), k4(:,:,:)
        
        integer :: nx, nz
        real(real64) :: dx, dz, dt
        
        ! Derivative arrays
        real(real64), allocatable :: dvx_dx(:,:), dvx_dz(:,:)
        real(real64), allocatable :: dvz_dx(:,:), dvz_dz(:,:)
        real(real64), allocatable :: dsxx_dx(:,:), dszz_dz(:,:)
        real(real64), allocatable :: dsxz_dx(:,:), dsxz_dz(:,:)
        
        logical :: block_output, legacy_output, SINGLE
        integer :: steps_per_block, block_count, ii, jj 
        character(len=256) :: current_block_file
        
        ! -------------------------------------------------------------------------
        nx = domain%nx
        nz = domain%nz
        dx = domain%dx
        dz = domain%dz
        dt = source%dt
        
        call init_spectral_grid(grid, nx, 1, nz, dx, 1.0_real64, dz)
        call init_sponge_layer(sponge, domain, 3.0_real64)
        
        ! Allocations
        allocate(Q(5, nx, nz), Q_stage(5, nx, nz)) 
        allocate( k1(5, nx, nz), k2(5, nx, nz), k3(5, nx, nz), k4(5, nx, nz))
        allocate(C(6, nx, nz), gamma(3, nx, nz), rho(nx, nz) )
        allocate(dvx_dx(nx, nz), dvx_dz(nx, nz), dvz_dx(nx, nz), dvz_dz(nx, nz))
        allocate(dsxx_dx(nx, nz), dszz_dz(nx, nz), dsxz_dx(nx, nz), dsxz_dz(nx, nz))
        allocate(kz(nx), kx(nx))
        
        ! Load coefficients 
        call material_rw2('c11.dat', C(1, :,:), .TRUE.)
        call material_rw2('c13.dat', C(2, :,:), .TRUE.)
        call material_rw2('c15.dat', C(3, :,:), .TRUE.)
        call material_rw2('c33.dat', C(4, :,:), .TRUE.)
        call material_rw2('c35.dat', C(5, :,:), .TRUE.)
        call material_rw2('c55.dat', C(6, :,:), .TRUE.)
        call material_rw2('gamma_x.dat',  gamma(1, :,:), .TRUE.)
        call material_rw2('gamma_z.dat',  gamma(2, :,:), .TRUE.)
        call material_rw2('gamma_xz.dat', gamma(3, :,:), .TRUE.)
        call material_rw2('density.dat', rho, .TRUE.)   
        
        Q = 0.0_real64
        
        ! Determine output precision flag
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        end if

        ! Setup I/O parameters and initialize block structures
        call setup_io_params_2d(nx, nz, domain%cpml, block_output, steps_per_block, legacy_output)
        if (block_output) call init_io_2d(nx, nz, domain%cpml, steps_per_block)
        block_count = 0
        
        ! 
        do it = 1, source%time_steps
            
            call evalute_rhs_2d(Q, k1, 0.0_real64)
            
            Q_stage = Q + 0.5_real64 * dt * k1
            call evalute_rhs_2d(Q_stage, k2, 0.5_real64) 
            
            Q_stage = Q + 0.5_real64 * dt * k2
            call evalute_rhs_2d(Q_stage, k3, 0.5_real64)
            
            Q_stage = Q + dt * k3
            call evalute_rhs_2d(Q_stage, k4, 1.0_real64)
            
            Q = Q + (dt / 6.0_real64) * (k1 + 2.0_real64 * k2 + 2.0_real64 * k3 + k4)
            
            ! Apply sponge damping
            call apply_sponge_damping_2d(sponge, domain, Q)
            
            if (mod(it, 100) == 0 .and. verbose) then
                print *, "Pseudospectral RK4 time step:", it, "of", source%time_steps
            end if
            
#ifdef SEIDART_OPENMP_GPU
        !$omp target update from(vx, vz)
#endif

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
        deallocate(dvx_dx, dvx_dz, dvz_dx, dvz_dz)
        deallocate(dsxx_dx, dszz_dz, dsxz_dx, dsxz_dz)
        deallocate(C, gamma, rho)
        
    end subroutine seismic2_pseudospectral_rk4
    
        ! ----------------------------------------------------------------------
    subroutine seismic25_pseudospectral_rk4(domain, source, SINGLE_OUTPUT)
        type(Domain), intent(in) :: domain
        type(Source), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT
        
        real(real64), allocatable :: C(:,:,:), rho(:,:), gamma(:,:,:)
        
        ! 2.5D State vector Q(9, nx, nz) complex wavefield:
        ! Q(1)=Vx, Q(2)=Vy, Q(3)=Vz, Q(4)=Sxx, Q(5)=Syy, Q(6)=Szz, Q(7)=Syz, Q(8)=Sxz, Q(9)=Sxy        
        complex(real64), allocatable :: Q(:,:,:)

        ! RK4 stage memory
        complex(real64), allocatable :: Q_stage(:,:,:), k1(:,:,:), k2(:,:,:), k3(:,:,:), k4(:,:,:)
        
        real(real64), allocatable :: kx(:), kz(:)
        integer :: nx, nz, it, istat
        real(real64) :: dx, dz, dt
        character(len=32) :: env_val
        
        ! Derivative arrays
        real(real64), allocatable :: dvx_dx(:,:), dvx_dz(:,:)
        real(real64), allocatable :: dvz_dx(:,:), dvz_dz(:,:)
        real(real64), allocatable :: dsxx_dx(:,:), dszz_dz(:,:)
        real(real64), allocatable :: dsxz_dx(:,:), dsxz_dz(:,:)
        
        logical :: block_output, legacy_output, SINGLE
        integer :: steps_per_block, block_count, ii, jj 
        character(len=256) :: current_block_file
        
        ! -------------------------------------------------------------------------
        nx = domain%nx
        nz = domain%nz
        dx = domain%dx
        dz = domain%dz
        dt = source%dt
        
        call init_spectral_grid(grid, nx, 1, nz, dx, 1.0_real64, dz)
        call init_sponge_layer(sponge, domain, 3.0_real64)
        
        ! Allocations
        allocate(Q(9, nx, nz), Q_stage(9, nx, nz)) 
        allocate( k1(9, nx, nz), k2(9, nx, nz), k3(9, nx, nz), k4(9, nx, nz))
        allocate(C(21, nx, nz), gamma(6, nx, nz), rho(nx, nz) )
        allocate(dvx_dx(nx, nz), dvx_dz(nx, nz), dvz_dx(nx, nz), dvz_dz(nx, nz))
        allocate(dsxx_dx(nx, nz), dszz_dz(nx, nz), dsxz_dx(nx, nz), dsxz_dz(nx, nz))
        allocate(kz(nx), kx(nx))
        
        ! Load coefficients 
        call material_rw2('c11.dat', C(1, :,:), .TRUE.)
        call material_rw2('c12.dat', C(2, :,:), .TRUE.)
        call material_rw2('c13.dat', C(3, :,:), .TRUE.)
        call material_rw2('c14.dat', C(4, :,:), .TRUE.)
        call material_rw2('c15.dat', C(5, :,:), .TRUE.)
        call material_rw2('c16.dat', C(6, :,:), .TRUE.)
        call material_rw2('c22.dat', C(7, :,:), .TRUE.)
        call material_rw2('c23.dat', C(8, :,:), .TRUE.)
        call material_rw2('c24.dat', C(9, :,:), .TRUE.)
        call material_rw2('c25.dat', C(10, :,:), .TRUE.)
        call material_rw2('c26.dat', C(11, :,:), .TRUE.)
        call material_rw2('c33.dat', C(12, :,:), .TRUE.)
        call material_rw2('c34.dat', C(13, :,:), .TRUE.)
        call material_rw2('c35.dat', C(14, :,:), .TRUE.)
        call material_rw2('c36.dat', C(15, :,:), .TRUE.)
        call material_rw2('c44.dat', C(16, :,:), .TRUE.)
        call material_rw2('c45.dat', C(17, :,:), .TRUE.)
        call material_rw2('c46.dat', C(18, :,:), .TRUE.)
        call material_rw2('c55.dat', C(19, :,:), .TRUE.)
        call material_rw2('c56.dat', C(20, :,:), .TRUE.)
        call material_rw2('c66.dat', C(21, :,:), .TRUE.)
        call material_rw2('gamma_x.dat',  gamma(1, :,:), .TRUE.)
        call material_rw2('gamma_y.dat',  gamma(2, :,:), .TRUE.)
        call material_rw2('gamma_z.dat',  gamma(3, :,:), .TRUE.)
        call material_rw2('gamma_yz.dat', gamma(4, :,:), .TRUE.)
        call material_rw2('gamma_xz.dat', gamma(5, :,:), .TRUE.)
        call material_rw2('gamma_xy.dat', gamma(6, :,:), .TRUE.)
        call material_rw2('density.dat', rho, .TRUE.)   
        
        Q = (0.0_real64, 0.0_real64)
        
        ! Determine output precision flag
        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        end if

        ! Setup I/O parameters and initialize block structures
        call setup_io_params_2d(nx, nz, domain%cpml, block_output, steps_per_block, legacy_output)
        if (block_output) call init_io_2d(nx, nz, domain%cpml, steps_per_block)
        block_count = 0
        
        ! 
        do it = 1, source%time_steps
            
            call evalute_rhs_2d(Q, k1, 0.0_real64)
            
            Q_stage = Q + 0.5_real64 * dt * k1
            call evalute_rhs_2d(Q_stage, k2, 0.5_real64) 
            
            Q_stage = Q + 0.5_real64 * dt * k2
            call evalute_rhs_2d(Q_stage, k3, 0.5_real64)
            
            Q_stage = Q + dt * k3
            call evalute_rhs_2d(Q_stage, k4, 1.0_real64)
            
            Q = Q + (dt / 6.0_real64) * (k1 + 2.0_real64 * k2 + 2.0_real64 * k3 + k4)
            
            ! Apply sponge damping
            call apply_sponge_damping_2d(sponge, domain, Q)
            
            if (mod(it, 100) == 0 .and. verbose) then
                print *, "Pseudospectral RK4 time step:", it, "of", source%time_steps
            end if
            
#ifdef SEIDART_OPENMP_GPU
        !$omp target update from(vx, vz)
#endif

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
        deallocate(dvx_dx, dvx_dz, dvz_dx, dvz_dz)
        deallocate(dsxx_dx, dszz_dz, dsxz_dx, dsxz_dz)
        deallocate(C, gamma, rho)
        
    end subroutine seismic25_pseudospectral_rk4
    
    
end module pseudospectralfdtd