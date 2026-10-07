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
    public :: seismic2_pseudospectral_rk4, seismic3_pseudospectral_rk4, seismic25_pseudospectral_rk4
    
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
    subroutine seismic2_pseudospectral_rk4(domain, source, density_method, verbose)
        type(Domain), intent(in) :: domain
        type(Source), intent(in) :: source
        character(len=*), intent(in) :: density_method
        logical, intent(in) :: verbose
        
        real(real64), allocatable :: c11(:,:), c13(:,:), c15(:,:)
        real(real64), allocatable :: c33(:,:), c35(:,:), c55(:,:)
        real(real64), allocatable :: rho(:,:)
        real(real64), allocatable :: gamma_x(:,:), gamma_z(:,:), gamma_xz(:,:)
        
        ! State vector Q(5, nx, nz)
        real(real64), allocatable :: Q(:,:,:)

        ! RK4 stage memory
        real(real64), allocatable :: Q_stage(:,:,:), k1(:,:,:), k2(:,:,:), k3(:,:,:), k4(:,:,:)
        
        ! Derivative arrays
        real(real64), allocatable :: dvx_dx(:,:), dvx_dz(:,:)
        real(real64), allocatable :: dvz_dx(:,:), dvz_dz(:,:)
        real(real64), allocatable :: dsxx_dx(:,:), dszz_dz(:,:)
        real(real64), allocatable :: dsxz_dx(:,:), dsxz_dz(:,:)
        
        ! -------------------------------------------------------------------------
        nx = domain%nx
        nz = domain%nz
        dx = domain%dx
        dz = domain%dz
        dt = source%dt
        
        call init_spectral_grid(grid, nx, 1, nz, dx, 1.0_real64, dz)
        
        ! Allocations
        allocate(Q(5, nx, nz), Q_stage(5, nx, nz), k1(5, nx, nz), k2(5, nx, nz), k3(5, nx, nz), k4(5, nx, nz))
        allocate(c11(nx, nz), c13(nx, nz), c15(nx, nz), c33(nx, nz), c35(nx, nz), c55(nx, nz))
        allocate(rho(nx, nz))
        allocate(gamma_x(nx, nz), gamma_z(nx, nz), gamma_xz(nx, nz))
        allocate(dvx_dx(nx, nz), dvx_dz(nx, nz), dvz_dx(nx, nz), dvz_dz(nx, nz))
        allocate(dsxx_dx(nx, nz), dszz_dz(nx, nz), dsxz_dx(nx, nz), dsxz_dz(nx, nz))
        allocate(kz(nx), kx(nx))
        
        ! Initialize Fourier wavenumbers kx, kz
        call init_wavenumbers(nx, nz, dx, dz, kx, kz)
        
        do it = 1, source%time_steps
            ! Compute derivatives in Fourier space
            call compute_derivatives(Q, kx, kz, dvx_dx, dvx_dz, dvz_dx, dvz_dz, dsxx_dx, dszz_dz, dsxz_dx, dsxz_dz)
            
            ! RK4 stages
            call rk4_step(Q, Q_stage, k1, k2, k3, k4, dt, c11, c13, c15, c33, c35, c55, rho)
            
            ! Apply source term
            call apply_source(Q_stage, source%source_time(it), source%source_position)
            
            ! Update state vector
            Q = Q_stage
            
            ! Output or visualization if needed
            if (mod(it, 100) == 0 .and. verbose) then
                print *, "Time step:", it
            end if
        end do
        
    end subroutine seismic2_pseudospectral_rk41
    
    
    
end module pseudospectralfdtd