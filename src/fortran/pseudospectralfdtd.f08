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
        ny = domain%ny
        nz = domain%nz
        dt = source%dt
        
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