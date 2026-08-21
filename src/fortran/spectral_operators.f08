module spectral_operators

    use iso_c_binding
    use iso_fortran_env, only: real64
    implicit none
    include 'fftw3.f03'
    public :: spectral_grid_t, init_spectral_grid, free_spectral_grid, grad3d 
    
    ! -------------------------------------------------------------------------
    type :: spectral_grid_t 
        integer :: nx, ny, nz, nkx
        real(real64) :: inv_n_total 
        real(real64), allocatable :: kx(:), ky(:), kz(:)
        type(c_ptr) :: plan_fwd = c_null_ptr 
        type(c_ptr) :: plan_bwd_x = c_null_ptr 
        type(c_ptr) :: plan_bwd_y = c_null_ptr 
        type(c_ptr) :: plan_bwd_z = c_null_ptr
        complex(real64), allocatable :: F_hat(:,:,:)
        complex(real64), allocatable :: dFx_hat(:,:,:), dFy_hat(:,:,:), &
                                        dFz_hat(:,:,:)
    end type spectral_grid_t
    
    contains
    
    subroutine init_spectral_grid(grid, nx, ny, nz, dx, dy, dz)
        type(spectral_grid_t), intent(inout) :: grid 
        integer, intent(in) :: nx, ny, nz
        real(real64), intent(in) :: dx, dy, dz
        
        real(real64), allocatable, target :: dummy_in(:,:,:), dummy_out(:,:,:)
        real(real64), parameter :: PI = 3.14159265358979323846_real64
        integer :: i,j,k
        
        grid%nx = nx; grid%ny = ny; grid%nz = nz
        grid%nkx = nx / 2 + 1
        grid%inv_n_total = 1.0_real64 / real(nx * ny * nz, real64) 
        
        allocate(grid%kx(grid%nkx), grid%ky(ny), grid%kz(nz))
        allocate(grid%F_hat(grid%nkx, ny, nz))
        allocate(grid%dFx_hat(grid%nkx, ny, nz))
        allocate(grid%dFy_hat(grid%nkx, ny, nz))
        allocate(grid%dFz_hat(grid%nkx, ny, nz))
        
        do i = 1,grid%nkx 
            grid%kx(i) = 2.0_real64 * PI * real(i - 1, real64) / &
                                    (real(nx, real64) * dx)
        end do
        do j = 1,ny
            grid%ky(j) = merge( 2.0_real64 * PI * real(j - 1, real64), &
                                2.0_real64 * PI * real(j - 1 - ny, real64), &
                                j <= ny/2 + 1) / &
                                    (real(ny, real64) * dy )
        end do
        do k = 1,nz
            grid%kz(k) = merge( 2.0_real64 * PI * real(k - 1, real64), &
                                2.0_real64 * PI * real(k - 1 - nz, real64), &
                                k <= nz/2 + 1 ) / &
                                    (real(nz, real64) * dz )
        end do
        
        allocate(dummy_in(nx, ny, nz), dummy_out(nx, ny, nz) )
        grid%plan_fwd = fftw_plan_dft_r2c_3d(nz, ny, nx, dummy_in, grid%F_hat, FFTW_ESTIMATE)
        grid%plan_bwd_x = fftw_plan_dft_c2r_3d(nz, ny, nx, grid%dFx_hat, dummy_out, FFTW_ESTIMATE)
        grid%plan_bwd_y = fftw_plan_dft_c2r_3d(nz, ny, nx, grid%dFy_hat, dummy_out, FFTW_ESTIMATE)
        grid%plan_bwd_z = fftw_plan_dft_c2r_3d(nz, ny, nx, grid%dFz_hat, dummy_out, FFTW_ESTIMATE)
        deallocate(dummy_in, dummy_out)
        
    end subroutine init_spectral_grid

    ! -------------------------------------------------------------------------
    subroutine grad3d(grid, F, dF_dx, dF_dy, dF_dz)
        type(spectral_grid_t), intent(inout) :: grid 
        real(real64), intent(inout) :: F(grid%nx, grid%ny, grid%nz)
        real(real64), intent(out) :: dF_dx(grid%nx, grid%ny, grid%nz)
        real(real64), intent(out) :: dF_dy(grid%nx, grid%ny, grid%nz)
        real(real64), intent(out) :: dF_dz(grid%nx, grid%ny, grid%nz)
        
        complex(real64), parameter :: imag_unit = (0.0_real64, 1.0_real64)
        integer :: i,j,k
        
        call fftw_execute_dft_r2c(grid%plan_fwd, F, grid%F_hat)
        
        !$omp parallel do collapse(3) private(i,j,k) 
        do k = 1, grid%nz
            do j = 1, grid%ny
                do i = 1, grid%nkx
                    grid%dFx_hat(i, j, k) = (imag_unit * grid%kx(i) * &
                            grid%inv_n_total) * grid%F_hat(i,j,k )
                    grid%dFy_hat(i, j, k) = (imag_unit * grid%ky(j) * &
                            grid%inv_n_total) * grid%F_hat(i, j, k)
                    grid%dFz_hat(i, j, k) = (imag_unit * grid%kz(k) * &
                            grid%inv_n_total) * grid%F_hat(i, j, k)
                end do
            end do
        end do
        !$omp end parallel do
        
        call fftw_execute_dft_c2r(grid%plan_bwd_x, grid%dFx_hat, dF_dx)
        call fftw_execute_dft_c2r(grid%plan_bwd_y, grid%dFy_hat, dF_dy)
        call fftw_execute_dft_c2r(grid%plan_bwd_z, grid%dFz_hat, dF_dz)
        
    end subroutine grad3d
    
    subroutine free_spectral_grid(grid)
        type(spectral_grid_t), intent(inout) :: grid
        if (c_associated(grid%plan_fwd))   call fftw_destroy_plan(grid%plan_fwd)
        if (c_associated(grid%plan_bwd_x)) call fftw_destroy_plan(grid%plan_bwd_x)
        if (c_associated(grid%plan_bwd_y)) call fftw_destroy_plan(grid%plan_bwd_y)
        if (c_associated(grid%plan_bwd_z)) call fftw_destroy_plan(grid%plan_bwd_z)
        if (allocated(grid%kx)) deallocate(grid%kx, grid%ky, grid%kz)
        if (allocated(grid%F_hat)) deallocate(grid%F_hat, grid%dFx_hat, grid%dFy_hat, grid%dFz_hat)
    end subroutine free_spectral_grid
    ! -------------------------------------------------------------------------

end module spectral_operators