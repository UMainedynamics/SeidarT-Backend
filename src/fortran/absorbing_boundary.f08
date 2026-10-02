module absorbing_boundary
    
    use iso_fortran_env, only: real64
    implicit none
    public :: sponge_layer_t, init_sponge_layer, free_sponge_layer, apply_sponge_damping

    type :: sponge_layer_t
        integer :: nx, ny, nz
        integer :: npml
        real(real64), allocatable :: damping_mask(:,:,:)
    end type sponge_layer_t

    contains
    
    ! --------------------------------------------------------------------------
    subroutine init_sponge_layer(sponge, domain, alpha)
        type(sponge_layer_t), intent(inout) :: sponge
        type(Domain_Type), intent(in) :: domain
        real(real64), intent(in) :: alpha
        
        real(real64), allocatable :: wx(:), wy(:), wz(:)
        real(real64) :: dist
        integer :: i, j, k
        
        if (allocated(sponge%damping_mask)) deallocate(sponge%damping_mask)
        allocate(sponge%damping_mask(domain%nx, domain%ny, domain%nz))
        allocate(wx(domain%nx), wy(domain%ny), wz(domain%nz))

        ! 1. Precompute 1D damping profiles (1.0 in interior, exp decay in layer)
        wx = 1.0_real64
        do i = 1, domain%npml
            dist = real(domain%npml - i + 1, real64)
            wx(i) = exp(-(alpha * dist / real(domain%npml, real64))**2)
            wx(domain%nx - i + 1) = wx(i)
        end do

        wy = 1.0_real64
        do j = 1, domain%npml
            dist = real(domain%npml - j + 1, real64)
            wy(j) = exp(-(alpha * dist / real(domain%npml, real64))**2)
            wy(domain%ny - j + 1) = wy(j)
        end do

        wz = 1.0_real64
        do k = 1, domain%npml
            dist = real(domain%npml - k + 1, real64)
            wz(k) = exp(-(alpha * dist / real(domain%npml, real64))**2)
            wz(domain%nz - k + 1) = wz(k)
        end do

        ! 2. Tensor product to build 3D damping mask
        do k = 1, domain%nz
            do j = 1, domain%ny
                do i = 1, domain%nx
                    sponge%damping_mask(i, j, k) = wx(i) * wy(j) * wz(k)
                end do
            end do
        end do

        deallocate(wx, wy, wz)
    end subroutine init_sponge_layer
    
    ! --------------------------------------------------------------------------
    ! Multiplies state vector and fluid pressure by damping mask at end of time step
    subroutine apply_sponge_damping(sponge, domain, Q, fluid_pressure)
        type(sponge_layer_t), intent(in) :: sponge
        type(Domain_Type), intent(in) :: domain
        real(real64), intent(inout) :: Q(21, domain%nx, domain%ny, domain%nz)
        real(real64), intent(inout) :: fluid_pressure(domain%nx, domain%ny, domain%nz)
        
        integer :: i, j, k, v
        integer :: nx, ny, nz
        real(real64) :: w

        nx = domain%nx
        ny = domain%ny
        nz = domain%nz

        !$omp target teams distribute parallel do collapse(3) &
        !$omp map(to: sponge%damping_mask) &
        !$omp map(tofrom: Q, fluid_pressure)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    w = sponge%damping_mask(i, j, k)
                    if (w < 0.999999_real64) then
                        do v = 1, 21
                            Q(v, i, j, k) = Q(v, i, j, k) * w
                        end do
                        fluid_pressure(i, j, k) = fluid_pressure(i, j, k) * w
                    end if
                end do
            end do
        end do
    end subroutine apply_sponge_damping

    ! --------------------------------------------------------------------------
    ! Destructor to release heap memory allocated for the sponge mask
    subroutine free_sponge_layer(sponge)
        type(sponge_layer_t), intent(inout) :: sponge
        if (allocated(sponge%damping_mask)) deallocate(sponge%damping_mask)
    end subroutine free_sponge_layer

end module absorbing_boundary