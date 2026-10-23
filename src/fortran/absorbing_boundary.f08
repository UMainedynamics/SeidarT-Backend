module absorbing_boundary
    
    use iso_fortran_env, only: real64
    use seidart_types
    implicit none
    
    
    public :: sponge_layer_t, init_sponge_layer2,init_sponge_layer3, free_sponge_layer
    public :: biot_sponge_damping3d, &
                elastic_sponge_damping_2d, elastic_sponge_damping_3d

    type :: sponge_layer_t
        integer :: nx, ny, nz
        integer :: npml
        real(real64), allocatable :: damping_mask(:,:,:)
    end type sponge_layer_t

    contains
    
    ! --------------------------------------------------------------------------
    subroutine init_sponge_layer2(sponge, domain, alpha)
        type(sponge_layer_t), intent(inout) :: sponge
        type(Domain_Type), intent(in)       :: domain
        real(real64), intent(in)           :: alpha
        
        real(real64), allocatable :: wx(:), wz(:)
        real(real64) :: dist
        integer :: i, k
        
        sponge%nx = domain%nx
        sponge%ny = 1
        sponge%nz = domain%nz
        sponge%npml = domain%npml

        if (allocated(sponge%damping_mask)) deallocate(sponge%damping_mask)
        ! Allocate a 2D slice in memory (ny = 1)
        allocate(sponge%damping_mask(domain%nx, 1, domain%nz))
        allocate(wx(domain%nx), wz(domain%nz))

        ! Precompute X-direction damping profile
        wx = 1.0_real64
        do i = 1, domain%npml
            dist = real(domain%npml - i + 1, real64)
            wx(i) = exp(-(alpha * dist / real(domain%npml, real64))**2)
            wx(domain%nx - i + 1) = wx(i)
        end do

        ! Precompute Z-direction damping profile
        wz = 1.0_real64
        do k = 1, domain%npml
            dist = real(domain%npml - k + 1, real64)
            wz(k) = exp(-(alpha * dist / real(domain%npml, real64))**2)
            wz(domain%nz - k + 1) = wz(k)
        end do

        ! Build 2D damping mask
        do k = 1, domain%nz
            do i = 1, domain%nx
                sponge%damping_mask(i, 1, k) = wx(i) * wz(k)
            end do
        end do

        deallocate(wx, wz)
    end subroutine init_sponge_layer2
    
    ! --------------------------------------------------------------------------
    subroutine init_sponge_layer3(sponge, domain, alpha)
        type(sponge_layer_t), intent(inout) :: sponge
        type(Domain_Type), intent(in)       :: domain
        real(real64), intent(in)           :: alpha
        
        real(real64), allocatable :: wx(:), wy(:), wz(:)
        real(real64) :: dist
        integer :: i, j, k
        
        ! Store domain dimensions and boundary width
        sponge%nx = domain%nx
        sponge%ny = domain%ny
        sponge%nz = domain%nz
        sponge%npml = domain%npml

        ! Allocate memory for 3D mask and temporary 1D profile vectors
        if (allocated(sponge%damping_mask)) deallocate(sponge%damping_mask)
        allocate(sponge%damping_mask(domain%nx, domain%ny, domain%nz))
        allocate(wx(domain%nx), wy(domain%ny), wz(domain%nz))

        ! 1. Precompute 1D damping profiles (1.0 in interior, quadratic exponential decay in sponge boundary)
        wx = 1.0_real64
        do i = 1, domain%npml
            dist = real(domain%npml - i + 1, real64)
            wx(i) = exp(-(alpha * dist / real(domain%npml, real64))**2)
            wx(domain%nx - i + 1) = wx(i)
        end do

        wy = 1.0_real64
        do j = 1, domain%ny
            if (domain%ny == 1) exit ! Guard for 2D setups where ny = 1
        end do
        wy = 1.0_real64
        do j = 1, domain%npml
            if (domain%ny <= 1) exit
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

        ! 2. Compute 3D tensor product mask
        do k = 1, domain%nz
            do j = 1, domain%ny
                do i = 1, domain%nx
                    sponge%damping_mask(i, j, k) = wx(i) * wy(j) * wz(k)
                end do
            end do
        end do

        deallocate(wx, wy, wz)
    end subroutine init_sponge_layer3
    ! --------------------------------------------------------------------------
    ! 2D adaptation of sponge damping
    subroutine elastic_sponge_damping_2d(sp, dom, Q_vec)
        type(sponge_layer_t), intent(in) :: sp
        type(Domain_Type), intent(in)   :: dom
        real(real64), intent(inout)     :: Q_vec(5, dom%nx, dom%nz)
        integer :: i, k, v
        real(real64) :: w

        do k = 1, dom%nz
            do i = 1, dom%nx
                w = sp%damping_mask(i, 1, k)
                if (w < 0.999999_real64) then
                    do v = 1, 5
                        Q_vec(v, i, k) = Q_vec(v, i, k) * w
                    end do
                end if
            end do
        end do
    end subroutine elastic_sponge_damping_2d
    
    ! --------------------------------------------------------------------------    
    ! Multiplies state vector and fluid pressure by damping mask at end of time step
    subroutine elastic_sponge_damping_3d(sponge, domain, Q)
        type(sponge_layer_t), intent(in) :: sponge
        type(Domain_Type), intent(in) :: domain
        real(real64), intent(inout) :: Q(9, domain%nx, domain%ny, domain%nz)
        
        integer :: i, j, k, v
        real(real64) :: w

        !$omp parallel do collapse(3) schedule(static)
        do k = 1, domain%nz
            do j = 1, domain%ny
                do i = 1, domain%nx
                    w = sponge%damping_mask(i, j, k)
                    if (w < 0.999999_real64) then
                        do v = 1, 21
                            Q(v, i, j, k) = Q(v, i, j, k) * w
                        end do                    
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine elastic_sponge_damping_3d

    ! --------------------------------------------------------------------------    
    ! Multiplies state vector and fluid pressure by damping mask at end of time step
    subroutine biot_sponge_damping_3d(sponge, domain, Q, fluid_pressure)
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
    end subroutine biot_sponge_damping_3d

    
    ! --------------------------------------------------------------------------
    ! Destructor to release heap memory allocated for the sponge mask
    subroutine free_sponge_layer(sponge)
        type(sponge_layer_t), intent(inout) :: sponge
        if (allocated(sponge%damping_mask)) deallocate(sponge%damping_mask)
    end subroutine free_sponge_layer

end module absorbing_boundary