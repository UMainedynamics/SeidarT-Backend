module seidart_types
    use iso_fortran_env, only: real64
    
    implicit none 
    
    public :: Domain_Type
    public :: Attenuation_Type, Stiffness_Type 
    public :: Source_Type, Permittivity_Type, Conductivity_Type
    
    ! -------------------------- Define Types ---------------------------------- 
    !I/O Types
    ! Domain parameters
    type :: Domain_Type 
        real(real64) :: dim 
        integer :: nx, ny, nz
        real(real64) :: dx, dy, dz 
        integer :: npml, nmats 
        character(len=:), allocatable :: image_file 
        ! character(len=256) :: image_file 
    end type Domain_Type
    
    ! Seismic and Electromagnetic Source
    type :: Source_Type
        real(real64) :: dt 
        integer :: time_steps
        real(real64) :: x, y, z
        integer :: xind, yind, zind
        integer :: i1, i2, j1, j2, k1, k2
        real(real64) :: source_frequency
        real(real64) :: dip, azimuth, strike, rake, plunge
        real(real64) :: amplitude
        integer :: half_span
        character(len=:), allocatable :: source_type
        character(len=:), allocatable :: source_wavelet
        real(real64), allocatable :: time_series(:)
        real(real64), allocatable :: spatial_kernel(:,:,:)
        ! Directional Force components (For AWD)
        real(real64) :: force_vec(3)
        
        ! Symmetric Moment Tensor (For Explosive, DC, CLVD)
        real(real64) :: moment_tensor(6) 
        
        ! Plane wave parameters
        real(real64) :: p_dir(3)        ! Propagation direction unit vector
        real(real64) :: e_pol(3)        ! Particle polarization unit vector
        real(real64) :: c_phase         ! Background wave phase velocity (m/s)
        real(real64) :: r0_ref(3)       ! Reference entry point for t=0 phase
        integer :: pml_thick            ! Thickness of sponge boundary layer
        real(real64), allocatable :: time_delay_3d(:,:,:) ! Spatial propagation delay tau(x,y,z)
        logical, allocatable :: injection_mask(:,:,:)     ! Active TFSF boundary cells
    end type Source_Type
    
    ! Seismic attenuation properties
    type :: Attenuation_Type
        integer :: id
        character(len=:), allocatable :: name
        real(real64) :: alpha_x, alpha_xy, alpha_xz, alpha_y, alpha_yz, alpha_z
        real(real64) :: reference_frequency
    end type Attenuation_Type
    
    ! Seismic stiffness coefficients
    type :: Stiffness_Type
        integer :: id
        real(real64) :: c11, c12, c13, c14, c15, c16
        real(real64) :: c22, c23, c24, c25, c26
        real(real64) :: c33, c34, c35, c36
        real(real64) :: c44, c45, c46
        real(real64) :: c55, c56
        real(real64) :: c66
        real(real64) :: density
    end type Stiffness_Type
    
    ! Electromagnetic permittivity properties
    type :: Permittivity_Type
        integer :: id
        real(real64) :: e11, e12, e13
        real(real64) :: e22, e23, e33
    end type Permittivity_Type
    
    ! Electromagnetic conductivity properties
    type :: Conductivity_Type
        integer :: id
        real(real64) :: s11, s12, s13
        real(real64) :: s22, s23, s33
    end type Conductivity_Type
    
        type :: spectral_grid_t 
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
    
    
end module seidart_types