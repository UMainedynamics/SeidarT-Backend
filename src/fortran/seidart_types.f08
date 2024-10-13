module seidart_types
    use iso_fortran_env, only: real64
    
    implicit none 
    
    public :: Domain_Type, Time_Parameters_Type
    public :: Attenuation_Type, Stiffness_Type 
    public :: Source_Type, Permittivity_Type, Conductivity_Type
    
    ! -------------------------- Define Types ---------------------------------- 
    !I/O Types
    ! Domain parameters
    type :: Domain_Type 
        real(real64) :: dim 
        integer :: nx, ny, nz
        real(real64) :: dx, dy, dz 
        integer :: cpml, nmats 
        character(len=:), allocatable :: image_file 
    end type Domain_Type
    
    ! Time parameters for Seismic and Electromagnetic sources
    type :: Time_Parameters_Type
        real(real64) :: dt
        integer :: time_steps
    end type Time_Parameters_Type
    
    ! Seismic Source
    type :: Source_Type
        real(real64) :: x, y, z
        integer :: xind, yind, zind
        real(real64) :: source_frequency
        real(real64) :: x_z_rotation, x_y_rotation
        real(real64) :: amplitude
        character(len=:), allocatable :: source_type
    end type Source_Type
    
    ! Seismic attenuation properties
    type :: Attenuation_Type
        integer :: id
        character(len=:), allocatable :: name
        real(real64) :: alpha_x, alpha_y, alpha_z
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
    
    
    ! CPMLFDTD Types 
    type :: Seismic2_Variables_Type
        real(real64), allocatable :: memdvx_dx(:,:)
        real(real64), allocatable :: memdvx_dz(:,:)
        real(real64), allocatable :: memdvz_dx(:,:)
        real(real64), allocatable :: memdvz_dz(:,:)
        real(real64), allocatable :: memdsigxx_dx (:,:)
        real(real64), allocatable :: memdsigzz_dz(:,:)
        real(real64), allocatable :: memdsigxz_dx(:,:)
        real(real64), allocatable :: memdsigxz_dz(:,:)
        real(real64), allocatable :: vx(:,:)
        real(real64), allocatable :: vz(:,:)
        real(real64), allocatable :: sigxx(:,:)
        real(real64), allocatable :: sigzz(:,:)
        real(real64), allocatable :: sigxz(:,:)
    end type Seismic2_Variables_Type

    type :: Seismic3_Variables_Type
        real(real64), allocatable :: memdvx_dx(:,:,:)
        real(real64), allocatable :: memdvx_dy(:,:,:)
        real(real64), allocatable :: memdvx_dz(:,:,:)
        real(real64), allocatable :: memdvy_dx(:,:,:)
        real(real64), allocatable :: memdvy_dy(:,:,:)
        real(real64), allocatable :: memdvy_dz(:,:,:)
        real(real64), allocatable :: memdvz_dx(:,:,:)
        real(real64), allocatable :: memdvz_dy(:,:,:)
        real(real64), allocatable :: memdvz_dz(:,:,:)
        real(real64), allocatable :: memdsigxx_dx (:,:,:)
        real(real64), allocatable :: memdsigyy_dy(:,:,:)
        real(real64), allocatable :: memdsigzz_dz(:,:,:)
        real(real64), allocatable :: memdsigxy_dx(:,:,:)
        real(real64), allocatable :: memdsigxy_dy(:,:,:)
        real(real64), allocatable :: memdsigxz_dx(:,:,:)
        real(real64), allocatable :: memdsigxz_dz(:,:,:)
        real(real64), allocatable :: memdsigyz_dy(:,:,:)
        real(real64), allocatable :: memdsigyz_dz(:,:,:)
        real(real64), allocatable :: vx(:,:,:)
        real(real64), allocatable :: vy(:,:,:) 
        real(real64), allocatable :: vz(:,:,:)
        real(real64), allocatable :: sigxx(:,:,:)
        real(real64), allocatable :: sigyy(:,:,:)
        real(real64), allocatable :: sigzz(:,:,:)
        real(real64), allocatable :: sigxy(:,:,:)
        real(real64), allocatable :: sigxz(:,:,:)
        real(real64), allocatable :: sigyz(:,:,:)
    end type Seismic3_Variables_Type

    type :: Electromagnetic2_Variables_Type
        real(real64), allocatable :: memdEz_dx(:,:), memdEx_dz(:,:)
        real(real64), allocatable :: memdHy_dx(:,:), memdHy_dz(:,:)
        real(real64), allocatable :: eps11(:,:), eps33(:,:)
        real(real64), allocatable :: eps13(:,:), sig11(:,:)
        real(real64), allocatable :: sig33(:,:), sig13(:,:)
        real(real64), allocatable :: Ex(:,:), Ez(:,:), Hy(:,:)
    end type Electromagnetic2_Variables_Type
    
    type :: Electromagnetic3_Variables_Type
        real(real64), allocatable :: memdEy_dx(:,:,:), memdEx_dy(:,:,:)
        real(real64), allocatable :: memdEz_dx(:,:,:), memdEx_dz(:,:,:)
        real(real64), allocatable :: memdEy_dz(:,:,:), memdEz_dy(:,:,:)        
        real(real64), allocatable :: memdHz_dx(:,:,:), memdHx_dz(:,:,:)
        real(real64), allocatable :: memdHy_dx(:,:,:), memdHx_dy(:,:,:)
        real(real64), allocatable :: memdHy_dz(:,:,:), memdHz_dy(:,:,:)        
        real(real64), allocatable :: eps11(:,:), eps22(:,:), eps33(:,:)
        real(real64), allocatable :: eps12(:,:), eps13(:,:), eps23(:,:)
        real(real64), allocatable :: sig11(:,:), sig22(:,:), sig33(:,:)
        real(real64), allocatable :: sig12(:,:), sig13(:,:), sig23(:,:)
        real(real64), allocatable :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
        real(real64), allocatable :: Hx(:,:,:), Hy(:,:,:), Hz(:,:,:)
    end type Electromagnetic3_Variables_Type

end module seidart_types