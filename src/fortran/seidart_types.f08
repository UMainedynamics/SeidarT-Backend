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
        integer :: cpml, nmats 
        character(len=:), allocatable :: image_file 
        ! character(len=256) :: image_file 
    end type Domain_Type
    
    ! Seismic and Electromagnetic Source
    type :: Source_Type
        real(real64) :: dt 
        integer :: time_steps
        real(real64) :: x, y, z
        integer :: xind, yind, zind
        real(real64) :: source_frequency
        real(real64) :: x_z_rotation, x_y_rotation
        real(real64) :: amplitude
        character(len=:), allocatable :: source_type
        ! character(len=5) :: source_type
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
    
end module seidart_types