module seidartio 

    use iso_fortran_env, only: real64
    use json_module 
    
    implicit none 
    
    public :: parse_json, Domain_Type, Material_Type, Time_Parameters_Type
    public :: Seismic_Source_Type, Attenuation_Type, Stiffness_Type 
    public :: Electromagnetic_Source_Type, Permittivity_Type, Conductivity_Type
    
    ! private :: 
    
    ! ---------------------------- Declarations --------------------------------
    ! Domain parameters
    type :: Domain_Type 
        integer :: dim 
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
        real(real64) :: source_frequency
        real(real64) :: x_z_rotation, x_y_rotation
        real(real64) :: amplitude
        character(len=:), allocatable :: source_type
    end type Seismic_Source_Type

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
    
    ! --------------------- Subroutine Definitions -----------------------------
    contains 
    
    subroutine parse_json(file_name, domain, seismic_source, 
                        seismic_stiffness, seismic_attenuation, &
                        electromagnetic_permittivity, electromagnetic_conductivity, &
                        electromagnetic_source
                        )
    
        implicit none 
        character(len=*), intent(in) :: file_name
        type(Domain_Type), intent(out) :: domain
        type(Time_Parameters_Type), intent(out) :: seismic_time_parameters 
        type(Time_Parameters_Type), intent(out) :: electromagnetic_time_parameters
        type(Source_Type), intent(out) :: seismic_source
        type(Source_Type), intent(out) :: electromagnetic_source
        type(Stiffness_Type), intent(out), allocatable :: seismic_stiffness(:)
        type(Attenuation_Type), intent(out), allocatable :: seismic_attenuation(:)
        type(Permittivity_Type), intent(out), allocatable :: electromagnetic_permittivity(:)
        type(Conductivity_Type), intent(out), allocatable :: electromagnetic_conductivity(:)
        

        type(json_file) :: json 
        integer :: i 
        
        ! --------------------- End Declarations -------------------------------
        
        call json%initialize()
        call json%load_file(trim(file_name))
        
        ! Parse the domain data
        call json%get('Domain/dim', domain%dim)
        call json%get('Domain/nx', domain%nx)
        call json%get('Domain/ny', domain%ny)
        call json%get('Domain/nz', domain%nz)
        call json%get('Domain/dx', domain%dx)
        call json%get('Domain/dy', domain%dy)
        call json%get('Domain/dz', domain%dz)
        call json%get('Domain/cpml', domain%cpml)
        call json%get('Domain/nmats', domain%nmats)
        call json%get('Domain/image_file', domain%image_file)
    
        
        ! Parse the Seismic Source
        call json%get('Seismic/Source/x', seismic_source%x)
        call json%get('Seismic/Source/y', seismic_source%y)
        call json%get('Seismic/Source/z', seismic_source%z)
        call json%get('Seismic/Source/source_frequency', seismic_source%source_frequency)
        call json%get('Seismic/Source/x-z_rotation', seismic_source%x_z_rotation)
        call json%get('Seismic/Source/x-y_rotation', seismic_source%x_y_rotation)
        call json%get('Seismic/Source/amplitude', seismic_source%amplitude)
        call json%get('Seismic/Source/source_type', seismic_source%source_type)

        ! Allocate and parse all tensor coefficients. This one's a doozy.
        allocate(seismic_stiffness(domain%nmats))
        allocate(seismic_attenuation(domain%nmats))
        allocate(electromagnetic_permittivity(domain%nmats))
        allocate(electromagnetic_conductivity(domain%nmats))
        
        do i = 1, domain%nmats
            ! Seismic
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/id', &
                            seismic_stiffness(i)%id)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c11', &
                            seismic_stiffness(i)%c11)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c12', &
                            seismic_stiffness(i)%c12)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c13', &
                            seismic_stiffness(i)%c13)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c14', &
                            seismic_stiffness(i)%c14)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c15', &
                            seismic_stiffness(i)%c15)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c16', &
                            seismic_stiffness(i)%c16)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c22', &
                            seismic_stiffness(i)%c22)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c23', &
                            seismic_stiffness(i)%c23)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c24', &
                            seismic_stiffness(i)%c24)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c25', &
                            seismic_stiffness(i)%c25)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c26', &
                            seismic_stiffness(i)%c26)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c33', &
                            seismic_stiffness(i)%c33)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c34', &
                            seismic_stiffness(i)%c34)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c35', &
                            seismic_stiffness(i)%c35)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c36', &
                            seismic_stiffness(i)%c36)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c44', &
                            seismic_stiffness(i)%c44)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c45', &
                            seismic_stiffness(i)%c45)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c46', &
                            seismic_stiffness(i)%c46)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c55', &
                            seismic_stiffness(i)%c55)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c56', &
                            seismic_stiffness(i)%c56)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i))//')/c66', &
                            seismic_stiffness(i)%c66) 
            ! EM      
            call json%get('Electromagnetic/Permittivity_Coefficients('//trim(adjustl(i))//')/e11',&
                            electromagnetic_permittivity(i)%e11)         
            call json%get('Electromagnetic/Permittivity_Coefficients('//trim(adjustl(i))//')/e12',&
                            electromagnetic_permittivity(i)%e12)         
            call json%get('Electromagnetic/Permittivity_Coefficients('//trim(adjustl(i))//')/e13',&
                            electromagnetic_permittivity(i)%e13)    
            call json%get('Electromagnetic/Permittivity_Coefficients('//trim(adjustl(i))//')/e22',&
                            electromagnetic_permittivity(i)%e22)   
            call json%get('Electromagnetic/Permittivity_Coefficients('//trim(adjustl(i))//')/e23',&
                            electromagnetic_permittivity(i)%e23)   
            call json%get('Electromagnetic/Permittivity_Coefficients('//trim(adjustl(i))//')/e33',&
                            electromagnetic_permittivity(i)%e33)   
            call json%get('Electromagnetic/Conductivity_Coefficients('//trim(adjustl(i))//')/s11',&
                            electromagnetic_conductivity(i)%s11)   
            call json%get('Electromagnetic/Conductivity_Coefficients('//trim(adjustl(i))//')/s12',&
                            electromagnetic_conductivity(i)%s12)   
            call json%get('Electromagnetic/Conductivity_Coefficients('//trim(adjustl(i))//')/s13',&
                            electromagnetic_conductivity(i)%s13)
            call json%get('Electromagnetic/Conductivity_Coefficients('//trim(adjustl(i))//')/s22',&
                            electromagnetic_conductivity(i)%s22)   
            call json%get('Electromagnetic/Conductivity_Coefficients('//trim(adjustl(i))//')/s23',&
                            electromagnetic_conductivity(i)%s23)   
            call json%get('Electromagnetic/Conductivity_Coefficients('//trim(adjustl(i))//')/s33',&
                            electromagnetic_conductivity(i)%s33)   
            
        end do

        call json%finalize()
    end subroutine parse_json
    
    ! --------------------------------------------------------------------------
    subroutine read_geometry(file_name, domain, geometry)
    
        implicit none 
        
        character(len=*), intent(in) :: file_name 
        type(Domain_Type), intent(in) :: domain 
        integer, allocatable, intent(out) :: geometry(:,:)
        
        integer :: i,j 
        integer :: nx, nz 
        integer : io_status 
        integer :: unit_number 
        
        ! Extract dimensions from the domain structure
        nx = domain%nx 
        nz = domain%nz 
        
        ! Allocate the 2D array based on nx and nz
        allocate(geometry(nx, nz))
        
        ! Open the geometry.dat file for reading
        open(nuewunit=unit_number, file=file_name, status = 'old', action = 'read', &
            iostat = io_status)
        
        if (io_status /= 0) then 
            print *, "Error opening file: ", file_name 
            stop 
        end if
        
        ! Read the nx-by-nz array from the file
        do j = 1, nz
            do i = 1, nx
            read(unit_number, *, iostat=io_status) geometry(i, j)
            if (io_status /= 0) then
                print *, "Error reading geometry data"
                stop
            end if
            end do
        end do
        
        ! Close the file 
        close(unit_number)
    end subroutine read_geometry

    ! --------------------------------------------------------------------------
    subroutine seismic_parameter_write(domain, geometry, stiffness, attenuation)
        implicit none 
        
        type(Domain_Type), intent(in) :: domain 
        integer, intent(in) :: geometry(:,:)
        type(Stiffness_Type), intent(in) :: stiffness 
        type(Attenuation_Type), intent(in) :: attenuation 
        
        real(real64), allocatable :: extended_stiffness(:,:,:)
        real(real64), allocatable :: extended_attenuation(:,:,:)
        integer :; nx, nz, cpml 
        integer :: i, j, id_value 
        integer :; nx_ext, nz_ext 
        integer :: unit_number, io_status 
        
        character(len=7) :: stiffness_fn
        character(len=10) :: attenuation_fn
        character(len=3), dimension(21) :: scoef_names = (/ &
                                            'c11', 'c12', 'c13', 'c14', 'c15', &
                                            'c16', 'c22', 'c23', 'c24', 'c25', &
                                            'c26', 'c33', 'c34', 'c35', 'c36', &
                                            'c44', 'c45', 'c46', 'c55', 'c56', &
                                            'c66' /)
        
        character(len=6), dimension(3) :: acoef_names = (/ 'gammax', 'gammay', 'gammaz'/)
        
        ! Extract domain dimensions and cpml
        nx = domain%nx
        nz = domain%nz
        cpml = domain%cpml

        ! Calculate extended dimensions
        nx_ext = 2 * cpml + nx
        nz_ext = 2 * cpml + nz

        ! Allocate the extended stiffness array with extra CPML padding
        allocate(extended_stiffness(nx_ext, nz_ext, 21))
        allocate(extended_attenuation(nx_ext, nz_ext, 3))  
        
        extended_stiffness(:,:,:) = 0.0 
        extended_attenuation(:,:,:) = 0.0

         ! Loop over the geometry array and assign stiffness coefficients
        do j = 1, nz
            do i = 1, nx
                id_value = geometry(i, j)  ! Get the material id from geometry

                ! Assign stiffness coefficients based on id_value from geometry
                extended_stiffness(i + cpml, j + cpml, 1) = stiffness%c11(id_value)
                extended_stiffness(i + cpml, j + cpml, 2) = stiffness%c12(id_value)
                extended_stiffness(i + cpml, j + cpml, 3) = stiffness%c13(id_value)
                extended_stiffness(i + cpml, j + cpml, 4) = stiffness%c14(id_value)
                extended_stiffness(i + cpml, j + cpml, 5) = stiffness%c15(id_value)
                extended_stiffness(i + cpml, j + cpml, 6) = stiffness%c16(id_value)
                extended_stiffness(i + cpml, j + cpml, 7) = stiffness%c22(id_value)
                extended_stiffness(i + cpml, j + cpml, 8) = stiffness%c23(id_value)
                extended_stiffness(i + cpml, j + cpml, 9) = stiffness%c24(id_value)
                extended_stiffness(i + cpml, j + cpml, 10) = stiffness%c25(id_value)
                extended_stiffness(i + cpml, j + cpml, 11) = stiffness%c26(id_value)
                extended_stiffness(i + cpml, j + cpml, 12) = stiffness%c33(id_value)
                extended_stiffness(i + cpml, j + cpml, 13) = stiffness%c34(id_value)
                extended_stiffness(i + cpml, j + cpml, 14) = stiffness%c35(id_value)
                extended_stiffness(i + cpml, j + cpml, 15) = stiffness%c36(id_value)
                extended_stiffness(i + cpml, j + cpml, 16) = stiffness%c44(id_value)
                extended_stiffness(i + cpml, j + cpml, 17) = stiffness%c45(id_value)
                extended_stiffness(i + cpml, j + cpml, 18) = stiffness%c46(id_value)
                extended_stiffness(i + cpml, j + cpml, 19) = stiffness%c55(id_value)
                extended_stiffness(i + cpml, j + cpml, 20) = stiffness%c56(id_value)
                extended_stiffness(i + cpml, j + cpml, 21) = stiffness%c66(id_value)
                ! same for attenuation
                extended_attenuation(i + cpml, j + cpml, 1) = attenuation%alpha_x(id_value)
                extended_attenuation(i + cpml, j + cpml, 2) = attenuation%alpha_y(id_value)
                extended_attenuation(i + cpml, j + cpml, 3) = attenuation%alpha_z(id_value)
            end do
        end do
        
        ! Extend the values at 1 + cpml to the 1 to cpml indices
        do k = 1, 21
          ! Extend along x-direction (rows)
          do j = 1, nz_ext
            do i = 1, cpml
              extended_stiffness(i, j, k) = extended_stiffness(1 + cpml, j, k)
              extended_attenuation(i, j, k) = extended_attenuation(1 + cpml, j, k)
              extended_stiffness(nx_ext - i + 1, j, k) = extended_stiffness(nx_ext - cpml, j, k)
              extended_attenuation(nx_ext - i + 1, j, k) = extended_attenuation(nx_ext - cpml, j, k)
            end do
          end do
        
          ! Extend along z-direction (columns)
          do i = 1, nx_ext
            do j = 1, cpml
              extended_stiffness(i, j, k) = extended_stiffness(i, 1 + cpml, k)
              extended_attenuation(i, j, k) = extended_attenuation(i, 1 + cpml, k)
              extended_stiffness(i, nz_ext - j + 1, k) = extended_stiffness(i, nz_ext - cpml, k)
              extended_attenuation(i, nz_ext - j + 1, k) = extended_attenuation(i, nz_ext - cpml, k)
            end do
          end do
        end do
        
        ! Write each 2D array (stiffness coefficients) to separate files
        do k = 1, 21
            ! Create filename based on coefficient name (e.g., 'c11.dat')
            write(stiffness_fn, '(A3, ".dat")') scoef_names(k)

            ! Open the file to write the 2D array
            open(newunit=unit_number, file=stiffness_fn, form='unformatted', &
                access='stream', status='replace', iostat=io_status)
            if (io_status /= 0) then
                print *, "Error opening file for writing: ", stiffness_fn
                stop
            end if

            ! Write the 2D array to the binary file
            write(unit_number) extended_stiffness(:,:,k)

            ! Close the file
            close(unit_number)
        end do
        
        do k = 1, 3
            ! Create filename based on coefficient name (e.g., 'c11.dat')
            write(attenuation_fn, '(A3, ".dat")') acoef_names(k)

            ! Open the file to write the 2D array
            open(newunit=unit_number, file=attenuation_fn, form='unformatted', &
                access='stream', status='replace', iostat=io_status)
            if (io_status /= 0) then
                print *, "Error opening file for writing: ", attenuation_fn
                stop
            end if

            ! Write the 2D array to the binary file
            write(unit_number) extended_attenuation(:,:,k)

            ! Close the file
            close(unit_number)
        end do
        
        deallocate(extended_stiffness)
        deallocate(extended_attenuation)
    end subroutine seismic_parameter_write
    
    ! --------------------------------------------------------------------------
    subroutine electromagnetic_parameter_write(domain, geometry, permittivity, conductivity)
        implicit none
        
        type(Domain_Type), intent(in) :: domain 
        integer, intent(in) :: geometry(:,:) 
        type(Permittivity_Type), intent(in) :: permittivity 
        type(Conductivity_Type), intent(in) :: conductivity 
        
        real(real64), allocatable :: extended_permittivity(:,:,:)
        real(real64), allocatable :: extended_conductivity(:,:,:)
        integer :: nx, nz, cpml 
        integer :: i, j, id_value 
        integer nx_ext, nz_ext 
        integer :: unit_number, io_status 
        
        character(len=9) :: permittivity_fn, conductivity_fn 
        character(len=5), dimension(6) :: pcoef_names = (/ &
                                'eps11','eps12','eps13','eps22','eps23','eps33')
        character(len=5), dimension(6) :: scoef_names = (/ &
                                'sig11','sig12','sig13','sig22','sig23','sig33')
        
        ! Extract domain dimensions and cpml
        nx = domain%nx
        nz = domain%nz
        cpml = domain%cpml

        ! Calculate extended dimensions
        nx_ext = 2 * cpml + nx
        nz_ext = 2 * cpml + nz
        
        ! Allocate the extended permittivity array with extra CPML padding
        allocate(extended_permittivity(nx_ext, nz_ext, 21))
        allocate(extended_conductivity(nx_ext, nz_ext, 3))  
        
        extended_permittivity(:,:,:) = 0.0 
        extended_conductivity(:,:,:) = 0.0
        
        do j = 1,nz 
            do i = 1,nx 
                id_value = geometry(i,j) 
                
                extended_permittivity(i + cpml, j + cpml, 1) = permittivity%e11(id_value)
                extended_permittivity(i + cpml, j + cpml, 2) = permittivity%e12(id_value)
                extended_permittivity(i + cpml, j + cpml, 3) = permittivity%e13(id_value)
                extended_permittivity(i + cpml, j + cpml, 4) = permittivity%e22(id_value)
                extended_permittivity(i + cpml, j + cpml, 5) = permittivity%e23(id_value)
                extended_permittivity(i + cpml, j + cpml, 6) = permittivity%e33(id_value)

                extended_conductivity(i + cpml, j + cpml, 1) = conductivity%s11(id_value)
                extended_conductivity(i + cpml, j + cpml, 2) = conductivity%s12(id_value)
                extended_conductivity(i + cpml, j + cpml, 3) = conductivity%s13(id_value)
                extended_conductivity(i + cpml, j + cpml, 4) = conductivity%s22(id_value)
                extended_conductivity(i + cpml, j + cpml, 5) = conductivity%s23(id_value)
                extended_conductivity(i + cpml, j + cpml, 6) = conductivity%s33(id_value)
            end do 
        end do 
        
        ! Extend the values at 1 + cpml to the 1 to cpml indices
        do k = 1, 6
          ! Extend along x-direction (rows)
          do j = 1, nz_ext
            do i = 1, cpml
              extended_permittivity(i, j, k) = extended_permittivity(1 + cpml, j, k)
              extended_conductivity(i, j, k) = extended_conductivity(1 + cpml, j, k)
              extended_permittivity(nx_ext - i + 1, j, k) = extended_permittivity(nx_ext - cpml, j, k)
              extended_conductivity(nx_ext - i + 1, j, k) = extended_conductivity(nx_ext - cpml, j, k)
            end do
          end do
        
          ! Extend along z-direction (columns)
          do i = 1, nx_ext
            do j = 1, cpml
              extended_permittivity(i, j, k) = extended_permittivity(i, 1 + cpml, k)
              extended_conductivity(i, j, k) = extended_conductivity(i, 1 + cpml, k)
              extended_permittivity(i, nz_ext - j + 1, k) = extended_permittivity(i, nz_ext - cpml, k)
              extended_conductivity(i, nz_ext - j + 1, k) = extended_conductivity(i, nz_ext - cpml, k)
            end do
          end do
        end do
        
        do k = 1, 6
            ! Create filename based on coefficient name (e.g., 'c11.dat')
            write(permittivity_fn, '(A3, ".dat")') pcoef_names(k)

            ! Open the file to write the 2D array
            open(newunit=unit_number, file=permittivity_fn, form='unformatted', &
                access='stream', status='replace', iostat=io_status)
            if (io_status /= 0) then
                print *, "Error opening file for writing: ", permittivity_fn
                stop
            end if

            ! Write the 2D array to the binary file
            write(unit_number) extended_permittivity(:,:,k)

            ! Close the file
            close(unit_number)
        end do
         do k = 1, 3
            ! Create filename based on coefficient name (e.g., 'c11.dat')
            write(conductivity_fn, '(A3, ".dat")') scoef_names(k)

            ! Open the file to write the 2D array
            open(newunit=unit_number, file=conductivity_fn, form='unformatted', &
                access='stream', status='replace', iostat=io_status)
            if (io_status /= 0) then
                print *, "Error opening file for writing: ", conductivity_fn
                stop
            end if

            ! Write the 2D array to the binary file
            write(unit_number) extended_conductivity(:,:,k)

            ! Close the file
            close(unit_number)
        end do
        
        deallocate(extended_permittivity)
        deallocate(extended_conductivity)
        
    end subroutine electromagnetic_parameter_write
    
    ! --------------------------------------------------------------------------
    subroutine loadsource(time_steps, filename, srcfn)
        
        implicit none
        
        character(len=*) :: filename
        integer :: time_steps
        real(real64),dimension(time_steps) :: srcfn
        integer :: unit_number, io_status
        
        
        open(newunit = unit_number, file = trim(filename), status="old", &
            form="unformatted", action = "read", iostat=io_status)
        
        if (io_status /= 0) then 
            print *, 'Error opening file: ', filename
            stop 
        end if
        
        read(unit_number) srcfn 
        
        close(unit_number)

    end subroutine loadsource
    
    ! --------------------------------------------------------------------------
    subroutine material_rw(filename, image_data, readfile)

        implicit none
        
        ! Arguments
        character(len=*), intent(in) :: filename
        real(real64), dimension(:,:), intent(inout) :: image_data
        logical, intent(in) :: readfile
        
        ! Local variables
        integer :: unit_number
        integer :: io_status

        ! Open the file with a dynamic unit number
        open(newunit=unit_number, file=trim(filename), form="unformatted", &
            status="old", action="readwrite", iostat=io_status)
        
        ! Check if the file opened successfully
        if (io_status /= 0) then
            print *, 'Error opening file: ', filename
            stop
        end if
        
        ! Read or write based on the value of readfile
        if (readfile) then
            read(unit_number) image_data
        else
            write(unit_number) image_data
        end if
        
        ! Close the file
        close(unit_number)

    end subroutine material_rw
    
    ! --------------------------------------------------------------------------
    subroutine write_image(image_data, nx, ny, nz, src, it, channel, SINGLE)

        implicit none

        ! Arguments
        real(real64), intent(in) :: image_data(..)
        integer, intent(in) :: nx, ny, nz, it
        integer, dimension(:), intent(in) :: src
        character(len=2), intent(in) :: channel
        logical, intent(in) :: SINGLE

        ! Local variables
        character(len=100) :: filename
        integer :: unit_number

        ! Write the filename based on the rank of the src array
        if (size(src) == 2) then
            WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a1, a4)" ) &
                    channel, '.', it, '.', src(1), '.', src(2), '.dat'
        else if (size(src) == 3) then
            WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a1, i0, a4)" ) &
                    channel, '.', it, '.', src(1), '.', src(2), '.', src(3), '.dat'
        else
            print *, "Error: src array must have 2 or 3 elements."
            stop
        end if

        ! Open the file with a dynamic unit number
        open(newunit=unit_number, form="unformatted", file=trim(filename))

        ! Write based on the SINGLE flag and rank of the image_data
        if (SINGLE) then
            if (rank(image_data) == 2) then
                write(unit_number) sngl(image_data(1:nx, 1:nz))
            else if (rank(image_data) == 3) then
                write(unit_number) sngl(image_data(1:nx, 1:ny, 1:nz))
            endif
        else
            if (rank(image_data) == 2) then
                write(unit_number) image_data(1:nx, 1:nz)
            else if (rank(image_data) == 3) then
                write(unit_number) image_data(1:nx, 1:ny, 1:nz)
            endif
        end if

        ! Close the file
        close(unit_number)

    end subroutine write_image

    ! --------------------------------------------------------------------------
    subroutine loadcpml(filename, image_data)

        use iso_fortran_env, only: real64
        implicit none

        ! Arguments
        character(len=*), intent(in) :: filename
        real(real64), dimension(:), intent(out) :: image_data

        ! Local variables
        integer :: unit_number
        integer :: io_status

        ! Open the file with a dynamic unit number and access the file as a stream
        open(newunit=unit_number, file=trim(filename), form="unformatted", access='stream', &
            status="old", action="read", iostat=io_status)

        ! Check if the file opened successfully
        if (io_status /= 0) then
            print *, 'Error opening file: ', filename
            stop
        end if

        ! Read the image data from the file
        read(unit_number) image_data

        ! Close the file
        close(unit_number)

    end subroutine loadcpml

end module seidartio
    
    