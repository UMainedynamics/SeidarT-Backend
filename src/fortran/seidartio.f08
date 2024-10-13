module seidartio 

    use iso_fortran_env, only: real64, real32
    use json_module
    use seidart_types 
    
    implicit none 
    
    public :: parse_json
    
    ! private :: 
    
    ! ---------------------------- Declarations --------------------------------
    
    
    ! --------------------- Subroutine Definitions -----------------------------
    contains 
    
    subroutine parse_json(file_name, domain, seismic_source, &
                        seismic_stiffness, seismic_attenuation, &
                        electromagnetic_permittivity, electromagnetic_conductivity, &
                        electromagnetic_source)
    
        implicit none 
        character(len=*), intent(in) :: file_name
        type(Domain_Type), intent(out) :: domain
        ! type(Time_Parameters_Type), intent(out) :: seismic_time_parameters 
        ! type(Time_Parameters_Type), intent(out) :: electromagnetic_time_parameters
        type(Source_Type), intent(out) :: seismic_source
        type(Source_Type), intent(out) :: electromagnetic_source
        type(Stiffness_Type), intent(out), allocatable :: seismic_stiffness(:)
        type(Attenuation_Type), intent(out), allocatable :: seismic_attenuation(:)
        type(Permittivity_Type), intent(out), allocatable :: electromagnetic_permittivity(:)
        type(Conductivity_Type), intent(out), allocatable :: electromagnetic_conductivity(:)
        

        type(json_file) :: json 
        integer :: i 
        character(len=10) :: i_str
        
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
        
        ! Parse the EM Source
        call json%get('Electromagnetic/Source/x', electromagnetic_source%x)
        call json%get('Electromagnetic/Source/y', electromagnetic_source%y)
        call json%get('Electromagnetic/Source/z', electromagnetic_source%z)
        call json%get('Electromagnetic/Source/source_frequency', electromagnetic_source%source_frequency)
        call json%get('Electromagnetic/Source/x-z_rotation', electromagnetic_source%x_z_rotation)
        call json%get('Electromagnetic/Source/x-y_rotation', electromagnetic_source%x_y_rotation)
        call json%get('Electromagnetic/Source/amplitude', electromagnetic_source%amplitude)
        call json%get('Electromagnetic/Source/source_type', electromagnetic_source%source_type)

        ! Allocate and parse all tensor coefficients. This one's a doozy.
        allocate(seismic_stiffness(domain%nmats))
        allocate(seismic_attenuation(domain%nmats))
        allocate(electromagnetic_permittivity(domain%nmats))
        allocate(electromagnetic_conductivity(domain%nmats))
        
        do i = 1, domain%nmats
             write(i_str, '(I0)') i  ! Convert integer i to a string without leading zeros
             
            ! Seismic
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/id', &
                            seismic_stiffness(i)%id)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c11', &
                            seismic_stiffness(i)%c11)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c12', &
                            seismic_stiffness(i)%c12)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c13', &
                            seismic_stiffness(i)%c13)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c14', &
                            seismic_stiffness(i)%c14)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c15', &
                            seismic_stiffness(i)%c15)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c16', &
                            seismic_stiffness(i)%c16)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c22', &
                            seismic_stiffness(i)%c22)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c23', &
                            seismic_stiffness(i)%c23)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c24', &
                            seismic_stiffness(i)%c24)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c25', &
                            seismic_stiffness(i)%c25)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c26', &
                            seismic_stiffness(i)%c26)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c33', &
                            seismic_stiffness(i)%c33)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c34', &
                            seismic_stiffness(i)%c34)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c35', &
                            seismic_stiffness(i)%c35)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c36', &
                            seismic_stiffness(i)%c36)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c44', &
                            seismic_stiffness(i)%c44)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c45', &
                            seismic_stiffness(i)%c45)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c46', &
                            seismic_stiffness(i)%c46)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c55', &
                            seismic_stiffness(i)%c55)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c56', &
                            seismic_stiffness(i)%c56)
            call json%get('Seismic/Stiffness_Coefficients('//trim(adjustl(i_str))//')/c66', &
                            seismic_stiffness(i)%c66) 
            ! EM      
            call json%get('Electromagnetic/Permittivity_Coefficients('//trim(adjustl(i_str))//')/e11',&
                            electromagnetic_permittivity(i)%e11)         
            call json%get('Electromagnetic/Permittivity_Coefficients('//trim(adjustl(i_str))//')/e12',&
                            electromagnetic_permittivity(i)%e12)         
            call json%get('Electromagnetic/Permittivity_Coefficients('//trim(adjustl(i_str))//')/e13',&
                            electromagnetic_permittivity(i)%e13)    
            call json%get('Electromagnetic/Permittivity_Coefficients('//trim(adjustl(i_str))//')/e22',&
                            electromagnetic_permittivity(i)%e22)   
            call json%get('Electromagnetic/Permittivity_Coefficients('//trim(adjustl(i_str))//')/e23',&
                            electromagnetic_permittivity(i)%e23)   
            call json%get('Electromagnetic/Permittivity_Coefficients('//trim(adjustl(i_str))//')/e33',&
                            electromagnetic_permittivity(i)%e33)   
            call json%get('Electromagnetic/Conductivity_Coefficients('//trim(adjustl(i_str))//')/s11',&
                            electromagnetic_conductivity(i)%s11)   
            call json%get('Electromagnetic/Conductivity_Coefficients('//trim(adjustl(i_str))//')/s12',&
                            electromagnetic_conductivity(i)%s12)   
            call json%get('Electromagnetic/Conductivity_Coefficients('//trim(adjustl(i_str))//')/s13',&
                            electromagnetic_conductivity(i)%s13)
            call json%get('Electromagnetic/Conductivity_Coefficients('//trim(adjustl(i_str))//')/s22',&
                            electromagnetic_conductivity(i)%s22)   
            call json%get('Electromagnetic/Conductivity_Coefficients('//trim(adjustl(i_str))//')/s23',&
                            electromagnetic_conductivity(i)%s23)   
            call json%get('Electromagnetic/Conductivity_Coefficients('//trim(adjustl(i_str))//')/s33',&
                            electromagnetic_conductivity(i)%s33)   
            
        end do

        call json%destroy()
    end subroutine parse_json
    
    ! --------------------------------------------------------------------------
    subroutine read_geometry(file_name, domain, geometry)
    
        implicit none 
        
        character(len=*), intent(in) :: file_name 
        type(Domain_Type), intent(in) :: domain 
        integer, allocatable, intent(out) :: geometry(:,:)
        
        integer :: i,j 
        integer :: nx, nz 
        integer :: io_status 
        integer :: unit_number 
        
        ! Extract dimensions from the domain structure
        nx = domain%nx 
        nz = domain%nz 
        
        ! Allocate the 2D array based on nx and nz
        allocate(geometry(nx, nz))
        
        ! Open the geometry.dat file for reading
        open(newunit=unit_number, file=file_name, status = 'old', action = 'read', &
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
        type(Stiffness_Type), intent(in) :: stiffness(:)
        type(Attenuation_Type), intent(in) :: attenuation(:)
        
        real(real64), allocatable :: extended_stiffness(:,:,:)
        real(real64), allocatable :: extended_attenuation(:,:,:)
        integer :: nx, nz, cpml 
        integer :: i, j, k, id_value 
        integer :: nx_ext, nz_ext 
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
                extended_stiffness(i + cpml, j + cpml, 1) = stiffness(id_value)%c11
                extended_stiffness(i + cpml, j + cpml, 2) = stiffness(id_value)%c12
                extended_stiffness(i + cpml, j + cpml, 3) = stiffness(id_value)%c13
                extended_stiffness(i + cpml, j + cpml, 4) = stiffness(id_value)%c14
                extended_stiffness(i + cpml, j + cpml, 5) = stiffness(id_value)%c15
                extended_stiffness(i + cpml, j + cpml, 6) = stiffness(id_value)%c16
                extended_stiffness(i + cpml, j + cpml, 7) = stiffness(id_value)%c22
                extended_stiffness(i + cpml, j + cpml, 8) = stiffness(id_value)%c23
                extended_stiffness(i + cpml, j + cpml, 9) = stiffness(id_value)%c24
                extended_stiffness(i + cpml, j + cpml, 10) = stiffness(id_value)%c25
                extended_stiffness(i + cpml, j + cpml, 11) = stiffness(id_value)%c26
                extended_stiffness(i + cpml, j + cpml, 12) = stiffness(id_value)%c33
                extended_stiffness(i + cpml, j + cpml, 13) = stiffness(id_value)%c34
                extended_stiffness(i + cpml, j + cpml, 14) = stiffness(id_value)%c35
                extended_stiffness(i + cpml, j + cpml, 15) = stiffness(id_value)%c36
                extended_stiffness(i + cpml, j + cpml, 16) = stiffness(id_value)%c44
                extended_stiffness(i + cpml, j + cpml, 17) = stiffness(id_value)%c45
                extended_stiffness(i + cpml, j + cpml, 18) = stiffness(id_value)%c46
                extended_stiffness(i + cpml, j + cpml, 19) = stiffness(id_value)%c55
                extended_stiffness(i + cpml, j + cpml, 20) = stiffness(id_value)%c56
                extended_stiffness(i + cpml, j + cpml, 21) = stiffness(id_value)%c66
                ! same for attenuation
                extended_attenuation(i + cpml, j + cpml, 1) = attenuation(id_value)%alpha_x
                extended_attenuation(i + cpml, j + cpml, 2) = attenuation(id_value)%alpha_y
                extended_attenuation(i + cpml, j + cpml, 3) = attenuation(id_value)%alpha_z
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
        type(Permittivity_Type), intent(in) :: permittivity(:) 
        type(Conductivity_Type), intent(in) :: conductivity(:) 
        
        real(real64), allocatable :: extended_permittivity(:,:,:)
        real(real64), allocatable :: extended_conductivity(:,:,:)
        integer :: nx, nz, cpml 
        integer :: i, j, k, id_value 
        integer nx_ext, nz_ext 
        integer :: unit_number, io_status 
        
        character(len=9) :: permittivity_fn, conductivity_fn 
        character(len=5), dimension(6) :: pcoef_names = (/ &
                                'eps11','eps12','eps13','eps22','eps23','eps33'/)
        character(len=5), dimension(6) :: scoef_names = (/ &
                                'sig11','sig12','sig13','sig22','sig23','sig33'/)
        
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
                
                extended_permittivity(i + cpml, j + cpml, 1) = permittivity(id_value)%e11
                extended_permittivity(i + cpml, j + cpml, 2) = permittivity(id_value)%e12
                extended_permittivity(i + cpml, j + cpml, 3) = permittivity(id_value)%e13
                extended_permittivity(i + cpml, j + cpml, 4) = permittivity(id_value)%e22
                extended_permittivity(i + cpml, j + cpml, 5) = permittivity(id_value)%e23
                extended_permittivity(i + cpml, j + cpml, 6) = permittivity(id_value)%e33

                extended_conductivity(i + cpml, j + cpml, 1) = conductivity(id_value)%s11
                extended_conductivity(i + cpml, j + cpml, 2) = conductivity(id_value)%s12
                extended_conductivity(i + cpml, j + cpml, 3) = conductivity(id_value)%s13
                extended_conductivity(i + cpml, j + cpml, 4) = conductivity(id_value)%s22
                extended_conductivity(i + cpml, j + cpml, 5) = conductivity(id_value)%s23
                extended_conductivity(i + cpml, j + cpml, 6) = conductivity(id_value)%s33
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
    subroutine loadsource(filename, time_steps, srcfn)
        
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
        real(real64), dimension(..), intent(inout) :: image_data
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
        select rank(image_data)
        rank(2)
            if (readfile) then
                read(unit_number) image_data
            else
                write(unit_number) image_data
            end if
        rank(3)
            if (readfile) then
                read(unit_number) image_data
            else
                write(unit_number) image_data
            end if
        end select
        
        ! Close the file
        close(unit_number)

    end subroutine material_rw
    
    ! --------------------------------------------------------------------------
    subroutine write_image(image_data, domain, source, it, channel, SINGLE)

        implicit none

        ! Arguments
        real(real64), intent(in) :: image_data(..)
        type(Domain_Type), intent(in) :: domain 
        type(Source_Type), intent(in) :: source
        integer, intent(in) :: it
        character(len=2), intent(in) :: channel
        logical, intent(in) :: SINGLE

        ! Local variables
        character(len=100) :: filename
        integer :: unit_number
        real(real32), allocatable :: img_single_2d(:,:)
        real(real32), allocatable :: img_single_3d(:,:,:)


        ! Write the filename based on the rank of the src array
        if (domain%dim == 2) then
            WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a1, a4)" ) &
                    channel, '.', it, '.', source%xind, '.', source%zind, '.dat'
        else if (domain%dim == 2.5) then
            WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a1, i0, a4)" ) &
                    channel, '.', it, '.', source%xind, '.', source%yind, '.', source%zind, '.dat'
        else
            print *, "Error: src array must have 2 or 3 elements."
            stop
        end if

        ! Open the file with a dynamic unit number
        open(newunit=unit_number, form="unformatted", file=trim(filename))
        
        ! Handle image_data conversion and writing based on its rank
        
        ! Write based on the SINGLE flag and rank of the image_data
        if (SINGLE) then
            ! Handle the single precision case
            select rank(image_data)
            rank (2)
                ! Allocate a 2D array for single precision
                allocate(img_single_2d(domain%nx, domain%nz))
                ! Convert image_data from double to single precision
                img_single_2d = real(image_data, kind=real32)
                write(unit_number) img_single_2d
            rank (3)
                ! Allocate a 3D array for single precision
                allocate(img_single_3d(domain%nx, domain%ny, domain%nz))
                ! Convert image_data from double to single precision
                img_single_3d = real(image_data, kind=real32)
                write(unit_number) img_single_3d

            end select
        else
            ! Handle the double precision case
            select rank(image_data)
            rank(2)
                write(unit_number) image_data
            rank(3)
                write(unit_number) image_data

            end select
        endif

        ! Close the file
        close(unit_number)

    end subroutine write_image
    
    ! ! --------------------------------------------------------------------------
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
    
    