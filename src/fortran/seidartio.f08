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
    
    subroutine parse_json(file_name, domain, seismic_source, electromagnetic_source)
    
        implicit none 
        character(len=*), intent(in) :: file_name
        type(Domain_Type), intent(out) :: domain
        type(Source_Type), intent(out) :: seismic_source
        type(Source_Type), intent(out) :: electromagnetic_source

        type(json_file) :: json

        ! --------------------- End Declarations -------------------------------
        call json%initialize()
        call json%load_file(trim(file_name))

        ! Parse the domain data
        call json%get('Domain.dim', domain%dim)
        call json%get('Domain.nx', domain%nx)
        call json%get('Domain.ny', domain%ny)
        call json%get('Domain.nz', domain%nz)
        call json%get('Domain.dx', domain%dx)
        call json%get('Domain.dy', domain%dy)
        call json%get('Domain.dz', domain%dz)
        call json%get('Domain.cpml', domain%cpml)
        call json%get('Domain.nmats', domain%nmats)
        call json%get('Domain.image_file', domain%image_file)
    
        ! Parse the Seismic Source
        call json%get('Seismic.Source.x', seismic_source%x)
        call json%get('Seismic.Source.y', seismic_source%y)
        call json%get('Seismic.Source.z', seismic_source%z)
        call json%get('Seismic.Source.xind', seismic_source%xind)
        call json%get('Seismic.Source.yind', seismic_source%yind)
        call json%get('Seismic.Source.zind', seismic_source%zind)
        call json%get('Seismic.Source.source_frequency', seismic_source%source_frequency)
        call json%get('Seismic.Source.x-z_rotation', seismic_source%x_z_rotation)
        call json%get('Seismic.Source.x-y_rotation', seismic_source%x_y_rotation)
        call json%get('Seismic.Source.amplitude', seismic_source%amplitude)
        call json%get('Seismic.Source.source_type', seismic_source%source_type)
        call json%get('Seismic.Source.dt', seismic_source%dt)
        call json%get('Seismic.Source.time_steps', seismic_source%time_steps)
        
        ! Parse the EM Source
        call json%get('Electromagnetic.Source.x', electromagnetic_source%x)
        call json%get('Electromagnetic.Source.y', electromagnetic_source%y)
        call json%get('Electromagnetic.Source.z', electromagnetic_source%z)
        call json%get('Electromagnetic.Source.xind', electromagnetic_source%xind)
        call json%get('Electromagnetic.Source.yind', electromagnetic_source%yind)
        call json%get('Electromagnetic.Source.zind', electromagnetic_source%zind)
        call json%get('Electromagnetic.Source.source_frequency', electromagnetic_source%source_frequency)
        call json%get('Electromagnetic.Source.x-z_rotation', electromagnetic_source%x_z_rotation)
        call json%get('Electromagnetic.Source.x-y_rotation', electromagnetic_source%x_y_rotation)
        call json%get('Electromagnetic.Source.amplitude', electromagnetic_source%amplitude)
        call json%get('Electromagnetic.Source.source_type', electromagnetic_source%source_type)
        call json%get('Electromagnetic.Source.dt', electromagnetic_source%dt)
        call json%get('Electromagnetic.Source.time_steps', electromagnetic_source%time_steps)

        call json%destroy()
    end subroutine parse_json
    
    ! --------------------------------------------------------------------------
    subroutine read_geometry(file_name, domain, geometry)
    
        implicit none 
        
        character(len=*), intent(in) :: file_name 
        type(Domain_Type), intent(in) :: domain 
        integer, allocatable, intent(out) :: geometry(:,:)
        
        integer :: nx, nz 
        integer :: io_status 
        integer :: unit_number 
        
        ! Extract dimensions from the domain structure
        nx = domain%nx 
        nz = domain%nz 
        
        ! Allocate the 2D array based on nx and nz
        allocate(geometry(nx, nz))
        
        ! Open the geometry.dat file for reading
        open(newunit=unit_number, file=trim(file_name), form='unformatted', &
                status='old', iostat=io_status)
        
        if (io_status /= 0) then 
            print *, "Error opening file: ", file_name 
            stop 
        end if
        
        ! Read the nx-by-nz array from the file
        read(unit_number) geometry
        if (io_status /= 0) then
            print *, "Error reading geometry data"
            stop
        end if
        
        ! Close the file 
        close(unit_number)
    end subroutine read_geometry

    ! --------------------------------------------------------------------------
    ! subroutine seismic_parameter_write(domain, geometry, stiffness, attenuation)
    !     implicit none 
        
    !     type(Domain_Type), intent(in) :: domain 
    !     integer, intent(in) :: geometry(:,:)
    !     type(Stiffness_Type), intent(in) :: stiffness(:)
    !     type(Attenuation_Type), intent(in) :: attenuation(:)
        
    !     real(real64), allocatable :: extended_stiffness(:,:,:)
    !     real(real64), allocatable :: extended_attenuation(:,:,:)
    !     real(real64), allocatable :: density_gradient(:,:)
        
    !     integer :: nx, nz, cpml 
    !     integer :: i, j, k, id_value 
    !     integer :: nx_ext, nz_ext 
    !     integer :: unit_number, io_status 
        
    !     character(len=7) :: stiffness_fn
    !     ! character(len=10) :: attenuation_fn
    !     character(len=3), dimension(22) :: scoef_names = (/ &
    !                                         'c11', 'c12', 'c13', 'c14', 'c15', &
    !                                         'c16', 'c22', 'c23', 'c24', 'c25', &
    !                                         'c26', 'c33', 'c34', 'c35', 'c36', &
    !                                         'c44', 'c45', 'c46', 'c55', 'c56', &
    !                                         'c66', 'rho' /)
        
    !     character(len=10), dimension(3) :: attenuation_fn = (/ &
    !                                 'gammax.dat', 'gammay.dat', 'gammaz.dat'/)
        
    !     ! Extract domain dimensions and cpml
    !     nx = domain%nx
    !     nz = domain%nz
    !     cpml = domain%cpml

    !     ! Calculate extended dimensions
    !     nx_ext = 2 * cpml + nx
    !     nz_ext = 2 * cpml + nz
        
    !     ! Allocate the extended stiffness array with extra CPML padding
    !     allocate(extended_stiffness(nx_ext, nz_ext, 22))
    !     allocate(extended_attenuation(nx_ext, nz_ext, 3))  
    !     allocate(density_gradient(nx, nz))
        
    !     ! Load the density gradient 
    !     call material_rw2('density_gradient.dat', density_gradient, .TRUE. )
        
    !     extended_stiffness(:,:,:) = 0.0 
    !     extended_attenuation(:,:,:) = 0.0

    !      ! Loop over the geometry array and assign stiffness coefficients
    !     do j = 1, nz
    !         do i = 1, nx
    !             id_value = geometry(i, j) + 1  ! Get the material id from geometry

    !             ! Assign stiffness coefficients based on id_value from geometry
    !             extended_stiffness(i + cpml, j + cpml, 1) = stiffness(id_value)%c11
    !             extended_stiffness(i + cpml, j + cpml, 2) = stiffness(id_value)%c12
    !             extended_stiffness(i + cpml, j + cpml, 3) = stiffness(id_value)%c13
    !             extended_stiffness(i + cpml, j + cpml, 4) = stiffness(id_value)%c14
    !             extended_stiffness(i + cpml, j + cpml, 5) = stiffness(id_value)%c15
    !             extended_stiffness(i + cpml, j + cpml, 6) = stiffness(id_value)%c16
    !             extended_stiffness(i + cpml, j + cpml, 7) = stiffness(id_value)%c22
    !             extended_stiffness(i + cpml, j + cpml, 8) = stiffness(id_value)%c23
    !             extended_stiffness(i + cpml, j + cpml, 9) = stiffness(id_value)%c24
    !             extended_stiffness(i + cpml, j + cpml, 10) = stiffness(id_value)%c25
    !             extended_stiffness(i + cpml, j + cpml, 11) = stiffness(id_value)%c26
    !             extended_stiffness(i + cpml, j + cpml, 12) = stiffness(id_value)%c33
    !             extended_stiffness(i + cpml, j + cpml, 13) = stiffness(id_value)%c34
    !             extended_stiffness(i + cpml, j + cpml, 14) = stiffness(id_value)%c35
    !             extended_stiffness(i + cpml, j + cpml, 15) = stiffness(id_value)%c36
    !             extended_stiffness(i + cpml, j + cpml, 16) = stiffness(id_value)%c44
    !             extended_stiffness(i + cpml, j + cpml, 17) = stiffness(id_value)%c45
    !             extended_stiffness(i + cpml, j + cpml, 18) = stiffness(id_value)%c46
    !             extended_stiffness(i + cpml, j + cpml, 19) = stiffness(id_value)%c55
    !             extended_stiffness(i + cpml, j + cpml, 20) = stiffness(id_value)%c56
    !             extended_stiffness(i + cpml, j + cpml, 21) = stiffness(id_value)%c66
    !             extended_stiffness(i + cpml, j + cpml, 22) = stiffness(id_value)%density
    !             ! same for attenuation
    !             extended_attenuation(i + cpml, j + cpml, 1) = attenuation(id_value)%alpha_x
    !             extended_attenuation(i + cpml, j + cpml, 2) = attenuation(id_value)%alpha_y
    !             extended_attenuation(i + cpml, j + cpml, 3) = attenuation(id_value)%alpha_z
    !         end do
    !     end do
        
    !     ! Scale the density by the gradient that was calculated around air boundaries
    !     extended_stiffness(cpml+1:nx+cpml,cpml+1:nz+cpml,22) = extended_stiffness(cpml+1:nx+cpml,cpml+1:nz+cpml,22) * density_gradient
        
    !     ! Extend the values at 1 + cpml to the 1 to cpml indices
    !     do k = 1, 21
    !         ! Extend along x-direction (rows)
    !         do j = 1, nz_ext
    !             do i = 1, cpml
    !                 extended_stiffness(i, j, k) = extended_stiffness(1 + cpml, j, k)
    !                 extended_stiffness(nx_ext - i + 1, j, k) = extended_stiffness(nx_ext - cpml, j, k)
    !             end do
    !         end do
            
    !         ! Extend along z-direction (columns)
    !         do i = 1, nx_ext
    !             do j = 1, cpml
    !                 extended_stiffness(i, j, k) = extended_stiffness(i, 1 + cpml, k)
    !                 extended_stiffness(i, nz_ext - j + 1, k) = extended_stiffness(i, nz_ext - cpml, k)
    !             end do
    !         end do
    !     end do
        
    !     do k = 1, 3
    !         ! Extend along x-direction (rows)
    !         do j = 1, nz_ext
    !             do i = 1, cpml
    !                 extended_attenuation(i, j, k) = extended_attenuation(1 + cpml, j, k)
    !                 extended_attenuation(nx_ext - i + 1, j, k) = extended_attenuation(nx_ext - cpml, j, k)
    !             end do
    !         end do
            
    !         ! Extend along z-direction (columns)
    !         do i = 1, nx_ext
    !             do j = 1, cpml
    !                 extended_attenuation(i, j, k) = extended_attenuation(i, 1 + cpml, k)
    !                 extended_attenuation(i, nz_ext - j + 1, k) = extended_attenuation(i, nz_ext - cpml, k)
    !             end do
    !         end do
    !     end do
        
    !     ! Write each 2D array (stiffness coefficients) to separate files
    !     do k = 1, 22
    !         ! Create filename based on coefficient name (e.g., 'c11.dat')
    !         write(stiffness_fn, '(A3, ".dat")') scoef_names(k)
    !         call material_rw2(stiffness_fn, extended_stiffness(:,:,k), .FALSE.)
    !     end do
        
    !     do k = 1, 3
    !         ! Create filename based on coefficient name (e.g., 'c11.dat')
    !         ! write(attenuation_fn, '(A3, ".dat")') acoef_names(k)

    !         ! Open the file to write the 2D array
    !         open(newunit=unit_number, file=attenuation_fn(k), form='unformatted', &
    !             access='stream', status='replace', iostat=io_status)
    !         if (io_status /= 0) then
    !             print *, "Error opening file for writing: ", attenuation_fn(k)
    !             stop
    !         end if

    !         ! Write the 2D array to the binary file
    !         write(unit_number) extended_attenuation(:,:,k)

    !         ! Close the file
    !         close(unit_number)
    !     end do
        
    !     deallocate(extended_stiffness)
    !     deallocate(extended_attenuation)
    ! end subroutine seismic_parameter_write
    
    ! ! --------------------------------------------------------------------------
    ! subroutine electromagnetic_parameter_write(domain, geometry, permittivity, conductivity)
    !     implicit none
        
    !     type(Domain_Type), intent(in) :: domain 
    !     integer, intent(in) :: geometry(:,:) 
    !     type(Permittivity_Type), intent(in) :: permittivity(:) 
    !     type(Conductivity_Type), intent(in) :: conductivity(:) 
        
    !     real(real64), allocatable :: extended_permittivity(:,:,:)
    !     real(real64), allocatable :: extended_conductivity(:,:,:)
    !     integer :: nx, nz, cpml 
    !     integer :: i, j, k, id_value 
    !     integer nx_ext, nz_ext 
    !     integer :: unit_number, io_status 
        
    !     character(len=9) :: permittivity_fn, conductivity_fn 
    !     character(len=5), dimension(6) :: pcoef_names = (/ &
    !                             'eps11','eps12','eps13','eps22','eps23','eps33'/)
    !     character(len=5), dimension(6) :: scoef_names = (/ &
    !                             'sig11','sig12','sig13','sig22','sig23','sig33'/)
        
    !     ! Extract domain dimensions and cpml
    !     nx = domain%nx
    !     nz = domain%nz
    !     cpml = domain%cpml

    !     ! Calculate extended dimensions
    !     nx_ext = 2 * cpml + nx
    !     nz_ext = 2 * cpml + nz
        
    !     ! Allocate the extended permittivity array with extra CPML padding
    !     allocate(extended_permittivity(nx_ext, nz_ext, 21))
    !     allocate(extended_conductivity(nx_ext, nz_ext, 3))  
        
    !     extended_permittivity(:,:,:) = 0.0 
    !     extended_conductivity(:,:,:) = 0.0
        
    !     do j = 1,nz 
    !         do i = 1,nx 
    !             id_value = geometry(i,j) 
                
    !             extended_permittivity(i + cpml, j + cpml, 1) = permittivity(id_value)%e11
    !             extended_permittivity(i + cpml, j + cpml, 2) = permittivity(id_value)%e12
    !             extended_permittivity(i + cpml, j + cpml, 3) = permittivity(id_value)%e13
    !             extended_permittivity(i + cpml, j + cpml, 4) = permittivity(id_value)%e22
    !             extended_permittivity(i + cpml, j + cpml, 5) = permittivity(id_value)%e23
    !             extended_permittivity(i + cpml, j + cpml, 6) = permittivity(id_value)%e33

    !             extended_conductivity(i + cpml, j + cpml, 1) = conductivity(id_value)%s11
    !             extended_conductivity(i + cpml, j + cpml, 2) = conductivity(id_value)%s12
    !             extended_conductivity(i + cpml, j + cpml, 3) = conductivity(id_value)%s13
    !             extended_conductivity(i + cpml, j + cpml, 4) = conductivity(id_value)%s22
    !             extended_conductivity(i + cpml, j + cpml, 5) = conductivity(id_value)%s23
    !             extended_conductivity(i + cpml, j + cpml, 6) = conductivity(id_value)%s33
    !         end do 
    !     end do 
        
    !     ! Extend the values at 1 + cpml to the 1 to cpml indices
    !     do k = 1, 6
    !       ! Extend along x-direction (rows)
    !       do j = 1, nz_ext
    !         do i = 1, cpml
    !           extended_permittivity(i, j, k) = extended_permittivity(1 + cpml, j, k)
    !           extended_conductivity(i, j, k) = extended_conductivity(1 + cpml, j, k)
    !           extended_permittivity(nx_ext - i + 1, j, k) = extended_permittivity(nx_ext - cpml, j, k)
    !           extended_conductivity(nx_ext - i + 1, j, k) = extended_conductivity(nx_ext - cpml, j, k)
    !         end do
    !       end do
        
    !       ! Extend along z-direction (columns)
    !       do i = 1, nx_ext
    !         do j = 1, cpml
    !           extended_permittivity(i, j, k) = extended_permittivity(i, 1 + cpml, k)
    !           extended_conductivity(i, j, k) = extended_conductivity(i, 1 + cpml, k)
    !           extended_permittivity(i, nz_ext - j + 1, k) = extended_permittivity(i, nz_ext - cpml, k)
    !           extended_conductivity(i, nz_ext - j + 1, k) = extended_conductivity(i, nz_ext - cpml, k)
    !         end do
    !       end do
    !     end do
        
    !     do k = 1, 6
    !         ! Create filename based on coefficient name (e.g., 'c11.dat')
    !         write(permittivity_fn, '(A5, ".dat")') pcoef_names(k)

    !         ! Open the file to write the 2D array
    !         open(newunit=unit_number, file=permittivity_fn, form='unformatted', &
    !             access='stream', status='replace', iostat=io_status)
    !         if (io_status /= 0) then
    !             print *, "Error opening file for writing: ", permittivity_fn
    !             stop
    !         end if

    !         ! Write the 2D array to the binary file
    !         write(unit_number) extended_permittivity(:,:,k)

    !         ! Close the file
    !         close(unit_number)
    !     end do
    !      do k = 1, 3
    !         ! Create filename based on coefficient name (e.g., 'c11.dat')
    !         write(conductivity_fn, '(A3, ".dat")') scoef_names(k)

    !         ! Open the file to write the 2D array
    !         open(newunit=unit_number, file=conductivity_fn, form='unformatted', &
    !             access='stream', status='replace', iostat=io_status)
    !         if (io_status /= 0) then
    !             print *, "Error opening file for writing: ", conductivity_fn
    !             stop
    !         end if

    !         ! Write the 2D array to the binary file
    !         write(unit_number) extended_conductivity(:,:,k)

    !         ! Close the file
    !         close(unit_number)
    !     end do
        
    !     deallocate(extended_permittivity)
    !     deallocate(extended_conductivity)
        
    ! end subroutine electromagnetic_parameter_write
    
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
        
        integer :: unit_number
        character(len=*) :: filename
        real(real64),dimension(:,:) :: image_data
        logical :: readfile
        
        open(newunit=unit_number, form="unformatted", file = trim(filename))
        
        if ( readfile ) then
            read(unit_number) image_data
        else
            write(unit_number) image_data
        endif
        
        close(unit = unit_number)

    end subroutine material_rw
    
    ! --------------------------------------------------------------------------
    subroutine material_rw2(filename, image_data, readfile)

        implicit none
        
        character(len=*) :: filename
        real(real64),dimension(:,:) :: image_data
        logical :: readfile
        integer :: unit_number
        
        
        open(newunit = unit_number, form="unformatted", file = trim(filename) )
        
        if ( readfile ) then
            read(unit_number) image_data
        else
            write(unit_number) image_data
        endif
        
        close(unit_number)

    end subroutine material_rw2
    
    ! --------------------------------------------------------------------------
    subroutine material_rw3(filename, image_data, readfile)

        implicit none
        
        character(len=*) :: filename
        real(real64),dimension(:,:,:) :: image_data
        logical :: readfile
        integer :: unit_number
        
        open(newunit = unit_number, form="unformatted", file = trim(filename))
        
        if ( readfile ) then
            read(unit_number) image_data
        else
            write(unit_number) image_data
        endif
        
        close(unit_number)

    end subroutine material_rw3
    
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
            WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a4)" ) &
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
    
    ! --------------------------------------------------------------------------
    subroutine write_image2(image_data, nx, nz, source, it, channel, SINGLE)
    
        implicit none

        integer :: unit_number, io_status
        integer :: nx, nz, it
        type(Source_Type), intent(in) :: source
        real(real64), intent(in) :: image_data(nx, nz)
        character(len=2) :: channel
        character(len=100) :: filename
        logical :: SINGLE

        ! WRITE (filename, "(a2, i6.6, '.dat')" ) channel, it
        WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a4)" ) &
            channel, '.', it, '.', source%xind, '.', source%zind, '.dat'
        open(newunit = unit_number, form = 'unformatted', action="write", &
                file = trim(filename), access="sequential", iostat=io_status)

        if (io_status /= 0) then
            print *, "Error opening file"
            stop
        end if
        
        if (SINGLE) then
            write(unit_number) sngl(image_data)
        else
            write(unit_number) image_data 
        end if 
        
        flush(unit_number)
        close(unit_number)

    end subroutine write_image2
    
    ! subroutine write_image2(image_data, nx, nz, src, it, channel, SINGLE)
    
    !     implicit none

    !     integer, parameter :: dp = kind(0.d0)
    !     integer :: unit_number
    !     integer :: nx, nz, it
    !     integer,dimension(2) :: src
    !     real(kind=dp) :: image_data(nx, nz)
    !     character(len=2) :: channel
    !     character(len=100) :: filename
    !     logical :: SINGLE

    !     ! WRITE (filename, "(a2, i6.6, '.dat')" ) channel, it
    !     WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a1, a4)" ) &
    !                 channel,'.', it,'.', src(1),'.', src(2), '.','.dat'
    !     open(newunit = unit_number, form = 'unformatted', file = trim(filename) )

    !     if (SINGLE) then
    !         write(unit_number) sngl(image_data)
    !     else
    !         write(unit_number) image_data 
    !     end if 


    !     close(unit = unit_number)

    ! end subroutine write_image2
    
    ! --------------------------------------------------------------------------
    subroutine write_image3(image_data, nx, ny, nz, source, it, channel, SINGLE)
    
        implicit none
    
        integer :: nx, ny, nz, it, unit_number
        type(Source_Type), intent(in) :: source
        real(real64), intent(in) :: image_data(nx, ny, nz)
        character(len=2) :: channel
        character(len=80) :: filename
        logical :: SINGLE
        
        WRITE (filename, "(a2, a1, i6.6, a1, i0, a1, i0, a1, i0, a4)" ) &
            channel, '.', it, '.', source%xind, '.', source%yind, '.', source%zind, '.dat'
        
        open(newunit = unit_number, form = 'unformatted', file = trim(filename) )
        
        if (SINGLE) then
            write(unit_number) sngl(image_data)
        else
            write(unit_number) image_data 
        end if 
        
        close(unit = unit_number)

    end subroutine write_image3
    
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
    
    