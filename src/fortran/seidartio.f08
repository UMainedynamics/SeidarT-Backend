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
        character(len=:), allocatable :: temp_string
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
        ! call json%get('Domain.image_file', domain%image_file)
        call json%get('Domain.image_file', temp_string)
        domain%image_file = trim(adjustl(temp_string))
        
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
        call json%get('Seismic.Source.y-z_rotation', seismic_source%y_z_rotation)
        call json%get('Seismic.Source.amplitude', seismic_source%amplitude)
        ! call json%get('Seismic.Source.source_type', seismic_source%source_type)
        call json%get('Seismic.Source.dt', seismic_source%dt)
        call json%get('Seismic.Source.time_steps', seismic_source%time_steps)
        
        call json%get('Seismic.Source.source_type', temp_string)
        seismic_source%source_type = trim(adjustl(temp_string))
        call json%get('Seismic.Source.source_wavelet', temp_string)
        seismic_source%source_wavelet = trim(adjustl(temp_string))
        
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
        call json%get('Electromagnetic.Source.y-z_rotation', electromagnetic_source%y_z_rotation)
        call json%get('Electromagnetic.Source.amplitude', electromagnetic_source%amplitude)
        ! call json%get('Electromagnetic.Source.source_type', electromagnetic_source%source_type)
        call json%get('Electromagnetic.Source.dt', electromagnetic_source%dt)
        call json%get('Electromagnetic.Source.time_steps', electromagnetic_source%time_steps)
        
        call json%get('Electromagnetic.Source.source_type', temp_string)
        electromagnetic_source%source_type = trim(adjustl(temp_string))
        
        call json%destroy()
    
    end subroutine parse_json
    
    ! --------------------------------------------------------------------------
    subroutine read_geometry(file_name, geometry)
    
        implicit none 
        
        character(len=*), intent(in) :: file_name 
        ! type(Domain_Type), intent(in) :: domain 
        integer, dimension(:,:) :: geometry(:,:)
        
        ! integer :: nx, nz 
        integer :: io_status 
        integer :: unit_number 
        
        ! ! Extract dimensions from the domain structure
        ! nx = domain%nx 
        ! nz = domain%nz 
        
        ! ! Allocate the 2D array based on nx and nz
        ! allocate(geometry(nx, nz))
        
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
    
    