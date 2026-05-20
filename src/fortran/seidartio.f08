module seidartio 

    use iso_fortran_env, only: real64, real32
    use, intrinsic :: iso_c_binding
    use json_module 
    use seidart_types 
    
    implicit none 
    
    public :: parse_json
    
    real(c_float), allocatable, save, target :: block_buffer(:,:,:,:,:)
    integer, save :: total_block_limit = 0
    integer, save :: current_step_in_block = 0
    integer, save :: output_cpml = 0

    interface
        integer(c_size_t) function ZSTD_compressBound(srcSize) bind(C, name="ZSTD_compressBound")
            import :: c_size_t
            integer(c_size_t), value :: srcSize
        end function ZSTD_compressBound

        integer(c_size_t) function ZSTD_compress(dst, dstCap, src, srcSize, compressionLevel) &
            bind(C, name="ZSTD_compress")
            import :: c_ptr, c_size_t, c_int
            type(c_ptr), value :: dst
            integer(c_size_t), value :: dstCap
            type(c_ptr), value :: src
            integer(c_size_t), value :: srcSize
            integer(c_int), value :: compressionLevel
        end function ZSTD_compress

        integer(c_int) function ZSTD_isError(code) bind(C, name="ZSTD_isError")
            import :: c_size_t, c_int
            integer(c_size_t), value :: code
        end function ZSTD_isError
    end interface
    
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
        call json%get('Electromagnetic.Source.source_wavelet', temp_string)
        electromagnetic_source%source_wavelet = trim(adjustl(temp_string))

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
    subroutine write_array(filename, length, array)
        
        implicit none
        
        character(len=*) :: filename
        integer :: length
        real(real64),dimension(length) :: array
        integer :: unit_number
        
        open(newunit=unit_number, form="unformatted", file = trim(filename))
        
        write(unit_number) array
        close(unit_number)

    end subroutine write_array
    
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
            print *, "Reading file: ", trim(filename)
            read(unit_number) image_data
        else
            write(unit_number) image_data
        endif
        
        close(unit_number)

    end subroutine material_rw2
    
        ! --------------------------------------------------------------------------
    subroutine material_rw2c(filename, image_data, readfile)

        implicit none
        
        character(len=*) :: filename
        complex(real64),dimension(:,:) :: image_data
        logical :: readfile
        integer :: unit_number
        
        
        open(newunit = unit_number, form="unformatted", file = trim(filename) )
        
        if ( readfile ) then
            print *, "Reading file: ", trim(filename)
            read(unit_number) image_data
        else
            write(unit_number) image_data
        endif
        
        close(unit_number)

    end subroutine material_rw2c
    
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
    
    ! ==========================================================================
    ! =========================== ZSTD block writing ===========================
    subroutine init_io(nx, ny, nz, cpml, n_steps_per_block)
        integer, intent(in) :: nx, ny, nz, cpml, n_steps_per_block
        integer :: output_nx, output_ny, output_nz

        if (n_steps_per_block < 1) then
            print *, 'Error: steps per zstd block must be at least 1'
            stop
        end if
        if (cpml < 0) then
            print *, 'Error: cpml must be non-negative for zstd block output'
            stop
        end if

        output_nx = nx - 2 * cpml
        output_ny = ny - 2 * cpml
        output_nz = nz - 2 * cpml
        if (output_nx < 1 .or. output_ny < 1 .or. output_nz < 1) then
            print *, 'Error: cpml removes the entire zstd block output domain'
            stop
        end if

        if (allocated(block_buffer)) deallocate(block_buffer)
        output_cpml = cpml
        total_block_limit = n_steps_per_block 
        ! shape: (x, y, z, component, step) 
        allocate(block_buffer(output_nx, output_ny, output_nz, 3, total_block_limit) )
        current_step_in_block = 0 

    end subroutine init_io  
    
    ! --------------------------------------------------------------------------
    subroutine add_step_to_block(filename, Ex, Ey, Ez) 
        character(len=*), intent(in) :: filename
        real(c_double), intent(in) :: Ex(:,:,:)
        real(c_double), intent(in) :: Ey(:,:,:)
        real(c_double), intent(in) :: Ez(:,:,:)

        if (.not. allocated(block_buffer)) then
            print *, 'Error: zstd block I/O has not been initialized'
            stop
        end if
        
        current_step_in_block = current_step_in_block + 1
        
        block_buffer(:,:,:,1,current_step_in_block) = &
            real(Ex(output_cpml+1:size(Ex,1)-output_cpml, &
                    output_cpml+1:size(Ex,2)-output_cpml, &
                    output_cpml+1:size(Ex,3)-output_cpml), c_float)
        block_buffer(:,:,:,2,current_step_in_block) = &
            real(Ey(output_cpml+1:size(Ey,1)-output_cpml, &
                    output_cpml+1:size(Ey,2)-output_cpml, &
                    output_cpml+1:size(Ey,3)-output_cpml), c_float)
        block_buffer(:,:,:,3,current_step_in_block) = &
            real(Ez(output_cpml+1:size(Ez,1)-output_cpml, &
                    output_cpml+1:size(Ez,2)-output_cpml, &
                    output_cpml+1:size(Ez,3)-output_cpml), c_float)
        
        ! if the block is full, write it and reset 
        if (current_step_in_block == total_block_limit) then
            call write_zstd_block(filename)
        end if
    end subroutine add_step_to_block 
    
    ! --------------------------------------------------------------------------
    subroutine write_zstd_block(filename) 
        character(len=*), intent(in) :: filename
        
        integer(c_size_t) :: src_size, max_dst_size, compressed_size
        integer(c_size_t) :: n_values
        integer(c_int) :: compression_level 
        type(c_ptr) :: src_ptr, dst_ptr 
        
        ! For file I/O
        integer :: i_unit, io_status
        character(kind=c_char), allocatable, target :: comp_buffer(:)

        if (.not. allocated(block_buffer)) then
            print *, 'Error: zstd block I/O has not been initialized'
            stop
        end if

        if (current_step_in_block < 1) return

        ! Compress only the populated leading portion of the Fortran-order buffer.
        n_values = int(size(block_buffer, 1), c_size_t) * &
                   int(size(block_buffer, 2), c_size_t) * &
                   int(size(block_buffer, 3), c_size_t) * &
                   int(size(block_buffer, 4), c_size_t) * &
                   int(current_step_in_block, c_size_t)
        src_size = n_values * c_sizeof(block_buffer(1,1,1,1,1))
        max_dst_size = ZSTD_compressBound(src_size)
        allocate(comp_buffer(max_dst_size))
        
        src_ptr = c_loc(block_buffer(1,1,1,1,1))
        dst_ptr = c_loc(comp_buffer(1))
        
        ! Level 3 is a good balance. Use 1 for speed or up to 22 for max compression.
        compression_level = 3 
        compressed_size = ZSTD_compress(dst_ptr, max_dst_size, src_ptr, src_size, compression_level)
        
        if (ZSTD_isError(compressed_size) /= 0) then
            print *, 'Error: ZSTD compression failed for file: ', trim(filename)
            deallocate(comp_buffer)
            stop
        end if

        open(newunit=i_unit, file=trim(filename), access='stream', form='unformatted', &
             status='replace', action='write', iostat=io_status)
        if (io_status /= 0) then
            print *, 'Error opening zstd output file: ', trim(filename)
            deallocate(comp_buffer)
            stop
        end if
        write(i_unit) comp_buffer(1:compressed_size)
        close(i_unit)
        
        deallocate(comp_buffer) 
        print*, "Successfully wrote zstd block to file: ", trim(filename), &
                " steps=", current_step_in_block, &
                " raw_bytes=", src_size, &
                " compressed_bytes=", compressed_size
        current_step_in_block = 0
    
    end subroutine write_zstd_block

    ! --------------------------------------------------------------------------
    subroutine finalize_io(filename)
        character(len=*), intent(in) :: filename

        if (allocated(block_buffer)) then
            if (current_step_in_block > 0) call write_zstd_block(filename)
            deallocate(block_buffer)
        end if
        total_block_limit = 0
        current_step_in_block = 0
        output_cpml = 0

    end subroutine finalize_io

    ! ==========================================================================
    ! ======================== 2D ZSTD block writing ===========================
    ! These are 2D companions to the 3D routines above. The trick is that
    ! init_io_2d allocates block_buffer with y-dimension = 1, so
    ! write_zstd_block and finalize_io work unchanged.
    ! ==========================================================================

    ! --------------------------------------------------------------------------
    !> Initialize the block buffer for 2D field output.
    !!
    !! Allocates block_buffer as (output_nx, 1, output_nz, 3, steps_per_block).
    !! The singleton y-dimension makes write_zstd_block and finalize_io
    !! work identically for 2D and 3D cases.
    subroutine init_io_2d(nx, nz, cpml, n_steps_per_block)
        integer, intent(in) :: nx, nz, cpml, n_steps_per_block
        integer :: output_nx, output_nz

        if (n_steps_per_block < 1) then
            print *, 'Error: steps per zstd block must be at least 1'
            stop
        end if
        if (cpml < 0) then
            print *, 'Error: cpml must be non-negative for zstd block output'
            stop
        end if

        output_nx = nx - 2 * cpml
        output_nz = nz - 2 * cpml
        if (output_nx < 1 .or. output_nz < 1) then
            print *, 'Error: cpml removes the entire zstd block output domain'
            stop
        end if

        if (allocated(block_buffer)) deallocate(block_buffer)
        output_cpml = cpml
        total_block_limit = n_steps_per_block
        ! shape: (x, 1, z, component, step) — y-dim is 1 for 2D
        allocate(block_buffer(output_nx, 1, output_nz, 3, total_block_limit))
        current_step_in_block = 0

    end subroutine init_io_2d

    ! --------------------------------------------------------------------------
    !> Accumulate one time step of 2D field data into the block buffer.
    !!
    !! Accepts rank-2 field arrays Ex(nx,nz), Ey(nx,nz), Ez(nx,nz),
    !! strips the CPML padding, converts to real32, and stores in the
    !! block_buffer(:,1,:,component,step) slice.
    !! When the block is full, automatically calls write_zstd_block.
    subroutine add_step_to_block_2d(filename, Ex, Ey, Ez)
        character(len=*), intent(in) :: filename
        real(c_double), intent(in) :: Ex(:,:)
        real(c_double), intent(in) :: Ey(:,:)
        real(c_double), intent(in) :: Ez(:,:)

        if (.not. allocated(block_buffer)) then
            print *, 'Error: zstd block I/O has not been initialized'
            stop
        end if

        current_step_in_block = current_step_in_block + 1

        block_buffer(:,1,:,1,current_step_in_block) = &
            real(Ex(output_cpml+1:size(Ex,1)-output_cpml, &
                    output_cpml+1:size(Ex,2)-output_cpml), c_float)
        block_buffer(:,1,:,2,current_step_in_block) = &
            real(Ey(output_cpml+1:size(Ey,1)-output_cpml, &
                    output_cpml+1:size(Ey,2)-output_cpml), c_float)
        block_buffer(:,1,:,3,current_step_in_block) = &
            real(Ez(output_cpml+1:size(Ez,1)-output_cpml, &
                    output_cpml+1:size(Ez,2)-output_cpml), c_float)

        ! if the block is full, write it and reset
        if (current_step_in_block == total_block_limit) then
            call write_zstd_block(filename)
        end if
    end subroutine add_step_to_block_2d

    ! --------------------------------------------------------------------------
    !> Configure block I/O parameters for a 2D simulation.
    !!
    !! Reads FDTD_LEGACY_OUTPUT and FDTD_STEPS_PER_BLOCK from the environment,
    !! computes output dimensions after CPML stripping, and auto-sizes the
    !! block to target ~2 GB of buffer memory.
    subroutine setup_io_params_2d(nx, nz, cpml, block_output, steps_per_block, legacy_output)
        integer, intent(in) :: nx, nz, cpml
        logical, intent(out) :: block_output, legacy_output
        integer, intent(out) :: steps_per_block

        character(len=32) :: env_val
        integer :: stat
        real(8) :: bytes_per_step
        integer :: output_nx, output_nz

        ! Block output default is true, legacy is false
        block_output = .true.
        legacy_output = .false.
        ! Check for legacy override
        call get_environment_variable("FDTD_LEGACY_OUTPUT", env_val, status = stat)
        if (stat == 0 .and. trim(env_val) == "TRUE") then
            legacy_output = .true.
            block_output = .false.
        endif

        output_nx = nx - 2 * cpml
        output_nz = nz - 2 * cpml
        if (output_nx < 1 .or. output_nz < 1) then
            print *, 'Error: cpml removes the entire zstd block output domain'
            stop
        end if

        ! check for steps per block override
        call get_environment_variable("FDTD_STEPS_PER_BLOCK", env_val, status = stat)
        if (stat == 0) then
            read(env_val, *) steps_per_block
        else
            ! auto calculate with target ~2GB buffer
            ! bytes = interior cells * 3 components * 4 bytes per real32
            ! (no y dimension for 2D)
            bytes_per_step = real(output_nx, 8) * output_nz * 3 * 4
            ! steps = 2GB (2^31 bytes) / bytes_per_step
            steps_per_block = max(1, int(2147483648.0_8 / bytes_per_step))
            if (steps_per_block > 500) steps_per_block = 500
        endif

        print*, "---- I/O Configuration (2D) ----"
        print*, "Block output: ", block_output
        print*, "Legacy output: ", legacy_output
        if (block_output) then
            print*, "Block output dimensions: ", output_nx, output_nz
            print*, "Block output precision: real32"
            print*, "Steps per block: ", steps_per_block
        end if
    end subroutine setup_io_params_2d

    
    ! --------------------------------------------------------------------------
    subroutine setup_io_params(nx, ny, nz, cpml, block_output, steps_per_block, legacy_output)
        integer, intent(in) :: nx, ny, nz, cpml
        logical, intent(out) :: block_output, legacy_output
        integer, intent(out) :: steps_per_block
        
        character(len=32) :: env_val 
        integer :: stat 
        real(8) :: bytes_per_step 
        integer :: output_nx, output_ny, output_nz
        
        ! Block output defaul is true, legacy is false 
        block_output = .true. 
        legacy_output = .false.
        ! Check for legacy override 
        call get_environment_variable("FDTD_LEGACY_OUTPUT", env_val, status = stat) 
        if (stat == 0 .and. trim(env_val) == "TRUE") then
            legacy_output = .true.
            block_output = .false.
        endif

        output_nx = nx - 2 * cpml
        output_ny = ny - 2 * cpml
        output_nz = nz - 2 * cpml
        if (output_nx < 1 .or. output_ny < 1 .or. output_nz < 1) then
            print *, 'Error: cpml removes the entire zstd block output domain'
            stop
        end if
        
        ! check for steps per block override 
        call get_environment_variable("FDTD_STEPS_PER_BLOCK", env_val, status = stat) 
        if (stat == 0) then
            read(env_val, *) steps_per_block
        else 
            ! auto calculate with target ~2GB buffer 
            ! bytes = interior cells * 3 components * 4 bytes per real32 
            bytes_per_step = real(output_nx, 8) * output_ny * output_nz * 3 * 4 
            
            ! steps = 2GB (2^31 bytes) / bytes_per_step 
            steps_per_block = max(1, int(2147483648.0_8 / bytes_per_step))
            if (steps_per_block > 500 ) steps_per_block = 500 
        endif
        
        print*, "---- I/O Configuration ----"
        print*, "Block output: ", block_output
        print*, "Legacy output: ", legacy_output
        print*, "Block output dimensions: ", output_nx, output_ny, output_nz
        print*, "Block output precision: real32"
        print*, "Steps per block: ", steps_per_block
        
    end subroutine setup_io_params
    

end module seidartio
