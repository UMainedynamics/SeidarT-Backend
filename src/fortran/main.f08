program main 
    ! ==========================================================================
    ! 
    !
    ! Refer to the Makefile in the root directory for compile directions and 
    ! options. To run use the command:
    ! 
    !       seidartfdtd <input.json> [seismic=<true|false>]
    ! ==========================================================================
    
    ! Types are defined in seidartio
    use seidartio 
    use cpmlfdtd 
    use seidart_types
    
    implicit none 
    
    type(Domain_Type) :: domain 
    type(Source_Type) :: seismic_source
    type(Source_Type) :: electromagnetic_source
    
    ! integer, allocatable :: geometry(:,:) 
    character(len=256) :: input_json_file 
    logical :: seismic
    integer :: argc 
    character(len=256) :: seismic_arg 
    
    ! ---------------------------- Parse CLI -----------------------------------
    ! Get the number of command line arguments
    argc = command_argument_count() 
    
     ! Check if the right number of arguments are provided
    if (argc < 1) then
        print *, 'Usage: main <input.json> [seismic=<true|false>]'
        stop
    end if
     
    ! Get the first command-line argument (input.json)
    call get_command_argument(1, input_json_file)
    
    ! Set seismic flag to true by default
    seismic = .true.
    
    ! Optional: Get the seismic flag from the second argument
    if (argc > 1) then
        call get_command_argument(2, seismic_arg)
        
        ! Check if the argument is 'seismic=false' or 'seismic=true'
        if (trim(seismic_arg) == 'seismic=false') then
            seismic = .false.
        else if (trim(seismic_arg) == 'seismic=true') then
            seismic = .true.
        else
            print *, 'Invalid argument for seismic flag. Use seismic=<true|false>'
        stop
        end if
    end if

    
    ! --------------------------------------------------------------------------
    ! Get going 
    call parse_json(trim(input_json_file), domain, seismic_source, electromagnetic_source)
    
    domain%nx = domain%nx + 2*domain%cpml
    domain%ny = domain%ny + 2*domain%cpml
    domain%nz = domain%nz + 2*domain%cpml
    
    ! Read the geometry.dat file into memory
    ! call read_geometry('geometry.dat', domain, geometry)
    
    ! if (seismic) then 
    !     print *, "Writing seismic model parameters to Fortran unformatted binary files."
    !     call seismic_parameter_write(domain, geometry, stiffness, attenuation)
    ! else 
    !     print *, "Writing electromagnetic model parameters to Fortran unformatted binary files."
    !     call electromagnetic_parameter_write(domain, geometry, permittivity, conductivity)
    ! end if
    
    
    ! --------------------------------------------------------------------------
    if (domain%dim == 2.5) then
        if (seismic) then
            print *, "Running 2.5D seismic model with ", seismic_source%time_steps, " time steps"
            call seismic25(domain, seismic_source, .TRUE.)
        else 
            print *, "Running 2.5D electromagnetic model with", electromagnetic_source%time_steps, "time steps"
            call electromag25(domain, electromagnetic_source, .TRUE.)
        endif
    else
        if (seismic) then  
            print *, "Running 2D seismic model with", seismic_source%time_steps, " time steps"
            call seismic2(domain, seismic_source, .TRUE.)
        else 
            print *, "Running 2D electromagnetic model with", electromagnetic_source%time_steps, " time steps"
            call electromag2(domain, electromagnetic_source, .TRUE.)
        endif
    endif 
    
    ! --------------------------------------------------------------------------
end program main 
