program main 
    
    ! Types are defined in seidartio
    use seidartio 
    use cpmlfdtd 
    use seidart_types
    
    implicit none 
    
    type(Domain_Type) :: domain 
    type(Stiffness_Type), allocatable :: stiffness(:)
    type(Attenuation_Type), allocatable :: attenuation(:)
    type(Permittivity_Type), allocatable :: permittivity(:)
    type(Conductivity_Type), allocatable :: conductivity(:)
    type(Source_Type) :: seismic_source
    type(Source_Type) :: electromagnetic_source
    
    integer, allocatable :: geometry(:,:) 
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

    ! Print out to confirm the seismic flag status
    if (seismic) then
        print *, "Seismic flag is TRUE"
    else
        print *, "Seismic flag is FALSE"
    end if
    
    ! --------------------------------------------------------------------------
    ! Get going 
    call parse_json(trim(input_json_file), domain, &
            seismic_source, stiffness, attenuation, &
            permittivity, conductivity, electromagnetic_source)
    
    ! Read the geometry.dat file into memory
    call read_geometry('geometry.dat', domain, geometry)
    
    if (seismic) then 
        print *, "Writing seismic model parameters to Fortran unformatted binary files."
        call seismic_parameter_write(domain, geometry, stiffness, attenuation)
    else 
        print *, "Writing electromagnetic model parameters to Fortran unformatted binary files."
        call electromagnetic_parameter_write(domain, geometry, permittivity, conductivity)
    end if
    
    
    ! --------------------------------------------------------------------------
    
    ! --------------------------------------------------------------------------
end program main 
