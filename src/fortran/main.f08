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
    character(len=256) :: arg
    character(len=256) :: input_json_file 
    logical :: seismic
    logical :: seismoacoustic
    integer :: argc, i
    character(len=256) :: key, value, density_method
    
    ! ---------------------------- Parse CLI -----------------------------------
    ! Get the number of command line arguments
    argc = command_argument_count() 

     ! Check if the right number of arguments are provided
    if (argc < 1) then
        print *, 'Usage: main <input.json> [seismic=<true|false>] [seismoacoustic=<true|false>]'
        stop
    end if
     
    ! Get the first command-line argument (input.json)
    call get_command_argument(1, input_json_file)
    
    ! Set defaults
    ! seismic flag 
    seismic = .true.
    density_method = 'arithmetic'
        
    ! loop over any further key=value flags
    do i = 2, argc
        call get_command_argument(i, arg)
        if (index(arg,'=') > 0) then
        key = adjustl(  trim(arg(1:index(arg,'=')-1))  )
        value = adjustl(  trim(arg(index(arg,'=')+1:))    )
        select case (key)
        case ('seismic')
            if (value == 'true') then
            seismic = .true.
            else if (value == 'false') then
            seismic = .false.
            end if
        case ('density_method')
            density_method = adjustl(value)
        case default
            ! ignore other flags
        end select
        end if
    end do

    ! Validate density_method
    select case (trim(density_method))
    case ('harmonic','geometric','arithmetic','none')
        ! OK
    case default
        print *, 'Error: invalid density_method "', density_method, '".'
        print *, '       Must be one of: harmonic, geometric, arithmetic, none'
        stop
    end select
    
    ! --------------------------------------------------------------------------
    ! Get going 
    call parse_json(trim(input_json_file), domain, seismic_source, electromagnetic_source)
    
    domain%nx = domain%nx + 2*domain%cpml
    domain%ny = domain%ny + 2*domain%cpml
    domain%nz = domain%nz + 2*domain%cpml
    
    !----------------------------------------------------------------------
    ! dispatch to the correct solver
    if (domain%dim == 2.5) then

        ! if (seismoacoustic) then
        ! print *, "Running 2.5D seismo‑acoustic model with", seismic_source%time_steps, "time steps"
        ! call seismoacoustic25(domain, seismic_source, .TRUE.)

        ! else if (seismic) then
        if (seismic) then
            print *, "Running 2.5D seismic model with", seismic_source%time_steps, "time steps"
            call seismic25(domain, seismic_source, density_method, .TRUE.)
        else
            print *, "Running 2.5D electromagnetic model with", electromagnetic_source%time_steps, "time steps"
            call electromag25(domain, electromagnetic_source, .TRUE.)
        end if
    else    ! dim == 2.0

        if (seismic) then
            print *, "Running 2D seismic model with", seismic_source%time_steps, "time steps"
            call seismic2(domain, seismic_source, density_method, .TRUE.)
        else
            print *, "Running 2D electromagnetic model with", electromagnetic_source%time_steps, "time steps"
            call electromag2(domain, electromagnetic_source, .TRUE.)
        end if

  end if

    
    ! --------------------------------------------------------------------------
end program main 
