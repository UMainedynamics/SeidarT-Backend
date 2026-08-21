module imex
    use iso_c_binding
    use iso_fortran_env, only: real64
    
    use spectral_operators
    use biot_physics
    use absorbing_boundary
    
    use seidartio
    use seidart_types
    use constants 
    use omp_lib
    
    implicit none
    include 'fftw3.f03'
    
    private
    public :: load_coefficients, ars_coefficients, biot_poroviscoelasticity3
    
contains
    
    ! --------------------------------------------------------------------------
    subroutine load_coefficients(nx, ny, nz, &
                                 C, drag_tensor, gamma_visco, &
                                 density_s, density_f, tau_relax, lwc, phi)
        integer, intent(in) :: nx, ny, nz
        real(real64), intent(out) :: C(21, nx, ny, nz)
        real(real64), intent(out) :: drag_tensor(6, nx, ny, nz), &
                                     gamma_visco(6, nx, ny, nz), &
                                     tau_relax(6, nx, ny, nz)
        real(real64), intent(out) :: density_s(nx, ny, nz), density_f(nx, ny, nz), &
                                     lwc(nx, ny, nz), phi(nx, ny, nz)
        
        ! Load stiffness coefficients (Upper-triangular Voigt packing)
        call material_rw3('c11.dat', C(1,:,:,:), .TRUE.)
        call material_rw3('c12.dat', C(2,:,:,:), .TRUE.)
        call material_rw3('c13.dat', C(3,:,:,:), .TRUE.)
        call material_rw3('c14.dat', C(4,:,:,:), .TRUE.)
        call material_rw3('c15.dat', C(5,:,:,:), .TRUE.)
        call material_rw3('c16.dat', C(6,:,:,:), .TRUE.)
        call material_rw3('c22.dat', C(7,:,:,:), .TRUE.)
        call material_rw3('c23.dat', C(8,:,:,:), .TRUE.)
        call material_rw3('c24.dat', C(9,:,:,:), .TRUE.)
        call material_rw3('c25.dat', C(10,:,:,:), .TRUE.)
        call material_rw3('c26.dat', C(11,:,:,:), .TRUE.)
        call material_rw3('c33.dat', C(12,:,:,:), .TRUE.)
        call material_rw3('c34.dat', C(13,:,:,:), .TRUE.)
        call material_rw3('c35.dat', C(14,:,:,:), .TRUE.)
        call material_rw3('c36.dat', C(15,:,:,:), .TRUE.)
        call material_rw3('c44.dat', C(16,:,:,:), .TRUE.)
        call material_rw3('c45.dat', C(17,:,:,:), .TRUE.)
        call material_rw3('c46.dat', C(18,:,:,:), .TRUE.)
        call material_rw3('c55.dat', C(19,:,:,:), .TRUE.)
        call material_rw3('c56.dat', C(20,:,:,:), .TRUE.)
        call material_rw3('c66.dat', C(21,:,:,:), .TRUE.)
        
        ! Load viscoelastic coefficients
        call material_rw3('gamma_x.dat',  gamma_visco(1,:,:,:), .TRUE.)
        call material_rw3('gamma_y.dat',  gamma_visco(2,:,:,:), .TRUE.)
        call material_rw3('gamma_z.dat',  gamma_visco(3,:,:,:), .TRUE.)
        call material_rw3('gamma_yz.dat', gamma_visco(4,:,:,:), .TRUE.)
        call material_rw3('gamma_xz.dat', gamma_visco(5,:,:,:), .TRUE.)
        call material_rw3('gamma_xy.dat', gamma_visco(6,:,:,:), .TRUE.)
        
        ! Load drag coefficients
        call material_rw3('drag_x.dat',  drag_tensor(1,:,:,:), .TRUE.)
        call material_rw3('drag_y.dat',  drag_tensor(2,:,:,:), .TRUE.)
        call material_rw3('drag_z.dat',  drag_tensor(3,:,:,:), .TRUE.)
        call material_rw3('drag_yz.dat', drag_tensor(4,:,:,:), .TRUE.)
        call material_rw3('drag_xz.dat', drag_tensor(5,:,:,:), .TRUE.)
        call material_rw3('drag_xy.dat', drag_tensor(6,:,:,:), .TRUE.)
        
        ! Load tau relaxation tensor
        call material_rw3('tau_x.dat',  tau_relax(1,:,:,:), .TRUE.)
        call material_rw3('tau_y.dat',  tau_relax(2,:,:,:), .TRUE.)
        call material_rw3('tau_z.dat',  tau_relax(3,:,:,:), .TRUE.)
        call material_rw3('tau_yz.dat', tau_relax(4,:,:,:), .TRUE.)
        call material_rw3('tau_xz.dat', tau_relax(5,:,:,:), .TRUE.)
        call material_rw3('tau_xy.dat', tau_relax(6,:,:,:), .TRUE.)
        
        call material_rw3('density_solid.dat',        density_s, .TRUE.)
        call material_rw3('density_fluid.dat',        density_f, .TRUE.)
        call material_rw3('porosity.dat',             phi,       .TRUE.)
        call material_rw3('liquid_water_content.dat', lwc,       .TRUE.)
        
    end subroutine load_coefficients
    
    ! --------------------------------------------------------------------------
    subroutine ars_coefficients(method, num_stages, A_exp, A_imp, b_exp, b_imp, c_vec)
        implicit none
        
        character(len=*), intent(in) :: method
        integer, intent(out) :: num_stages
        real(real64), allocatable, intent(out) :: A_exp(:,:), A_imp(:,:)
        real(real64), allocatable, intent(out) :: b_exp(:), b_imp(:), c_vec(:)
        
        real(real64) :: gamma, delta
        
        select case (trim(method))
        case ('ARS222', 'default')
            num_stages = 3
            allocate(A_exp(3,3), A_imp(3,3), b_exp(3), b_imp(3), c_vec(3))
            
            gamma = 1.0_real64 - 1.0_real64 / sqrt(2.0_real64)
            delta = 1.0_real64 - 1.0_real64 / (2.0_real64 * gamma)
            
            ! Explicit tableau (A_exp)
            A_exp = 0.0_real64 
            A_exp(2,1) = gamma 
            A_exp(3,1) = delta 
            A_exp(3,2) = 1.0_real64 - delta
            
            ! Implicit tableau (A_imp) - DIRK
            A_imp = 0.0_real64
            A_imp(2,2) = gamma
            A_imp(3,2) = 1.0_real64 - gamma
            A_imp(3,3) = gamma
            
            ! Stage weights
            b_exp = (/ 0.0_real64, 1.0_real64 - gamma, gamma /)
            b_imp = (/ 0.0_real64, 1.0_real64 - gamma, gamma /)
            c_vec = (/ 0.0_real64, gamma, 1.0_real64 /)
        case default
            error stop "Unknown IMEX method specified."
        end select
        
    end subroutine ars_coefficients
    
    ! -------------------------------------------------------------------------
    subroutine biot_poroviscoelasticity3(domain, source, method)
        implicit none 
        
        ! Input arguments 
        type(Domain_Type), intent(in) :: domain 
        type(Source_Type), intent(in) :: source
        character(len=*), intent(in)  :: method 
        
        type(spectral_grid_t) :: grid
        type(sponge_layer_t)  :: sponge
                
        ! Tableau allocations populated by ars_coefficients
        real(real64), allocatable :: A_exp(:,:), A_imp(:,:)
        real(real64), allocatable :: b_exp(:), b_imp(:), c_vec(:)
        
        ! Local variables
        integer :: nx, ny, nz, it, s, stage, num_stages
        real(real64) :: dt 
        
        real(real64), allocatable :: Q(:,:,:,:), Q_star(:,:,:,:), Q_stage(:,:,:,:)
        real(real64), allocatable :: T(:,:,:,:,:), H(:,:,:,:,:)
        
        real(real64), allocatable :: C(:,:,:,:), gamma_visco(:,:,:,:), &
                                     drag_tensor(:,:,:,:), tau_relax(:,:,:,:)
        real(real64), allocatable :: fluid_pressure(:,:,:), fluid_pressure_dot(:,:,:), &
                                     density_s(:,:,:), density_f(:,:,:), lwc(:,:,:), phi(:,:,:)
        
        ! Spatial gradient work arrays
        real(real64), allocatable :: dvx_dx(:,:,:), dvx_dy(:,:,:), dvx_dz(:,:,:)
        real(real64), allocatable :: dvy_dx(:,:,:), dvy_dy(:,:,:), dvy_dz(:,:,:)
        real(real64), allocatable :: dvz_dx(:,:,:), dvz_dy(:,:,:), dvz_dz(:,:,:)
        real(real64), allocatable :: dsxx_dx(:,:,:), dsxx_dy(:,:,:), dsxx_dz(:,:,:)
        real(real64), allocatable :: dsyy_dx(:,:,:), dsyy_dy(:,:,:), dsyy_dz(:,:,:)
        real(real64), allocatable :: dszz_dx(:,:,:), dszz_dy(:,:,:), dszz_dz(:,:,:)
        real(real64), allocatable :: dsyz_dx(:,:,:), dsyz_dy(:,:,:), dsyz_dz(:,:,:)
        real(real64), allocatable :: dsxz_dx(:,:,:), dsxz_dy(:,:,:), dsxz_dz(:,:,:)
        real(real64), allocatable :: dsxy_dx(:,:,:), dsxy_dy(:,:,:), dsxy_dz(:,:,:)
        real(real64), allocatable :: dp_dx(:,:,:),   dp_dy(:,:,:),   dp_dz(:,:,:)

        ! Plane wave source variables
        real(real64) :: XMIN, XMAX, XMID, YMIN, YMAX, YMID, ZMIN, ZMAX, ZMID
        real(real64) :: XLOC, YLOC, ZLOC, cbackground
        real(real64) :: r0(3), p(3), ehat(3)
        logical :: active(6)
        integer :: i_min, i_max, j_min, j_max, k_min, k_max
        real(real64), allocatable :: eig_array(:,:,:)
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:)
        real(real64), allocatable :: srcxx(:), srcyy(:), srczz(:), srcyz(:), srcxz(:), srcxy(:)
        
        ! -------------------------------------------------------------------------
        nx = domain%nx
        ny = domain%ny
        nz = domain%nz
        dt = domain%dt
        
        call ars_coefficients(method, num_stages, A_exp, A_imp, b_exp, b_imp, c_vec)
        
        ! Allocations 
        allocate(Q(21, nx, ny, nz), Q_star(21, nx, ny, nz), Q_stage(21, nx, ny, nz))
        allocate(T(21, nx, ny, nz, num_stages), H(21, nx, ny, nz, num_stages))
        
        allocate(C(21, nx, ny, nz))
        allocate(gamma_visco(6, nx, ny, nz))
        allocate(drag_tensor(6, nx, ny, nz), tau_relax(6, nx, ny, nz))
        allocate(fluid_pressure(nx, ny, nz), fluid_pressure_dot(nx, ny, nz))
        allocate(density_s(nx, ny, nz), density_f(nx, ny, nz))
        allocate(lwc(nx, ny, nz), phi(nx, ny, nz))
        
        ! Derivative allocations
        allocate(dvx_dx(nx,ny,nz), dvx_dy(nx,ny,nz), dvx_dz(nx,ny,nz))
        allocate(dvy_dx(nx,ny,nz), dvy_dy(nx,ny,nz), dvy_dz(nx,ny,nz))
        allocate(dvz_dx(nx,ny,nz), dvz_dy(nx,ny,nz), dvz_dz(nx,ny,nz))
        allocate(dsxx_dx(nx,ny,nz), dsxx_dy(nx,ny,nz), dsxx_dz(nx,ny,nz))
        allocate(dsyy_dx(nx,ny,nz), dsyy_dy(nx,ny,nz), dsyy_dz(nx,ny,nz))
        allocate(dszz_dx(nx,ny,nz), dszz_dy(nx,ny,nz), dszz_dz(nx,ny,nz))
        allocate(dsyz_dx(nx,ny,nz), dsyz_dy(nx,ny,nz), dsyz_dz(nx,ny,nz))
        allocate(dsxz_dx(nx,ny,nz), dsxz_dy(nx,ny,nz), dsxz_dz(nx,ny,nz))
        allocate(dsxy_dx(nx,ny,nz), dsxy_dy(nx,ny,nz), dsxy_dz(nx,ny,nz))
        allocate(dp_dx(nx,ny,nz),   dp_dy(nx,ny,nz),   dp_dz(nx,ny,nz))
        
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps))
        allocate(srcxx(source%time_steps), srcyy(source%time_steps), srczz(source%time_steps))
        allocate(srcyz(source%time_steps), srcxz(source%time_steps), srcxy(source%time_steps))
        allocate(eig_array(nx, ny, nz))
        
        ! Initialize FFTW grid
        call init_spectral_grid(grid, nx, ny, nz, domain%dx, domain%dy, domain%dz)
        
        ! Load material parameters
        call load_coefficients(nx, ny, nz, C, drag_tensor, gamma_visco, &
                               density_s, density_f, tau_relax, lwc, phi)
        
        Q = 0.0_real64
        fluid_pressure = 0.0_real64
        
        ! ------------------------------------------------------------------------
        ! Initialize the source 
        select case (trim(source%type))
        case('ac')
            init_source_weight_drop()
        case('tnt')
            init_source_explosive()     
        case('dc')
            init_source_double_couple() 
        case('clvd')
            init_source_clvd()          
        case('pw')
            init_source_plane_wave()
        end select 
        
        ! ------------------------------------------------------------------------
        ! Initialize absorbing boundary
        call init_sponge_layer(sponge, nx, ny, nz, domain%cpml, 0.025_real64)
        
        ! ------------------------------------------------------------------------
        ! ----- IMEX Time Loop -----
        do it = 1, source%time_steps
    
            do stage = 1, num_stages
                
                ! 1. Assemble intermediate explicit state
                Q_star = Q 
                do s = 1, stage - 1
                    if (abs(A_exp(stage, s)) > 1.0e-14_real64) then 
                        Q_star = Q_star + (dt * A_exp(stage, s)) * T(:,:,:,:,s) 
                    end if
                    if (abs(A_imp(stage, s)) > 1.0e-14_real64) then 
                        Q_star = Q_star + (dt * A_imp(stage, s)) * H(:,:,:,:,s)
                    end if 
                end do 
                
                ! 2. Solve local implicit stage for Q_stage 
                if (abs(A_imp(stage, stage)) > 1.0e-14_real64) then 
                    call implicit_drag_kernel(nx, ny, nz, Q_star, Q_stage, &
                                              drag_tensor, tau_relax, &
                                              density_s, density_f, lwc, &
                                              dt * A_imp(stage, stage))
                else
                    Q_stage = Q_star 
                end if
                
                ! 3. Compute spatial gradients via FFTW
                call grad3d(grid, Q_stage(1,:,:,:), dvx_dx,  dvx_dy,  dvx_dz)
                call grad3d(grid, Q_stage(2,:,:,:), dvy_dx,  dvy_dy,  dvy_dz)
                call grad3d(grid, Q_stage(3,:,:,:), dvz_dx,  dvz_dy,  dvz_dz)
                call grad3d(grid, Q_stage(7,:,:,:), dsxx_dx, dsxx_dy, dsxx_dz)
                call grad3d(grid, Q_stage(8,:,:,:), dsyy_dx, dsyy_dy, dsyy_dz)
                call grad3d(grid, Q_stage(9,:,:,:), dszz_dx, dszz_dy, dszz_dz)
                call grad3d(grid, Q_stage(10,:,:,:), dsyz_dx, dsyz_dy, dsyz_dz)
                call grad3d(grid, Q_stage(11,:,:,:), dsxz_dx, dsxz_dy, dsxz_dz)
                call grad3d(grid, Q_stage(12,:,:,:), dsxy_dx, dsxy_dy, dsxy_dz)
                call grad3d(grid, fluid_pressure,   dp_dx,   dp_dy,   dp_dz)
                
                ! 4. Evaluate explicit spatial RHS (T)
                call explicit_constitutive_kernel(nx, ny, nz, Q_stage, &
                                                  dvx_dx, dvx_dy, dvx_dz, &
                                                  dvy_dx, dvy_dy, dvy_dz, &
                                                  dvz_dx, dvz_dy, dvz_dz, &
                                                  dsxx_dx, dsyy_dy, dszz_dz, &
                                                  dsyz_dy, dsyz_dz, &
                                                  dsxz_dx, dsxz_dz, &
                                                  dsxy_dx, dsxy_dy, &
                                                  dp_dx,   dp_dy,   dp_dz, &
                                                  C, gamma_visco, &
                                                  density_s, density_f, phi, lwc, &
                                                  domain%bulk_mod_fluid, &
                                                  fluid_pressure_dot, T(:,:,:,:,stage))
                
                ! 5. Evaluate stiff evaluation H
                if (abs(A_imp(stage, stage)) > 1.0e-14_real64) then 
                    H(:,:,:,:,stage) = (Q_stage - Q_star) / (dt * A_imp(stage, stage))
                else
                    H(:,:,:,:,stage) = 0.0_real64
                end if
            end do
            
            ! Advance solution to next time step
            do s = 1, num_stages
                if (abs(b_exp(s)) > 1.0e-14_real64) then 
                    Q = Q + (dt * b_exp(s)) * T(:,:,:,:,s)
                end if 
                if (abs(b_imp(s)) > 1.0e-14_real64) then 
                    Q = Q + (dt * b_imp(s)) * H(:,:,:,:,s)
                end if 
            end do
            
            fluid_pressure = fluid_pressure + dt * fluid_pressure_dot 
            
            ! Attenuate waves at absorbing boundary
            call apply_sponge_damping(sponge, Q, fluid_pressure)
            
        end do
        
        call free_spectral_grid(grid)
        call free_sponge_layer(sponge)
        
    end subroutine biot_poroviscoelasticity3
    
end module imex