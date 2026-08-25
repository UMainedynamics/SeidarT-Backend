module source_module
    use iso_c_binding
    use iso_fortran_env, only: real64
    use seidart_types
    use constants
    
    implicit none
    include 'fftw3.f03'

    private
    public :: init_source_weight_drop,   &
              init_source_explosive,     &
              init_source_double_couple, &
              init_source_clvd,          &
              init_source_plane_wave,    &
              inject_source_explicit_stage, &
              free_source

contains

    ! --------------------------------------------------------------------------
    !> 1. ACCELERATED WEIGHT DROP (AWD)
    subroutine init_source_weight_drop(src, domain, F_spec, freq_spec, n_spec, grid)
        type(Source_Type), intent(inout) :: src
        type(Domain_Type), intent(in) :: domain
        complex(real64), intent(in) :: F_spec(n_spec)
        real(real64), intent(in) :: freq_spec(n_spec)
        integer, intent(in) :: n_spec,
        type(spectral_grid_t), intent(inout) :: grid

        
        real(real64) :: az_rad, dip_rad

        az_rad  = src%azimuth * DEG2RAD
        dip_rad = src%dip * DEG2RAD

        src%force_vec(1) = src%amplitude * cos(dip_rad) * sin(az_rad)
        src%force_vec(2) = src%amplitude * cos(dip_rad) * cos(az_rad)
        src%force_vec(3) = src%amplitude * sin(dip_rad)

        call synthesize_from_spectrum(F_spec, freq_spec, n_spec, n_steps, dt, src%time_series)
        call precompute_spatial_kernel(src, grid, domain)
    end subroutine init_source_weight_drop

    ! --------------------------------------------------------------------------
    !> 2. EXPLOSIVE (ISOTROPIC MOMENT TENSOR)
    subroutine init_source_explosive(src, M0, &
                                     F_spec, freq_spec, n_spec, &
                                     n_steps, dt, grid, domain_dx, domain_dy, domain_dz)
        type(Source_Type), intent(inout) :: src
        type(Domain_Type), intent(in) :: domain
        
        real(real64), intent(in) :: M0
        complex(real64), intent(in) :: F_spec(n_spec)
        real(real64), intent(in) :: freq_spec(n_spec)
        integer, intent(in) :: n_spec, n_steps
        real(real64), intent(in) :: dt
        type(spectral_grid_t), intent(inout) :: grid
        
        
        src%moment_tensor = 0.0_real64
        src%moment_tensor(1) = M0
        src%moment_tensor(2) = M0
        src%moment_tensor(3) = M0

        call synthesize_from_spectrum(F_spec, freq_spec, n_spec, n_steps, dt, src%time_series)
        call precompute_spatial_kernel(src, grid, domain)
    end subroutine init_source_explosive

    ! --------------------------------------------------------------------------
    !> 3. DOUBLE COUPLE (FAULT RUPTURE)
    subroutine init_source_double_couple(src, domain, M0, strike_deg, dip_deg, rake_deg, &
                                         F_spec, freq_spec, n_spec, &
                                         n_steps, dt, grid, domain_dx, domain_dy, domain_dz)
        type(Source_Type), intent(inout) :: src
        type(Domain_Type)
        real(real64), intent(in) :: M0
        real(real64), intent(in) :: strike_deg, dip_deg, rake_deg
        complex(real64), intent(in) :: F_spec(n_spec)
        real(real64), intent(in) :: freq_spec(n_spec)
        integer, intent(in) :: n_spec
        real(real64), intent(in) :: dt
        type(spectral_grid_t), intent(inout) :: grid
        real(real64), intent(in) :: domain_dx, domain_dy, domain_dz

        real(real64), parameter :: DEG2RAD = 3.14159265358979323846_real64 / 180.0_real64
        real(real64) :: phi, del, lam
        real(real64) :: s_phi, c_phi, s_2phi, c_2phi
        real(real64) :: s_del, c_del, s_2del, c_2del
        real(real64) :: s_lam, c_lam

        phi = strike_deg * DEG2RAD
        del = dip_deg    * DEG2RAD
        lam = rake_deg   * DEG2RAD

        s_phi = sin(phi); c_phi = cos(phi)
        s_2phi = sin(2.0_real64 * phi); c_2phi = cos(2.0_real64 * phi)
        s_del = sin(del); c_del = cos(del)
        s_2del = sin(2.0_real64 * del); c_2del = cos(2.0_real64 * del)
        s_lam = sin(lam); c_lam = cos(lam)

        src%moment_tensor(1) = -M0 * (s_del * c_lam * s_2phi + s_2del * s_lam * s_phi**2)
        src%moment_tensor(2) =  M0 * (s_del * c_lam * s_2phi - s_2del * s_lam * c_phi**2)
        src%moment_tensor(3) =  M0 * (s_2del * s_lam)
        src%moment_tensor(4) = -M0 * (c_del * c_lam * s_phi - c_2del * s_lam * c_phi)
        src%moment_tensor(5) = -M0 * (c_del * c_lam * c_phi + c_2del * s_lam * s_phi)
        src%moment_tensor(6) =  M0 * (s_del * c_lam * c_2phi + 0.5_real64 * s_2del * s_lam * s_2phi)

        call synthesize_from_spectrum(F_spec, freq_spec, n_spec, n_steps, dt, src%time_series)
        call precompute_spatial_kernel(src, grid, domain_dx, domain_dy, domain_dz)
    end subroutine init_source_double_couple

    ! --------------------------------------------------------------------------
    !> 4. COMPENSATED LINEAR VECTOR DIPOLE (CLVD)
    subroutine init_source_clvd(src, domain, M0, axis_azimuth_deg, axis_plunge_deg, &
                                F_spec, freq_spec, n_spec, &
                                n_steps, dt, grid, domain_dx, domain_dy, domain_dz)
        type(source_t), intent(inout) :: src
        real(real64), intent(in) :: xs, ys, zs
        real(real64), intent(in) :: M0
        real(real64), intent(in) :: axis_azimuth_deg, axis_plunge_deg
        complex(real64), intent(in) :: F_spec(n_spec)
        real(real64), intent(in) :: freq_spec(n_spec)
        integer, intent(in) :: n_spec, n_steps
        real(real64), intent(in) :: dt
        type(spectral_grid_t), intent(inout) :: grid
        real(real64), intent(in) :: domain_dx, domain_dy, domain_dz

        real(real64), parameter :: DEG2RAD = 3.14159265358979323846_real64 / 180.0_real64
        real(real64) :: az_rad, pl_rad, ex, ey, ez

        src%xs = xs; src%ys = ys; src%zs = zs
        src%source_type = 2
        src%n_steps = n_steps
        src%dt = dt

        az_rad = axis_azimuth_deg * DEG2RAD
        pl_rad = axis_plunge_deg  * DEG2RAD

        ex = cos(pl_rad) * sin(az_rad)
        ey = cos(pl_rad) * cos(az_rad)
        ez = sin(pl_rad)

        src%moment_tensor(1) = M0 * (3.0_real64 * ex * ex - 1.0_real64)
        src%moment_tensor(2) = M0 * (3.0_real64 * ey * ey - 1.0_real64)
        src%moment_tensor(3) = M0 * (3.0_real64 * ez * ez - 1.0_real64)
        src%moment_tensor(4) = M0 * (3.0_real64 * ey * ez)
        src%moment_tensor(5) = M0 * (3.0_real64 * ex * ez)
        src%moment_tensor(6) = M0 * (3.0_real64 * ex * ey)

        call synthesize_from_spectrum(F_spec, freq_spec, n_spec, n_steps, dt, src%time_series)
        call precompute_spatial_kernel(src, grid, domain_dx, domain_dy, domain_dz)
    end subroutine init_source_clvd

    ! --------------------------------------------------------------------------
    !> 5. INCIDENT PLANE WAVE (BOUNDARY TFSF INJECTION)
    !> Injects a coherent plane wave across the inner boundary of the sponge layer.
    !> prop_azimuth_deg: 0=North (+y), 90=East (+x)
    !> prop_dip_deg: 0=Horizontal, 90=Directly downward (+z)
    !> pol_type: 'P', 'SV', 'SH'
    subroutine init_source_plane_wave(src, nx, ny, nz, dx, dy, dz, cpml_nodes, &
                                      c_background, prop_azimuth_deg, prop_dip_deg, pol_type, &
                                      F_spec, freq_spec, n_spec, n_steps, dt)
        type(source_t), intent(inout) :: src
        integer, intent(in) :: nx, ny, nz, cpml_nodes
        real(real64), intent(in) :: dx, dy, dz, c_background
        real(real64), intent(in) :: prop_azimuth_deg, prop_dip_deg
        character(len=*), intent(in) :: pol_type
        complex(real64), intent(in) :: F_spec(n_spec)
        real(real64), intent(in) :: freq_spec(n_spec)
        integer, intent(in) :: n_spec, n_steps
        real(real64), intent(in) :: dt

        real(real64), parameter :: DEG2RAD = 3.14159265358979323846_real64 / 180.0_real64
        real(real64) :: az_rad, dip_rad
        real(real64) :: px, py, pz, sv_x, sv_y, sv_z, sh_x, sh_y, sh_z
        real(real64) :: rx, ry, rz, dot_val
        integer :: i, j, k, i_min, i_max, j_min, j_max, k_min, k_max

        src%source_type = 3
        src%n_steps = n_steps
        src%dt = dt
        src%c_phase = c_background
        src%pml_thick = cpml_nodes

        az_rad  = prop_azimuth_deg * DEG2RAD
        dip_rad = prop_dip_deg     * DEG2RAD

        ! Direction unit vector p (x=East, y=North, z=Down)
        px = cos(dip_rad) * sin(az_rad)
        py = cos(dip_rad) * cos(az_rad)
        pz = sin(dip_rad)
        src%p_dir = (/ px, py, pz /)

        ! Polarization vectors
        select case (trim(pol_type))
        case ('P', 'COMPRESSIONAL', 'ACOUSTIC')
            src%e_pol = src%p_dir
        case ('SH', 'SHEAR_HORIZONTAL')
            sh_x =  cos(az_rad)
            sh_y = -sin(az_rad)
            sh_z =  0.0_real64
            src%e_pol = (/ sh_x, sh_y, sh_z /)
        case ('SV', 'SHEAR_VERTICAL')
            sv_x = -sin(dip_rad) * sin(az_rad)
            sv_y = -sin(dip_rad) * cos(az_rad)
            sv_z =  cos(dip_rad)
            src%e_pol = (/ sv_x, sv_y, sv_z /)
        case default
            src%e_pol = src%p_dir
        end select

        ! Reference plane entry corner
        src%r0_ref(1) = merge(real(cpml_nodes+1, real64)*dx, real(nx-cpml_nodes, real64)*dx, px >= 0.0_real64)
        src%r0_ref(2) = merge(real(cpml_nodes+1, real64)*dy, real(ny-cpml_nodes, real64)*dy, py >= 0.0_real64)
        src%r0_ref(3) = merge(real(cpml_nodes+1, real64)*dz, real(nz-cpml_nodes, real64)*dz, pz >= 0.0_real64)

        allocate(src%time_delay_3d(nx, ny, nz))
        allocate(src%injection_mask(nx, ny, nz))
        src%injection_mask = .false.

        i_min = cpml_nodes + 2; i_max = nx - cpml_nodes - 1
        j_min = cpml_nodes + 2; j_max = ny - cpml_nodes - 1
        k_min = cpml_nodes + 2; k_max = nz - cpml_nodes - 1

        ! Mark TFSF boundary faces
        do k = 1, nz
            rz = real(k - 1, real64) * dz
            do j = 1, ny
                ry = real(j - 1, real64) * dy
                do i = 1, nx
                    rx = real(i - 1, real64) * dx

                    ! Compute geometric phase travel-time delay
                    dot_val = (rx - src%r0_ref(1))*px + (ry - src%r0_ref(2))*py + (rz - src%r0_ref(3))*pz
                    src%time_delay_3d(i, j, k) = dot_val / src%c_phase

                    ! Active on the boundary boundary layer of the inner box
                    if ((i == i_min .or. i == i_max .or. &
                         j == j_min .or. j == j_max .or. &
                         k == k_min .or. k == k_max) .and. &
                        (i >= i_min .and. i <= i_max .and. &
                         j >= j_min .and. j <= j_max .and. &
                         k >= k_min .and. k <= k_max)) then
                        src%injection_mask(i, j, k) = .true.
                    end if
                end do
            end do
        end do

        call synthesize_from_spectrum(F_spec, freq_spec, n_spec, n_steps, dt, src%time_series)

    end subroutine init_source_plane_wave

    ! --------------------------------------------------------------------------
    !> Transforms 1D frequency spectrum to real time domain
    subroutine synthesize_from_spectrum(F_in, freq_in, n_in, n_steps, dt, s_time)
        complex(real64), intent(in)  :: F_in(n_in)
        real(real64), intent(in)     :: freq_in(n_in)
        integer, intent(in)          :: n_in, n_steps
        real(real64), intent(in)     :: dt
        real(real64), allocatable, intent(out) :: s_time(:)

        complex(real64), allocatable :: H_unif(:)
        type(C_PTR) :: plan_1d
        real(real64) :: df_target, target_freq, f_frac, norm_factor
        integer :: n_freq, k, j

        if (allocated(s_time)) deallocate(s_time)
        allocate(s_time(n_steps))

        n_freq = n_steps / 2 + 1
        df_target = 1.0_real64 / (real(n_steps, real64) * dt)
        allocate(H_unif(n_freq))
        H_unif = cmplx(0.0_real64, 0.0_real64, kind=real64)

        do k = 1, n_freq
            target_freq = real(k - 1, real64) * df_target

            if (target_freq >= freq_in(1) .and. target_freq <= freq_in(n_in)) then
                do j = 1, n_in - 1
                    if (target_freq >= freq_in(j) .and. target_freq <= freq_in(j+1)) then
                        f_frac = (target_freq - freq_in(j)) / max(freq_in(j+1) - freq_in(j), 1.0e-14_real64)
                        H_unif(k) = (1.0_real64 - f_frac) * F_in(j) + f_frac * F_in(j+1)
                        exit
                    end if
                end do
            else if (target_freq < freq_in(1)) then
                H_unif(k) = F_in(1)
            else
                H_unif(k) = cmplx(0.0_real64, 0.0_real64, kind=real64)
            end if
        end do

        norm_factor = 1.0_real64 / real(n_steps, real64)
        H_unif = H_unif * norm_factor

        plan_1d = fftw_plan_dft_c2r_1d(n_steps, H_unif, s_time, FFTW_ESTIMATE)
        call fftw_execute_dft_c2r(plan_1d, H_unif, s_time)
        call fftw_destroy_plan(plan_1d)

        deallocate(H_unif)
    end subroutine synthesize_from_spectrum

    ! --------------------------------------------------------------------------
    !> Builds compact 3D spatial subgrid kernel using exact k-space phase shifts
    subroutine precompute_spatial_kernel(src, grid, domain)
        type(source_t), intent(inout) :: src
        type(spectral_grid_t), intent(inout) :: grid
        real(Domain_Type), intent(in) :: domain

        real(real64), allocatable :: full_spatial(:,:,:)
        complex(real64), parameter :: imag_unit = (0.0_real64, 1.0_real64)
        real(real64) :: v_cell, phase
        integer :: i, j, k, di, dj, dk

        src%half_span = 4
        allocate(full_spatial(domain%nx, domain%ny, domain%nz))
        allocate(src%spatial_kernel(-src%half_span:src%half_span, &
                                    -src%half_span:src%half_span, &
                                    -src%half_span:src%half_span))

        src%isrc = nint(src%xs / dx) + 1
        src%jsrc = nint(src%ys / dy) + 1
        src%ksrc = nint(src%zs / dz) + 1

        v_cell = dx * dy * dz

        !$omp parallel do collapse(3) private(i,j,k,phase)
        do k = 1, domain%nz
            do j = 1, domain%ny
                do i = 1, domain%nkx
                    phase = -(grid%kx(i)*src%xs + grid%ky(j)*src%ys + grid%kz(k)*src%zs)
                    grid%F_hat(i, j, k) = (exp(imag_unit * phase) / v_cell) * grid%inv_n_total
                end do
            end do
        end do
        !$omp end parallel do

        call fftw_execute_dft_c2r(grid%plan_bwd_x, grid%F_hat, full_spatial)

        src%i1 = max(1, src%isrc - src%half_span)
        src%i2 = min(grid%nx, src%isrc + src%half_span)
        src%j1 = max(1, src%jsrc - src%half_span)
        src%j2 = min(grid%ny, src%jsrc + src%half_span)
        src%k1 = max(1, src%ksrc - src%half_span)
        src%k2 = min(grid%nz, src%ksrc + src%half_span)

        src%spatial_kernel = 0.0_real64
        do k = src%k1, src%k2
            dk = k - src%ksrc
            do j = src%j1, src%j2
                dj = j - src%jsrc
                do i = src%i1, src%i2
                    di = i - src%isrc
                    src%spatial_kernel(di, dj, dk) = full_spatial(i, j, k)
                end do
            end do
        end do

        deallocate(full_spatial)
    end subroutine precompute_spatial_kernel

    ! --------------------------------------------------------------------------
    !> Injects source amplitude into explicit stage RHS tensor T at stage time
    subroutine inject_source_explicit_stage(src, nx, ny, nz, T_stage, it, stage_time_frac, density_s)
        type(Source_Type), intent(in) :: src
        integer, intent(in) :: nx, ny, nz
        real(real64), intent(inout) :: T_stage(21, nx, ny, nz)
        integer, intent(in) :: it
        real(real64), intent(in) :: stage_time_frac
        real(real64), intent(in) :: density_s(nx, ny, nz)

        real(real64) :: s_val, t_curr, t_delayed, t_idx, w_ker, frac
        real(real64) :: z_impedance, rho_val
        integer :: idx_floor, idx_ceil, i, j, k, di, dj, dk

        t_curr = (real(it - 1, real64) + stage_time_frac) * src%dt
        
        if (src%source_type == "pw") then
            
        ! ---------------------------------------------
        ! Incident Boundary Plane Wave (TFSF Injection)
        ! ---------------------------------------------
    
            !$omp parallel do collapse(3) private(i,j,k,t_delayed,t_idx,idx_floor,idx_ceil,frac,s_val,rho_val,z_impedance)
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        if (src%injection_mask(i, j, k)) then
                            t_delayed = t_curr - src%time_delay_3d(i, j, k)

                            if (t_delayed >= 0.0_real64) then
                                t_idx = t_delayed / src%dt
                                idx_floor = int(t_idx) + 1
                                idx_ceil  = min(src%n_steps, idx_floor + 1)
                                frac = t_idx - real(idx_floor - 1, real64)

                                if (idx_floor >= 1 .and. idx_floor <= src%n_steps) then
                                    s_val = (1.0_real64 - frac) * src%time_series(idx_floor) + frac * src%time_series(idx_ceil)
                                else
                                    s_val = 0.0_real64
                                end if

                                rho_val = density_s(i, j, k)
                                z_impedance = rho_val * src%c_phase

                                ! 1. Particle Velocity Injection (1:3)
                                T_stage(1, i, j, k) = T_stage(1, i, j, k) + src%e_pol(1) * s_val
                                T_stage(2, i, j, k) = T_stage(2, i, j, k) + src%e_pol(2) * s_val
                                T_stage(3, i, j, k) = T_stage(3, i, j, k) + src%e_pol(3) * s_val

                                ! 2. Characteristic Plane Wave Stress Injection (7:12)
                                T_stage(7,  i, j, k) = T_stage(7,  i, j, k) - z_impedance * (src%p_dir(1) * src%e_pol(1)) * s_val
                                T_stage(8,  i, j, k) = T_stage(8,  i, j, k) - z_impedance * (src%p_dir(2) * src%e_pol(2)) * s_val
                                T_stage(9,  i, j, k) = T_stage(9,  i, j, k) - z_impedance * (src%p_dir(3) * src%e_pol(3)) * s_val
                                T_stage(10, i, j, k) = T_stage(10, i, j, k) - z_impedance * 0.5_real64 * &
                                    (src%p_dir(2)*src%e_pol(3) + src%p_dir(3)*src%e_pol(2)) * s_val
                                T_stage(11, i, j, k) = T_stage(11, i, j, k) - z_impedance * 0.5_real64 * &
                                    (src%p_dir(1)*src%e_pol(3) + src%p_dir(3)*src%e_pol(1)) * s_val
                                T_stage(12, i, j, k) = T_stage(12, i, j, k) - z_impedance * 0.5_real64 * &
                                    (src%p_dir(1)*src%e_pol(2) + src%p_dir(2)*src%e_pol(1)) * s_val
                            end if
                        end if
                    end do
                end do
            end do
            !$omp end parallel do
        else
            ! ------------------------------------------
            ! Local Point Sources (AWD & Moment Tensors)
            ! ------------------------------------------
            t_idx = t_curr / src%dt
            idx_floor = int(t_idx) + 1
            idx_ceil  = min(src%n_steps, idx_floor + 1)
            frac = t_idx - real(idx_floor - 1, real64)

            if (idx_floor >= 1 .and. idx_floor <= src%n_steps) then
                s_val = (1.0_real64 - frac) * src%time_series(idx_floor) + frac * src%time_series(idx_ceil)
            else
                return
            end if

            do k = src%k1, src%k2
                dk = k - src%ksrc
                do j = src%j1, src%j2
                    dj = j - src%jsrc
                    do i = src%i1, src%i2
                        di = i - src%isrc
                        w_ker = src%spatial_kernel(di, dj, dk)

                        if (src%source_type == "ac") then
                            ! Force -> Solid Accelerations (1:3)
                            T_stage(1, i, j, k) = T_stage(1, i, j, k) + (s_val * src%force_vec(1) * w_ker) / density_s(i, j, k)
                            T_stage(2, i, j, k) = T_stage(2, i, j, k) + (s_val * src%force_vec(2) * w_ker) / density_s(i, j, k)
                            T_stage(3, i, j, k) = T_stage(3, i, j, k) + (s_val * src%force_vec(3) * w_ker) / density_s(i, j, k)
                        else
                            ! Moment Tensor -> Stress Rates (7:12)
                            T_stage(7,  i, j, k) = T_stage(7,  i, j, k) + (s_val * src%moment_tensor(1) * w_ker)
                            T_stage(8,  i, j, k) = T_stage(8,  i, j, k) + (s_val * src%moment_tensor(2) * w_ker)
                            T_stage(9,  i, j, k) = T_stage(9,  i, j, k) + (s_val * src%moment_tensor(3) * w_ker)
                            T_stage(10, i, j, k) = T_stage(10, i, j, k) + (s_val * src%moment_tensor(4) * w_ker)
                            T_stage(11, i, j, k) = T_stage(11, i, j, k) + (s_val * src%moment_tensor(5) * w_ker)
                            T_stage(12, i, j, k) = T_stage(12, i, j, k) + (s_val * src%moment_tensor(6) * w_ker)
                        end if
                    end do
                end do
            end do

        end if

    end subroutine inject_source_explicit_stage

    ! --------------------------------------------------------------------------
    subroutine free_source(src)
        type(source_t), intent(inout) :: src
        if (allocated(src%time_series))     deallocate(src%time_series)
        if (allocated(src%spatial_kernel))  deallocate(src%spatial_kernel)
        if (allocated(src%time_delay_3d))   deallocate(src%time_delay_3d)
        if (allocated(src%injection_mask))  deallocate(src%injection_mask)
    end subroutine free_source

end module source_module