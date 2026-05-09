module discontinuous_galerkin_methods

    use iso_fortran_env, only: real64

    implicit none

    integer, parameter :: DG_P = 4
    integer, parameter :: DG_NP = DG_P + 1

contains

    pure logical function is_fluid_injection_source(source_type) result(is_fluid)
        character(len=*), intent(in) :: source_type
        character(len=:), allocatable :: st

        st = trim(adjustl(source_type))
        is_fluid = st == 'fluid' .or. st == 'fluid_injection' .or. &
            st == 'fluid-injection' .or. st == 'fi' .or. st == 'q'
    end function is_fluid_injection_source

    pure logical function is_pressure_injection_source(source_type) result(is_pressure)
        character(len=*), intent(in) :: source_type
        character(len=:), allocatable :: st

        st = trim(adjustl(source_type))
        is_pressure = st == 'pressure' .or. st == 'pressure_injection' .or. &
            st == 'pressure-injection' .or. st == 'pi' .or. st == 'p'
    end function is_pressure_injection_source

    pure real(real64) function positive(value, floor_value)
        real(real64), intent(in) :: value, floor_value

        positive = max(value, floor_value)
    end function positive

    subroutine initialize_gll4(nodes, weights, deriv)
        real(real64), intent(out) :: nodes(0:DG_P), weights(0:DG_P)
        real(real64), intent(out) :: deriv(0:DG_P,0:DG_P)

        real(real64) :: bary(0:DG_P), prod
        integer :: i, j, m

        nodes(0) = -1.0_real64
        nodes(1) = -sqrt(3.0_real64 / 7.0_real64)
        nodes(2) = 0.0_real64
        nodes(3) = sqrt(3.0_real64 / 7.0_real64)
        nodes(4) = 1.0_real64

        weights(0) = 1.0_real64 / 10.0_real64
        weights(1) = 49.0_real64 / 90.0_real64
        weights(2) = 32.0_real64 / 45.0_real64
        weights(3) = 49.0_real64 / 90.0_real64
        weights(4) = 1.0_real64 / 10.0_real64

        do i = 0, DG_P
            prod = 1.0_real64
            do m = 0, DG_P
                if (m /= i) prod = prod * (nodes(i) - nodes(m))
            enddo
            bary(i) = 1.0_real64 / prod
        enddo

        deriv(:,:) = 0.0_real64
        do i = 0, DG_P
            do j = 0, DG_P
                if (i /= j) deriv(i,j) = bary(j) / (bary(i) * (nodes(i) - nodes(j)))
            enddo
            deriv(i,i) = -sum(deriv(i,:))
        enddo
    end subroutine initialize_gll4

    subroutine fill_modal_from_cell(cell, nodal)
        real(real64), intent(in) :: cell(:,:)
        real(real64), intent(out) :: nodal(0:,0:,:,:)

        integer :: a, b, i, k

        !$omp parallel do collapse(2) private(a,b) schedule(static)
        do k = 1, size(cell, 2)
            do i = 1, size(cell, 1)
                do b = 0, DG_P
                    do a = 0, DG_P
                        nodal(a,b,i,k) = cell(i,k)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine fill_modal_from_cell

    subroutine cell_average(nodal, weights, cell)
        real(real64), intent(in) :: nodal(0:,0:,:,:)
        real(real64), intent(in) :: weights(0:DG_P)
        real(real64), intent(out) :: cell(:,:)

        integer :: a, b, i, k
        real(real64) :: wsum

        wsum = 4.0_real64
        cell(:,:) = 0.0_real64
        !$omp parallel do collapse(2) private(a,b) schedule(static)
        do k = 1, size(cell, 2)
            do i = 1, size(cell, 1)
                do b = 0, DG_P
                    do a = 0, DG_P
                        cell(i,k) = cell(i,k) + weights(a) * weights(b) * nodal(a,b,i,k)
                    enddo
                enddo
                cell(i,k) = cell(i,k) / wsum
            enddo
        enddo
        !$omp end parallel do
    end subroutine cell_average

    subroutine dg_derivatives(field, deriv, dx, dz, dfdx, dfdz)
        real(real64), intent(in) :: field(0:,0:,:,:)
        real(real64), intent(in) :: deriv(0:DG_P,0:DG_P)
        real(real64), intent(in) :: dx, dz
        real(real64), intent(out) :: dfdx(0:,0:,:,:)
        real(real64), intent(out) :: dfdz(0:,0:,:,:)

        integer :: a, b, m, i, k

        dfdx(:,:,:,:) = 0.0_real64
        dfdz(:,:,:,:) = 0.0_real64
        !$omp parallel do collapse(2) private(a,b,m) schedule(static)
        do k = 1, size(field, 4)
            do i = 1, size(field, 3)
                do b = 0, DG_P
                    do a = 0, DG_P
                        do m = 0, DG_P
                            dfdx(a,b,i,k) = dfdx(a,b,i,k) + deriv(a,m) * field(m,b,i,k)
                            dfdz(a,b,i,k) = dfdz(a,b,i,k) + deriv(b,m) * field(a,m,i,k)
                        enddo
                        dfdx(a,b,i,k) = 2.0_real64 * dfdx(a,b,i,k) / dx
                        dfdz(a,b,i,k) = 2.0_real64 * dfdz(a,b,i,k) / dz
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine dg_derivatives

    subroutine apply_rusanov_penalty(u, rhs, weights, dx, dz, speed)
        real(real64), intent(in) :: u(0:,0:,:,:)
        real(real64), intent(inout) :: rhs(0:,0:,:,:)
        real(real64), intent(in) :: weights(0:DG_P)
        real(real64), intent(in) :: dx, dz, speed

        integer :: b, i, k
        integer :: nx, nz
        real(real64) :: left_state, right_state, jump

        nx = size(u, 3)
        nz = size(u, 4)

        !$omp parallel do collapse(2) private(b,left_state,right_state,jump) schedule(static)
        do k = 1, nz
            do i = 1, nx
                do b = 0, DG_P
                    if (i > 1) then
                        left_state = u(DG_P,b,i-1,k)
                    else
                        left_state = u(0,b,i,k)
                    endif
                    jump = u(0,b,i,k) - left_state
                    rhs(0,b,i,k) = rhs(0,b,i,k) - speed * jump * 2.0_real64 / (dx * weights(0))

                    if (i < nx) then
                        right_state = u(0,b,i+1,k)
                    else
                        right_state = u(DG_P,b,i,k)
                    endif
                    jump = u(DG_P,b,i,k) - right_state
                    rhs(DG_P,b,i,k) = rhs(DG_P,b,i,k) - speed * jump * 2.0_real64 / (dx * weights(DG_P))
                enddo
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do collapse(2) private(b,left_state,right_state,jump) schedule(static)
        do k = 1, nz
            do i = 1, nx
                do b = 0, DG_P
                    if (k > 1) then
                        left_state = u(b,DG_P,i,k-1)
                    else
                        left_state = u(b,0,i,k)
                    endif
                    jump = u(b,0,i,k) - left_state
                    rhs(b,0,i,k) = rhs(b,0,i,k) - speed * jump * 2.0_real64 / (dz * weights(0))

                    if (k < nz) then
                        right_state = u(b,0,i,k+1)
                    else
                        right_state = u(b,DG_P,i,k)
                    endif
                    jump = u(b,DG_P,i,k) - right_state
                    rhs(b,DG_P,i,k) = rhs(b,DG_P,i,k) - speed * jump * 2.0_real64 / (dz * weights(DG_P))
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine apply_rusanov_penalty

    subroutine apply_rusanov_penalty_masked(u, rhs, active, weights, dx, dz, speed)
        real(real64), intent(in) :: u(0:,0:,:,:)
        real(real64), intent(inout) :: rhs(0:,0:,:,:)
        real(real64), intent(in) :: active(:,:)
        real(real64), intent(in) :: weights(0:DG_P)
        real(real64), intent(in) :: dx, dz, speed

        integer :: b, i, k
        integer :: nx, nz
        real(real64) :: left_state, right_state, jump

        nx = size(u, 3)
        nz = size(u, 4)

        !$omp parallel do collapse(2) private(b,left_state,right_state,jump) schedule(static)
        do k = 1, nz
            do i = 1, nx
                if (active(i,k) <= 0.5_real64) cycle
                do b = 0, DG_P
                    if (i > 1) then
                        if (active(i-1,k) > 0.5_real64) then
                            left_state = u(DG_P,b,i-1,k)
                        else
                            left_state = u(0,b,i,k)
                        endif
                    else
                        left_state = u(0,b,i,k)
                    endif
                    jump = u(0,b,i,k) - left_state
                    rhs(0,b,i,k) = rhs(0,b,i,k) - speed * jump * 2.0_real64 / (dx * weights(0))

                    if (i < nx) then
                        if (active(i+1,k) > 0.5_real64) then
                            right_state = u(0,b,i+1,k)
                        else
                            right_state = u(DG_P,b,i,k)
                        endif
                    else
                        right_state = u(DG_P,b,i,k)
                    endif
                    jump = u(DG_P,b,i,k) - right_state
                    rhs(DG_P,b,i,k) = rhs(DG_P,b,i,k) - speed * jump * 2.0_real64 / (dx * weights(DG_P))
                enddo
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do collapse(2) private(b,left_state,right_state,jump) schedule(static)
        do k = 1, nz
            do i = 1, nx
                if (active(i,k) <= 0.5_real64) cycle
                do b = 0, DG_P
                    if (k > 1) then
                        if (active(i,k-1) > 0.5_real64) then
                            left_state = u(b,DG_P,i,k-1)
                        else
                            left_state = u(b,0,i,k)
                        endif
                    else
                        left_state = u(b,0,i,k)
                    endif
                    jump = u(b,0,i,k) - left_state
                    rhs(b,0,i,k) = rhs(b,0,i,k) - speed * jump * 2.0_real64 / (dz * weights(0))

                    if (k < nz) then
                        if (active(i,k+1) > 0.5_real64) then
                            right_state = u(b,0,i,k+1)
                        else
                            right_state = u(b,DG_P,i,k)
                        endif
                    else
                        right_state = u(b,DG_P,i,k)
                    endif
                    jump = u(b,DG_P,i,k) - right_state
                    rhs(b,DG_P,i,k) = rhs(b,DG_P,i,k) - speed * jump * 2.0_real64 / (dz * weights(DG_P))
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine apply_rusanov_penalty_masked

    subroutine apply_rusanov_penalty_masked_local(u, rhs, active, weights, dx, dz, speed)
        real(real64), intent(in) :: u(0:,0:,:,:)
        real(real64), intent(inout) :: rhs(0:,0:,:,:)
        real(real64), intent(in) :: active(:,:), speed(:,:)
        real(real64), intent(in) :: weights(0:DG_P)
        real(real64), intent(in) :: dx, dz

        integer :: b, i, k
        integer :: nx, nz
        real(real64) :: left_state, right_state, jump, face_speed

        nx = size(u, 3)
        nz = size(u, 4)

        !$omp parallel do collapse(2) private(b,left_state,right_state,jump,face_speed) schedule(static)
        do k = 1, nz
            do i = 1, nx
                if (active(i,k) <= 0.5_real64) cycle
                do b = 0, DG_P
                    if (i > 1) then
                        if (active(i-1,k) > 0.5_real64) then
                            left_state = u(DG_P,b,i-1,k)
                            face_speed = max(speed(i,k), speed(i-1,k))
                        else
                            left_state = u(0,b,i,k)
                            face_speed = 0.0_real64
                        endif
                    else
                        left_state = u(0,b,i,k)
                        face_speed = 0.0_real64
                    endif
                    jump = u(0,b,i,k) - left_state
                    rhs(0,b,i,k) = rhs(0,b,i,k) - face_speed * jump * 2.0_real64 / (dx * weights(0))

                    if (i < nx) then
                        if (active(i+1,k) > 0.5_real64) then
                            right_state = u(0,b,i+1,k)
                            face_speed = max(speed(i,k), speed(i+1,k))
                        else
                            right_state = u(DG_P,b,i,k)
                            face_speed = 0.0_real64
                        endif
                    else
                        right_state = u(DG_P,b,i,k)
                        face_speed = 0.0_real64
                    endif
                    jump = u(DG_P,b,i,k) - right_state
                    rhs(DG_P,b,i,k) = rhs(DG_P,b,i,k) - face_speed * jump * 2.0_real64 / (dx * weights(DG_P))
                enddo
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do collapse(2) private(b,left_state,right_state,jump,face_speed) schedule(static)
        do k = 1, nz
            do i = 1, nx
                if (active(i,k) <= 0.5_real64) cycle
                do b = 0, DG_P
                    if (k > 1) then
                        if (active(i,k-1) > 0.5_real64) then
                            left_state = u(b,DG_P,i,k-1)
                            face_speed = max(speed(i,k), speed(i,k-1))
                        else
                            left_state = u(b,0,i,k)
                            face_speed = 0.0_real64
                        endif
                    else
                        left_state = u(b,0,i,k)
                        face_speed = 0.0_real64
                    endif
                    jump = u(b,0,i,k) - left_state
                    rhs(b,0,i,k) = rhs(b,0,i,k) - face_speed * jump * 2.0_real64 / (dz * weights(0))

                    if (k < nz) then
                        if (active(i,k+1) > 0.5_real64) then
                            right_state = u(b,0,i,k+1)
                            face_speed = max(speed(i,k), speed(i,k+1))
                        else
                            right_state = u(b,DG_P,i,k)
                            face_speed = 0.0_real64
                        endif
                    else
                        right_state = u(b,DG_P,i,k)
                        face_speed = 0.0_real64
                    endif
                    jump = u(b,DG_P,i,k) - right_state
                    rhs(b,DG_P,i,k) = rhs(b,DG_P,i,k) - face_speed * jump * 2.0_real64 / (dz * weights(DG_P))
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine apply_rusanov_penalty_masked_local

    subroutine damp_biot_darcy_fields(qx, qz, active, viscosity, permeability_x, permeability_z, rho_fluid, dt_step)
        real(real64), intent(inout) :: qx(0:,0:,:,:), qz(0:,0:,:,:)
        real(real64), intent(in) :: active(:,:), viscosity(:,:), permeability_x(:,:), permeability_z(:,:), rho_fluid(:,:)
        real(real64), intent(in) :: dt_step

        integer :: i, k
        real(real64) :: rate_x, rate_z, rho_f

        !$omp parallel do collapse(2) private(rate_x,rate_z,rho_f) schedule(static)
        do k = 1, size(qx, 4)
            do i = 1, size(qx, 3)
                if (active(i,k) > 0.5_real64) then
                    rho_f = positive(rho_fluid(i,k), 1.0_real64)
                    if (permeability_x(i,k) > 0.0_real64) then
                        rate_x = positive(viscosity(i,k), 0.0_real64) / (permeability_x(i,k) * rho_f)
                        qx(:,:,i,k) = qx(:,:,i,k) * exp(-rate_x * dt_step)
                    else
                        qx(:,:,i,k) = 0.0_real64
                    endif
                    if (permeability_z(i,k) > 0.0_real64) then
                        rate_z = positive(viscosity(i,k), 0.0_real64) / (permeability_z(i,k) * rho_f)
                        qz(:,:,i,k) = qz(:,:,i,k) * exp(-rate_z * dt_step)
                    else
                        qz(:,:,i,k) = 0.0_real64
                    endif
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine damp_biot_darcy_fields

    subroutine zero_inactive_biot_fields(qx, qz, pressure, active)
        real(real64), intent(inout) :: qx(0:,0:,:,:), qz(0:,0:,:,:), pressure(0:,0:,:,:)
        real(real64), intent(in) :: active(:,:)

        integer :: i, k

        !$omp parallel do collapse(2) schedule(static)
        do k = 1, size(qx, 4)
            do i = 1, size(qx, 3)
                if (active(i,k) <= 0.5_real64) then
                    qx(:,:,i,k) = 0.0_real64
                    qz(:,:,i,k) = 0.0_real64
                    pressure(:,:,i,k) = 0.0_real64
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine zero_inactive_biot_fields

    subroutine zero_inactive_mechanical_fields(vx, vz, sigmaxx, sigmazz, sigmaxz, active)
        real(real64), intent(inout) :: vx(0:,0:,:,:), vz(0:,0:,:,:)
        real(real64), intent(inout) :: sigmaxx(0:,0:,:,:), sigmazz(0:,0:,:,:), sigmaxz(0:,0:,:,:)
        real(real64), intent(in) :: active(:,:)

        integer :: i, k

        !$omp parallel do collapse(2) schedule(static)
        do k = 1, size(vx, 4)
            do i = 1, size(vx, 3)
                if (active(i,k) <= 0.5_real64) then
                    vx(:,:,i,k) = 0.0_real64
                    vz(:,:,i,k) = 0.0_real64
                    sigmaxx(:,:,i,k) = 0.0_real64
                    sigmazz(:,:,i,k) = 0.0_real64
                    sigmaxz(:,:,i,k) = 0.0_real64
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine zero_inactive_mechanical_fields

    subroutine apply_spectral_filter(q)
        real(real64), intent(inout) :: q(0:DG_P, 0:DG_P)
        integer :: i, j
        real(real64) :: filter_weight, alpha, s, Nc

        alpha = 36.0_real64 ! Strength of filter
        Nc = 0.0_real64    ! Cutoff (0 means filter all modes)
        s = 12.0_real64    ! Order of filter

        ! This dampens the highest order modes (the ringing)
        ! while leaving the physical wave (low modes) alone.
        do j = 0, DG_P
            do i = 0, DG_P
                filter_weight = exp(-alpha * ((real(i+j, real64)/real(2*DG_P, real64))**s))
                q(i,j) = q(i,j) * filter_weight
            enddo
        enddo
    end subroutine apply_spectral_filter

end module discontinuous_galerkin_methods
