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

end module discontinuous_galerkin_methods
