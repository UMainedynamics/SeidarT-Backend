module biotdg

    use iso_fortran_env, only: real64
    use seidartio
    use seidart_types
    use constants

    implicit none

    integer, parameter :: DG_P = 4
    integer, parameter :: DG_NP = DG_P + 1

contains

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

        do k = 1, size(cell, 2)
            do i = 1, size(cell, 1)
                do b = 0, DG_P
                    do a = 0, DG_P
                        nodal(a,b,i,k) = cell(i,k)
                    enddo
                enddo
            enddo
        enddo
    end subroutine fill_modal_from_cell

    subroutine cell_average(nodal, weights, cell)
        real(real64), intent(in) :: nodal(0:,0:,:,:)
        real(real64), intent(in) :: weights(0:DG_P)
        real(real64), intent(out) :: cell(:,:)

        integer :: a, b, i, k
        real(real64) :: wsum

        wsum = 4.0_real64
        cell(:,:) = 0.0_real64
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

        do k = 1, nz
            do i = 1, nx
                do b = 0, DG_P
                    if (i > 1) then
                        left_state = u(DG_P,b,i-1,k)
                    else
                        left_state = -u(0,b,i,k)
                    endif
                    jump = u(0,b,i,k) - left_state
                    rhs(0,b,i,k) = rhs(0,b,i,k) - speed * jump * 2.0_real64 / (dx * weights(0))

                    if (i < nx) then
                        right_state = u(0,b,i+1,k)
                    else
                        right_state = -u(DG_P,b,i,k)
                    endif
                    jump = u(DG_P,b,i,k) - right_state
                    rhs(DG_P,b,i,k) = rhs(DG_P,b,i,k) - speed * jump * 2.0_real64 / (dx * weights(DG_P))
                enddo
            enddo
        enddo

        do k = 1, nz
            do i = 1, nx
                do b = 0, DG_P
                    if (k > 1) then
                        left_state = u(b,DG_P,i,k-1)
                    else
                        left_state = -u(b,0,i,k)
                    endif
                    jump = u(b,0,i,k) - left_state
                    rhs(b,0,i,k) = rhs(b,0,i,k) - speed * jump * 2.0_real64 / (dz * weights(0))

                    if (k < nz) then
                        right_state = u(b,0,i,k+1)
                    else
                        right_state = -u(b,DG_P,i,k)
                    endif
                    jump = u(b,DG_P,i,k) - right_state
                    rhs(b,DG_P,i,k) = rhs(b,DG_P,i,k) - speed * jump * 2.0_real64 / (dz * weights(DG_P))
                enddo
            enddo
        enddo
    end subroutine apply_rusanov_penalty

    subroutine rhs_biot2(vx, vz, qx, qz, pressure, sigmaxx, sigmazz, sigmaxz, &
                         c11, c13, c15, c33, c35, c55, rho, rho_fluid, viscosity, &
                         alpha_x, alpha_z, biot_m, permeability_x, permeability_z, &
                         deriv, weights, dx, dz, max_speed, &
                         rvx, rvz, rqx, rqz, rp, rsxx, rszz, rsxz)
        real(real64), intent(in) :: vx(0:,0:,:,:), vz(0:,0:,:,:)
        real(real64), intent(in) :: qx(0:,0:,:,:), qz(0:,0:,:,:)
        real(real64), intent(in) :: pressure(0:,0:,:,:)
        real(real64), intent(in) :: sigmaxx(0:,0:,:,:), sigmazz(0:,0:,:,:)
        real(real64), intent(in) :: sigmaxz(0:,0:,:,:)
        real(real64), intent(in) :: c11(:,:), c13(:,:), c15(:,:), c33(:,:), c35(:,:), c55(:,:)
        real(real64), intent(in) :: rho(:,:), rho_fluid(:,:), viscosity(:,:)
        real(real64), intent(in) :: alpha_x(:,:), alpha_z(:,:), biot_m(:,:)
        real(real64), intent(in) :: permeability_x(:,:), permeability_z(:,:)
        real(real64), intent(in) :: deriv(0:DG_P,0:DG_P), weights(0:DG_P)
        real(real64), intent(in) :: dx, dz, max_speed
        real(real64), intent(out) :: rvx(0:,0:,:,:), rvz(0:,0:,:,:)
        real(real64), intent(out) :: rqx(0:,0:,:,:), rqz(0:,0:,:,:)
        real(real64), intent(out) :: rp(0:,0:,:,:), rsxx(0:,0:,:,:)
        real(real64), intent(out) :: rszz(0:,0:,:,:), rsxz(0:,0:,:,:)

        real(real64), allocatable :: dvx_dx(:,:,:,:), dvx_dz(:,:,:,:), dvz_dx(:,:,:,:), dvz_dz(:,:,:,:)
        real(real64), allocatable :: dqx_dx(:,:,:,:), dqx_dz(:,:,:,:), dqz_dx(:,:,:,:), dqz_dz(:,:,:,:)
        real(real64), allocatable :: dp_dx(:,:,:,:), dp_dz(:,:,:,:)
        real(real64), allocatable :: dsxx_dx(:,:,:,:), dsxx_dz(:,:,:,:), dszz_dx(:,:,:,:), dszz_dz(:,:,:,:)
        real(real64), allocatable :: dsxz_dx(:,:,:,:), dsxz_dz(:,:,:,:)
        integer :: a, b, i, k
        real(real64) :: div_vs, div_q, p_rate, drag_x, drag_z

        allocate(dvx_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)), dvx_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dvz_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)), dvz_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dqx_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)), dqx_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dqz_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)), dqz_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dp_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)), dp_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dsxx_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)), dsxx_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dszz_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)), dszz_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dsxz_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)), dsxz_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))

        call dg_derivatives(vx, deriv, dx, dz, dvx_dx, dvx_dz)
        call dg_derivatives(vz, deriv, dx, dz, dvz_dx, dvz_dz)
        call dg_derivatives(qx, deriv, dx, dz, dqx_dx, dqx_dz)
        call dg_derivatives(qz, deriv, dx, dz, dqz_dx, dqz_dz)
        call dg_derivatives(pressure, deriv, dx, dz, dp_dx, dp_dz)
        call dg_derivatives(sigmaxx, deriv, dx, dz, dsxx_dx, dsxx_dz)
        call dg_derivatives(sigmazz, deriv, dx, dz, dszz_dx, dszz_dz)
        call dg_derivatives(sigmaxz, deriv, dx, dz, dsxz_dx, dsxz_dz)

        do k = 1, size(vx, 4)
            do i = 1, size(vx, 3)
                do b = 0, DG_P
                    do a = 0, DG_P
                        div_vs = dvx_dx(a,b,i,k) + dvz_dz(a,b,i,k)
                        div_q = dqx_dx(a,b,i,k) + dqz_dz(a,b,i,k)
                        p_rate = -biot_m(i,k) * &
                            (0.5_real64 * (alpha_x(i,k) + alpha_z(i,k)) * div_vs + div_q)
                        drag_x = permeability_x(i,k) / &
                            (positive(viscosity(i,k), 1.0e-6_real64) * positive(rho_fluid(i,k), 1.0_real64))
                        drag_z = permeability_z(i,k) / &
                            (positive(viscosity(i,k), 1.0e-6_real64) * positive(rho_fluid(i,k), 1.0_real64))

                        rp(a,b,i,k) = p_rate
                        rsxx(a,b,i,k) = c11(i,k)*dvx_dx(a,b,i,k) + c13(i,k)*dvz_dz(a,b,i,k) + &
                            c15(i,k)*(dvz_dx(a,b,i,k) + dvx_dz(a,b,i,k)) - alpha_x(i,k) * p_rate
                        rszz(a,b,i,k) = c13(i,k)*dvx_dx(a,b,i,k) + c33(i,k)*dvz_dz(a,b,i,k) + &
                            c35(i,k)*(dvz_dx(a,b,i,k) + dvx_dz(a,b,i,k)) - alpha_z(i,k) * p_rate
                        rsxz(a,b,i,k) = c15(i,k)*dvx_dx(a,b,i,k) + c35(i,k)*dvz_dz(a,b,i,k) + &
                            c55(i,k)*(dvz_dx(a,b,i,k) + dvx_dz(a,b,i,k))
                        rvx(a,b,i,k) = (dsxx_dx(a,b,i,k) + dsxz_dz(a,b,i,k)) / positive(rho(i,k), 1.0_real64)
                        rvz(a,b,i,k) = (dsxz_dx(a,b,i,k) + dszz_dz(a,b,i,k)) / positive(rho(i,k), 1.0_real64)
                        rqx(a,b,i,k) = -drag_x * dp_dx(a,b,i,k)
                        rqz(a,b,i,k) = -drag_z * dp_dz(a,b,i,k)
                    enddo
                enddo
            enddo
        enddo

        call apply_rusanov_penalty(vx, rvx, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(vz, rvz, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(qx, rqx, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(qz, rqz, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(pressure, rp, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(sigmaxx, rsxx, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(sigmazz, rszz, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(sigmaxz, rsxz, weights, dx, dz, max_speed)

        deallocate(dvx_dx, dvx_dz, dvz_dx, dvz_dz, dqx_dx, dqx_dz, dqz_dx, dqz_dz)
        deallocate(dp_dx, dp_dz, dsxx_dx, dsxx_dz, dszz_dx, dszz_dz, dsxz_dx, dsxz_dz)
    end subroutine rhs_biot2

    subroutine biotdg2(domain, source, SINGLE_OUTPUT)
        implicit none

        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        integer :: nx, nz
        integer :: it, isource, ksource
        real(real64) :: dx, dz, dt, max_speed
        logical :: SINGLE
        real(real64) :: nodes(0:DG_P), weights(0:DG_P), deriv(0:DG_P,0:DG_P)

        real(real64), allocatable :: c11(:,:), c13(:,:), c15(:,:), c33(:,:), c35(:,:), c55(:,:)
        real(real64), allocatable :: rho(:,:), rho_fluid(:,:), viscosity(:,:)
        real(real64), allocatable :: alpha_x(:,:), alpha_z(:,:), biot_m(:,:)
        real(real64), allocatable :: permeability_x(:,:), permeability_z(:,:)
        real(real64), allocatable :: cell_vx(:,:), cell_vz(:,:), cell_out(:,:)
        real(real64), allocatable :: vx(:,:,:,:), vz(:,:,:,:), qx(:,:,:,:), qz(:,:,:,:), pressure(:,:,:,:)
        real(real64), allocatable :: sigmaxx(:,:,:,:), sigmazz(:,:,:,:), sigmaxz(:,:,:,:)
        real(real64), allocatable :: rvx(:,:,:,:), rvz(:,:,:,:), rqx(:,:,:,:), rqz(:,:,:,:), rp(:,:,:,:)
        real(real64), allocatable :: rsxx(:,:,:,:), rszz(:,:,:,:), rsxz(:,:,:,:)
        real(real64), allocatable :: tvx(:,:,:,:), tvz(:,:,:,:), tqx(:,:,:,:), tqz(:,:,:,:), tp(:,:,:,:)
        real(real64), allocatable :: tsxx(:,:,:,:), tszz(:,:,:,:), tsxz(:,:,:,:)
        real(real64), allocatable :: k1vx(:,:,:,:), k1vz(:,:,:,:), k1qx(:,:,:,:), k1qz(:,:,:,:), k1p(:,:,:,:)
        real(real64), allocatable :: k1sxx(:,:,:,:), k1szz(:,:,:,:), k1sxz(:,:,:,:)
        real(real64), allocatable :: k2vx(:,:,:,:), k2vz(:,:,:,:), k2qx(:,:,:,:), k2qz(:,:,:,:), k2p(:,:,:,:)
        real(real64), allocatable :: k2sxx(:,:,:,:), k2szz(:,:,:,:), k2sxz(:,:,:,:)
        real(real64), allocatable :: k3vx(:,:,:,:), k3vz(:,:,:,:), k3qx(:,:,:,:), k3qz(:,:,:,:), k3p(:,:,:,:)
        real(real64), allocatable :: k3sxx(:,:,:,:), k3szz(:,:,:,:), k3sxz(:,:,:,:)
        real(real64), allocatable :: k4vx(:,:,:,:), k4vz(:,:,:,:), k4qx(:,:,:,:), k4qz(:,:,:,:), k4p(:,:,:,:)
        real(real64), allocatable :: k4sxx(:,:,:,:), k4szz(:,:,:,:), k4sxz(:,:,:,:)
        real(real64), allocatable :: srcx(:), srcz(:), srcxx(:), srcxz(:), srczz(:)
        real(real64), allocatable :: velocnorm(:), pressurenorm(:)

        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        endif

        nx = domain%nx
        nz = domain%nz
        dx = domain%dx
        dz = domain%dz
        dt = source%dt

        call initialize_gll4(nodes, weights, deriv)

        allocate(c11(nx,nz), c13(nx,nz), c15(nx,nz), c33(nx,nz), c35(nx,nz), c55(nx,nz))
        allocate(rho(nx,nz), rho_fluid(nx,nz), viscosity(nx,nz))
        allocate(alpha_x(nx,nz), alpha_z(nx,nz), biot_m(nx,nz))
        allocate(permeability_x(nx,nz), permeability_z(nx,nz))
        allocate(cell_vx(nx,nz), cell_vz(nx,nz), cell_out(nx,nz))

        allocate(vx(0:DG_P,0:DG_P,nx,nz), vz(0:DG_P,0:DG_P,nx,nz))
        allocate(qx(0:DG_P,0:DG_P,nx,nz), qz(0:DG_P,0:DG_P,nx,nz), pressure(0:DG_P,0:DG_P,nx,nz))
        allocate(sigmaxx(0:DG_P,0:DG_P,nx,nz), sigmazz(0:DG_P,0:DG_P,nx,nz), sigmaxz(0:DG_P,0:DG_P,nx,nz))
        allocate(rvx(0:DG_P,0:DG_P,nx,nz), rvz(0:DG_P,0:DG_P,nx,nz), rqx(0:DG_P,0:DG_P,nx,nz))
        allocate(rqz(0:DG_P,0:DG_P,nx,nz), rp(0:DG_P,0:DG_P,nx,nz))
        allocate(rsxx(0:DG_P,0:DG_P,nx,nz), rszz(0:DG_P,0:DG_P,nx,nz), rsxz(0:DG_P,0:DG_P,nx,nz))
        allocate(tvx(0:DG_P,0:DG_P,nx,nz), tvz(0:DG_P,0:DG_P,nx,nz), tqx(0:DG_P,0:DG_P,nx,nz))
        allocate(tqz(0:DG_P,0:DG_P,nx,nz), tp(0:DG_P,0:DG_P,nx,nz))
        allocate(tsxx(0:DG_P,0:DG_P,nx,nz), tszz(0:DG_P,0:DG_P,nx,nz), tsxz(0:DG_P,0:DG_P,nx,nz))

        allocate(k1vx(0:DG_P,0:DG_P,nx,nz), k1vz(0:DG_P,0:DG_P,nx,nz), k1qx(0:DG_P,0:DG_P,nx,nz))
        allocate(k1qz(0:DG_P,0:DG_P,nx,nz), k1p(0:DG_P,0:DG_P,nx,nz))
        allocate(k1sxx(0:DG_P,0:DG_P,nx,nz), k1szz(0:DG_P,0:DG_P,nx,nz), k1sxz(0:DG_P,0:DG_P,nx,nz))
        allocate(k2vx, source=k1vx); allocate(k2vz, source=k1vz); allocate(k2qx, source=k1qx)
        allocate(k2qz, source=k1qz); allocate(k2p, source=k1p)
        allocate(k2sxx, source=k1sxx); allocate(k2szz, source=k1szz); allocate(k2sxz, source=k1sxz)
        allocate(k3vx, source=k1vx); allocate(k3vz, source=k1vz); allocate(k3qx, source=k1qx)
        allocate(k3qz, source=k1qz); allocate(k3p, source=k1p)
        allocate(k3sxx, source=k1sxx); allocate(k3szz, source=k1szz); allocate(k3sxz, source=k1sxz)
        allocate(k4vx, source=k1vx); allocate(k4vz, source=k1vz); allocate(k4qx, source=k1qx)
        allocate(k4qz, source=k1qz); allocate(k4p, source=k1p)
        allocate(k4sxx, source=k1sxx); allocate(k4szz, source=k1szz); allocate(k4sxz, source=k1sxz)

        allocate(srcx(source%time_steps), srcz(source%time_steps))
        allocate(srcxx(source%time_steps), srcxz(source%time_steps), srczz(source%time_steps))
        allocate(velocnorm(source%time_steps), pressurenorm(source%time_steps))

        call material_rw2('c11.dat', c11, .TRUE.)
        call material_rw2('c13.dat', c13, .TRUE.)
        call material_rw2('c15.dat', c15, .TRUE.)
        call material_rw2('c33.dat', c33, .TRUE.)
        call material_rw2('c35.dat', c35, .TRUE.)
        call material_rw2('c55.dat', c55, .TRUE.)
        call material_rw2('rho.dat', rho, .TRUE.)
        call material_rw2('rho_fluid.dat', rho_fluid, .TRUE.)
        call material_rw2('viscosity.dat', viscosity, .TRUE.)
        call material_rw2('alpha_x.dat', alpha_x, .TRUE.)
        call material_rw2('alpha_z.dat', alpha_z, .TRUE.)
        call material_rw2('biot_m.dat', biot_m, .TRUE.)
        call material_rw2('permeability_x.dat', permeability_x, .TRUE.)
        call material_rw2('permeability_z.dat', permeability_z, .TRUE.)
        call material_rw2('initialconditionVx.dat', cell_vx, .TRUE.)
        call material_rw2('initialconditionVz.dat', cell_vz, .TRUE.)

        call fill_modal_from_cell(cell_vx, vx)
        call fill_modal_from_cell(cell_vz, vz)
        qx(:,:,:,:) = 0.0_real64
        qz(:,:,:,:) = 0.0_real64
        pressure(:,:,:,:) = 0.0_real64
        sigmaxx(:,:,:,:) = 0.0_real64
        sigmazz(:,:,:,:) = 0.0_real64
        sigmaxz(:,:,:,:) = 0.0_real64
        srcx(:) = 0.0_real64
        srcz(:) = 0.0_real64
        srcxx(:) = 0.0_real64
        srcxz(:) = 0.0_real64
        srczz(:) = 0.0_real64
        velocnorm(:) = 0.0_real64
        pressurenorm(:) = 0.0_real64

        isource = source%xind + domain%cpml
        ksource = source%zind + domain%cpml
        max_speed = sqrt(maxval(max(c11, c33)) / positive(minval(rho), 1.0_real64))

        if (source%source_type == 'ac') then
            call loadsource('seismicsourcex.dat', source%time_steps, srcx)
            call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        else
            call loadsource('seismicsourcexx.dat', source%time_steps, srcxx)
            call loadsource('seismicsourcexz.dat', source%time_steps, srcxz)
            call loadsource('seismicsourcezz.dat', source%time_steps, srczz)
        endif

        do it = 1, source%time_steps
            call rhs_biot2(vx, vz, qx, qz, pressure, sigmaxx, sigmazz, sigmaxz, &
                           c11, c13, c15, c33, c35, c55, rho, rho_fluid, viscosity, &
                           alpha_x, alpha_z, biot_m, permeability_x, permeability_z, &
                           deriv, weights, dx, dz, max_speed, &
                           k1vx, k1vz, k1qx, k1qz, k1p, k1sxx, k1szz, k1sxz)
            tvx = vx + 0.5_real64 * dt * k1vx
            tvz = vz + 0.5_real64 * dt * k1vz
            tqx = qx + 0.5_real64 * dt * k1qx
            tqz = qz + 0.5_real64 * dt * k1qz
            tp = pressure + 0.5_real64 * dt * k1p
            tsxx = sigmaxx + 0.5_real64 * dt * k1sxx
            tszz = sigmazz + 0.5_real64 * dt * k1szz
            tsxz = sigmaxz + 0.5_real64 * dt * k1sxz

            call rhs_biot2(tvx, tvz, tqx, tqz, tp, tsxx, tszz, tsxz, &
                           c11, c13, c15, c33, c35, c55, rho, rho_fluid, viscosity, &
                           alpha_x, alpha_z, biot_m, permeability_x, permeability_z, &
                           deriv, weights, dx, dz, max_speed, &
                           k2vx, k2vz, k2qx, k2qz, k2p, k2sxx, k2szz, k2sxz)
            tvx = vx + 0.5_real64 * dt * k2vx
            tvz = vz + 0.5_real64 * dt * k2vz
            tqx = qx + 0.5_real64 * dt * k2qx
            tqz = qz + 0.5_real64 * dt * k2qz
            tp = pressure + 0.5_real64 * dt * k2p
            tsxx = sigmaxx + 0.5_real64 * dt * k2sxx
            tszz = sigmazz + 0.5_real64 * dt * k2szz
            tsxz = sigmaxz + 0.5_real64 * dt * k2sxz

            call rhs_biot2(tvx, tvz, tqx, tqz, tp, tsxx, tszz, tsxz, &
                           c11, c13, c15, c33, c35, c55, rho, rho_fluid, viscosity, &
                           alpha_x, alpha_z, biot_m, permeability_x, permeability_z, &
                           deriv, weights, dx, dz, max_speed, &
                           k3vx, k3vz, k3qx, k3qz, k3p, k3sxx, k3szz, k3sxz)
            tvx = vx + dt * k3vx
            tvz = vz + dt * k3vz
            tqx = qx + dt * k3qx
            tqz = qz + dt * k3qz
            tp = pressure + dt * k3p
            tsxx = sigmaxx + dt * k3sxx
            tszz = sigmazz + dt * k3szz
            tsxz = sigmaxz + dt * k3sxz

            call rhs_biot2(tvx, tvz, tqx, tqz, tp, tsxx, tszz, tsxz, &
                           c11, c13, c15, c33, c35, c55, rho, rho_fluid, viscosity, &
                           alpha_x, alpha_z, biot_m, permeability_x, permeability_z, &
                           deriv, weights, dx, dz, max_speed, &
                           k4vx, k4vz, k4qx, k4qz, k4p, k4sxx, k4szz, k4sxz)

            vx = vx + dt * (k1vx + 2.0_real64*k2vx + 2.0_real64*k3vx + k4vx) / 6.0_real64
            vz = vz + dt * (k1vz + 2.0_real64*k2vz + 2.0_real64*k3vz + k4vz) / 6.0_real64
            qx = qx + dt * (k1qx + 2.0_real64*k2qx + 2.0_real64*k3qx + k4qx) / 6.0_real64
            qz = qz + dt * (k1qz + 2.0_real64*k2qz + 2.0_real64*k3qz + k4qz) / 6.0_real64
            pressure = pressure + dt * (k1p + 2.0_real64*k2p + 2.0_real64*k3p + k4p) / 6.0_real64
            sigmaxx = sigmaxx + dt * (k1sxx + 2.0_real64*k2sxx + 2.0_real64*k3sxx + k4sxx) / 6.0_real64
            sigmazz = sigmazz + dt * (k1szz + 2.0_real64*k2szz + 2.0_real64*k3szz + k4szz) / 6.0_real64
            sigmaxz = sigmaxz + dt * (k1sxz + 2.0_real64*k2sxz + 2.0_real64*k3sxz + k4sxz) / 6.0_real64

            vx(0,0,isource,ksource) = vx(0,0,isource,ksource) + srcx(it) * dt / positive(rho(isource,ksource), 1.0_real64)
            vz(0,0,isource,ksource) = vz(0,0,isource,ksource) + srcz(it) * dt / positive(rho(isource,ksource), 1.0_real64)
            sigmaxx(0,0,isource,ksource) = sigmaxx(0,0,isource,ksource) + srcxx(it) / positive(rho(isource,ksource), 1.0_real64)
            sigmazz(0,0,isource,ksource) = sigmazz(0,0,isource,ksource) + srczz(it) / positive(rho(isource,ksource), 1.0_real64)
            sigmaxz(0,0,isource,ksource) = sigmaxz(0,0,isource,ksource) + srcxz(it) / positive(rho(isource,ksource), 1.0_real64)

            call cell_average(vx, weights, cell_out)
            call write_image2(cell_out, nx, nz, source, it, 'Vx', SINGLE)
            velocnorm(it) = maxval(abs(cell_out))
            call cell_average(vz, weights, cell_out)
            call write_image2(cell_out, nx, nz, source, it, 'Vz', SINGLE)
            velocnorm(it) = max(velocnorm(it), maxval(abs(cell_out)))
            call cell_average(qx, weights, cell_out)
            call write_image2(cell_out, nx, nz, source, it, 'Qx', SINGLE)
            call cell_average(qz, weights, cell_out)
            call write_image2(cell_out, nx, nz, source, it, 'Qz', SINGLE)
            call cell_average(pressure, weights, cell_out)
            call write_image2(cell_out, nx, nz, source, it, 'Pp', SINGLE)
            pressurenorm(it) = maxval(abs(cell_out))

            if (velocnorm(it) > stability_threshold) stop 'Biot DG 2D solution became unstable'
        enddo

        call write_array('velocity_norm.dat', source%time_steps, velocnorm)
        call write_array('pore_pressure_norm.dat', source%time_steps, pressurenorm)
    end subroutine biotdg2

    subroutine biotdg25(domain, source, SINGLE_OUTPUT)
        implicit none

        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        integer :: nx, ny, nz, it, isource, jsource, ksource
        real(real64) :: dx, dy, dz, dt, max_speed, rho0
        logical :: SINGLE

        real(real64), allocatable :: c11(:,:), c12(:,:), c13(:,:), c22(:,:), c23(:,:), c33(:,:), rho(:,:)
        real(real64), allocatable :: alpha_x(:,:), alpha_y(:,:), alpha_z(:,:), biot_m(:,:)
        real(real64), allocatable :: rho_fluid(:,:), viscosity(:,:), permeability_x(:,:), permeability_y(:,:), permeability_z(:,:)
        real(real64), allocatable :: vx(:,:,:), vy(:,:,:), vz(:,:,:), qx(:,:,:), qy(:,:,:), qz(:,:,:), pressure(:,:,:)
        real(real64), allocatable :: sxx(:,:,:), syy(:,:,:), szz(:,:,:)
        real(real64), allocatable :: k1vx(:,:,:), k1vy(:,:,:), k1vz(:,:,:), k1qx(:,:,:), k1qy(:,:,:), k1qz(:,:,:), k1p(:,:,:)
        real(real64), allocatable :: k1sxx(:,:,:), k1syy(:,:,:), k1szz(:,:,:)
        real(real64), allocatable :: k2vx(:,:,:), k2vy(:,:,:), k2vz(:,:,:), k2qx(:,:,:), k2qy(:,:,:), k2qz(:,:,:), k2p(:,:,:)
        real(real64), allocatable :: k2sxx(:,:,:), k2syy(:,:,:), k2szz(:,:,:)
        real(real64), allocatable :: k3vx(:,:,:), k3vy(:,:,:), k3vz(:,:,:), k3qx(:,:,:), k3qy(:,:,:), k3qz(:,:,:), k3p(:,:,:)
        real(real64), allocatable :: k3sxx(:,:,:), k3syy(:,:,:), k3szz(:,:,:)
        real(real64), allocatable :: k4vx(:,:,:), k4vy(:,:,:), k4vz(:,:,:), k4qx(:,:,:), k4qy(:,:,:), k4qz(:,:,:), k4p(:,:,:)
        real(real64), allocatable :: k4sxx(:,:,:), k4syy(:,:,:), k4szz(:,:,:)
        real(real64), allocatable :: tvx(:,:,:), tvy(:,:,:), tvz(:,:,:), tqx(:,:,:), tqy(:,:,:), tqz(:,:,:), tp(:,:,:)
        real(real64), allocatable :: tsxx(:,:,:), tsyy(:,:,:), tszz(:,:,:)
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:), srcxx(:), srcyy(:), srczz(:)
        real(real64), allocatable :: velocnorm(:), pressurenorm(:)

        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        endif

        nx = domain%nx; ny = domain%ny; nz = domain%nz
        dx = domain%dx; dy = domain%dy; dz = domain%dz; dt = source%dt

        allocate(c11(nx,nz), c12(nx,nz), c13(nx,nz), c22(nx,nz), c23(nx,nz), c33(nx,nz), rho(nx,nz))
        allocate(alpha_x(nx,nz), alpha_y(nx,nz), alpha_z(nx,nz), biot_m(nx,nz))
        allocate(rho_fluid(nx,nz), viscosity(nx,nz), permeability_x(nx,nz), permeability_y(nx,nz), permeability_z(nx,nz))
        allocate(vx(nx,ny,nz), vy(nx,ny,nz), vz(nx,ny,nz), qx(nx,ny,nz), qy(nx,ny,nz), qz(nx,ny,nz), pressure(nx,ny,nz))
        allocate(sxx(nx,ny,nz), syy(nx,ny,nz), szz(nx,ny,nz))
        allocate(k1vx(nx,ny,nz), k1vy(nx,ny,nz), k1vz(nx,ny,nz), k1qx(nx,ny,nz), k1qy(nx,ny,nz), k1qz(nx,ny,nz), k1p(nx,ny,nz))
        allocate(k1sxx(nx,ny,nz), k1syy(nx,ny,nz), k1szz(nx,ny,nz))
        allocate(k2vx, source=k1vx); allocate(k2vy, source=k1vy); allocate(k2vz, source=k1vz)
        allocate(k2qx, source=k1qx); allocate(k2qy, source=k1qy); allocate(k2qz, source=k1qz); allocate(k2p, source=k1p)
        allocate(k2sxx, source=k1sxx); allocate(k2syy, source=k1syy); allocate(k2szz, source=k1szz)
        allocate(k3vx, source=k1vx); allocate(k3vy, source=k1vy); allocate(k3vz, source=k1vz)
        allocate(k3qx, source=k1qx); allocate(k3qy, source=k1qy); allocate(k3qz, source=k1qz); allocate(k3p, source=k1p)
        allocate(k3sxx, source=k1sxx); allocate(k3syy, source=k1syy); allocate(k3szz, source=k1szz)
        allocate(k4vx, source=k1vx); allocate(k4vy, source=k1vy); allocate(k4vz, source=k1vz)
        allocate(k4qx, source=k1qx); allocate(k4qy, source=k1qy); allocate(k4qz, source=k1qz); allocate(k4p, source=k1p)
        allocate(k4sxx, source=k1sxx); allocate(k4syy, source=k1syy); allocate(k4szz, source=k1szz)
        allocate(tvx, source=k1vx); allocate(tvy, source=k1vy); allocate(tvz, source=k1vz)
        allocate(tqx, source=k1qx); allocate(tqy, source=k1qy); allocate(tqz, source=k1qz); allocate(tp, source=k1p)
        allocate(tsxx, source=k1sxx); allocate(tsyy, source=k1syy); allocate(tszz, source=k1szz)
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps))
        allocate(srcxx(source%time_steps), srcyy(source%time_steps), srczz(source%time_steps))
        allocate(velocnorm(source%time_steps), pressurenorm(source%time_steps))

        call material_rw2('c11.dat', c11, .TRUE.); call material_rw2('c12.dat', c12, .TRUE.)
        call material_rw2('c13.dat', c13, .TRUE.); call material_rw2('c22.dat', c22, .TRUE.)
        call material_rw2('c23.dat', c23, .TRUE.); call material_rw2('c33.dat', c33, .TRUE.)
        call material_rw2('rho.dat', rho, .TRUE.); call material_rw2('alpha_x.dat', alpha_x, .TRUE.)
        call material_rw2('alpha_y.dat', alpha_y, .TRUE.); call material_rw2('alpha_z.dat', alpha_z, .TRUE.)
        call material_rw2('biot_m.dat', biot_m, .TRUE.); call material_rw2('rho_fluid.dat', rho_fluid, .TRUE.)
        call material_rw2('viscosity.dat', viscosity, .TRUE.); call material_rw2('permeability_x.dat', permeability_x, .TRUE.)
        call material_rw2('permeability_y.dat', permeability_y, .TRUE.); call material_rw2('permeability_z.dat', permeability_z, .TRUE.)
        call material_rw3('initialconditionVx.dat', vx, .TRUE.)
        call material_rw3('initialconditionVy.dat', vy, .TRUE.)
        call material_rw3('initialconditionVz.dat', vz, .TRUE.)

        qx = 0.0_real64; qy = 0.0_real64; qz = 0.0_real64; pressure = 0.0_real64
        sxx = 0.0_real64; syy = 0.0_real64; szz = 0.0_real64
        srcx = 0.0_real64; srcy = 0.0_real64; srcz = 0.0_real64
        srcxx = 0.0_real64; srcyy = 0.0_real64; srczz = 0.0_real64
        velocnorm = 0.0_real64; pressurenorm = 0.0_real64

        isource = source%xind + domain%cpml
        jsource = source%yind + domain%cpml
        ksource = source%zind + domain%cpml
        max_speed = sqrt(maxval(max(c11, c33)) / positive(minval(rho), 1.0_real64))

        if (source%source_type == 'ac') then
            call loadsource('seismicsourcex.dat', source%time_steps, srcx)
            call loadsource('seismicsourcey.dat', source%time_steps, srcy)
            call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        else
            call loadsource('seismicsourcexx.dat', source%time_steps, srcxx)
            call loadsource('seismicsourceyy.dat', source%time_steps, srcyy)
            call loadsource('seismicsourcezz.dat', source%time_steps, srczz)
        endif

        do it = 1, source%time_steps
            call rhs_biot25_cell(vx,vy,vz,qx,qy,qz,pressure,sxx,syy,szz,c11,c12,c13,c22,c23,c33,rho, &
                rho_fluid,viscosity,alpha_x,alpha_y,alpha_z,biot_m,permeability_x,permeability_y,permeability_z, &
                dx,dy,dz,max_speed,k1vx,k1vy,k1vz,k1qx,k1qy,k1qz,k1p,k1sxx,k1syy,k1szz)
            tvx=vx+0.5_real64*dt*k1vx; tvy=vy+0.5_real64*dt*k1vy; tvz=vz+0.5_real64*dt*k1vz
            tqx=qx+0.5_real64*dt*k1qx; tqy=qy+0.5_real64*dt*k1qy; tqz=qz+0.5_real64*dt*k1qz
            tp=pressure+0.5_real64*dt*k1p; tsxx=sxx+0.5_real64*dt*k1sxx; tsyy=syy+0.5_real64*dt*k1syy; tszz=szz+0.5_real64*dt*k1szz
            call rhs_biot25_cell(tvx,tvy,tvz,tqx,tqy,tqz,tp,tsxx,tsyy,tszz,c11,c12,c13,c22,c23,c33,rho, &
                rho_fluid,viscosity,alpha_x,alpha_y,alpha_z,biot_m,permeability_x,permeability_y,permeability_z, &
                dx,dy,dz,max_speed,k2vx,k2vy,k2vz,k2qx,k2qy,k2qz,k2p,k2sxx,k2syy,k2szz)
            tvx=vx+0.5_real64*dt*k2vx; tvy=vy+0.5_real64*dt*k2vy; tvz=vz+0.5_real64*dt*k2vz
            tqx=qx+0.5_real64*dt*k2qx; tqy=qy+0.5_real64*dt*k2qy; tqz=qz+0.5_real64*dt*k2qz
            tp=pressure+0.5_real64*dt*k2p; tsxx=sxx+0.5_real64*dt*k2sxx; tsyy=syy+0.5_real64*dt*k2syy; tszz=szz+0.5_real64*dt*k2szz
            call rhs_biot25_cell(tvx,tvy,tvz,tqx,tqy,tqz,tp,tsxx,tsyy,tszz,c11,c12,c13,c22,c23,c33,rho, &
                rho_fluid,viscosity,alpha_x,alpha_y,alpha_z,biot_m,permeability_x,permeability_y,permeability_z, &
                dx,dy,dz,max_speed,k3vx,k3vy,k3vz,k3qx,k3qy,k3qz,k3p,k3sxx,k3syy,k3szz)
            tvx=vx+dt*k3vx; tvy=vy+dt*k3vy; tvz=vz+dt*k3vz; tqx=qx+dt*k3qx; tqy=qy+dt*k3qy; tqz=qz+dt*k3qz
            tp=pressure+dt*k3p; tsxx=sxx+dt*k3sxx; tsyy=syy+dt*k3syy; tszz=szz+dt*k3szz
            call rhs_biot25_cell(tvx,tvy,tvz,tqx,tqy,tqz,tp,tsxx,tsyy,tszz,c11,c12,c13,c22,c23,c33,rho, &
                rho_fluid,viscosity,alpha_x,alpha_y,alpha_z,biot_m,permeability_x,permeability_y,permeability_z, &
                dx,dy,dz,max_speed,k4vx,k4vy,k4vz,k4qx,k4qy,k4qz,k4p,k4sxx,k4syy,k4szz)

            vx=vx+dt*(k1vx+2.0_real64*k2vx+2.0_real64*k3vx+k4vx)/6.0_real64
            vy=vy+dt*(k1vy+2.0_real64*k2vy+2.0_real64*k3vy+k4vy)/6.0_real64
            vz=vz+dt*(k1vz+2.0_real64*k2vz+2.0_real64*k3vz+k4vz)/6.0_real64
            qx=qx+dt*(k1qx+2.0_real64*k2qx+2.0_real64*k3qx+k4qx)/6.0_real64
            qy=qy+dt*(k1qy+2.0_real64*k2qy+2.0_real64*k3qy+k4qy)/6.0_real64
            qz=qz+dt*(k1qz+2.0_real64*k2qz+2.0_real64*k3qz+k4qz)/6.0_real64
            pressure=pressure+dt*(k1p+2.0_real64*k2p+2.0_real64*k3p+k4p)/6.0_real64
            sxx=sxx+dt*(k1sxx+2.0_real64*k2sxx+2.0_real64*k3sxx+k4sxx)/6.0_real64
            syy=syy+dt*(k1syy+2.0_real64*k2syy+2.0_real64*k3syy+k4syy)/6.0_real64
            szz=szz+dt*(k1szz+2.0_real64*k2szz+2.0_real64*k3szz+k4szz)/6.0_real64

            rho0 = positive(rho(isource,ksource), 1.0_real64)
            vx(isource,jsource,ksource)=vx(isource,jsource,ksource)+srcx(it)*dt/rho0
            vy(isource,jsource,ksource)=vy(isource,jsource,ksource)+srcy(it)*dt/rho0
            vz(isource,jsource,ksource)=vz(isource,jsource,ksource)+srcz(it)*dt/rho0
            sxx(isource,jsource,ksource)=sxx(isource,jsource,ksource)+srcxx(it)/rho0
            syy(isource,jsource,ksource)=syy(isource,jsource,ksource)+srcyy(it)/rho0
            szz(isource,jsource,ksource)=szz(isource,jsource,ksource)+srczz(it)/rho0

            vx(1,:,:)=0.0_real64; vy(:,1,:)=0.0_real64; vz(:,:,1)=0.0_real64
            vx(nx,:,:)=0.0_real64; vy(:,ny,:)=0.0_real64; vz(:,:,nz)=0.0_real64

            velocnorm(it)=maxval(sqrt(vx**2+vy**2+vz**2)); pressurenorm(it)=maxval(abs(pressure))
            if (velocnorm(it) > stability_threshold) stop 'Biot DG 2.5D solution became unstable'
            call write_image3(vx,nx,ny,nz,source,it,'Vx',SINGLE); call write_image3(vy,nx,ny,nz,source,it,'Vy',SINGLE)
            call write_image3(vz,nx,ny,nz,source,it,'Vz',SINGLE); call write_image3(qx,nx,ny,nz,source,it,'Qx',SINGLE)
            call write_image3(qy,nx,ny,nz,source,it,'Qy',SINGLE); call write_image3(qz,nx,ny,nz,source,it,'Qz',SINGLE)
            call write_image3(pressure,nx,ny,nz,source,it,'Pp',SINGLE)
        enddo

        call write_array('velocity_norm.dat', source%time_steps, velocnorm)
        call write_array('pore_pressure_norm.dat', source%time_steps, pressurenorm)
    end subroutine biotdg25

    subroutine biotdg3(domain, source, SINGLE_OUTPUT)
        implicit none

        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        integer :: nx, ny, nz, i, j, k, it, isource, jsource, ksource
        real(real64) :: dx, dy, dz, dt, rho0, div_v, div_q, px, py, pz
        logical :: SINGLE

        real(real64), allocatable :: c11(:,:,:), c12(:,:,:), c13(:,:,:), c22(:,:,:), c23(:,:,:), c33(:,:,:), rho(:,:,:)
        real(real64), allocatable :: ax(:,:,:), ay(:,:,:), az(:,:,:), m(:,:,:), rho_f(:,:,:), eta(:,:,:), kx(:,:,:), ky(:,:,:), kz(:,:,:)
        real(real64), allocatable :: vx(:,:,:), vy(:,:,:), vz(:,:,:), qx(:,:,:), qy(:,:,:), qz(:,:,:), p(:,:,:)
        real(real64), allocatable :: sxx(:,:,:), syy(:,:,:), szz(:,:,:)
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:), srcxx(:), srcyy(:), srczz(:)
        real(real64), allocatable :: velocnorm(:), pressurenorm(:)

        if (present(SINGLE_OUTPUT)) then
            SINGLE = SINGLE_OUTPUT
        else
            SINGLE = .TRUE.
        endif

        nx=domain%nx; ny=domain%ny; nz=domain%nz
        dx=domain%dx; dy=domain%dy; dz=domain%dz; dt=source%dt

        allocate(c11(nx,ny,nz), c12(nx,ny,nz), c13(nx,ny,nz), c22(nx,ny,nz), c23(nx,ny,nz), c33(nx,ny,nz), rho(nx,ny,nz))
        allocate(ax(nx,ny,nz), ay(nx,ny,nz), az(nx,ny,nz), m(nx,ny,nz), rho_f(nx,ny,nz), eta(nx,ny,nz))
        allocate(kx(nx,ny,nz), ky(nx,ny,nz), kz(nx,ny,nz))
        allocate(vx(nx,ny,nz), vy(nx,ny,nz), vz(nx,ny,nz), qx(nx,ny,nz), qy(nx,ny,nz), qz(nx,ny,nz), p(nx,ny,nz))
        allocate(sxx(nx,ny,nz), syy(nx,ny,nz), szz(nx,ny,nz))
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps))
        allocate(srcxx(source%time_steps), srcyy(source%time_steps), srczz(source%time_steps))
        allocate(velocnorm(source%time_steps), pressurenorm(source%time_steps))

        call material_rw3('c11.dat', c11, .TRUE.); call material_rw3('c12.dat', c12, .TRUE.)
        call material_rw3('c13.dat', c13, .TRUE.); call material_rw3('c22.dat', c22, .TRUE.)
        call material_rw3('c23.dat', c23, .TRUE.); call material_rw3('c33.dat', c33, .TRUE.)
        call material_rw3('rho.dat', rho, .TRUE.); call material_rw3('alpha_x.dat', ax, .TRUE.)
        call material_rw3('alpha_y.dat', ay, .TRUE.); call material_rw3('alpha_z.dat', az, .TRUE.)
        call material_rw3('biot_m.dat', m, .TRUE.); call material_rw3('rho_fluid.dat', rho_f, .TRUE.)
        call material_rw3('viscosity.dat', eta, .TRUE.); call material_rw3('permeability_x.dat', kx, .TRUE.)
        call material_rw3('permeability_y.dat', ky, .TRUE.); call material_rw3('permeability_z.dat', kz, .TRUE.)
        call material_rw3('initialconditionVx.dat', vx, .TRUE.)
        call material_rw3('initialconditionVy.dat', vy, .TRUE.)
        call material_rw3('initialconditionVz.dat', vz, .TRUE.)

        qx=0.0_real64; qy=0.0_real64; qz=0.0_real64; p=0.0_real64
        sxx=0.0_real64; syy=0.0_real64; szz=0.0_real64
        srcx=0.0_real64; srcy=0.0_real64; srcz=0.0_real64
        srcxx=0.0_real64; srcyy=0.0_real64; srczz=0.0_real64
        velocnorm=0.0_real64; pressurenorm=0.0_real64

        isource=source%xind+domain%cpml; jsource=source%yind+domain%cpml; ksource=source%zind+domain%cpml

        if (source%source_type == 'ac') then
            call loadsource('seismicsourcex.dat', source%time_steps, srcx)
            call loadsource('seismicsourcey.dat', source%time_steps, srcy)
            call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        else
            call loadsource('seismicsourcexx.dat', source%time_steps, srcxx)
            call loadsource('seismicsourceyy.dat', source%time_steps, srcyy)
            call loadsource('seismicsourcezz.dat', source%time_steps, srczz)
        endif

        do it=1,source%time_steps
            do k=2,nz-1
                do j=2,ny-1
                    do i=2,nx-1
                        div_v=(vx(i+1,j,k)-vx(i-1,j,k))/(2.0_real64*dx)+ &
                              (vy(i,j+1,k)-vy(i,j-1,k))/(2.0_real64*dy)+ &
                              (vz(i,j,k+1)-vz(i,j,k-1))/(2.0_real64*dz)
                        div_q=(qx(i+1,j,k)-qx(i-1,j,k))/(2.0_real64*dx)+ &
                              (qy(i,j+1,k)-qy(i,j-1,k))/(2.0_real64*dy)+ &
                              (qz(i,j,k+1)-qz(i,j,k-1))/(2.0_real64*dz)
                        p(i,j,k)=p(i,j,k)-dt*m(i,j,k)*(((ax(i,j,k)+ay(i,j,k)+az(i,j,k))/3.0_real64)*div_v+div_q)

                        sxx(i,j,k)=sxx(i,j,k)+dt*(c11(i,j,k)*(vx(i+1,j,k)-vx(i-1,j,k))/(2.0_real64*dx)+ &
                            c12(i,j,k)*(vy(i,j+1,k)-vy(i,j-1,k))/(2.0_real64*dy)+ &
                            c13(i,j,k)*(vz(i,j,k+1)-vz(i,j,k-1))/(2.0_real64*dz)-ax(i,j,k)*p(i,j,k))
                        syy(i,j,k)=syy(i,j,k)+dt*(c12(i,j,k)*(vx(i+1,j,k)-vx(i-1,j,k))/(2.0_real64*dx)+ &
                            c22(i,j,k)*(vy(i,j+1,k)-vy(i,j-1,k))/(2.0_real64*dy)+ &
                            c23(i,j,k)*(vz(i,j,k+1)-vz(i,j,k-1))/(2.0_real64*dz)-ay(i,j,k)*p(i,j,k))
                        szz(i,j,k)=szz(i,j,k)+dt*(c13(i,j,k)*(vx(i+1,j,k)-vx(i-1,j,k))/(2.0_real64*dx)+ &
                            c23(i,j,k)*(vy(i,j+1,k)-vy(i,j-1,k))/(2.0_real64*dy)+ &
                            c33(i,j,k)*(vz(i,j,k+1)-vz(i,j,k-1))/(2.0_real64*dz)-az(i,j,k)*p(i,j,k))
                    enddo
                enddo
            enddo

            do k=2,nz-1
                do j=2,ny-1
                    do i=2,nx-1
                        px=(p(i+1,j,k)-p(i-1,j,k))/(2.0_real64*dx)
                        py=(p(i,j+1,k)-p(i,j-1,k))/(2.0_real64*dy)
                        pz=(p(i,j,k+1)-p(i,j,k-1))/(2.0_real64*dz)
                        vx(i,j,k)=vx(i,j,k)+dt*(sxx(i+1,j,k)-sxx(i-1,j,k))/(2.0_real64*dx)/positive(rho(i,j,k),1.0_real64)
                        vy(i,j,k)=vy(i,j,k)+dt*(syy(i,j+1,k)-syy(i,j-1,k))/(2.0_real64*dy)/positive(rho(i,j,k),1.0_real64)
                        vz(i,j,k)=vz(i,j,k)+dt*(szz(i,j,k+1)-szz(i,j,k-1))/(2.0_real64*dz)/positive(rho(i,j,k),1.0_real64)
                        qx(i,j,k)=qx(i,j,k)-dt*kx(i,j,k)*px/(positive(eta(i,j,k),1.0e-6_real64)*positive(rho_f(i,j,k),1.0_real64))
                        qy(i,j,k)=qy(i,j,k)-dt*ky(i,j,k)*py/(positive(eta(i,j,k),1.0e-6_real64)*positive(rho_f(i,j,k),1.0_real64))
                        qz(i,j,k)=qz(i,j,k)-dt*kz(i,j,k)*pz/(positive(eta(i,j,k),1.0e-6_real64)*positive(rho_f(i,j,k),1.0_real64))
                    enddo
                enddo
            enddo

            rho0=positive(rho(isource,jsource,ksource),1.0_real64)
            vx(isource,jsource,ksource)=vx(isource,jsource,ksource)+srcx(it)*dt/rho0
            vy(isource,jsource,ksource)=vy(isource,jsource,ksource)+srcy(it)*dt/rho0
            vz(isource,jsource,ksource)=vz(isource,jsource,ksource)+srcz(it)*dt/rho0
            sxx(isource,jsource,ksource)=sxx(isource,jsource,ksource)+srcxx(it)/rho0
            syy(isource,jsource,ksource)=syy(isource,jsource,ksource)+srcyy(it)/rho0
            szz(isource,jsource,ksource)=szz(isource,jsource,ksource)+srczz(it)/rho0

            vx(1,:,:)=0.0_real64; vy(:,1,:)=0.0_real64; vz(:,:,1)=0.0_real64
            vx(nx,:,:)=0.0_real64; vy(:,ny,:)=0.0_real64; vz(:,:,nz)=0.0_real64
            velocnorm(it)=maxval(sqrt(vx**2+vy**2+vz**2)); pressurenorm(it)=maxval(abs(p))
            if (velocnorm(it) > stability_threshold) stop 'Biot DG 3D solution became unstable'

            call write_image3(vx,nx,ny,nz,source,it,'Vx',SINGLE); call write_image3(vy,nx,ny,nz,source,it,'Vy',SINGLE)
            call write_image3(vz,nx,ny,nz,source,it,'Vz',SINGLE); call write_image3(qx,nx,ny,nz,source,it,'Qx',SINGLE)
            call write_image3(qy,nx,ny,nz,source,it,'Qy',SINGLE); call write_image3(qz,nx,ny,nz,source,it,'Qz',SINGLE)
            call write_image3(p,nx,ny,nz,source,it,'Pp',SINGLE)
        enddo

        call write_array('velocity_norm.dat', source%time_steps, velocnorm)
        call write_array('pore_pressure_norm.dat', source%time_steps, pressurenorm)
    end subroutine biotdg3

    subroutine rhs_biot25_cell(vx,vy,vz,qx,qy,qz,p,sxx,syy,szz,c11,c12,c13,c22,c23,c33,rho, &
                               rho_f,eta,ax,ay,az,m,kx,ky,kz,dx,dy,dz,speed, &
                               rvx,rvy,rvz,rqx,rqy,rqz,rp,rsxx,rsyy,rszz)
        implicit none

        real(real64), intent(in) :: vx(:,:,:), vy(:,:,:), vz(:,:,:), qx(:,:,:), qy(:,:,:), qz(:,:,:)
        real(real64), intent(in) :: p(:,:,:), sxx(:,:,:), syy(:,:,:), szz(:,:,:)
        real(real64), intent(in) :: c11(:,:), c12(:,:), c13(:,:), c22(:,:), c23(:,:), c33(:,:), rho(:,:)
        real(real64), intent(in) :: rho_f(:,:), eta(:,:), ax(:,:), ay(:,:), az(:,:), m(:,:), kx(:,:), ky(:,:), kz(:,:)
        real(real64), intent(in) :: dx, dy, dz, speed
        real(real64), intent(out) :: rvx(:,:,:), rvy(:,:,:), rvz(:,:,:), rqx(:,:,:), rqy(:,:,:), rqz(:,:,:)
        real(real64), intent(out) :: rp(:,:,:), rsxx(:,:,:), rsyy(:,:,:), rszz(:,:,:)

        integer :: i, j, k, nx, ny, nz
        real(real64) :: div_v, div_q, px, py, pz, dampx, dampy, dampz

        nx=size(vx,1); ny=size(vx,2); nz=size(vx,3)
        rvx=0.0_real64; rvy=0.0_real64; rvz=0.0_real64; rqx=0.0_real64; rqy=0.0_real64; rqz=0.0_real64
        rp=0.0_real64; rsxx=0.0_real64; rsyy=0.0_real64; rszz=0.0_real64

        do k=2,nz-1
            do j=2,ny-1
                do i=2,nx-1
                    div_v = (vx(i+1,j,k)-vx(i-1,j,k))/(2.0_real64*dx) + &
                            (vy(i,j+1,k)-vy(i,j-1,k))/(2.0_real64*dy) + &
                            (vz(i,j,k+1)-vz(i,j,k-1))/(2.0_real64*dz)
                    div_q = (qx(i+1,j,k)-qx(i-1,j,k))/(2.0_real64*dx) + &
                            (qy(i,j+1,k)-qy(i,j-1,k))/(2.0_real64*dy) + &
                            (qz(i,j,k+1)-qz(i,j,k-1))/(2.0_real64*dz)
                    rp(i,j,k) = -m(i,k) * (((ax(i,k)+ay(i,k)+az(i,k))/3.0_real64)*div_v + div_q)

                    rsxx(i,j,k) = c11(i,k)*(vx(i+1,j,k)-vx(i-1,j,k))/(2.0_real64*dx) + &
                                  c12(i,k)*(vy(i,j+1,k)-vy(i,j-1,k))/(2.0_real64*dy) + &
                                  c13(i,k)*(vz(i,j,k+1)-vz(i,j,k-1))/(2.0_real64*dz) - ax(i,k)*rp(i,j,k)
                    rsyy(i,j,k) = c12(i,k)*(vx(i+1,j,k)-vx(i-1,j,k))/(2.0_real64*dx) + &
                                  c22(i,k)*(vy(i,j+1,k)-vy(i,j-1,k))/(2.0_real64*dy) + &
                                  c23(i,k)*(vz(i,j,k+1)-vz(i,j,k-1))/(2.0_real64*dz) - ay(i,k)*rp(i,j,k)
                    rszz(i,j,k) = c13(i,k)*(vx(i+1,j,k)-vx(i-1,j,k))/(2.0_real64*dx) + &
                                  c23(i,k)*(vy(i,j+1,k)-vy(i,j-1,k))/(2.0_real64*dy) + &
                                  c33(i,k)*(vz(i,j,k+1)-vz(i,j,k-1))/(2.0_real64*dz) - az(i,k)*rp(i,j,k)

                    px=(p(i+1,j,k)-p(i-1,j,k))/(2.0_real64*dx); py=(p(i,j+1,k)-p(i,j-1,k))/(2.0_real64*dy)
                    pz=(p(i,j,k+1)-p(i,j,k-1))/(2.0_real64*dz)
                    rvx(i,j,k)=(sxx(i+1,j,k)-sxx(i-1,j,k))/(2.0_real64*dx)/positive(rho(i,k),1.0_real64)
                    rvy(i,j,k)=(syy(i,j+1,k)-syy(i,j-1,k))/(2.0_real64*dy)/positive(rho(i,k),1.0_real64)
                    rvz(i,j,k)=(szz(i,j,k+1)-szz(i,j,k-1))/(2.0_real64*dz)/positive(rho(i,k),1.0_real64)
                    rqx(i,j,k)=-kx(i,k)*px/(positive(eta(i,k),1.0e-6_real64)*positive(rho_f(i,k),1.0_real64))
                    rqy(i,j,k)=-ky(i,k)*py/(positive(eta(i,k),1.0e-6_real64)*positive(rho_f(i,k),1.0_real64))
                    rqz(i,j,k)=-kz(i,k)*pz/(positive(eta(i,k),1.0e-6_real64)*positive(rho_f(i,k),1.0_real64))

                    dampx = speed*((vx(i+1,j,k)-2.0_real64*vx(i,j,k)+vx(i-1,j,k))/(dx*dx))
                    dampy = speed*((vy(i,j+1,k)-2.0_real64*vy(i,j,k)+vy(i,j-1,k))/(dy*dy))
                    dampz = speed*((vz(i,j,k+1)-2.0_real64*vz(i,j,k)+vz(i,j,k-1))/(dz*dz))
                    rvx(i,j,k)=rvx(i,j,k)+dampx; rvy(i,j,k)=rvy(i,j,k)+dampy; rvz(i,j,k)=rvz(i,j,k)+dampz
                enddo
            enddo
        enddo
    end subroutine rhs_biot25_cell

end module biotdg
