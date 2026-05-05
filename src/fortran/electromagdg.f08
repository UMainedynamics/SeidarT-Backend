module electromagdg

    use iso_fortran_env, only: real64
    use seidartio
    use seidart_types
    use constants
    use discontinuous_galerkin_methods, only: DG_P, positive, initialize_gll4, &
        fill_modal_from_cell, cell_average, dg_derivatives, apply_rusanov_penalty

    implicit none

contains

    subroutine rhs_electromagdg2(ex, ez, hy, eps11, eps13, eps33, sig11, sig13, sig33, &
                                 deriv, weights, dx, dz, max_speed, rex, rez, rhy)
        real(real64), intent(in) :: ex(0:,0:,:,:), ez(0:,0:,:,:), hy(0:,0:,:,:)
        real(real64), intent(in) :: eps11(:,:), eps13(:,:), eps33(:,:)
        real(real64), intent(in) :: sig11(:,:), sig13(:,:), sig33(:,:)
        real(real64), intent(in) :: deriv(0:DG_P,0:DG_P), weights(0:DG_P)
        real(real64), intent(in) :: dx, dz, max_speed
        real(real64), intent(out) :: rex(0:,0:,:,:), rez(0:,0:,:,:), rhy(0:,0:,:,:)

        real(real64), allocatable :: dex_dx(:,:,:,:), dex_dz(:,:,:,:)
        real(real64), allocatable :: dez_dx(:,:,:,:), dez_dz(:,:,:,:)
        real(real64), allocatable :: dhy_dx(:,:,:,:), dhy_dz(:,:,:,:)
        integer :: a, b, i, k
        real(real64) :: det, inv11, inv13, inv33, rhs_x, rhs_z

        allocate(dex_dx(0:DG_P,0:DG_P,size(ex,3),size(ex,4)))
        allocate(dex_dz(0:DG_P,0:DG_P,size(ex,3),size(ex,4)))
        allocate(dez_dx(0:DG_P,0:DG_P,size(ex,3),size(ex,4)))
        allocate(dez_dz(0:DG_P,0:DG_P,size(ex,3),size(ex,4)))
        allocate(dhy_dx(0:DG_P,0:DG_P,size(ex,3),size(ex,4)))
        allocate(dhy_dz(0:DG_P,0:DG_P,size(ex,3),size(ex,4)))

        call dg_derivatives(ex, deriv, dx, dz, dex_dx, dex_dz)
        call dg_derivatives(ez, deriv, dx, dz, dez_dx, dez_dz)
        call dg_derivatives(hy, deriv, dx, dz, dhy_dx, dhy_dz)

        !$omp parallel do collapse(2) private(a,b,det,inv11,inv13,inv33,rhs_x,rhs_z) schedule(static)
        do k = 1, size(ex, 4)
            do i = 1, size(ex, 3)
                det = eps11(i,k)*eps33(i,k) - eps13(i,k)*eps13(i,k)
                det = positive(det, 1.0e-24_real64)
                inv11 = eps33(i,k) / det
                inv13 = -eps13(i,k) / det
                inv33 = eps11(i,k) / det
                do b = 0, DG_P
                    do a = 0, DG_P
                        rhs_x = -dhy_dz(a,b,i,k) - sig11(i,k)*ex(a,b,i,k) - &
                            sig13(i,k)*ez(a,b,i,k)
                        rhs_z = dhy_dx(a,b,i,k) - sig13(i,k)*ex(a,b,i,k) - &
                            sig33(i,k)*ez(a,b,i,k)
                        rex(a,b,i,k) = inv11*rhs_x + inv13*rhs_z
                        rez(a,b,i,k) = inv13*rhs_x + inv33*rhs_z
                        rhy(a,b,i,k) = (dez_dx(a,b,i,k) - dex_dz(a,b,i,k)) / mu0
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do

        call apply_rusanov_penalty(ex, rex, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(ez, rez, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(hy, rhy, weights, dx, dz, max_speed)

        deallocate(dex_dx, dex_dz, dez_dx, dez_dz, dhy_dx, dhy_dz)
    end subroutine rhs_electromagdg2

    subroutine electromagdg2(domain, source, SINGLE_OUTPUT)
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        integer :: nx, nz, it, isource, ksource
        real(real64) :: dx, dz, dt, max_speed, eps_floor
        logical :: SINGLE
        real(real64) :: nodes(0:DG_P), weights(0:DG_P), deriv(0:DG_P,0:DG_P)
        real(real64), allocatable :: eps11(:,:), eps13(:,:), eps33(:,:), sig11(:,:), sig13(:,:), sig33(:,:)
        real(real64), allocatable :: cell_ex(:,:), cell_ez(:,:), cell_out(:,:)
        real(real64), allocatable :: ex(:,:,:,:), ez(:,:,:,:), hy(:,:,:,:)
        real(real64), allocatable :: k1ex(:,:,:,:), k1ez(:,:,:,:), k1hy(:,:,:,:)
        real(real64), allocatable :: k2ex(:,:,:,:), k2ez(:,:,:,:), k2hy(:,:,:,:)
        real(real64), allocatable :: k3ex(:,:,:,:), k3ez(:,:,:,:), k3hy(:,:,:,:)
        real(real64), allocatable :: k4ex(:,:,:,:), k4ez(:,:,:,:), k4hy(:,:,:,:)
        real(real64), allocatable :: tex(:,:,:,:), tez(:,:,:,:), thy(:,:,:,:)
        real(real64), allocatable :: srcx(:), srcz(:), fieldnorm(:)

        SINGLE = .TRUE.; if (present(SINGLE_OUTPUT)) SINGLE = SINGLE_OUTPUT
        nx = domain%nx; nz = domain%nz; dx = domain%dx; dz = domain%dz; dt = source%dt
        call initialize_gll4(nodes, weights, deriv)

        allocate(eps11(nx,nz), eps13(nx,nz), eps33(nx,nz), sig11(nx,nz), sig13(nx,nz), sig33(nx,nz))
        allocate(cell_ex(nx,nz), cell_ez(nx,nz), cell_out(nx,nz))
        allocate(ex(0:DG_P,0:DG_P,nx,nz), ez(0:DG_P,0:DG_P,nx,nz), hy(0:DG_P,0:DG_P,nx,nz))
        allocate(k1ex(0:DG_P,0:DG_P,nx,nz), k1ez(0:DG_P,0:DG_P,nx,nz), k1hy(0:DG_P,0:DG_P,nx,nz))
        allocate(k2ex, source=k1ex); allocate(k2ez, source=k1ez); allocate(k2hy, source=k1hy)
        allocate(k3ex, source=k1ex); allocate(k3ez, source=k1ez); allocate(k3hy, source=k1hy)
        allocate(k4ex, source=k1ex); allocate(k4ez, source=k1ez); allocate(k4hy, source=k1hy)
        allocate(tex, source=k1ex); allocate(tez, source=k1ez); allocate(thy, source=k1hy)
        allocate(srcx(source%time_steps), srcz(source%time_steps), fieldnorm(source%time_steps))

        call material_rw2('e11.dat', eps11, .TRUE.)
        call material_rw2('e13.dat', eps13, .TRUE.)
        call material_rw2('e33.dat', eps33, .TRUE.)
        call material_rw2('s11.dat', sig11, .TRUE.)
        call material_rw2('s13.dat', sig13, .TRUE.)
        call material_rw2('s33.dat', sig33, .TRUE.)
        call material_rw2('initialconditionEx.dat', cell_ex, .TRUE.)
        call material_rw2('initialconditionEz.dat', cell_ez, .TRUE.)

        eps11 = eps11 * eps0; eps13 = eps13 * eps0; eps33 = eps33 * eps0
        call fill_modal_from_cell(cell_ex, ex); call fill_modal_from_cell(cell_ez, ez)
        hy = 0.0_real64; srcx = 0.0_real64; srcz = 0.0_real64; fieldnorm = 0.0_real64
        isource = source%xind + domain%cpml; ksource = source%zind + domain%cpml
        eps_floor = positive(min(minval(eps11), minval(eps33)), eps0)
        max_speed = sqrt(1.0_real64 / (mu0 * eps_floor))

        if (source%source_type /= 'pw') then
            call loadsource('electromagneticsourcex.dat', source%time_steps, srcx)
            call loadsource('electromagneticsourcez.dat', source%time_steps, srcz)
        endif

        do it = 1, source%time_steps
            call rhs_electromagdg2(ex,ez,hy,eps11,eps13,eps33,sig11,sig13,sig33, &
                deriv,weights,dx,dz,max_speed,k1ex,k1ez,k1hy)
            tex=ex+0.5_real64*dt*k1ex; tez=ez+0.5_real64*dt*k1ez; thy=hy+0.5_real64*dt*k1hy
            call rhs_electromagdg2(tex,tez,thy,eps11,eps13,eps33,sig11,sig13,sig33, &
                deriv,weights,dx,dz,max_speed,k2ex,k2ez,k2hy)
            tex=ex+0.5_real64*dt*k2ex; tez=ez+0.5_real64*dt*k2ez; thy=hy+0.5_real64*dt*k2hy
            call rhs_electromagdg2(tex,tez,thy,eps11,eps13,eps33,sig11,sig13,sig33, &
                deriv,weights,dx,dz,max_speed,k3ex,k3ez,k3hy)
            tex=ex+dt*k3ex; tez=ez+dt*k3ez; thy=hy+dt*k3hy
            call rhs_electromagdg2(tex,tez,thy,eps11,eps13,eps33,sig11,sig13,sig33, &
                deriv,weights,dx,dz,max_speed,k4ex,k4ez,k4hy)

            ex = ex + dt*(k1ex+2.0_real64*k2ex+2.0_real64*k3ex+k4ex)/6.0_real64
            ez = ez + dt*(k1ez+2.0_real64*k2ez+2.0_real64*k3ez+k4ez)/6.0_real64
            hy = hy + dt*(k1hy+2.0_real64*k2hy+2.0_real64*k3hy+k4hy)/6.0_real64

            ex(0,0,isource,ksource) = ex(0,0,isource,ksource) + srcx(it)*dt/positive(eps11(isource,ksource), eps0)
            ez(0,0,isource,ksource) = ez(0,0,isource,ksource) + srcz(it)*dt/positive(eps33(isource,ksource), eps0)

            call cell_average(ex, weights, cell_out); call write_image2(cell_out, nx, nz, source, it, 'Ex', SINGLE)
            fieldnorm(it) = maxval(abs(cell_out))
            call cell_average(ez, weights, cell_out); call write_image2(cell_out, nx, nz, source, it, 'Ez', SINGLE)
            fieldnorm(it) = max(fieldnorm(it), maxval(abs(cell_out)))
            if (fieldnorm(it) > stability_threshold) stop 'Electromagnetic DG 2D solution became unstable'
        enddo

        call write_array('em_field_norm.dat', source%time_steps, fieldnorm)
    end subroutine electromagdg2

end module electromagdg
