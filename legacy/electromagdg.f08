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

        integer :: a, b, i, k
        real(real64) :: det, inv11, inv13, inv33, rhs_x, rhs_z
        real(real64), allocatable :: dex_dx(:,:,:,:), dex_dz(:,:,:,:)
        real(real64), allocatable :: dez_dx(:,:,:,:), dez_dz(:,:,:,:)
        real(real64), allocatable :: dhy_dx(:,:,:,:), dhy_dz(:,:,:,:)

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
                        rhs_x = -dhy_dz(a,b,i,k) - sig11(i,k)*ex(a,b,i,k) - sig13(i,k)*ez(a,b,i,k)
                        rhs_z =  dhy_dx(a,b,i,k) - sig13(i,k)*ex(a,b,i,k) - sig33(i,k)*ez(a,b,i,k)
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

    subroutine rhs_electromagdg25_cell(ex, ey, ez, hx, hy, hz, eps11, eps12, eps13, &
                                       eps22, eps23, eps33, sig11, sig12, sig13, &
                                       sig22, sig23, sig33, dx, dy, dz, &
                                       rex, rey, rez, rhx, rhy, rhz)
        real(real64), intent(in) :: ex(:,:,:), ey(:,:,:), ez(:,:,:)
        real(real64), intent(in) :: hx(:,:,:), hy(:,:,:), hz(:,:,:)
        real(real64), intent(in) :: eps11(:,:), eps12(:,:), eps13(:,:), eps22(:,:), eps23(:,:), eps33(:,:)
        real(real64), intent(in) :: sig11(:,:), sig12(:,:), sig13(:,:), sig22(:,:), sig23(:,:), sig33(:,:)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), intent(out) :: rex(:,:,:), rey(:,:,:), rez(:,:,:)
        real(real64), intent(out) :: rhx(:,:,:), rhy(:,:,:), rhz(:,:,:)

        integer :: i, j, k, nx, ny, nz
        real(real64) :: dHz_dy, dHy_dz, dHx_dz, dHz_dx, dHy_dx, dHx_dy
        real(real64) :: dEy_dz, dEz_dy, dEz_dx, dEx_dz, dEx_dy, dEy_dx
        real(real64) :: det, inv11, inv12, inv13, inv22, inv23, inv33
        real(real64) :: rhs_x, rhs_y, rhs_z

        nx = size(ex, 1); ny = size(ex, 2); nz = size(ex, 3)
        rex = 0.0_real64; rey = 0.0_real64; rez = 0.0_real64
        rhx = 0.0_real64; rhy = 0.0_real64; rhz = 0.0_real64

        !$omp parallel do collapse(3) private(dEy_dz,dEz_dy,dEz_dx,dEx_dz,dEx_dy,dEy_dx) schedule(static)
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1
                    dEy_dz = (ey(i,j,k+1) - ey(i,j,k-1)) / (2.0_real64*dz)
                    dEz_dy = (ez(i,j+1,k) - ez(i,j-1,k)) / (2.0_real64*dy)
                    dEz_dx = (ez(i+1,j,k) - ez(i-1,j,k)) / (2.0_real64*dx)
                    dEx_dz = (ex(i,j,k+1) - ex(i,j,k-1)) / (2.0_real64*dz)
                    dEx_dy = (ex(i,j+1,k) - ex(i,j-1,k)) / (2.0_real64*dy)
                    dEy_dx = (ey(i+1,j,k) - ey(i-1,j,k)) / (2.0_real64*dx)
                    rhx(i,j,k) = (dEy_dz - dEz_dy) / mu0
                    rhy(i,j,k) = (dEz_dx - dEx_dz) / mu0
                    rhz(i,j,k) = (dEx_dy - dEy_dx) / mu0
                enddo
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do collapse(3) private(dHz_dy,dHy_dz,dHx_dz,dHz_dx,dHy_dx,dHx_dy,det, &
        !$omp& inv11,inv12,inv13,inv22,inv23,inv33,rhs_x,rhs_y,rhs_z) schedule(static)
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1
                    dHz_dy = (hz(i,j+1,k) - hz(i,j-1,k)) / (2.0_real64*dy)
                    dHy_dz = (hy(i,j,k+1) - hy(i,j,k-1)) / (2.0_real64*dz)
                    dHx_dz = (hx(i,j,k+1) - hx(i,j,k-1)) / (2.0_real64*dz)
                    dHz_dx = (hz(i+1,j,k) - hz(i-1,j,k)) / (2.0_real64*dx)
                    dHy_dx = (hy(i+1,j,k) - hy(i-1,j,k)) / (2.0_real64*dx)
                    dHx_dy = (hx(i,j+1,k) - hx(i,j-1,k)) / (2.0_real64*dy)

                    det = eps11(i,k)*(eps22(i,k)*eps33(i,k) - eps23(i,k)*eps23(i,k)) - &
                        eps12(i,k)*(eps12(i,k)*eps33(i,k) - eps23(i,k)*eps13(i,k)) + &
                        eps13(i,k)*(eps12(i,k)*eps23(i,k) - eps22(i,k)*eps13(i,k))
                    det = positive(det, 1.0e-24_real64)
                    inv11 = (eps22(i,k)*eps33(i,k) - eps23(i,k)*eps23(i,k)) / det
                    inv12 = -(eps12(i,k)*eps33(i,k) - eps23(i,k)*eps13(i,k)) / det
                    inv13 = (eps12(i,k)*eps23(i,k) - eps22(i,k)*eps13(i,k)) / det
                    inv22 = (eps11(i,k)*eps33(i,k) - eps13(i,k)*eps13(i,k)) / det
                    inv23 = -(eps11(i,k)*eps23(i,k) - eps12(i,k)*eps13(i,k)) / det
                    inv33 = (eps11(i,k)*eps22(i,k) - eps12(i,k)*eps12(i,k)) / det

                    rhs_x = dHz_dy - dHy_dz - sig11(i,k)*ex(i,j,k) - sig12(i,k)*ey(i,j,k) - sig13(i,k)*ez(i,j,k)
                    rhs_y = dHx_dz - dHz_dx - sig12(i,k)*ex(i,j,k) - sig22(i,k)*ey(i,j,k) - sig23(i,k)*ez(i,j,k)
                    rhs_z = dHy_dx - dHx_dy - sig13(i,k)*ex(i,j,k) - sig23(i,k)*ey(i,j,k) - sig33(i,k)*ez(i,j,k)
                    rex(i,j,k) = inv11*rhs_x + inv12*rhs_y + inv13*rhs_z
                    rey(i,j,k) = inv12*rhs_x + inv22*rhs_y + inv23*rhs_z
                    rez(i,j,k) = inv13*rhs_x + inv23*rhs_y + inv33*rhs_z
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine rhs_electromagdg25_cell

    subroutine rhs_electromagdg3_cell(ex, ey, ez, hx, hy, hz, eps11, eps12, eps13, &
                                      eps22, eps23, eps33, sig11, sig12, sig13, &
                                      sig22, sig23, sig33, dx, dy, dz, &
                                      rex, rey, rez, rhx, rhy, rhz)
        real(real64), intent(in) :: ex(:,:,:), ey(:,:,:), ez(:,:,:)
        real(real64), intent(in) :: hx(:,:,:), hy(:,:,:), hz(:,:,:)
        real(real64), intent(in) :: eps11(:,:,:), eps12(:,:,:), eps13(:,:,:)
        real(real64), intent(in) :: eps22(:,:,:), eps23(:,:,:), eps33(:,:,:)
        real(real64), intent(in) :: sig11(:,:,:), sig12(:,:,:), sig13(:,:,:)
        real(real64), intent(in) :: sig22(:,:,:), sig23(:,:,:), sig33(:,:,:)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), intent(out) :: rex(:,:,:), rey(:,:,:), rez(:,:,:)
        real(real64), intent(out) :: rhx(:,:,:), rhy(:,:,:), rhz(:,:,:)

        integer :: i, j, k, nx, ny, nz
        real(real64) :: dHz_dy, dHy_dz, dHx_dz, dHz_dx, dHy_dx, dHx_dy
        real(real64) :: dEy_dz, dEz_dy, dEz_dx, dEx_dz, dEx_dy, dEy_dx
        real(real64) :: det, inv11, inv12, inv13, inv22, inv23, inv33
        real(real64) :: rhs_x, rhs_y, rhs_z

        nx = size(ex, 1); ny = size(ex, 2); nz = size(ex, 3)
        rex = 0.0_real64; rey = 0.0_real64; rez = 0.0_real64
        rhx = 0.0_real64; rhy = 0.0_real64; rhz = 0.0_real64

        !$omp parallel do collapse(3) private(dEy_dz,dEz_dy,dEz_dx,dEx_dz,dEx_dy,dEy_dx) schedule(static)
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1
                    dEy_dz = (ey(i,j,k+1) - ey(i,j,k-1)) / (2.0_real64*dz)
                    dEz_dy = (ez(i,j+1,k) - ez(i,j-1,k)) / (2.0_real64*dy)
                    dEz_dx = (ez(i+1,j,k) - ez(i-1,j,k)) / (2.0_real64*dx)
                    dEx_dz = (ex(i,j,k+1) - ex(i,j,k-1)) / (2.0_real64*dz)
                    dEx_dy = (ex(i,j+1,k) - ex(i,j-1,k)) / (2.0_real64*dy)
                    dEy_dx = (ey(i+1,j,k) - ey(i-1,j,k)) / (2.0_real64*dx)
                    rhx(i,j,k) = (dEy_dz - dEz_dy) / mu0
                    rhy(i,j,k) = (dEz_dx - dEx_dz) / mu0
                    rhz(i,j,k) = (dEx_dy - dEy_dx) / mu0
                enddo
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do collapse(3) private(dHz_dy,dHy_dz,dHx_dz,dHz_dx,dHy_dx,dHx_dy,det, &
        !$omp& inv11,inv12,inv13,inv22,inv23,inv33,rhs_x,rhs_y,rhs_z) schedule(static)
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1
                    dHz_dy = (hz(i,j+1,k) - hz(i,j-1,k)) / (2.0_real64*dy)
                    dHy_dz = (hy(i,j,k+1) - hy(i,j,k-1)) / (2.0_real64*dz)
                    dHx_dz = (hx(i,j,k+1) - hx(i,j,k-1)) / (2.0_real64*dz)
                    dHz_dx = (hz(i+1,j,k) - hz(i-1,j,k)) / (2.0_real64*dx)
                    dHy_dx = (hy(i+1,j,k) - hy(i-1,j,k)) / (2.0_real64*dx)
                    dHx_dy = (hx(i,j+1,k) - hx(i,j-1,k)) / (2.0_real64*dy)

                    det = eps11(i,j,k)*(eps22(i,j,k)*eps33(i,j,k) - eps23(i,j,k)*eps23(i,j,k)) - &
                        eps12(i,j,k)*(eps12(i,j,k)*eps33(i,j,k) - eps23(i,j,k)*eps13(i,j,k)) + &
                        eps13(i,j,k)*(eps12(i,j,k)*eps23(i,j,k) - eps22(i,j,k)*eps13(i,j,k))
                    det = positive(det, 1.0e-24_real64)
                    inv11 = (eps22(i,j,k)*eps33(i,j,k) - eps23(i,j,k)*eps23(i,j,k)) / det
                    inv12 = -(eps12(i,j,k)*eps33(i,j,k) - eps23(i,j,k)*eps13(i,j,k)) / det
                    inv13 = (eps12(i,j,k)*eps23(i,j,k) - eps22(i,j,k)*eps13(i,j,k)) / det
                    inv22 = (eps11(i,j,k)*eps33(i,j,k) - eps13(i,j,k)*eps13(i,j,k)) / det
                    inv23 = -(eps11(i,j,k)*eps23(i,j,k) - eps12(i,j,k)*eps13(i,j,k)) / det
                    inv33 = (eps11(i,j,k)*eps22(i,j,k) - eps12(i,j,k)*eps12(i,j,k)) / det

                    rhs_x = dHz_dy - dHy_dz - sig11(i,j,k)*ex(i,j,k) - sig12(i,j,k)*ey(i,j,k) - sig13(i,j,k)*ez(i,j,k)
                    rhs_y = dHx_dz - dHz_dx - sig12(i,j,k)*ex(i,j,k) - sig22(i,j,k)*ey(i,j,k) - sig23(i,j,k)*ez(i,j,k)
                    rhs_z = dHy_dx - dHx_dy - sig13(i,j,k)*ex(i,j,k) - sig23(i,j,k)*ey(i,j,k) - sig33(i,j,k)*ez(i,j,k)
                    rex(i,j,k) = inv11*rhs_x + inv12*rhs_y + inv13*rhs_z
                    rey(i,j,k) = inv12*rhs_x + inv22*rhs_y + inv23*rhs_z
                    rez(i,j,k) = inv13*rhs_x + inv23*rhs_y + inv33*rhs_z
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine rhs_electromagdg3_cell

    subroutine zero_em_boundaries(ex, ey, ez, hx, hy, hz)
        real(real64), intent(inout) :: ex(:,:,:), ey(:,:,:), ez(:,:,:)
        real(real64), intent(inout) :: hx(:,:,:), hy(:,:,:), hz(:,:,:)
        integer :: nx, ny, nz

        nx = size(ex, 1); ny = size(ex, 2); nz = size(ex, 3)
        ex(1,:,:) = 0.0_real64; ex(nx,:,:) = 0.0_real64
        ex(:,1,:) = 0.0_real64; ex(:,ny,:) = 0.0_real64
        ex(:,:,1) = 0.0_real64; ex(:,:,nz) = 0.0_real64
        ey(1,:,:) = 0.0_real64; ey(nx,:,:) = 0.0_real64
        ey(:,1,:) = 0.0_real64; ey(:,ny,:) = 0.0_real64
        ey(:,:,1) = 0.0_real64; ey(:,:,nz) = 0.0_real64
        ez(1,:,:) = 0.0_real64; ez(nx,:,:) = 0.0_real64
        ez(:,1,:) = 0.0_real64; ez(:,ny,:) = 0.0_real64
        ez(:,:,1) = 0.0_real64; ez(:,:,nz) = 0.0_real64
        hx(1,:,:) = 0.0_real64; hx(nx,:,:) = 0.0_real64
        hx(:,1,:) = 0.0_real64; hx(:,ny,:) = 0.0_real64
        hx(:,:,1) = 0.0_real64; hx(:,:,nz) = 0.0_real64
        hy(1,:,:) = 0.0_real64; hy(nx,:,:) = 0.0_real64
        hy(:,1,:) = 0.0_real64; hy(:,ny,:) = 0.0_real64
        hy(:,:,1) = 0.0_real64; hy(:,:,nz) = 0.0_real64
        hz(1,:,:) = 0.0_real64; hz(nx,:,:) = 0.0_real64
        hz(:,1,:) = 0.0_real64; hz(:,ny,:) = 0.0_real64
        hz(:,:,1) = 0.0_real64; hz(:,:,nz) = 0.0_real64
    end subroutine zero_em_boundaries

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
        call fill_modal_from_cell(cell_ex, ex)
        call fill_modal_from_cell(cell_ez, ez)
        hy = 0.0_real64; srcx = 0.0_real64; srcz = 0.0_real64; fieldnorm = 0.0_real64
        isource = source%xind + domain%cpml; ksource = source%zind + domain%cpml
        eps_floor = positive(min(minval(eps11), minval(eps33)), eps0)
        max_speed = sqrt(1.0_real64 / (mu0 * eps_floor))

        if (source%source_type == 'pw') stop 'Electromagnetic DG 2D plane-wave source is not implemented'
        call loadsource('electromagneticsourcex.dat', source%time_steps, srcx)
        call loadsource('electromagneticsourcez.dat', source%time_steps, srcz)

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
            
            


            call cell_average(ex, weights, cell_out); call write_image2(cell_out, nx, nz, source, it, 'Ex', SINGLE)
            fieldnorm(it) = maxval(abs(cell_out))
            call cell_average(ez, weights, cell_out); call write_image2(cell_out, nx, nz, source, it, 'Ez', SINGLE)
            fieldnorm(it) = max(fieldnorm(it), maxval(abs(cell_out)))
            if (fieldnorm(it) > stability_threshold) stop 'Electromagnetic DG 2D solution became unstable'
        enddo

        call write_array('em_field_norm.dat', source%time_steps, fieldnorm)
        deallocate(eps11, eps13, eps33, sig11, sig13, sig33)
        deallocate(cell_ex, cell_ez, cell_out, ex, ez, hy)
        deallocate(k1ex, k1ez, k1hy, k2ex, k2ez, k2hy, k3ex, k3ez, k3hy, k4ex, k4ez, k4hy)
        deallocate(tex, tez, thy, srcx, srcz, fieldnorm)
    end subroutine electromagdg2

    subroutine electromagdg25(domain, source, SINGLE_OUTPUT)
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        integer :: nx, ny, nz, it, isource, jsource, ksource
        real(real64) :: dx, dy, dz, dt, fieldnorm, eps_floor
        logical :: SINGLE
        real(real64), allocatable :: eps11(:,:), eps12(:,:), eps13(:,:), eps22(:,:), eps23(:,:), eps33(:,:)
        real(real64), allocatable :: sig11(:,:), sig12(:,:), sig13(:,:), sig22(:,:), sig23(:,:), sig33(:,:)
        real(real64), allocatable :: ex(:,:,:), ey(:,:,:), ez(:,:,:), hx(:,:,:), hy(:,:,:), hz(:,:,:)
        real(real64), allocatable :: k1ex(:,:,:), k1ey(:,:,:), k1ez(:,:,:), k1hx(:,:,:), k1hy(:,:,:), k1hz(:,:,:)
        real(real64), allocatable :: k2ex(:,:,:), k2ey(:,:,:), k2ez(:,:,:), k2hx(:,:,:), k2hy(:,:,:), k2hz(:,:,:)
        real(real64), allocatable :: k3ex(:,:,:), k3ey(:,:,:), k3ez(:,:,:), k3hx(:,:,:), k3hy(:,:,:), k3hz(:,:,:)
        real(real64), allocatable :: k4ex(:,:,:), k4ey(:,:,:), k4ez(:,:,:), k4hx(:,:,:), k4hy(:,:,:), k4hz(:,:,:)
        real(real64), allocatable :: tex(:,:,:), tey(:,:,:), tez(:,:,:), thx(:,:,:), thy(:,:,:), thz(:,:,:)
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:), norm(:)

        SINGLE = .TRUE.; if (present(SINGLE_OUTPUT)) SINGLE = SINGLE_OUTPUT
        nx=domain%nx; ny=domain%ny; nz=domain%nz; dx=domain%dx; dy=domain%dy; dz=domain%dz; dt=source%dt
        allocate(eps11(nx,nz), eps12(nx,nz), eps13(nx,nz), eps22(nx,nz), eps23(nx,nz), eps33(nx,nz))
        allocate(sig11(nx,nz), sig12(nx,nz), sig13(nx,nz), sig22(nx,nz), sig23(nx,nz), sig33(nx,nz))
        allocate(ex(nx,ny,nz), ey(nx,ny,nz), ez(nx,ny,nz), hx(nx,ny,nz), hy(nx,ny,nz), hz(nx,ny,nz))
        allocate(k1ex(nx,ny,nz), k1ey(nx,ny,nz), k1ez(nx,ny,nz), k1hx(nx,ny,nz), k1hy(nx,ny,nz), k1hz(nx,ny,nz))
        allocate(k2ex, source=k1ex); allocate(k2ey, source=k1ey); allocate(k2ez, source=k1ez)
        allocate(k2hx, source=k1hx); allocate(k2hy, source=k1hy); allocate(k2hz, source=k1hz)
        allocate(k3ex, source=k1ex); allocate(k3ey, source=k1ey); allocate(k3ez, source=k1ez)
        allocate(k3hx, source=k1hx); allocate(k3hy, source=k1hy); allocate(k3hz, source=k1hz)
        allocate(k4ex, source=k1ex); allocate(k4ey, source=k1ey); allocate(k4ez, source=k1ez)
        allocate(k4hx, source=k1hx); allocate(k4hy, source=k1hy); allocate(k4hz, source=k1hz)
        allocate(tex, source=k1ex); allocate(tey, source=k1ey); allocate(tez, source=k1ez)
        allocate(thx, source=k1hx); allocate(thy, source=k1hy); allocate(thz, source=k1hz)
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps), norm(source%time_steps))

        call material_rw2('e11.dat', eps11, .TRUE.); call material_rw2('e12.dat', eps12, .TRUE.)
        call material_rw2('e13.dat', eps13, .TRUE.); call material_rw2('e22.dat', eps22, .TRUE.)
        call material_rw2('e23.dat', eps23, .TRUE.); call material_rw2('e33.dat', eps33, .TRUE.)
        call material_rw2('s11.dat', sig11, .TRUE.); call material_rw2('s12.dat', sig12, .TRUE.)
        call material_rw2('s13.dat', sig13, .TRUE.); call material_rw2('s22.dat', sig22, .TRUE.)
        call material_rw2('s23.dat', sig23, .TRUE.); call material_rw2('s33.dat', sig33, .TRUE.)
        call material_rw3('initialconditionEx.dat', ex, .TRUE.)
        call material_rw3('initialconditionEy.dat', ey, .TRUE.)
        call material_rw3('initialconditionEz.dat', ez, .TRUE.)

        eps11=eps11*eps0; eps12=eps12*eps0; eps13=eps13*eps0
        eps22=eps22*eps0; eps23=eps23*eps0; eps33=eps33*eps0
        hx=0.0_real64; hy=0.0_real64; hz=0.0_real64
        srcx=0.0_real64; srcy=0.0_real64; srcz=0.0_real64; norm=0.0_real64
        isource=source%xind+domain%cpml; jsource=source%yind+domain%cpml; ksource=source%zind+domain%cpml
        eps_floor = positive(min(minval(eps11), min(minval(eps22), minval(eps33))), eps0)
        if (source%source_type == 'pw') stop 'Electromagnetic DG 2.5D plane-wave source is not implemented'
        call loadsource('electromagneticsourcex.dat', source%time_steps, srcx)
        call loadsource('electromagneticsourcey.dat', source%time_steps, srcy)
        call loadsource('electromagneticsourcez.dat', source%time_steps, srcz)

        do it=1,source%time_steps
            call rhs_electromagdg25_cell(ex,ey,ez,hx,hy,hz,eps11,eps12,eps13,eps22,eps23,eps33, &
                sig11,sig12,sig13,sig22,sig23,sig33,dx,dy,dz,k1ex,k1ey,k1ez,k1hx,k1hy,k1hz)
            tex=ex+0.5_real64*dt*k1ex; tey=ey+0.5_real64*dt*k1ey; tez=ez+0.5_real64*dt*k1ez
            thx=hx+0.5_real64*dt*k1hx; thy=hy+0.5_real64*dt*k1hy; thz=hz+0.5_real64*dt*k1hz
            call rhs_electromagdg25_cell(tex,tey,tez,thx,thy,thz,eps11,eps12,eps13,eps22,eps23,eps33, &
                sig11,sig12,sig13,sig22,sig23,sig33,dx,dy,dz,k2ex,k2ey,k2ez,k2hx,k2hy,k2hz)
            tex=ex+0.5_real64*dt*k2ex; tey=ey+0.5_real64*dt*k2ey; tez=ez+0.5_real64*dt*k2ez
            thx=hx+0.5_real64*dt*k2hx; thy=hy+0.5_real64*dt*k2hy; thz=hz+0.5_real64*dt*k2hz
            call rhs_electromagdg25_cell(tex,tey,tez,thx,thy,thz,eps11,eps12,eps13,eps22,eps23,eps33, &
                sig11,sig12,sig13,sig22,sig23,sig33,dx,dy,dz,k3ex,k3ey,k3ez,k3hx,k3hy,k3hz)
            tex=ex+dt*k3ex; tey=ey+dt*k3ey; tez=ez+dt*k3ez
            thx=hx+dt*k3hx; thy=hy+dt*k3hy; thz=hz+dt*k3hz
            call rhs_electromagdg25_cell(tex,tey,tez,thx,thy,thz,eps11,eps12,eps13,eps22,eps23,eps33, &
                sig11,sig12,sig13,sig22,sig23,sig33,dx,dy,dz,k4ex,k4ey,k4ez,k4hx,k4hy,k4hz)

            ex=ex+dt*(k1ex+2.0_real64*k2ex+2.0_real64*k3ex+k4ex)/6.0_real64
            ey=ey+dt*(k1ey+2.0_real64*k2ey+2.0_real64*k3ey+k4ey)/6.0_real64
            ez=ez+dt*(k1ez+2.0_real64*k2ez+2.0_real64*k3ez+k4ez)/6.0_real64
            hx=hx+dt*(k1hx+2.0_real64*k2hx+2.0_real64*k3hx+k4hx)/6.0_real64
            hy=hy+dt*(k1hy+2.0_real64*k2hy+2.0_real64*k3hy+k4hy)/6.0_real64
            hz=hz+dt*(k1hz+2.0_real64*k2hz+2.0_real64*k3hz+k4hz)/6.0_real64
            
            
            if ( source%source_type == 'pw' ) then 
                t = it * source%dt
                if ( active(1) .or. active(2) ) then 
                    ii = nint(XLOC / domain%dx)
                    do jj = j_min,j_max
                        do kk = k_min, k_max
                            r = (/ XLOC, jj * domain%dy, kk * domain%dz /)
                            Ev = (/ Ex(ii,jj,kk), Ey(ii,jj,kk), Ez(ii,jj,kk) /)
                            Hv = (/ Hx(ii,jj,kk), Hy(ii,jj,kk), Hz(ii,jj,kk) /)
                            call boundary_em_field(r, r0, t, &
                                            source%source_frequency, &
                                            source%x_z_rotation, &
                                            source%x_y_rotation, &
                                            source%y_z_rotation, &
                                            vbackground, eta, &
                                            source%amplitude, &
                                            source%source_wavelet, &
                                            Ev, Hv )
                            Ex(ii,jj,kk) = Ex(ii,jj,kk) + Ev(1) 
                            Ey(ii,jj,kk) = Ey(ii,jj,kk) + Ev(2) 
                            Ez(ii,jj,kk) = Ez(ii,jj,kk) + Ev(3) 
                            Hx(ii,jj,kk) = Hx(ii,jj,kk) + Hv(1)
                            Hy(ii,jj,kk) = Hy(ii,jj,kk) + Hv(2)
                            Hz(ii,jj,kk) = Hz(ii,jj,kk) + Hv(3)
                        enddo
                    enddo
                endif 
                if ( active(3) .or. active(4)) then 
                    jj = nint(YLOC / domain%dy)
                    do ii = i_min,i_max
                        do kk = k_min, k_max
                            r = (/ ii * domain%dx, YLOC, kk * domain%dz /)
                            Ev = (/ Ex(ii,jj,kk), Ey(ii,jj,kk), Ez(ii,jj,kk) /)
                            Hv = (/ Hx(ii,jj,kk), Hy(ii,jj,kk), Hz(ii,jj,kk) /)
                            call boundary_em_field(r, r0, t, &
                                            source%source_frequency, &
                                            source%x_z_rotation, &
                                            source%x_y_rotation, &
                                            source%y_z_rotation, &
                                            vbackground, eta, &
                                            source%amplitude, &
                                            source%source_wavelet, &
                                            Ev, Hv )
                            Ex(ii,jj,kk) = Ex(ii,jj,kk) + Ev(1) 
                            Ey(ii,jj,kk) = Ey(ii,jj,kk) + Ev(2) 
                            Ez(ii,jj,kk) = Ez(ii,jj,kk) + Ev(3) 
                            Hx(ii,jj,kk) = Hx(ii,jj,kk) + Hv(1)
                            Hy(ii,jj,kk) = Hy(ii,jj,kk) + Hv(2)
                            Hz(ii,jj,kk) = Hz(ii,jj,kk) + Hv(3)
                        enddo
                    enddo
                endif 
                if ( active(5) .or. active(6)) then 
                    kk = nint(ZLOC / domain%dz)
                    do ii = i_min,i_max
                        do jj = j_min, j_max
                            r = (/ ii * domain%dx, jj * domain%dy, ZLOC/)
                            Ev = (/ Ex(ii,jj,kk), Ey(ii,jj,kk), Ez(ii,jj,kk) /)
                            Hv = (/ Hx(ii,jj,kk), Hy(ii,jj,kk), Hz(ii,jj,kk) /)
                            call boundary_em_field(r, r0, t, &
                                            source%source_frequency, &
                                            source%x_z_rotation, &
                                            source%x_y_rotation, &
                                            source%y_z_rotation, &
                                            vbackground, eta, &
                                            source%amplitude, &
                                            source%source_wavelet, &
                                            Ev, Hv )
                            Ex(ii,jj,kk) = Ex(ii,jj,kk) + Ev(1) 
                            Ey(ii,jj,kk) = Ey(ii,jj,kk) + Ev(2) 
                            Ez(ii,jj,kk) = Ez(ii,jj,kk) + Ev(3) 
                            Hx(ii,jj,kk) = Hx(ii,jj,kk) + Hv(1)
                            Hy(ii,jj,kk) = Hy(ii,jj,kk) + Hv(2)
                            Hz(ii,jj,kk) = Hz(ii,jj,kk) + Hv(3)
                        enddo
                    enddo
                endif 
            else
                ! add the source (force vector located at a given grid point)
                Ex(isource,jsource,ksource) = Ex(isource,jsource,ksource) + & 
                            srcx(it) * dt / eps11(isource,ksource)
                Ey(isource,jsource,ksource) = Ey(isource,jsource,ksource) + & 
                            srcy(it) * dt / eps22(isource,ksource) 
                Ez(isource,jsource,ksource) = Ez(isource,jsource,ksource) + & 
                            srcz(it) * dt / eps33(isource,ksource)
            endif 

            call zero_em_boundaries(ex, ey, ez, hx, hy, hz)
            fieldnorm=maxval(sqrt(ex**2+ey**2+ez**2)); norm(it)=fieldnorm
            if (fieldnorm > stability_threshold) stop 'Electromagnetic DG 2.5D solution became unstable'
            call write_image3(ex,nx,ny,nz,source,it,'Ex',SINGLE)
            call write_image3(ey,nx,ny,nz,source,it,'Ey',SINGLE)
            call write_image3(ez,nx,ny,nz,source,it,'Ez',SINGLE)
        enddo

        call write_array('em_field_norm.dat', source%time_steps, norm)
        deallocate(eps11, eps12, eps13, eps22, eps23, eps33, sig11, sig12, sig13, sig22, sig23, sig33)
        deallocate(ex, ey, ez, hx, hy, hz, k1ex, k1ey, k1ez, k1hx, k1hy, k1hz)
        deallocate(k2ex, k2ey, k2ez, k2hx, k2hy, k2hz, k3ex, k3ey, k3ez, k3hx, k3hy, k3hz)
        deallocate(k4ex, k4ey, k4ez, k4hx, k4hy, k4hz, tex, tey, tez, thx, thy, thz)
        deallocate(srcx, srcy, srcz, norm)
    end subroutine electromagdg25

    subroutine electromagdg3(domain, source, SINGLE_OUTPUT)
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        integer :: nx, ny, nz, it, isource, jsource, ksource
        real(real64) :: dx, dy, dz, dt, fieldnorm, eps_floor
        logical :: SINGLE
        real(real64), allocatable :: eps11(:,:,:), eps12(:,:,:), eps13(:,:,:), eps22(:,:,:), eps23(:,:,:), eps33(:,:,:)
        real(real64), allocatable :: sig11(:,:,:), sig12(:,:,:), sig13(:,:,:), sig22(:,:,:), sig23(:,:,:), sig33(:,:,:)
        real(real64), allocatable :: ex(:,:,:), ey(:,:,:), ez(:,:,:), hx(:,:,:), hy(:,:,:), hz(:,:,:)
        real(real64), allocatable :: k1ex(:,:,:), k1ey(:,:,:), k1ez(:,:,:), k1hx(:,:,:), k1hy(:,:,:), k1hz(:,:,:)
        real(real64), allocatable :: k2ex(:,:,:), k2ey(:,:,:), k2ez(:,:,:), k2hx(:,:,:), k2hy(:,:,:), k2hz(:,:,:)
        real(real64), allocatable :: k3ex(:,:,:), k3ey(:,:,:), k3ez(:,:,:), k3hx(:,:,:), k3hy(:,:,:), k3hz(:,:,:)
        real(real64), allocatable :: k4ex(:,:,:), k4ey(:,:,:), k4ez(:,:,:), k4hx(:,:,:), k4hy(:,:,:), k4hz(:,:,:)
        real(real64), allocatable :: tex(:,:,:), tey(:,:,:), tez(:,:,:), thx(:,:,:), thy(:,:,:), thz(:,:,:)
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:), norm(:)

        SINGLE = .TRUE.; if (present(SINGLE_OUTPUT)) SINGLE = SINGLE_OUTPUT
        nx=domain%nx; ny=domain%ny; nz=domain%nz; dx=domain%dx; dy=domain%dy; dz=domain%dz; dt=source%dt
        allocate(eps11(nx,ny,nz), eps12(nx,ny,nz), eps13(nx,ny,nz), eps22(nx,ny,nz), eps23(nx,ny,nz), eps33(nx,ny,nz))
        allocate(sig11(nx,ny,nz), sig12(nx,ny,nz), sig13(nx,ny,nz), sig22(nx,ny,nz), sig23(nx,ny,nz), sig33(nx,ny,nz))
        allocate(ex(nx,ny,nz), ey(nx,ny,nz), ez(nx,ny,nz), hx(nx,ny,nz), hy(nx,ny,nz), hz(nx,ny,nz))
        allocate(k1ex(nx,ny,nz), k1ey(nx,ny,nz), k1ez(nx,ny,nz), k1hx(nx,ny,nz), k1hy(nx,ny,nz), k1hz(nx,ny,nz))
        allocate(k2ex, source=k1ex); allocate(k2ey, source=k1ey); allocate(k2ez, source=k1ez)
        allocate(k2hx, source=k1hx); allocate(k2hy, source=k1hy); allocate(k2hz, source=k1hz)
        allocate(k3ex, source=k1ex); allocate(k3ey, source=k1ey); allocate(k3ez, source=k1ez)
        allocate(k3hx, source=k1hx); allocate(k3hy, source=k1hy); allocate(k3hz, source=k1hz)
        allocate(k4ex, source=k1ex); allocate(k4ey, source=k1ey); allocate(k4ez, source=k1ez)
        allocate(k4hx, source=k1hx); allocate(k4hy, source=k1hy); allocate(k4hz, source=k1hz)
        allocate(tex, source=k1ex); allocate(tey, source=k1ey); allocate(tez, source=k1ez)
        allocate(thx, source=k1hx); allocate(thy, source=k1hy); allocate(thz, source=k1hz)
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps), norm(source%time_steps))

        call material_rw3('e11.dat', eps11, .TRUE.); call material_rw3('e12.dat', eps12, .TRUE.)
        call material_rw3('e13.dat', eps13, .TRUE.); call material_rw3('e22.dat', eps22, .TRUE.)
        call material_rw3('e23.dat', eps23, .TRUE.); call material_rw3('e33.dat', eps33, .TRUE.)
        call material_rw3('s11.dat', sig11, .TRUE.); call material_rw3('s12.dat', sig12, .TRUE.)
        call material_rw3('s13.dat', sig13, .TRUE.); call material_rw3('s22.dat', sig22, .TRUE.)
        call material_rw3('s23.dat', sig23, .TRUE.); call material_rw3('s33.dat', sig33, .TRUE.)
        call material_rw3('initialconditionEx.dat', ex, .TRUE.)
        call material_rw3('initialconditionEy.dat', ey, .TRUE.)
        call material_rw3('initialconditionEz.dat', ez, .TRUE.)

        eps11=eps11*eps0; eps12=eps12*eps0; eps13=eps13*eps0
        eps22=eps22*eps0; eps23=eps23*eps0; eps33=eps33*eps0
        hx=0.0_real64; hy=0.0_real64; hz=0.0_real64
        srcx=0.0_real64; srcy=0.0_real64; srcz=0.0_real64; norm=0.0_real64
        isource=source%xind+domain%cpml; jsource=source%yind+domain%cpml; ksource=source%zind+domain%cpml
        eps_floor = positive(min(minval(eps11), min(minval(eps22), minval(eps33))), eps0)
        if (source%source_type == 'pw') stop 'Electromagnetic DG 3D plane-wave source is not implemented'
        call loadsource('electromagneticsourcex.dat', source%time_steps, srcx)
        call loadsource('electromagneticsourcey.dat', source%time_steps, srcy)
        call loadsource('electromagneticsourcez.dat', source%time_steps, srcz)

        do it=1,source%time_steps
            call rhs_electromagdg3_cell(ex,ey,ez,hx,hy,hz,eps11,eps12,eps13,eps22,eps23,eps33, &
                sig11,sig12,sig13,sig22,sig23,sig33,dx,dy,dz,k1ex,k1ey,k1ez,k1hx,k1hy,k1hz)
            tex=ex+0.5_real64*dt*k1ex; tey=ey+0.5_real64*dt*k1ey; tez=ez+0.5_real64*dt*k1ez
            thx=hx+0.5_real64*dt*k1hx; thy=hy+0.5_real64*dt*k1hy; thz=hz+0.5_real64*dt*k1hz
            call rhs_electromagdg3_cell(tex,tey,tez,thx,thy,thz,eps11,eps12,eps13,eps22,eps23,eps33, &
                sig11,sig12,sig13,sig22,sig23,sig33,dx,dy,dz,k2ex,k2ey,k2ez,k2hx,k2hy,k2hz)
            tex=ex+0.5_real64*dt*k2ex; tey=ey+0.5_real64*dt*k2ey; tez=ez+0.5_real64*dt*k2ez
            thx=hx+0.5_real64*dt*k2hx; thy=hy+0.5_real64*dt*k2hy; thz=hz+0.5_real64*dt*k2hz
            call rhs_electromagdg3_cell(tex,tey,tez,thx,thy,thz,eps11,eps12,eps13,eps22,eps23,eps33, &
                sig11,sig12,sig13,sig22,sig23,sig33,dx,dy,dz,k3ex,k3ey,k3ez,k3hx,k3hy,k3hz)
            tex=ex+dt*k3ex; tey=ey+dt*k3ey; tez=ez+dt*k3ez
            thx=hx+dt*k3hx; thy=hy+dt*k3hy; thz=hz+dt*k3hz
            call rhs_electromagdg3_cell(tex,tey,tez,thx,thy,thz,eps11,eps12,eps13,eps22,eps23,eps33, &
                sig11,sig12,sig13,sig22,sig23,sig33,dx,dy,dz,k4ex,k4ey,k4ez,k4hx,k4hy,k4hz)

            ex=ex+dt*(k1ex+2.0_real64*k2ex+2.0_real64*k3ex+k4ex)/6.0_real64
            ey=ey+dt*(k1ey+2.0_real64*k2ey+2.0_real64*k3ey+k4ey)/6.0_real64
            ez=ez+dt*(k1ez+2.0_real64*k2ez+2.0_real64*k3ez+k4ez)/6.0_real64
            hx=hx+dt*(k1hx+2.0_real64*k2hx+2.0_real64*k3hx+k4hx)/6.0_real64
            hy=hy+dt*(k1hy+2.0_real64*k2hy+2.0_real64*k3hy+k4hy)/6.0_real64
            hz=hz+dt*(k1hz+2.0_real64*k2hz+2.0_real64*k3hz+k4hz)/6.0_real64
            
            
            ! Add the source terms
            if ( source%source_type == 'pw' ) then
                t = it * source%dt
                if ( active(1) .or. active(2) ) then
                    ii = nint(XLOC / domain%dx)
                    do jj = j_min,j_max
                        do kk = k_min,k_max
                            r = (/ XLOC, jj * domain%dy, kk * domain%dz /)
                            call boundary_em_field(r, r0, t, &
                                            source%source_frequency, &
                                            source%x_z_rotation, &
                                            source%x_y_rotation, &
                                            source%y_z_rotation, &
                                            vbackground, eta, &
                                            source%amplitude, &
                                            source%source_wavelet, &
                                            Ev, Hv )
                            Ex(ii,jj,kk) = Ex(ii,jj,kk) + Ev(1)
                            Ey(ii,jj,kk) = Ey(ii,jj,kk) + Ev(2)
                            Ez(ii,jj,kk) = Ez(ii,jj,kk) + Ev(3)
                            Hx(ii,jj,kk) = Hx(ii,jj,kk) + Hv(1)
                            Hy(ii,jj,kk) = Hy(ii,jj,kk) + Hv(2)
                            Hz(ii,jj,kk) = Hz(ii,jj,kk) + Hv(3)
                        enddo
                    enddo
                endif
                if ( active(3) .or. active(4) ) then
                    jj = nint(YLOC / domain%dy)
                    do ii = i_min,i_max
                        do kk = k_min,k_max
                            r = (/ ii * domain%dx, YLOC, kk * domain%dz /)
                            call boundary_em_field(r, r0, t, &
                                            source%source_frequency, &
                                            source%x_z_rotation, &
                                            source%x_y_rotation, &
                                            source%y_z_rotation, &
                                            vbackground, eta, &
                                            source%amplitude, &
                                            source%source_wavelet, &
                                            Ev, Hv )
                            Ex(ii,jj,kk) = Ex(ii,jj,kk) + Ev(1)
                            Ey(ii,jj,kk) = Ey(ii,jj,kk) + Ev(2)
                            Ez(ii,jj,kk) = Ez(ii,jj,kk) + Ev(3)
                            Hx(ii,jj,kk) = Hx(ii,jj,kk) + Hv(1)
                            Hy(ii,jj,kk) = Hy(ii,jj,kk) + Hv(2)
                            Hz(ii,jj,kk) = Hz(ii,jj,kk) + Hv(3)
                        enddo
                    enddo
                endif
                if ( active(5) .or. active(6) ) then
                    kk = nint(ZLOC / domain%dz)
                    do ii = i_min,i_max
                        do jj = j_min,j_max
                            r = (/ ii * domain%dx, jj * domain%dy, ZLOC /)
                            call boundary_em_field(r, r0, t, &
                                            source%source_frequency, &
                                            source%x_z_rotation, &
                                            source%x_y_rotation, &
                                            source%y_z_rotation, &
                                            vbackground, eta, &
                                            source%amplitude, &
                                            source%source_wavelet, &
                                            Ev, Hv )
                            Ex(ii,jj,kk) = Ex(ii,jj,kk) + Ev(1)
                            Ey(ii,jj,kk) = Ey(ii,jj,kk) + Ev(2)
                            Ez(ii,jj,kk) = Ez(ii,jj,kk) + Ev(3)
                            Hx(ii,jj,kk) = Hx(ii,jj,kk) + Hv(1)
                            Hy(ii,jj,kk) = Hy(ii,jj,kk) + Hv(2)
                            Hz(ii,jj,kk) = Hz(ii,jj,kk) + Hv(3)
                        enddo
                    enddo
                endif
            else
                ! add the source (force vector located at a given grid point)
                Ex(isource,jsource,ksource) = Ex(isource,jsource,ksource) + &
                            srcx(it) * dt / eps11(isource,jsource,ksource)
                Ey(isource,jsource,ksource) = Ey(isource,jsource,ksource) + &
                            srcy(it) * dt / eps22(isource,jsource,ksource)
                Ez(isource,jsource,ksource) = Ez(isource,jsource,ksource) + &
                            srcz(it) * dt / eps33(isource,jsource,ksource)
            endif
            
            call zero_em_boundaries(ex, ey, ez, hx, hy, hz)
            fieldnorm=maxval(sqrt(ex**2+ey**2+ez**2)); norm(it)=fieldnorm
            if (fieldnorm > stability_threshold) stop 'Electromagnetic DG 3D solution became unstable'
            call write_image3(ex,nx,ny,nz,source,it,'Ex',SINGLE)
            call write_image3(ey,nx,ny,nz,source,it,'Ey',SINGLE)
            call write_image3(ez,nx,ny,nz,source,it,'Ez',SINGLE)
        enddo

        call write_array('em_field_norm.dat', source%time_steps, norm)
        deallocate(eps11, eps12, eps13, eps22, eps23, eps33, sig11, sig12, sig13, sig22, sig23, sig33)
        deallocate(ex, ey, ez, hx, hy, hz, k1ex, k1ey, k1ez, k1hx, k1hy, k1hz)
        deallocate(k2ex, k2ey, k2ez, k2hx, k2hy, k2hz, k3ex, k3ey, k3ez, k3hx, k3hy, k3hz)
        deallocate(k4ex, k4ey, k4ez, k4hx, k4hy, k4hz, tex, tey, tez, thx, thy, thz)
        deallocate(srcx, srcy, srcz, norm)
    end subroutine electromagdg3

end module electromagdg
