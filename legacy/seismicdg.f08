module seismicdg

    use iso_fortran_env, only: real64
    use seidartio
    use seidart_types
    use constants
    use discontinuous_galerkin_methods, only: DG_P, positive, initialize_gll4, &
        fill_modal_from_cell, cell_average, dg_derivatives, apply_rusanov_penalty

    implicit none

contains

    subroutine rhs_seismicdg2(vx, vz, sigmaxx, sigmazz, sigmaxz, &
                              c11, c13, c15, c33, c35, c55, rho, &
                              gamma_x, gamma_z, gamma_xz, &
                              deriv, weights, dx, dz, max_speed, &
                              rvx, rvz, rsxx, rszz, rsxz)
        real(real64), intent(in) :: vx(0:,0:,:,:), vz(0:,0:,:,:)
        real(real64), intent(in) :: sigmaxx(0:,0:,:,:), sigmazz(0:,0:,:,:)
        real(real64), intent(in) :: sigmaxz(0:,0:,:,:)
        real(real64), intent(in) :: c11(:,:), c13(:,:), c15(:,:)
        real(real64), intent(in) :: c33(:,:), c35(:,:), c55(:,:), rho(:,:)
        real(real64), intent(in) :: gamma_x(:,:), gamma_z(:,:), gamma_xz(:,:)
        real(real64), intent(in) :: deriv(0:DG_P,0:DG_P), weights(0:DG_P)
        real(real64), intent(in) :: dx, dz, max_speed
        real(real64), intent(out) :: rvx(0:,0:,:,:), rvz(0:,0:,:,:)
        real(real64), intent(out) :: rsxx(0:,0:,:,:), rszz(0:,0:,:,:), rsxz(0:,0:,:,:)

        real(real64), allocatable :: dvx_dx(:,:,:,:), dvx_dz(:,:,:,:)
        real(real64), allocatable :: dvz_dx(:,:,:,:), dvz_dz(:,:,:,:)
        real(real64), allocatable :: dsxx_dx(:,:,:,:), dsxx_dz(:,:,:,:)
        real(real64), allocatable :: dszz_dx(:,:,:,:), dszz_dz(:,:,:,:)
        real(real64), allocatable :: dsxz_dx(:,:,:,:), dsxz_dz(:,:,:,:)
        integer :: a, b, i, k

        allocate(dvx_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dvx_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dvz_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dvz_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dsxx_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dsxx_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dszz_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dszz_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dsxz_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dsxz_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))

        call dg_derivatives(vx, deriv, dx, dz, dvx_dx, dvx_dz)
        call dg_derivatives(vz, deriv, dx, dz, dvz_dx, dvz_dz)
        call dg_derivatives(sigmaxx, deriv, dx, dz, dsxx_dx, dsxx_dz)
        call dg_derivatives(sigmazz, deriv, dx, dz, dszz_dx, dszz_dz)
        call dg_derivatives(sigmaxz, deriv, dx, dz, dsxz_dx, dsxz_dz)

        !$omp parallel do collapse(2) private(a,b) schedule(static)
        do k = 1, size(vx, 4)
            do i = 1, size(vx, 3)
                do b = 0, DG_P
                    do a = 0, DG_P
                        rsxx(a,b,i,k) = c11(i,k)*dvx_dx(a,b,i,k) + &
                            c13(i,k)*dvz_dz(a,b,i,k) + &
                            c15(i,k)*(dvz_dx(a,b,i,k) + dvx_dz(a,b,i,k)) - &
                            gamma_x(i,k)*sigmaxx(a,b,i,k)
                        rszz(a,b,i,k) = c13(i,k)*dvx_dx(a,b,i,k) + &
                            c33(i,k)*dvz_dz(a,b,i,k) + &
                            c35(i,k)*(dvz_dx(a,b,i,k) + dvx_dz(a,b,i,k)) - &
                            gamma_z(i,k)*sigmazz(a,b,i,k)
                        rsxz(a,b,i,k) = c15(i,k)*dvx_dx(a,b,i,k) + &
                            c35(i,k)*dvz_dz(a,b,i,k) + &
                            c55(i,k)*(dvz_dx(a,b,i,k) + dvx_dz(a,b,i,k)) - &
                            gamma_xz(i,k)*sigmaxz(a,b,i,k)
                        rvx(a,b,i,k) = (dsxx_dx(a,b,i,k) + dsxz_dz(a,b,i,k)) / &
                            positive(rho(i,k), 1.0_real64)
                        rvz(a,b,i,k) = (dsxz_dx(a,b,i,k) + dszz_dz(a,b,i,k)) / &
                            positive(rho(i,k), 1.0_real64)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do

        call apply_rusanov_penalty(vx, rvx, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(vz, rvz, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(sigmaxx, rsxx, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(sigmazz, rszz, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(sigmaxz, rsxz, weights, dx, dz, max_speed)

        deallocate(dvx_dx, dvx_dz, dvz_dx, dvz_dz)
        deallocate(dsxx_dx, dsxx_dz, dszz_dx, dszz_dz, dsxz_dx, dsxz_dz)
    end subroutine rhs_seismicdg2

    subroutine rhs_seismicdg25_cell(vx, vy, vz, sxx, syy, szz, sxy, syz, sxz, &
                                    c11, c12, c13, c14, c15, c16, &
                                    c22, c23, c24, c25, c26, c33, &
                                    c34, c35, c36, c44, c45, c46, &
                                    c55, c56, c66, rho, gamma_x, gamma_y, gamma_z, &
                                    gamma_xy, gamma_yz, gamma_xz, dx, dy, dz, &
                                    rvx, rvy, rvz, rsxx, rsyy, rszz, rsxy, rsyz, rsxz)
        real(real64), intent(in) :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
        real(real64), intent(in) :: sxx(:,:,:), syy(:,:,:), szz(:,:,:)
        real(real64), intent(in) :: sxy(:,:,:), syz(:,:,:), sxz(:,:,:)
        real(real64), intent(in) :: c11(:,:), c12(:,:), c13(:,:), c14(:,:), c15(:,:), c16(:,:)
        real(real64), intent(in) :: c22(:,:), c23(:,:), c24(:,:), c25(:,:), c26(:,:), c33(:,:)
        real(real64), intent(in) :: c34(:,:), c35(:,:), c36(:,:), c44(:,:), c45(:,:), c46(:,:)
        real(real64), intent(in) :: c55(:,:), c56(:,:), c66(:,:), rho(:,:)
        real(real64), intent(in) :: gamma_x(:,:), gamma_y(:,:), gamma_z(:,:)
        real(real64), intent(in) :: gamma_xy(:,:), gamma_yz(:,:), gamma_xz(:,:)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), intent(out) :: rvx(:,:,:), rvy(:,:,:), rvz(:,:,:)
        real(real64), intent(out) :: rsxx(:,:,:), rsyy(:,:,:), rszz(:,:,:)
        real(real64), intent(out) :: rsxy(:,:,:), rsyz(:,:,:), rsxz(:,:,:)

        integer :: i, j, k, nx, ny, nz
        real(real64) :: dvx_dx, dvx_dy, dvx_dz, dvy_dx, dvy_dy, dvy_dz
        real(real64) :: dvz_dx, dvz_dy, dvz_dz, exx, eyy, ezz, eyz, exz, exy
        real(real64) :: dsxx_dx, dsxy_dy, dsxz_dz, dsxy_dx, dsyy_dy, dsyz_dz
        real(real64) :: dsxz_dx, dsyz_dy, dszz_dz, inv_rho

        nx = size(vx, 1); ny = size(vx, 2); nz = size(vx, 3)
        rvx = 0.0_real64; rvy = 0.0_real64; rvz = 0.0_real64
        rsxx = 0.0_real64; rsyy = 0.0_real64; rszz = 0.0_real64
        rsxy = 0.0_real64; rsyz = 0.0_real64; rsxz = 0.0_real64

        !$omp parallel do collapse(3) private(dvx_dx,dvx_dy,dvx_dz,dvy_dx,dvy_dy,dvy_dz, &
        !$omp& dvz_dx,dvz_dy,dvz_dz,exx,eyy,ezz,eyz,exz,exy) schedule(static)
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1
                    dvx_dx = (vx(i+1,j,k) - vx(i-1,j,k)) / (2.0_real64*dx)
                    dvx_dy = (vx(i,j+1,k) - vx(i,j-1,k)) / (2.0_real64*dy)
                    dvx_dz = (vx(i,j,k+1) - vx(i,j,k-1)) / (2.0_real64*dz)
                    dvy_dx = (vy(i+1,j,k) - vy(i-1,j,k)) / (2.0_real64*dx)
                    dvy_dy = (vy(i,j+1,k) - vy(i,j-1,k)) / (2.0_real64*dy)
                    dvy_dz = (vy(i,j,k+1) - vy(i,j,k-1)) / (2.0_real64*dz)
                    dvz_dx = (vz(i+1,j,k) - vz(i-1,j,k)) / (2.0_real64*dx)
                    dvz_dy = (vz(i,j+1,k) - vz(i,j-1,k)) / (2.0_real64*dy)
                    dvz_dz = (vz(i,j,k+1) - vz(i,j,k-1)) / (2.0_real64*dz)
                    exx = dvx_dx; eyy = dvy_dy; ezz = dvz_dz
                    eyz = dvy_dz + dvz_dy
                    exz = dvx_dz + dvz_dx
                    exy = dvx_dy + dvy_dx
                    rsxx(i,j,k) = c11(i,k)*exx + c12(i,k)*eyy + c13(i,k)*ezz + &
                        c14(i,k)*eyz + c15(i,k)*exz + c16(i,k)*exy - gamma_x(i,k)*sxx(i,j,k)
                    rsyy(i,j,k) = c12(i,k)*exx + c22(i,k)*eyy + c23(i,k)*ezz + &
                        c24(i,k)*eyz + c25(i,k)*exz + c26(i,k)*exy - gamma_y(i,k)*syy(i,j,k)
                    rszz(i,j,k) = c13(i,k)*exx + c23(i,k)*eyy + c33(i,k)*ezz + &
                        c34(i,k)*eyz + c35(i,k)*exz + c36(i,k)*exy - gamma_z(i,k)*szz(i,j,k)
                    rsyz(i,j,k) = c14(i,k)*exx + c24(i,k)*eyy + c34(i,k)*ezz + &
                        c44(i,k)*eyz + c45(i,k)*exz + c46(i,k)*exy - gamma_yz(i,k)*syz(i,j,k)
                    rsxz(i,j,k) = c15(i,k)*exx + c25(i,k)*eyy + c35(i,k)*ezz + &
                        c45(i,k)*eyz + c55(i,k)*exz + c56(i,k)*exy - gamma_xz(i,k)*sxz(i,j,k)
                    rsxy(i,j,k) = c16(i,k)*exx + c26(i,k)*eyy + c36(i,k)*ezz + &
                        c46(i,k)*eyz + c56(i,k)*exz + c66(i,k)*exy - gamma_xy(i,k)*sxy(i,j,k)
                enddo
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do collapse(3) private(dsxx_dx,dsxy_dy,dsxz_dz,dsxy_dx,dsyy_dy, &
        !$omp& dsyz_dz,dsxz_dx,dsyz_dy,dszz_dz,inv_rho) schedule(static)
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1
                    dsxx_dx = (sxx(i+1,j,k) - sxx(i-1,j,k)) / (2.0_real64*dx)
                    dsxy_dy = (sxy(i,j+1,k) - sxy(i,j-1,k)) / (2.0_real64*dy)
                    dsxz_dz = (sxz(i,j,k+1) - sxz(i,j,k-1)) / (2.0_real64*dz)
                    dsxy_dx = (sxy(i+1,j,k) - sxy(i-1,j,k)) / (2.0_real64*dx)
                    dsyy_dy = (syy(i,j+1,k) - syy(i,j-1,k)) / (2.0_real64*dy)
                    dsyz_dz = (syz(i,j,k+1) - syz(i,j,k-1)) / (2.0_real64*dz)
                    dsxz_dx = (sxz(i+1,j,k) - sxz(i-1,j,k)) / (2.0_real64*dx)
                    dsyz_dy = (syz(i,j+1,k) - syz(i,j-1,k)) / (2.0_real64*dy)
                    dszz_dz = (szz(i,j,k+1) - szz(i,j,k-1)) / (2.0_real64*dz)
                    inv_rho = 1.0_real64 / positive(rho(i,k), 1.0_real64)
                    rvx(i,j,k) = (dsxx_dx + dsxy_dy + dsxz_dz) * inv_rho
                    rvy(i,j,k) = (dsxy_dx + dsyy_dy + dsyz_dz) * inv_rho
                    rvz(i,j,k) = (dsxz_dx + dsyz_dy + dszz_dz) * inv_rho
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine rhs_seismicdg25_cell

    subroutine rhs_seismicdg3_cell(vx, vy, vz, sxx, syy, szz, sxy, syz, sxz, &
                                   c11, c12, c13, c14, c15, c16, &
                                   c22, c23, c24, c25, c26, c33, &
                                   c34, c35, c36, c44, c45, c46, &
                                   c55, c56, c66, rho, gamma_x, gamma_y, gamma_z, &
                                   gamma_xy, gamma_yz, gamma_xz, dx, dy, dz, &
                                   rvx, rvy, rvz, rsxx, rsyy, rszz, rsxy, rsyz, rsxz)
        real(real64), intent(in) :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
        real(real64), intent(in) :: sxx(:,:,:), syy(:,:,:), szz(:,:,:)
        real(real64), intent(in) :: sxy(:,:,:), syz(:,:,:), sxz(:,:,:)
        real(real64), intent(in) :: c11(:,:,:), c12(:,:,:), c13(:,:,:), c14(:,:,:)
        real(real64), intent(in) :: c15(:,:,:), c16(:,:,:), c22(:,:,:), c23(:,:,:)
        real(real64), intent(in) :: c24(:,:,:), c25(:,:,:), c26(:,:,:), c33(:,:,:)
        real(real64), intent(in) :: c34(:,:,:), c35(:,:,:), c36(:,:,:), c44(:,:,:)
        real(real64), intent(in) :: c45(:,:,:), c46(:,:,:), c55(:,:,:), c56(:,:,:), c66(:,:,:)
        real(real64), intent(in) :: rho(:,:,:), gamma_x(:,:,:), gamma_y(:,:,:), gamma_z(:,:,:)
        real(real64), intent(in) :: gamma_xy(:,:,:), gamma_yz(:,:,:), gamma_xz(:,:,:)
        real(real64), intent(in) :: dx, dy, dz
        real(real64), intent(out) :: rvx(:,:,:), rvy(:,:,:), rvz(:,:,:)
        real(real64), intent(out) :: rsxx(:,:,:), rsyy(:,:,:), rszz(:,:,:)
        real(real64), intent(out) :: rsxy(:,:,:), rsyz(:,:,:), rsxz(:,:,:)

        integer :: i, j, k, nx, ny, nz
        real(real64) :: dvx_dx, dvx_dy, dvx_dz, dvy_dx, dvy_dy, dvy_dz
        real(real64) :: dvz_dx, dvz_dy, dvz_dz, exx, eyy, ezz, eyz, exz, exy
        real(real64) :: dsxx_dx, dsxy_dy, dsxz_dz, dsxy_dx, dsyy_dy, dsyz_dz
        real(real64) :: dsxz_dx, dsyz_dy, dszz_dz, inv_rho

        nx = size(vx, 1); ny = size(vx, 2); nz = size(vx, 3)
        rvx = 0.0_real64; rvy = 0.0_real64; rvz = 0.0_real64
        rsxx = 0.0_real64; rsyy = 0.0_real64; rszz = 0.0_real64
        rsxy = 0.0_real64; rsyz = 0.0_real64; rsxz = 0.0_real64

        !$omp parallel do collapse(3) private(dvx_dx,dvx_dy,dvx_dz,dvy_dx,dvy_dy,dvy_dz, &
        !$omp& dvz_dx,dvz_dy,dvz_dz,exx,eyy,ezz,eyz,exz,exy) schedule(static)
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1
                    dvx_dx = (vx(i+1,j,k) - vx(i-1,j,k)) / (2.0_real64*dx)
                    dvx_dy = (vx(i,j+1,k) - vx(i,j-1,k)) / (2.0_real64*dy)
                    dvx_dz = (vx(i,j,k+1) - vx(i,j,k-1)) / (2.0_real64*dz)
                    dvy_dx = (vy(i+1,j,k) - vy(i-1,j,k)) / (2.0_real64*dx)
                    dvy_dy = (vy(i,j+1,k) - vy(i,j-1,k)) / (2.0_real64*dy)
                    dvy_dz = (vy(i,j,k+1) - vy(i,j,k-1)) / (2.0_real64*dz)
                    dvz_dx = (vz(i+1,j,k) - vz(i-1,j,k)) / (2.0_real64*dx)
                    dvz_dy = (vz(i,j+1,k) - vz(i,j-1,k)) / (2.0_real64*dy)
                    dvz_dz = (vz(i,j,k+1) - vz(i,j,k-1)) / (2.0_real64*dz)
                    exx = dvx_dx; eyy = dvy_dy; ezz = dvz_dz
                    eyz = dvy_dz + dvz_dy
                    exz = dvx_dz + dvz_dx
                    exy = dvx_dy + dvy_dx
                    rsxx(i,j,k) = c11(i,j,k)*exx + c12(i,j,k)*eyy + c13(i,j,k)*ezz + &
                        c14(i,j,k)*eyz + c15(i,j,k)*exz + c16(i,j,k)*exy - gamma_x(i,j,k)*sxx(i,j,k)
                    rsyy(i,j,k) = c12(i,j,k)*exx + c22(i,j,k)*eyy + c23(i,j,k)*ezz + &
                        c24(i,j,k)*eyz + c25(i,j,k)*exz + c26(i,j,k)*exy - gamma_y(i,j,k)*syy(i,j,k)
                    rszz(i,j,k) = c13(i,j,k)*exx + c23(i,j,k)*eyy + c33(i,j,k)*ezz + &
                        c34(i,j,k)*eyz + c35(i,j,k)*exz + c36(i,j,k)*exy - gamma_z(i,j,k)*szz(i,j,k)
                    rsyz(i,j,k) = c14(i,j,k)*exx + c24(i,j,k)*eyy + c34(i,j,k)*ezz + &
                        c44(i,j,k)*eyz + c45(i,j,k)*exz + c46(i,j,k)*exy - gamma_yz(i,j,k)*syz(i,j,k)
                    rsxz(i,j,k) = c15(i,j,k)*exx + c25(i,j,k)*eyy + c35(i,j,k)*ezz + &
                        c45(i,j,k)*eyz + c55(i,j,k)*exz + c56(i,j,k)*exy - gamma_xz(i,j,k)*sxz(i,j,k)
                    rsxy(i,j,k) = c16(i,j,k)*exx + c26(i,j,k)*eyy + c36(i,j,k)*ezz + &
                        c46(i,j,k)*eyz + c56(i,j,k)*exz + c66(i,j,k)*exy - gamma_xy(i,j,k)*sxy(i,j,k)
                enddo
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do collapse(3) private(dsxx_dx,dsxy_dy,dsxz_dz,dsxy_dx,dsyy_dy, &
        !$omp& dsyz_dz,dsxz_dx,dsyz_dy,dszz_dz,inv_rho) schedule(static)
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1
                    dsxx_dx = (sxx(i+1,j,k) - sxx(i-1,j,k)) / (2.0_real64*dx)
                    dsxy_dy = (sxy(i,j+1,k) - sxy(i,j-1,k)) / (2.0_real64*dy)
                    dsxz_dz = (sxz(i,j,k+1) - sxz(i,j,k-1)) / (2.0_real64*dz)
                    dsxy_dx = (sxy(i+1,j,k) - sxy(i-1,j,k)) / (2.0_real64*dx)
                    dsyy_dy = (syy(i,j+1,k) - syy(i,j-1,k)) / (2.0_real64*dy)
                    dsyz_dz = (syz(i,j,k+1) - syz(i,j,k-1)) / (2.0_real64*dz)
                    dsxz_dx = (sxz(i+1,j,k) - sxz(i-1,j,k)) / (2.0_real64*dx)
                    dsyz_dy = (syz(i,j+1,k) - syz(i,j-1,k)) / (2.0_real64*dy)
                    dszz_dz = (szz(i,j,k+1) - szz(i,j,k-1)) / (2.0_real64*dz)
                    inv_rho = 1.0_real64 / positive(rho(i,j,k), 1.0_real64)
                    rvx(i,j,k) = (dsxx_dx + dsxy_dy + dsxz_dz) * inv_rho
                    rvy(i,j,k) = (dsxy_dx + dsyy_dy + dsyz_dz) * inv_rho
                    rvz(i,j,k) = (dsxz_dx + dsyz_dy + dszz_dz) * inv_rho
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine rhs_seismicdg3_cell

    subroutine zero_seismic_boundaries(vx, vy, vz, sxx, syy, szz, sxy, syz, sxz)
        real(real64), intent(inout) :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
        real(real64), intent(inout) :: sxx(:,:,:), syy(:,:,:), szz(:,:,:)
        real(real64), intent(inout) :: sxy(:,:,:), syz(:,:,:), sxz(:,:,:)
        integer :: nx, ny, nz

        nx = size(vx, 1); ny = size(vx, 2); nz = size(vx, 3)
        vx(1,:,:) = 0.0_real64; vx(nx,:,:) = 0.0_real64
        vx(:,1,:) = 0.0_real64; vx(:,ny,:) = 0.0_real64
        vx(:,:,1) = 0.0_real64; vx(:,:,nz) = 0.0_real64
        vy(1,:,:) = 0.0_real64; vy(nx,:,:) = 0.0_real64
        vy(:,1,:) = 0.0_real64; vy(:,ny,:) = 0.0_real64
        vy(:,:,1) = 0.0_real64; vy(:,:,nz) = 0.0_real64
        vz(1,:,:) = 0.0_real64; vz(nx,:,:) = 0.0_real64
        vz(:,1,:) = 0.0_real64; vz(:,ny,:) = 0.0_real64
        vz(:,:,1) = 0.0_real64; vz(:,:,nz) = 0.0_real64
        sxx(1,:,:) = 0.0_real64; sxx(nx,:,:) = 0.0_real64
        sxx(:,1,:) = 0.0_real64; sxx(:,ny,:) = 0.0_real64
        sxx(:,:,1) = 0.0_real64; sxx(:,:,nz) = 0.0_real64
        syy(1,:,:) = 0.0_real64; syy(nx,:,:) = 0.0_real64
        syy(:,1,:) = 0.0_real64; syy(:,ny,:) = 0.0_real64
        syy(:,:,1) = 0.0_real64; syy(:,:,nz) = 0.0_real64
        szz(1,:,:) = 0.0_real64; szz(nx,:,:) = 0.0_real64
        szz(:,1,:) = 0.0_real64; szz(:,ny,:) = 0.0_real64
        szz(:,:,1) = 0.0_real64; szz(:,:,nz) = 0.0_real64
        sxy(1,:,:) = 0.0_real64; sxy(nx,:,:) = 0.0_real64
        sxy(:,1,:) = 0.0_real64; sxy(:,ny,:) = 0.0_real64
        sxy(:,:,1) = 0.0_real64; sxy(:,:,nz) = 0.0_real64
        syz(1,:,:) = 0.0_real64; syz(nx,:,:) = 0.0_real64
        syz(:,1,:) = 0.0_real64; syz(:,ny,:) = 0.0_real64
        syz(:,:,1) = 0.0_real64; syz(:,:,nz) = 0.0_real64
        sxz(1,:,:) = 0.0_real64; sxz(nx,:,:) = 0.0_real64
        sxz(:,1,:) = 0.0_real64; sxz(:,ny,:) = 0.0_real64
        sxz(:,:,1) = 0.0_real64; sxz(:,:,nz) = 0.0_real64
    end subroutine zero_seismic_boundaries

    subroutine seismicdg2(domain, source, SINGLE_OUTPUT)
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        integer :: nx, nz, it, isource, ksource
        real(real64) :: dx, dz, dt, max_speed, rho0
        logical :: SINGLE
        real(real64) :: nodes(0:DG_P), weights(0:DG_P), deriv(0:DG_P,0:DG_P)
        real(real64), allocatable :: c11(:,:), c13(:,:), c15(:,:), c33(:,:), c35(:,:), c55(:,:)
        real(real64), allocatable :: rho(:,:), gamma_x(:,:), gamma_z(:,:), gamma_xz(:,:)
        real(real64), allocatable :: cell_vx(:,:), cell_vz(:,:), cell_out(:,:)
        real(real64), allocatable :: vx(:,:,:,:), vz(:,:,:,:), sxx(:,:,:,:), szz(:,:,:,:), sxz(:,:,:,:)
        real(real64), allocatable :: k1vx(:,:,:,:), k1vz(:,:,:,:), k1sxx(:,:,:,:), k1szz(:,:,:,:), k1sxz(:,:,:,:)
        real(real64), allocatable :: k2vx(:,:,:,:), k2vz(:,:,:,:), k2sxx(:,:,:,:), k2szz(:,:,:,:), k2sxz(:,:,:,:)
        real(real64), allocatable :: k3vx(:,:,:,:), k3vz(:,:,:,:), k3sxx(:,:,:,:), k3szz(:,:,:,:), k3sxz(:,:,:,:)
        real(real64), allocatable :: k4vx(:,:,:,:), k4vz(:,:,:,:), k4sxx(:,:,:,:), k4szz(:,:,:,:), k4sxz(:,:,:,:)
        real(real64), allocatable :: tvx(:,:,:,:), tvz(:,:,:,:), tsxx(:,:,:,:), tszz(:,:,:,:), tsxz(:,:,:,:)
        real(real64), allocatable :: srcx(:), srcz(:), srcxx(:), srcxz(:), srczz(:), velocnorm(:)

        SINGLE = .TRUE.; if (present(SINGLE_OUTPUT)) SINGLE = SINGLE_OUTPUT
        nx = domain%nx; nz = domain%nz; dx = domain%dx; dz = domain%dz; dt = source%dt
        call initialize_gll4(nodes, weights, deriv)

        allocate(c11(nx,nz), c13(nx,nz), c15(nx,nz), c33(nx,nz), c35(nx,nz), c55(nx,nz))
        allocate(rho(nx,nz), gamma_x(nx,nz), gamma_z(nx,nz), gamma_xz(nx,nz))
        allocate(cell_vx(nx,nz), cell_vz(nx,nz), cell_out(nx,nz))
        allocate(vx(0:DG_P,0:DG_P,nx,nz), vz(0:DG_P,0:DG_P,nx,nz))
        allocate(sxx(0:DG_P,0:DG_P,nx,nz), szz(0:DG_P,0:DG_P,nx,nz), sxz(0:DG_P,0:DG_P,nx,nz))
        allocate(k1vx(0:DG_P,0:DG_P,nx,nz), k1vz(0:DG_P,0:DG_P,nx,nz))
        allocate(k1sxx(0:DG_P,0:DG_P,nx,nz), k1szz(0:DG_P,0:DG_P,nx,nz), k1sxz(0:DG_P,0:DG_P,nx,nz))
        allocate(k2vx, source=k1vx); allocate(k2vz, source=k1vz)
        allocate(k2sxx, source=k1sxx); allocate(k2szz, source=k1szz); allocate(k2sxz, source=k1sxz)
        allocate(k3vx, source=k1vx); allocate(k3vz, source=k1vz)
        allocate(k3sxx, source=k1sxx); allocate(k3szz, source=k1szz); allocate(k3sxz, source=k1sxz)
        allocate(k4vx, source=k1vx); allocate(k4vz, source=k1vz)
        allocate(k4sxx, source=k1sxx); allocate(k4szz, source=k1szz); allocate(k4sxz, source=k1sxz)
        allocate(tvx, source=k1vx); allocate(tvz, source=k1vz)
        allocate(tsxx, source=k1sxx); allocate(tszz, source=k1szz); allocate(tsxz, source=k1sxz)
        allocate(srcx(source%time_steps), srcz(source%time_steps), srcxx(source%time_steps))
        allocate(srcxz(source%time_steps), srczz(source%time_steps), velocnorm(source%time_steps))

        call material_rw2('c11.dat', c11, .TRUE.); call material_rw2('c13.dat', c13, .TRUE.)
        call material_rw2('c15.dat', c15, .TRUE.); call material_rw2('c33.dat', c33, .TRUE.)
        call material_rw2('c35.dat', c35, .TRUE.); call material_rw2('c55.dat', c55, .TRUE.)
        call material_rw2('rho.dat', rho, .TRUE.)
        call material_rw2('gamma_x.dat', gamma_x, .TRUE.); call material_rw2('gamma_z.dat', gamma_z, .TRUE.)
        call material_rw2('gamma_xz.dat', gamma_xz, .TRUE.)
        call material_rw2('initialconditionVx.dat', cell_vx, .TRUE.)
        call material_rw2('initialconditionVz.dat', cell_vz, .TRUE.)

        call fill_modal_from_cell(cell_vx, vx)
        call fill_modal_from_cell(cell_vz, vz)
        sxx = 0.0_real64; szz = 0.0_real64; sxz = 0.0_real64
        srcx = 0.0_real64; srcz = 0.0_real64; srcxx = 0.0_real64
        srcxz = 0.0_real64; srczz = 0.0_real64; velocnorm = 0.0_real64
        isource = source%xind + domain%cpml
        ksource = source%zind + domain%cpml
        rho0 = positive(rho(isource,ksource), 1.0_real64)
        max_speed = sqrt(maxval(max(c11, c33)) / positive(minval(rho), 1.0_real64))

        if (source%source_type == 'pw') then
            stop 'Seismic DG 2D plane-wave source is not implemented'
        else if (source%source_type == 'ac') then
            call loadsource('seismicsourcex.dat', source%time_steps, srcx)
            call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        else
            call loadsource('seismicsourcexx.dat', source%time_steps, srcxx)
            call loadsource('seismicsourcexz.dat', source%time_steps, srcxz)
            call loadsource('seismicsourcezz.dat', source%time_steps, srczz)
        endif

        do it = 1, source%time_steps
            call rhs_seismicdg2(vx,vz,sxx,szz,sxz,c11,c13,c15,c33,c35,c55,rho, &
                gamma_x,gamma_z,gamma_xz,deriv,weights,dx,dz,max_speed,k1vx,k1vz,k1sxx,k1szz,k1sxz)
            tvx=vx+0.5_real64*dt*k1vx; tvz=vz+0.5_real64*dt*k1vz
            tsxx=sxx+0.5_real64*dt*k1sxx; tszz=szz+0.5_real64*dt*k1szz; tsxz=sxz+0.5_real64*dt*k1sxz
            call rhs_seismicdg2(tvx,tvz,tsxx,tszz,tsxz,c11,c13,c15,c33,c35,c55,rho, &
                gamma_x,gamma_z,gamma_xz,deriv,weights,dx,dz,max_speed,k2vx,k2vz,k2sxx,k2szz,k2sxz)
            tvx=vx+0.5_real64*dt*k2vx; tvz=vz+0.5_real64*dt*k2vz
            tsxx=sxx+0.5_real64*dt*k2sxx; tszz=szz+0.5_real64*dt*k2szz; tsxz=sxz+0.5_real64*dt*k2sxz
            call rhs_seismicdg2(tvx,tvz,tsxx,tszz,tsxz,c11,c13,c15,c33,c35,c55,rho, &
                gamma_x,gamma_z,gamma_xz,deriv,weights,dx,dz,max_speed,k3vx,k3vz,k3sxx,k3szz,k3sxz)
            tvx=vx+dt*k3vx; tvz=vz+dt*k3vz; tsxx=sxx+dt*k3sxx; tszz=szz+dt*k3szz; tsxz=sxz+dt*k3sxz
            call rhs_seismicdg2(tvx,tvz,tsxx,tszz,tsxz,c11,c13,c15,c33,c35,c55,rho, &
                gamma_x,gamma_z,gamma_xz,deriv,weights,dx,dz,max_speed,k4vx,k4vz,k4sxx,k4szz,k4sxz)

            vx = vx + dt*(k1vx + 2.0_real64*k2vx + 2.0_real64*k3vx + k4vx)/6.0_real64
            vz = vz + dt*(k1vz + 2.0_real64*k2vz + 2.0_real64*k3vz + k4vz)/6.0_real64
            sxx = sxx + dt*(k1sxx + 2.0_real64*k2sxx + 2.0_real64*k3sxx + k4sxx)/6.0_real64
            szz = szz + dt*(k1szz + 2.0_real64*k2szz + 2.0_real64*k3szz + k4szz)/6.0_real64
            sxz = sxz + dt*(k1sxz + 2.0_real64*k2sxz + 2.0_real64*k3sxz + k4sxz)/6.0_real64

            sxx(:,:,isource,ksource) = sxx(:,:,isource,ksource) + srcxx(it) / rho0
            sxz(:,:,isource,ksource) = sxz(:,:,isource,ksource) + srcxz(it) / rho0
            szz(:,:,isource,ksource) = szz(:,:,isource,ksource) + srczz(it) / rho0
            vx(:,:,isource,ksource) = vx(:,:,isource,ksource) + dt*srcx(it) / rho0
            vz(:,:,isource,ksource) = vz(:,:,isource,ksource) + dt*srcz(it) / rho0

            call cell_average(vx, weights, cell_out)
            call write_image2(cell_out, nx, nz, source, it, 'Vx', SINGLE)
            velocnorm(it) = maxval(abs(cell_out))
            call cell_average(vz, weights, cell_out)
            call write_image2(cell_out, nx, nz, source, it, 'Vz', SINGLE)
            velocnorm(it) = max(velocnorm(it), maxval(abs(cell_out)))
            if (velocnorm(it) > stability_threshold) stop 'Seismic DG 2D solution became unstable'
        enddo

        call write_array('velocity_norm.dat', source%time_steps, velocnorm)
        deallocate(c11, c13, c15, c33, c35, c55, rho, gamma_x, gamma_z, gamma_xz)
        deallocate(cell_vx, cell_vz, cell_out, vx, vz, sxx, szz, sxz)
        deallocate(k1vx, k1vz, k1sxx, k1szz, k1sxz, k2vx, k2vz, k2sxx, k2szz, k2sxz)
        deallocate(k3vx, k3vz, k3sxx, k3szz, k3sxz, k4vx, k4vz, k4sxx, k4szz, k4sxz)
        deallocate(tvx, tvz, tsxx, tszz, tsxz, srcx, srcz, srcxx, srcxz, srczz, velocnorm)
    end subroutine seismicdg2

    subroutine seismicdg25(domain, source, SINGLE_OUTPUT)
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        real(real64), allocatable :: c11(:,:), c12(:,:), c13(:,:), c14(:,:), c15(:,:), c16(:,:)
        real(real64), allocatable :: c22(:,:), c23(:,:), c24(:,:), c25(:,:), c26(:,:), c33(:,:)
        real(real64), allocatable :: c34(:,:), c35(:,:), c36(:,:), c44(:,:), c45(:,:), c46(:,:)
        real(real64), allocatable :: c55(:,:), c56(:,:), c66(:,:), rho(:,:)
        real(real64), allocatable :: gx(:,:), gy(:,:), gz(:,:), gxy(:,:), gyz(:,:), gxz(:,:)
        real(real64), allocatable :: vx(:,:,:), vy(:,:,:), vz(:,:,:), sxx(:,:,:), syy(:,:,:), szz(:,:,:)
        real(real64), allocatable :: sxy(:,:,:), syz(:,:,:), sxz(:,:,:), k1(:,:,:,:), k2(:,:,:,:), k3(:,:,:,:), k4(:,:,:,:)
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:), srcxx(:), srcxy(:), srcxz(:)
        real(real64), allocatable :: srcyy(:), srcyz(:), srczz(:), velocnorm(:)
        integer :: nx, ny, nz, it, isource, jsource, ksource
        real(real64) :: dx, dy, dz, dt, rho0
        logical :: SINGLE

        SINGLE = .TRUE.; if (present(SINGLE_OUTPUT)) SINGLE = SINGLE_OUTPUT
        nx = domain%nx; ny = domain%ny; nz = domain%nz
        dx = domain%dx; dy = domain%dy; dz = domain%dz; dt = source%dt
        allocate(c11(nx,nz), c12(nx,nz), c13(nx,nz), c14(nx,nz), c15(nx,nz), c16(nx,nz))
        allocate(c22(nx,nz), c23(nx,nz), c24(nx,nz), c25(nx,nz), c26(nx,nz), c33(nx,nz))
        allocate(c34(nx,nz), c35(nx,nz), c36(nx,nz), c44(nx,nz), c45(nx,nz), c46(nx,nz))
        allocate(c55(nx,nz), c56(nx,nz), c66(nx,nz), rho(nx,nz))
        allocate(gx(nx,nz), gy(nx,nz), gz(nx,nz), gxy(nx,nz), gyz(nx,nz), gxz(nx,nz))
        allocate(vx(nx,ny,nz), vy(nx,ny,nz), vz(nx,ny,nz), sxx(nx,ny,nz), syy(nx,ny,nz), szz(nx,ny,nz))
        allocate(sxy(nx,ny,nz), syz(nx,ny,nz), sxz(nx,ny,nz), k1(nx,ny,nz,9), k2(nx,ny,nz,9))
        allocate(k3(nx,ny,nz,9), k4(nx,ny,nz,9))
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps))
        allocate(srcxx(source%time_steps), srcxy(source%time_steps), srcxz(source%time_steps))
        allocate(srcyy(source%time_steps), srcyz(source%time_steps), srczz(source%time_steps))
        allocate(velocnorm(source%time_steps))
        call read_seismic_materials2(c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
            c34,c35,c36,c44,c45,c46,c55,c56,c66,rho,gx,gy,gz,gxy,gyz,gxz)
        call material_rw3('initialconditionVx.dat', vx, .TRUE.)
        call material_rw3('initialconditionVy.dat', vy, .TRUE.)
        call material_rw3('initialconditionVz.dat', vz, .TRUE.)
        sxx = 0.0_real64; syy = 0.0_real64; szz = 0.0_real64
        sxy = 0.0_real64; syz = 0.0_real64; sxz = 0.0_real64; velocnorm = 0.0_real64
        call zero_source_arrays(srcx,srcy,srcz,srcxx,srcxy,srcxz,srcyy,srcyz,srczz)
        if (source%source_type == 'pw') stop 'Seismic DG 2.5D plane-wave source is not implemented'
        call load_seismic_sources(source, srcx,srcy,srcz,srcxx,srcxy,srcxz,srcyy,srcyz,srczz)
        isource = source%xind + domain%cpml; jsource = source%yind + domain%cpml
        ksource = source%zind + domain%cpml
        rho0 = positive(rho(isource,ksource), 1.0_real64)

        do it = 1, source%time_steps
            call rhs_seismicdg25_cell(vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,c11,c12,c13,c14,c15,c16, &
                c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,rho,gx,gy,gz, &
                gxy,gyz,gxz,dx,dy,dz,k1(:,:,:,1),k1(:,:,:,2),k1(:,:,:,3),k1(:,:,:,4), &
                k1(:,:,:,5),k1(:,:,:,6),k1(:,:,:,7),k1(:,:,:,8),k1(:,:,:,9))
            call rhs_seismicdg25_cell(vx+0.5_real64*dt*k1(:,:,:,1),vy+0.5_real64*dt*k1(:,:,:,2), &
                vz+0.5_real64*dt*k1(:,:,:,3),sxx+0.5_real64*dt*k1(:,:,:,4), &
                syy+0.5_real64*dt*k1(:,:,:,5),szz+0.5_real64*dt*k1(:,:,:,6), &
                sxy+0.5_real64*dt*k1(:,:,:,7),syz+0.5_real64*dt*k1(:,:,:,8), &
                sxz+0.5_real64*dt*k1(:,:,:,9),c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                c34,c35,c36,c44,c45,c46,c55,c56,c66,rho,gx,gy,gz,gxy,gyz,gxz,dx,dy,dz, &
                k2(:,:,:,1),k2(:,:,:,2),k2(:,:,:,3),k2(:,:,:,4),k2(:,:,:,5),k2(:,:,:,6), &
                k2(:,:,:,7),k2(:,:,:,8),k2(:,:,:,9))
            call rhs_seismicdg25_cell(vx+0.5_real64*dt*k2(:,:,:,1),vy+0.5_real64*dt*k2(:,:,:,2), &
                vz+0.5_real64*dt*k2(:,:,:,3),sxx+0.5_real64*dt*k2(:,:,:,4), &
                syy+0.5_real64*dt*k2(:,:,:,5),szz+0.5_real64*dt*k2(:,:,:,6), &
                sxy+0.5_real64*dt*k2(:,:,:,7),syz+0.5_real64*dt*k2(:,:,:,8), &
                sxz+0.5_real64*dt*k2(:,:,:,9),c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                c34,c35,c36,c44,c45,c46,c55,c56,c66,rho,gx,gy,gz,gxy,gyz,gxz,dx,dy,dz, &
                k3(:,:,:,1),k3(:,:,:,2),k3(:,:,:,3),k3(:,:,:,4),k3(:,:,:,5),k3(:,:,:,6), &
                k3(:,:,:,7),k3(:,:,:,8),k3(:,:,:,9))
            call rhs_seismicdg25_cell(vx+dt*k3(:,:,:,1),vy+dt*k3(:,:,:,2),vz+dt*k3(:,:,:,3), &
                sxx+dt*k3(:,:,:,4),syy+dt*k3(:,:,:,5),szz+dt*k3(:,:,:,6),sxy+dt*k3(:,:,:,7), &
                syz+dt*k3(:,:,:,8),sxz+dt*k3(:,:,:,9),c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,rho,gx,gy,gz,gxy,gyz,gxz,dx,dy,dz, &
                k4(:,:,:,1),k4(:,:,:,2),k4(:,:,:,3),k4(:,:,:,4),k4(:,:,:,5),k4(:,:,:,6), &
                k4(:,:,:,7),k4(:,:,:,8),k4(:,:,:,9))
            call update_seismic3(vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,k1,k2,k3,k4,dt)
            call inject_seismic3(vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,isource,jsource,ksource, &
                rho0,dt,it,srcx,srcy,srcz,srcxx,srcxy,srcxz,srcyy,srcyz,srczz)
            call zero_seismic_boundaries(vx,vy,vz,sxx,syy,szz,sxy,syz,sxz)
            call write_image3(vx, nx, ny, nz, source, it, 'Vx', SINGLE)
            call write_image3(vy, nx, ny, nz, source, it, 'Vy', SINGLE)
            call write_image3(vz, nx, ny, nz, source, it, 'Vz', SINGLE)
            velocnorm(it) = max(maxval(abs(vx)), max(maxval(abs(vy)), maxval(abs(vz))))
            if (velocnorm(it) > stability_threshold) stop 'Seismic DG 2.5D solution became unstable'
        enddo
        call write_array('velocity_norm.dat', source%time_steps, velocnorm)
        deallocate(c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
        deallocate(rho,gx,gy,gz,gxy,gyz,gxz,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,k1,k2,k3,k4)
        deallocate(srcx,srcy,srcz,srcxx,srcxy,srcxz,srcyy,srcyz,srczz,velocnorm)
    end subroutine seismicdg25

    subroutine seismicdg3(domain, source, SINGLE_OUTPUT)
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        real(real64), allocatable :: c11(:,:,:), c12(:,:,:), c13(:,:,:), c14(:,:,:), c15(:,:,:), c16(:,:,:)
        real(real64), allocatable :: c22(:,:,:), c23(:,:,:), c24(:,:,:), c25(:,:,:), c26(:,:,:), c33(:,:,:)
        real(real64), allocatable :: c34(:,:,:), c35(:,:,:), c36(:,:,:), c44(:,:,:), c45(:,:,:), c46(:,:,:)
        real(real64), allocatable :: c55(:,:,:), c56(:,:,:), c66(:,:,:), rho(:,:,:)
        real(real64), allocatable :: gx(:,:,:), gy(:,:,:), gz(:,:,:), gxy(:,:,:), gyz(:,:,:), gxz(:,:,:)
        real(real64), allocatable :: vx(:,:,:), vy(:,:,:), vz(:,:,:), sxx(:,:,:), syy(:,:,:), szz(:,:,:)
        real(real64), allocatable :: sxy(:,:,:), syz(:,:,:), sxz(:,:,:), k1(:,:,:,:), k2(:,:,:,:), k3(:,:,:,:), k4(:,:,:,:)
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:), srcxx(:), srcxy(:), srcxz(:)
        real(real64), allocatable :: srcyy(:), srcyz(:), srczz(:), velocnorm(:)
        integer :: nx, ny, nz, it, isource, jsource, ksource
        real(real64) :: dx, dy, dz, dt, rho0
        logical :: SINGLE

        SINGLE = .TRUE.; if (present(SINGLE_OUTPUT)) SINGLE = SINGLE_OUTPUT
        nx = domain%nx; ny = domain%ny; nz = domain%nz
        dx = domain%dx; dy = domain%dy; dz = domain%dz; dt = source%dt
        allocate(c11(nx,ny,nz), c12(nx,ny,nz), c13(nx,ny,nz), c14(nx,ny,nz), c15(nx,ny,nz), c16(nx,ny,nz))
        allocate(c22(nx,ny,nz), c23(nx,ny,nz), c24(nx,ny,nz), c25(nx,ny,nz), c26(nx,ny,nz), c33(nx,ny,nz))
        allocate(c34(nx,ny,nz), c35(nx,ny,nz), c36(nx,ny,nz), c44(nx,ny,nz), c45(nx,ny,nz), c46(nx,ny,nz))
        allocate(c55(nx,ny,nz), c56(nx,ny,nz), c66(nx,ny,nz), rho(nx,ny,nz))
        allocate(gx(nx,ny,nz), gy(nx,ny,nz), gz(nx,ny,nz), gxy(nx,ny,nz), gyz(nx,ny,nz), gxz(nx,ny,nz))
        allocate(vx(nx,ny,nz), vy(nx,ny,nz), vz(nx,ny,nz), sxx(nx,ny,nz), syy(nx,ny,nz), szz(nx,ny,nz))
        allocate(sxy(nx,ny,nz), syz(nx,ny,nz), sxz(nx,ny,nz), k1(nx,ny,nz,9), k2(nx,ny,nz,9))
        allocate(k3(nx,ny,nz,9), k4(nx,ny,nz,9))
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps))
        allocate(srcxx(source%time_steps), srcxy(source%time_steps), srcxz(source%time_steps))
        allocate(srcyy(source%time_steps), srcyz(source%time_steps), srczz(source%time_steps))
        allocate(velocnorm(source%time_steps))
        call read_seismic_materials3(c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
            c34,c35,c36,c44,c45,c46,c55,c56,c66,rho,gx,gy,gz,gxy,gyz,gxz)
        call material_rw3('initialconditionVx.dat', vx, .TRUE.)
        call material_rw3('initialconditionVy.dat', vy, .TRUE.)
        call material_rw3('initialconditionVz.dat', vz, .TRUE.)
        sxx = 0.0_real64; syy = 0.0_real64; szz = 0.0_real64
        sxy = 0.0_real64; syz = 0.0_real64; sxz = 0.0_real64; velocnorm = 0.0_real64
        call zero_source_arrays(srcx,srcy,srcz,srcxx,srcxy,srcxz,srcyy,srcyz,srczz)
        if (source%source_type == 'pw') stop 'Seismic DG 3D plane-wave source is not implemented'
        call load_seismic_sources(source, srcx,srcy,srcz,srcxx,srcxy,srcxz,srcyy,srcyz,srczz)
        isource = source%xind + domain%cpml; jsource = source%yind + domain%cpml
        ksource = source%zind + domain%cpml
        rho0 = positive(rho(isource,jsource,ksource), 1.0_real64)

        do it = 1, source%time_steps
            call rhs_seismicdg3_cell(vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,c11,c12,c13,c14,c15,c16, &
                c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,rho,gx,gy,gz, &
                gxy,gyz,gxz,dx,dy,dz,k1(:,:,:,1),k1(:,:,:,2),k1(:,:,:,3),k1(:,:,:,4), &
                k1(:,:,:,5),k1(:,:,:,6),k1(:,:,:,7),k1(:,:,:,8),k1(:,:,:,9))
            call rhs_seismicdg3_cell(vx+0.5_real64*dt*k1(:,:,:,1),vy+0.5_real64*dt*k1(:,:,:,2), &
                vz+0.5_real64*dt*k1(:,:,:,3),sxx+0.5_real64*dt*k1(:,:,:,4), &
                syy+0.5_real64*dt*k1(:,:,:,5),szz+0.5_real64*dt*k1(:,:,:,6), &
                sxy+0.5_real64*dt*k1(:,:,:,7),syz+0.5_real64*dt*k1(:,:,:,8), &
                sxz+0.5_real64*dt*k1(:,:,:,9),c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                c34,c35,c36,c44,c45,c46,c55,c56,c66,rho,gx,gy,gz,gxy,gyz,gxz,dx,dy,dz, &
                k2(:,:,:,1),k2(:,:,:,2),k2(:,:,:,3),k2(:,:,:,4),k2(:,:,:,5),k2(:,:,:,6), &
                k2(:,:,:,7),k2(:,:,:,8),k2(:,:,:,9))
            call rhs_seismicdg3_cell(vx+0.5_real64*dt*k2(:,:,:,1),vy+0.5_real64*dt*k2(:,:,:,2), &
                vz+0.5_real64*dt*k2(:,:,:,3),sxx+0.5_real64*dt*k2(:,:,:,4), &
                syy+0.5_real64*dt*k2(:,:,:,5),szz+0.5_real64*dt*k2(:,:,:,6), &
                sxy+0.5_real64*dt*k2(:,:,:,7),syz+0.5_real64*dt*k2(:,:,:,8), &
                sxz+0.5_real64*dt*k2(:,:,:,9),c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                c34,c35,c36,c44,c45,c46,c55,c56,c66,rho,gx,gy,gz,gxy,gyz,gxz,dx,dy,dz, &
                k3(:,:,:,1),k3(:,:,:,2),k3(:,:,:,3),k3(:,:,:,4),k3(:,:,:,5),k3(:,:,:,6), &
                k3(:,:,:,7),k3(:,:,:,8),k3(:,:,:,9))
            call rhs_seismicdg3_cell(vx+dt*k3(:,:,:,1),vy+dt*k3(:,:,:,2),vz+dt*k3(:,:,:,3), &
                sxx+dt*k3(:,:,:,4),syy+dt*k3(:,:,:,5),szz+dt*k3(:,:,:,6),sxy+dt*k3(:,:,:,7), &
                syz+dt*k3(:,:,:,8),sxz+dt*k3(:,:,:,9),c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,rho,gx,gy,gz,gxy,gyz,gxz,dx,dy,dz, &
                k4(:,:,:,1),k4(:,:,:,2),k4(:,:,:,3),k4(:,:,:,4),k4(:,:,:,5),k4(:,:,:,6), &
                k4(:,:,:,7),k4(:,:,:,8),k4(:,:,:,9))
            call update_seismic3(vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,k1,k2,k3,k4,dt)
            call inject_seismic3(vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,isource,jsource,ksource, &
                rho0,dt,it,srcx,srcy,srcz,srcxx,srcxy,srcxz,srcyy,srcyz,srczz)
            call zero_seismic_boundaries(vx,vy,vz,sxx,syy,szz,sxy,syz,sxz)
            call write_image3(vx, nx, ny, nz, source, it, 'Vx', SINGLE)
            call write_image3(vy, nx, ny, nz, source, it, 'Vy', SINGLE)
            call write_image3(vz, nx, ny, nz, source, it, 'Vz', SINGLE)
            velocnorm(it) = max(maxval(abs(vx)), max(maxval(abs(vy)), maxval(abs(vz))))
            if (velocnorm(it) > stability_threshold) stop 'Seismic DG 3D solution became unstable'
        enddo
        call write_array('velocity_norm.dat', source%time_steps, velocnorm)
        deallocate(c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
        deallocate(rho,gx,gy,gz,gxy,gyz,gxz,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,k1,k2,k3,k4)
        deallocate(srcx,srcy,srcz,srcxx,srcxy,srcxz,srcyy,srcyz,srczz,velocnorm)
    end subroutine seismicdg3

    subroutine zero_source_arrays(srcx, srcy, srcz, srcxx, srcxy, srcxz, srcyy, srcyz, srczz)
        real(real64), intent(out) :: srcx(:), srcy(:), srcz(:), srcxx(:), srcxy(:), srcxz(:)
        real(real64), intent(out) :: srcyy(:), srcyz(:), srczz(:)

        srcx = 0.0_real64; srcy = 0.0_real64; srcz = 0.0_real64
        srcxx = 0.0_real64; srcxy = 0.0_real64; srcxz = 0.0_real64
        srcyy = 0.0_real64; srcyz = 0.0_real64; srczz = 0.0_real64
    end subroutine zero_source_arrays

    subroutine load_seismic_sources(source, srcx, srcy, srcz, srcxx, srcxy, srcxz, srcyy, srcyz, srczz)
        type(Source_Type), intent(in) :: source
        real(real64), intent(inout) :: srcx(:), srcy(:), srcz(:), srcxx(:), srcxy(:), srcxz(:)
        real(real64), intent(inout) :: srcyy(:), srcyz(:), srczz(:)

        if (source%source_type == 'ac') then
            call loadsource('seismicsourcex.dat', source%time_steps, srcx)
            call loadsource('seismicsourcey.dat', source%time_steps, srcy)
            call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        else
            call loadsource('seismicsourcexx.dat', source%time_steps, srcxx)
            call loadsource('seismicsourcexy.dat', source%time_steps, srcxy)
            call loadsource('seismicsourcexz.dat', source%time_steps, srcxz)
            call loadsource('seismicsourceyy.dat', source%time_steps, srcyy)
            call loadsource('seismicsourceyz.dat', source%time_steps, srcyz)
            call loadsource('seismicsourcezz.dat', source%time_steps, srczz)
        endif
    end subroutine load_seismic_sources

    subroutine update_seismic3(vx, vy, vz, sxx, syy, szz, sxy, syz, sxz, k1, k2, k3, k4, dt)
        real(real64), intent(inout) :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
        real(real64), intent(inout) :: sxx(:,:,:), syy(:,:,:), szz(:,:,:), sxy(:,:,:), syz(:,:,:), sxz(:,:,:)
        real(real64), intent(in) :: k1(:,:,:,:), k2(:,:,:,:), k3(:,:,:,:), k4(:,:,:,:), dt

        vx = vx + dt*(k1(:,:,:,1) + 2.0_real64*k2(:,:,:,1) + 2.0_real64*k3(:,:,:,1) + k4(:,:,:,1))/6.0_real64
        vy = vy + dt*(k1(:,:,:,2) + 2.0_real64*k2(:,:,:,2) + 2.0_real64*k3(:,:,:,2) + k4(:,:,:,2))/6.0_real64
        vz = vz + dt*(k1(:,:,:,3) + 2.0_real64*k2(:,:,:,3) + 2.0_real64*k3(:,:,:,3) + k4(:,:,:,3))/6.0_real64
        sxx = sxx + dt*(k1(:,:,:,4) + 2.0_real64*k2(:,:,:,4) + 2.0_real64*k3(:,:,:,4) + k4(:,:,:,4))/6.0_real64
        syy = syy + dt*(k1(:,:,:,5) + 2.0_real64*k2(:,:,:,5) + 2.0_real64*k3(:,:,:,5) + k4(:,:,:,5))/6.0_real64
        szz = szz + dt*(k1(:,:,:,6) + 2.0_real64*k2(:,:,:,6) + 2.0_real64*k3(:,:,:,6) + k4(:,:,:,6))/6.0_real64
        sxy = sxy + dt*(k1(:,:,:,7) + 2.0_real64*k2(:,:,:,7) + 2.0_real64*k3(:,:,:,7) + k4(:,:,:,7))/6.0_real64
        syz = syz + dt*(k1(:,:,:,8) + 2.0_real64*k2(:,:,:,8) + 2.0_real64*k3(:,:,:,8) + k4(:,:,:,8))/6.0_real64
        sxz = sxz + dt*(k1(:,:,:,9) + 2.0_real64*k2(:,:,:,9) + 2.0_real64*k3(:,:,:,9) + k4(:,:,:,9))/6.0_real64
    end subroutine update_seismic3

    subroutine inject_seismic3(vx, vy, vz, sxx, syy, szz, sxy, syz, sxz, isource, jsource, ksource, &
                               rho0, dt, it, srcx, srcy, srcz, srcxx, srcxy, srcxz, srcyy, srcyz, srczz)
        real(real64), intent(inout) :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
        real(real64), intent(inout) :: sxx(:,:,:), syy(:,:,:), szz(:,:,:), sxy(:,:,:), syz(:,:,:), sxz(:,:,:)
        integer, intent(in) :: isource, jsource, ksource, it
        real(real64), intent(in) :: rho0, dt
        real(real64), intent(in) :: srcx(:), srcy(:), srcz(:), srcxx(:), srcxy(:), srcxz(:)
        real(real64), intent(in) :: srcyy(:), srcyz(:), srczz(:)

        vx(isource,jsource,ksource) = vx(isource,jsource,ksource) + dt*srcx(it) / rho0
        vy(isource,jsource,ksource) = vy(isource,jsource,ksource) + dt*srcy(it) / rho0
        vz(isource,jsource,ksource) = vz(isource,jsource,ksource) + dt*srcz(it) / rho0
        sxx(isource,jsource,ksource) = sxx(isource,jsource,ksource) + srcxx(it) / rho0
        sxy(isource,jsource,ksource) = sxy(isource,jsource,ksource) + srcxy(it) / rho0
        sxz(isource,jsource,ksource) = sxz(isource,jsource,ksource) + srcxz(it) / rho0
        syy(isource,jsource,ksource) = syy(isource,jsource,ksource) + srcyy(it) / rho0
        syz(isource,jsource,ksource) = syz(isource,jsource,ksource) + srcyz(it) / rho0
        szz(isource,jsource,ksource) = szz(isource,jsource,ksource) + srczz(it) / rho0
    end subroutine inject_seismic3

    subroutine read_seismic_materials2(c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                                       c34,c35,c36,c44,c45,c46,c55,c56,c66,rho,gx,gy,gz,gxy,gyz,gxz)
        real(real64), intent(out) :: c11(:,:), c12(:,:), c13(:,:), c14(:,:), c15(:,:), c16(:,:)
        real(real64), intent(out) :: c22(:,:), c23(:,:), c24(:,:), c25(:,:), c26(:,:), c33(:,:)
        real(real64), intent(out) :: c34(:,:), c35(:,:), c36(:,:), c44(:,:), c45(:,:), c46(:,:)
        real(real64), intent(out) :: c55(:,:), c56(:,:), c66(:,:), rho(:,:)
        real(real64), intent(out) :: gx(:,:), gy(:,:), gz(:,:), gxy(:,:), gyz(:,:), gxz(:,:)

        call material_rw2('c11.dat', c11, .TRUE.); call material_rw2('c12.dat', c12, .TRUE.)
        call material_rw2('c13.dat', c13, .TRUE.); call material_rw2('c14.dat', c14, .TRUE.)
        call material_rw2('c15.dat', c15, .TRUE.); call material_rw2('c16.dat', c16, .TRUE.)
        call material_rw2('c22.dat', c22, .TRUE.); call material_rw2('c23.dat', c23, .TRUE.)
        call material_rw2('c24.dat', c24, .TRUE.); call material_rw2('c25.dat', c25, .TRUE.)
        call material_rw2('c26.dat', c26, .TRUE.); call material_rw2('c33.dat', c33, .TRUE.)
        call material_rw2('c34.dat', c34, .TRUE.); call material_rw2('c35.dat', c35, .TRUE.)
        call material_rw2('c36.dat', c36, .TRUE.); call material_rw2('c44.dat', c44, .TRUE.)
        call material_rw2('c45.dat', c45, .TRUE.); call material_rw2('c46.dat', c46, .TRUE.)
        call material_rw2('c55.dat', c55, .TRUE.); call material_rw2('c56.dat', c56, .TRUE.)
        call material_rw2('c66.dat', c66, .TRUE.); call material_rw2('rho.dat', rho, .TRUE.)
        call material_rw2('gamma_x.dat', gx, .TRUE.); call material_rw2('gamma_y.dat', gy, .TRUE.)
        call material_rw2('gamma_z.dat', gz, .TRUE.); call material_rw2('gamma_xy.dat', gxy, .TRUE.)
        call material_rw2('gamma_yz.dat', gyz, .TRUE.); call material_rw2('gamma_xz.dat', gxz, .TRUE.)
    end subroutine read_seismic_materials2

    subroutine read_seismic_materials3(c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                                       c34,c35,c36,c44,c45,c46,c55,c56,c66,rho,gx,gy,gz,gxy,gyz,gxz)
        real(real64), intent(out) :: c11(:,:,:), c12(:,:,:), c13(:,:,:), c14(:,:,:), c15(:,:,:), c16(:,:,:)
        real(real64), intent(out) :: c22(:,:,:), c23(:,:,:), c24(:,:,:), c25(:,:,:), c26(:,:,:), c33(:,:,:)
        real(real64), intent(out) :: c34(:,:,:), c35(:,:,:), c36(:,:,:), c44(:,:,:), c45(:,:,:), c46(:,:,:)
        real(real64), intent(out) :: c55(:,:,:), c56(:,:,:), c66(:,:,:), rho(:,:,:)
        real(real64), intent(out) :: gx(:,:,:), gy(:,:,:), gz(:,:,:), gxy(:,:,:), gyz(:,:,:), gxz(:,:,:)

        call material_rw3('c11.dat', c11, .TRUE.); call material_rw3('c12.dat', c12, .TRUE.)
        call material_rw3('c13.dat', c13, .TRUE.); call material_rw3('c14.dat', c14, .TRUE.)
        call material_rw3('c15.dat', c15, .TRUE.); call material_rw3('c16.dat', c16, .TRUE.)
        call material_rw3('c22.dat', c22, .TRUE.); call material_rw3('c23.dat', c23, .TRUE.)
        call material_rw3('c24.dat', c24, .TRUE.); call material_rw3('c25.dat', c25, .TRUE.)
        call material_rw3('c26.dat', c26, .TRUE.); call material_rw3('c33.dat', c33, .TRUE.)
        call material_rw3('c34.dat', c34, .TRUE.); call material_rw3('c35.dat', c35, .TRUE.)
        call material_rw3('c36.dat', c36, .TRUE.); call material_rw3('c44.dat', c44, .TRUE.)
        call material_rw3('c45.dat', c45, .TRUE.); call material_rw3('c46.dat', c46, .TRUE.)
        call material_rw3('c55.dat', c55, .TRUE.); call material_rw3('c56.dat', c56, .TRUE.)
        call material_rw3('c66.dat', c66, .TRUE.); call material_rw3('rho.dat', rho, .TRUE.)
        call material_rw3('gamma_x.dat', gx, .TRUE.); call material_rw3('gamma_y.dat', gy, .TRUE.)
        call material_rw3('gamma_z.dat', gz, .TRUE.); call material_rw3('gamma_xy.dat', gxy, .TRUE.)
        call material_rw3('gamma_yz.dat', gyz, .TRUE.); call material_rw3('gamma_xz.dat', gxz, .TRUE.)
    end subroutine read_seismic_materials3

end module seismicdg
