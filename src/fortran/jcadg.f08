module jcadg

    use iso_fortran_env, only: real64
    use seidartio
    use seidart_types
    use constants
    use discontinuous_galerkin_methods, only: DG_P, positive, initialize_gll4, &
        fill_modal_from_cell, cell_average, dg_derivatives, apply_rusanov_penalty, &
        is_pressure_injection_source

    implicit none

contains

    subroutine rhs_jca2(vx, vz, pressure, rho_x, rho_z, bulk, attenuation, &
                        deriv, weights, dx, dz, max_speed, rvx, rvz, rp)
        real(real64), intent(in) :: vx(0:,0:,:,:), vz(0:,0:,:,:), pressure(0:,0:,:,:)
        real(real64), intent(in) :: rho_x(:,:), rho_z(:,:), bulk(:,:), attenuation(:,:)
        real(real64), intent(in) :: deriv(0:DG_P,0:DG_P), weights(0:DG_P)
        real(real64), intent(in) :: dx, dz, max_speed
        real(real64), intent(out) :: rvx(0:,0:,:,:), rvz(0:,0:,:,:), rp(0:,0:,:,:)

        real(real64), allocatable :: dp_dx(:,:,:,:), dp_dz(:,:,:,:)
        real(real64), allocatable :: dvx_dx(:,:,:,:), dvz_dz(:,:,:,:)
        real(real64), allocatable :: scratch(:,:,:,:)
        integer :: a, b, i, k

        allocate(dp_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)), dp_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(dvx_dx(0:DG_P,0:DG_P,size(vx,3),size(vx,4)), dvz_dz(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))
        allocate(scratch(0:DG_P,0:DG_P,size(vx,3),size(vx,4)))

        call dg_derivatives(pressure, deriv, dx, dz, dp_dx, dp_dz)
        call dg_derivatives(vx, deriv, dx, dz, dvx_dx, scratch)
        call dg_derivatives(vz, deriv, dx, dz, scratch, dvz_dz)

        !$omp parallel do collapse(2) private(a,b) schedule(static)
        do k = 1, size(vx, 4)
            do i = 1, size(vx, 3)
                do b = 0, DG_P
                    do a = 0, DG_P
                        rvx(a,b,i,k) = -dp_dx(a,b,i,k) / positive(rho_x(i,k), 1.0e-6_real64)
                        rvz(a,b,i,k) = -dp_dz(a,b,i,k) / positive(rho_z(i,k), 1.0e-6_real64)
                        rp(a,b,i,k) = -bulk(i,k) * (dvx_dx(a,b,i,k) + dvz_dz(a,b,i,k)) - &
                                      attenuation(i,k) * pressure(a,b,i,k)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do

        call apply_rusanov_penalty(vx, rvx, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(vz, rvz, weights, dx, dz, max_speed)
        call apply_rusanov_penalty(pressure, rp, weights, dx, dz, max_speed)
    end subroutine rhs_jca2

    subroutine jcadg2(domain, source, SINGLE_OUTPUT)
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        integer :: nx, nz, it, isource, ksource
        real(real64) :: dx, dz, dt, max_speed
        logical :: SINGLE
        real(real64) :: nodes(0:DG_P), weights(0:DG_P), deriv(0:DG_P,0:DG_P)
        real(real64), allocatable :: rho_x(:,:), rho_z(:,:), bulk(:,:), attenuation(:,:)
        real(real64), allocatable :: cell_vx(:,:), cell_vz(:,:), cell_out(:,:)
        real(real64), allocatable :: vx(:,:,:,:), vz(:,:,:,:), pressure(:,:,:,:)
        real(real64), allocatable :: k1vx(:,:,:,:), k1vz(:,:,:,:), k1p(:,:,:,:)
        real(real64), allocatable :: k2vx(:,:,:,:), k2vz(:,:,:,:), k2p(:,:,:,:)
        real(real64), allocatable :: k3vx(:,:,:,:), k3vz(:,:,:,:), k3p(:,:,:,:)
        real(real64), allocatable :: k4vx(:,:,:,:), k4vz(:,:,:,:), k4p(:,:,:,:)
        real(real64), allocatable :: tvx(:,:,:,:), tvz(:,:,:,:), tp(:,:,:,:)
        real(real64), allocatable :: srcx(:), srcz(:), srcxx(:), srczz(:), srcxz(:), srcp(:)
        real(real64), allocatable :: velocnorm(:), pressurenorm(:)

        SINGLE = .TRUE.; if (present(SINGLE_OUTPUT)) SINGLE = SINGLE_OUTPUT
        nx = domain%nx; nz = domain%nz; dx = domain%dx; dz = domain%dz; dt = source%dt
        call initialize_gll4(nodes, weights, deriv)

        allocate(rho_x(nx,nz), rho_z(nx,nz), bulk(nx,nz), attenuation(nx,nz))
        allocate(cell_vx(nx,nz), cell_vz(nx,nz), cell_out(nx,nz))
        allocate(vx(0:DG_P,0:DG_P,nx,nz), vz(0:DG_P,0:DG_P,nx,nz), pressure(0:DG_P,0:DG_P,nx,nz))
        allocate(k1vx(0:DG_P,0:DG_P,nx,nz), k1vz(0:DG_P,0:DG_P,nx,nz), k1p(0:DG_P,0:DG_P,nx,nz))
        allocate(k2vx, source=k1vx); allocate(k2vz, source=k1vz); allocate(k2p, source=k1p)
        allocate(k3vx, source=k1vx); allocate(k3vz, source=k1vz); allocate(k3p, source=k1p)
        allocate(k4vx, source=k1vx); allocate(k4vz, source=k1vz); allocate(k4p, source=k1p)
        allocate(tvx, source=k1vx); allocate(tvz, source=k1vz); allocate(tp, source=k1p)
        allocate(srcx(source%time_steps), srcz(source%time_steps), srcxx(source%time_steps))
        allocate(srczz(source%time_steps), srcxz(source%time_steps))
        allocate(srcp(source%time_steps))
        allocate(velocnorm(source%time_steps), pressurenorm(source%time_steps))

        call material_rw2('jca_rho_x.dat', rho_x, .TRUE.)
        call material_rw2('jca_rho_z.dat', rho_z, .TRUE.)
        call material_rw2('jca_bulk_modulus.dat', bulk, .TRUE.)
        call material_rw2('jca_attenuation.dat', attenuation, .TRUE.)
        call material_rw2('initialconditionVx.dat', cell_vx, .TRUE.)
        call material_rw2('initialconditionVz.dat', cell_vz, .TRUE.)
        call fill_modal_from_cell(cell_vx, vx)
        call fill_modal_from_cell(cell_vz, vz)
        pressure = 0.0_real64
        srcx = 0.0_real64; srcz = 0.0_real64; srcxx = 0.0_real64; srczz = 0.0_real64; srcxz = 0.0_real64
        srcp = 0.0_real64
        velocnorm = 0.0_real64; pressurenorm = 0.0_real64
        isource = source%xind + domain%cpml; ksource = source%zind + domain%cpml
        max_speed = sqrt(maxval(bulk / positive(minval(rho_x), 1.0e-6_real64)))

        if (source%source_type == 'ac') then
            call loadsource('seismicsourcex.dat', source%time_steps, srcx)
            call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        else if (is_pressure_injection_source(source%source_type)) then
            call loadsource('seismicsourcep.dat', source%time_steps, srcp)
        else
            call loadsource('seismicsourcexx.dat', source%time_steps, srcxx)
            call loadsource('seismicsourcexz.dat', source%time_steps, srcxz)
            call loadsource('seismicsourcezz.dat', source%time_steps, srczz)
        endif

        do it = 1, source%time_steps
            call rhs_jca2(vx,vz,pressure,rho_x,rho_z,bulk,attenuation,deriv,weights,dx,dz,max_speed,k1vx,k1vz,k1p)
            tvx=vx+0.5_real64*dt*k1vx; tvz=vz+0.5_real64*dt*k1vz; tp=pressure+0.5_real64*dt*k1p
            call rhs_jca2(tvx,tvz,tp,rho_x,rho_z,bulk,attenuation,deriv,weights,dx,dz,max_speed,k2vx,k2vz,k2p)
            tvx=vx+0.5_real64*dt*k2vx; tvz=vz+0.5_real64*dt*k2vz; tp=pressure+0.5_real64*dt*k2p
            call rhs_jca2(tvx,tvz,tp,rho_x,rho_z,bulk,attenuation,deriv,weights,dx,dz,max_speed,k3vx,k3vz,k3p)
            tvx=vx+dt*k3vx; tvz=vz+dt*k3vz; tp=pressure+dt*k3p
            call rhs_jca2(tvx,tvz,tp,rho_x,rho_z,bulk,attenuation,deriv,weights,dx,dz,max_speed,k4vx,k4vz,k4p)

            vx = vx + dt*(k1vx+2.0_real64*k2vx+2.0_real64*k3vx+k4vx)/6.0_real64
            vz = vz + dt*(k1vz+2.0_real64*k2vz+2.0_real64*k3vz+k4vz)/6.0_real64
            pressure = pressure + dt*(k1p+2.0_real64*k2p+2.0_real64*k3p+k4p)/6.0_real64

            vx(0,0,isource,ksource) = vx(0,0,isource,ksource) + srcx(it) * dt / positive(rho_x(isource,ksource), 1.0e-6_real64)
            vz(0,0,isource,ksource) = vz(0,0,isource,ksource) + srcz(it) * dt / positive(rho_z(isource,ksource), 1.0e-6_real64)
            pressure(0,0,isource,ksource) = pressure(0,0,isource,ksource) + (srcxx(it) + srczz(it) + srcxz(it)) * dt
            pressure(0,0,isource,ksource) = pressure(0,0,isource,ksource) + srcp(it) * dt

            call cell_average(vx, weights, cell_out); call write_image2(cell_out, nx, nz, source, it, 'Vx', SINGLE)
            velocnorm(it) = maxval(abs(cell_out))
            call cell_average(vz, weights, cell_out); call write_image2(cell_out, nx, nz, source, it, 'Vz', SINGLE)
            velocnorm(it) = max(velocnorm(it), maxval(abs(cell_out)))
            call cell_average(pressure, weights, cell_out); call write_image2(cell_out, nx, nz, source, it, 'Pp', SINGLE)
            pressurenorm(it) = maxval(abs(cell_out))
            if (velocnorm(it) > stability_threshold) stop 'JCA DG 2D solution became unstable'
        enddo

        call write_array('velocity_norm.dat', source%time_steps, velocnorm)
        call write_array('pore_pressure_norm.dat', source%time_steps, pressurenorm)
    end subroutine jcadg2

    subroutine jcadg25(domain, source, SINGLE_OUTPUT)
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        integer :: nx, ny, nz, i, j, k, it, isource, jsource, ksource
        real(real64) :: dx, dy, dz, dt, div_v, px, py, pz
        logical :: SINGLE
        real(real64), allocatable :: rho_x(:,:), rho_y(:,:), rho_z(:,:), bulk(:,:), attenuation(:,:)
        real(real64), allocatable :: vx(:,:,:), vy(:,:,:), vz(:,:,:), pressure(:,:,:)
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:), srcxx(:), srcyy(:), srczz(:), srcp(:)
        real(real64), allocatable :: velocnorm(:), pressurenorm(:)

        SINGLE = .TRUE.; if (present(SINGLE_OUTPUT)) SINGLE = SINGLE_OUTPUT
        nx=domain%nx; ny=domain%ny; nz=domain%nz; dx=domain%dx; dy=domain%dy; dz=domain%dz; dt=source%dt
        allocate(rho_x(nx,nz), rho_y(nx,nz), rho_z(nx,nz), bulk(nx,nz), attenuation(nx,nz))
        allocate(vx(nx,ny,nz), vy(nx,ny,nz), vz(nx,ny,nz), pressure(nx,ny,nz))
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps))
        allocate(srcxx(source%time_steps), srcyy(source%time_steps), srczz(source%time_steps))
        allocate(srcp(source%time_steps))
        allocate(velocnorm(source%time_steps), pressurenorm(source%time_steps))

        call material_rw2('jca_rho_x.dat', rho_x, .TRUE.)
        call material_rw2('jca_rho_y.dat', rho_y, .TRUE.)
        call material_rw2('jca_rho_z.dat', rho_z, .TRUE.)
        call material_rw2('jca_bulk_modulus.dat', bulk, .TRUE.)
        call material_rw2('jca_attenuation.dat', attenuation, .TRUE.)
        call material_rw3('initialconditionVx.dat', vx, .TRUE.)
        call material_rw3('initialconditionVy.dat', vy, .TRUE.)
        call material_rw3('initialconditionVz.dat', vz, .TRUE.)
        pressure = 0.0_real64
        srcx=0.0_real64; srcy=0.0_real64; srcz=0.0_real64; srcxx=0.0_real64; srcyy=0.0_real64; srczz=0.0_real64
        srcp=0.0_real64
        velocnorm=0.0_real64; pressurenorm=0.0_real64
        isource=source%xind+domain%cpml; jsource=source%yind+domain%cpml; ksource=source%zind+domain%cpml

        if (source%source_type == 'ac') then
            call loadsource('seismicsourcex.dat', source%time_steps, srcx)
            call loadsource('seismicsourcey.dat', source%time_steps, srcy)
            call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        else if (is_pressure_injection_source(source%source_type)) then
            call loadsource('seismicsourcep.dat', source%time_steps, srcp)
        else
            call loadsource('seismicsourcexx.dat', source%time_steps, srcxx)
            call loadsource('seismicsourceyy.dat', source%time_steps, srcyy)
            call loadsource('seismicsourcezz.dat', source%time_steps, srczz)
        endif

        do it=1,source%time_steps
            !$omp parallel do collapse(3) private(div_v) schedule(static)
            do k=2,nz-1
                do j=2,ny-1
                    do i=2,nx-1
                        div_v=(vx(i+1,j,k)-vx(i-1,j,k))/(2.0_real64*dx)+ &
                              (vy(i,j+1,k)-vy(i,j-1,k))/(2.0_real64*dy)+ &
                              (vz(i,j,k+1)-vz(i,j,k-1))/(2.0_real64*dz)
                        pressure(i,j,k)=pressure(i,j,k)-dt*(bulk(i,k)*div_v+attenuation(i,k)*pressure(i,j,k))
                    enddo
                enddo
            enddo
            !$omp end parallel do

            !$omp parallel do collapse(3) private(px,py,pz) schedule(static)
            do k=2,nz-1
                do j=2,ny-1
                    do i=2,nx-1
                        px=(pressure(i+1,j,k)-pressure(i-1,j,k))/(2.0_real64*dx)
                        py=(pressure(i,j+1,k)-pressure(i,j-1,k))/(2.0_real64*dy)
                        pz=(pressure(i,j,k+1)-pressure(i,j,k-1))/(2.0_real64*dz)
                        vx(i,j,k)=vx(i,j,k)-dt*px/positive(rho_x(i,k),1.0e-6_real64)
                        vy(i,j,k)=vy(i,j,k)-dt*py/positive(rho_y(i,k),1.0e-6_real64)
                        vz(i,j,k)=vz(i,j,k)-dt*pz/positive(rho_z(i,k),1.0e-6_real64)
                    enddo
                enddo
            enddo
            !$omp end parallel do

            vx(isource,jsource,ksource)=vx(isource,jsource,ksource)+srcx(it)*dt/positive(rho_x(isource,ksource),1.0e-6_real64)
            vy(isource,jsource,ksource)=vy(isource,jsource,ksource)+srcy(it)*dt/positive(rho_y(isource,ksource),1.0e-6_real64)
            vz(isource,jsource,ksource)=vz(isource,jsource,ksource)+srcz(it)*dt/positive(rho_z(isource,ksource),1.0e-6_real64)
            pressure(isource,jsource,ksource)=pressure(isource,jsource,ksource)+(srcxx(it)+srcyy(it)+srczz(it))*dt
            pressure(isource,jsource,ksource)=pressure(isource,jsource,ksource)+srcp(it)*dt
            vx(1,:,:)=0.0_real64; vy(:,1,:)=0.0_real64; vz(:,:,1)=0.0_real64
            vx(nx,:,:)=0.0_real64; vy(:,ny,:)=0.0_real64; vz(:,:,nz)=0.0_real64

            velocnorm(it)=maxval(sqrt(vx**2+vy**2+vz**2)); pressurenorm(it)=maxval(abs(pressure))
            if (velocnorm(it) > stability_threshold) stop 'JCA 2.5D solution became unstable'
            call write_image3(vx,nx,ny,nz,source,it,'Vx',SINGLE); call write_image3(vy,nx,ny,nz,source,it,'Vy',SINGLE)
            call write_image3(vz,nx,ny,nz,source,it,'Vz',SINGLE); call write_image3(pressure,nx,ny,nz,source,it,'Pp',SINGLE)
        enddo

        call write_array('velocity_norm.dat', source%time_steps, velocnorm)
        call write_array('pore_pressure_norm.dat', source%time_steps, pressurenorm)
    end subroutine jcadg25

    subroutine jcadg3(domain, source, SINGLE_OUTPUT)
        type(Domain_Type), intent(in) :: domain
        type(Source_Type), intent(in) :: source
        logical, intent(in), optional :: SINGLE_OUTPUT

        integer :: nx, ny, nz, i, j, k, it, isource, jsource, ksource
        real(real64) :: dx, dy, dz, dt, div_v, px, py, pz
        logical :: SINGLE
        real(real64), allocatable :: rho_x(:,:,:), rho_y(:,:,:), rho_z(:,:,:), bulk(:,:,:), attenuation(:,:,:)
        real(real64), allocatable :: vx(:,:,:), vy(:,:,:), vz(:,:,:), pressure(:,:,:)
        real(real64), allocatable :: srcx(:), srcy(:), srcz(:), srcxx(:), srcyy(:), srczz(:), srcp(:)
        real(real64), allocatable :: velocnorm(:), pressurenorm(:)

        SINGLE = .TRUE.; if (present(SINGLE_OUTPUT)) SINGLE = SINGLE_OUTPUT
        nx=domain%nx; ny=domain%ny; nz=domain%nz; dx=domain%dx; dy=domain%dy; dz=domain%dz; dt=source%dt
        allocate(rho_x(nx,ny,nz), rho_y(nx,ny,nz), rho_z(nx,ny,nz), bulk(nx,ny,nz), attenuation(nx,ny,nz))
        allocate(vx(nx,ny,nz), vy(nx,ny,nz), vz(nx,ny,nz), pressure(nx,ny,nz))
        allocate(srcx(source%time_steps), srcy(source%time_steps), srcz(source%time_steps))
        allocate(srcxx(source%time_steps), srcyy(source%time_steps), srczz(source%time_steps))
        allocate(srcp(source%time_steps))
        allocate(velocnorm(source%time_steps), pressurenorm(source%time_steps))

        call material_rw3('jca_rho_x.dat', rho_x, .TRUE.)
        call material_rw3('jca_rho_y.dat', rho_y, .TRUE.)
        call material_rw3('jca_rho_z.dat', rho_z, .TRUE.)
        call material_rw3('jca_bulk_modulus.dat', bulk, .TRUE.)
        call material_rw3('jca_attenuation.dat', attenuation, .TRUE.)
        call material_rw3('initialconditionVx.dat', vx, .TRUE.)
        call material_rw3('initialconditionVy.dat', vy, .TRUE.)
        call material_rw3('initialconditionVz.dat', vz, .TRUE.)
        pressure = 0.0_real64
        srcx=0.0_real64; srcy=0.0_real64; srcz=0.0_real64; srcxx=0.0_real64; srcyy=0.0_real64; srczz=0.0_real64
        srcp=0.0_real64
        velocnorm=0.0_real64; pressurenorm=0.0_real64
        isource=source%xind+domain%cpml; jsource=source%yind+domain%cpml; ksource=source%zind+domain%cpml

        if (source%source_type == 'ac') then
            call loadsource('seismicsourcex.dat', source%time_steps, srcx)
            call loadsource('seismicsourcey.dat', source%time_steps, srcy)
            call loadsource('seismicsourcez.dat', source%time_steps, srcz)
        else if (is_pressure_injection_source(source%source_type)) then
            call loadsource('seismicsourcep.dat', source%time_steps, srcp)
        else
            call loadsource('seismicsourcexx.dat', source%time_steps, srcxx)
            call loadsource('seismicsourceyy.dat', source%time_steps, srcyy)
            call loadsource('seismicsourcezz.dat', source%time_steps, srczz)
        endif

        do it=1,source%time_steps
            !$omp parallel do collapse(3) private(div_v) schedule(static)
            do k=2,nz-1
                do j=2,ny-1
                    do i=2,nx-1
                        div_v=(vx(i+1,j,k)-vx(i-1,j,k))/(2.0_real64*dx)+ &
                              (vy(i,j+1,k)-vy(i,j-1,k))/(2.0_real64*dy)+ &
                              (vz(i,j,k+1)-vz(i,j,k-1))/(2.0_real64*dz)
                        pressure(i,j,k)=pressure(i,j,k)-dt*(bulk(i,j,k)*div_v+attenuation(i,j,k)*pressure(i,j,k))
                    enddo
                enddo
            enddo
            !$omp end parallel do

            !$omp parallel do collapse(3) private(px,py,pz) schedule(static)
            do k=2,nz-1
                do j=2,ny-1
                    do i=2,nx-1
                        px=(pressure(i+1,j,k)-pressure(i-1,j,k))/(2.0_real64*dx)
                        py=(pressure(i,j+1,k)-pressure(i,j-1,k))/(2.0_real64*dy)
                        pz=(pressure(i,j,k+1)-pressure(i,j,k-1))/(2.0_real64*dz)
                        vx(i,j,k)=vx(i,j,k)-dt*px/positive(rho_x(i,j,k),1.0e-6_real64)
                        vy(i,j,k)=vy(i,j,k)-dt*py/positive(rho_y(i,j,k),1.0e-6_real64)
                        vz(i,j,k)=vz(i,j,k)-dt*pz/positive(rho_z(i,j,k),1.0e-6_real64)
                    enddo
                enddo
            enddo
            !$omp end parallel do

            vx(isource,jsource,ksource)=vx(isource,jsource,ksource)+srcx(it)*dt/positive(rho_x(isource,jsource,ksource),1.0e-6_real64)
            vy(isource,jsource,ksource)=vy(isource,jsource,ksource)+srcy(it)*dt/positive(rho_y(isource,jsource,ksource),1.0e-6_real64)
            vz(isource,jsource,ksource)=vz(isource,jsource,ksource)+srcz(it)*dt/positive(rho_z(isource,jsource,ksource),1.0e-6_real64)
            pressure(isource,jsource,ksource)=pressure(isource,jsource,ksource)+(srcxx(it)+srcyy(it)+srczz(it))*dt
            pressure(isource,jsource,ksource)=pressure(isource,jsource,ksource)+srcp(it)*dt
            vx(1,:,:)=0.0_real64; vy(:,1,:)=0.0_real64; vz(:,:,1)=0.0_real64
            vx(nx,:,:)=0.0_real64; vy(:,ny,:)=0.0_real64; vz(:,:,nz)=0.0_real64

            velocnorm(it)=maxval(sqrt(vx**2+vy**2+vz**2)); pressurenorm(it)=maxval(abs(pressure))
            if (velocnorm(it) > stability_threshold) stop 'JCA 3D solution became unstable'
            call write_image3(vx,nx,ny,nz,source,it,'Vx',SINGLE); call write_image3(vy,nx,ny,nz,source,it,'Vy',SINGLE)
            call write_image3(vz,nx,ny,nz,source,it,'Vz',SINGLE); call write_image3(pressure,nx,ny,nz,source,it,'Pp',SINGLE)
        enddo

        call write_array('velocity_norm.dat', source%time_steps, velocnorm)
        call write_array('pore_pressure_norm.dat', source%time_steps, pressurenorm)
    end subroutine jcadg3

end module jcadg
