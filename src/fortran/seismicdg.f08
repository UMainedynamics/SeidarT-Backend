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

        call fill_modal_from_cell(cell_vx, vx); call fill_modal_from_cell(cell_vz, vz)
        sxx = 0.0_real64; szz = 0.0_real64; sxz = 0.0_real64
        srcx = 0.0_real64; srcz = 0.0_real64; srcxx = 0.0_real64
        srcxz = 0.0_real64; srczz = 0.0_real64; velocnorm = 0.0_real64
        isource = source%xind 
        ksource = source%zind
        max_speed = sqrt(maxval(max(c11, c33)) / positive(minval(rho), 1.0_real64))

        if ( is_plane_wave_source(source%source_type) ) then 
            XMIN = domain%dx
            XMAX = domain%dx * domain%nx
            XMID = 0.50_real64 * (XMIN + XMAX)
            ZMIN = domain%dz 
            ZMAX = domain%dz * domain%nz
            ZMID = 0.50_real64 * (ZMIN + ZMAX)

            call direction_polarization_vector( source%x_z_rotation, &
                                                source%x_y_rotation, &
                                                source%y_z_rotation, &
                                                p, ehat)
            call select_injection_faces(p, active)

            r0 = (/ XMID, 0.0_real64, ZMID /)
            XLOC = XMID
            ZLOC = ZMID
            if ( active(1) ) then
                r0(1) = XMIN
                XLOC = XMIN
            endif
            if ( active(2) ) then
                r0(1) = XMAX
                XLOC = XMAX
            endif
            if ( active(5) ) then
                r0(3) = ZMIN
                ZLOC = ZMIN
            endif
            if ( active(6) ) then
                r0(3) = ZMAX
                ZLOC = ZMAX
            endif

            i_min =  1
            i_max = domain%nx 
            j_min = 1
            j_max = domain%nz 

            call array_eigenvalues4_2(c11, c13, c15, c33, c35, c55, &
                                      eig_array, nx, nz)
            cbackground = sqrt(maxval(eig_array) / minval(rho))
        else if ( source%source_type == 'ac') then 
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

            ! Add source terms
            if ( is_plane_wave_source(source%source_type) ) then
                t = it * source%dt
                if ( active(1) .or. active(2) ) then
                    xi = nint(XLOC / domain%dx)
                    do ip = j_min,j_max
                        r = (/ XLOC, 0.0_real64, ip * domain%dz /)
                        call boundary_seismic_field(r, r0, t, &
                                        source%source_frequency, &
                                        source%x_z_rotation, &
                                        source%x_y_rotation, &
                                        source%y_z_rotation, &
                                        cbackground, source%amplitude, &
                                        source%source_wavelet, &
                                        Vvec, strain_vec, source%source_type)
                        vx(xi,ip) = vx(xi,ip) + Vvec(1)
                        vz(xi,ip) = vz(xi,ip) + Vvec(3)
                        sigmaxx(xi,ip) = sigmaxx(xi,ip) + &
                            c11(xi,ip) * strain_vec(1) + c13(xi,ip) * strain_vec(3) + &
                            c15(xi,ip) * strain_vec(5)
                        sigmazz(xi,ip) = sigmazz(xi,ip) + &
                            c13(xi,ip) * strain_vec(1) + c33(xi,ip) * strain_vec(3) + &
                            c35(xi,ip) * strain_vec(5)
                        sigmaxz(xi,ip) = sigmaxz(xi,ip) + &
                            c15(xi,ip) * strain_vec(1) + c35(xi,ip) * strain_vec(3) + &
                            c55(xi,ip) * strain_vec(5)
                    enddo
                endif
                if ( active(5) .or. active(6) ) then
                    zi = nint(ZLOC / domain%dz)
                    do ip = i_min,i_max
                        r = (/ ip * domain%dx, 0.0_real64, ZLOC /)
                        call boundary_seismic_field(r, r0, t, &
                                        source%source_frequency, &
                                        source%x_z_rotation, &
                                        source%x_y_rotation, &
                                        source%y_z_rotation, &
                                        cbackground, source%amplitude, &
                                        source%source_wavelet, &
                                        Vvec, strain_vec, source%source_type)
                        vx(ip,zi) = vx(ip,zi) + Vvec(1)
                        vz(ip,zi) = vz(ip,zi) + Vvec(3)
                        sigmaxx(ip,zi) = sigmaxx(ip,zi) + &
                            c11(ip,zi) * strain_vec(1) + c13(ip,zi) * strain_vec(3) + &
                            c15(ip,zi) * strain_vec(5)
                        sigmazz(ip,zi) = sigmazz(ip,zi) + &
                            c13(ip,zi) * strain_vec(1) + c33(ip,zi) * strain_vec(3) + &
                            c35(ip,zi) * strain_vec(5)
                        sigmaxz(ip,zi) = sigmaxz(ip,zi) + &
                            c15(ip,zi) * strain_vec(1) + c35(ip,zi) * strain_vec(3) + &
                            c55(ip,zi) * strain_vec(5)
                    enddo
                endif
            else
                ! Add the source term. If it is an accelerated weight drop the src_ij terms are zero
                ! If it is any other source, the src_i terms are zero
                sigmaxx(isource,jsource) = sigmaxx(isource, jsource) + srcxx(it) / rho(isource,jsource)
                sigmaxz(isource+1,jsource+1) = sigmaxz(isource+1, jsource+1) + srcxz(it) / rho(isource+1,jsource+1)  
                sigmazz(isource,jsource) = sigmazz(isource, jsource) + srczz(it) / rho(isource,jsource) 
                vx(isource,jsource) = vx(isource,jsource) + srcx(it) / rho(isource,jsource)
                vz(isource,jsource) = vz(isource,jsource) + srcz(it) / rho(isource,jsource)
            endif


            call cell_average(vx, weights, cell_out); call write_image2(cell_out, nx, nz, source, it, 'Vx', SINGLE)
            velocnorm(it) = maxval(abs(cell_out))
            call cell_average(vz, weights, cell_out); call write_image2(cell_out, nx, nz, source, it, 'Vz', SINGLE)
            velocnorm(it) = max(velocnorm(it), maxval(abs(cell_out)))
            if (velocnorm(it) > stability_threshold) stop 'Seismic DG 2D solution became unstable'
        enddo

        call write_array('velocity_norm.dat', source%time_steps, velocnorm)
    end subroutine seismicdg2

end module seismicdg
