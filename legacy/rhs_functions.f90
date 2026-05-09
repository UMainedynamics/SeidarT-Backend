
    subroutine rhs_biot25_cell(vx,vy,vz,qx,qy,qz,p,sxx,syy,szz,sxy,sxz,syz, &
        c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
        rho,rho_f,eta,ax,ay,az,m,kx,ky,kz,gx,gy,gz,dx,dy,dz,speed, &
        rvx,rvy,rvz,rqx,rqy,rqz,rp,rsxx,rsyy,rszz,rsxy,rsxz,rsyz)
        implicit none

        real(real64), intent(in) :: vx(:,:,:), vy(:,:,:), vz(:,:,:), qx(:,:,:), qy(:,:,:), qz(:,:,:)
real(real64), intent(in) :: p(:,:,:), sxx(:,:,:), syy(:,:,:), szz(:,:,:), sxy(:,:,:), sxz(:,:,:), syz(:,:,:)
real(real64), intent(in) :: c11(:,:), c12(:,:), c13(:,:), c14(:,:), c15(:,:), c16(:,:)
real(real64), intent(in) :: c22(:,:), c23(:,:), c24(:,:), c25(:,:), c26(:,:)
real(real64), intent(in) :: c33(:,:), c34(:,:), c35(:,:), c36(:,:)
real(real64), intent(in) :: c44(:,:), c45(:,:), c46(:,:), c55(:,:), c56(:,:), c66(:,:)
real(real64), intent(in) :: rho(:,:), rho_f(:,:), eta(:,:), ax(:,:), ay(:,:), az(:,:)
real(real64), intent(in) :: m(:,:), kx(:,:), ky(:,:), kz(:,:), gx(:,:), gy(:,:), gz(:,:)
real(real64), intent(out) :: rvx(:,:,:), rvy(:,:,:), rvz(:,:,:), rqx(:,:,:), rqy(:,:,:), rqz(:,:,:)
real(real64), intent(out) :: rp(:,:,:), rsxx(:,:,:), rsyy(:,:,:), rszz(:,:,:), rsxy(:,:,:), rsxz(:,:,:), rsyz(:,:,:)

        real(real64), intent(in) :: dx, dy, dz, speed

        integer :: i, j, k, nx, ny, nz
        real(real64) :: exx, eyy, ezz, eyz, exz, exy
        real(real64) :: div_q, px, py, pz, dampx, dampy, dampz

        nx=size(vx,1); ny=size(vx,2); nz=size(vx,3)
        rvx=0.0_real64; rvy=0.0_real64; rvz=0.0_real64; rqx=0.0_real64; rqy=0.0_real64; rqz=0.0_real64
        rp=0.0_real64; rsxx=0.0_real64; rsyy=0.0_real64; rszz=0.0_real64; rsxy=0.0_real64; rsxz=0.0_real64; rsyz=0.0_real64

        !$omp parallel do collapse(3) private(exx,eyy,ezz,eyz,exz,exy,div_q,px,py,pz,dampx,dampy,dampz) schedule(static)
        do k=2,nz-1
            do j=2,ny-1
                do i=2,nx-1
                    exx = (vx(i+1,j,k)-vx(i-1,j,k))/(2.0_real64*dx)
                    eyy = (vy(i,j+1,k)-vy(i,j-1,k))/(2.0_real64*dy)
                    ezz = (vz(i,j,k+1)-vz(i,j,k-1))/(2.0_real64*dz)
                    eyz = (vy(i,j,k+1)-vy(i,j,k-1))/(2.0_real64*dz) + (vz(i,j+1,k)-vz(i,j-1,k))/(2.0_real64*dy)
                    exz = (vx(i,j,k+1)-vx(i,j,k-1))/(2.0_real64*dz) + (vz(i+1,j,k)-vz(i-1,j,k))/(2.0_real64*dx)
                    exy = (vx(i,j+1,k)-vx(i,j-1,k))/(2.0_real64*dy) + (vy(i+1,j,k)-vy(i-1,j,k))/(2.0_real64*dx)
                    
                    div_q = (qx(i+1,j,k)-qx(i-1,j,k))/(2.0_real64*dx) + &
                            (qy(i,j+1,k)-qy(i,j-1,k))/(2.0_real64*dy) + &
                            (qz(i,j,k+1)-qz(i,j,k-1))/(2.0_real64*dz)
                            
                    rp(i,j,k) = -m(i,k) * (((ax(i,k)+ay(i,k)+az(i,k))/3.0_real64)*(exx+eyy+ezz) + div_q)

                    rsxx(i,j,k) = c11(i,k)*exx + c12(i,k)*eyy + c13(i,k)*ezz + c14(i,k)*eyz + c15(i,k)*exz + c16(i,k)*exy - ax(i,k)*rp(i,j,k) - gx(i,k)*sxx(i,j,k)
                    rsyy(i,j,k) = c12(i,k)*exx + c22(i,k)*eyy + c23(i,k)*ezz + c24(i,k)*eyz + c25(i,k)*exz + c26(i,k)*exy - ay(i,k)*rp(i,j,k) - gy(i,k)*syy(i,j,k)
                    rszz(i,j,k) = c13(i,k)*exx + c23(i,k)*eyy + c33(i,k)*ezz + c34(i,k)*eyz + c35(i,k)*exz + c36(i,k)*exy - az(i,k)*rp(i,j,k) - gz(i,k)*szz(i,j,k)
                    rsyz(i,j,k) = c14(i,k)*exx + c24(i,k)*eyy + c34(i,k)*ezz + c44(i,k)*eyz + c45(i,k)*exz + c46(i,k)*exy - gy(i,k)*syz(i,j,k)
                    rsxz(i,j,k) = c15(i,k)*exx + c25(i,k)*eyy + c35(i,k)*ezz + c45(i,k)*eyz + c55(i,k)*exz + c56(i,k)*exy - gx(i,k)*sxz(i,j,k)
                    rsxy(i,j,k) = c16(i,k)*exx + c26(i,k)*eyy + c36(i,k)*ezz + c46(i,k)*eyz + c56(i,k)*exz + c66(i,k)*exy - gx(i,k)*sxy(i,j,k)

                    px=(p(i+1,j,k)-p(i-1,j,k))/(2.0_real64*dx)
                    py=(p(i,j+1,k)-p(i,j-1,k))/(2.0_real64*dy)
                    pz=(p(i,j,k+1)-p(i,j,k-1))/(2.0_real64*dz)
                    
                    rvx(i,j,k)=((sxx(i+1,j,k)-sxx(i-1,j,k))/(2.0_real64*dx) + &
                                (sxy(i,j+1,k)-sxy(i,j-1,k))/(2.0_real64*dy) + &
                                (sxz(i,j,k+1)-sxz(i,j,k-1))/(2.0_real64*dz))/positive(rho(i,k),1.0_real64)
                    rvy(i,j,k)=((sxy(i+1,j,k)-sxy(i-1,j,k))/(2.0_real64*dx) + &
                                (syy(i,j+1,k)-syy(i,j-1,k))/(2.0_real64*dy) + &
                                (syz(i,j,k+1)-syz(i,j,k-1))/(2.0_real64*dz))/positive(rho(i,k),1.0_real64)
                    rvz(i,j,k)=((sxz(i+1,j,k)-sxz(i-1,j,k))/(2.0_real64*dx) + &
                                (syz(i,j+1,k)-syz(i,j-1,k))/(2.0_real64*dy) + &
                                (szz(i,j,k+1)-szz(i,j,k-1))/(2.0_real64*dz))/positive(rho(i,k),1.0_real64)
                                
                    rqx(i,j,k)=-kx(i,k)*px/(positive(eta(i,k),1.0e-6_real64)*positive(rho_f(i,k),1.0_real64))
                    rqy(i,j,k)=-ky(i,k)*py/(positive(eta(i,k),1.0e-6_real64)*positive(rho_f(i,k),1.0_real64))
                    rqz(i,j,k)=-kz(i,k)*pz/(positive(eta(i,k),1.0e-6_real64)*positive(rho_f(i,k),1.0_real64))

                    dampx = speed*((vx(i+1,j,k)-2.0_real64*vx(i,j,k)+vx(i-1,j,k))/(dx*dx)) * dx * 0.1_real64
                    dampy = speed*((vy(i,j+1,k)-2.0_real64*vy(i,j,k)+vy(i,j-1,k))/(dy*dy)) * dy * 0.1_real64
                    dampz = speed*((vz(i,j,k+1)-2.0_real64*vz(i,j,k)+vz(i,j,k-1))/(dz*dz)) * dz * 0.1_real64
                    rvx(i,j,k)=rvx(i,j,k)+dampx; rvy(i,j,k)=rvy(i,j,k)+dampy; rvz(i,j,k)=rvz(i,j,k)+dampz
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine rhs_biot25_cell
    

    subroutine rhs_biot3_cell(vx,vy,vz,qx,qy,qz,p,sxx,syy,szz,sxy,sxz,syz, &
        c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
        rho,rho_f,eta,ax,ay,az,m,kx,ky,kz,gx,gy,gz,dx,dy,dz,speed, &
        rvx,rvy,rvz,rqx,rqy,rqz,rp,rsxx,rsyy,rszz,rsxy,rsxz,rsyz)
        implicit none

        real(real64), intent(in) :: vx(:,:,:), vy(:,:,:), vz(:,:,:), qx(:,:,:), qy(:,:,:), qz(:,:,:)
real(real64), intent(in) :: p(:,:,:), sxx(:,:,:), syy(:,:,:), szz(:,:,:), sxy(:,:,:), sxz(:,:,:), syz(:,:,:)
real(real64), intent(in) :: c11(:,:,:), c12(:,:,:), c13(:,:,:), c14(:,:,:), c15(:,:,:), c16(:,:,:)
real(real64), intent(in) :: c22(:,:,:), c23(:,:,:), c24(:,:,:), c25(:,:,:), c26(:,:,:)
real(real64), intent(in) :: c33(:,:,:), c34(:,:,:), c35(:,:,:), c36(:,:,:)
real(real64), intent(in) :: c44(:,:,:), c45(:,:,:), c46(:,:,:), c55(:,:,:), c56(:,:,:), c66(:,:,:)
real(real64), intent(in) :: rho(:,:,:), rho_f(:,:,:), eta(:,:,:), ax(:,:,:), ay(:,:,:), az(:,:,:)
real(real64), intent(in) :: m(:,:,:), kx(:,:,:), ky(:,:,:), kz(:,:,:), gx(:,:,:), gy(:,:,:), gz(:,:,:)
real(real64), intent(out) :: rvx(:,:,:), rvy(:,:,:), rvz(:,:,:), rqx(:,:,:), rqy(:,:,:), rqz(:,:,:)
real(real64), intent(out) :: rp(:,:,:), rsxx(:,:,:), rsyy(:,:,:), rszz(:,:,:), rsxy(:,:,:), rsxz(:,:,:), rsyz(:,:,:)

        real(real64), intent(in) :: dx, dy, dz, speed

        integer :: i, j, k, nx, ny, nz
        real(real64) :: exx, eyy, ezz, eyz, exz, exy
        real(real64) :: div_q, px, py, pz, dampx, dampy, dampz

        nx=size(vx,1); ny=size(vx,2); nz=size(vx,3)
        rvx=0.0_real64; rvy=0.0_real64; rvz=0.0_real64; rqx=0.0_real64; rqy=0.0_real64; rqz=0.0_real64
        rp=0.0_real64; rsxx=0.0_real64; rsyy=0.0_real64; rszz=0.0_real64; rsxy=0.0_real64; rsxz=0.0_real64; rsyz=0.0_real64

        !$omp parallel do collapse(3) private(exx,eyy,ezz,eyz,exz,exy,div_q,px,py,pz,dampx,dampy,dampz) schedule(static)
        do k=2,nz-1
            do j=2,ny-1
                do i=2,nx-1
                    exx = (vx(i+1,j,k)-vx(i-1,j,k))/(2.0_real64*dx)
                    eyy = (vy(i,j+1,k)-vy(i,j-1,k))/(2.0_real64*dy)
                    ezz = (vz(i,j,k+1)-vz(i,j,k-1))/(2.0_real64*dz)
                    eyz = (vy(i,j,k+1)-vy(i,j,k-1))/(2.0_real64*dz) + (vz(i,j+1,k)-vz(i,j-1,k))/(2.0_real64*dy)
                    exz = (vx(i,j,k+1)-vx(i,j,k-1))/(2.0_real64*dz) + (vz(i+1,j,k)-vz(i-1,j,k))/(2.0_real64*dx)
                    exy = (vx(i,j+1,k)-vx(i,j-1,k))/(2.0_real64*dy) + (vy(i+1,j,k)-vy(i-1,j,k))/(2.0_real64*dx)
                    
                    div_q = (qx(i+1,j,k)-qx(i-1,j,k))/(2.0_real64*dx) + &
                            (qy(i,j+1,k)-qy(i,j-1,k))/(2.0_real64*dy) + &
                            (qz(i,j,k+1)-qz(i,j,k-1))/(2.0_real64*dz)
                            
                    rp(i,j,k) = -m(i,j,k) * (((ax(i,j,k)+ay(i,j,k)+az(i,j,k))/3.0_real64)*(exx+eyy+ezz) + div_q)

                    rsxx(i,j,k) = c11(i,j,k)*exx + c12(i,j,k)*eyy + c13(i,j,k)*ezz + c14(i,j,k)*eyz + c15(i,j,k)*exz + c16(i,j,k)*exy - ax(i,j,k)*rp(i,j,k) - gx(i,j,k)*sxx(i,j,k)
                    rsyy(i,j,k) = c12(i,j,k)*exx + c22(i,j,k)*eyy + c23(i,j,k)*ezz + c24(i,j,k)*eyz + c25(i,j,k)*exz + c26(i,j,k)*exy - ay(i,j,k)*rp(i,j,k) - gy(i,j,k)*syy(i,j,k)
                    rszz(i,j,k) = c13(i,j,k)*exx + c23(i,j,k)*eyy + c33(i,j,k)*ezz + c34(i,j,k)*eyz + c35(i,j,k)*exz + c36(i,j,k)*exy - az(i,j,k)*rp(i,j,k) - gz(i,j,k)*szz(i,j,k)
                    rsyz(i,j,k) = c14(i,j,k)*exx + c24(i,j,k)*eyy + c34(i,j,k)*ezz + c44(i,j,k)*eyz + c45(i,j,k)*exz + c46(i,j,k)*exy - gy(i,j,k)*syz(i,j,k)
                    rsxz(i,j,k) = c15(i,j,k)*exx + c25(i,j,k)*eyy + c35(i,j,k)*ezz + c45(i,j,k)*eyz + c55(i,j,k)*exz + c56(i,j,k)*exy - gx(i,j,k)*sxz(i,j,k)
                    rsxy(i,j,k) = c16(i,j,k)*exx + c26(i,j,k)*eyy + c36(i,j,k)*ezz + c46(i,j,k)*eyz + c56(i,j,k)*exz + c66(i,j,k)*exy - gx(i,j,k)*sxy(i,j,k)

                    px=(p(i+1,j,k)-p(i-1,j,k))/(2.0_real64*dx)
                    py=(p(i,j+1,k)-p(i,j-1,k))/(2.0_real64*dy)
                    pz=(p(i,j,k+1)-p(i,j,k-1))/(2.0_real64*dz)
                    
                    rvx(i,j,k)=((sxx(i+1,j,k)-sxx(i-1,j,k))/(2.0_real64*dx) + &
                                (sxy(i,j+1,k)-sxy(i,j-1,k))/(2.0_real64*dy) + &
                                (sxz(i,j,k+1)-sxz(i,j,k-1))/(2.0_real64*dz))/positive(rho(i,j,k),1.0_real64)
                    rvy(i,j,k)=((sxy(i+1,j,k)-sxy(i-1,j,k))/(2.0_real64*dx) + &
                                (syy(i,j+1,k)-syy(i,j-1,k))/(2.0_real64*dy) + &
                                (syz(i,j,k+1)-syz(i,j,k-1))/(2.0_real64*dz))/positive(rho(i,j,k),1.0_real64)
                    rvz(i,j,k)=((sxz(i+1,j,k)-sxz(i-1,j,k))/(2.0_real64*dx) + &
                                (syz(i,j+1,k)-syz(i,j-1,k))/(2.0_real64*dy) + &
                                (szz(i,j,k+1)-szz(i,j,k-1))/(2.0_real64*dz))/positive(rho(i,j,k),1.0_real64)
                                
                    rqx(i,j,k)=-kx(i,j,k)*px/(positive(eta(i,j,k),1.0e-6_real64)*positive(rho_f(i,j,k),1.0_real64))
                    rqy(i,j,k)=-ky(i,j,k)*py/(positive(eta(i,j,k),1.0e-6_real64)*positive(rho_f(i,j,k),1.0_real64))
                    rqz(i,j,k)=-kz(i,j,k)*pz/(positive(eta(i,j,k),1.0e-6_real64)*positive(rho_f(i,j,k),1.0_real64))

                    dampx = speed*((vx(i+1,j,k)-2.0_real64*vx(i,j,k)+vx(i-1,j,k))/(dx*dx)) * dx * 0.1_real64
                    dampy = speed*((vy(i,j+1,k)-2.0_real64*vy(i,j,k)+vy(i,j-1,k))/(dy*dy)) * dy * 0.1_real64
                    dampz = speed*((vz(i,j,k+1)-2.0_real64*vz(i,j,k)+vz(i,j,k-1))/(dz*dz)) * dz * 0.1_real64
                    rvx(i,j,k)=rvx(i,j,k)+dampx; rvy(i,j,k)=rvy(i,j,k)+dampy; rvz(i,j,k)=rvz(i,j,k)+dampz
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine rhs_biot3_cell
    
