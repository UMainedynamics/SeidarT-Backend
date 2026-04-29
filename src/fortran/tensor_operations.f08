module tensor_operations 

use iso_fortran_env, only: real64 

implicit none 

private 

public :: eigen2, eigen4, & 
            array_eigenvalues22, array_eigenvalues42, &
            array_eigenvalues23, array_eigenvalues43

contains 

    ! --------------------------------------------------------------------------
    subroutine eigen2(a11, a12, a13, a22, a23, a33, lambda)
        !! Eigenvalues of a 2nd order (3x3) real symmetric matrix 
        real(real64), intent(in) :: a11, a12, a13, a22, a23, a33
        real(real64), intent(out) :: lambda(3) 

        real(real64) :: A(3,3)
        real(real64) :: work_query(1) 
        real(real64), allocatable :: work(:) 
        integer :: lwork, info 

        A(:,:) = 0.0_real64
        A(1,1) = a11 
        A(2,2) = a22 
        A(3,3) = a33 
        A(2,1) = a12; A(1,2) = a12 
        A(3,1) = a13; A(1,3) = a13 
        A(3,2) = a23; A(2,3) = a23

        lwork = -1 
        call dsyev('N', 'U', 3, A, 3, lambda, work_query, lwork, info) 
        if (info /= 0) then 
            write(*,*) 'eigen2: DSYEV failed, info = ', info 
            error stop 'DSYEV failure in eigen2'
        endif

        lwork = max(1, int(work_query(1)))
        allocate(work(lwork))

        call dsyev('N', 'U', 3, A, 3, lambda, work, lwork, info)
        deallocate(work)

        if (info /= 0) then
            write(*,*) 'eigen2: DSYEV failed, info = ', info
            error stop 'DSYEV failure in eigen2'
        end if

    end subroutine eigen2 

    ! --------------------------------------------------------------------------
    subroutine eigen4(a11, a12,    a13,    a14,    a15,    a16, & 
                                a22,    a23,    a24,    a25,    a26, &
                                        a33,    a34,    a35,    a36, &
                                                a44,    a45,    a46, &
                                                        a55,    a56, &
                                                                a66, lambda)
        !! Eigenvalues of a voigt matrix
        real(real64), intent(in) :: a11, a12, a13, a14, a15, a16
        real(real64), intent(in) :: a22, a23, a24, a25, a26
        real(real64), intent(in) :: a33, a34, a35, a36
        real(real64), intent(in) :: a44, a45, a46
        real(real64), intent(in) :: a55, a56, a66
        real(real64), intent(out) :: lambda(6) 

        real(real64) :: A(6,6)
        real(real64) :: work_query(1)
        real(real64), allocatable :: work(:) 
        integer :: lwork, info
        real(real64) :: s(6)
        integer :: ii,jj

        A(:,:) = 0.0_real64
        A(1,1) = a11 
        A(2,2) = a22 
        A(3,3) = a33 
        A(4,4) = a44 
        A(5,5) = a55 
        A(6,6) = a66 

        A(2,1) = a12; A(1,2) = a12 
        A(3,1) = a13; A(1,3) = a13 
        A(4,1) = a14; A(1,4) = a14
        A(5,1) = a15; A(1,5) = a15 
        A(6,1) = a16; A(1,6) = a16 

        A(3,2) = a23; A(2,3) = a23 
        A(4,2) = a24; A(2,4) = a24
        A(5,2) = a25; A(2,5) = a25 
        A(6,2) = a26; A(2,6) = a26 

        A(4,3) = a34; A(3,4) = a34
        A(5,3) = a35; A(3,5) = a35 
        A(6,3) = a36; A(3,6) = a36 

        A(5,4) = a45; A(4,5) = a45 
        A(6,4) = a46; A(4,6) = a46 

        A(6,5) = a56; A(5,6) = a56 

        s(1:3) = 1.0_real64
        s(4:6) = sqrt(2.0_real64)


        ! Convert the Voigt 6x6 to the Mandel/Kelvin 6x6 
        do concurrent (ii = 1:6, jj = 1:6)
            A(ii,jj) = s(ii) * s(jj) * A(ii,jj) 
        end do 


        lwork = -1 
        call dsyev('N', 'U', 6, A, 6, lambda, work_query, lwork, info) 
        if (info /= 0) then 
            write(*,*) 'eigen4: DSYEV failed, info = ', info 
            error stop 'DSYEV failure in eigen4'
        endif
        
        lwork = max(1, int(work_query(1)))
        allocate(work(lwork))

        call dsyev('N', 'U', 6, A, 6, lambda, work, lwork, info)
        deallocate(work)

        if (info /= 0) then
            write(*,*) 'eigen4: DSYEV failed, info = ', info
            error stop 'DSYEV failure in eigen4'
        end if
    end subroutine eigen4
    ! --------------------------------------------------------------------------

    subroutine array_eigenvalues2_2(a11, a12, a22, eig_array, nx, ny)
        !! Second order tensor, 2D case
        real(real64), intent(in) :: a11(:,:), a12(:,:), a22(:,:)
        
        integer, intent(in) :: nx, ny 
        real(real64), intent(out) :: eig_array(nx, ny, 2)
        integer :: i,j
        real(real64) :: lambda(3) 

        eig_array(:,:,:) = 0.0_real64 
        
        do j = 1,ny 
            do i=1,nx 
                call eigen2(a11(i,j), a12(i,j), 0.0_real64, &
                            a22(i,j), 0.0_real64, 0.0_real64, lambda )
                eig_array(i,j,:) = lambda(2:3)
            enddo 
        enddo 

    end subroutine array_eigenvalues2_2

    ! --------------------------------------------------------------------------
    subroutine array_eigenvalues4_2(a11, a13, a15, a33, a35, a55, &
                                    eig_array, nx, ny)
        !! 4th order, 2D
        real(real64), intent(in) :: a11(:,:), a13(:,:), a15(:,:), a33(:,:), a35(:,:), a55(:,:)
        integer, intent(in) :: nx, ny

        real(real64), intent(out) :: eig_array(nx, ny, 4)
        
        integer :: i,j
        real(real64) :: lambda(6) 

        eig_array(:,:,:) = 0.0_real64 
        
        do j = 1,ny 
            do i=1,nx 
                call eigen4(a11(i,j), 0.0_real64, a13(i,j), 0.0_real64, a15(i,j), 0.0_real64, &
                            0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
                            a33(i,j), 0.0_real64, a35(i,j), 0.0_real64, &
                            0.0_real64, 0.0_real64, 0.0_real64, &
                            a55(i,j), 0.0_real64, 0.0_real64, lambda )
                eig_array(i,j,:) = lambda(3:6)
            enddo 
        enddo 

    end subroutine array_eigenvalues4_2

    ! --------------------------------------------------------------------------

    subroutine array_eigenvalues2_25(a11, a12, a13, a22, a23, a33, eig_array, nx, ny)
        !! Second order, 3D
        real(real64), intent(in) :: a11(:,:), a12(:,:), a13(:,:)
        real(real64), intent(in) :: a22(:,:), a23(:,:), a33(:,:) 
        integer, intent(in) :: nx, ny 
        real(real64), intent(out) :: eig_array(nx, ny, 3)
        integer :: i,j
        real(real64) :: lambda(3) 

        eig_array(:,:,:) = 0.0_real64 
        
        do j = 1,ny 
            do i=1,nx 
                call eigen2(a11(i,j), a12(i,j), a13(i,j), &
                            a22(i,j), a23(i,j), a33(i,j), lambda )
                eig_array(i,j,:) = lambda(:)
            enddo 
        enddo 

    end subroutine array_eigenvalues2_25

    ! --------------------------------------------------------------------------
    subroutine array_eigenvalues4_25(a11, a12, a13, a14, a15, a16, &
                                    a22, a23, a24, a25, a26, &
                                    a33, a34, a35, a36, a44, a45, a46, &
                                    a55, a56, a66, eig_array, nx, ny)
        !! 4th order 3D
        real(real64), intent(in) :: a11(:,:), a12(:,:), a13(:,:), a14(:,:), a15(:,:), a16(:,:)
        real(real64), intent(in) :: a22(:,:), a23(:,:), a24(:,:), a25(:,:), a26(:,:) 
        real(real64), intent(in) :: a33(:,:), a34(:,:), a35(:,:), a36(:,:), a44(:,:) 
        real(real64), intent(in) :: a45(:,:), a46(:,:), a55(:,:), a56(:,:), a66(:,:)
        integer, intent(in) :: nx, ny

        real(real64), intent(out) :: eig_array(nx, ny, 6)
        
        integer :: i,j
        real(real64) :: lambda(6) 

        eig_array(:,:,:) = 0.0_real64 
        
        do j = 1,ny 
            do i=1,nx 
                call eigen4(a11(i,j), a12(i,j), a13(i,j), a14(i,j), a15(i,j), a16(i,j), &
                            a22(i,j), a23(i,j), a24(i,j), a25(i,j), a26(i,j), &
                            a33(i,j), a34(i,j), a35(i,j), a36(i,j), &
                            a44(i,j), a45(i,j), a46(i,j), &
                            a55(i,j), a56(i,j), a66(i,j), lambda )
                eig_array(i,j,:) = lambda(:)
            enddo 
        enddo 

    end subroutine array_eigenvalues4_25

    ! --------------------------------------------------------------------------

    subroutine array_eigenvalues2_3(a11, a12, a13, a22, a23, a33, &
                                    eig_array, nx, ny, nz)
        !! Second order, 3D
        real(real64), intent(in) :: a11(:,:,:), a12(:,:,:), a13(:,:,:)
        real(real64), intent(in) :: a22(:,:,:), a23(:,:,:), a33(:,:,:) 
        integer, intent(in) :: nx, ny, nz
        real(real64), intent(out) :: eig_array(nx, ny, nz, 3)
        integer :: i,j,k
        real(real64) :: lambda(3) 

        eig_array(:,:,:) = 0.0_real64 
        
        do k = 1,nz 
            do j = 1,ny 
                do i=1,nx 
                    call eigen2(a11(i,j,k), a12(i,j,k), a13(i,j,k), &
                                a22(i,j,k), a23(i,j,k), a33(i,j,k), lambda )
                    eig_array(i,j,k,:) = lambda(:)
                enddo 
            enddo
        enddo  

    end subroutine array_eigenvalues2_3

    ! --------------------------------------------------------------------------
    subroutine array_eigenvalues4_3(a11, a12, a13, a14, a15, a16, &
                                    a22, a23, a24, a25, a26, &
                                    a33, a34, a35, a36, a44, a45, a46, &
                                    a55, a56, a66, eig_array, nx, ny, nz)
        !! 4th order 3D
        real(real64), intent(in) :: a11(:,:,:), a12(:,:,:), a13(:,:,:), a14(:,:,:), a15(:,:,:), a16(:,:,:)
        real(real64), intent(in) :: a22(:,:,:), a23(:,:,:), a24(:,:,:), a25(:,:,:), a26(:,:,:) 
        real(real64), intent(in) :: a33(:,:,:), a34(:,:,:), a35(:,:,:), a36(:,:,:), a44(:,:,:) 
        real(real64), intent(in) :: a45(:,:,:), a46(:,:,:), a55(:,:,:), a56(:,:,:), a66(:,:,:)
        integer, intent(in) :: nx, ny, nz

        real(real64), intent(out) :: eig_array(nx, ny, nz, 6)
        
        integer :: i,j
        real(real64) :: lambda(6) 

        eig_array(:,:,:,:) = 0.0_real64 
        
        do k = 1,nz
            do j = 1,ny 
                do i=1,nx 
                    call eigen4(a11(i,j,k), a12(i,j,k), a13(i,j,k), a14(i,j,k), a15(i,j,k), a16(i,j,k), &
                                a22(i,j,k), a23(i,j,k), a24(i,j,k), a25(i,j,k), a26(i,j,k), &
                                a33(i,j,k), a34(i,j,k), a35(i,j,k), a36(i,j,k), &
                                a44(i,j,k), a45(i,j,k), a46(i,j,k), &
                                a55(i,j,k), a56(i,j,k), a66(i,j,k), lambda )
                    eig_array(i,j,k,:) = lambda(:)
                enddo 
            enddo 
        enddo

    end subroutine array_eigenvalues4_25
end module tensor_operations