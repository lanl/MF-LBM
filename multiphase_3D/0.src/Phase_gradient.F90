#include "./preprocessor.h"
!=====================================================================================================================================
!---------------------- calculate color gradient ----------------------
!=====================================================================================================================================
subroutine color_gradient
    use Fluid_multiphase, only: cn_x, cn_y, cn_z, phi, ISO4, c_norm, curv
    use Misc_module
    IMPLICIT NONE
    integer :: i, j, k, num, n, ie, iex, iey, iez, overlap_temp
    real(kind=8) :: kxx, kxy, kxz, kyx, kyy, kyz, kzx, kzy, kzz, tmp

    ! ~~~~~~~~~~~~~~~~~~~~~~~ extrapolate phi values to solid boundary nodes ~~~~~~~~~~~~~~~~~~
    !$omp parallel do private(i,j,k,n,iex,iey,iez,ie)
    !$acc kernels present(phi,solid_boundary_nodes)
    !$acc loop independent device_type(NVIDIA)
    do num = 1, num_solid_boundary                        !cover MPI domain with three layer of ghost nodes, to be used in color gradient calculation
        i = solid_boundary_nodes(num)%ix                  !also cover additional layer for the mirror BC at y=1
        j = solid_boundary_nodes(num)%iy
        k = solid_boundary_nodes(num)%iz
        phi(i, j, k) = 0d0
        do n = 1, solid_boundary_nodes(num)%i_fluid_num
            ie = solid_boundary_nodes(num)%neighbor_list(n)
            iex = i + ex(ie)
            iey = j + ey(ie)
            iez = k + ez(ie)
            phi(i, j, k) = phi(i, j, k) + phi(iex, iey, iez)*w_equ(ie)
        end do
        phi(i, j, k) = phi(i, j, k)/solid_boundary_nodes(num)%la_weight
    end do
    !$acc end kernels

    ! ~~~~~~~~~~~~~~~~~~ calculate normal directions of interfaces from phi gradient ~~~~~~~~~~~~~~~~~~
    overlap_temp = 2
    !$omp parallel do private(i,j) collapse(2)
    !$acc kernels present(phi,cn_x,cn_y,cn_z,c_norm)
    do k = 1 - overlap_temp, nz + overlap_temp
        do j = 1 - overlap_temp, ny + overlap_temp
            do i = 1 - overlap_temp, nx + overlap_temp

                cn_x(i, j, k) = &
                    ISO4(1)*(phi(i + 1, j, k) - phi(i - 1, j, k)) + &
                    ISO4(2)*( &
                    phi(i + 1, j + 1, k) - phi(i - 1, j - 1, k) + &
                    phi(i + 1, j - 1, k) - phi(i - 1, j + 1, k) + &
                    phi(i + 1, j, k + 1) - phi(i - 1, j, k - 1) + &
                    phi(i + 1, j, k - 1) - phi(i - 1, j, k + 1))

                cn_y(i, j, k) = &
                    ISO4(1)*(phi(i, j + 1, k) - phi(i, j - 1, k)) + &
                    ISO4(2)*( &
                    phi(i + 1, j + 1, k) - phi(i - 1, j - 1, k) + &
                    phi(i - 1, j + 1, k) - phi(i + 1, j - 1, k) + &
                    phi(i, j + 1, k + 1) - phi(i, j - 1, k - 1) + &
                    phi(i, j + 1, k - 1) - phi(i, j - 1, k + 1))

                cn_z(i, j, k) = &
                    ISO4(1)*(phi(i, j, k + 1) - phi(i, j, k - 1)) + &
                    ISO4(2)*( &
                    phi(i + 1, j, k + 1) - phi(i - 1, j, k - 1) + &
                    phi(i - 1, j, k + 1) - phi(i + 1, j, k - 1) + &
                    phi(i, j + 1, k + 1) - phi(i, j - 1, k - 1) + &
                    phi(i, j - 1, k + 1) - phi(i, j + 1, k - 1))

                c_norm(i, j, k) = dsqrt(cn_x(i, j, k)*cn_x(i, j, k) + cn_y(i, j, k)*cn_y(i, j, k) + cn_z(i, j, k)*cn_z(i, j, k))
                if (c_norm(i, j, k) < 1d-6 .or. walls(i, j, k) == 1) then
                    cn_x(i, j, k) = 0d0
                    cn_y(i, j, k) = 0d0
                    cn_z(i, j, k) = 0d0
                    c_norm(i, j, k) = 0d0
                else
                    cn_x(i, j, k) = cn_x(i, j, k)/c_norm(i, j, k)
                    cn_y(i, j, k) = cn_y(i, j, k)/c_norm(i, j, k)
                    cn_z(i, j, k) = cn_z(i, j, k)/c_norm(i, j, k)    !normalized color gradient - interface normal direction
                end if

            end do
        end do
    end do
    !$acc end kernels

    ! ~~~~~~~~~~~~~~~ numerically alter the normal directions of interfaces to desired contact angle ~~~~~~~~~~~~~~~~
    call alter_color_gradient_solid_surface

    ! ~~~~~~~~~~~~~~ extrapolate normal direction info to solid boundary nodes, to minimize unbalanced forces ~~~~~~~~~~~~~~
    !$omp parallel do private(i,j,k,n,iex,iey,iez,ie)
    !$acc kernels present(cn_x,cn_y,cn_z,solid_boundary_nodes)
    !$acc loop independent device_type(NVIDIA)
    do num = 1, num_solid_boundary              !cover MPI domain with one layer of ghost nodes, to be used in curvature calculation
        i = solid_boundary_nodes(num)%ix
        j = solid_boundary_nodes(num)%iy
        k = solid_boundary_nodes(num)%iz
        if (i >= 0 .and. i <= nx + 1 .and. j >= 0 .and. j <= ny + 1 .and. k >= 0 .and. k <= nz + 1) then
            cn_x(i, j, k) = 0d0
            cn_y(i, j, k) = 0d0
            cn_z(i, j, k) = 0d0
            do n = 1, solid_boundary_nodes(num)%i_fluid_num
                ie = solid_boundary_nodes(num)%neighbor_list(n)
                iex = i + ex(ie)
                iey = j + ey(ie)
                iez = k + ez(ie)
                cn_x(i, j, k) = cn_x(i, j, k) + cn_x(iex, iey, iez)*w_equ(ie)
                cn_y(i, j, k) = cn_y(i, j, k) + cn_y(iex, iey, iez)*w_equ(ie)
                cn_z(i, j, k) = cn_z(i, j, k) + cn_z(iex, iey, iez)*w_equ(ie)
            end do
            cn_x(i, j, k) = cn_x(i, j, k)/solid_boundary_nodes(num)%la_weight
            cn_y(i, j, k) = cn_y(i, j, k)/solid_boundary_nodes(num)%la_weight
            cn_z(i, j, k) = cn_z(i, j, k)/solid_boundary_nodes(num)%la_weight
        end if
    end do
    !$acc end kernels

    ! ~~~~~~~~~~~~~~~~~~ calculate CSF forces based on interace curvature  ~~~~~~~~~~~~~~~~~~
    !$omp parallel do private(i,j,kxx,kxy,kxz,kyx,kyy,kyz,kzx,kzy,kzz) collapse(2)
    !$acc kernels present(curv,cn_x,cn_y,cn_z)
    !$acc loop device_type(NVIDIA)
    do k = 1, nz
        ! !$acc loop device_type(NVIDIA)
        do j = 1, ny
            !    !$acc loop device_type(NVIDIA) vector
            do i = 1, nx
                !  if(walls(i,j,k)==0)then
                kxx = &
                    ISO4(1)*(cn_x(i + 1, j, k) - cn_x(i - 1, j, k)) + &
                    ISO4(2)*( &
                    cn_x(i + 1, j + 1, k) - cn_x(i - 1, j - 1, k) + &
                    cn_x(i + 1, j - 1, k) - cn_x(i - 1, j + 1, k) + &
                    cn_x(i + 1, j, k + 1) - cn_x(i - 1, j, k - 1) + &
                    cn_x(i + 1, j, k - 1) - cn_x(i - 1, j, k + 1))

                kyy = &
                    ISO4(1)*(cn_y(i, j + 1, k) - cn_y(i, j - 1, k)) + &
                    ISO4(2)*( &
                    cn_y(i + 1, j + 1, k) - cn_y(i - 1, j - 1, k) + &
                    cn_y(i - 1, j + 1, k) - cn_y(i + 1, j - 1, k) + &
                    cn_y(i, j + 1, k + 1) - cn_y(i, j - 1, k - 1) + &
                    cn_y(i, j + 1, k - 1) - cn_y(i, j - 1, k + 1))

                kzz = &
                    ISO4(1)*(cn_z(i, j, k + 1) - cn_z(i, j, k - 1)) + &
                    ISO4(2)*( &
                    cn_z(i + 1, j, k + 1) - cn_z(i - 1, j, k - 1) + &
                    cn_z(i - 1, j, k + 1) - cn_z(i + 1, j, k - 1) + &
                    cn_z(i, j + 1, k + 1) - cn_z(i, j - 1, k - 1) + &
                    cn_z(i, j - 1, k + 1) - cn_z(i, j + 1, k - 1))

                kxy = &
                    ISO4(1)*(cn_x(i, j + 1, k) - cn_x(i, j - 1, k)) + &
                    ISO4(2)*( &
                    cn_x(i + 1, j + 1, k) - cn_x(i - 1, j - 1, k) + &
                    cn_x(i - 1, j + 1, k) - cn_x(i + 1, j - 1, k) + &
                    cn_x(i, j + 1, k + 1) - cn_x(i, j - 1, k - 1) + &
                    cn_x(i, j + 1, k - 1) - cn_x(i, j - 1, k + 1))

                kxz = &
                    ISO4(1)*(cn_x(i, j, k + 1) - cn_x(i, j, k - 1)) + &
                    ISO4(2)*( &
                    cn_x(i + 1, j, k + 1) - cn_x(i - 1, j, k - 1) + &
                    cn_x(i - 1, j, k + 1) - cn_x(i + 1, j, k - 1) + &
                    cn_x(i, j + 1, k + 1) - cn_x(i, j - 1, k - 1) + &
                    cn_x(i, j - 1, k + 1) - cn_x(i, j + 1, k - 1))

                kyx = &
                    ISO4(1)*(cn_y(i + 1, j, k) - cn_y(i - 1, j, k)) + &
                    ISO4(2)*( &
                    cn_y(i + 1, j + 1, k) - cn_y(i - 1, j - 1, k) + &
                    cn_y(i + 1, j - 1, k) - cn_y(i - 1, j + 1, k) + &
                    cn_y(i + 1, j, k + 1) - cn_y(i - 1, j, k - 1) + &
                    cn_y(i + 1, j, k - 1) - cn_y(i - 1, j, k + 1))

                kyz = &
                    ISO4(1)*(cn_y(i, j, k + 1) - cn_y(i, j, k - 1)) + &
                    ISO4(2)*( &
                    cn_y(i + 1, j, k + 1) - cn_y(i - 1, j, k - 1) + &
                    cn_y(i - 1, j, k + 1) - cn_y(i + 1, j, k - 1) + &
                    cn_y(i, j + 1, k + 1) - cn_y(i, j - 1, k - 1) + &
                    cn_y(i, j - 1, k + 1) - cn_y(i, j + 1, k - 1))

                kzx = &
                    ISO4(1)*(cn_z(i + 1, j, k) - cn_z(i - 1, j, k)) + &
                    ISO4(2)*( &
                    cn_z(i + 1, j + 1, k) - cn_z(i - 1, j - 1, k) + &
                    cn_z(i + 1, j - 1, k) - cn_z(i - 1, j + 1, k) + &
                    cn_z(i + 1, j, k + 1) - cn_z(i - 1, j, k - 1) + &
                    cn_z(i + 1, j, k - 1) - cn_z(i - 1, j, k + 1))

                kzy = &
                    ISO4(1)*(cn_z(i, j + 1, k) - cn_z(i, j - 1, k)) + &
                    ISO4(2)*( &
                    cn_z(i + 1, j + 1, k) - cn_z(i - 1, j - 1, k) + &
                    cn_z(i - 1, j + 1, k) - cn_z(i + 1, j - 1, k) + &
                    cn_z(i, j + 1, k + 1) - cn_z(i, j - 1, k - 1) + &
                    cn_z(i, j + 1, k - 1) - cn_z(i, j - 1, k + 1))

                curv(i, j, k) = (cn_x(i, j, k)**2 - 1d0)*kxx + (cn_y(i, j, k)**2 - 1d0)*kyy + (cn_z(i, j, k)**2 - 1d0)*kzz + &
         cn_x(i, j, k)*cn_y(i, j, k)*(kxy + kyx) + cn_x(i, j, k)*cn_z(i, j, k)*(kxz + kzx) + cn_y(i, j, k)*cn_z(i, j, k)*(kzy + kyz)

                !curv(i,j,k) = - kxx - kyy - kzz
            end do
        end do
    end do
    !$acc end kernels

    return
end subroutine color_gradient

!=====================================================================================================================================
!------------------ alter solid surface normal directions to control wettability ----------------------
!=====================================================================================================================================
! scheme from Akai et al., 2018
subroutine alter_color_gradient_solid_surface
    use Fluid_multiphase
    use Misc_module
    IMPLICIT NONE
    integer :: i, j, k, num
    real(kind=8) :: nwx, nwy, nwz, theta_local, local_eps
    real(kind=8) :: cnx0, cny0, cnz0, cnxp, cnyp, cnzp, cnxm, cnym, cnzm, coe1, coe2, tmp1, tmp2, theta_tmp, distP, distM

    local_eps = 1d-6

    !$omp parallel do default(none)&
    !$omp & private(i,j,k,nwx,nwy,nwz,cnx0, cny0, cnz0, cnxp, cnyp, cnzp, cnxm, cnym, cnzm, coe1, coe2, tmp1,tmp2,theta_tmp, distP, distM, theta_local) &
    !$omp & shared(fluid_boundary_nodes, cn_x, cn_y, cn_z, c_norm, local_eps,num_fluid_boundary)
    !$acc kernels present(cn_x,cn_y,cn_z,c_norm,fluid_boundary_nodes)
    !$acc loop independent device_type(NVIDIA)
    do num = 1, num_fluid_boundary
        i = fluid_boundary_nodes(num)%ix
        j = fluid_boundary_nodes(num)%iy
        k = fluid_boundary_nodes(num)%iz

        if (c_norm(i, j, k) > local_eps) then
            theta_local = fluid_boundary_nodes(num)%theta
            nwx = fluid_boundary_nodes(num)%nwx
            nwy = fluid_boundary_nodes(num)%nwy
            nwz = fluid_boundary_nodes(num)%nwz
            cnx0 = cn_x(i, j, k)
            cny0 = cn_y(i, j, k)
            cnz0 = cn_z(i, j, k)
            theta_tmp = dacos(nwx*cnx0 + nwy*cny0 + nwz*cnz0)
            tmp1 = 1.0/dsin(theta_tmp)
            tmp2 = dcos(theta_local)
            coe1 = dsin(theta_local)*dcos(theta_tmp)*tmp1
            coe2 = dsin(theta_tmp)*tmp1
            cnxp = (tmp2 - coe1)*nwx + coe2*cnx0
            cnyp = (tmp2 - coe1)*nwy + coe2*cny0
            cnzp = (tmp2 - coe1)*nwz + coe2*cnz0
            cnxm = (tmp2 + coe1)*nwx - coe2*cnx0
            cnym = (tmp2 + coe1)*nwy - coe2*cny0
            cnzm = (tmp2 + coe1)*nwz - coe2*cnz0
            distP = (cnxp - cnx0)*(cnxp - cnx0) + (cnyp - cny0)*(cnyp - cny0) + (cnzp - cnz0)*(cnzp - cnz0)
            distM = (cnxm - cnx0)*(cnxm - cnx0) + (cnym - cny0)*(cnym - cny0) + (cnzm - cnz0)*(cnzm - cnz0)
            if (distP <= distM) then
                cn_x(i, j, k) = cnxp
                cn_y(i, j, k) = cnyp
                cn_z(i, j, k) = cnzp
            else
                cn_x(i, j, k) = cnxm
                cn_y(i, j, k) = cnym
                cn_z(i, j, k) = cnzm
            end if
        end if
    end do
    !$acc end kernels

    return
end subroutine alter_color_gradient_solid_surface

! old iteration scheme from Leclaire et al., 2017
subroutine alter_color_gradient_solid_surface2
  use Fluid_multiphase
  use Misc_module
  IMPLICIT NONE
  integer :: i,j,k,num,iteration,iteration_max
  real(kind=8) :: nwx,nwy,nwz,cos_theta,lambda,local_eps
  real(kind=8) :: vcx0,vcy0,vcz0,vcx1,vcy1,vcz1,vcx2,vcy2,vcz2,err0,err1,err2,tmp

  lambda = 0.5d0
  local_eps = 1d-6
  iteration_max=4
  !$omp parallel do default(none)&
  !$omp & private(i,j,k,nwx,nwy,nwz,vcx0,vcy0,vcz0,vcx1,vcy1,vcz1,vcx2,vcy2,vcz2,tmp,err1,err0,err2,cos_theta,iteration) &
  !$omp & shared(fluid_boundary_nodes, cn_x, cn_y, cn_z, c_norm, local_eps, lambda, iteration_max,num_fluid_boundary)
  !$acc kernels present(cn_x,cn_y,cn_z,c_norm,fluid_boundary_nodes)
  !$acc loop independent device_type(NVIDIA)
  do num=1, num_fluid_boundary
      i=fluid_boundary_nodes(num)%ix
      j=fluid_boundary_nodes(num)%iy
      k=fluid_boundary_nodes(num)%iz

      if(c_norm(i,j,k)>local_eps)then

          cos_theta=dcos(fluid_boundary_nodes(num)%theta)
          nwx = fluid_boundary_nodes(num)%nwx
          nwy = fluid_boundary_nodes(num)%nwy
          nwz = fluid_boundary_nodes(num)%nwz
          vcx0=cn_x(i,j,k)
          vcy0=cn_y(i,j,k)
          vcz0=cn_z(i,j,k)
          vcx1 = vcx0 - lambda*(vcx0+nwx)
          vcy1 = vcy0 - lambda*(vcy0+nwy)
          vcz1 = vcz0 - lambda*(vcz0+nwz)

          err0 = (nwx*vcx0+nwy*vcy0+nwz*vcz0) - cos_theta

          if((dabs(vcx0+nwx)+dabs(vcy0+nwy)+dabs(vcz0+nwz)>local_eps.or.dabs(vcx0-nwx)+dabs(vcy0-nwy)+dabs(vcz0-nwz)>local_eps).and.err0>local_eps)then
              !do not perform alteration when the normal direction of the solid surface aligned with the fluid interface direction, or the initial fluid direction is already the desired direction
              err1 = (nwx*vcx1+nwy*vcy1+nwz*vcz1) - dsqrt(vcx1*vcx1+vcy1*vcy1+vcz1*vcz1)*cos_theta
              tmp = 1d0/(err1-err0)
              vcx2 = tmp*(vcx0*err1 - vcx1*err0)
              vcy2 = tmp*(vcy0*err1 - vcy1*err0)
              vcz2 = tmp*(vcz0*err1 - vcz1*err0)

              err2 =  (nwx*vcx2+nwy*vcy2+nwz*vcz2) - dsqrt(vcx2*vcx2+vcy2*vcy2+vcz2*vcz2)*cos_theta

              if(err2>local_eps)then
                  do iteration = 2, iteration_max
                      vcx0 = vcx1
                      vcy0 = vcy1
                      vcz0 = vcz1
                      vcx1 = vcx2
                      vcy1 = vcy2
                      vcz1 = vcz2
                      err0 = (nwx*vcx0+nwy*vcy0+nwz*vcz0) - dsqrt(vcx0*vcx0+vcy0*vcy0+vcz0*vcz0)*cos_theta
                      err1 = (nwx*vcx1+nwy*vcy1+nwz*vcz1) - dsqrt(vcx1*vcx1+vcy1*vcy1+vcz1*vcz1)*cos_theta
                      tmp = 1d0/(err1-err0)
                      vcx2 = tmp*(vcx0*err1 - vcx1*err0)
                      vcy2 = tmp*(vcy0*err1 - vcy1*err0)
                      vcz2 = tmp*(vcz0*err1 - vcz1*err0)
                      err2 =  (nwx*vcx2+nwy*vcy2+nwz*vcz2) - dsqrt(vcx2*vcx2+vcy2*vcy2+vcz2*vcz2)*cos_theta
                  ! if(err2<1d-6)exit
                  enddo
                  ! if(iteration>=iteration_max)print*,'after',iteration_max,' iterations, theta=', dacos((nwx*vcx2+nwy*vcy2+nwz*vcz2)/dsqrt(vcx2*vcx2+vcy2*vcy2+vcz2*vcz2))/pi*180
              endif
              tmp = 1d0/((1d-30)+dsqrt(vcx2*vcx2+vcy2*vcy2+vcz2*vcz2))
              cn_x(i,j,k)  = vcx2*tmp
              cn_y(i,j,k)  = vcy2*tmp
              cn_z(i,j,k)  = vcz2*tmp
          endif
      endif
  enddo
  !$acc end kernels

  return
end subroutine alter_color_gradient_solid_surface2
!*************alter solid surface color gradient to control wettability**********************************************************

!#ifdef multirange_gradient
!subroutine color_gradient_multirange
!
!    use Fluid_multiphase,  only : c_x,c_y,c_z, phi,ISO8
!    use Misc_module, only : nx,ny,nz
!    IMPLICIT NONE
!    integer :: i,j,k
!
!    !$omp parallel do private(i,j) collapse(2)
!    do k=1,nz
!        do j=1,ny
!            do i=1,nx
!                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cx~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                c_x(i,j,k) = &
!                    ISO8(1)*(phi(i+1,j  ,k  ) - phi(i-1,j  ,k  )) + &
!
!                    ISO8(2)*(&
!                    phi(i+1,j+1,k  ) - phi(i-1,j-1,k  ) + &
!                    phi(i+1,j-1,k  ) - phi(i-1,j+1,k  ) + &
!                    phi(i+1,j  ,k+1) - phi(i-1,j  ,k-1) + &
!                    phi(i+1,j  ,k-1) - phi(i-1,j  ,k+1)) + &
!
!                    ISO8(3)*(&
!                    phi(i+1,j+1,k+1) - phi(i-1,j-1,k-1) + &
!                    phi(i+1,j+1,k-1) - phi(i-1,j-1,k+1) + &
!                    phi(i+1,j-1,k+1) - phi(i-1,j+1,k-1) + &
!                    phi(i+1,j-1,k-1) - phi(i-1,j+1,k+1)) + &
!
!                    ISO8(4)*(phi(i+2,j  ,k  ) - phi(i-2,j  ,k  )) + &
!
!                    ISO8(5)*(&
!                    phi(i+2,j+1,k  ) - phi(i-2,j-1,k  ) + &
!                    phi(i+2,j-1,k  ) - phi(i-2,j+1,k  ) + &
!                    phi(i+2,j  ,k+1) - phi(i-2,j  ,k-1) + &
!                    phi(i+2,j  ,k-1) - phi(i-2,j  ,k+1) + &
!                    phi(i+1,j+2,k  ) - phi(i-1,j-2,k  ) + &
!                    phi(i+1,j-2,k  ) - phi(i-1,j+2,k  ) + &
!                    phi(i+1,j  ,k+2) - phi(i-1,j  ,k-2) + &
!                    phi(i+1,j  ,k-2) - phi(i-1,j  ,k+2)) + &
!
!                    ISO8(6)*(&
!                    phi(i+2,j+1,k+1) - phi(i-2,j-1,k-1) + &
!                    phi(i+2,j+1,k-1) - phi(i-2,j-1,k+1) + &
!                    phi(i+2,j-1,k+1) - phi(i-2,j+1,k-1) + &
!                    phi(i+2,j-1,k-1) - phi(i-2,j+1,k+1) + &
!                    phi(i+1,j+2,k+1) - phi(i-1,j-2,k-1) + &
!                    phi(i+1,j+2,k-1) - phi(i-1,j-2,k+1) + &
!                    phi(i+1,j-2,k+1) - phi(i-1,j+2,k-1) + &
!                    phi(i+1,j-2,k-1) - phi(i-1,j+2,k+1) + &
!                    phi(i+1,j+1,k+2) - phi(i-1,j-1,k-2) + &
!                    phi(i+1,j+1,k-2) - phi(i-1,j-1,k+2) + &
!                    phi(i+1,j-1,k+2) - phi(i-1,j+1,k-2) + &
!                    phi(i+1,j-1,k-2) - phi(i-1,j+1,k+2)) + &
!
!                    ISO8(7)*(&
!                    phi(i+2,j+2,k  ) - phi(i-2,j-2,k  ) + &
!                    phi(i+2,j-2,k  ) - phi(i-2,j+2,k  ) + &
!                    phi(i+2,j  ,k+2) - phi(i-2,j  ,k-2) + &
!                    phi(i+2,j  ,k-2) - phi(i-2,j  ,k+2))
!
!                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cy~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                c_y(i,j,k) = &
!                    ISO8(1)*(phi(i  ,j+1,k  ) - phi(i  ,j-1,k  )) + &
!
!                    ISO8(2)*(&
!                    phi(i+1,j+1,k  ) - phi(i-1,j-1,k  ) + &
!                    phi(i-1,j+1,k  ) - phi(i+1,j-1,k  ) + &
!                    phi(i  ,j+1,k+1) - phi(i  ,j-1,k-1) + &
!                    phi(i  ,j+1,k-1) - phi(i  ,j-1,k+1)) + &
!
!                    ISO8(3)*(&
!                    phi(i+1,j+1,k+1) - phi(i-1,j-1,k-1) + &
!                    phi(i+1,j+1,k-1) - phi(i-1,j-1,k+1) + &
!                    phi(i-1,j+1,k-1) - phi(i+1,j-1,k+1) + &
!                    phi(i-1,j+1,k+1) - phi(i+1,j-1,k-1)) + &
!
!                    ISO8(4)*(phi(i  ,j+2,k  ) - phi(i  ,j-2,k  )) + &
!
!                    ISO8(5)*(&
!                    phi(i+2,j+1,k  ) - phi(i-2,j-1,k  ) + &
!                    phi(i-2,j+1,k  ) - phi(i+2,j-1,k  ) + &
!                    phi(i  ,j+2,k+1) - phi(i  ,j-2,k-1) + &
!                    phi(i  ,j+2,k-1) - phi(i  ,j-2,k+1) + &
!                    phi(i+1,j+2,k  ) - phi(i-1,j-2,k  ) + &
!                    phi(i-1,j+2,k  ) - phi(i+1,j-2,k  ) + &
!                    phi(i  ,j+1,k+2) - phi(i  ,j-1,k-2) + &
!                    phi(i  ,j+1,k-2) - phi(i  ,j-1,k+2)) + &
!
!                    ISO8(6)*(&
!                    phi(i+2,j+1,k+1) - phi(i-2,j-1,k-1) + &
!                    phi(i+2,j+1,k-1) - phi(i-2,j-1,k+1) + &
!                    phi(i-2,j+1,k+1) - phi(i+2,j-1,k-1) + &
!                    phi(i-2,j+1,k-1) - phi(i+2,j-1,k+1) + &
!                    phi(i+1,j+2,k+1) - phi(i-1,j-2,k-1) + &
!                    phi(i+1,j+2,k-1) - phi(i-1,j-2,k+1) + &
!                    phi(i-1,j+2,k+1) - phi(i+1,j-2,k-1) + &
!                    phi(i-1,j+2,k-1) - phi(i+1,j-2,k+1) + &
!                    phi(i+1,j+1,k+2) - phi(i-1,j-1,k-2) + &
!                    phi(i+1,j+1,k-2) - phi(i-1,j-1,k+2) + &
!                    phi(i-1,j+1,k+2) - phi(i+1,j-1,k-2) + &
!                    phi(i-1,j+1,k-2) - phi(i+1,j-1,k+2)) + &
!
!                    ISO8(7)*(&
!                    phi(i+2,j+2,k  ) - phi(i-2,j-2,k  ) + &
!                    phi(i-2,j+2,k  ) - phi(i+2,j-2,k  ) + &
!                    phi(i  ,j+2,k+2) - phi(i  ,j-2,k-2) + &
!                    phi(i  ,j+2,k-2) - phi(i  ,j-2,k+2))
!                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cz~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                c_z(i,j,k) = &
!                    ISO8(1)*(phi(i  ,j  ,k+1) -phi(i  ,j  ,k-1)) + &
!
!                    ISO8(2)*(&
!                    phi(i  ,j+1,k+1) - phi(i  ,j-1,k-1) + &
!                    phi(i  ,j-1,k+1) - phi(i  ,j+1,k-1) + &
!                    phi(i+1,j  ,k+1) - phi(i-1,j  ,k-1) + &
!                    phi(i-1,j  ,k+1) - phi(i+1,j  ,k-1)) + &
!
!                    +ISO8(3)*(&
!                    phi(i+1,j+1,k+1) - phi(i-1,j-1,k-1) + &
!                    phi(i+1,j-1,k+1) - phi(i-1,j+1,k-1) + &
!                    phi(i-1,j+1,k+1) - phi(i+1,j-1,k-1) + &
!                    phi(i-1,j-1,k+1) - phi(i+1,j+1,k-1)) + &
!
!                    ISO8(4)*(phi(i  ,j  ,k+2) - phi(i  ,j  ,k-2)) + &
!
!                    ISO8(5)*(&
!                    phi(i  ,j+1,k+2) - phi(i  ,j-1,k-2) + &
!                    phi(i  ,j-1,k+2) - phi(i  ,j+1,k-2) + &
!                    phi(i+2,j  ,k+1) - phi(i-2,j  ,k-1) + &
!                    phi(i-2,j  ,k+1) - phi(i+2,j  ,k-1) + &
!                    phi(i  ,j+2,k+1) - phi(i  ,j-2,k-1) + &
!                    phi(i  ,j-2,k+1) - phi(i  ,j+2,k-1) + &
!                    phi(i+1,j  ,k+2) - phi(i-1,j  ,k-2) + &
!                    phi(i-1,j  ,k+2) - phi(i+1,j  ,k-2)) + &
!
!                    ISO8(6)*(&
!                    phi(i+2,j+1,k+1) - phi(i-2,j-1,k-1) + &
!                    phi(i+2,j-1,k+1) - phi(i-2,j+1,k-1) + &
!                    phi(i-2,j+1,k+1) - phi(i+2,j-1,k-1) + &
!                    phi(i-2,j-1,k+1) - phi(i+2,j+1,k-1) + &
!                    phi(i+1,j+2,k+1) - phi(i-1,j-2,k-1) + &
!                    phi(i+1,j-2,k+1) - phi(i-1,j+2,k-1) + &
!                    phi(i-1,j+2,k+1) - phi(i+1,j-2,k-1) + &
!                    phi(i-1,j-2,k+1) - phi(i+1,j+2,k-1) + &
!                    phi(i+1,j+1,k+2) - phi(i-1,j-1,k-2) + &
!                    phi(i+1,j-1,k+2) - phi(i-1,j+1,k-2) + &
!                    phi(i-1,j+1,k+2) - phi(i+1,j-1,k-2) + &
!                    phi(i-1,j-1,k+2) - phi(i+1,j+1,k-2)) + &
!
!                    ISO8(7)*(&
!                    phi(i  ,j+2,k+2) - phi(i  ,j-2,k-2) + &
!                    phi(i  ,j-2,k+2) - phi(i  ,j+2,k-2) + &
!                    phi(i+2,j  ,k+2) - phi(i-2,j  ,k-2) + &
!                    phi(i-2,j  ,k+2) - phi(i+2,j  ,k-2))
!
!            enddo
!        enddo
!    enddo
!
!    return
!endsubroutine color_gradient_multirange
!#endif