!=====================================================================================================================================
!---------------------- odd step kernel ----------------------
! complete two streaming steps
!=====================================================================================================================================
subroutine kernel_odd(ixmin, ixmax, iymin, iymax, izmin, izmax, async_label)
    use Misc_module
    use Fluid_singlephase
    IMPLICIT NONE
    integer :: i, j, k, m, ixmin, ixmax, iymin, iymax, izmin, izmax, async_label
    real(kind=8) :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9
    real(kind=8) :: m_rho,m_e, m_e2, m_jx, m_qx, m_jy, m_qy, m_jz, m_qz, m_3pxx, m_3pixx, m_pww, m_piww, m_pxy, m_pyz, m_pzx, m_tx, m_ty, m_tz
    real(kind=8) :: ft0, ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8, ft9, ft10, ft11, ft12, ft13, ft14, ft15, ft16, ft17, ft18
    real(kind=8) :: ux, uy, uz, den, u2, fx, fy, fz

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~openmp
    !$omp PARALLEL DO default(none) &
    !$omp & SHARED(&
    !$omp & walls,f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18, &
    !$omp & ixmin,ixmax,iymin,iymax,izmin,izmax,force_Z,s_e,s_e2,s_q,s_nu,s_pi,s_t) &
    !$omp & PRIVATE(&
    !$omp & i,j,&
    !$omp & ux,uy,uz,u2,den,fx,fy,fz,&
    !$omp & ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft17,ft18,&
    !$omp & sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, &
    !$omp & m_rho, m_e, m_e2, m_jx, m_qx, m_jy, m_qy, m_jz, m_qz, m_3pxx, m_3pixx, m_pww, m_piww, m_pxy, m_pyz, m_pzx, m_tx, m_ty, m_tz) collapse(2)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
    !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,walls)async(async_label)
    !$acc loop collapse(3) device_type(NVIDIA)
    do k = izmin, izmax
        DO j = iymin, iymax
            !DIR$ IVDEP
            do i = ixmin, ixmax
                if(walls(i, j, k)==0)then
                  !+++++++++- AA pattern pull step++++++++++++
                  ft0 = f0(i, j, k)
                  ft1 = f1(i - 1, j, k)
                  ft2 = f2(i + 1, j, k)
                  ft3 = f3(i, j - 1, k)
                  ft4 = f4(i, j + 1, k)
                  ft5 = f5(i, j, k - 1)
                  ft6 = f6(i, j, k + 1)
                  ft7 = f7(i - 1, j - 1, k)
                  ft8 = f8(i + 1, j - 1, k)
                  ft9 = f9(i - 1, j + 1, k)
                  ft10 = f10(i + 1, j + 1, k)
                  ft11 = f11(i - 1, j, k - 1)
                  ft12 = f12(i + 1, j, k - 1)
                  ft13 = f13(i - 1, j, k + 1)
                  ft14 = f14(i + 1, j, k + 1)
                  ft15 = f15(i, j - 1, k - 1)
                  ft16 = f16(i, j + 1, k - 1)
                  ft17 = f17(i, j - 1, k + 1)
                  ft18 = f18(i, j + 1, k + 1)

                  fx = 0d0
                  fy = 0d0
                  fz = force_Z   !body force force_z along flow direction

                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MRT kernel, repeated part~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  !+++++++++- !calculate macroscopic variables++++++++++++
        den = ft0 + ft1 + ft2 + ft3 + ft4 + ft5 + ft6 + ft7 + ft8 + ft9 + ft10 + ft11 + ft12 + ft13 + ft14 + ft15 + ft16 + ft17 + ft18
                  ux = ft1 - ft2 + ft7 - ft8 + ft9 - ft10 + ft11 - ft12 + ft13 - ft14 + 0.5d0*fx
                  uy = ft3 - ft4 + ft7 + ft8 - ft9 - ft10 + ft15 - ft16 + ft17 - ft18 + 0.5d0*fy
                  uz = ft5 - ft6 + ft11 + ft12 - ft13 - ft14 + ft15 + ft16 - ft17 - ft18 + 0.5d0*fz

                  u2 = ux*ux + uy*uy + uz*uz

                  !PDFs summations for computation efficiency purpose
                  sum1 = ft1 + ft2 + ft3 + ft4 + ft5 + ft6
                  sum2 = ft7 + ft8 + ft9 + ft10 + ft11 + ft12 + ft13 + ft14 + ft15 + ft16 + ft17 + ft18
                  sum3 = ft7 - ft8 + ft9 - ft10 + ft11 - ft12 + ft13 - ft14
                  sum4 = ft7 + ft8 - ft9 - ft10 + ft15 - ft16 + ft17 - ft18
                  sum5 = ft11 + ft12 - ft13 - ft14 + ft15 + ft16 - ft17 - ft18
                  sum6 = 2d0*(ft1 + ft2) - ft3 - ft4 - ft5 - ft6
                  sum7 = ft7 + ft8 + ft9 + ft10 + ft11 + ft12 + ft13 + ft14 - 2d0*(ft15 + ft16 + ft17 + ft18)
                  sum8 = ft3 + ft4 - ft5 - ft6
                  sum9 = ft7 + ft8 + ft9 + ft10 - ft11 - ft12 - ft13 - ft14
                  !PDF to moment
                  m_rho = den
                  m_e = -30d0*ft0 - 11d0*sum1 + 8d0*sum2
                  m_e2 = 12d0*ft0 - 4d0*sum1 + sum2
                  m_jx = ft1 - ft2 + sum3
                  m_qx = -4d0*(ft1 - ft2) + sum3
                  m_jy = ft3 - ft4 + sum4
                  m_qy = -4d0*(ft3 - ft4) + sum4
                  m_jz = ft5 - ft6 + sum5
                  m_qz = -4d0*(ft5 - ft6) + sum5
                  m_3pxx = sum6 + sum7
                  m_3pixx = -2d0*sum6 + sum7
                  m_pww = sum8 + sum9
                  m_piww = -2d0*sum8 + sum9
                  m_pxy = ft7 - ft8 - ft9 + ft10
                  m_pyz = ft15 - ft16 - ft17 + ft18
                  m_pzx = ft11 - ft12 - ft13 + ft14
                  m_tx = ft7 - ft8 + ft9 - ft10 - ft11 + ft12 - ft13 + ft14
                  m_ty = -ft7 - ft8 + ft9 + ft10 + ft15 - ft16 + ft17 - ft18
                  m_tz = ft11 + ft12 - ft13 - ft14 - ft15 - ft16 + ft17 + ft18

                  !relaxtion in moment space
          m_e = m_e - s_e*(m_e - (-11.0d0*den + 19.0d0*u2)) + (38d0 - 19d0*s_e)*(fx*ux + fy*uy + fz*uz)                           !m1
  m_e2 = m_e2 - s_e2*(m_e2 - (mrt_e2_coef1*den + mrt_e2_coef2*u2)) + (-11d0 + 5.5d0*s_e2)*(fx*ux + fy*uy + fz*uz)                     !m2
                  m_jx = m_jx + fx                                                                                                   !m3
                  m_qx = m_qx - s_q*(m_qx - (-0.666666666666666667d0*ux)) + (-0.666666666666666667d0 + 0.333333333333333333d0*s_q)*fx !m4
                  m_jy = m_jy + fy                                                                                                   !m5
                  m_qy = m_qy - s_q*(m_qy - (-0.666666666666666667d0*uy)) + (-0.666666666666666667d0 + 0.333333333333333333d0*s_q)*fy !m6
                  m_jz = m_jz + fz                                                                                                   !m7
                  m_qz = m_qz - s_q*(m_qz - (-0.666666666666666667d0*uz)) + (-0.666666666666666667d0 + 0.333333333333333333d0*s_q)*fz !m8

          m_3pxx = m_3pxx - s_nu*(m_3pxx - (3d0*ux*ux - u2)) + (2d0 - s_nu)*(2d0*fx*ux - fy*uy - fz*uz)                            !m9
    m_3pixx = m_3pixx - s_pi*(m_3pixx - mrt_omega_xx*(3d0*ux*ux - u2)) + (1d0 - 0.5d0*s_pi)*(-2d0*fx*ux + fy*uy + fz*uz)         !m10
                  m_pww = m_pww - s_nu*(m_pww - (uy*uy - uz*uz)) + (2d0 - s_nu)*(fy*uy - fz*uz)                                      !m11
          m_piww = m_piww - s_pi*(m_piww - mrt_omega_xx*(uy*uy - uz*uz)) + (1d0 - 0.5d0*s_pi)*(-fy*uy + fz*uz)                   !m12
                  m_pxy = m_pxy - s_nu*(m_pxy - (ux*uy)) + (1d0 - 0.5d0*s_nu)*(fx*uy + fy*ux)                                        !m13
                  m_pyz = m_pyz - s_nu*(m_pyz - (uy*uz)) + (1d0 - 0.5d0*s_nu)*(fy*uz + fz*uy)                                        !m14
                  m_pzx = m_pzx - s_nu*(m_pzx - (ux*uz)) + (1d0 - 0.5d0*s_nu)*(fx*uz + fz*ux)                                        !m15
                  m_tx = m_tx - s_t*(m_tx)                                                                                         !m16
                  m_ty = m_ty - s_t*(m_ty)                                                                                         !m17
                  m_tz = m_tz - s_t*(m_tz)                                                                                         !m18

                  !transform back to PDFs
                  !coeffcients for performance
                  m_rho = mrt_coef1*m_rho            ! 1/19
                  m_e = mrt_coef2*m_e              ! 1/2394
                  m_e2 = mrt_coef3*m_e2              ! 1/252
                  m_jx = 0.1d0*m_jx
                  m_qx = 0.025d0*m_qx
                  m_jy = 0.1d0*m_jy
                  m_qy = 0.025d0*m_qy
                  m_jz = 0.1d0*m_jz
                  m_qz = 0.025d0*m_qz
                  m_3pxx = 2d0*mrt_coef4*m_3pxx      ! 1/36
                  m_3pixx = mrt_coef4*m_3pixx        ! 1/72
                  m_pww = 6d0*mrt_coef4*m_pww       ! 1/12
                  m_piww = 3d0*mrt_coef4*m_piww      ! 1/24
                  m_pxy = 0.25d0*m_pxy
                  m_pyz = 0.25d0*m_pyz
                  m_pzx = 0.25d0*m_pzx
                  m_tx = 0.125d0*m_tx
                  m_ty = 0.125d0*m_ty
                  m_tz = 0.125d0*m_tz
                  sum1 = m_rho - 11d0*m_e - 4d0*m_e2
                  sum2 = 2d0*m_3pxx - 4d0*m_3pixx
                  sum3 = m_pww - 2d0*m_piww
                  sum4 = m_rho + 8d0*m_e + m_e2
                  sum5 = m_jx + m_qx
                  sum6 = m_jy + m_qy
                  sum7 = m_jz + m_qz
                  sum8 = m_3pxx + m_3pixx
                  sum9 = m_pww + m_piww

                  ft0 = m_rho - 30d0*m_e + 12d0*m_e2
                  ft1 = sum1 + m_jx - 4d0*m_qx + sum2
                  ft2 = sum1 - m_jx + 4d0*m_qx + sum2 
                  ft3 = sum1 + m_jy - 4d0*m_qy - 0.5d0*sum2 + sum3 
                  ft4 = sum1 - m_jy + 4d0*m_qy - 0.5d0*sum2 + sum3 
                  ft5 = sum1 + m_jz - 4d0*m_qz - 0.5d0*sum2 - sum3 
                  ft6 = sum1 - m_jz + 4d0*m_qz - 0.5d0*sum2 - sum3 
            ft7 = sum4 + sum5 + sum6 + sum8 + sum9 + m_pxy + m_tx - m_ty 
              ft8 = sum4 - sum5 + sum6 + sum8 + sum9 - m_pxy - m_tx - m_ty 
              ft9 = sum4 + sum5 - sum6 + sum8 + sum9 - m_pxy + m_tx + m_ty 
            ft10 = sum4 - sum5 - sum6 + sum8 + sum9 + m_pxy - m_tx + m_ty 
            ft11 = sum4 + sum5 + sum7 + sum8 - sum9 + m_pzx - m_tx + m_tz 
            ft12 = sum4 - sum5 + sum7 + sum8 - sum9 - m_pzx + m_tx + m_tz 
            ft13 = sum4 + sum5 - sum7 + sum8 - sum9 - m_pzx - m_tx - m_tz 
            ft14 = sum4 - sum5 - sum7 + sum8 - sum9 + m_pzx + m_tx - m_tz 
              ft15 = sum4 + sum6 + sum7 - sum8*2d0 + m_pyz + m_ty - m_tz 
              ft16 = sum4 - sum6 + sum7 - sum8*2d0 - m_pyz - m_ty - m_tz 
              ft17 = sum4 + sum6 - sum7 - sum8*2d0 - m_pyz + m_ty + m_tz 
              ft18 = sum4 - sum6 - sum7 - sum8*2d0 + m_pyz - m_ty + m_tz 
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MRT kernel, repeated part~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                  !+++++++++- AA pattern push step++++++++++++
                  f0(i, j, k) = ft0
                  f2(i + 1, j, k) = ft1
                  f1(i - 1, j, k) = ft2
                  f4(i, j + 1, k) = ft3
                  f3(i, j - 1, k) = ft4
                  f6(i, j, k + 1) = ft5
                  f5(i, j, k - 1) = ft6
                  f10(i + 1, j + 1, k) = ft7
                  f9(i - 1, j + 1, k) = ft8
                  f8(i + 1, j - 1, k) = ft9
                  f7(i - 1, j - 1, k) = ft10
                  f14(i + 1, j, k + 1) = ft11
                  f13(i - 1, j, k + 1) = ft12
                  f12(i + 1, j, k - 1) = ft13
                  f11(i - 1, j, k - 1) = ft14
                  f18(i, j + 1, k + 1) = ft15
                  f17(i, j - 1, k + 1) = ft16
                  f16(i, j + 1, k - 1) = ft17
                  f15(i, j - 1, k - 1) = ft18
                endif

            END DO
        END DO
    end do
    !$acc end kernels

    RETURN
end subroutine kernel_odd

!=====================================================================================================================================
!---------------------- even step kernel ----------------------
! no steaming steps
!=====================================================================================================================================
subroutine kernel_even(ixmin, ixmax, iymin, iymax, izmin, izmax, async_label)
    use Misc_module
    use Fluid_singlephase
    IMPLICIT NONE
    integer :: i, j, k, m, ixmin, ixmax, iymin, iymax, izmin, izmax, async_label
    integer(kind=1) :: wall_indicator
    real(kind=8) :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9
    real(kind=8) :: m_rho,m_e, m_e2, m_jx, m_qx, m_jy, m_qy, m_jz, m_qz, m_3pxx, m_3pixx, m_pww, m_piww, m_pxy, m_pyz, m_pzx, m_tx, m_ty, m_tz
    real(kind=8) :: ft0, ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8, ft9, ft10, ft11, ft12, ft13, ft14, ft15, ft16, ft17, ft18
    real(kind=8) :: ux, uy, uz, den, u2, fx, fy, fz

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~openmp
    !$omp PARALLEL DO default(none) &
    !$omp & SHARED(&
    !$omp & walls,f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18, &
    !$omp & ixmin,ixmax,iymin,iymax,izmin,izmax,force_Z,s_e,s_e2,s_q,s_nu,s_pi,s_t) &
    !$omp & PRIVATE(&
    !$omp & i,j,wall_indicator,&
    !$omp & ux,uy,uz,u2,den,fx,fy,fz,&
    !$omp & ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft17,ft18,&
    !$omp & sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, &
    !$omp & m_rho, m_e, m_e2, m_jx, m_qx, m_jy, m_qy, m_jz, m_qz, m_3pxx, m_3pixx, m_pww, m_piww, m_pxy, m_pyz, m_pzx, m_tx, m_ty, m_tz) collapse(2)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
    !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,walls)async(async_label)
    !$acc loop collapse(3) device_type(NVIDIA)
    do k = izmin, izmax
        DO j = iymin, iymax
            !DIR$ IVDEP
            do i = ixmin, ixmax
                if(walls(i,j,k)==0)then
                  !+++++++++- AA pattern pull step++++++++++++
                  ft0 = f0(i, j, k)
                  ft1 = f2(i, j, k)
                  ft2 = f1(i, j, k)
                  ft3 = f4(i, j, k)
                  ft4 = f3(i, j, k)
                  ft5 = f6(i, j, k)
                  ft6 = f5(i, j, k)
                  ft7 = f10(i, j, k)
                  ft8 = f9(i, j, k)
                  ft9 = f8(i, j, k)
                  ft10 = f7(i, j, k)
                  ft11 = f14(i, j, k)
                  ft12 = f13(i, j, k)
                  ft13 = f12(i, j, k)
                  ft14 = f11(i, j, k)
                  ft15 = f18(i, j, k)
                  ft16 = f17(i, j, k)
                  ft17 = f16(i, j, k)
                  ft18 = f15(i, j, k)

                  fx = 0d0
                  fy = 0d0
                  fz = force_Z   !body force force_z along flow direction

                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MRT kernel, repeated part~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  !+++++++++- !calculate macroscopic variables++++++++++++
        den = ft0 + ft1 + ft2 + ft3 + ft4 + ft5 + ft6 + ft7 + ft8 + ft9 + ft10 + ft11 + ft12 + ft13 + ft14 + ft15 + ft16 + ft17 + ft18
                  ux = ft1 - ft2 + ft7 - ft8 + ft9 - ft10 + ft11 - ft12 + ft13 - ft14 + 0.5d0*fx
                  uy = ft3 - ft4 + ft7 + ft8 - ft9 - ft10 + ft15 - ft16 + ft17 - ft18 + 0.5d0*fy
                  uz = ft5 - ft6 + ft11 + ft12 - ft13 - ft14 + ft15 + ft16 - ft17 - ft18 + 0.5d0*fz
                  u2 = ux*ux + uy*uy + uz*uz

                  !PDFs summations for computation efficiency purpose
                  sum1 = ft1 + ft2 + ft3 + ft4 + ft5 + ft6
                  sum2 = ft7 + ft8 + ft9 + ft10 + ft11 + ft12 + ft13 + ft14 + ft15 + ft16 + ft17 + ft18
                  sum3 = ft7 - ft8 + ft9 - ft10 + ft11 - ft12 + ft13 - ft14
                  sum4 = ft7 + ft8 - ft9 - ft10 + ft15 - ft16 + ft17 - ft18
                  sum5 = ft11 + ft12 - ft13 - ft14 + ft15 + ft16 - ft17 - ft18
                  sum6 = 2d0*(ft1 + ft2) - ft3 - ft4 - ft5 - ft6
                  sum7 = ft7 + ft8 + ft9 + ft10 + ft11 + ft12 + ft13 + ft14 - 2d0*(ft15 + ft16 + ft17 + ft18)
                  sum8 = ft3 + ft4 - ft5 - ft6
                  sum9 = ft7 + ft8 + ft9 + ft10 - ft11 - ft12 - ft13 - ft14
                  !PDF to moment
                  m_rho = den
                  m_e = -30d0*ft0 - 11d0*sum1 + 8d0*sum2
                  m_e2 = 12d0*ft0 - 4d0*sum1 + sum2
                  m_jx = ft1 - ft2 + sum3
                  m_qx = -4d0*(ft1 - ft2) + sum3
                  m_jy = ft3 - ft4 + sum4
                  m_qy = -4d0*(ft3 - ft4) + sum4
                  m_jz = ft5 - ft6 + sum5
                  m_qz = -4d0*(ft5 - ft6) + sum5
                  m_3pxx = sum6 + sum7
                  m_3pixx = -2d0*sum6 + sum7
                  m_pww = sum8 + sum9
                  m_piww = -2d0*sum8 + sum9
                  m_pxy = ft7 - ft8 - ft9 + ft10
                  m_pyz = ft15 - ft16 - ft17 + ft18
                  m_pzx = ft11 - ft12 - ft13 + ft14
                  m_tx = ft7 - ft8 + ft9 - ft10 - ft11 + ft12 - ft13 + ft14
                  m_ty = -ft7 - ft8 + ft9 + ft10 + ft15 - ft16 + ft17 - ft18
                  m_tz = ft11 + ft12 - ft13 - ft14 - ft15 - ft16 + ft17 + ft18

                  !relaxtion in moment space
          m_e = m_e - s_e*(m_e - (-11.0d0*den + 19.0d0*u2)) + (38d0 - 19d0*s_e)*(fx*ux + fy*uy + fz*uz)                           !m1
  m_e2 = m_e2 - s_e2*(m_e2 - (mrt_e2_coef1*den + mrt_e2_coef2*u2)) + (-11d0 + 5.5d0*s_e2)*(fx*ux + fy*uy + fz*uz)                     !m2
                  m_jx = m_jx + fx                                                                                                   !m3
                  m_qx = m_qx - s_q*(m_qx - (-0.666666666666666667d0*ux)) + (-0.666666666666666667d0 + 0.333333333333333333d0*s_q)*fx !m4
                  m_jy = m_jy + fy                                                                                                   !m5
                  m_qy = m_qy - s_q*(m_qy - (-0.666666666666666667d0*uy)) + (-0.666666666666666667d0 + 0.333333333333333333d0*s_q)*fy !m6
                  m_jz = m_jz + fz                                                                                                   !m7
                  m_qz = m_qz - s_q*(m_qz - (-0.666666666666666667d0*uz)) + (-0.666666666666666667d0 + 0.333333333333333333d0*s_q)*fz !m8

          m_3pxx = m_3pxx - s_nu*(m_3pxx - (3d0*ux*ux - u2)) + (2d0 - s_nu)*(2d0*fx*ux - fy*uy - fz*uz)                            !m9
    m_3pixx = m_3pixx - s_pi*(m_3pixx - mrt_omega_xx*(3d0*ux*ux - u2)) + (1d0 - 0.5d0*s_pi)*(-2d0*fx*ux + fy*uy + fz*uz)         !m10
                  m_pww = m_pww - s_nu*(m_pww - (uy*uy - uz*uz)) + (2d0 - s_nu)*(fy*uy - fz*uz)                                      !m11
          m_piww = m_piww - s_pi*(m_piww - mrt_omega_xx*(uy*uy - uz*uz)) + (1d0 - 0.5d0*s_pi)*(-fy*uy + fz*uz)                   !m12
                  m_pxy = m_pxy - s_nu*(m_pxy - (ux*uy)) + (1d0 - 0.5d0*s_nu)*(fx*uy + fy*ux)                                        !m13
                  m_pyz = m_pyz - s_nu*(m_pyz - (uy*uz)) + (1d0 - 0.5d0*s_nu)*(fy*uz + fz*uy)                                        !m14
                  m_pzx = m_pzx - s_nu*(m_pzx - (ux*uz)) + (1d0 - 0.5d0*s_nu)*(fx*uz + fz*ux)                                        !m15
                  m_tx = m_tx - s_t*(m_tx)                                                                                         !m16
                  m_ty = m_ty - s_t*(m_ty)                                                                                         !m17
                  m_tz = m_tz - s_t*(m_tz)                                                                                         !m18

                  !transform back to PDFs
                  !coeffcients for performance
                  m_rho = mrt_coef1*m_rho            ! 1/19
                  m_e = mrt_coef2*m_e              ! 1/2394
                  m_e2 = mrt_coef3*m_e2              ! 1/252
                  m_jx = 0.1d0*m_jx
                  m_qx = 0.025d0*m_qx
                  m_jy = 0.1d0*m_jy
                  m_qy = 0.025d0*m_qy
                  m_jz = 0.1d0*m_jz
                  m_qz = 0.025d0*m_qz
                  m_3pxx = 2d0*mrt_coef4*m_3pxx      ! 1/36
                  m_3pixx = mrt_coef4*m_3pixx        ! 1/72
                  m_pww = 6d0*mrt_coef4*m_pww       ! 1/12
                  m_piww = 3d0*mrt_coef4*m_piww      ! 1/24
                  m_pxy = 0.25d0*m_pxy
                  m_pyz = 0.25d0*m_pyz
                  m_pzx = 0.25d0*m_pzx
                  m_tx = 0.125d0*m_tx
                  m_ty = 0.125d0*m_ty
                  m_tz = 0.125d0*m_tz
                  sum1 = m_rho - 11d0*m_e - 4d0*m_e2
                  sum2 = 2d0*m_3pxx - 4d0*m_3pixx
                  sum3 = m_pww - 2d0*m_piww
                  sum4 = m_rho + 8d0*m_e + m_e2
                  sum5 = m_jx + m_qx
                  sum6 = m_jy + m_qy
                  sum7 = m_jz + m_qz
                  sum8 = m_3pxx + m_3pixx
                  sum9 = m_pww + m_piww

                  ft0 = m_rho - 30d0*m_e + 12d0*m_e2
                  ft1 = sum1 + m_jx - 4d0*m_qx + sum2 
                  ft2 = sum1 - m_jx + 4d0*m_qx + sum2 
                  ft3 = sum1 + m_jy - 4d0*m_qy - 0.5d0*sum2 + sum3 
                  ft4 = sum1 - m_jy + 4d0*m_qy - 0.5d0*sum2 + sum3 
                  ft5 = sum1 + m_jz - 4d0*m_qz - 0.5d0*sum2 - sum3 
                  ft6 = sum1 - m_jz + 4d0*m_qz - 0.5d0*sum2 - sum3 
                  ft7 = sum4 + sum5 + sum6 + sum8 + sum9 + m_pxy + m_tx - m_ty 
                  ft8 = sum4 - sum5 + sum6 + sum8 + sum9 - m_pxy - m_tx - m_ty 
                  ft9 = sum4 + sum5 - sum6 + sum8 + sum9 - m_pxy + m_tx + m_ty 
                  ft10 = sum4 - sum5 - sum6 + sum8 + sum9 + m_pxy - m_tx + m_ty 
                  ft11 = sum4 + sum5 + sum7 + sum8 - sum9 + m_pzx - m_tx + m_tz 
                  ft12 = sum4 - sum5 + sum7 + sum8 - sum9 - m_pzx + m_tx + m_tz 
                  ft13 = sum4 + sum5 - sum7 + sum8 - sum9 - m_pzx - m_tx - m_tz 
                  ft14 = sum4 - sum5 - sum7 + sum8 - sum9 + m_pzx + m_tx - m_tz 
                  ft15 = sum4 + sum6 + sum7 - sum8*2d0 + m_pyz + m_ty - m_tz 
                  ft16 = sum4 - sum6 + sum7 - sum8*2d0 - m_pyz - m_ty - m_tz 
                  ft17 = sum4 + sum6 - sum7 - sum8*2d0 - m_pyz + m_ty + m_tz 
                  ft18 = sum4 - sum6 - sum7 - sum8*2d0 + m_pyz - m_ty + m_tz 
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MRT kernel, repeated part~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                  !+++++++++- AA pattern++++++++++++
                  f0(i, j, k) = ft0
                  f1(i, j, k) = ft1
                  f2(i, j, k) = ft2
                  f3(i, j, k) = ft3
                  f4(i, j, k) = ft4
                  f5(i, j, k) = ft5
                  f6(i, j, k) = ft6
                  f7(i, j, k) = ft7
                  f8(i, j, k) = ft8
                  f9(i, j, k) = ft9
                  f10(i, j, k) = ft10
                  f11(i, j, k) = ft11
                  f12(i, j, k) = ft12
                  f13(i, j, k) = ft13
                  f14(i, j, k) = ft14
                  f15(i, j, k) = ft15
                  f16(i, j, k) = ft16
                  f17(i, j, k) = ft17
                  f18(i, j, k) = ft18
                endif
            END DO
        END DO
    end do
    !$acc end kernels

    RETURN
end subroutine kernel_even
