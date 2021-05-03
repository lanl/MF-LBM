#include "./preprocessor.h"
!=====================================================================================================================================
!---------------------- odd step kernel ----------------------
! complete two streaming steps
!=====================================================================================================================================
subroutine kernel_odd_color(ixmin,ixmax,iymin,iymax,izmin,izmax,async_label)
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    IMPLICIT NONE
    integer :: i,j,k,m, ixmin,ixmax,iymin,iymax,izmin,izmax,async_label
    integer(kind=1) :: wall_indicator
    real(kind=8) :: fx,fy,fz,omega,cnx,cny,cnz,fluid_indicator   
    real(kind=8) :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9
    real(kind=8) :: m_rho,m_e, m_e2, m_jx, m_qx, m_jy, m_qy, m_jz, m_qz, m_3pxx, m_3pixx, m_pww, m_piww, m_pxy, m_pyz, m_pzx, m_tx, m_ty, m_tz
    real(kind=8) :: ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft17,ft18
    real(kind=8) :: g1t0,g1t1,g1t2,g1t3,g1t4,g1t5,g1t6,g1t7,g1t8,g1t9,g1t10,g1t11,g1t12,g1t13,g1t14,g1t15,g1t16,g1t17,g1t18
    real(kind=8) :: g2t0,g2t1,g2t2,g2t3,g2t4,g2t5,g2t6,g2t7,g2t8,g2t9,g2t10,g2t11,g2t12,g2t13,g2t14,g2t15,g2t16,g2t17,g2t18
    real(kind=8) :: s_e,s_e2,s_q,s_nu,s_pi,s_t  !relaxation parameters
    real(kind=8) :: ux1,uy1,uz1,den,u2,rho1,rho2,rho1rho2denbeta,tmp,tmp1

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~openmp
    !$OMP PARALLEL DO default(none) &
    !$OMP & SHARED(&
    !$OMP & walls,f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18, &
    !$OMP & g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,phi,&
    !$OMP & cn_x,cn_y,cn_z,gamma,la_nui1,la_nui2,beta,ixmin,ixmax,iymin,iymax,izmin,izmax,curv,c_norm,force_Z) &
    !$OMP & PRIVATE(&
    !$OMP & i,j,wall_indicator,omega,s_e,s_e2,s_q,s_nu,s_pi,s_t,&
    !$OMP & cnx,cny,cnz,ux1,uy1,uz1,u2,rho1,rho2,den,tmp,tmp1,&
    !$OMP & ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft17,ft18,&
    !$OMP & g1t0,g1t1,g1t2,g1t3,g1t4,g1t5,g1t6,g1t7,g1t8,g1t9,g1t10,g1t11,g1t12,g1t13,g1t14,g1t15,g1t16,g1t17,g1t18,&
    !$OMP & g2t0,g2t1,g2t2,g2t3,g2t4,g2t5,g2t6,g2t7,g2t8,g2t9,g2t10,g2t11,g2t12,g2t13,g2t14,g2t15,g2t16,g2t17,g2t18,&
    !$OMP & sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9,fx,fy,fz, &
    !$OMP & m_rho, m_e, m_e2, m_jx, m_qx, m_jy, m_qy, m_jz, m_qz, m_3pxx, m_3pixx, m_pww, m_piww, m_pxy, m_pyz, m_pzx, m_tx, m_ty, m_tz) collapse(2)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
    !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
    !$acc & g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,phi,walls,c_norm,cn_x,cn_y,cn_z,curv)async(async_label)
    !$acc loop collapse(3) device_type(NVIDIA)
    do k=izmin,izmax
        DO j=iymin,iymax
            !DIR$ IVDEP
            do i=ixmin,ixmax
                !branch treatment, do not perform collision on wall nodes
                wall_indicator = walls(i,j,k)
        
                !+++++++++- AA pattern pull step++++++++++++
                g1t0=  f0(i  ,j  ,k  )
                g1t1=  f1(i-1,j  ,k  )
                g1t2=  f2(i+1,j  ,k  )
                g1t3=  f3(i  ,j-1,k  )
                g1t4=  f4(i  ,j+1,k  )
                g1t5=  f5(i  ,j  ,k-1)
                g1t6=  f6(i  ,j  ,k+1)
                g1t7=  f7(i-1,j-1,k  )
                g1t8=  f8(i+1,j-1,k  )
                g1t9=  f9(i-1,j+1,k  )
                g1t10= f10(i+1,j+1,k  )
                g1t11= f11(i-1,j  ,k-1)
                g1t12= f12(i+1,j  ,k-1)
                g1t13= f13(i-1,j  ,k+1)
                g1t14= f14(i+1,j  ,k+1)
                g1t15= f15(i  ,j-1,k-1)
                g1t16= f16(i  ,j+1,k-1)
                g1t17= f17(i  ,j-1,k+1)
                g1t18= f18(i  ,j+1,k+1)

                g2t0=  g0(i  ,j  ,k  )
                g2t1=  g1(i-1,j  ,k  )
                g2t2=  g2(i+1,j  ,k  )
                g2t3=  g3(i  ,j-1,k  )
                g2t4=  g4(i  ,j+1,k  )
                g2t5=  g5(i  ,j  ,k-1)
                g2t6=  g6(i  ,j  ,k+1)
                g2t7=  g7(i-1,j-1,k  )
                g2t8=  g8(i+1,j-1,k  )
                g2t9=  g9(i-1,j+1,k  )
                g2t10= g10(i+1,j+1,k  )
                g2t11= g11(i-1,j  ,k-1)
                g2t12= g12(i+1,j  ,k-1)
                g2t13= g13(i-1,j  ,k+1)
                g2t14= g14(i+1,j  ,k+1)
                g2t15= g15(i  ,j-1,k-1)
                g2t16= g16(i  ,j+1,k-1)
                g2t17= g17(i  ,j-1,k+1)
                g2t18= g18(i  ,j+1,k+1)

                !let ft be the bulk PDF
                ft0 = g1t0 + g2t0
                ft1  = g1t1 + g2t1
                ft2  = g1t2 + g2t2
                ft3  = g1t3 + g2t3
                ft4  = g1t4 + g2t4
                ft5  = g1t5 + g2t5
                ft6  = g1t6 + g2t6
                ft7  = g1t7 + g2t7
                ft8  = g1t8 + g2t8
                ft9  = g1t9 + g2t9
                ft10  = g1t10 + g2t10
                ft11  = g1t11 + g2t11
                ft12  = g1t12 + g2t12
                ft13  = g1t13 + g2t13
                ft14  = g1t14 + g2t14
                ft15  = g1t15 + g2t15
                ft16  = g1t16 + g2t16
                ft17  = g1t17 + g2t17
                ft18  = g1t18 + g2t18

                !order parameter
                rho1=g1t0+g1t1+g1t2+g1t3+g1t4+g1t5+g1t6+g1t7+g1t8+g1t9+g1t10+g1t11+g1t12+g1t13+g1t14+g1t15+g1t16+g1t17+g1t18
                rho2=g2t0+g2t1+g2t2+g2t3+g2t4+g2t5+g2t6+g2t7+g2t8+g2t9+g2t10+g2t11+g2t12+g2t13+g2t14+g2t15+g2t16+g2t17+g2t18
                phi(i,j,k) = (rho1-rho2)/(rho1+rho2)* (1-wall_indicator) + phi(i,j,k)*wall_indicator

                cnx=cn_x(i,j,k)
                cny=cn_y(i,j,k)
                cnz=cn_z(i,j,k)
                
                tmp = 0.5d0*gamma*curv(i,j,k)*c_norm(i,j,k)
                fx = tmp*cnx
                fy = tmp*cny
                fz = tmp*cnz + force_Z   !body force force_z along flow direction

                
                !++++++++++++-MRT COLLISION+++++++++++++-
                !select viscosity+++++++++++++++++++++++-
                omega = 1d0/(6d0/( (1.0d0 + phi(i,j,k))*la_nui1 + (1.0d0 - phi(i,j,k))*la_nui2 ) + 0.5d0)
                !MRT PARAMETERS
                s_nu =  omega
#if mrt==1
                !************bounceback opt************
                s_e =  omega
                s_e2 = omega
                s_pi = omega
                s_q =  8.0d0*(2.0d0-omega)/(8.0d0-omega)
                s_t = s_q
#elif mrt==2
                !************original************
                s_e =  1.19d0
                s_e2 = 1.4d0
                s_pi = 1.4d0
                s_q= 1.2d0
                s_t = 1.98d0

#elif mrt==3            
                !************SRT************
                s_e =  omega
                s_e2 = omega
                s_pi = omega
                s_q =  omega
                s_t = omega
#elif mrt==4            
                !************advection opt************
                s_e =  omega
                s_e2 = omega
                s_pi = omega
                s_q = (6d0-3d0*omega)/(3d0-omega)
                s_t = omega
#endif 
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MRT kernel, repeated part~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                !+++++++++- !calculate macroscopic variables++++++++++++
                den=   rho1 + rho2
                ux1=  (ft1 - ft2 + ft7 - ft8 + ft9 - ft10 + ft11 - ft12 + ft13 - ft14 + 0.5d0*fx)* (1-wall_indicator)
                uy1=  (ft3 - ft4 + ft7 + ft8 - ft9 - ft10 + ft15 - ft16 + ft17 - ft18 + 0.5d0*fy)* (1-wall_indicator)
                uz1=  (ft5 - ft6 + ft11+ ft12- ft13- ft14 + ft15 + ft16 - ft17 - ft18 + 0.5d0*fz)* (1-wall_indicator)

                u2 = ux1*ux1+uy1*uy1+uz1*uz1

                !PDFs summations for computation efficiency purpose
                sum1 = ft1 + ft2 + ft3 + ft4 + ft5 + ft6
                sum2 = ft7 + ft8 + ft9 + ft10 + ft11 + ft12 + ft13 + ft14 + ft15 + ft16 + ft17 + ft18
                sum3 = ft7 - ft8 + ft9 - ft10 + ft11 - ft12 + ft13 - ft14
                sum4 = ft7 + ft8 - ft9 - ft10 + ft15 - ft16 + ft17 - ft18
                sum5 = ft11+ ft12 -ft13- ft14 + ft15 + ft16 - ft17 - ft18
                sum6 = 2d0*( ft1 + ft2) - ft3 - ft4 - ft5 - ft6
                sum7 = ft7 + ft8 + ft9 + ft10 + ft11 + ft12 + ft13 + ft14 - 2d0*(ft15 + ft16 + ft17 + ft18)
                sum8 = ft3 + ft4 - ft5 - ft6
                sum9 = ft7 + ft8 + ft9 + ft10 - ft11 - ft12 - ft13 - ft14
                !PDF to moment
                m_rho = den
                m_e =  -30d0*ft0 - 11d0*sum1  + 8d0*sum2
                m_e2 =  12d0*ft0 -  4d0*sum1  + sum2
                m_jx = ft1 - ft2 + sum3
                m_qx = -4d0*(ft1 - ft2) + sum3
                m_jy = ft3 - ft4 + sum4
                m_qy = -4d0*(ft3 - ft4) + sum4
                m_jz = ft5 - ft6 + sum5
                m_qz = -4d0*(ft5 - ft6) + sum5
                m_3pxx = sum6 + sum7
                m_3pixx = -2d0*sum6 + sum7
                m_pww = sum8 + sum9
                m_piww =-2d0*sum8 + sum9
                m_pxy = ft7  - ft8  - ft9  + ft10
                m_pyz = ft15 - ft16 - ft17 + ft18
                m_pzx = ft11 - ft12 - ft13 + ft14
                m_tx =   ft7  - ft8  + ft9  - ft10 - ft11 + ft12 - ft13 + ft14
                m_ty = - ft7  - ft8  + ft9  + ft10 + ft15 - ft16 + ft17 - ft18
                m_tz =   ft11 + ft12 - ft13 - ft14 - ft15 - ft16 + ft17 + ft18
         
                !relaxtion in moment space                
                m_e = m_e - s_e*(m_e - (-11.0d0*den+19.0d0*u2)) + (38d0-19d0*s_e)*(fx*ux1+fy*uy1+fz*uz1)                           !m1
                m_e2 = m_e2 - s_e2*(m_e2 - (3.0d0*den - 5.5d0*u2)) + (-11d0+5.5d0*s_e2)*(fx*ux1+fy*uy1+fz*uz1)                     !m2               
                m_jx = m_jx + fx                                                                                                   !m3
                m_qx = m_qx - s_q*(m_qx - (-0.666666666666666667d0*ux1)) + (-0.666666666666666667d0+0.333333333333333333d0*s_q)*fx !m4
                m_jy = m_jy + fy                                                                                                   !m5
                m_qy = m_qy - s_q*(m_qy - (-0.666666666666666667d0*uy1)) + (-0.666666666666666667d0+0.333333333333333333d0*s_q)*fy !m6
                m_jz = m_jz + fz                                                                                                   !m7
                m_qz = m_qz - s_q*(m_qz - (-0.666666666666666667d0*uz1)) + (-0.666666666666666667d0+0.333333333333333333d0*s_q)*fz !m8

                m_3pxx = m_3pxx - s_nu*(m_3pxx - (3d0*ux1*ux1-u2)) + (2d0-s_nu)*(2d0*fx*ux1-fy*uy1-fz*uz1)                            !m9                  
                m_3pixx = m_3pixx - s_pi*(m_3pixx - (-1.5d0*ux1*ux1+0.5d0*u2)) + (1d0-0.5d0*s_pi)*(-2d0*fx*ux1+fy*uy1+fz*uz1)         !m10
                m_pww = m_pww - s_nu*(m_pww - (uy1*uy1-uz1*uz1)) + (2d0-s_nu)*(fy*uy1-fz*uz1)                                      !m11
                m_piww = m_piww - s_pi*(m_piww - (-0.5d0)*(uy1*uy1-uz1*uz1)) + (1d0-0.5d0*s_pi)*(-fy*uy1+fz*uz1)                   !m12         
                m_pxy = m_pxy - s_nu*(m_pxy - (ux1*uy1)) + (1d0-0.5d0*s_nu)*(fx*uy1+fy*ux1)                                        !m13
                m_pyz = m_pyz - s_nu*(m_pyz - (uy1*uz1)) + (1d0-0.5d0*s_nu)*(fy*uz1+fz*uy1)                                        !m14
                m_pzx = m_pzx - s_nu*(m_pzx - (ux1*uz1)) + (1d0-0.5d0*s_nu)*(fx*uz1+fz*ux1)                                        !m15        
                m_tx =  m_tx - s_t*(m_tx )                                                                                         !m16
                m_ty =  m_ty - s_t*(m_ty )                                                                                         !m17
                m_tz =  m_tz - s_t*(m_tz )                                                                                         !m18

                !transform back to PDFs
                !coeffcients for performance
                m_rho = mrt_coef1*m_rho            ! 1/19
                m_e =   mrt_coef2*m_e              ! 1/2394
                m_e2 =  mrt_coef3*m_e2              ! 1/252
                m_jx = 0.1d0*m_jx
                m_qx = 0.025d0*m_qx
                m_jy = 0.1d0*m_jy
                m_qy = 0.025d0*m_qy
                m_jz = 0.1d0*m_jz
                m_qz = 0.025d0*m_qz
                m_3pxx = 2d0*mrt_coef4*m_3pxx      ! 1/36
                m_3pixx = mrt_coef4*m_3pixx        ! 1/72
                m_pww =  6d0*mrt_coef4*m_pww       ! 1/12
                m_piww = 3d0*mrt_coef4*m_piww      ! 1/24
                m_pxy = 0.25d0*m_pxy
                m_pyz = 0.25d0*m_pyz
                m_pzx = 0.25d0*m_pzx
                m_tx = 0.125d0* m_tx
                m_ty = 0.125d0* m_ty
                m_tz = 0.125d0* m_tz
                sum1 = m_rho - 11d0*m_e - 4d0*m_e2
                sum2 = 2d0*m_3pxx - 4d0*m_3pixx
                sum3 = m_pww - 2d0*m_piww
                sum4 = m_rho + 8d0*m_e + m_e2
                sum5 = m_jx + m_qx
                sum6 = m_jy + m_qy
                sum7 = m_jz + m_qz
                sum8 = m_3pxx + m_3pixx
                sum9 = m_pww + m_piww

                ft0 =  m_rho -  30d0*m_e + 12d0*m_e2
                ft1 =  sum1 +  m_jx - 4d0*m_qx + sum2
                ft2 =  sum1 -  m_jx + 4d0*m_qx + sum2
                ft3 =  sum1 +  m_jy - 4d0*m_qy - 0.5d0*sum2 + sum3
                ft4 =  sum1 -  m_jy + 4d0*m_qy - 0.5d0*sum2 + sum3
                ft5 =  sum1 +  m_jz - 4d0*m_qz - 0.5d0*sum2 - sum3
                ft6 =  sum1 -  m_jz + 4d0*m_qz - 0.5d0*sum2 - sum3
                ft7 =  sum4 +  sum5 + sum6 + sum8 + sum9  + m_pxy + m_tx - m_ty
                ft8 =  sum4 -  sum5 + sum6 + sum8 + sum9  - m_pxy - m_tx - m_ty
                ft9 =  sum4 +  sum5 - sum6 + sum8 + sum9  - m_pxy + m_tx + m_ty
                ft10 =  sum4 -  sum5 - sum6 + sum8 + sum9  + m_pxy - m_tx + m_ty
                ft11 =  sum4 +  sum5 + sum7 + sum8 - sum9  + m_pzx - m_tx + m_tz
                ft12 =  sum4 -  sum5 + sum7 + sum8 - sum9  - m_pzx + m_tx + m_tz
                ft13 =  sum4 +  sum5 - sum7 + sum8 - sum9  - m_pzx - m_tx - m_tz
                ft14 =  sum4 -  sum5 - sum7 + sum8 - sum9  + m_pzx + m_tx - m_tz
                ft15 =  sum4 +  sum6 + sum7 - sum8*2d0     + m_pyz + m_ty - m_tz
                ft16 =  sum4 -  sum6 + sum7 - sum8*2d0     - m_pyz - m_ty - m_tz
                ft17 =  sum4 +  sum6 - sum7 - sum8*2d0     - m_pyz + m_ty + m_tz
                ft18 =  sum4 -  sum6 - sum7 - sum8*2d0     + m_pyz - m_ty + m_tz
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MRT kernel, repeated part~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                !++++++++++++-recoloring & streaming to opposite direction++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                tmp1 = rho1/den
                f0(i,j,k) = tmp1*ft0
                g0(i,j,k) = ft0*(1d0-tmp1)

                !R-K method
                tmp =  rho1*rho2*beta/den

                g1t1  =  (tmp1*ft1  + w_equ_1*tmp * ( cnx)) * (1-wall_indicator) +  f2(i+1,j  ,k  )*wall_indicator
                g1t2  =  (tmp1*ft2  + w_equ_1*tmp * (-cnx)) * (1-wall_indicator) +  f1(i-1,j  ,k  )*wall_indicator
                g1t3  =  (tmp1*ft3  + w_equ_1*tmp * ( cny)) * (1-wall_indicator) +  f4(i  ,j+1,k  )*wall_indicator
                g1t4  =  (tmp1*ft4  + w_equ_1*tmp * (-cny)) * (1-wall_indicator) +  f3(i  ,j-1,k  )*wall_indicator
                g1t5  =  (tmp1*ft5  + w_equ_1*tmp * ( cnz)) * (1-wall_indicator) +  f6(i  ,j  ,k+1)*wall_indicator
                g1t6  =  (tmp1*ft6  + w_equ_1*tmp * (-cnz)) * (1-wall_indicator) +  f5(i  ,j  ,k-1)*wall_indicator
                g1t7 =  (tmp1*ft7  + rk_weight2 * tmp * ( cnx + cny))* (1-wall_indicator) +f10(i+1,j+1,k  )*wall_indicator
                g1t8 =  (tmp1*ft8  + rk_weight2 * tmp * (-cnx + cny))* (1-wall_indicator) +  f9(i-1,j+1,k  )*wall_indicator
                g1t9 =  (tmp1*ft9  + rk_weight2 * tmp * ( cnx - cny))* (1-wall_indicator) +  f8(i+1,j-1,k  )*wall_indicator
                g1t10  =  (tmp1*ft10 + rk_weight2 * tmp * (-cnx - cny))* (1-wall_indicator) +  f7(i-1,j-1,k  )*wall_indicator
                g1t11 =  (tmp1*ft11 + rk_weight2 * tmp * ( cnx + cnz))* (1-wall_indicator) +  f14(i+1,j  ,k+1)*wall_indicator
                g1t12 =  (tmp1*ft12 + rk_weight2 * tmp * (-cnx + cnz))* (1-wall_indicator) +  f13(i-1,j  ,k+1)*wall_indicator
                g1t13 =  (tmp1*ft13 + rk_weight2 * tmp * ( cnx - cnz))* (1-wall_indicator) +  f12(i+1,j  ,k-1)*wall_indicator
                g1t14 =  (tmp1*ft14 + rk_weight2 * tmp * (-cnx - cnz))* (1-wall_indicator) +  f11(i-1,j  ,k-1)*wall_indicator
                g1t15 =  (tmp1*ft15 + rk_weight2 * tmp * ( cny + cnz))* (1-wall_indicator) +  f18(i  ,j+1,k+1)*wall_indicator
                g1t16 =  (tmp1*ft16 + rk_weight2 * tmp * (-cny + cnz))* (1-wall_indicator) +  f17(i  ,j-1,k+1)*wall_indicator
                g1t17 =  (tmp1*ft17 + rk_weight2 * tmp * ( cny - cnz))* (1-wall_indicator) +  f16(i  ,j+1,k-1)*wall_indicator
                g1t18 =  (tmp1*ft18 + rk_weight2 * tmp * (-cny - cnz))* (1-wall_indicator) +  f15(i  ,j-1,k-1)*wall_indicator

                g2t1  =  (ft1  - g1t1) * (1-wall_indicator) +  g2(i+1,j  ,k  )*wall_indicator
                g2t2  =  (ft2  - g1t2) * (1-wall_indicator) +  g1(i-1,j  ,k  )*wall_indicator
                g2t3  =  (ft3  - g1t3) * (1-wall_indicator) +  g4(i  ,j+1,k  )*wall_indicator
                g2t4  =  (ft4  - g1t4) * (1-wall_indicator) +  g3(i  ,j-1,k  )*wall_indicator
                g2t5  =  (ft5  - g1t5) * (1-wall_indicator) +  g6(i  ,j  ,k+1)*wall_indicator
                g2t6  =  (ft6  - g1t6) * (1-wall_indicator) +  g5(i  ,j  ,k-1)*wall_indicator
                g2t7 =  (ft7  - g1t7)* (1-wall_indicator) +    g10(i+1,j+1,k  )*wall_indicator
                g2t8 =  (ft8  - g1t8)* (1-wall_indicator) +    g9(i-1,j+1,k  )*wall_indicator
                g2t9 =  (ft9  - g1t9)* (1-wall_indicator) +    g8(i+1,j-1,k  )*wall_indicator
                g2t10  =  (ft10 - g1t10)* (1-wall_indicator) + g7(i-1,j-1,k  )*wall_indicator
                g2t11 =  (ft11 - g1t11)* (1-wall_indicator) + g14(i+1,j  ,k+1)*wall_indicator
                g2t12 =  (ft12 - g1t12)* (1-wall_indicator) + g13(i-1,j  ,k+1)*wall_indicator
                g2t13 =  (ft13 - g1t13)* (1-wall_indicator) + g12(i+1,j  ,k-1)*wall_indicator
                g2t14 =  (ft14 - g1t14)* (1-wall_indicator) + g11(i-1,j  ,k-1)*wall_indicator
                g2t15 =  (ft15 - g1t15)* (1-wall_indicator) + g18(i  ,j+1,k+1)*wall_indicator
                g2t16 =  (ft16 - g1t16)* (1-wall_indicator) + g17(i  ,j-1,k+1)*wall_indicator
                g2t17 =  (ft17 - g1t17)* (1-wall_indicator) + g16(i  ,j+1,k-1)*wall_indicator
                g2t18 =  (ft18 - g1t18)* (1-wall_indicator) + g15(i  ,j-1,k-1)*wall_indicator

                !+++++++++- AA pattern push step++++++++++++
                f2(i+1,j  ,k  ) = g1t1
                f1(i-1,j  ,k  ) = g1t2
                f4(i  ,j+1,k  ) = g1t3
                f3(i  ,j-1,k  ) = g1t4
                f6(i  ,j  ,k+1) = g1t5
                f5(i  ,j  ,k-1) = g1t6
                f10(i+1,j+1,k  ) = g1t7
                f9(i-1,j+1,k  ) = g1t8
                f8(i+1,j-1,k  ) = g1t9
                f7(i-1,j-1,k  ) = g1t10
                f14(i+1,j  ,k+1) = g1t11
                f13(i-1,j  ,k+1) = g1t12
                f12(i+1,j  ,k-1) = g1t13
                f11(i-1,j  ,k-1) = g1t14
                f18(i  ,j+1,k+1) = g1t15
                f17(i  ,j-1,k+1) = g1t16
                f16(i  ,j+1,k-1) = g1t17
                f15(i  ,j-1,k-1) = g1t18

                g2(i+1,j  ,k  ) = g2t1
                g1(i-1,j  ,k  ) = g2t2
                g4(i  ,j+1,k  ) = g2t3
                g3(i  ,j-1,k  ) = g2t4
                g6(i  ,j  ,k+1) = g2t5
                g5(i  ,j  ,k-1) = g2t6
                g10(i+1,j+1,k  ) = g2t7
                g9(i-1,j+1,k  ) = g2t8
                g8(i+1,j-1,k  ) = g2t9
                g7(i-1,j-1,k  ) = g2t10
                g14(i+1,j  ,k+1) = g2t11
                g13(i-1,j  ,k+1) = g2t12
                g12(i+1,j  ,k-1) = g2t13
                g11(i-1,j  ,k-1) = g2t14
                g18(i  ,j+1,k+1) = g2t15
                g17(i  ,j-1,k+1) = g2t16
                g16(i  ,j+1,k-1) = g2t17
                g15(i  ,j-1,k-1) = g2t18
            ENDDO
        ENDDO
    enddo
    !$acc end kernels

    RETURN
endsubroutine kernel_odd_color




!=====================================================================================================================================
!---------------------- even step kernel ----------------------
! no steaming steps
!=====================================================================================================================================
subroutine kernel_even_color(ixmin,ixmax,iymin,iymax,izmin,izmax,async_label)
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    IMPLICIT NONE
    integer :: i,j,k,m, ixmin,ixmax,iymin,iymax,izmin,izmax,async_label
    integer(kind=1) :: wall_indicator
    real(kind=8) :: f_force,cnx,cny,cnz,rho1,rho2,tmp1   !color model
    real(kind=8) :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9
    real(kind=8) :: m_rho,m_e, m_e2, m_jx, m_qx, m_jy, m_qy, m_jz, m_qz, m_3pxx, m_3pixx, m_pww, m_piww, m_pxy, m_pyz, m_pzx, m_tx, m_ty, m_tz
    real(kind=8) :: ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft17,ft18
    real(kind=8) :: g1t0,g1t1,g1t2,g1t3,g1t4,g1t5,g1t6,g1t7,g1t8,g1t9,g1t10,g1t11,g1t12,g1t13,g1t14,g1t15,g1t16,g1t17,g1t18
    real(kind=8) :: g2t0,g2t1,g2t2,g2t3,g2t4,g2t5,g2t6,g2t7,g2t8,g2t9,g2t10,g2t11,g2t12,g2t13,g2t14,g2t15,g2t16,g2t17,g2t18
    real(kind=8) :: s_e,s_e2,s_q,s_nu,s_pi,s_t  !relaxation parameters
    real(kind=8) :: fx,fy,fz, ux1,uy1,uz1,den,tmp,omega,u2

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~openmp
    !$OMP PARALLEL DO default(none) &
    !$OMP & SHARED(&
    !$OMP & walls,f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18, &
    !$OMP & g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,phi,&
    !$OMP & cn_x,cn_y,cn_z,gamma,la_nui1,la_nui2,beta,ixmin,ixmax,iymin,iymax,izmin,izmax,curv,c_norm,force_Z) &
    !$OMP & PRIVATE(&
    !$OMP & i,j,wall_indicator,omega,s_e,s_e2,s_q,s_nu,s_pi,s_t,&
    !$OMP & cnx,cny,cnz,ux1,uy1,uz1,u2,rho1,rho2,den,tmp1,tmp,&
    !$OMP & ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft17,ft18,&
    !$OMP & g1t0,g1t1,g1t2,g1t3,g1t4,g1t5,g1t6,g1t7,g1t8,g1t9,g1t10,g1t11,g1t12,g1t13,g1t14,g1t15,g1t16,g1t17,g1t18,&
    !$OMP & g2t0,g2t1,g2t2,g2t3,g2t4,g2t5,g2t6,g2t7,g2t8,g2t9,g2t10,g2t11,g2t12,g2t13,g2t14,g2t15,g2t16,g2t17,g2t18,&
    !$OMP & sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9,fx,fy,fz,&
    !$OMP & m_rho, m_e, m_e2, m_jx, m_qx, m_jy, m_qy, m_jz, m_qz, m_3pxx, m_3pixx, m_pww, m_piww, m_pxy, m_pyz, m_pzx, m_tx, m_ty, m_tz) collapse(2)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
    !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
    !$acc & g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,phi,walls,c_norm,cn_x,cn_y,cn_z,curv)async(async_label)
    !$acc loop collapse(3) device_type(NVIDIA)
    do k=izmin,izmax
        DO j=iymin,iymax
            !DIR$ IVDEP
            do i=ixmin,ixmax
                !branch treatment, do not perform collision on wall nodes
                wall_indicator = walls(i,j,k)
                !+++++++++- AA pattern pull step++++++++++++
                g1t0 = f0(i,j,k)
                g1t1  = f2(i,j,k)
                g1t2  = f1(i,j,k)
                g1t3  = f4(i,j,k)
                g1t4  = f3(i,j,k)
                g1t5  = f6(i,j,k)
                g1t6  = f5(i,j,k)
                g1t7  = f10(i,j,k)
                g1t8  = f9(i,j,k)
                g1t9  = f8(i,j,k)
                g1t10  = f7(i,j,k)
                g1t11  = f14(i,j,k)
                g1t12  = f13(i,j,k)
                g1t13  = f12(i,j,k)
                g1t14  = f11(i,j,k)
                g1t15  = f18(i,j,k)
                g1t16  = f17(i,j,k)
                g1t17  = f16(i,j,k)
                g1t18  = f15(i,j,k)

                g2t0 = g0(i,j,k)
                g2t1  = g2(i,j,k)
                g2t2  = g1(i,j,k)
                g2t3  = g4(i,j,k)
                g2t4  = g3(i,j,k)
                g2t5  = g6(i,j,k)
                g2t6  = g5(i,j,k)
                g2t7  = g10(i,j,k)
                g2t8  = g9(i,j,k)
                g2t9  = g8(i,j,k)
                g2t10  = g7(i,j,k)
                g2t11  = g14(i,j,k)
                g2t12  = g13(i,j,k)
                g2t13  = g12(i,j,k)
                g2t14  = g11(i,j,k)
                g2t15  = g18(i,j,k)
                g2t16  = g17(i,j,k)
                g2t17  = g16(i,j,k)
                g2t18  = g15(i,j,k)

                !ft: bulk PDF
                ft0 = g1t0 + g2t0
                ft1  = g1t1 + g2t1
                ft2  = g1t2 + g2t2
                ft3  = g1t3 + g2t3
                ft4  = g1t4 + g2t4
                ft5  = g1t5 + g2t5
                ft6  = g1t6 + g2t6
                ft7  = g1t7 + g2t7
                ft8  = g1t8 + g2t8
                ft9  = g1t9 + g2t9
                ft10  = g1t10 + g2t10
                ft11  = g1t11 + g2t11
                ft12  = g1t12 + g2t12
                ft13  = g1t13 + g2t13
                ft14  = g1t14 + g2t14
                ft15  = g1t15 + g2t15
                ft16  = g1t16 + g2t16
                ft17  = g1t17 + g2t17
                ft18  = g1t18 + g2t18

                !order parameter
                rho1=g1t0+g1t1+g1t2+g1t3+g1t4+g1t5+g1t6+g1t7+g1t8+g1t9+g1t10+g1t11+g1t12+g1t13+g1t14+g1t15+g1t16+g1t17+g1t18
                rho2=g2t0+g2t1+g2t2+g2t3+g2t4+g2t5+g2t6+g2t7+g2t8+g2t9+g2t10+g2t11+g2t12+g2t13+g2t14+g2t15+g2t16+g2t17+g2t18
                phi(i,j,k) = (rho1-rho2)/(rho1+rho2)* (1-wall_indicator) + phi(i,j,k)*wall_indicator

                cnx=cn_x(i,j,k)
                cny=cn_y(i,j,k)
                cnz=cn_z(i,j,k)

                tmp = 0.5d0*gamma*curv(i,j,k)*c_norm(i,j,k)
                fx = tmp*cnx
                fy = tmp*cny
                fz = tmp*cnz + force_Z   !body force force_z along flow direction

                !++++++++++++-MRT COLLISION+++++++++++++-
                !select viscosity+++++++++++++++++++++++-
                omega = 1d0/(6d0/( (1.0d0 + phi(i,j,k))*la_nui1 + (1.0d0 - phi(i,j,k))*la_nui2 ) + 0.5d0)
                !MRT PARAMETERS
                s_nu =  omega
#if mrt==1
                !************bounceback opt************
                s_e =  omega
                s_e2 = omega
                s_pi = omega
                s_q =  8.0d0*(2.0d0-omega)/(8.0d0-omega)
                s_t = s_q
#elif mrt==2
                !************original opt************
                s_e =  1.19d0
                s_e2 = 1.4d0
                s_pi = 1.4d0
                s_q= 1.2d0
                s_t = 1.98d0               
#elif mrt==3
                !************SRT************
                s_e =  omega
                s_e2 = omega
                s_pi = omega
                s_q =  omega
                s_t = omega
#elif mrt==4            
                !************advection opt************
                s_e =  omega
                s_e2 = omega
                s_pi = omega
                s_q = (6d0-3d0*omega)/(3d0-omega)
                s_t = omega
#endif

                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MRT kernel, repeated part~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                !+++++++++- !calculate macroscopic variables++++++++++++
                den=   rho1 + rho2

                ux1=  (ft1 - ft2 + ft7 - ft8 + ft9 - ft10 + ft11 - ft12 + ft13 - ft14 + 0.5d0*fx)* (1-wall_indicator)
                uy1=  (ft3 - ft4 + ft7 + ft8 - ft9 - ft10 + ft15 - ft16 + ft17 - ft18 + 0.5d0*fy)* (1-wall_indicator)
                uz1=  (ft5 - ft6 + ft11+ ft12- ft13- ft14 + ft15 + ft16 - ft17 - ft18 + 0.5d0*fz)* (1-wall_indicator)
                u2 = ux1*ux1+uy1*uy1+uz1*uz1

                !PDFs summations for computation efficiency purpose
                sum1 = ft1 + ft2 + ft3 + ft4 + ft5 + ft6
                sum2 = ft7 + ft8 + ft9 + ft10 + ft11 + ft12 + ft13 + ft14 + ft15 + ft16 + ft17 + ft18
                sum3 = ft7 - ft8 + ft9 - ft10 + ft11 - ft12 + ft13 - ft14
                sum4 = ft7 + ft8 - ft9 - ft10 + ft15 - ft16 + ft17 - ft18
                sum5 = ft11+ ft12 -ft13- ft14 + ft15 + ft16 - ft17 - ft18
                sum6 = 2d0*( ft1 + ft2) - ft3 - ft4 - ft5 - ft6
                sum7 = ft7 + ft8 + ft9 + ft10 + ft11 + ft12 + ft13 + ft14 - 2d0*(ft15 + ft16 + ft17 + ft18)
                sum8 = ft3 + ft4 - ft5 - ft6
                sum9 = ft7 + ft8 + ft9 + ft10 - ft11 - ft12 - ft13 - ft14
                !PDF to moment
                m_rho = den
                m_e =  -30d0*ft0 - 11d0*sum1  + 8d0*sum2
                m_e2 =  12d0*ft0 -  4d0*sum1  + sum2
                m_jx = ft1 - ft2 + sum3
                m_qx = -4d0*(ft1 - ft2) + sum3
                m_jy = ft3 - ft4 + sum4
                m_qy = -4d0*(ft3 - ft4) + sum4
                m_jz = ft5 - ft6 + sum5
                m_qz = -4d0*(ft5 - ft6) + sum5
                m_3pxx = sum6 + sum7
                m_3pixx = -2d0*sum6 + sum7
                m_pww = sum8 + sum9
                m_piww =-2d0*sum8 + sum9
                m_pxy = ft7  - ft8  - ft9  + ft10
                m_pyz = ft15 - ft16 - ft17 + ft18
                m_pzx = ft11 - ft12 - ft13 + ft14
                m_tx =   ft7  - ft8  + ft9  - ft10 - ft11 + ft12 - ft13 + ft14
                m_ty = - ft7  - ft8  + ft9  + ft10 + ft15 - ft16 + ft17 - ft18
                m_tz =   ft11 + ft12 - ft13 - ft14 - ft15 - ft16 + ft17 + ft18

                !relaxtion in moment space                
                m_e = m_e - s_e*(m_e - (-11.0d0*den+19.0d0*u2)) + (38d0-19d0*s_e)*(fx*ux1+fy*uy1+fz*uz1)                           !m1
                m_e2 = m_e2 - s_e2*(m_e2 - (3.0d0*den - 5.5d0*u2)) + (-11d0+5.5d0*s_e2)*(fx*ux1+fy*uy1+fz*uz1)                     !m2               
                m_jx = m_jx + fx                                                                                                   !m3
                m_qx = m_qx - s_q*(m_qx - (-0.666666666666666667d0*ux1)) + (-0.666666666666666667d0+0.333333333333333333d0*s_q)*fx !m4
                m_jy = m_jy + fy                                                                                                   !m5
                m_qy = m_qy - s_q*(m_qy - (-0.666666666666666667d0*uy1)) + (-0.666666666666666667d0+0.333333333333333333d0*s_q)*fy !m6
                m_jz = m_jz + fz                                                                                                   !m7
                m_qz = m_qz - s_q*(m_qz - (-0.666666666666666667d0*uz1)) + (-0.666666666666666667d0+0.333333333333333333d0*s_q)*fz !m8

                m_3pxx = m_3pxx - s_nu*(m_3pxx - (3d0*ux1*ux1-u2)) + (2d0-s_nu)*(2d0*fx*ux1-fy*uy1-fz*uz1)                         !m9                  
                m_3pixx = m_3pixx - s_pi*(m_3pixx - (-1.5d0*ux1*ux1+0.5d0*u2)) + (1d0-0.5d0*s_pi)*(-2d0*fx*ux1+fy*uy1+fz*uz1)      !m10
                m_pww = m_pww - s_nu*(m_pww - (uy1*uy1-uz1*uz1)) + (2d0-s_nu)*(fy*uy1-fz*uz1)                                      !m11
                m_piww = m_piww - s_pi*(m_piww - (-0.5d0)*(uy1*uy1-uz1*uz1)) + (1d0-0.5d0*s_pi)*(-fy*uy1+fz*uz1)                   !m12         
                m_pxy = m_pxy - s_nu*(m_pxy - (ux1*uy1)) + (1d0-0.5d0*s_nu)*(fx*uy1+fy*ux1)                                        !m13
                m_pyz = m_pyz - s_nu*(m_pyz - (uy1*uz1)) + (1d0-0.5d0*s_nu)*(fy*uz1+fz*uy1)                                        !m14
                m_pzx = m_pzx - s_nu*(m_pzx - (ux1*uz1)) + (1d0-0.5d0*s_nu)*(fx*uz1+fz*ux1)                                        !m15        
                m_tx =  m_tx - s_t*(m_tx )                                                                                         !m16
                m_ty =  m_ty - s_t*(m_ty )                                                                                         !m17
                m_tz =  m_tz - s_t*(m_tz )                                                                                         !m18


                !transform back to PDFs
                !coeffcients for performance
                m_rho = mrt_coef1*m_rho            ! 1/19
                m_e =   mrt_coef2*m_e              ! 1/2394
                m_e2 =  mrt_coef3*m_e2              ! 1/252
                m_jx = 0.1d0*m_jx
                m_qx = 0.025d0*m_qx
                m_jy = 0.1d0*m_jy
                m_qy = 0.025d0*m_qy
                m_jz = 0.1d0*m_jz
                m_qz = 0.025d0*m_qz
                m_3pxx = 2d0*mrt_coef4*m_3pxx      ! 1/36
                m_3pixx = mrt_coef4*m_3pixx        ! 1/72
                m_pww =  6d0*mrt_coef4*m_pww       ! 1/12
                m_piww = 3d0*mrt_coef4*m_piww      ! 1/24
                m_pxy = 0.25d0*m_pxy
                m_pyz = 0.25d0*m_pyz
                m_pzx = 0.25d0*m_pzx
                m_tx = 0.125d0* m_tx
                m_ty = 0.125d0* m_ty
                m_tz = 0.125d0* m_tz
                sum1 = m_rho - 11d0*m_e - 4d0*m_e2
                sum2 = 2d0*m_3pxx - 4d0*m_3pixx
                sum3 = m_pww - 2d0*m_piww
                sum4 = m_rho + 8d0*m_e + m_e2
                sum5 = m_jx + m_qx
                sum6 = m_jy + m_qy
                sum7 = m_jz + m_qz
                sum8 = m_3pxx + m_3pixx
                sum9 = m_pww + m_piww

                ft0 =  m_rho -  30d0*m_e + 12d0*m_e2
                ft1 =  sum1 +  m_jx - 4d0*m_qx + sum2
                ft2 =  sum1 -  m_jx + 4d0*m_qx + sum2
                ft3 =  sum1 +  m_jy - 4d0*m_qy - 0.5d0*sum2 + sum3
                ft4 =  sum1 -  m_jy + 4d0*m_qy - 0.5d0*sum2 + sum3
                ft5 =  sum1 +  m_jz - 4d0*m_qz - 0.5d0*sum2 - sum3
                ft6 =  sum1 -  m_jz + 4d0*m_qz - 0.5d0*sum2 - sum3
                ft7 =  sum4 +  sum5 + sum6 + sum8 + sum9  + m_pxy + m_tx - m_ty
                ft8 =  sum4 -  sum5 + sum6 + sum8 + sum9  - m_pxy - m_tx - m_ty
                ft9 =  sum4 +  sum5 - sum6 + sum8 + sum9  - m_pxy + m_tx + m_ty
                ft10 =  sum4 -  sum5 - sum6 + sum8 + sum9  + m_pxy - m_tx + m_ty
                ft11 =  sum4 +  sum5 + sum7 + sum8 - sum9  + m_pzx - m_tx + m_tz
                ft12 =  sum4 -  sum5 + sum7 + sum8 - sum9  - m_pzx + m_tx + m_tz
                ft13 =  sum4 +  sum5 - sum7 + sum8 - sum9  - m_pzx - m_tx - m_tz
                ft14 =  sum4 -  sum5 - sum7 + sum8 - sum9  + m_pzx + m_tx - m_tz
                ft15 =  sum4 +  sum6 + sum7 - sum8*2d0     + m_pyz + m_ty - m_tz
                ft16 =  sum4 -  sum6 + sum7 - sum8*2d0     - m_pyz - m_ty - m_tz
                ft17 =  sum4 +  sum6 - sum7 - sum8*2d0     - m_pyz + m_ty + m_tz
                ft18 =  sum4 -  sum6 - sum7 - sum8*2d0     + m_pyz - m_ty + m_tz
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MRT kernel, repeated part~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                !++++++++++++-recoloring++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                tmp1 = rho1/den
                f0(i,j,k) = tmp1*ft0
                g0(i,j,k) = ft0*(1d0-tmp1)

                !R-K method
                tmp = rho1*rho2*beta/den

                g1t1 = (tmp1*ft1 + w_equ_1*tmp * ( cnx)) * (1-wall_indicator) + f1(i,j,k)*wall_indicator
                g1t2 = (tmp1*ft2 + w_equ_1*tmp * (-cnx)) * (1-wall_indicator) + f2(i,j,k)*wall_indicator
                g1t3 = (tmp1*ft3 + w_equ_1*tmp * ( cny)) * (1-wall_indicator) + f3(i,j,k)*wall_indicator
                g1t4 = (tmp1*ft4 + w_equ_1*tmp * (-cny)) * (1-wall_indicator) + f4(i,j,k)*wall_indicator
                g1t5 = (tmp1*ft5 + w_equ_1*tmp * ( cnz)) * (1-wall_indicator) + f5(i,j,k)*wall_indicator
                g1t6 = (tmp1*ft6 + w_equ_1*tmp * (-cnz)) * (1-wall_indicator) + f6(i,j,k)*wall_indicator
                g1t7  = (tmp1*ft7  + rk_weight2 * tmp * ( cnx + cny)) * (1-wall_indicator) + f7(i,j,k)*wall_indicator
                g1t8  = (tmp1*ft8  + rk_weight2 * tmp * (-cnx + cny)) * (1-wall_indicator) + f8(i,j,k)*wall_indicator
                g1t9  = (tmp1*ft9  + rk_weight2 * tmp * ( cnx - cny)) * (1-wall_indicator) + f9(i,j,k)*wall_indicator
                g1t10 = (tmp1*ft10 + rk_weight2 * tmp * (-cnx - cny)) * (1-wall_indicator) + f10(i,j,k)*wall_indicator
                g1t11 = (tmp1*ft11 + rk_weight2 * tmp * ( cnx + cnz)) * (1-wall_indicator) + f11(i,j,k)*wall_indicator
                g1t12 = (tmp1*ft12 + rk_weight2 * tmp * (-cnx + cnz)) * (1-wall_indicator) + f12(i,j,k)*wall_indicator
                g1t13 = (tmp1*ft13 + rk_weight2 * tmp * ( cnx - cnz)) * (1-wall_indicator) + f13(i,j,k)*wall_indicator
                g1t14 = (tmp1*ft14 + rk_weight2 * tmp * (-cnx - cnz)) * (1-wall_indicator) + f14(i,j,k)*wall_indicator
                g1t15 = (tmp1*ft15 + rk_weight2 * tmp * ( cny + cnz)) * (1-wall_indicator) + f15(i,j,k)*wall_indicator
                g1t16 = (tmp1*ft16 + rk_weight2 * tmp * (-cny + cnz)) * (1-wall_indicator) + f16(i,j,k)*wall_indicator
                g1t17 = (tmp1*ft17 + rk_weight2 * tmp * ( cny - cnz)) * (1-wall_indicator) + f17(i,j,k)*wall_indicator
                g1t18 = (tmp1*ft18 + rk_weight2 * tmp * (-cny - cnz)) * (1-wall_indicator) + f18(i,j,k)*wall_indicator

                g2t1= ( ft1 - g1t1) * (1-wall_indicator) +g1(i,j,k)*wall_indicator
                g2t2= ( ft2 - g1t2) * (1-wall_indicator) +g2(i,j,k)*wall_indicator
                g2t3= ( ft3 - g1t3) * (1-wall_indicator) +g3(i,j,k)*wall_indicator
                g2t4= ( ft4 - g1t4) * (1-wall_indicator) +g4(i,j,k)*wall_indicator
                g2t5= ( ft5 - g1t5) * (1-wall_indicator) +g5(i,j,k)*wall_indicator
                g2t6= ( ft6 - g1t6) * (1-wall_indicator) +g6(i,j,k)*wall_indicator
                g2t7= ( ft7 - g1t7) * (1-wall_indicator) +g7(i,j,k)*wall_indicator
                g2t8= ( ft8 - g1t8) * (1-wall_indicator) +g8(i,j,k)*wall_indicator
                g2t9= ( ft9 - g1t9) * (1-wall_indicator) +g9(i,j,k)*wall_indicator
                g2t10= ( ft10 - g1t10) * (1-wall_indicator) +g10(i,j,k)*wall_indicator
                g2t11= ( ft11 - g1t11) * (1-wall_indicator) +g11(i,j,k)*wall_indicator
                g2t12= ( ft12 - g1t12) * (1-wall_indicator) +g12(i,j,k)*wall_indicator
                g2t13= ( ft13 - g1t13) * (1-wall_indicator) +g13(i,j,k)*wall_indicator
                g2t14= ( ft14 - g1t14) * (1-wall_indicator) +g14(i,j,k)*wall_indicator
                g2t15= ( ft15 - g1t15) * (1-wall_indicator) +g15(i,j,k)*wall_indicator
                g2t16= ( ft16 - g1t16) * (1-wall_indicator) +g16(i,j,k)*wall_indicator
                g2t17= ( ft17 - g1t17) * (1-wall_indicator) +g17(i,j,k)*wall_indicator
                g2t18= ( ft18 - g1t18) * (1-wall_indicator) +g18(i,j,k)*wall_indicator

                !+++++++++- AA pattern++++++++++++
                f1(i,j,k)  = g1t1
                f2(i,j,k)  = g1t2
                f3(i,j,k)  = g1t3
                f4(i,j,k)  = g1t4
                f5(i,j,k)  = g1t5
                f6(i,j,k)  = g1t6
                f7(i,j,k)  = g1t7
                f8(i,j,k)  = g1t8
                f9(i,j,k)  = g1t9
                f10(i,j,k)  = g1t10
                f11(i,j,k)  = g1t11
                f12(i,j,k)  = g1t12
                f13(i,j,k)  = g1t13
                f14(i,j,k)  = g1t14
                f15(i,j,k)  = g1t15
                f16(i,j,k)  = g1t16
                f17(i,j,k)  = g1t17
                f18(i,j,k)  = g1t18

                g1(i,j,k)  = g2t1
                g2(i,j,k)  = g2t2
                g3(i,j,k)  = g2t3
                g4(i,j,k)  = g2t4
                g5(i,j,k)  = g2t5
                g6(i,j,k)  = g2t6
                g7(i,j,k)  = g2t7
                g8(i,j,k)  = g2t8
                g9(i,j,k)  = g2t9
                g10(i,j,k)  = g2t10
                g11(i,j,k)  = g2t11
                g12(i,j,k)  = g2t12
                g13(i,j,k)  = g2t13
                g14(i,j,k)  = g2t14
                g15(i,j,k)  = g2t15
                g16(i,j,k)  = g2t16
                g17(i,j,k)  = g2t17
                g18(i,j,k)  = g2t18

            ENDDO
        ENDDO
    enddo
    !$acc end kernels

    RETURN
end subroutine kernel_even_color



