#include "./preprocessor.h"
!=======================================================================================================================================================================================================================
!---------------------- multiphase outlet convective BC ----------------------
! convective bc applies to unknown PDFs only
!=======================================================================================================================================================================================================================
!************************** before odd step kernel *****************************************
subroutine outlet_convective_BC_before_odd    !before streaming type BC
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    real(kind=8) :: temp,u_convec,ft5,ft6,ft11,ft12,ft13,ft14,ft15,ft16,ft17,ft18
    integer(kind=1) :: wall_indicator

    if(idz==npz-1)then
        !inlet average velocity
        u_convec = uin_avg
        temp = 1d0/(1d0+u_convec)

        !$omp parallel do private(wall_indicator,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18, &
        !$acc & g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,phi,walls,u,v,w,f_convec_bc,g_convec_bc,phi_convec_bc)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,nz)
                
                phi(i,j,nz+1) = ((phi_convec_bc(i,j) + u_convec*phi(i,j,nz))*temp) * (1-wall_indicator) + phi(i,j,nz+1)*wall_indicator
                phi_convec_bc(i,j) = phi(i,j,nz+1)   !store PDF for next step
                phi(i,j,nz+2) = phi(i,j,nz+1)
                phi(i,j,nz+3) = phi(i,j,nz+1)
                phi(i,j,nz+4) = phi(i,j,nz+1)   !overlap_phi=4

                !if outlet convective BC
                f6(i    , j    , nz + 1) = ((f_convec_bc(i,j,6) + u_convec * f6(i    , j    , nz + 1 - 1)) * temp) * (1 - wall_indicator) + f6(i    , j    , nz + 1)*wall_indicator
                f13(i - 1, j    , nz + 1) = ((f_convec_bc(i,j,13) + u_convec * f13(i - 1, j    , nz + 1 - 1)) * temp) * (1 - wall_indicator) + f13(i - 1, j    , nz + 1)*wall_indicator
                f14(i + 1, j    , nz + 1) = ((f_convec_bc(i,j,14) + u_convec * f14(i + 1, j    , nz + 1 - 1)) * temp) * (1 - wall_indicator) + f14(i + 1, j    , nz + 1)*wall_indicator
                f17(i    , j - 1, nz + 1) = ((f_convec_bc(i,j,17) + u_convec * f17(i    , j - 1, nz + 1 - 1)) * temp) * (1 - wall_indicator) + f17(i    , j - 1, nz + 1)*wall_indicator
                f18(i    , j + 1, nz + 1) = ((f_convec_bc(i,j,18) + u_convec * f18(i    , j + 1, nz + 1 - 1)) * temp) * (1 - wall_indicator) + f18(i    , j + 1, nz + 1)*wall_indicator

                g6(i    , j    , nz + 1) = ((g_convec_bc(i,j,6) + u_convec * g6(i    , j    , nz + 1 - 1)) * temp) * (1 - wall_indicator) + g6(i    , j    , nz + 1)*wall_indicator
                g13(i - 1, j    , nz + 1) = ((g_convec_bc(i,j,13) + u_convec * g13(i - 1, j    , nz + 1 - 1)) * temp) * (1 - wall_indicator) + g13(i - 1, j    , nz + 1)*wall_indicator
                g14(i + 1, j    , nz + 1) = ((g_convec_bc(i,j,14) + u_convec * g14(i + 1, j    , nz + 1 - 1)) * temp) * (1 - wall_indicator) + g14(i + 1, j    , nz + 1)*wall_indicator
                g17(i    , j - 1, nz + 1) = ((g_convec_bc(i,j,17) + u_convec * g17(i    , j - 1, nz + 1 - 1)) * temp) * (1 - wall_indicator) + g17(i    , j - 1, nz + 1)*wall_indicator
                g18(i    , j + 1, nz + 1) = ((g_convec_bc(i,j,18) + u_convec * g18(i    , j + 1, nz + 1 - 1)) * temp) * (1 - wall_indicator) + g18(i    , j + 1, nz + 1)*wall_indicator

                f_convec_bc(i,j,6) = f6(i    , j    , nz + 1)            
                f_convec_bc(i,j,13) = f13(i - 1, j    , nz + 1)
                f_convec_bc(i,j,14) = f14(i + 1, j    , nz + 1)             
                f_convec_bc(i,j,17) = f17(i    , j - 1, nz + 1)
                f_convec_bc(i,j,18) = f18(i    , j + 1, nz + 1)

                g_convec_bc(i,j,6) = g6(i    , j    , nz + 1)            
                g_convec_bc(i,j,13) = g13(i - 1, j    , nz + 1)
                g_convec_bc(i,j,14) = g14(i + 1, j    , nz + 1)             
                g_convec_bc(i,j,17) = g17(i    , j - 1, nz + 1)
                g_convec_bc(i,j,18) = g18(i    , j + 1, nz + 1)
            enddo
        enddo
        !$acc end kernels
    endif
end subroutine outlet_convective_BC_before_odd

!************************** after odd step kernel *****************************************
subroutine outlet_convective_BC_after_odd   !after streaming type BC
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    real(kind=8) :: temp,u_convec,ft5,ft6,ft11,ft12,ft13,ft14,ft15,ft16,ft17,ft18
    integer(kind=1) :: wall_indicator

    if(idz==npz-1)then
        !inlet average velocity
        u_convec = uin_avg
        temp = 1d0/(1d0+u_convec)

        !$omp parallel do private(wall_indicator,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc & g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,phi,walls,u,v,w,f_convec_bc,g_convec_bc,phi_convec_bc)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,nz)

                phi(i,j,nz+1) = ((phi_convec_bc(i,j) + u_convec*phi(i,j,nz))*temp) * (1-wall_indicator) + phi(i,j,nz+1)*wall_indicator
                phi_convec_bc(i,j) = phi(i,j,nz+1)   !store PDF for next step
                phi(i,j,nz+2) = phi(i,j,nz+1)
                phi(i,j,nz+3) = phi(i,j,nz+1)
                phi(i,j,nz+4) = phi(i,j,nz+1)   !overlap_phi=4
                !outlet convective bc
              
                f5(i,j,nz) = (( f_convec_bc(i,j,6) + u_convec*f5(i,j,nz-1))*temp  )* (1 - wall_indicator)   + f5(i,j,nz)*wall_indicator              
                f11(i,j,nz) = (( f_convec_bc(i,j,14) + u_convec*f11(i,j,nz-1))*temp  )* (1 - wall_indicator)   + f11(i,j,nz)*wall_indicator
                f12(i,j,nz) = (( f_convec_bc(i,j,13) + u_convec*f12(i,j,nz-1))*temp  )* (1 - wall_indicator)   + f12(i,j,nz)*wall_indicator             
                f15(i,j,nz) = (( f_convec_bc(i,j,18) + u_convec*f15(i,j,nz-1))*temp  )* (1 - wall_indicator)   + f15(i,j,nz)*wall_indicator
                f16(i,j,nz) = (( f_convec_bc(i,j,17) + u_convec*f16(i,j,nz-1))*temp  )* (1 - wall_indicator)   + f16(i,j,nz)*wall_indicator
                
                g5(i,j,nz) = (( g_convec_bc(i,j,6) + u_convec*g5(i,j,nz-1))*temp  )* (1 - wall_indicator)   + g5(i,j,nz)*wall_indicator            
                g11(i,j,nz) = (( g_convec_bc(i,j,14) + u_convec*g11(i,j,nz-1))*temp  )* (1 - wall_indicator)   + g11(i,j,nz)*wall_indicator
                g12(i,j,nz) = (( g_convec_bc(i,j,13) + u_convec*g12(i,j,nz-1))*temp  )* (1 - wall_indicator)   + g12(i,j,nz)*wall_indicator                
                g15(i,j,nz) = (( g_convec_bc(i,j,18) + u_convec*g15(i,j,nz-1))*temp  )* (1 - wall_indicator)   + g15(i,j,nz)*wall_indicator
                g16(i,j,nz) = (( g_convec_bc(i,j,17) + u_convec*g16(i,j,nz-1))*temp  )* (1 - wall_indicator)   + g16(i,j,nz)*wall_indicator
        
                f_convec_bc(i,j,6) = f5(i,j,nz)                        
                f_convec_bc(i,j,14) = f11(i,j,nz)  
                f_convec_bc(i,j,13) = f12(i,j,nz)               
                f_convec_bc(i,j,18) = f15(i,j,nz)
                f_convec_bc(i,j,17) = f16(i,j,nz)
        
                g_convec_bc(i,j,6) = g5(i,j,nz)                        
                g_convec_bc(i,j,14) = g11(i,j,nz)
                g_convec_bc(i,j,13) = g12(i,j,nz)               
                g_convec_bc(i,j,18) = g15(i,j,nz)
                g_convec_bc(i,j,17) = g16(i,j,nz)
            enddo
        enddo
        !$acc end kernels
    endif

end subroutine outlet_convective_BC_after_odd







!=======================================================================================================================================================================================================================
!---------------------- Zou-He type pressure open outlet boundary conditions ----------------------
!=======================================================================================================================================================================================================================
!************************** before odd step kernel *****************************************
subroutine outlet_Zou_He_pressure_BC_before_odd    !before streaming type BC
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    integer(kind=1) :: wall_indicator
    real(kind=8) :: tmp1,tmp2,tnx,tny,rho1,rho2,ux1,uy1,uz1

    if(idz==npz-1)then
        !$omp parallel
        !$omp do private(wall_indicator,tmp1,tmp2,tnx,tny,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc &g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,phi,walls)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,nz)
                phi(i,j,nz+1) = phi(i,j,nz)
                phi(i,j,nz+2) = phi(i,j,nz)
                phi(i,j,nz+3) = phi(i,j,nz)
                phi(i,j,nz+4) = phi(i,j,nz)   !overlap_phi=4

                !outlet pressure BC    k=1a
                tmp1 = (&
                    f0(i,j,nz)+&
                    f1(i-1,j,nz)+&
                    f2(i+1,j,nz)+&
                    f3(i,j-1,nz)+&
                    f4(i,j+1,nz)+& 
                    f7(i-1,j-1,nz)+&
                    f8(i+1,j-1,nz)+&
                    f9(i-1,j+1,nz)+&
                    f10(i+1,j+1,nz)+ 2d0*(&
                    f5(i,j,nz-1)+&
                    f11(i-1,j,nz-1)+&
                    f12(i+1,j,nz-1)+&
                    f15(i,j-1,nz-1)+&
                    f16(i,j+1,nz-1))+&
                    g0(i,j,nz)+&
                    g1(i-1,j,nz)+&
                    g2(i+1,j,nz)+&
                    g3(i,j-1,nz)+&
                    g4(i,j+1,nz)+&
                    g7(i-1,j-1,nz)+&
                    g8(i+1,j-1,nz)+&
                    g9(i-1,j+1,nz)+&
                    g10(i+1,j+1,nz)+ 2d0*(&
                    g5(i,j,nz-1)+&
                    g11(i-1,j,nz-1)+&
                    g12(i+1,j,nz-1)+&
                    g15(i,j-1,nz-1)+&
                    g16(i,j+1,nz-1))) - rho_out

                tmp2 = tmp1*0.5d0*(1d0-phi(i,j,nz))    !fluid 2 net flux
                tmp1 = tmp1 - tmp2                   !fluid 1 net flux

                tnx = 0.5d0*(&
                    f1(i-1,j,nz)+f7(i-1,j-1,nz)+f9(i-1,j+1,nz)-(&
                    f2(i+1,j,nz)+f8(i+1,j-1,nz)+f10(i+1,j+1,nz)))
                tny = 0.5d0*(&
                    f3(i,j-1,nz)+f7(i-1,j-1,nz)+f8(i+1,j-1,nz)-(&
                    f4(i,j+1,nz)+f10(i+1,j+1,nz)+f9(i-1,j+1,nz)))
                f6(i,j,nz+1)     = (f5(i,j,nz-1)-0.333333333333333333d0*tmp1) * (1d0-wall_indicator)             + f6(i,j,nz+1)*wall_indicator
                f13(i-1,j,nz+1) = (f12(i+1,j,nz-1)-0.166666666666666667d0*tmp1 - tnx) * (1d0-wall_indicator) + f13(i-1,j,nz+1)*wall_indicator
                f14(i+1,j,nz+1) = (f11(i-1,j,nz-1)-0.166666666666666667d0*tmp1 + tnx) * (1d0-wall_indicator) + f14(i+1,j,nz+1)*wall_indicator
                f17(i,j-1,nz+1) = (f16(i,j+1,nz-1)-0.166666666666666667d0*tmp1 - tny) * (1d0-wall_indicator) + f17(i,j-1,nz+1)*wall_indicator
                f18(i,j+1,nz+1) = (f15(i,j-1,nz-1)-0.166666666666666667d0*tmp1 + tny) * (1d0-wall_indicator) + f18(i,j+1,nz+1)*wall_indicator

                tnx = 0.5d0*(&
                    g1(i-1,j,nz)+g7(i-1,j-1,nz)+g9(i-1,j+1,nz)-(&
                    g2(i+1,j,nz)+g8(i+1,j-1,nz)+g10(i+1,j+1,nz)))
                tny = 0.5d0*(&
                    g3(i,j-1,nz)+g7(i-1,j-1,nz)+g8(i+1,j-1,nz)-(&
                    g4(i,j+1,nz)+g10(i+1,j+1,nz)+g9(i-1,j+1,nz)))
                g6(i,j,nz+1)     = (g5(i,j,nz-1)-0.333333333333333333d0*tmp2) * (1d0-wall_indicator)             + g6(i,j,nz+1)*wall_indicator
                g13(i-1,j,nz+1) = (g12(i+1,j,nz-1)-0.166666666666666667d0*tmp2 - tnx) * (1d0-wall_indicator) + g13(i-1,j,nz+1)*wall_indicator
                g14(i+1,j,nz+1) = (g11(i-1,j,nz-1)-0.166666666666666667d0*tmp2 + tnx) * (1d0-wall_indicator) + g14(i+1,j,nz+1)*wall_indicator
                g17(i,j-1,nz+1) = (g16(i,j+1,nz-1)-0.166666666666666667d0*tmp2 - tny) * (1d0-wall_indicator) + g17(i,j-1,nz+1)*wall_indicator
                g18(i,j+1,nz+1) = (g15(i,j-1,nz-1)-0.166666666666666667d0*tmp2 + tny) * (1d0-wall_indicator) + g18(i,j+1,nz+1)*wall_indicator
            enddo
        enddo
        !$acc end kernels
        !$omp end parallel        
    endif

    return
end subroutine outlet_Zou_He_pressure_BC_before_odd

!************************** after odd step kernel *****************************************
subroutine outlet_Zou_He_pressure_BC_after_odd     !after streaming type BC
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    integer(kind=1) :: wall_indicator
    real(kind=8) :: tmp1,tmp2,tnx,tny

    if(idz==npz-1)then
        !$omp parallel
        !$omp do private(wall_indicator,tmp1,tmp2,tnx,tny,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc &g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,phi,walls)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,nz)
                phi(i,j,nz+1) = phi(i,j,nz)
                phi(i,j,nz+2) = phi(i,j,nz)
                phi(i,j,nz+3) = phi(i,j,nz)
                phi(i,j,nz+4) = phi(i,j,nz)   !overlap_phi=4
                
                !inlet pressure BC    k=1
                tmp1 = (&
                    f0(i,j,nz)+&
                    f2(i,j,nz)+&
                    f1(i,j,nz)+&
                    f4(i,j,nz)+&
                    f3(i,j,nz)+&
                    f8(i,j,nz)+&
                    f7(i,j,nz)+&
                    f10(i,j,nz)+&
                    f9(i,j,nz)+2d0*(&
                    f6(i,j,nz)+&
                    f14(i,j,nz)+&
                    f13(i,j,nz)+&
                    f18(i,j,nz)+&
                    f17(i,j,nz))+&
                    g0(i,j,nz)+&
                    g2(i,j,nz)+&
                    g1(i,j,nz)+&
                    g4(i,j,nz)+&
                    g3(i,j,nz)+&
                    g8(i,j,nz)+&
                    g7(i,j,nz)+&
                    g10(i,j,nz)+&
                    g9(i,j,nz)+2d0*(&
                    g6(i,j,nz)+&
                    g14(i,j,nz)+&
                    g13(i,j,nz)+&
                    g18(i,j,nz)+&
                    g17(i,j,nz)) ) - rho_out

                tmp2 = tmp1*0.5d0*(1d0-phi(i,j,nz))    !fluid 2 net flux
                tmp1 = tmp1 - tmp2                   !fluid 1 net flux

                tnx = 0.5d0*( f2(i,j,nz)+f8(i,j,nz)+f10(i,j,nz)-(f1(i,j,nz)+f7(i,j,nz)+f9(i,j,nz)) )
                tny = 0.5d0*( f4(i,j,nz)+f10(i,j,nz)+f9(i,j,nz)-(f3(i,j,nz)+f7(i,j,nz)+f8(i,j,nz)) )
                f5(i,j,nz) = (f6(i,j,nz)-0.333333333333333333d0*tmp1)* (1d0-wall_indicator)   + f5(i,j,nz)*wall_indicator
                f11(i,j,nz) = (f14(i,j,nz)-0.166666666666666667d0*tmp1 + tnx)* (1d0-wall_indicator)   + f11(i,j,nz)*wall_indicator
                f12(i,j,nz) = (f13(i,j,nz)-0.166666666666666667d0*tmp1 - tnx)* (1d0-wall_indicator)   + f12(i,j,nz)*wall_indicator
                f15(i,j,nz) = (f18(i,j,nz)-0.166666666666666667d0*tmp1 + tny)* (1d0-wall_indicator)   + f15(i,j,nz)*wall_indicator
                f16(i,j,nz) = (f17(i,j,nz)-0.166666666666666667d0*tmp1 - tny)* (1d0-wall_indicator)   + f16(i,j,nz)*wall_indicator

                tnx = 0.5d0*( g2(i,j,nz)+g8(i,j,nz)+g10(i,j,nz)-(g1(i,j,nz)+g7(i,j,nz)+g9(i,j,nz)) )
                tny = 0.5d0*( g4(i,j,nz)+g10(i,j,nz)+g9(i,j,nz)-(g3(i,j,nz)+g7(i,j,nz)+g8(i,j,nz)) )
                g5(i,j,nz) = (g6(i,j,nz)-0.333333333333333333d0*tmp2)* (1d0-wall_indicator)   + g5(i,j,nz)*wall_indicator
                g11(i,j,nz) = (g14(i,j,nz)-0.166666666666666667d0*tmp2 + tnx)* (1d0-wall_indicator)   + g11(i,j,nz)*wall_indicator
                g12(i,j,nz) = (g13(i,j,nz)-0.166666666666666667d0*tmp2 - tnx)* (1d0-wall_indicator)   + g12(i,j,nz)*wall_indicator
                g15(i,j,nz) = (g18(i,j,nz)-0.166666666666666667d0*tmp2 + tny)* (1d0-wall_indicator)   + g15(i,j,nz)*wall_indicator
                g16(i,j,nz) = (g17(i,j,nz)-0.166666666666666667d0*tmp2 - tny)* (1d0-wall_indicator)   + g16(i,j,nz)*wall_indicator
            enddo
        enddo
        !$acc end kernels
        !$omp end parallel
    endif
    return
end subroutine outlet_Zou_He_pressure_BC_after_odd




