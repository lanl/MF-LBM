#include "./preprocessor.h"
!============================================================================================================================================================================
!---------------------- Bounce-back type velocity open inlet boundary conditions (default) ----------------------
!============================================================================================================================================================================
!************************** before odd step kernel *****************************************
subroutine inlet_bounce_back_velocity_BC_before_odd    !before streaming type BC
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    real(kind=8) :: tmp1,tmp2
    integer(kind=1) :: wall_indicator

    if(idz==0)then
        !$omp parallel do private(wall_indicator,tmp1,tmp2,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc &g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,walls,w_in)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,1)

                phi(i,j,0) = phi_inlet*(1-wall_indicator)+  phi(i,j,0)*wall_indicator
                phi(i,j,-1) =  phi(i,j,0)
                phi(i,j,-2) =  phi(i,j,0)
                phi(i,j,-3) =  phi(i,j,0)   !overlap_phi = 4
                
                !inlet velocity BC    k=1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                tmp2 = w_in(i,j)*relaxation
                tmp1 = tmp2*sa_inject              !fluid 1 injection
                tmp2 = tmp2-tmp1                   !fluid 2 injection

                f5(i,j,0) = (f6(i, j, 1)+ 6.0d0*w_equ_1*tmp1 ) * (1-wall_indicator) + f5(i,j,0)*wall_indicator
                f11(i-1,j,0) = (f14(i,j,1)+ 6.0d0*w_equ_2*tmp1 ) * (1-wall_indicator) + f11(i-1,j,0)*wall_indicator
                f12(i+1,j,0) = (f13(i,j,1)+ 6.0d0*w_equ_2*tmp1 ) * (1-wall_indicator) + f12(i+1,j,0)*wall_indicator
                f15(i,j-1,0) = (f18(i,j,1)+ 6.0d0*w_equ_2*tmp1 ) * (1-wall_indicator) + f15(i,j-1,0)*wall_indicator
                f16(i,j+1,0) = (f17(i,j,1)+ 6.0d0*w_equ_2*tmp1 ) * (1-wall_indicator) + f16(i,j+1,0)*wall_indicator

                g5(i,j,0) = (g6(i, j, 1)+ 6.0d0*w_equ_1*tmp2 ) * (1-wall_indicator) + g5(i,j,0)*wall_indicator
                g11(i-1,j,0) = (g14(i,j,1)+ 6.0d0*w_equ_2*tmp2 ) * (1-wall_indicator) + g11(i-1,j,0)*wall_indicator
                g12(i+1,j,0) = (g13(i,j,1)+ 6.0d0*w_equ_2*tmp2 ) * (1-wall_indicator) + g12(i+1,j,0)*wall_indicator
                g15(i,j-1,0) = (g18(i,j,1)+ 6.0d0*w_equ_2*tmp2 ) * (1-wall_indicator) + g15(i,j-1,0)*wall_indicator
                g16(i,j+1,0) = (g17(i,j,1)+ 6.0d0*w_equ_2*tmp2 ) * (1-wall_indicator) + g16(i,j+1,0)*wall_indicator
            enddo
        enddo
        !$acc end kernels
    endif

    return
end subroutine inlet_bounce_back_velocity_BC_before_odd

!************************** after odd step kernel *****************************************
subroutine inlet_bounce_back_velocity_BC_after_odd   !after streaming type BC
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    real(kind=8) :: tmp1,tmp2
    integer(kind=1) :: wall_indicator

    if(idz==0)then
        !$omp parallel do private(wall_indicator,tmp1,tmp2,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc &g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,walls,w_in)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,1)

                phi(i,j,0) = phi_inlet*(1-wall_indicator)+  phi(i,j,0)*wall_indicator
                phi(i,j,-1) =  phi(i,j,0)
                phi(i,j,-2) =  phi(i,j,0)
                phi(i,j,-3) =  phi(i,j,0)   !overlap_phi = 4
                
                !inlet velocity BC    k=1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                tmp2 = w_in(i,j)*relaxation
                tmp1 = tmp2*sa_inject              !fluid 1 injection
                tmp2 = tmp2-tmp1                   !fluid 2 injection

                f6(i,j,1) = (f5(i,j,0)+ 6.0d0*w_equ_1*tmp1 )* (1-wall_indicator)   + f6(i,j,1)*wall_indicator
                f13(i,j,1) = (f12(i+1,j,0)+ 6.0d0*w_equ_2*tmp1 )* (1-wall_indicator)   + f13(i,j,1)*wall_indicator
                f14(i,j,1) = (f11(i-1,j,0)+ 6.0d0*w_equ_2*tmp1 )* (1-wall_indicator)   + f14(i,j,1)*wall_indicator
                f17(i,j,1) = (f16(i,j+1,0)+ 6.0d0*w_equ_2*tmp1 )* (1-wall_indicator)   + f17(i,j,1)*wall_indicator
                f18(i,j,1) = (f15(i,j-1,0)+ 6.0d0*w_equ_2*tmp1 )* (1-wall_indicator)   + f18(i,j,1)*wall_indicator

                g6(i,j,1) = (g5(i,j,0)+ 6.0d0*w_equ_1*tmp2 )* (1-wall_indicator)   + g6(i,j,1)*wall_indicator
                g13(i,j,1) = (g12(i+1,j,0)+ 6.0d0*w_equ_2*tmp2 )* (1-wall_indicator)   + g13(i,j,1)*wall_indicator
                g14(i,j,1) = (g11(i-1,j,0)+ 6.0d0*w_equ_2*tmp2 )* (1-wall_indicator)   + g14(i,j,1)*wall_indicator
                g17(i,j,1) = (g16(i,j+1,0)+ 6.0d0*w_equ_2*tmp2 )* (1-wall_indicator)   + g17(i,j,1)*wall_indicator
                g18(i,j,1) = (g15(i,j-1,0)+ 6.0d0*w_equ_2*tmp2 )* (1-wall_indicator)   + g18(i,j,1)*wall_indicator
            enddo
        enddo
        !$acc end kernels
    endif

end subroutine inlet_bounce_back_velocity_BC_after_odd






!==================================================================================================================================================================
!---------------------- Zou-He type pressure/velocity open inlet boundary conditions ----------------------
! currently, there should be only one dominant phase at the inlet boundary nodes, 
! otherwise recolor scheme conflicts with zou-he BC for individual fluid component (momentumn for individual fluid is not conserved)
!==================================================================================================================================================================
!************************** before odd step kernel *****************************************
subroutine inlet_Zou_He_pressure_BC_before_odd    !before streaming type BC
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    integer(kind=1) :: wall_indicator
    real(kind=8) :: tmp1,tmp2,tnx,tny,ux1,uy1,uz1

    if(idz==0)then
        !$omp parallel do private(wall_indicator,tmp1,tmp2,tnx,tny,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc &g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,walls,w_in)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,1)

                phi(i,j,0) = phi_inlet*(1-wall_indicator)+  phi(i,j,0)*wall_indicator
                phi(i,j,-1) =  phi(i,j,0)
                phi(i,j,-2) =  phi(i,j,0)
                phi(i,j,-3) =  phi(i,j,0)   !overlap_phi = 4

                !inlet pressure BC    k=1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                !Zou-He pressure BC applied to the bulk PDF
                ! bulk fluid injection
                tmp1 = (rho_in - &
                    (f0(i,j,1)+&
                    f1(i-1,j,1)+&
                    f2(i+1,j,1)+&
                    f3(i,j-1,1)+&
                    f4(i,j+1,1)+&
                    f7(i-1,j-1,1)+&
                    f8(i+1,j-1,1)+&
                    f9(i-1,j+1,1)+&
                    f10(i+1,j+1,1)+ 2d0*(&
                    f6(i,j,2)+&
                    f14(i+1,j,2)+&
                    f13(i-1,j,2)+&
                    f18(i,j+1,2)+&
                    f17(i,j-1,2)) )  )*relaxation

                tnx = 0.5d0*(&
                    f1(i-1,j,1)+f7(i-1,j-1,1)+f9(i-1,j+1,1)-(&
                    f2(i+1,j,1)+f8(i+1,j-1,1)+f10(i+1,j+1,1)))
                tny = 0.5d0*(&
                    f3(i,j-1,1)+f7(i-1,j-1,1)+f8(i+1,j-1,1)-(&
                    f4(i,j+1,1)+f10(i+1,j+1,1)+f9(i-1,j+1,1)))
                f5(i,j,0)     = (f6(i, j, 2)+0.333333333333333333d0*tmp1) * (1-wall_indicator)        + f5(i,j,0)*wall_indicator
                f11(i-1,j,0) = (f14(i+1,j,2)+0.166666666666666667d0*tmp1 - tnx) * (1-wall_indicator) + f11(i-1,j,0)*wall_indicator
                f12(i+1,j,0) = (f13(i-1,j,2)+0.166666666666666667d0*tmp1 + tnx) * (1-wall_indicator) + f12(i+1,j,0)*wall_indicator
                f15(i,j-1,0) = (f18(i,j+1,2)+0.166666666666666667d0*tmp1 - tny) * (1-wall_indicator) + f15(i,j-1,0)*wall_indicator
                f16(i,j+1,0) = (f17(i,j-1,2)+0.166666666666666667d0*tmp1 + tny) * (1-wall_indicator) + f16(i,j+1,0)*wall_indicator

                !bounce-back for the other phase   
                g5(i,j,0) = g6(i,j,1) 
                g11(i-1,j,0) = g14(i,j,1) 
                g12(i+1,j,0) = g13(i,j,1) 
                g15(i,j-1,0) = g18(i,j,1) 
                g16(i,j+1,0) = g17(i,j,1)                                
            enddo
        enddo
        !$acc end kernels

    endif

    return
end subroutine inlet_Zou_He_pressure_BC_before_odd

!************************** after odd step kernel *****************************************
subroutine inlet_Zou_He_pressure_BC_after_odd   !after streaming type BC
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    integer(kind=1) :: wall_indicator
    real(kind=8) :: tmp1,tmp2,tnx,tny,ux1,uy1,uz1

    if(idz==0)then
        !$omp parallel do private(wall_indicator,tmp1,tmp2,tnx,tny,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc &g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,walls,w_in)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,1)

                phi(i,j,0) = phi_inlet*(1-wall_indicator)+  phi(i,j,0)*wall_indicator
                phi(i,j,-1) =  phi(i,j,0)
                phi(i,j,-2) =  phi(i,j,0)
                phi(i,j,-3) =  phi(i,j,0)   !overlap_phi = 4

                !inlet pressure BC    k=1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                !Zou-He pressure BC applied to the bulk PDF
                ! bulk fluid injection
                tmp1 = (rho_in - (&
                    f0(i,j,1)+&
                    f2(i,j,1)+&
                    f1(i,j,1)+&
                    f4(i,j,1)+&
                    f3(i,j,1)+&
                    f8(i,j,1)+&
                    f7(i,j,1)+&
                    f10(i,j,1)+&
                    f9(i,j,1)+2d0*(&
                    f5(i,j,1)+&
                    f11(i,j,1)+&
                    f12(i,j,1)+&
                    f15(i,j,1)+&
                    f16(i,j,1))) )*relaxation

                tnx = 0.5d0*( f2(i,j,1)+f8(i,j,1)+f10(i,j,1)-(f1(i,j,1)+f7(i,j,1)+f9(i,j,1)) )
                tny = 0.5d0*( f4(i,j,1)+f9(i,j,1)+f10(i,j,1)-(f3(i,j,1)+f8(i,j,1)+f7(i,j,1)) )
                f6(i,j,1) = (f5(i,j,1)+0.333333333333333333d0*tmp1)* (1-wall_indicator)   + f6(i,j,1)*wall_indicator
                f13(i,j,1) = (f12(i,j,1)+0.166666666666666667d0*tmp1 + tnx)* (1-wall_indicator)   + f13(i,j,1)*wall_indicator
                f14(i,j,1) = (f11(i,j,1)+0.166666666666666667d0*tmp1 - tnx)* (1-wall_indicator)   + f14(i,j,1)*wall_indicator
                f17(i,j,1) = (f16(i,j,1)+0.166666666666666667d0*tmp1 + tny)* (1-wall_indicator)   + f17(i,j,1)*wall_indicator
                f18(i,j,1) = (f15(i,j,1)+0.166666666666666667d0*tmp1 - tny)* (1-wall_indicator)   + f18(i,j,1)*wall_indicator

                !bounce-back for the other phase   
                g6(i,j,1) = g5(i,j,0) 
                g13(i,j,1) = g12(i+1,j,0) 
                g14(i,j,1) = g11(i-1,j,0)    
                g17(i,j,1) = g16(i,j+1,0)    
                g18(i,j,1) = g15(i,j-1,0)                                       
            enddo
        enddo
        !$acc end kernels

    endif
end subroutine inlet_Zou_He_pressure_BC_after_odd



