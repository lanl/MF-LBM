#include "./preprocessor.h"
!============================================================================================================================================================================
!---------------------- Bounce-back type velocity open inlet boundary conditions (default) ----------------------
!============================================================================================================================================================================
!************************** before odd step kernel *****************************************
subroutine inlet_bounce_back_velocity_BC_before_odd    !before streaming type BC
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    real(kind=8) :: tmp1,tmp2
    integer(kind=1) :: wall_indicator

    if(idz==0)then
        !$omp parallel do private(wall_indicator,tmp1,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,walls,w_in)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,1)
                
                !inlet velocity BC    k=1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                tmp1 = w_in(i,j)*relaxation

                f5(i,j,0) = (f6(i, j, 1)+ 6.0d0*w_equ_1*tmp1 ) * (1-wall_indicator) + f5(i,j,0)*wall_indicator
                f11(i-1,j,0) = (f14(i,j,1)+ 6.0d0*w_equ_2*tmp1 ) * (1-wall_indicator) + f11(i-1,j,0)*wall_indicator
                f12(i+1,j,0) = (f13(i,j,1)+ 6.0d0*w_equ_2*tmp1 ) * (1-wall_indicator) + f12(i+1,j,0)*wall_indicator
                f15(i,j-1,0) = (f18(i,j,1)+ 6.0d0*w_equ_2*tmp1 ) * (1-wall_indicator) + f15(i,j-1,0)*wall_indicator
                f16(i,j+1,0) = (f17(i,j,1)+ 6.0d0*w_equ_2*tmp1 ) * (1-wall_indicator) + f16(i,j+1,0)*wall_indicator

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
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    real(kind=8) :: tmp1
    integer(kind=1) :: wall_indicator

    if(idz==0)then
        !$omp parallel do private(wall_indicator,tmp1,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,walls,w_in)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,1)
                
                !inlet velocity BC    k=1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                tmp1 = w_in(i,j)*relaxation

                f6(i,j,1) = (f5(i,j,0)+ 6.0d0*w_equ_1*tmp1 )* (1-wall_indicator)   + f6(i,j,1)*wall_indicator
                f13(i,j,1) = (f12(i+1,j,0)+ 6.0d0*w_equ_2*tmp1 )* (1-wall_indicator)   + f13(i,j,1)*wall_indicator
                f14(i,j,1) = (f11(i-1,j,0)+ 6.0d0*w_equ_2*tmp1 )* (1-wall_indicator)   + f14(i,j,1)*wall_indicator
                f17(i,j,1) = (f16(i,j+1,0)+ 6.0d0*w_equ_2*tmp1 )* (1-wall_indicator)   + f17(i,j,1)*wall_indicator
                f18(i,j,1) = (f15(i,j-1,0)+ 6.0d0*w_equ_2*tmp1 )* (1-wall_indicator)   + f18(i,j,1)*wall_indicator

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
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    integer(kind=1) :: wall_indicator
    real(kind=8) :: tmp1,tmp2,tnx,tny,ux1,uy1,uz1

    if(idz==0)then
        !$omp parallel do private(wall_indicator,tmp1,tnx,tny,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,walls,w_in)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,1)

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
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    integer(kind=1) :: wall_indicator
    real(kind=8) :: tmp1,tmp2,tnx,tny,ux1,uy1,uz1

    if(idz==0)then
        !$omp parallel do private(wall_indicator,tmp1,tnx,tny,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,walls,w_in)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,1)

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
            enddo
        enddo
        !$acc end kernels

    endif
end subroutine inlet_Zou_He_pressure_BC_after_odd







!=======================================================================================================================================================================================================================
!---------------------- outlet convective BC ----------------------
! convective bc applies to unknown PDFs only
!=======================================================================================================================================================================================================================
!************************** before odd step kernel *****************************************
subroutine outlet_convective_BC_before_odd    !before streaming type BC
    use Misc_module
    use Fluid_singlephase
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
        !$acc & walls,u,v,w,f_convec_bc)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,nz)
            
                !if outlet convective BC
                f6(i    , j    , nz + 1) = ((f_convec_bc(i,j,6) + u_convec * f6(i    , j    , nz + 1 - 1)) * temp) * (1 - wall_indicator) + f6(i    , j    , nz + 1)*wall_indicator
                f13(i - 1, j    , nz + 1) = ((f_convec_bc(i,j,13) + u_convec * f13(i - 1, j    , nz + 1 - 1)) * temp) * (1 - wall_indicator) + f13(i - 1, j    , nz + 1)*wall_indicator
                f14(i + 1, j    , nz + 1) = ((f_convec_bc(i,j,14) + u_convec * f14(i + 1, j    , nz + 1 - 1)) * temp) * (1 - wall_indicator) + f14(i + 1, j    , nz + 1)*wall_indicator
                f17(i    , j - 1, nz + 1) = ((f_convec_bc(i,j,17) + u_convec * f17(i    , j - 1, nz + 1 - 1)) * temp) * (1 - wall_indicator) + f17(i    , j - 1, nz + 1)*wall_indicator
                f18(i    , j + 1, nz + 1) = ((f_convec_bc(i,j,18) + u_convec * f18(i    , j + 1, nz + 1 - 1)) * temp) * (1 - wall_indicator) + f18(i    , j + 1, nz + 1)*wall_indicator

                f_convec_bc(i,j,6) = f6(i    , j    , nz + 1)            
                f_convec_bc(i,j,13) = f13(i - 1, j    , nz + 1)
                f_convec_bc(i,j,14) = f14(i + 1, j    , nz + 1)             
                f_convec_bc(i,j,17) = f17(i    , j - 1, nz + 1)
                f_convec_bc(i,j,18) = f18(i    , j + 1, nz + 1)
            enddo
        enddo
        !$acc end kernels
    endif
end subroutine outlet_convective_BC_before_odd

!************************** after odd step kernel *****************************************
subroutine outlet_convective_BC_after_odd   !after streaming type BC
    use Misc_module
    use Fluid_singlephase
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
        !$acc & walls,u,v,w,f_convec_bc)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,nz)
              
                f5(i,j,nz) = (( f_convec_bc(i,j,6) + u_convec*f5(i,j,nz-1))*temp  )* (1 - wall_indicator)   + f5(i,j,nz)*wall_indicator              
                f11(i,j,nz) = (( f_convec_bc(i,j,14) + u_convec*f11(i,j,nz-1))*temp  )* (1 - wall_indicator)   + f11(i,j,nz)*wall_indicator
                f12(i,j,nz) = (( f_convec_bc(i,j,13) + u_convec*f12(i,j,nz-1))*temp  )* (1 - wall_indicator)   + f12(i,j,nz)*wall_indicator             
                f15(i,j,nz) = (( f_convec_bc(i,j,18) + u_convec*f15(i,j,nz-1))*temp  )* (1 - wall_indicator)   + f15(i,j,nz)*wall_indicator
                f16(i,j,nz) = (( f_convec_bc(i,j,17) + u_convec*f16(i,j,nz-1))*temp  )* (1 - wall_indicator)   + f16(i,j,nz)*wall_indicator
        
                f_convec_bc(i,j,6) = f5(i,j,nz)                        
                f_convec_bc(i,j,14) = f11(i,j,nz)  
                f_convec_bc(i,j,13) = f12(i,j,nz)               
                f_convec_bc(i,j,18) = f15(i,j,nz)
                f_convec_bc(i,j,17) = f16(i,j,nz)
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
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    integer(kind=1) :: wall_indicator
    real(kind=8) :: tmp1,tmp2,tnx,tny,rho1,ux1,uy1,uz1

    if(idz==npz-1)then
        !$omp parallel
        !$omp do private(wall_indicator,tmp1,tnx,tny,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,walls)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,nz)

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
                    f16(i,j+1,nz-1))) - rho_out

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
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k
    integer(kind=1) :: wall_indicator
    real(kind=8) :: tmp1,tnx,tny

    if(idz==npz-1)then
        !$omp parallel
        !$omp do private(wall_indicator,tmp1,tnx,tny,i)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,walls)
        !$acc loop collapse(2) device_type(NVIDIA)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,nz)
                
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
                    f17(i,j,nz))) - rho_out

                tnx = 0.5d0*( f2(i,j,nz)+f8(i,j,nz)+f10(i,j,nz)-(f1(i,j,nz)+f7(i,j,nz)+f9(i,j,nz)) )
                tny = 0.5d0*( f4(i,j,nz)+f10(i,j,nz)+f9(i,j,nz)-(f3(i,j,nz)+f7(i,j,nz)+f8(i,j,nz)) )
                f5(i,j,nz) = (f6(i,j,nz)-0.333333333333333333d0*tmp1)* (1d0-wall_indicator)   + f5(i,j,nz)*wall_indicator
                f11(i,j,nz) = (f14(i,j,nz)-0.166666666666666667d0*tmp1 + tnx)* (1d0-wall_indicator)   + f11(i,j,nz)*wall_indicator
                f12(i,j,nz) = (f13(i,j,nz)-0.166666666666666667d0*tmp1 - tnx)* (1d0-wall_indicator)   + f12(i,j,nz)*wall_indicator
                f15(i,j,nz) = (f18(i,j,nz)-0.166666666666666667d0*tmp1 + tny)* (1d0-wall_indicator)   + f15(i,j,nz)*wall_indicator
                f16(i,j,nz) = (f17(i,j,nz)-0.166666666666666667d0*tmp1 - tny)* (1d0-wall_indicator)   + f16(i,j,nz)*wall_indicator
            enddo
        enddo
        !$acc end kernels
        !$omp end parallel
    endif
    return
end subroutine outlet_Zou_He_pressure_BC_after_odd