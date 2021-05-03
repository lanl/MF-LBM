#include "./preprocessor.h"
!===================================================================================================================================================================================
!---------------------- place a porous plate in the domain: 0 - no; 1 - block fluid 1; 2 - block fluid 2  ----------------------
! example1: inject fluid 2 (wetting) during imbibition cycle and block fluid 1 from exiting the inlet
! example2: inject fluid 1 (nonwetting) during drainage cycle and block fluid 1 from exiting the outlet
!===================================================================================================================================================================================
!************************** before odd step kernel *****************************************
subroutine porous_plate_BC_before_odd    !before streaming type BC
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k,zmin,zmax, z_porous_plate_local
    
    zmin = idz*nz + 1
    zmax = idz*nz + nz
    if(z_porous_plate>=zmin.and.z_porous_plate<=zmax)then
        z_porous_plate_local = z_porous_plate - idz*nz
        if(porous_plate_cmd==1)then  !block fluid 1, default (assuming fluid 1 is the nonwetting phase)
            !$omp parallel do private(i)
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
            !$acc kernels present(&
            !$acc &f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
            !$acc &g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18)
            !$acc loop collapse(2) device_type(NVIDIA)
            do j=1,ny
                do i=1,nx          			
                    f6(i,j,z_porous_plate_local) = f5(i, j, z_porous_plate_local-1)
                    f13(i-1,j,z_porous_plate_local) = f12(i,j,z_porous_plate_local-1)
                    f14(i+1,j,z_porous_plate_local) = f11(i,j,z_porous_plate_local-1)
                    f17(i,j-1,z_porous_plate_local) = f16(i,j,z_porous_plate_local-1)
                    f18(i,j+1,z_porous_plate_local) = f15(i,j,z_porous_plate_local-1)
    
                    f5(i,j,z_porous_plate_local) = f6(i, j, z_porous_plate_local+1)
                    f12(i+1,j,z_porous_plate_local) = f13(i,j,z_porous_plate_local+1)
                    f11(i-1,j,z_porous_plate_local) = f14(i,j,z_porous_plate_local+1)
                    f16(i,j+1,z_porous_plate_local) = f17(i,j,z_porous_plate_local+1)
                    f15(i,j-1,z_porous_plate_local) = f18(i,j,z_porous_plate_local+1)
    
                    g6(i,j,z_porous_plate_local) = g6(i, j, z_porous_plate_local+1)
                    g13(i,j,z_porous_plate_local) = g13(i,j,z_porous_plate_local+1)
                    g14(i,j,z_porous_plate_local) = g14(i,j,z_porous_plate_local+1)
                    g17(i,j,z_porous_plate_local) = g17(i,j,z_porous_plate_local+1)
                    g18(i,j,z_porous_plate_local) = g18(i,j,z_porous_plate_local+1)
    
                    g5(i,j,z_porous_plate_local) = g5(i, j, z_porous_plate_local-1)
                    g12(i,j,z_porous_plate_local) = g12(i,j,z_porous_plate_local-1)
                    g11(i,j,z_porous_plate_local) = g11(i,j,z_porous_plate_local-1)
                    g16(i,j,z_porous_plate_local) = g16(i,j,z_porous_plate_local-1)
                    g15(i,j,z_porous_plate_local) = g15(i,j,z_porous_plate_local-1)	                                   
                enddo
            enddo
            !$acc end kernels               
        elseif(porous_plate_cmd==2)then !block fluid 2
            !$omp parallel do private(i)
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
            !$acc kernels present(&
            !$acc &f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
            !$acc &g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18)
            !$acc loop collapse(2) device_type(NVIDIA)
            do j=1,ny
                do i=1,nx          			
                    g6(i,j,z_porous_plate_local) = g5(i, j, z_porous_plate_local-1)
                    g13(i-1,j,z_porous_plate_local) = g12(i,j,z_porous_plate_local-1)
                    g14(i+1,j,z_porous_plate_local) = g11(i,j,z_porous_plate_local-1)
                    g17(i,j-1,z_porous_plate_local) = g16(i,j,z_porous_plate_local-1)
                    g18(i,j+1,z_porous_plate_local) = g15(i,j,z_porous_plate_local-1)
    
                    g5(i,j,z_porous_plate_local) = g6(i, j, z_porous_plate_local+1)
                    g12(i+1,j,z_porous_plate_local) = g13(i,j,z_porous_plate_local+1)
                    g11(i-1,j,z_porous_plate_local) = g14(i,j,z_porous_plate_local+1)
                    g16(i,j+1,z_porous_plate_local) = g17(i,j,z_porous_plate_local+1)
                    g15(i,j-1,z_porous_plate_local) = g18(i,j,z_porous_plate_local+1)
    
                    f6(i,j,z_porous_plate_local) = f6(i, j, z_porous_plate_local+1)
                    f13(i,j,z_porous_plate_local) = f13(i,j,z_porous_plate_local+1)
                    f14(i,j,z_porous_plate_local) = f14(i,j,z_porous_plate_local+1)
                    f17(i,j,z_porous_plate_local) = f17(i,j,z_porous_plate_local+1)
                    f18(i,j,z_porous_plate_local) = f18(i,j,z_porous_plate_local+1)
    
                    f5(i,j,z_porous_plate_local) = f5(i, j, z_porous_plate_local-1)
                    f12(i,j,z_porous_plate_local) = f12(i,j,z_porous_plate_local-1)
                    f11(i,j,z_porous_plate_local) = f11(i,j,z_porous_plate_local-1)
                    f16(i,j,z_porous_plate_local) = f16(i,j,z_porous_plate_local-1)
                    f15(i,j,z_porous_plate_local) = f15(i,j,z_porous_plate_local-1)	                                   
                enddo
            enddo
            !$acc end kernels               
        endif
    endif

    return
end subroutine porous_plate_BC_before_odd

!************************** after odd step kernel *****************************************
subroutine porous_plate_BC_after_odd    !after streaming type BC
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    integer :: i,j,k,zmin,zmax,z_porous_plate_local
    
    zmin = idz*nz + 1
    zmax = idz*nz + nz

    if(z_porous_plate>=zmin.and.z_porous_plate<=zmax)then
        z_porous_plate_local = z_porous_plate - idz*nz  
        if(porous_plate_cmd==1)then  !block fluid 1, default (assuming fluid 1 is the nonwetting phase)
            !$omp parallel do private(i)
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
            !$acc kernels present(&
            !$acc &f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
            !$acc &g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18)
            !$acc loop collapse(2) device_type(NVIDIA)
            do j=1,ny
                do i=1,nx
                    f5(i,j,z_porous_plate_local-1) = f6(i,j,z_porous_plate_local)
                    f11(i,j,z_porous_plate_local-1) = f14(i+1,j,z_porous_plate_local)
                    f12(i,j,z_porous_plate_local-1) = f13(i-1,j,z_porous_plate_local)
                    f15(i,j,z_porous_plate_local-1) = f18(i,j+1,z_porous_plate_local)
                    f16(i,j,z_porous_plate_local-1) = f17(i,j-1,z_porous_plate_local)

                    f6(i,j,z_porous_plate_local+1) = f5(i,j,z_porous_plate_local)
                    f14(i,j,z_porous_plate_local+1) = f11(i-1,j,z_porous_plate_local)
                    f13(i,j,z_porous_plate_local+1) = f12(i+1,j,z_porous_plate_local)
                    f18(i,j,z_porous_plate_local+1) = f15(i,j-1,z_porous_plate_local)
                    f17(i,j,z_porous_plate_local+1) = f16(i,j+1,z_porous_plate_local)

                    g5(i,j,z_porous_plate_local-1) = g5(i,j,z_porous_plate_local)
                    g11(i,j,z_porous_plate_local-1) = g11(i,j,z_porous_plate_local)
                    g12(i,j,z_porous_plate_local-1) = g12(i,j,z_porous_plate_local)
                    g15(i,j,z_porous_plate_local-1) = g15(i,j,z_porous_plate_local)
                    g16(i,j,z_porous_plate_local-1) = g16(i,j,z_porous_plate_local)

                    g6(i,j,z_porous_plate_local+1) = g6(i,j,z_porous_plate_local)
                    g14(i,j,z_porous_plate_local+1) = g14(i,j,z_porous_plate_local)
                    g13(i,j,z_porous_plate_local+1) = g13(i,j,z_porous_plate_local)
                    g18(i,j,z_porous_plate_local+1) = g18(i,j,z_porous_plate_local)
                    g17(i,j,z_porous_plate_local+1) = g17(i,j,z_porous_plate_local)
                enddo
            enddo
            !$acc end kernels
        elseif(porous_plate_cmd==2)then !block fluid 2
            !$omp parallel do private(i)
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_openacc
            !$acc kernels present(&
            !$acc &f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
            !$acc &g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18)
            !$acc loop collapse(2) device_type(NVIDIA)
            do j=1,ny
                do i=1,nx
                    g5(i,j,z_porous_plate_local-1) = g6(i,j,z_porous_plate_local)
                    g11(i,j,z_porous_plate_local-1) = g14(i+1,j,z_porous_plate_local)
                    g12(i,j,z_porous_plate_local-1) = g13(i-1,j,z_porous_plate_local)
                    g15(i,j,z_porous_plate_local-1) = g18(i,j+1,z_porous_plate_local)
                    g16(i,j,z_porous_plate_local-1) = g17(i,j-1,z_porous_plate_local)

                    g6(i,j,z_porous_plate_local+1) = g5(i,j,z_porous_plate_local)
                    g14(i,j,z_porous_plate_local+1) = g11(i-1,j,z_porous_plate_local)
                    g13(i,j,z_porous_plate_local+1) = g12(i+1,j,z_porous_plate_local)
                    g18(i,j,z_porous_plate_local+1) = g15(i,j-1,z_porous_plate_local)
                    g17(i,j,z_porous_plate_local+1) = g16(i,j+1,z_porous_plate_local)

                    f5(i,j,z_porous_plate_local-1) = f5(i,j,z_porous_plate_local)
                    f11(i,j,z_porous_plate_local-1) = f11(i,j,z_porous_plate_local)
                    f12(i,j,z_porous_plate_local-1) = f12(i,j,z_porous_plate_local)
                    f15(i,j,z_porous_plate_local-1) = f15(i,j,z_porous_plate_local)
                    f16(i,j,z_porous_plate_local-1) = f16(i,j,z_porous_plate_local)

                    f6(i,j,z_porous_plate_local+1) = f6(i,j,z_porous_plate_local)
                    f14(i,j,z_porous_plate_local+1) = f14(i,j,z_porous_plate_local)
                    f13(i,j,z_porous_plate_local+1) = f13(i,j,z_porous_plate_local)
                    f18(i,j,z_porous_plate_local+1) = f18(i,j,z_porous_plate_local)
                    f17(i,j,z_porous_plate_local+1) = f17(i,j,z_porous_plate_local)
                enddo
            enddo
            !$acc end kernels
        endif
    endif

    return
end subroutine porous_plate_BC_after_odd
