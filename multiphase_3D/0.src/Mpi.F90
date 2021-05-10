#include "./preprocessor.h"
!=====================================================================================================================================
!------- MPI receive initialization --------
!=====================================================================================================================================
subroutine mpi_irecv_initialization
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    INTEGER :: MPI_ERR

    !faces z
    if(mpi_z)then
        CALL MPI_IRECV( recv_pdf_zP, isize_pdf_z, MPI_DOUBLE_PRECISION, idzP, tag5, MPI_COMM_VGRID, MPI_REQ_z(1), MPI_ERR )
        CALL MPI_IRECV( recv_pdf_zM, isize_pdf_z, MPI_DOUBLE_PRECISION, idzM, tag6, MPI_COMM_VGRID, MPI_REQ_z(2), MPI_ERR )

        CALL MPI_IRECV( recv_phi_zP, isize_phi_z, MPI_DOUBLE_PRECISION, idzP, TAG_phi_zP, MPI_COMM_VGRID, MPI_REQ_phi_Z(1), MPI_ERR )
        CALL MPI_IRECV( recv_phi_zM, isize_phi_z, MPI_DOUBLE_PRECISION, idzM, TAG_phi_zM, MPI_COMM_VGRID, MPI_REQ_phi_Z(2), MPI_ERR )
    endif

    !faces y
    if(mpi_y)then
        CALL MPI_IRECV( recv_pdf_yP, isize_pdf_y, MPI_DOUBLE_PRECISION, idyP, tag3, MPI_COMM_VGRID, MPI_REQ_y(1), MPI_ERR )
        CALL MPI_IRECV( recv_pdf_yM, isize_pdf_y, MPI_DOUBLE_PRECISION, idyM, tag4, MPI_COMM_VGRID, MPI_REQ_y(2), MPI_ERR )

        CALL MPI_IRECV( recv_phi_yP, isize_phi_y, MPI_DOUBLE_PRECISION, idyP, TAG_phi_yP, MPI_COMM_VGRID, MPI_REQ_phi_Y(1), MPI_ERR )
        CALL MPI_IRECV( recv_phi_yM, isize_phi_y, MPI_DOUBLE_PRECISION, idyM, TAG_phi_yM, MPI_COMM_VGRID, MPI_REQ_phi_Y(2), MPI_ERR )
    endif 


    !edges x
    if(mpi_y.and.mpi_z)then  
        CALL MPI_IRECV( recv_pdf_yPzP, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyPzP, TAG7, MPI_COMM_VGRID, MPI_REQ_EX(1), MPI_ERR )
        CALL MPI_IRECV( recv_pdf_yPzM, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyPzM, TAG8, MPI_COMM_VGRID, MPI_REQ_EX(2), MPI_ERR )
        CALL MPI_IRECV( recv_pdf_yMzP, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyMzP, TAG9, MPI_COMM_VGRID, MPI_REQ_EX(3), MPI_ERR )
        CALL MPI_IRECV( recv_pdf_yMzM, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyMzM, TAG10, MPI_COMM_VGRID, MPI_REQ_EX(4), MPI_ERR )

        CALL MPI_IRECV( recv_phi_yPzP, isize_phi_ex, MPI_DOUBLE_PRECISION,idyPzP, TAG_phi_yP_zP, MPI_COMM_VGRID, MPI_REQ_phi_EX(1), MPI_ERR )
        CALL MPI_IRECV( recv_phi_yMzP, isize_phi_ex, MPI_DOUBLE_PRECISION,idyMzP, TAG_phi_yM_zP, MPI_COMM_VGRID, MPI_REQ_phi_EX(2), MPI_ERR )
        CALL MPI_IRECV( recv_phi_yPzM, isize_phi_ex, MPI_DOUBLE_PRECISION,idyPzM, TAG_phi_yP_zM, MPI_COMM_VGRID, MPI_REQ_phi_EX(3), MPI_ERR )
        CALL MPI_IRECV( recv_phi_yMzM, isize_phi_ex, MPI_DOUBLE_PRECISION,idyMzM, TAG_phi_yM_zM, MPI_COMM_VGRID, MPI_REQ_phi_EX(4), MPI_ERR )
    endif

    return
end subroutine mpi_irecv_initialization

!=====================================================================================================================================
!------- MPI send request --------
!=====================================================================================================================================
subroutine mpi_send_req
    use mpi_variable
    use Misc_module
    IMPLICIT NONE
    include 'mpif.h'
    INTEGER :: MPI_ERR

    !face z
    if(mpi_z)then
        !$acc wait(z_pdf_update_async)
        CALL MPI_Isend( send_pdf_zM, isize_pdf_z, MPI_DOUBLE_PRECISION, idzM, tag5, MPI_COMM_VGRID, MPI_REQ_z(3), MPI_ERR )
        CALL MPI_Isend( send_pdf_zP, isize_pdf_z, MPI_DOUBLE_PRECISION, idzP, tag6, MPI_COMM_VGRID, MPI_REQ_z(4), MPI_ERR )

        !$acc wait(z_phi_update_async)
        CALL MPI_Isend( send_phi_zM, isize_phi_z, MPI_DOUBLE_PRECISION, idzM, TAG_phi_zP, MPI_COMM_VGRID, MPI_REQ_phi_Z(3), MPI_ERR )
        CALL MPI_Isend( send_phi_zP, isize_phi_z, MPI_DOUBLE_PRECISION, idzP, TAG_phi_zM, MPI_COMM_VGRID, MPI_REQ_phi_Z(4), MPI_ERR )
    endif

    !face y
    if(mpi_y)then
        !$acc wait(y_pdf_update_async)
        CALL MPI_Isend( send_pdf_yM, isize_pdf_y, MPI_DOUBLE_PRECISION, idyM, tag3, MPI_COMM_VGRID, MPI_REQ_y(3), MPI_ERR )
        CALL MPI_Isend( send_pdf_yP, isize_pdf_y, MPI_DOUBLE_PRECISION, idyP, tag4, MPI_COMM_VGRID, MPI_REQ_y(4), MPI_ERR )

        !$acc wait(y_phi_update_async)
        CALL MPI_Isend( send_phi_yM, isize_phi_y, MPI_DOUBLE_PRECISION, idyM, TAG_phi_yP, MPI_COMM_VGRID, MPI_REQ_phi_Y(3), MPI_ERR )
        CALL MPI_Isend( send_phi_yP, isize_phi_y, MPI_DOUBLE_PRECISION, idyP, TAG_phi_yM, MPI_COMM_VGRID, MPI_REQ_phi_Y(4), MPI_ERR )
    endif

    !$acc wait(edge_pdf_update_async)
    !$acc wait(edge_phi_update_async)

    !edge x
    if(mpi_y.and.mpi_z)then 
        CALL MPI_Isend( send_pdf_yMzM, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyMzM, TAG7, MPI_COMM_VGRID, MPI_REQ_EX(5), MPI_ERR )
        CALL MPI_Isend( send_pdf_yMzP, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyMzP, TAG8, MPI_COMM_VGRID, MPI_REQ_EX(6), MPI_ERR )
        CALL MPI_Isend( send_pdf_yPzM, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyPzM, TAG9, MPI_COMM_VGRID, MPI_REQ_EX(7), MPI_ERR )
        CALL MPI_Isend( send_pdf_yPzP, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyPzP, TAG10, MPI_COMM_VGRID, MPI_REQ_EX(8), MPI_ERR )

        CALL MPI_Isend( send_phi_yMzM, isize_phi_ex, MPI_DOUBLE_PRECISION,idyMzM, TAG_phi_yP_zP, MPI_COMM_VGRID, MPI_REQ_phi_EX(5), MPI_ERR )
        CALL MPI_Isend( send_phi_yPzM, isize_phi_ex, MPI_DOUBLE_PRECISION,idyPzM, TAG_phi_yM_zP, MPI_COMM_VGRID, MPI_REQ_phi_EX(6), MPI_ERR )
        CALL MPI_Isend( send_phi_yMzP, isize_phi_ex, MPI_DOUBLE_PRECISION,idyMzP, TAG_phi_yP_zM, MPI_COMM_VGRID, MPI_REQ_phi_EX(7), MPI_ERR )
        CALL MPI_Isend( send_phi_yPzP, isize_phi_ex, MPI_DOUBLE_PRECISION,idyPzP, TAG_phi_yM_zM, MPI_COMM_VGRID, MPI_REQ_phi_EX(8), MPI_ERR )
    endif 

    return
end subroutine mpi_send_req


!=====================================================================================================================================
!------- pdf data packing (pull type) --------
!=====================================================================================================================================
subroutine mpi_pdf_pull_async_pack
    use Misc_module
    use Fluid_multiphase
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k
    INTEGER :: MPI_ERR

    !$OMP PARALLEL
    
    if(mpi_z)then
        !****************************z direction******************************************
        !$OMP DO private(i)
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc & g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,send_pdf_zP,send_pdf_zM)async(LBM_async_z)
        !$acc loop
        do j = 1, ny
            do i = 1,  nx
                send_pdf_zM(i,j,1) = f6(i,j,1)
                send_pdf_zM(i,j,2) = f14(i,j,1)
                send_pdf_zM(i,j,3) = f13(i,j,1)
                send_pdf_zM(i,j,4) = f18(i,j,1)
                send_pdf_zM(i,j,5) = f17(i,j,1)

                send_pdf_zP(i,j,1) = f5(i,j,nz)
                send_pdf_zP(i,j,2) = f11(i,j,nz)
                send_pdf_zP(i,j,3) = f12(i,j,nz)
                send_pdf_zP(i,j,4) = f15(i,j,nz)
                send_pdf_zP(i,j,5) = f16(i,j,nz)

                send_pdf_zM(i,j,6) = g6(i,j,1)
                send_pdf_zM(i,j,7) = g14(i,j,1)
                send_pdf_zM(i,j,8) = g13(i,j,1)
                send_pdf_zM(i,j,9) = g18(i,j,1)
                send_pdf_zM(i,j,10) = g17(i,j,1)

                send_pdf_zP(i,j,6) = g5(i,j,nz)
                send_pdf_zP(i,j,7) = g11(i,j,nz)
                send_pdf_zP(i,j,8) = g12(i,j,nz)
                send_pdf_zP(i,j,9) = g15(i,j,nz)
                send_pdf_zP(i,j,10) = g16(i,j,nz)
            enddo
        enddo
        !$acc end kernels
        !$acc update host(send_pdf_zP,send_pdf_zM)async(z_pdf_update_async)wait(LBM_async_z)
    endif

    if(mpi_y)then
        !****************************y direction******************************************
        !$OMP DO private(i)
        !the reason that acc_wait async_z is due to y halo calculation is not complete calculation. z-halo caculation is complete halo layer calulation
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,send_pdf_yP,send_pdf_yM)async(LBM_async_y)wait(LBM_async_z)
        !$acc loop
        do k = 1, nz
            do i = 1, nx
                send_pdf_yM(i,k,1) = f4(i,1,k)
                send_pdf_yM(i,k,2) = f10(i,1,k)
                send_pdf_yM(i,k,3) = f9(i,1,k)
                send_pdf_yM(i,k,4) = f16(i,1,k)
                send_pdf_yM(i,k,5) = f18(i,1,k)

                send_pdf_yP(i,k,1) = f3(i,ny,k)
                send_pdf_yP(i,k,2) = f7(i,ny,k)
                send_pdf_yP(i,k,3) = f8(i,ny,k)
                send_pdf_yP(i,k,4) = f15(i,ny,k)
                send_pdf_yP(i,k,5) = f17(i,ny,k)

                send_pdf_yM(i,k,6) = g4(i,1,k)
                send_pdf_yM(i,k,7) = g10(i,1,k)
                send_pdf_yM(i,k,8) = g9(i,1,k)
                send_pdf_yM(i,k,9) = g16(i,1,k)
                send_pdf_yM(i,k,10) = g18(i,1,k)

                send_pdf_yP(i,k,6) = g3(i,ny,k)
                send_pdf_yP(i,k,7) = g7(i,ny,k)
                send_pdf_yP(i,k,8) = g8(i,ny,k)
                send_pdf_yP(i,k,9) = g15(i,ny,k)
                send_pdf_yP(i,k,10) = g17(i,ny,k)
            enddo
        enddo
        !$acc end kernels
        !$acc update host(send_pdf_yP,send_pdf_yM)async(y_pdf_update_async)wait(LBM_async_y)
    endif

    !****************************edges******************************************
    if(mpi_y.and.mpi_z)then
        ! Edges along X
        !$OMP DO
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,&
        !$acc & send_pdf_yPzP,send_pdf_yPzM,send_pdf_yMzP,send_pdf_yMzM)async(edge_pdf_async)wait(LBM_async_z,LBM_async_y)
        !$acc loop
        DO i = 1, NX
            send_pdf_yPzP(i) = f15(i,ny,nz)
            send_pdf_yPzM(i) = f17(i,ny,1 )
            send_pdf_yMzP(i) = f16(i,1 ,nz)
            send_pdf_yMzM(i) = f18(i,1 ,1 )

            send_pdf_yPzP(i+nx) = g15(i,ny,nz)
            send_pdf_yPzM(i+nx) = g17(i,ny,1 )
            send_pdf_yMzP(i+nx) = g16(i,1 ,nz)
            send_pdf_yMzM(i+nx) = g18(i,1 ,1 )
        END DO
        !$acc end kernels
        !$acc update host(send_pdf_yPzP,send_pdf_yPzM,send_pdf_yMzP,send_pdf_yMzM)async(edge_pdf_update_async)wait(edge_pdf_async)
    endif

    !$OMP END PARALLEL

    return
end subroutine mpi_pdf_pull_async_pack


!=====================================================================================================================================
!------- pdf data updating (pull type) --------
!=====================================================================================================================================
subroutine mpi_pdf_pull_sync_update
    use Misc_module
    use Fluid_multiphase
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k
    INTEGER :: MPI_ERR

    !$OMP PARALLEL

    if(mpi_z)then
        !$OMP SINGLE 
        CALL MPI_WAITALL(4, MPI_REQ_Z, MPI_STAT, MPI_ERR)
        !$OMP END SINGLE
        !$OMP BARRIER

        ! Z direction 
        !$OMP DO private(i)
        !$acc update device(recv_pdf_zM,recv_pdf_zP)async(z_pdf_update_async)
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc & g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,send_pdf_zP,send_pdf_zM)async(z_pdf_update_async)
        !$acc loop
        do j = 1, ny
            do i = 1, nx
                f5(i,j,1-1) = recv_pdf_zM(i,j,1)
                f11(i,j,1-1) = recv_pdf_zM(i,j,2)
                f12(i,j,1-1) = recv_pdf_zM(i,j,3)
                f15(i,j,1-1) = recv_pdf_zM(i,j,4)
                f16(i,j,1-1) = recv_pdf_zM(i,j,5)

                f6(i,j,1+nz) = recv_pdf_zP(i,j,1)
                f14(i,j,1+nz) = recv_pdf_zP(i,j,2)
                f13(i,j,1+nz) = recv_pdf_zP(i,j,3)
                f18(i,j,1+nz) = recv_pdf_zP(i,j,4)
                f17(i,j,1+nz) = recv_pdf_zP(i,j,5)

                g5(i,j,1-1) = recv_pdf_zM(i,j,6)
                g11(i,j,1-1) = recv_pdf_zM(i,j,7)
                g12(i,j,1-1) = recv_pdf_zM(i,j,8)
                g15(i,j,1-1) = recv_pdf_zM(i,j,9)
                g16(i,j,1-1) = recv_pdf_zM(i,j,10)

                g6(i,j,1+nz) = recv_pdf_zP(i,j,6)
                g14(i,j,1+nz) = recv_pdf_zP(i,j,7)
                g13(i,j,1+nz) = recv_pdf_zP(i,j,8)
                g18(i,j,1+nz) = recv_pdf_zP(i,j,9)
                g17(i,j,1+nz) = recv_pdf_zP(i,j,10)
            enddo
        enddo
        !$acc end kernels
    endif

    if(mpi_y)then
        !$OMP SINGLE
        CALL MPI_WAITALL(4, MPI_REQ_Y, MPI_STAT, MPI_ERR) 
        !$OMP END SINGLE
        !$OMP BARRIER

        ! Y direction 
        !$OMP DO private(i)
        !$acc update device(recv_pdf_yM,recv_pdf_yP)async(y_pdf_update_async)
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc & g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,send_pdf_yP,send_pdf_yM)async(y_pdf_update_async)
        !$acc loop
        do k = 1, nz
            do i = 1, nx
                f3(i,1-1,k) = recv_pdf_yM(i,k,1)
                f7(i,1-1,k) = recv_pdf_yM(i,k,2)
                f8(i,1-1,k) = recv_pdf_yM(i,k,3)
                f15(i,1-1,k) = recv_pdf_yM(i,k,4)
                f17(i,1-1,k) = recv_pdf_yM(i,k,5)

                f4(i,1+ny,k) = recv_pdf_yP(i,k,1)
                f10(i,1+ny,k) = recv_pdf_yP(i,k,2)
                f9(i,1+ny,k) = recv_pdf_yP(i,k,3)
                f16(i,1+ny,k) = recv_pdf_yP(i,k,4)
                f18(i,1+ny,k) = recv_pdf_yP(i,k,5)

                g3(i,1-1,k) = recv_pdf_yM(i,k,6)
                g7(i,1-1,k) = recv_pdf_yM(i,k,7)
                g8(i,1-1,k) = recv_pdf_yM(i,k,8)
                g15(i,1-1,k) = recv_pdf_yM(i,k,9)
                g17(i,1-1,k) = recv_pdf_yM(i,k,10)

                g4(i,1+ny,k) = recv_pdf_yP(i,k,6)
                g10(i,1+ny,k) = recv_pdf_yP(i,k,7)
                g9(i,1+ny,k) = recv_pdf_yP(i,k,8)
                g16(i,1+ny,k) = recv_pdf_yP(i,k,9)
                g18(i,1+ny,k) = recv_pdf_yP(i,k,10)
            enddo
        enddo
        !$acc end kernels
    endif



    if(mpi_y.and.mpi_z)then 
        !$OMP SINGLE
        CALL MPI_WAITALL(8, MPI_REQ_EX, MPI_ESTAT, MPI_ERR)
        !$OMP END SINGLE
        !$OMP BARRIER
        
        !$OMP DO
        !$acc update device(recv_pdf_yPzP,recv_pdf_yPzM,recv_pdf_yMzP,recv_pdf_yMzM)async(edge_pdf_update_async)
        ! Edges along x, acc wait is required because face communication will contaminate edge communication !!!
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc & g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,&
        !$acc & recv_pdf_yPzP,recv_pdf_yPzM,recv_pdf_yMzP,recv_pdf_yMzM)wait(z_pdf_update_async,y_pdf_update_async)
        !$acc loop
        DO i = 1, NX
            f18(i,ny+1,nz+1) = recv_pdf_yPzP(i)
            f16(i,ny+1,1-1 ) = recv_pdf_yPzM(i)
            f17(i,1-1 ,nz+1) = recv_pdf_yMzP(i)
            f15(i,1-1 ,1-1 ) = recv_pdf_yMzM(i)

            g18(i,ny+1,nz+1) = recv_pdf_yPzP(i+nx)
            g16(i,ny+1,1-1 ) = recv_pdf_yPzM(i+nx)
            g17(i,1-1 ,nz+1) = recv_pdf_yMzP(i+nx)
            g15(i,1-1 ,1-1 ) = recv_pdf_yMzM(i+nx)
        END DO
        !$acc end kernels
    endif

    !$OMP END PARALLEL

    return
end subroutine mpi_pdf_pull_sync_update
!**********************MPI PDF pull*****************************



!=====================================================================================================================================
!------- pdf data packing (push type) --------
!=====================================================================================================================================
subroutine mpi_pdf_push_async_pack
    use Misc_module
    use Fluid_multiphase
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k
    INTEGER :: MPI_ERR

    !$OMP PARALLEL
    
    if(mpi_z)then
        !****************************z direction******************************************
        !$OMP DO private(i)
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,send_pdf_zP,send_pdf_zM)async(LBM_async_z)
        !$acc loop
        do j = 1, ny
            !DIR$ vector always
            do i = 1,  nx
                send_pdf_zM(i,j,1) = f5(i,j,1-1)
                send_pdf_zM(i,j,2) = f11(i,j,1-1)
                send_pdf_zM(i,j,3) = f12(i,j,1-1)
                send_pdf_zM(i,j,4) = f15(i,j,1-1)
                send_pdf_zM(i,j,5) = f16(i,j,1-1)

                send_pdf_zP(i,j,1) = f6(i,j,nz+1)
                send_pdf_zP(i,j,2) = f14(i,j,nz+1)
                send_pdf_zP(i,j,3) = f13(i,j,nz+1)
                send_pdf_zP(i,j,4) = f18(i,j,nz+1)
                send_pdf_zP(i,j,5) = f17(i,j,nz+1)

                send_pdf_zM(i,j,6) = g5(i,j,1-1)
                send_pdf_zM(i,j,7) = g11(i,j,1-1)
                send_pdf_zM(i,j,8) = g12(i,j,1-1)
                send_pdf_zM(i,j,9) = g15(i,j,1-1)
                send_pdf_zM(i,j,10) = g16(i,j,1-1)

                send_pdf_zP(i,j,6) = g6(i,j,nz+1)
                send_pdf_zP(i,j,7) = g14(i,j,nz+1)
                send_pdf_zP(i,j,8) = g13(i,j,nz+1)
                send_pdf_zP(i,j,9) = g18(i,j,nz+1)
                send_pdf_zP(i,j,10) = g17(i,j,nz+1)
            enddo
        enddo
        !$acc end kernels
        !$acc update host(send_pdf_zP,send_pdf_zM)async(z_pdf_update_async)wait(LBM_async_z)
    endif

    if(mpi_y)then
        !****************************y direction******************************************
        !$OMP DO private(i)
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,send_pdf_yP,send_pdf_yM)async(LBM_async_y)wait(LBM_async_z)
        !$acc loop
        do k = 1, nz
            !DIR$ vector always
            do i = 1, nx
                send_pdf_yM(i,k,1) = f3(i,1-1,k)
                send_pdf_yM(i,k,2) = f7(i,1-1,k)
                send_pdf_yM(i,k,3) = f8(i,1-1,k)
                send_pdf_yM(i,k,4) = f15(i,1-1,k)
                send_pdf_yM(i,k,5) = f17(i,1-1,k)

                send_pdf_yP(i,k,1) = f4(i,ny+1,k)
                send_pdf_yP(i,k,2) = f10(i,ny+1,k)
                send_pdf_yP(i,k,3) = f9(i,ny+1,k)
                send_pdf_yP(i,k,4) = f16(i,ny+1,k)
                send_pdf_yP(i,k,5) = f18(i,ny+1,k)

                send_pdf_yM(i,k,6) = g3(i,1-1,k)
                send_pdf_yM(i,k,7) = g7(i,1-1,k)
                send_pdf_yM(i,k,8) = g8(i,1-1,k)
                send_pdf_yM(i,k,9) = g15(i,1-1,k)
                send_pdf_yM(i,k,10) = g17(i,1-1,k)

                send_pdf_yP(i,k,6) = g4(i,ny+1,k)
                send_pdf_yP(i,k,7) = g10(i,ny+1,k)
                send_pdf_yP(i,k,8) = g9(i,ny+1,k)
                send_pdf_yP(i,k,9) = g16(i,ny+1,k)
                send_pdf_yP(i,k,10) = g18(i,ny+1,k)
            enddo
        enddo
        !$acc end kernels
        !$acc update host(send_pdf_yP,send_pdf_yM)async(y_pdf_update_async)wait(LBM_async_y)
    endif
    

    !****************************edges******************************************

    if(mpi_y.and.mpi_z)then
        ! Edges along X
        !$OMP DO 
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,&
        !$acc & send_pdf_yPzP,send_pdf_yPzM,send_pdf_yMzP,send_pdf_yMzM)async(edge_pdf_async)wait(LBM_async_z,LBM_async_y)
        !$acc loop
        DO i = 1, NX
            send_pdf_yPzP(i) = f18(i,ny+1,nz+1)
            send_pdf_yPzM(i) = f16(i,ny+1,1 -1)
            send_pdf_yMzP(i) = f17(i,1 -1,nz+1)
            send_pdf_yMzM(i) = f15(i,1 -1,1 -1)

            send_pdf_yPzP(i+nx) = g18(i,ny+1,nz+1)
            send_pdf_yPzM(i+nx) = g16(i,ny+1,1 -1)
            send_pdf_yMzP(i+nx) = g17(i,1 -1,nz+1)
            send_pdf_yMzM(i+nx) = g15(i,1 -1,1 -1)
        END DO
        !$acc end kernels
        !$acc update host(send_pdf_yPzP,send_pdf_yPzM,send_pdf_yMzP,send_pdf_yMzM)async(edge_pdf_update_async)wait(edge_pdf_async)
    endif

    !$OMP END PARALLEL

    return
end subroutine mpi_pdf_push_async_pack

!=====================================================================================================================================
!------- pdf data updating (push type) --------
!=====================================================================================================================================
subroutine mpi_pdf_push_sync_update
    use Misc_module
    use Fluid_multiphase
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k
    INTEGER :: MPI_ERR

    !$OMP PARALLEL

    if(mpi_z)then
        ! Z direction
        !$OMP SINGLE 
        CALL MPI_WAITALL(4, MPI_REQ_Z, MPI_STAT, MPI_ERR)
        !$OMP END SINGLE
        !$OMP BARRIER

        !$OMP DO private(i)
        !$acc update device(recv_pdf_zM,recv_pdf_zP)async(z_pdf_update_async)
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,send_pdf_zP,send_pdf_zM)async(z_pdf_update_async)
        !$acc loop
        do j = 1, ny
            do i = 1, nx
                f6(i,j,1) = recv_pdf_zM(i,j,1)
                f14(i,j,1) = recv_pdf_zM(i,j,2)
                f13(i,j,1) = recv_pdf_zM(i,j,3)
                f18(i,j,1) = recv_pdf_zM(i,j,4)
                f17(i,j,1) = recv_pdf_zM(i,j,5)

                f5(i,j,nz) = recv_pdf_zP(i,j,1)
                f11(i,j,nz) = recv_pdf_zP(i,j,2)
                f12(i,j,nz) = recv_pdf_zP(i,j,3)
                f15(i,j,nz) = recv_pdf_zP(i,j,4)
                f16(i,j,nz) = recv_pdf_zP(i,j,5)

                g6(i,j,1) = recv_pdf_zM(i,j,6)
                g14(i,j,1) = recv_pdf_zM(i,j,7)
                g13(i,j,1) = recv_pdf_zM(i,j,8)
                g18(i,j,1) = recv_pdf_zM(i,j,9)
                g17(i,j,1) = recv_pdf_zM(i,j,10)

                g5(i,j,nz) = recv_pdf_zP(i,j,6)
                g11(i,j,nz) = recv_pdf_zP(i,j,7)
                g12(i,j,nz) = recv_pdf_zP(i,j,8)
                g15(i,j,nz) = recv_pdf_zP(i,j,9)
                g16(i,j,nz) = recv_pdf_zP(i,j,10)
            enddo
        enddo
        !$acc end kernels
    endif

    if(mpi_y)then
        ! Y direction

        !$OMP SINGLE 
        CALL MPI_WAITALL(4, MPI_REQ_Y, MPI_STAT, MPI_ERR) 
        !$OMP END SINGLE
        !$OMP BARRIER
        
        !$OMP DO private(i)
        !$acc update device(recv_pdf_yM,recv_pdf_yP)async(y_pdf_update_async)
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,send_pdf_yP,send_pdf_yM)async(y_pdf_update_async)
        !$acc loop
        do k = 1, nz
            do i = 1, nx
                f4(i,1,k) = recv_pdf_yM(i,k,1)
                f10(i,1,k) = recv_pdf_yM(i,k,2)
                f9(i,1,k) = recv_pdf_yM(i,k,3)
                f16(i,1,k) = recv_pdf_yM(i,k,4)
                f18(i,1,k) = recv_pdf_yM(i,k,5)

                f3(i,ny,k) = recv_pdf_yP(i,k,1)
                f7(i,ny,k) = recv_pdf_yP(i,k,2)
                f8(i,ny,k) = recv_pdf_yP(i,k,3)
                f15(i,ny,k) = recv_pdf_yP(i,k,4)
                f17(i,ny,k) = recv_pdf_yP(i,k,5)

                g4(i,1,k) = recv_pdf_yM(i,k,6)
                g10(i,1,k) = recv_pdf_yM(i,k,7)
                g9(i,1,k) = recv_pdf_yM(i,k,8)
                g16(i,1,k) = recv_pdf_yM(i,k,9)
                g18(i,1,k) = recv_pdf_yM(i,k,10)

                g3(i,ny,k) = recv_pdf_yP(i,k,6)
                g7(i,ny,k) = recv_pdf_yP(i,k,7)
                g8(i,ny,k) = recv_pdf_yP(i,k,8)
                g15(i,ny,k) = recv_pdf_yP(i,k,9)
                g17(i,ny,k) = recv_pdf_yP(i,k,10)
            enddo
        enddo
        !$acc end kernels
    endif
    

    if(mpi_y.and.mpi_z)then  
        ! Edges along x
        !$OMP SINGLE 
        CALL MPI_WAITALL(8, MPI_REQ_EX, MPI_ESTAT, MPI_ERR)
        !$OMP END SINGLE
        !$OMP BARRIER
    
        !$OMP DO 
        !$acc update device(recv_pdf_yPzP,recv_pdf_yPzM,recv_pdf_yMzP,recv_pdf_yMzM)async(edge_pdf_update_async)
        ! Edges along x, acc wait is required because face communication will contaminate edge communication !!!
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,&
        !$acc & recv_pdf_yPzP,recv_pdf_yPzM,recv_pdf_yMzP,recv_pdf_yMzM)async(edge_pdf_update_async)wait(z_pdf_update_async,y_pdf_update_async)
        !$acc loop
        DO i = 1, NX
            f15(i,ny,nz) = recv_pdf_yPzP(i)
            f17(i,ny,1 ) = recv_pdf_yPzM(i)
            f16(i,1 ,nz) = recv_pdf_yMzP(i)
            f18(i,1 ,1 ) = recv_pdf_yMzM(i)

            g15(i,ny,nz) = recv_pdf_yPzP(i+nx)
            g17(i,ny,1 ) = recv_pdf_yPzM(i+nx)
            g16(i,1 ,nz) = recv_pdf_yMzP(i+nx)
            g18(i,1 ,1 ) = recv_pdf_yMzM(i+nx)
        END DO
        !$acc end kernels
    endif

    !$OMP END PARALLEL

    return
end subroutine mpi_pdf_push_sync_update






!=====================================================================================================================================
!------- phi data packing --------
!=====================================================================================================================================
subroutine MPI_phi_async_pack
    use Misc_module
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k

    !$OMP PARALLEL

    !**************************** faces ******************************************

    if(mpi_z)then
        !$OMP DO private(i,j)
        !$acc kernels present(phi,send_phi_zM,send_phi_zP)async(LBM_async_z)
        !$acc loop collapse(2)
        do k=1,overlap_phi
            do j = 1, ny
                do i = 1, nx
                    send_phi_zM(i,j,k) = phi(i,j,k)
                    send_phi_zP(i,j,k) = phi(i,j,nz+k-overlap_phi)
                enddo
            enddo
        enddo
        !$acc end kernels
        !$acc update host(send_phi_zM,send_phi_zP)async(z_phi_update_async)wait(LBM_async_z)
    endif
        
    if(mpi_y)then
        !$OMP DO private(i,j)
        !$acc kernels present(phi,send_phi_yM,send_phi_yP)async(LBM_async_y)wait(LBM_async_z)
        !$acc loop collapse(2)
        do k = 1, nz
            do j = 1,overlap_phi
                do i = 1, nx
                    send_phi_yM(i,j,k) = phi(i,j,k)
                    send_phi_yP(i,j,k) = phi(i,ny+j-overlap_phi,k)
                enddo
            enddo
        enddo
        !$acc end kernels
        !$acc update host(send_phi_yM,send_phi_yP)async(y_phi_update_async)wait(LBM_async_y)
    endif

    !**************************** edges ******************************************

    if(mpi_y.and.mpi_z)then  
        !$OMP DO private(i,j)
        !!$acc wait(LBM_kernels_z,LBM_kernels_y)
        !$acc kernels present(phi,send_phi_yMzM,send_phi_yPzP,send_phi_yMzP,send_phi_yPzM)async(edge_phi_async)wait(LBM_async_z,LBM_async_y)
        !$acc loop collapse(2)
        do k = 1, overlap_phi
            do j=1,overlap_phi
                do i = 1, nx
                    send_phi_yMzM(i,j,k) = phi(i,j,k)
                    send_phi_yPzP(i,j,k) = phi(i,ny+j-overlap_phi,nz+k-overlap_phi)
                    send_phi_yMzP(i,j,k) = phi(i,j,nz+k-overlap_phi)
                    send_phi_yPzM(i,j,k) = phi(i,ny+j-overlap_phi,k)
                enddo
            enddo
        enddo
        !$acc end kernels
        !$acc update host(send_phi_yMzM,send_phi_yPzP,send_phi_yMzP,send_phi_yPzM)async(edge_phi_update_async)wait(edge_phi_async)
    endif

    !$OMP END PARALLEL

    return
end subroutine MPI_phi_async_pack


!=====================================================================================================================================
!------- phi data updating --------
!=====================================================================================================================================
subroutine MPI_phi_sync_update
    use Misc_module
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k,lx,ly,lz
    INTEGER :: MPI_ERR
    
    !$OMP PARALLEL

    !**************************** face update ******************************************

    if(mpi_z)then
        !$OMP SINGLE 
        CALL MPI_WAITALL(4, MPI_REQ_phi_Z, MPI_STAT, MPI_ERR)   
        !$OMP END SINGLE
        !$OMP BARRIER   

        !$acc update device(recv_phi_zM,recv_phi_zP)async(z_phi_update_async)
        if(kper.eq.1 .or. idz.ne.0)then   !  iper indicates perodic condition
            !$OMP DO private(i,j)
            !$acc kernels present(phi,recv_phi_zM)async(z_phi_update_async)
            !$acc loop  collapse(2)
            do k=1,overlap_phi       
                do j = 1, ny
                    do i = 1, nx
                        phi(i,j,k-overlap_phi ) = recv_phi_zM(i,j,k)
                    enddo
                enddo
            enddo
            !$acc end kernels
        endif
        if(kper.eq.1 .or. idz.ne.npz-1)then
            !$OMP DO private(i,j)
            !$acc kernels present(phi,recv_phi_zP)async(z_phi_update_async)
            !$acc loop  collapse(2)
            do k=1,overlap_phi            
                do j = 1, ny
                    do i = 1, nx
                        phi(i,j,k+nz ) = recv_phi_zP(i,j,k)
                    enddo
                enddo
            enddo
            !$acc end kernels
        endif 
    endif

    if(mpi_y)then  
        !$OMP SINGLE
        CALL MPI_WAITALL(4, MPI_REQ_phi_Y, MPI_STAT, MPI_ERR)  
        !$OMP END SINGLE
        !$OMP BARRIER

        !$acc update device(recv_phi_yM,recv_phi_yP)async(y_phi_update_async)
        if(jper.eq.1 .or. idy.ne.0)then   !  iper indicates perodic condition      
            !$OMP DO private(i,j)
            !$acc kernels present(phi,recv_phi_yM)async(y_phi_update_async)
            !$acc loop collapse(2)
            do k = 1, nz
                do j = 1, overlap_phi
                    do i = 1, nx
                        phi(i,j-overlap_phi,k ) = recv_phi_yM(i,j,k)
                    enddo
                enddo
            enddo
            !$acc end kernels
        endif
        if(jper.eq.1 .or. idy.ne.npy-1)then
            !$OMP DO private(i,j)
            !$acc kernels present(phi,recv_phi_yP)async(y_phi_update_async)
            !$acc loop collapse(2)
            do k = 1, nz
                do j = 1, overlap_phi
                    do i = 1, nx
                        phi(i,j+ny,k ) = recv_phi_yP(i,j,k)
                    enddo
                enddo
            enddo
            !$acc end kernels
        endif
    endif


    !**************************** edge update ******************************************
    if(mpi_y.and.mpi_z)then 
        !$OMP SINGLE
        CALL MPI_WAITALL(8, MPI_REQ_phi_EX, MPI_ESTAT, MPI_ERR)
        !$OMP END SINGLE
        !$OMP BARRIER

        !$acc update device(recv_phi_yPzP,recv_phi_yPzM,recv_phi_yMzP,recv_phi_yMzM)async(edge_phi_update_async)
        if((jper.eq.1 .or. idy.ne.0).and.(kper.eq.1 .or. idz.ne.0))then 
            !$OMP DO private(i,j)  
            !$acc kernels present(phi,recv_phi_yMzM)async(edge_phi_update_async)
            !$acc loop collapse(2)
            do k = 1,overlap_phi
                do j=1,overlap_phi
                    do i = 1, nx
                        phi(i,j-overlap_phi,k-overlap_phi) = recv_phi_yMzM(i,j,k)
                    enddo
                enddo
            enddo
            !$acc end kernels
        endif
        if((jper.eq.1 .or. idy.ne.npy-1).and.(kper.eq.1 .or. idz.ne.0))then   
            !$OMP DO private(i,j)
            !$acc kernels present(phi,recv_phi_yPzM)async(edge_phi_update_async)
            !$acc loop collapse(2)
            do k = 1,overlap_phi
                do j=1,overlap_phi
                    do i = 1, nx
                        phi(i,j+ny,k-overlap_phi) = recv_phi_yPzM(i,j,k)
                    enddo
                enddo
            enddo
            !$acc end kernels
        endif
        if((jper.eq.1 .or. idy.ne.npy-1).and.(kper.eq.1 .or. idz.ne.npz-1))then
            !$OMP DO private(i,j)
            !$acc kernels present(phi,recv_phi_yPzP)async(edge_phi_update_async)
            !$acc loop collapse(2)
            do k = 1,overlap_phi
                do j=1,overlap_phi
                    do i = 1, nx
                        phi(i,j+ny,k+nz) = recv_phi_yPzP(i,j,k)
                    enddo
                enddo
            enddo
            !$acc end kernels
        endif
        if((jper.eq.1 .or. idy.ne.0).and.(kper.eq.1 .or. idz.ne.npz-1))then
            !$OMP DO private(i,j)
            !$acc kernels present(phi,recv_phi_yMzP)async(edge_phi_update_async)
            !$acc loop collapse(2)
            do k = 1,overlap_phi
                do j=1,overlap_phi
                    do i = 1, nx
                        phi(i,j-overlap_phi,k+nz) = recv_phi_yMzP(i,j,k)
                    enddo
                enddo
            enddo
        !$acc end kernels
        endif
    endif

    !#if defined(z_mpi)&&defined(x_mpi)
    !        CALL MPI_WAITALL(8, MPI_REQ_phi_EY, MPI_ESTAT, MPI_ERR)
    !    !$acc update device(recv_phi_xPzP,recv_phi_xPzM,recv_phi_xMzP,recv_phi_xMzM)async(y_edge_phi_recv_unpack_async)
    !    !$acc kernels present(phi,recv_phi_xPzP,recv_phi_xPzM,recv_phi_xMzP,recv_phi_xMzM)async(y_edge_phi_recv_unpack_async)
    !    !$acc loop collapse(3)
    !    do k = 1,overlap_phi
    !        do j=1,ny
    !            do i = 1,overlap_phi
    !                phi(i-overlap_phi,j,k-overlap_phi) = recv_phi_xMzM(i,j,k)
    !                phi(i+nx,j,k-overlap_phi) = recv_phi_xPzM(i,j,k)
    !                phi(i-overlap_phi,j,k+nz) = recv_phi_xMzP(i,j,k)
    !                phi(i+nx,j,k+nz) = recv_phi_xPzP(i,j,k)
    !            enddo
    !        enddo
    !    enddo
    !    !$acc end kernels
    !#endif
    !
    !#if defined(x_mpi)&&defined(y_mpi)
    !    CALL MPI_WAITALL(8, MPI_REQ_phi_EZ, MPI_ESTAT, MPI_ERR)
    !    !$acc update device(recv_phi_xPyP,recv_phi_xPyM,recv_phi_xMyP,recv_phi_xMyM)async(z_edge_phi_recv_unpack_async)
    !    !$acc kernels present(phi,recv_phi_xPyP,recv_phi_xPyM,recv_phi_xMyP,recv_phi_xMyM)async(z_edge_phi_recv_unpack_async)
    !    !$acc loop collapse(3)
    !    do k = 1,nz
    !        do j=1,overlap_phi
    !            do i = 1,overlap_phi
    !                phi(i-overlap_phi,j-overlap_phi,k) = recv_phi_xMyM(i,j,k)
    !                phi(i+nx,j-overlap_phi,k) = recv_phi_xPyM(i,j,k)
    !                phi(i-overlap_phi,j+ny,k) = recv_phi_xMyP(i,j,k)
    !                phi(i+nx,j+ny,k) = recv_phi_xPyP(i,j,k)
    !            enddo
    !        enddo
    !    enddo
    !    !$acc end kernels
    !#endif

    !$OMP END PARALLEL

    return
end subroutine MPI_phi_sync_update







