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
    endif

    !faces y
    if(mpi_y)then
        CALL MPI_IRECV( recv_pdf_yP, isize_pdf_y, MPI_DOUBLE_PRECISION, idyP, tag3, MPI_COMM_VGRID, MPI_REQ_y(1), MPI_ERR )
        CALL MPI_IRECV( recv_pdf_yM, isize_pdf_y, MPI_DOUBLE_PRECISION, idyM, tag4, MPI_COMM_VGRID, MPI_REQ_y(2), MPI_ERR )
    endif 


    !edges x
    if(mpi_y.and.mpi_z)then  
        CALL MPI_IRECV( recv_pdf_yPzP, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyPzP, TAG7, MPI_COMM_VGRID, MPI_REQ_EX(1), MPI_ERR )
        CALL MPI_IRECV( recv_pdf_yPzM, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyPzM, TAG8, MPI_COMM_VGRID, MPI_REQ_EX(2), MPI_ERR )
        CALL MPI_IRECV( recv_pdf_yMzP, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyMzP, TAG9, MPI_COMM_VGRID, MPI_REQ_EX(3), MPI_ERR )
        CALL MPI_IRECV( recv_pdf_yMzM, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyMzM, TAG10, MPI_COMM_VGRID, MPI_REQ_EX(4), MPI_ERR )
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
    endif

    !face y
    if(mpi_y)then
        !$acc wait(y_pdf_update_async)
        CALL MPI_Isend( send_pdf_yM, isize_pdf_y, MPI_DOUBLE_PRECISION, idyM, tag3, MPI_COMM_VGRID, MPI_REQ_y(3), MPI_ERR )
        CALL MPI_Isend( send_pdf_yP, isize_pdf_y, MPI_DOUBLE_PRECISION, idyP, tag4, MPI_COMM_VGRID, MPI_REQ_y(4), MPI_ERR )
    endif

    !$acc wait(edge_pdf_update_async)

    !edge x
    if(mpi_y.and.mpi_z)then 
        CALL MPI_Isend( send_pdf_yMzM, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyMzM, TAG7, MPI_COMM_VGRID, MPI_REQ_EX(5), MPI_ERR )
        CALL MPI_Isend( send_pdf_yMzP, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyMzP, TAG8, MPI_COMM_VGRID, MPI_REQ_EX(6), MPI_ERR )
        CALL MPI_Isend( send_pdf_yPzM, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyPzM, TAG9, MPI_COMM_VGRID, MPI_REQ_EX(7), MPI_ERR )
        CALL MPI_Isend( send_pdf_yPzP, isize_pdf_ex, MPI_DOUBLE_PRECISION,idyPzP, TAG10, MPI_COMM_VGRID, MPI_REQ_EX(8), MPI_ERR )
    endif 

    return
end subroutine mpi_send_req


!=====================================================================================================================================
!------- pdf data packing (pull type) --------
!=====================================================================================================================================
subroutine mpi_pdf_pull_async_pack
    use Misc_module
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
        !$acc & send_pdf_zP,send_pdf_zM)async(LBM_async_z)
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
            enddo
        enddo
        !$acc end kernels
        !$acc update host(send_pdf_zP,send_pdf_zM)async(z_pdf_update_async)wait(LBM_async_z)
    endif

    if(mpi_y)then
        !****************************y direction******************************************
        !$OMP DO private(i)
        !the reason that acc_wait async_z is due to y halo calculation is not complete calculation. z-halo caculation is complete halo layer calulation
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,send_pdf_yP,send_pdf_yM)async(LBM_async_y)wait(LBM_async_z)
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
            enddo
        enddo
        !$acc end kernels
        !$acc update host(send_pdf_yP,send_pdf_yM)async(y_pdf_update_async)wait(LBM_async_y)
    endif

    !****************************edges******************************************
    if(mpi_y.and.mpi_z)then
        ! Edges along X
        !$OMP DO
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc & send_pdf_yPzP,send_pdf_yPzM,send_pdf_yMzP,send_pdf_yMzM)async(edge_pdf_async)wait(LBM_async_z,LBM_async_y)
        !$acc loop
        DO i = 1, NX
            send_pdf_yPzP(i) = f15(i,ny,nz)
            send_pdf_yPzM(i) = f17(i,ny,1 )
            send_pdf_yMzP(i) = f16(i,1 ,nz)
            send_pdf_yMzM(i) = f18(i,1 ,1 )
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
        !$acc & ,send_pdf_zP,send_pdf_zM)async(z_pdf_update_async)
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
        !$acc & ,send_pdf_yP,send_pdf_yM)async(y_pdf_update_async)
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
        !$acc & recv_pdf_yPzP,recv_pdf_yPzM,recv_pdf_yMzP,recv_pdf_yMzM)wait(z_pdf_update_async,y_pdf_update_async)
        !$acc loop
        DO i = 1, NX
            f18(i,ny+1,nz+1) = recv_pdf_yPzP(i)
            f16(i,ny+1,1-1 ) = recv_pdf_yPzM(i)
            f17(i,1-1 ,nz+1) = recv_pdf_yMzP(i)
            f15(i,1-1 ,1-1 ) = recv_pdf_yMzM(i)
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
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,send_pdf_zP,send_pdf_zM)async(LBM_async_z)
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
            enddo
        enddo
        !$acc end kernels
        !$acc update host(send_pdf_zP,send_pdf_zM)async(z_pdf_update_async)wait(LBM_async_z)
    endif

    if(mpi_y)then
        !****************************y direction******************************************
        !$OMP DO private(i)
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,send_pdf_yP,send_pdf_yM)async(LBM_async_y)wait(LBM_async_z)
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
            enddo
        enddo
        !$acc end kernels
        !$acc update host(send_pdf_yP,send_pdf_yM)async(y_pdf_update_async)wait(LBM_async_y)
    endif
    

    !****************************edges******************************************

    if(mpi_y.and.mpi_z)then
        ! Edges along X
        !$OMP DO 
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc & send_pdf_yPzP,send_pdf_yPzM,send_pdf_yMzP,send_pdf_yMzM)async(edge_pdf_async)wait(LBM_async_z,LBM_async_y)
        !$acc loop
        DO i = 1, NX
            send_pdf_yPzP(i) = f18(i,ny+1,nz+1)
            send_pdf_yPzM(i) = f16(i,ny+1,1 -1)
            send_pdf_yMzP(i) = f17(i,1 -1,nz+1)
            send_pdf_yMzM(i) = f15(i,1 -1,1 -1)
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
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,send_pdf_zP,send_pdf_zM)async(z_pdf_update_async)
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
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,send_pdf_yP,send_pdf_yM)async(y_pdf_update_async)
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
        !$acc kernels present(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,&
        !$acc & recv_pdf_yPzP,recv_pdf_yPzM,recv_pdf_yMzP,recv_pdf_yMzM)async(edge_pdf_update_async)wait(z_pdf_update_async,y_pdf_update_async)
        !$acc loop
        DO i = 1, NX
            f15(i,ny,nz) = recv_pdf_yPzP(i)
            f17(i,ny,1 ) = recv_pdf_yPzM(i)
            f16(i,1 ,nz) = recv_pdf_yMzP(i)
            f18(i,1 ,1 ) = recv_pdf_yMzM(i)
        END DO
        !$acc end kernels
    endif

    !$OMP END PARALLEL

    return
end subroutine mpi_pdf_push_sync_update









