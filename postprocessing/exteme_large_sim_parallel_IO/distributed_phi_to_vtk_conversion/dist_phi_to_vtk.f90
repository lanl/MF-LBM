
!==============================================================================
! This code converts distributed phase filed data (save_phi function in MF-LBM)
! to vtk files for easy visualization. save_phi is activated when
! extreme_large_sim_cmd = 1 in control file (99% of chance that you don't need this)
!==============================================================================

MODULE misc
    IMPLICIT NONE
    SAVE
    real(kind=4),ALLOCATABLE,DIMENSION(:,:,:) :: phi
    integer :: np
    integer :: ntime0,ntime_interval,ntime_max
    integer :: nxglobal,nyglobal,nzglobal
END MODULE misc


!################################################################################################################################
!Main Start
!################################################################################################################################
program main_multiphase
    use misc
    IMPLICIT NONE
    integer :: i,j,k
    double precision :: temp
    character (len=50) :: flnm   !file name

    call system('mkdir ./output')

    open(5,file='./input_parameters.txt',status='old')
    read(5,*)np  ! total number of MPI processes in the original simulation
    read(5,*)nxglobal,nyglobal,nzglobal  !dimensions of the entire simulation domain
    read(5,*)ntime0,ntime_interval,ntime_max  !start time step, time interval, end time step
    close(5)

    allocate(phi(nxglobal,nyglobal,nzglobal))
    !$OMP parallel DO private(i,j)
    do k=1,nzglobal
        do j=1,nyglobal
            do i=1,nxglobal
                phi(i,j,k)=0.0
            enddo
        enddo
    enddo

    call read_save_phi

    deallocate(phi)

End program main_multiphase
!################################################################################################################################
!Main End
!################################################################################################################################

!************************************* read and save macro field files *********************************
subroutine read_save_phi
    use Misc
    implicit none
    integer :: i,j,k,l,m,n, nx,ny,nz
    integer :: n_step,ncount,num, nt,n_fluid_node_local, id,idx,idy,idz
    character (len=30) :: flnm   !file name
    real (kind=4):: v1
    type geometry
    integer :: n_fluid_node_local
        integer, allocatable, dimension(:) :: ix,iy,iz
    end type geometry

    type(geometry), allocatable :: fluid_array(:)

    allocate(fluid_array(0:np-1))

    do id=0,np-1
        write(flnm,"('./data/geometry_id',i5.5)")id
        open(unit=9+id, file=trim(flnm), FORM='unformatted', status='old',access='stream')
        read(9+id)fluid_array(id)%n_fluid_node_local  !total number of fluid nodes of local MPI subdomain
        allocate(fluid_array(id)%ix(fluid_array(id)%n_fluid_node_local),fluid_array(id)%iy(fluid_array(id)%n_fluid_node_local),fluid_array(id)%iz(fluid_array(id)%n_fluid_node_local))
        read(9+id)idx,idy,idz,nx,ny,nz
        do num=1,fluid_array(id)%n_fluid_node_local
            read(9+id)i,j,k
            l = idx*nx + i
            m = idy*ny + j
            n = idz*nz + k
            fluid_array(id)%ix(num) = l
            fluid_array(id)%iy(num) = m
            fluid_array(id)%iz(num) = n
        enddo
        close(9+id)
    enddo

    ncount = (ntime_max-ntime0)/ntime_interval + 1
    print*,'number of vtk files pending processing:', ncount

    do n_step=1,ncount
        nt = ntime0 + (n_step-1)*ntime_interval
        print*, 'Start processing phase field of ntime=', nt
        !$omp parallel do private(flnm,num,v1,i,j,k)
        do id=0,np-1
            write(flnm,"('phi_nt',i9.9,'_id',i5.5)")nt,id
            open(unit=9+id, file='./data/'//trim(flnm), FORM='unformatted', status='old',access='stream')
            do num=1,fluid_array(id)%n_fluid_node_local
                read(9+id)v1
                i = fluid_array(id)%ix(num)
                j = fluid_array(id)%iy(num)
                k = fluid_array(id)%iz(num)
                phi(i,j,k) = v1
            enddo
            close(9+id)
        enddo

        call VTK_small_bin(nt)
        call VTK_small_bin_half(nt)  
        print*, 'End processing phase field of ntime=', nt
    enddo

    deallocate(fluid_array)

    return
end subroutine read_save_phi
!************************************* read and save macro field files *********************************


!************************************* VTK legacy output **************************************
subroutine VTK_small_bin(nt)   !phase field, single precision
    use Misc
    IMPLICIT NONE
    character*30 :: flnm
    integer :: i,j,k, nt
    integer(kind=8) :: num
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10
    integer   :: ivtk = 9, int

    write(flnm,'(i10.10,".vtk")')nt

    open(unit=ivtk,file='./output/small_'//flnm,FORM='unformatted',access='stream',status='replace',convert='BIG_ENDIAN')

    lf = char(10) ! line feed character

    buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
    buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
    buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
    buffer = 'DATASET STRUCTURED_POINTS '//lf                                          ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')nxglobal
    write(str2(1:10),'(i10)')nyglobal
    write(str3(1:10),'(i10)')nzglobal
    buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'SPACING '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

    num=nxglobal*nyglobal*nzglobal

    write(str1(1:10),'(i10)')num
    buffer = 'POINT_DATA '//str1//lf                                               ; write(ivtk) trim(buffer)

    !scalar - phase field
    buffer = 'SCALARS phi float'//lf                                          ; write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf                                          ; write(ivtk) trim(buffer)
    write(ivtk)(((phi(i,j,k),i=1,nxglobal),j=1,nyglobal),k=1,nzglobal)

    close(ivtk)

    return
end

subroutine  VTK_small_bin_half(nt)   !solid geometry
    use Misc
    IMPLICIT NONE
    character*30 :: flnm
    integer :: i,j,k, nt
    integer(kind=8) :: num
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10
    integer   :: ivtk = 9, int

    write(flnm,'(i10.10,".vtk")')nt

    open(unit=ivtk,file='./output/small_half_'//flnm,FORM='unformatted',access='stream',status='replace',convert='BIG_ENDIAN')

    lf = char(10) ! line feed character

    buffer = '# vtk DataFile Version 3.0'//lf
    write(ivtk) trim(buffer)
    buffer = 'vtk output'//lf
    write(ivtk) trim(buffer)
    buffer = 'BINARY'//lf
    write(ivtk) trim(buffer)
    buffer = 'DATASET STRUCTURED_POINTS '//lf
    write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')nxglobal/2
    write(str2(1:10),'(i10)')nyglobal/2
    write(str3(1:10),'(i10)')nzglobal/2
    buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf
    write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf
    write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'SPACING '//str1//' '//str2//' '//str3//lf
    write(ivtk) trim(buffer)

    num = int(nxglobal,kind=8)*int(nyglobal,kind=8)*int(nzglobal,kind=8)/8
    write(str1(1:10),'(i10)')num
    buffer = 'POINT_DATA '//str1//lf
    write(ivtk) trim(buffer)

    !scalar - walls
    buffer = 'SCALARS phi float'//lf                                          ;
    write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write(ivtk) trim(buffer)
    write(ivtk)(((phi(i,j,k),i=1,nxglobal,2),j=1,nyglobal,2),k=1,nzglobal,2)

    close(ivtk)

    return
end

!************************************* VTK legacy output **************************************