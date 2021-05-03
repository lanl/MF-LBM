module wall
    implicit none
    integer :: nx,ny,nz,nx0,ny0,nz0,inlet_buffer,outlet_buffer,mirror_option,crop_option  
    double precision :: xc,yc,zc
    integer (kind=1), allocatable, dimension(:,:,:) ::  walls_raw, walls,walls_tmp
    double precision :: porosity
    character(len=300) :: wall_file_path
end module wall

!**********************main***************************
program Main
    use wall
    implicit none
    integer :: i,j,k
    integer (kind=8) :: icount, itemp
    character(len=300) :: dummy

    call system('mkdir ./output')

    open(5,file='./input_parameter.txt',status='old')
    read(5,'(A)') dummy
    read(5,'(A)') wall_file_path
    read(5,*) inlet_buffer,outlet_buffer    ! add inlet and outlet buffer layers to implement inlet/outlet boundary conditions
    read(5,*) crop_option  
    read(5,*) xc,yc,zc   !crop center coordinate
    read(5,*) nx,ny,nz   !cropped domain size (not including buffer layers)
    close(5)
 
    OPEN(UNIT=9,FILE=trim(adjustl(wall_file_path)),STATUS='OLD',FORM='UNFORMATTED',access='stream')
    read(9)nx0,ny0,nz0
    print*,'original porous media sample size: ', 'nx0=', nx0, 'ny0=', ny0,'nz0=', nz0
    allocate(walls_raw(1:nx0,1:ny0,1:nz0)) 
    read(9)(((walls_raw(i,j,k),i=1,nx0),j=1,ny0),k=1,nz0)
    close(9)


    if(crop_option==0)then
        nx = nx0
        ny = ny0
        nz = nz0 + inlet_buffer + outlet_buffer
        allocate(walls(1:nx,1:ny,1:nz))
        walls = 0
        icount = 0  
        do k=1+inlet_buffer,nz-outlet_buffer
            do j=1,ny
                do i=1,nx
                    walls(i,j,k) = walls_raw(i,j,k-inlet_buffer)
                enddo
            enddo
        enddo  
    else
        allocate(walls(1:nx,1:ny,1:nz+inlet_buffer+outlet_buffer))
        walls=0 
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(walls_raw(int(xc-dble(1+nx)/2d0 +0.5d0)+i,int(yc-dble(1+ny)/2d0 +0.5d0)+j,int(zc-dble(1+nz)/2d0 +0.5d0)+k)==1)then
                        walls(i,j,k+inlet_buffer) = 1
                    endif            
                enddo
            enddo
        enddo
        nz = nz +inlet_buffer + outlet_buffer  !final grid nodes
    endif

    !apply channel walls
    walls(1,:,:)=1
    walls(nx,:,:)=1
    walls(:,1,:)=1
    walls(:,ny,:)=1

    print*, 'final grid: nx= ', nx, ' ny= ', ny, ' nx= ', nz   

    call visual
    !call visual_skip_half

    open(11, file='./output/walls_processed.dat', status='replace',form='unformatted',access='stream')
    write(11)nx,ny,nz
    write(11)(((walls(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    close(11)

    deallocate(walls,walls_raw)

end program Main
!**********************main***************************


!*******************************************visualization**************************************************    
subroutine visual   !solid geometry
    use wall
    IMPLICIT NONE
    integer (kind=8) :: num
    integer :: i,j,k
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10
    integer   :: ivtk = 9, int

    open(unit=ivtk,file='./output/geometry.vtk',FORM='unformatted',access='stream',status='replace',convert='BIG_ENDIAN')

    lf = char(10) ! line feed character

    buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
    buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
    buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
    buffer = 'DATASET STRUCTURED_POINTS '//lf                                          ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')nx
    write(str2(1:10),'(i10)')ny
    write(str3(1:10),'(i10)')nz
    buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'SPACING '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

    num = int(nx,kind=8)*int(ny,kind=8)*int(nz,kind=8)
    write(str1(1:10),'(i10)')num
    buffer = 'POINT_DATA '//str1//lf                                               ; write(ivtk) trim(buffer)

    !scalar - walls
    buffer = 'SCALARS walls int'//lf                                          ; write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf                                          ; write(ivtk) trim(buffer)
    write(ivtk)(((int(walls(i,j,k)),i=1,nx),j=1,ny),k=1,nz)

    close(ivtk)

    return
end

subroutine visual_skip_half   !solid geometry
    use wall
    IMPLICIT NONE
    
    integer :: i,j,k
    integer (kind=8) :: num
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10
    integer   :: ivtk = 9, int

    open(unit=ivtk,file='./output/geometry_skip.vtk',FORM='unformatted',access='stream',status='replace',convert='BIG_ENDIAN')

    lf = char(10) ! line feed character

    buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
    buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
    buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
    buffer = 'DATASET STRUCTURED_POINTS '//lf                                          ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')nx/2
    write(str2(1:10),'(i10)')ny/2
    write(str3(1:10),'(i10)')nz/2
    buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'SPACING '//str1//' '//str2//' '//str3//lf                                               ; write(ivtk) trim(buffer)

    num = int(nx,kind=8)*int(ny,kind=8)*int(nz,kind=8)/8
    write(str1(1:10),'(i10)')num
    buffer = 'POINT_DATA '//str1//lf                                               ; write(ivtk) trim(buffer)

    !scalar - walls
    buffer = 'SCALARS walls int'//lf                                          ; write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf                                          ; write(ivtk) trim(buffer)
    write(ivtk)(((int(walls(i,j,k)),i=1,nx,2),j=1,ny,2),k=1,nz,2)

    close(ivtk)

    return
end   

    !*******************************************visualization**************************************************
