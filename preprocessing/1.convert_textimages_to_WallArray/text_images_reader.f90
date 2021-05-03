!==============================================================================
! This a sample code to convert text images to wall array and stored in binary 
! format. The modification program may be needed to add buffer layers and cropping.
!==============================================================================

module wall
    implicit none
    integer :: nx,ny,nz
    character(len=50) :: images_name
    character(len=300) :: data_location
    integer :: images_type, images_digit, images_start_number
    integer, allocatable, dimension(:,:,:) :: walls
    integer (kind=1), allocatable, dimension(:,:,:) :: walls_int1
    double precision :: porosity
end module wall

!**********************main***************************
program Main
    use wall
    implicit none
    integer :: i,j,k
    integer (kind=8) :: icount, itemp
    character(len=100) :: dummy

    call system('mkdir ./output')

    open(5,file='input_parameters.txt',status='old')
    read(5,'(A)') dummy
    read(5,'(A)') data_location
    read(5,'(A)') dummy
    read(5,'(A)') images_name
    read(5,*) images_digit   ! number of digit of the text image serial
    read(5,*) images_start_number ! starting slice number 
    read(5,*) images_type    ! type1 - solid node value <=127; type2 (2) - solid node value >127
    read(5,*) nx,ny,nz   !  dimensions
    close(5)

    print*,'original grid size: ', nx,ny,nz   ! dimensions

    allocate(walls(1:nx,1:ny,1:nz), walls_int1(1:nx,1:ny,1:nz))
    walls=0
    walls_int1=0

    call read_walls_text

    if(images_type==2)then  ! type2 - solid node value >127
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,i) 
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(walls(i,j,k)>127)then
                        walls_int1(i,j,k) = 1
                    endif            
                enddo
            enddo
        enddo 
    elseif(images_type==1)then   ! type1 - solid node value <=127
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,i) 
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if(walls(i,j,k)<=127)then
                        walls_int1(i,j,k) = 1
                    endif            
                enddo
            enddo
        enddo 
    else
        print*, 'input parameter images_type error!'
        stop
    endif

    icount = 0
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,i) reduction(+:icount)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                if(walls_int1(i,j,k)==0)then
                    icount = icount + 1
                endif            
            enddo
        enddo
    enddo 
    itemp =  int(nx,kind=8)*int(ny,kind=8)*int(nz,kind=8)      
    porosity = dble(icount)/dble(itemp)
    print*,'porosity is', porosity

    !apply channel walls
    walls_int1(1,:,:)=1
    walls_int1(nx,:,:)=1
    walls_int1(:,1,:)=1
    walls_int1(:,ny,:)=1
    
    call visual
    !call visual_skip_half     ! if the geometry is tooooooo large

    open(11, file='./output/walls_processed.dat', status='replace',form='unformatted',access='stream')   ! format used in the main simulation code
    write(11)nx,ny,nz
    write(11)(((walls_int1(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    close(11)

    deallocate(walls,walls_int1)

end program Main
!**********************main***************************

!**********************read walls**************************
subroutine read_walls_text
    use wall
    integer :: ix, iy, iz, iz_abs
    character(len=50) :: str

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iz,iz_abs,iy,ix,str)
    !$OMP DO
    do iz = 1, nz
        iz_abs = iz + images_start_number - 1 
        write(str,*) iz_abs 
        select case(images_digit)
            case (1)
                open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//trim(adjustl(str))//'.txt',status='unknown')
            case (2)
                if (iz_abs<=9) then
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//'0'//trim(adjustl(str))//'.txt',status='unknown')
                else
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//trim(adjustl(str))//'.txt',status='unknown')
                end if
            case (3)
                if (iz_abs<=9) then
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//'00'//trim(adjustl(str))//'.txt',status='unknown')
                else if (iz_abs<=99) then
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//'0'//trim(adjustl(str))//'.txt',status='unknown')
                else
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//trim(adjustl(str))//'.txt',status='unknown')
                end if
            case (4)
                if (iz_abs<=9) then
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//'000'//trim(adjustl(str))//'.txt',status='unknown')
                else if (iz_abs<=99) then
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//'00'//trim(adjustl(str))//'.txt',status='unknown')
                else if (iz_abs<=999) then
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//'0'//trim(adjustl(str))//'.txt',status='unknown')
                else
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//trim(adjustl(str))//'.txt',status='unknown')
                end if
            case (5)
                if (iz_abs<=9) then
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//'0000'//trim(adjustl(str))//'.txt',status='unknown')
                else if (iz_abs<=99) then
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//'000'//trim(adjustl(str))//'.txt',status='unknown')
                else if (iz_abs<=999) then
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//'00'//trim(adjustl(str))//'.txt',status='unknown')
                else if (iz_abs<=9999) then
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//'0'//trim(adjustl(str))//'.txt',status='unknown')
                else
                    open(unit=10+iz_abs, file=trim(adjustl(data_location))//'/'//trim(adjustl(images_name))//trim(adjustl(str))//'.txt',status='unknown')
                end if
        end select
        do iy = 1,ny
            read(10+iz_abs,*) (walls(ix,iy,iz), ix= 1,nx)
        end do
        close(10+iz_abs)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    return
end subroutine read_walls_text
!**********************read walls***************************



!*******************************************visualization**************************************************    
subroutine visual   !solid geometry
    use wall
    IMPLICIT NONE
    
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

    write(str1(1:10),'(i10)')nx*ny*nz
    buffer = 'POINT_DATA '//str1//lf                                               ; write(ivtk) trim(buffer)

    !scalar - walls
    buffer = 'SCALARS walls int'//lf                                          ; write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf                                          ; write(ivtk) trim(buffer)
    write(ivtk)(((int(walls_int1(i,j,k)),i=1,nx),j=1,ny),k=1,nz)

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
    write(ivtk)(((int(walls_int1(i,j,k)),i=1,nx,2),j=1,ny,2),k=1,nz,2)

    close(ivtk)

    return
end    
!*******************************************visualization**************************************************
