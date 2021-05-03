!==============================================================================
!   This a sample code to create simple geometries to be used in MF-LBM
!==============================================================================

module wall
    implicit none
    integer :: nx,ny,nz
    integer(kind=1), allocatable, dimension(:,:,:) :: walls    ! kind=1 to reduce file size
end module wall

!**********************main***************************
program Main
    use wall
    integer :: i,j,k, inlet_buffer, outlet_buffer
    double precision :: r1, r2, xc, yc, zc

    ! domain dimension
    nx = 60
    ny = 60
    nz = 80   ! flow direction

    ! create inlet/outlet reservoirs
    inlet_buffer = 5
    outlet_buffer = 5

    allocate(walls(1:nx,1:ny,1:nz))

    walls=0
    
    xc = 0.5d0 * (nx + 1)
    yc = 0.5d0 * (ny + 1)
    zc = 0.5d0 * (nz + 1)
    r1 = xc
    r2 = 0.3d0 * xc

    ! change the square duct to a tube and place a sphere obstacle in the center
    do k=1,nz
        do j=1,ny
            do i=1,nx
                if(( i - xc )**2+( j - yc )**2 > r1**2)then  ! tube
                    walls(i,j,k)=1
                endif
                if(( i - xc )**2 + ( j - yc )**2 + ( k - zc )**2 < r2**2)then   !sphere
                    walls(i,j,k)=1
                endif
                !create inlet/outlet reservoirs
                if( k<=inlet_buffer .or. k>=nz-outlet_buffer+1)then
                    walls(i,j,k)=0
                endif
            enddo
        enddo
    enddo

    !apply channel walls
    walls(1,:,:)=1
    walls(nx,:,:)=1
    walls(:,1,:)=1
    walls(:,ny,:)=1

    call VTK_walls_bin

    call write_walls('walls.dat')

    deallocate(walls)

end program Main
!**********************main***************************

!**********************write walls***************************
subroutine write_walls(filename)
    use wall
    integer :: i,j,k
    character(len=*) :: filename

    open(11, file=filename, status='replace',form='unformatted',access='stream')    ! default geometry file format used in MF-LBM
    write(11)nx,ny,nz
    write(11)(((walls(i,j,k),i=1,nx),j=1,ny),k=1,nz)   
    !!
    close(11)

    return
end subroutine write_walls
!**********************write walls***************************

!*******************************************visualization**************************************************
subroutine VTK_walls_bin   !solid geometry
    use wall
    IMPLICIT NONE
    integer :: i,j,k
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10
    integer   :: ivtk = 9, int

    open(unit=ivtk,file='./walls_bin.vtk',FORM='unformatted',access='stream',status='replace',convert='BIG_ENDIAN')

    lf = char(10) ! line feed character

    buffer = '# vtk DataFile Version 3.0'//lf                            ; write(ivtk) trim(buffer)
    buffer = 'vtk output'//lf                                             ; write(ivtk) trim(buffer)
    buffer = 'BINARY'//lf                                                 ; write(ivtk) trim(buffer)
    buffer = 'DATASET STRUCTURED_POINTS '//lf                              ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')nx
    write(str2(1:10),'(i10)')ny
    write(str3(1:10),'(i10)')nz
    buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf                 ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf                         ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'SPACING '//str1//' '//str2//' '//str3//lf                          ; write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')nx*ny*nz
    buffer = 'POINT_DATA '//str1//lf                                               ; write(ivtk) trim(buffer)

    !scalar - walls
    buffer = 'SCALARS walls int'//lf                                          ; write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf                                          ; write(ivtk) trim(buffer)
    write(ivtk)(((int(walls(i,j,k)),i=1,nx),j=1,ny),k=1,nz)

    close(ivtk)

    return
end subroutine VTK_walls_bin 
!*******************************************visualization**************************************************
