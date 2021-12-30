!==============================================================================
! This code converts distributed macro data (save_macro function in MF-LBM)
! to full domain field data for easy analysis. save_macro is activated when
! extreme_large_sim_cmd = 1 in control file (99% of chance that you don't need this)
!==============================================================================

MODULE misc
    IMPLICIT NONE
    SAVE
    real(kind=8), ALLOCATABLE, DIMENSION(:, :, :) :: u, v, w, rho
    real(kind=4), ALLOCATABLE, DIMENSION(:, :, :) :: usp, vsp, wsp, rhosp     !single precision version
    integer :: np, input_precision_choice, output_precision_choice, skip_point_choice
    integer :: ntime0, ntime_interval, ntime_max
    integer :: nxglobal, nyglobal, nzglobal
    character(len=300) :: data_location
END MODULE misc

!################################################################################################################################
!Main Start
!################################################################################################################################
program main_singlephase
    use misc
    IMPLICIT NONE
    integer :: i, j, k
    double precision :: temp
    character(len=300) :: flnm, syscmd, dummy   !file name

    open (5, file='input_parameters.txt', status='old')
    read (5, '(A)') dummy
    read (5, '(A)') data_location
    read (5, *) np  ! total number of MPI processes in the original simulation
    read (5, *) nxglobal, nyglobal, nzglobal  !dimensions of the entire simulation domain
    read (5, *) ntime0, ntime_interval, ntime_max  !start time step, time interval, end time step
    read (5, *) input_precision_choice
    read (5, *) output_precision_choice
    read (5, *) skip_point_choice
    close (5)

 allocate(u(nxglobal,nyglobal,nzglobal),v(nxglobal,nyglobal,nzglobal),w(nxglobal,nyglobal,nzglobal),rho(nxglobal,nyglobal,nzglobal))
    !$omp parallel DO private(i,j)
    do k = 1, nzglobal
        do j = 1, nyglobal
            do i = 1, nxglobal
                u(i, j, k) = 0d0
                v(i, j, k) = 0d0
                w(i, j, k) = 0d0
                rho(i, j, k) = 1d0
            end do
        end do
    end do

    if (input_precision_choice == 1) then
      allocate(usp(nxglobal,nyglobal,nzglobal),vsp(nxglobal,nyglobal,nzglobal),wsp(nxglobal,nyglobal,nzglobal),rhosp(nxglobal,nyglobal,nzglobal))
        !$omp parallel DO private(i,j)
        do k = 1, nzglobal
            do j = 1, nyglobal
                do i = 1, nxglobal
                    usp(i, j, k) = 0d0
                    vsp(i, j, k) = 0d0
                    wsp(i, j, k) = 0d0
                    rhosp(i, j, k) = 1d0
                end do
            end do
        end do
    elseif (input_precision_choice /= 2) then
        print *, 'Precision choice error! Aborting...'
        stop
    end if

    write (syscmd, "(A)") 'mkdir '//trim(data_location)//'/processed_field_data'
    call system(trim(syscmd))

    call read_save_macro

    deallocate (u, v, w, rho)
    if (input_precision_choice == 1) then
        deallocate (usp, vsp, wsp, rhosp)
    end if

End program main_singlephase
!################################################################################################################################
!Main End
!################################################################################################################################

!************************************* read and save macro field files *********************************
subroutine read_save_macro
    use Misc
    implicit none
    integer :: i, j, k, l, m, n, nx, ny, nz
    integer :: n_step, ncount, num, nt, id, idx, idy, idz
    character(len=30) :: flnm   !file name
    double precision :: v1, v2, v3, v4, v5

    ncount = (ntime_max - ntime0)/ntime_interval + 1
    print *, 'number of vtk files pending processing:', ncount

    do n_step = 1, ncount
        nt = ntime0 + (n_step - 1)*ntime_interval
        print *, 'Start processing macro field of ntime=', nt
        !$omp parallel do private(flnm,i,j,k,idx,idy,idz,nx,ny,nz)
        do id = 0, np - 1
            write (flnm, "('full_nt',i9.9,'_id',i5.5)") nt, id
            open (unit=9 + id, file=trim(data_location)//'/'//trim(flnm), FORM='unformatted', status='old', access='stream')
            if (input_precision_choice == 1) then
                read (9 + id) idx, idy, idz, nx, ny, nz
                read (9 + id) (((usp(i + idx*nx, j + idy*ny, k + idz*nz), i=1, nx), j=1, ny), k=1, nz)
                read (9 + id) (((vsp(i + idx*nx, j + idy*ny, k + idz*nz), i=1, nx), j=1, ny), k=1, nz)
                read (9 + id) (((wsp(i + idx*nx, j + idy*ny, k + idz*nz), i=1, nx), j=1, ny), k=1, nz)
                read (9 + id) (((rhosp(i + idx*nx, j + idy*ny, k + idz*nz), i=1, nx), j=1, ny), k=1, nz)
            else
                read (9 + id) idx, idy, idz, nx, ny, nz
                read (9 + id) (((u(i + idx*nx, j + idy*ny, k + idz*nz), i=1, nx), j=1, ny), k=1, nz)
                read (9 + id) (((v(i + idx*nx, j + idy*ny, k + idz*nz), i=1, nx), j=1, ny), k=1, nz)
                read (9 + id) (((w(i + idx*nx, j + idy*ny, k + idz*nz), i=1, nx), j=1, ny), k=1, nz)
                read (9 + id) (((rho(i + idx*nx, j + idy*ny, k + idz*nz), i=1, nx), j=1, ny), k=1, nz)
            end if
            close (9 + id)
        end do

        !call save_macro(nt)   
        if (input_precision_choice == 1) then
            !$omp parallel DO private(i,j)
            do k = 1, nzglobal
                do j = 1, nyglobal
                    do i = 1, nxglobal
                        u(i, j, k) = usp(i, j, k)
                        v(i, j, k) = vsp(i, j, k)
                        w(i, j, k) = wsp(i, j, k)
                        rho(i, j, k) = rhosp(i, j, k)
                    end do
                end do
            end do
        end if
        if (skip_point_choice == 1) then
            call VTK_detail_bin_half_sp(nt)  !single precision only
        else
            call VTK_detail_bin(nt)
        end if

        print *, 'End processing macro field of ntime=', nt
    end do

    return
end subroutine read_save_macro
!************************************* read and save macro field files *********************************

!************************************* write data *********************************
subroutine save_macro(nt)
    use Misc
    IMPLICIT NONE
    character*30 :: flnm
    integer :: nt

    if (input_precision_choice == 1) then
        write (flnm, '(i10.10,".dat")') nt
  open (unit=9, file=trim(data_location)//'/processed_field_data/rho_ntime_'//trim(flnm), FORM='unformatted', access='stream', status='replace', convert='BIG_ENDIAN')
        write (9) rhosp
        close (9)
  open (unit=9, file=trim(data_location)//'/processed_field_data/vel_ntime_'//trim(flnm), FORM='unformatted', access='stream', status='replace', convert='BIG_ENDIAN')
        write (9) usp
        write (9) vsp
        write (9) wsp
        close (9)
    else
        write (flnm, '(i10.10,".dat")') nt
  open (unit=9, file=trim(data_location)//'/processed_field_data/rho_ntime_'//trim(flnm), FORM='unformatted', access='stream', status='replace', convert='BIG_ENDIAN')
        write (9) rho
        close (9)
  open (unit=9, file=trim(data_location)//'/processed_field_data/vel_ntime_'//trim(flnm), FORM='unformatted', access='stream', status='replace', convert='BIG_ENDIAN')
        write (9) u
        write (9) v
        write (9) w
        close (9)
    end if

    return
end subroutine save_macro

!************************************* VTK legacy output **************************************
subroutine VTK_detail_bin(nt)
    use Misc
    IMPLICIT NONE
    character*30 :: flnm, fmt
    integer :: i, j, k, nt
    integer(kind=8) :: num
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10, str4*14
    integer   :: ivtk = 9, int

    write (flnm, '(i10.10,".vtk")') nt

    open (unit=ivtk, file=trim(data_location)//'/processed_field_data/detail_'//flnm, FORM='unformatted', access='stream', status='replace', convert='BIG_ENDIAN')

    lf = char(10) ! line feed character

    buffer = '# vtk DataFile Version 3.0'//lf; write (ivtk) trim(buffer)
    buffer = 'vtk output'//lf; write (ivtk) trim(buffer)
    buffer = 'BINARY'//lf; write (ivtk) trim(buffer)
    buffer = 'DATASET STRUCTURED_POINTS '//lf; write (ivtk) trim(buffer)

    write (str1(1:10), '(i10)') nxglobal
    write (str2(1:10), '(i10)') nyglobal
    write (str3(1:10), '(i10)') nzglobal
    buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf; write (ivtk) trim(buffer)

    write (str1(1:10), '(i10)') 1
    write (str2(1:10), '(i10)') 1
    write (str3(1:10), '(i10)') 1
    buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf; write (ivtk) trim(buffer)

    write (str1(1:10), '(i10)') 1
    write (str2(1:10), '(i10)') 1
    write (str3(1:10), '(i10)') 1
    buffer = 'SPACING '//str1//' '//str2//' '//str3//lf; write (ivtk) trim(buffer)

    num = int(nxGlobal, kind=8)*int(nyGlobal, kind=8)*int(nzGlobal, kind=8)

    write (str4(1:14), '(i14)') num
    buffer = 'POINT_DATA '//str4//lf; write (ivtk) trim(buffer)

    if (output_precision_choice == 1.or.input_precision_choice == 1) then   ! dp in sp out / sp in sp out
        fmt = 'float'
        buffer = 'SCALARS density '//fmt//lf
        write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf
        write (ivtk) trim(buffer)
        write (ivtk) (((real(rho(i, j, k)), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
        buffer = 'SCALARS velocity_X '//fmt//lf
        write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf
        write (ivtk) trim(buffer)
        write (ivtk) (((real(u(i, j, k)), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
        buffer = 'SCALARS velocity_Y '//fmt//lf
        write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf
        write (ivtk) trim(buffer)
        write (ivtk) (((real(v(i, j, k)), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
        buffer = 'SCALARS velocity_Z '//fmt//lf
        write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf
        write (ivtk) trim(buffer)
        write (ivtk) (((real(w(i, j, k)), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
    else
        fmt = 'double'
        buffer = 'SCALARS density '//fmt//lf
        write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf
        write (ivtk) trim(buffer)
        write (ivtk) (((rho(i, j, k), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
        buffer = 'SCALARS velocity_X '//fmt//lf
        write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf
        write (ivtk) trim(buffer)
        write (ivtk) (((u(i, j, k), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
        buffer = 'SCALARS velocity_Y '//fmt//lf
        write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf
        write (ivtk) trim(buffer)
        write (ivtk) (((v(i, j, k), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
        buffer = 'SCALARS velocity_Z '//fmt//lf
        write (ivtk) trim(buffer)
        buffer = 'LOOKUP_TABLE default'//lf
        write (ivtk) trim(buffer)
        write (ivtk) (((w(i, j, k), i=1, nxGlobal), j=1, nyGlobal), k=1, nzGlobal)
    end if

    close (ivtk)

    return
end

subroutine VTK_detail_bin_half_sp(nt)
    use Misc
    IMPLICIT NONE
    character*30 :: flnm
    integer :: i, j, k, nt
    integer(kind=8) :: num
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10, str4*14
    integer   :: ivtk = 9, int

    write (flnm, '(i10.10,".vtk")') nt

   open (unit=ivtk, file=trim(data_location)//'/processed_field_data/detail_half_'//flnm, FORM='unformatted', access='stream', status='replace', convert='BIG_ENDIAN')

    lf = char(10) ! line feed character

    buffer = '# vtk DataFile Version 3.0'//lf; write (ivtk) trim(buffer)
    buffer = 'vtk output'//lf; write (ivtk) trim(buffer)
    buffer = 'BINARY'//lf; write (ivtk) trim(buffer)
    buffer = 'DATASET STRUCTURED_POINTS '//lf; write (ivtk) trim(buffer)

    write (str1(1:10), '(i10)') nxglobal/2
    write (str2(1:10), '(i10)') nyglobal/2
    write (str3(1:10), '(i10)') nzglobal/2
    buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf; write (ivtk) trim(buffer)

    write (str1(1:10), '(i10)') 1
    write (str2(1:10), '(i10)') 1
    write (str3(1:10), '(i10)') 1
    buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf; write (ivtk) trim(buffer)

    write (str1(1:10), '(i10)') 1
    write (str2(1:10), '(i10)') 1
    write (str3(1:10), '(i10)') 1
    buffer = 'SPACING '//str1//' '//str2//' '//str3//lf; write (ivtk) trim(buffer)

    num = int(nxGlobal, kind=8)*int(nyGlobal, kind=8)*int(nzGlobal, kind=8)/8

    write (str4(1:14), '(i14)') num
    buffer = 'POINT_DATA '//str4//lf; write (ivtk) trim(buffer)

    !scalar - density
    buffer = 'SCALARS density float'//lf
    write (ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write (ivtk) trim(buffer)
    write (ivtk) (((real(rho(i, j, k)), i=1, nxGlobal, 2), j=1, nyGlobal, 2), k=1, nzGlobal, 2)

    buffer = 'SCALARS velocity_X float'//lf
    write (ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write (ivtk) trim(buffer)
    write (ivtk) (((real(u(i, j, k)), i=1, nxGlobal, 2), j=1, nyGlobal, 2), k=1, nzGlobal, 2)

    buffer = 'SCALARS velocity_Y float'//lf
    write (ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write (ivtk) trim(buffer)
    write (ivtk) (((real(v(i, j, k)), i=1, nxGlobal, 2), j=1, nyGlobal, 2), k=1, nzGlobal, 2)

    buffer = 'SCALARS velocity_Z float'//lf
    write (ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write (ivtk) trim(buffer)
    write (ivtk) (((real(w(i, j, k)), i=1, nxGlobal, 2), j=1, nyGlobal, 2), k=1, nzGlobal, 2)

    close (ivtk)

    return
end

!************************************* VTK legacy output **************************************
