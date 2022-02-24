!=====================================================================================================================================
!------- set up MPI basics (follow Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]'s code) --------
!=====================================================================================================================================
SUBROUTINE set_MPI
    use Misc_module
    use mpi_variable
    implicit none
    include 'mpif.h'
    INTEGER :: complete, direction, partial, shift
    INTEGER :: MPI_ERR, x, xeast, xwest, y, ynorth, ysouth, z, ztop, zbottom
    INTEGER :: xmin, xmax, xl, xu, xlg, xug
    INTEGER :: ymin, ymax, yl, yu, ylg, yug
    INTEGER :: zmin, zmax, zl, zu, zlg, zug
    INTEGER, DIMENSION(1:mpi_dim) :: dims
    LOGICAL, DIMENSION(1:mpi_dim) :: periodic
    LOGICAL :: reorder

    dims(1) = npx
    dims(2) = npy
    dims(3) = npz
    periodic(1) = .true.
    periodic(2) = .true.
    periodic(3) = .true.
    reorder = .false.
    xmin = 1
    xmax = nxGlobal
    ymin = 1
    ymax = nyGlobal
    zmin = 1
    zmax = nzGlobal

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, mpi_dim, dims, periodic, reorder, MPI_COMM_VGRID, MPI_ERR)
    CALL MPI_COMM_RANK(MPI_COMM_VGRID, id, MPI_ERR)
    CALL MPI_CART_COORDS(MPI_COMM_VGRID, id, mpi_dim, mpi_coords, MPI_ERR)
    idx = mpi_coords(1)
    idy = mpi_coords(2)
    idz = mpi_coords(3)

    !-------- Compute limits [ (xl,xu) x (yl,yu) x (zl,zu)] for this processor -----
    !  Partitioning in the x direction
    complete = (xmax - xmin)/dims(1)
    partial = (xmax - xmin) - complete*dims(1)
    IF (idx + 1 <= partial) THEN
        xl = xmin + (complete + 1)*idx
        xu = xmin + (complete + 1)*(idx + 1) - 1
    ELSE
        xl = xmin + complete*idx + partial
        xu = xmin + complete*(idx + 1) + partial - 1
    END IF
    IF (MOD(idx + 1, dims(1)) == 0) xu = xu + 1
    !  Partitioning in the y direction
    complete = (ymax - ymin)/dims(2)
    partial = (ymax - ymin) - complete*dims(2)
    IF (idy + 1 <= partial) THEN
        yl = ymin + (complete + 1)*idy
        yu = ymin + (complete + 1)*(idy + 1) - 1
    ELSE
        yl = ymin + complete*idy + partial
        yu = ymin + complete*(idy + 1) + partial - 1
    END IF
    IF (MOD(idy + 1, dims(2)) == 0) yu = yu + 1
    !  Partitioning in the y direction
    complete = (zmax - zmin)/dims(3)
    partial = (zmax - zmin) - complete*dims(3)
    IF (idz + 1 <= partial) THEN
        zl = zmin + (complete + 1)*idz
        zu = zmin + (complete + 1)*(idz + 1) - 1
    ELSE
        zl = zmin + complete*idz + partial
        zu = zmin + complete*(idz + 1) + partial - 1
    END IF
    IF (MOD(idz + 1, dims(3)) == 0) zu = zu + 1

    ! Set domain limits
    NX = xu - xl + 1
    NY = yu - yl + 1
    NZ = zu - zl + 1

    shift = 1
    direction = 0
    CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, idxM, idxP, MPI_ERR)
    direction = 1
    CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, idyM, idyP, MPI_ERR)
    direction = 2
    CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, idzM, idzP, MPI_ERR)

    CALL MPI_CART_COORDS(MPI_COMM_VGRID, idxP, mpi_dim, mpi_coords, MPI_ERR)
    xeast = mpi_coords(1)
    CALL MPI_CART_COORDS(MPI_COMM_VGRID, idxM, mpi_dim, mpi_coords, MPI_ERR)
    xwest = mpi_coords(1)
    CALL MPI_CART_COORDS(MPI_COMM_VGRID, idyP, mpi_dim, mpi_coords, MPI_ERR)
    ynorth = mpi_coords(2)
    CALL MPI_CART_COORDS(MPI_COMM_VGRID, idyM, mpi_dim, mpi_coords, MPI_ERR)
    ysouth = mpi_coords(2)
    CALL MPI_CART_COORDS(MPI_COMM_VGRID, idzP, mpi_dim, mpi_coords, MPI_ERR)
    ztop = mpi_coords(3)
    CALL MPI_CART_COORDS(MPI_COMM_VGRID, idzM, mpi_dim, mpi_coords, MPI_ERR)
    zbottom = mpi_coords(3)

    !  Neighbors along x edges
    mpi_coords(1) = idx      ! current node x
    mpi_coords(2) = ynorth ! y plus node y
    mpi_coords(3) = ztop   ! z plus node z
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idyPzP, MPI_ERR)
    mpi_coords(1) = idx
    mpi_coords(2) = ysouth
    mpi_coords(3) = ztop
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idyMzP, MPI_ERR)
    mpi_coords(1) = idx
    mpi_coords(2) = ysouth
    mpi_coords(3) = zbottom
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idyMzM, MPI_ERR)
    mpi_coords(1) = idx
    mpi_coords(2) = ynorth
    mpi_coords(3) = zbottom
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idyPzM, MPI_ERR)

    !  Neighbors along y edges
    mpi_coords(1) = xeast
    mpi_coords(2) = idy
    mpi_coords(3) = ztop
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxPzP, MPI_ERR)
    mpi_coords(1) = xwest
    mpi_coords(2) = idy
    mpi_coords(3) = ztop
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxMzP, MPI_ERR)
    mpi_coords(1) = xwest
    mpi_coords(2) = idy
    mpi_coords(3) = zbottom
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxMzM, MPI_ERR)
    mpi_coords(1) = xeast
    mpi_coords(2) = idy
    mpi_coords(3) = zbottom
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxPzM, MPI_ERR)

    !  Neighbors along z edges
    mpi_coords(1) = xeast
    mpi_coords(2) = ynorth
    mpi_coords(3) = idz
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxPyP, MPI_ERR)
    mpi_coords(1) = xwest
    mpi_coords(2) = ysouth
    mpi_coords(3) = idz
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxMyM, MPI_ERR)
    mpi_coords(1) = xwest
    mpi_coords(2) = ynorth
    mpi_coords(3) = idz
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxMyP, MPI_ERR)
    mpi_coords(1) = xeast
    mpi_coords(2) = ysouth
    mpi_coords(3) = idz
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxPyM, MPI_ERR)

    !  Neighbors corners (not used in this code as no domain decomposition along x direction)
    mpi_coords(1) = xeast
    mpi_coords(2) = ynorth
    mpi_coords(3) = ztop
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxPyPzP, MPI_ERR)
    mpi_coords(1) = xeast
    mpi_coords(2) = ysouth
    mpi_coords(3) = ztop
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxPyMzP, MPI_ERR)
    mpi_coords(1) = xeast
    mpi_coords(2) = ynorth
    mpi_coords(3) = zbottom
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxPyPzM, MPI_ERR)
    mpi_coords(1) = xeast
    mpi_coords(2) = ysouth
    mpi_coords(3) = zbottom
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxPyMzM, MPI_ERR)
    mpi_coords(1) = xwest
    mpi_coords(2) = ynorth
    mpi_coords(3) = ztop
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxMyPzP, MPI_ERR)
    mpi_coords(1) = xwest
    mpi_coords(2) = ysouth
    mpi_coords(3) = ztop
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxMyMzP, MPI_ERR)
    mpi_coords(1) = xwest
    mpi_coords(2) = ynorth
    mpi_coords(3) = zbottom
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxMyPzM, MPI_ERR)
    mpi_coords(1) = xwest
    mpi_coords(2) = ysouth
    mpi_coords(3) = zbottom
    CALL MPI_CART_RANK(MPI_COMM_VGRID, mpi_coords, idxMyMzM, MPI_ERR)

    RETURN
END SUBROUTINE set_MPI

!=====================================================================================================================================
!---------------------- MPI data gather/scatter ----------------------
!=====================================================================================================================================
subroutine AllGather(gg1, i1, i2, j1, j2, k1, k2, rg, ff, l1, l2, m1, m2, n1, n2, rf)
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i1, i2, j1, j2, k1, k2, rg, l1, l2, m1, m2, n1, n2, rf
    real(kind=8)  :: gg1(i1 - rg:i2 + rg, j1 - rg:j2 + rg, k1 - rg:k2 + rg)
    real(kind=8)  :: ff(l1 - rf:l2 + rf, m1 - rf:m2 + rf, n1 - rf:n2 + rf)
    real(kind=8), allocatable, dimension(:, :, :) :: tt
    integer :: i, j, k, L, M, N, out1, out2, out3, rank, isize, status(MPI_STATUS_SIZE)

    isize = (i2 - i1 + 1)*(j2 - j1 + 1)*(k2 - k1 + 1)
    allocate (tt(i1:i2, j1:j2, k1:k2))

    if (id .ne. 0) then
        !$omp PARALLEL DO private(i,j)
        do k = k1, k2
            do j = j1, j2
                do i = i1, i2
                    tt(i, j, k) = gg1(i, j, k)
                end do
            end do
        end do
        call mpi_send(tt, isize, MPI_DOUBLE_PRECISION, 0, 600 + id, MPI_COMM_VGRID, ierr)
    else
        do rank = 1, np - 1
            call mpi_recv(tt, isize, MPI_DOUBLE_PRECISION, rank, 600 + rank, MPI_COMM_VGRID, status, ierr)

            CALL MPI_CART_COORDS(MPI_COMM_VGRID, rank, mpi_dim, mpi_coords, ierr)
            out1 = mpi_coords(1)
            out2 = mpi_coords(2)
            out3 = mpi_coords(3)

            !$omp PARALLEL DO private(i,j,l,m,n)
            do k = k1, k2
                do j = j1, j2
                    do i = i1, i2
                        L = out1*nx + i
                        M = out2*ny + j
                        N = out3*nz + k
                        ff(L, M, N) = tt(i, j, k)   !for id=1 to np-1
                    end do
                end do
            end do
        end do

        CALL MPI_CART_COORDS(MPI_COMM_VGRID, id, mpi_dim, mpi_coords, ierr)
        out1 = mpi_coords(1)
        out2 = mpi_coords(2)
        out3 = mpi_coords(3)

        !$omp PARALLEL DO private(i,j,l,m,n)
        do k = k1, k2
            do j = j1, j2
                do i = i1, i2
                    L = out1*nx + i
                    M = out2*ny + j
                    N = out3*nz + k
                    ff(L, M, N) = gg1(i, j, k)   !for id=0
                end do
            end do
        end do
    end if

    deallocate (tt)
    return
end subroutine AllGather
!~~~~~~~~~~~~~~~~~~~~~~~~~ single precision version ~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine AllGather_SP(gg1, i1, i2, j1, j2, k1, k2, rg, ff, l1, l2, m1, m2, n1, n2, rf)
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i1, i2, j1, j2, k1, k2, rg, l1, l2, m1, m2, n1, n2, rf
    real(kind=8) :: gg1(i1 - rg:i2 + rg, j1 - rg:j2 + rg, k1 - rg:k2 + rg)   !keep original precision
    real(kind=4)  :: ff(l1 - rf:l2 + rf, m1 - rf:m2 + rf, n1 - rf:n2 + rf)
    real(kind=4), allocatable, dimension(:, :, :) :: tt
    integer :: i, j, k, L, M, N, out1, out2, out3, rank, isize, status(MPI_STATUS_SIZE)

    isize = (i2 - i1 + 1)*(j2 - j1 + 1)*(k2 - k1 + 1)
    allocate (tt(i1:i2, j1:j2, k1:k2))
    if (id .ne. 0) then
        !$omp PARALLEL DO private(i,j)
        do k = k1, k2
            do j = j1, j2
                do i = i1, i2
                    tt(i, j, k) = gg1(i, j, k)
                end do
            end do
        end do
        call mpi_send(tt, isize, MPI_FLOAT, 0, 600 + id, MPI_COMM_VGRID, ierr)
    else
        do rank = 1, np - 1
            call mpi_recv(tt, isize, MPI_FLOAT, rank, 600 + rank, MPI_COMM_VGRID, status, ierr)
            CALL MPI_CART_COORDS(MPI_COMM_VGRID, rank, mpi_dim, mpi_coords, ierr)
            out1 = mpi_coords(1)
            out2 = mpi_coords(2)
            out3 = mpi_coords(3)

            !$omp PARALLEL DO private(i,j,l,m,n)
            do k = k1, k2
                do j = j1, j2
                    do i = i1, i2
                        L = out1*nx + i
                        M = out2*ny + j
                        N = out3*nz + k
                        ff(L, M, N) = tt(i, j, k)   !for id=1 to np-1
                    end do
                end do
            end do
        end do

        CALL MPI_CART_COORDS(MPI_COMM_VGRID, id, mpi_dim, mpi_coords, ierr)
        out1 = mpi_coords(1)
        out2 = mpi_coords(2)
        out3 = mpi_coords(3)

        !$omp PARALLEL DO private(i,j,l,m,n)
        do k = k1, k2
            do j = j1, j2
                do i = i1, i2
                    L = out1*nx + i
                    M = out2*ny + j
                    N = out3*nz + k
                    ff(L, M, N) = gg1(i, j, k)   !for id=0
                end do
            end do
        end do
    end if

    deallocate (tt)
    return
end subroutine AllGather_SP

!=====================================================================================================================================
!------------------- Geometry info (int 1 data type) MPI transfer (not used in GPU computing region) ----------------------
!=====================================================================================================================================
!********************** walls transfer along x axis *****************************
subroutine xtransport_walls(lx, ly, lz)
    use Misc_module
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i, j, k, m, m1, m2, lx, ly, lz
    integer :: isize, status(MPI_STATUS_SIZE)
    integer, allocatable, dimension(:, :, :)::send1, send2, recv1, recv2

    allocate ( &
        send1(lx, 1 - ly:ny + ly, 1 - lz:nz + lz), &
        send2(lx, 1 - ly:ny + ly, 1 - lz:nz + lz), &
        recv1(lx, 1 - ly:ny + ly, 1 - lz:nz + lz), &
        recv2(lx, 1 - ly:ny + ly, 1 - lz:nz + lz))
    isize = (ny + ly + ly)*(nz + lz + lz)*lx

    !$omp PARALLEL DO private(i,j)
    do k = 1 - lz, nz + lz
        do j = 1 - ly, ny + ly
            do i = 1, lx
                send1(i, j, k) = walls(i, j, k)
                send2(i, j, k) = walls(nx + i - lx, j, k)
            end do
        end do
    end do

    call mpi_sendrecv(send1, isize, MPI_integer, idxM, 9000, recv2, &
                      isize, MPI_integer, idxP, 9000, MPI_COMM_VGRID, status, ierr)
    call mpi_sendrecv(send2, isize, MPI_integer, idxP, 8000, recv1, &
                      isize, MPI_integer, idxM, 8000, MPI_COMM_VGRID, status, ierr)
    if (iper .eq. 1 .or. idx .ne. 0) then   !  iper indicates perodic condition
        !$omp PARALLEL DO private(i,j)
        do k = 1 - lz, nz + lz
            do j = 1 - ly, ny + ly
                do i = 1, lx
                    walls(i - lx, j, k) = recv1(i, j, k)
                end do
            end do
        end do
    end if
    if (iper .eq. 1 .or. idx .ne. npx - 1) then
        !$omp PARALLEL DO private(i,j)
        do k = 1 - lz, nz + lz
            do j = 1 - ly, ny + ly
                do i = 1, lx
                    walls(i + nx, j, k) = recv2(i, j, k)
                end do
            end do
        end do
    end if

    deallocate (send1, send2, recv1, recv2)
    return
end subroutine xtransport_walls

!********************** walls transfer along y axis *****************************
subroutine ytransport_walls(lx, ly, lz)
    use Misc_module
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i, j, k, m, m1, m2, lx, ly, lz
    integer :: isize, status(MPI_STATUS_SIZE)
    integer, allocatable, dimension(:, :, :)::send1, send2, recv1, recv2

    allocate ( &
        send1(1 - lx:nx + lx, ly, 1 - lz:nz + lz), &
        send2(1 - lx:nx + lx, ly, 1 - lz:nz + lz), &
        recv1(1 - lx:nx + lx, ly, 1 - lz:nz + lz), &
        recv2(1 - lx:nx + lx, ly, 1 - lz:nz + lz))
    isize = (nx + lx + lx)*(nz + lz + lz)*ly

    !$omp PARALLEL DO private(i,j)
    do k = 1 - lz, nz + lz
        do j = 1, ly
            do i = 1 - lx, nx + lx
                send1(i, j, k) = walls(i, j, k)
                send2(i, j, k) = walls(i, ny + j - ly, k)
            end do
        end do
    end do

    call mpi_sendrecv(send1, isize, MPI_integer, idyM, 9000, recv2, &
                      isize, MPI_integer, idyP, 9000, MPI_COMM_VGRID, status, ierr)
    call mpi_sendrecv(send2, isize, MPI_integer, idyP, 8000, recv1, &
                      isize, MPI_integer, idyM, 8000, MPI_COMM_VGRID, status, ierr)
    if (jper .eq. 1 .or. idy .ne. 0) then   !  iper indicates perodic condition
        !$omp PARALLEL DO private(i,j)
        do k = 1 - lz, nz + lz
            do j = 1, ly
                do i = 1 - lx, nx + lx
                    walls(i, j - ly, k) = recv1(i, j, k)
                end do
            end do
        end do
    end if
    if (jper .eq. 1 .or. idy .ne. npy - 1) then
        !$omp PARALLEL DO private(i,j)
        do k = 1 - lz, nz + lz
            do j = 1, ly
                do i = 1 - lx, nx + lx
                    walls(i, j + ny, k) = recv2(i, j, k)
                end do
            end do
        end do
    end if

    deallocate (send1, send2, recv1, recv2)
    return
end subroutine ytransport_walls

!********************** walls transfer along z axis *****************************
subroutine ztransport_walls(lx, ly, lz)
    use Misc_module
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i, j, k, m, m1, m2, lx, ly, lz
    integer :: isize, status(MPI_STATUS_SIZE)
    integer, allocatable, dimension(:, :, :)::send1, send2, recv1, recv2

    allocate ( &
        send1(1 - lx:nx + lx, 1 - ly:ny + ly, lz), &
        send2(1 - lx:nx + lx, 1 - ly:ny + ly, lz), &
        recv1(1 - lx:nx + lx, 1 - ly:ny + ly, lz), &
        recv2(1 - lx:nx + lx, 1 - ly:ny + ly, lz))
    isize = (nx + lx + lx)*(ny + ly + ly)*lz

    !$omp PARALLEL DO private(i,j)
    do k = 1, lz
        do j = 1 - ly, ny + ly
            do i = 1 - lx, nx + lx
                send1(i, j, k) = walls(i, j, k)
                send2(i, j, k) = walls(i, j, nz + k - lz)
            end do
        end do
    end do

    call mpi_sendrecv(send1, isize, MPI_integer, idzM, 9000, recv2, &
                      isize, MPI_integer, idzP, 9000, MPI_COMM_VGRID, status, ierr)
    call mpi_sendrecv(send2, isize, MPI_integer, idzP, 8000, recv1, &
                      isize, MPI_integer, idzM, 8000, MPI_COMM_VGRID, status, ierr)
    if (kper .eq. 1 .or. idz .ne. 0) then   !  iper indicates perodic condition
        !$omp PARALLEL DO private(i,j)
        do k = 1, lz
            do j = 1 - ly, ny + ly
                do i = 1 - lx, nx + lx
                    walls(i, j, k - lz) = recv1(i, j, k)
                end do
            end do
        end do
    end if
    if (kper .eq. 1 .or. idz .ne. npz - 1) then
        !$omp PARALLEL DO private(i,j)
        do k = 1, lz
            do j = 1 - ly, ny + ly
                do i = 1 - lx, nx + lx
                    walls(i, j, k + nz) = recv2(i, j, k)
                end do
            end do
        end do
    end if

    deallocate (send1, send2, recv1, recv2)
    return
end subroutine ztransport_walls
!********************** walls transfer along z axis *****************************
