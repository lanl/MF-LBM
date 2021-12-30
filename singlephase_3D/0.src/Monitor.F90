!=======================================================================================================================================================
!---------------------- monitor_singlephase  unsteady flow  ----------------------
!=======================================================================================================================================================
subroutine monitor
    use Misc_module
    use Fluid_singlephase
    use mpi_variable
    USE,INTRINSIC :: IEEE_ARITHMETIC
    IMPLICIT NONE
    include 'mpif.h'
    real(kind=8),allocatable,dimension(:) :: fl_0,pre_0
    character (len=20) :: flnm,flnm1   !file name
    integer :: i,j,k ,L,M,N,rank,status(MPI_STATUS_SIZE),o1,o2,o3, icount, itemp
    integer (kind=1) :: wall_indicator
    real(kind=8) :: umax,temp,mass,temp1,temp2,temp3
    real(kind=8) :: prek,fl_avg,tmp,fl_avg_whole, rel_err

    !********************************* preperation ***********************************

    call compute_macro_vars

    !$omp parallel

    prek=0d0
    temp=0d0
    umax=0d0
    !$omp do private(i,j,prek,temp)
    !$acc kernels present(fl,pre,rho)
    !$acc loop
    do k=1,nz
        temp=0d0
        prek =0d0
        !$acc loop reduction(+:prek,temp) collapse(2)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,k)
                temp=temp+w(i,j,k) 
                prek = prek + rho(i,j,k)   !calculate average bulk pressure
            enddo
        enddo
        pre(k)=prek                       !pressure sum, further divided by wallz_prof to obtain averaged pressure profile along the flow direction z
        fl(k)=temp
    enddo
    !$acc end kernels

    !$omp do private(i,j,temp) reduction(max:umax)
    !$acc kernels present(u,v,w)
    !$acc loop reduction(max:umax)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                temp = u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)
                if(umax<temp)then
                    umax=temp
                endif
            enddo
        enddo
    enddo
    !$acc end kernels
    !$omp end parallel

    !******************************* packing and send/recv data **********************************
    if(id.ne.0)then
        !packing the data for transfer
        !$omp parallel do
        !$acc parallel present(fl,pre,tk)
        !$acc loop independent
        do k=1,nz
            tk(k)     =fl(k)
            tk(nz+k)  =pre(k)
        enddo
        !$acc end parallel
        !$acc update host (tk)
        tk(2*nz+1)=umax
        call mpi_send(tk,tk_isize,MPI_DOUBLE_PRECISION,0,600+id,MPI_COMM_VGRID,ierr)
    else
        allocate(fl_0(nzGlobal),pre_0(nzGlobal))
        !$acc update host (fl,pre)
        !$omp parallel
        !$omp do
        do k=1,nzGlobal
            fl_0(k)=0.0d0
            pre_0(k)=0.0d0
        enddo
        !$omp do
        do k=1,nz
            fl_0(idz*nz+k)= fl_0(idz*nz+k) + fl(k)
            pre_0(idz*nz+k)= pre_0(idz*nz+k) + pre(k)
        enddo
        !$omp end parallel
        do rank = 1, np-1
            call mpi_recv(tk,tk_isize,MPI_DOUBLE_PRECISION,rank,600+rank,MPI_COMM_VGRID,status,ierr)

            CALL MPI_CART_COORDS(MPI_COMM_VGRID, rank, mpi_dim, mpi_coords, ierr)
            o3=mpi_coords(3)

            !$omp parallel do
            do k=1,nz
                fl_0(o3*nz+k)= fl_0(o3*nz+k) + tk(k)
                pre_0(o3*nz+k)= pre_0(o3*nz+k) + tk(nz+k)
            enddo
            if(umax<tk(2*nz+1))then
                umax=tk(2*nz+1)
            endif
        enddo

        umax_global = dsqrt(umax)

        !********************************************* flowrate calculation **********************************************
        fl_avg_whole = 0d0
        
        !$omp parallel do reduction (+:fl_avg_whole)
        do k=1,nzGlobal   !average velocity of whole domain
            fl_avg_whole = fl_avg_whole + fl_0(k)
        enddo      
        fl_avg_whole = fl_avg_whole/dble(nzGlobal)
        flowrate = fl_avg_whole

        !********************************************* save data ************************************************************
        open(unit=13, file='out1.output/flowrate_time.dat' ,status='unknown',position='append')
        write(13,"(I10,2(1x,E14.6))")ntime,fl_avg_whole,umax
        close(13)      

        if(kper==0.and.domain_wall_status_z_min==0.and.domain_wall_status_z_max==0)then    !non-periodic BC along flow direction (z)    
            temp1 = pre_0(1+n_exclude_inlet)/pore_profile_z(1+n_exclude_inlet)   
            temp2 = pre_0(nzglobal-n_exclude_outlet)/pore_profile_z(nzglobal-n_exclude_outlet)
            temp3 = temp1 - temp2
            open(unit=14, file='out1.output/pre.dat' ,status='unknown',position='append')
            write(14,"(I10,3(1x,e14.7))")ntime,temp1,temp2,temp3
            close(14)
        endif
    
        if(MOD(ntime,ntime_monitor_profile)==0)then
            write(flnm,'(I10.8,".dat")')ntime
            flnm1 = trim(flnm)//trim(flnm1)   !?
            open(unit=8,file='out1.output/profile/monitor'//flnm,status='replace')
            do k=1,nzGlobal
                write(8,"(3(1x,e14.7))")dble(k),fl_0(k),pre_0(k)/(dble(pore_profile_z(k))+eps)
            enddo
            close(8)
        endif

        deallocate(fl_0,pre_0)

        if(steady_state_option==1)then !steady state simulation based on flowrate
          monitor_previous_value = monitor_current_value  ! store previous step
          monitor_current_value  = flowrate   
          rel_err = dabs(monitor_current_value-monitor_previous_value)/(dabs(monitor_current_value)+eps)
          open(unit=15, file='out1.output/steady_monitor_flowrate_error.dat' ,status='unknown',position='append')
          write(15,"(I10,(1x,e14.7))")ntime,rel_err
          close(15)
          if(rel_err<convergence_criteria.and.ntime>ntime_monitor)then
            print*, 'Simulation converged based on flowrate convergence! Relative error =', rel_err
            simulation_end_indicator = 1   !successfully end simulation
          endif
        endif

        ! check simulation status
        if(ieee_is_nan(flowrate))then
            print*, 'Simulation failed due to NAN of "flowrate"!'
            simulation_end_indicator = 3
        elseif(umax_global > 0.5)then
            print*, 'Simulation failed due to maximum velocity larger than 0.5!'
            simulation_end_indicator = 3
        endif 
        
    endif

    call MPI_Bcast(simulation_end_indicator,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)

    return
end subroutine monitor




