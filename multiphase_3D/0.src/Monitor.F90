#include "./preprocessor.h"
!=======================================================================================================================================================
!---------------------- monitor_multiphase  unsteady flow  ----------------------
!=======================================================================================================================================================
subroutine monitor
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    USE,INTRINSIC :: IEEE_ARITHMETIC
    IMPLICIT NONE
    include 'mpif.h'
    real(kind=8),allocatable,dimension(:) :: fl1_0,fl2_0,pre_0,mass1_0,mass2_0,vol1_0,vol2_0
    character (len=20) :: flnm,flnm1   !file name
    integer :: i,j,k ,L,M,N,rank,status(MPI_STATUS_SIZE),o1,o2,o3, icount, itemp
    integer (kind=1) :: wall_indicator
    real(kind=8) :: umax,temp,mass,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,usq1,usq2
    real(kind=8) :: prek,fl1_avg,fl2_avg,fl_avg,tmp,fx,fy,fz,fl1_avg_whole,fl2_avg_whole,fl_avg_whole
    real(kind=8) :: slope1, slope2,slope3, av1, av2 ,av3 ! least square
    real(kind=8), dimension(3) :: w_local

    !********************************* preperation ***********************************
    prek=0d0
    temp=0d0
    umax=0d0

    call compute_macro_vars


    usq1=0d0   !kinetic energy of each phase
    usq2=0d0
    !$omp parallel do private(i,j,temp,wall_indicator) reduction(max:umax)reduction(+:usq1,usq2)
    !$acc kernels present(u,v,w,walls,phi)
    !$acc loop reduction(max:umax)reduction(+:usq1,usq2)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,k)
                temp = (u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))* (1-wall_indicator)
                if(umax<temp)then
                    umax=temp
                endif
                if(phi(i,j,k)>0.999d0)then
                    usq1 = usq1 + temp
                elseif(phi(i,j,k)<-0.999d0)then
                    usq2 = usq2 + temp
                endif
            enddo
        enddo
    enddo
    !$acc end kernels

    !$omp parallel do private(i,j,prek,temp1,temp2,temp3,temp4,temp5,temp6,wall_indicator)
    !$acc kernels present(fl1,fl2,mass1,mass2,vol1,vol2,pre,rho,phi,walls)
    !$acc loop
    do k=1,nz
        temp1=0d0
        temp2=0d0
        prek =0d0
        temp3=0d0
        temp4=0d0
        temp5=0d0
        temp6=0d0
        !$acc loop reduction(+:prek,temp1,temp2,temp3,temp4,temp5,temp6) collapse(2)
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,k)
                temp3= temp3 + 0.5d0*(1d0+phi(i,j,k))* (1-wall_indicator)   !volume fraction
                temp4= temp4 + 0.5d0*(1d0-phi(i,j,k))* (1-wall_indicator)
                temp5 = temp5 + rho(i,j,k)*0.5d0*(1d0+phi(i,j,k))* (1-wall_indicator)  !mass fraction
                temp6 = temp6 + rho(i,j,k)*0.5d0*(1d0-phi(i,j,k))* (1-wall_indicator)   
                temp1=temp1+w(i,j,k)*0.5d0*(1d0+phi(i,j,k))* (1-wall_indicator)     !volume fractional flow fluid 1
                temp2=temp2+w(i,j,k)*0.5d0*(1d0-phi(i,j,k))* (1-wall_indicator)     !volume fractional flow fluid 2
                prek = prek + rho(i,j,k)* (1-wall_indicator)    !calculate average bulk pressure
            enddo
        enddo
        pre(k)=prek                       !pressure sum, further divided by wallz_prof to obtain averaged pressure profile along the flow direction z
        fl1(k)=temp1
        fl2(k)=temp2
        vol1(k)= temp3
        vol2(k)= temp4
        mass1(k)= temp5
        mass2(k)= temp6
    enddo
    !$acc end kernels

    !******************************* packing and send/recv data **********************************
    if(id.ne.0)then
        !packing the data for transfer
        !$omp parallel do
        !$acc parallel present(fl1,fl2,mass1,mass2,vol1,vol2,pre,tk)
        !$acc loop independent
        do k=1,nz
            tk(k)     =fl1(k)
            tk(nz+k)  =fl2(k)
            tk(2*nz+k)=vol1(k)
            tk(3*nz+k)=vol2(k)
            tk(4*nz+k)=mass1(k)
            tk(5*nz+k)=mass2(k)
            tk(6*nz+k)=Pre(k)
        enddo
        !$acc end parallel
        !$acc update host (tk)
        tk(7*nz+1)=umax
        tk(7*nz+2)=usq1
        tk(7*nz+3)=usq2
        call mpi_send(tk,tk_isize,MPI_DOUBLE_PRECISION,0,600+id,MPI_COMM_VGRID,ierr)
    else
        allocate(fl1_0(nzGlobal),fl2_0(nzGlobal),vol1_0(nzGlobal),vol2_0(nzGlobal),mass1_0(nzGlobal),mass2_0(nzGlobal),pre_0(nzGlobal))
        !$acc update host (fl1,fl2,mass1,mass2,vol1,vol2,pre)
        !$omp parallel do
        do k=1,nzGlobal
            fl1_0(k)=0.0d0
            fl2_0(k)=0.0d0
            vol1_0(k)=0.0d0
            vol2_0(k)=0.0d0
            mass1_0(k)=0.0d0
            mass2_0(k)=0.0d0
            pre_0(k)=0.0d0
        enddo
        !$omp do
        do k=1,nz
            fl1_0(idz*nz+k)= fl1_0(idz*nz+k) + fl1(k)
            fl2_0(idz*nz+k)= fl2_0(idz*nz+k) + fl2(k)
            vol1_0(idz*nz+k)= vol1_0(idz*nz+k) + vol1(k)
            vol2_0(idz*nz+k)= vol2_0(idz*nz+k) + vol2(k)
            mass1_0(idz*nz+k)= mass1_0(idz*nz+k) + mass1(k)
            mass2_0(idz*nz+k)= mass2_0(idz*nz+k) + mass2(k)
            pre_0(idz*nz+k)= pre_0(idz*nz+k) + pre(k)
        enddo

        do rank = 1, np-1
            call mpi_recv(tk,tk_isize,MPI_DOUBLE_PRECISION,rank,600+rank,MPI_COMM_VGRID,status,ierr)
            CALL MPI_CART_COORDS(MPI_COMM_VGRID, rank, mpi_dim, mpi_coords, ierr)
            o3=mpi_coords(3)
            !$omp parallel do
            do k=1,nz
                fl1_0(o3*nz+k)= fl1_0(o3*nz+k) + tk(k)
                fl2_0(o3*nz+k)= fl2_0(o3*nz+k) + tk(nz+k)
                vol1_0(o3*nz+k)= vol1_0(o3*nz+k) + tk(2*nz+k)
                vol2_0(o3*nz+k)= vol2_0(o3*nz+k) + tk(3*nz+k)
                mass1_0(o3*nz+k)= mass1_0(o3*nz+k) + tk(4*nz+k)
                mass2_0(o3*nz+k)= mass2_0(o3*nz+k) + tk(5*nz+k)
                pre_0(o3*nz+k)= pre_0(o3*nz+k) + tk(6*nz+k)
            enddo
            if(umax<tk(7*nz+1))then
                umax=tk(7*nz+1)
            endif
            usq1 = usq1 + tk(7*nz+2)
            usq2 = usq2 + tk(7*nz+3)
        enddo

        kinetic_energy(1) = 0.5d0*usq1
        kinetic_energy(2) = 0.5d0*usq2

        umax_global = dsqrt(umax)

        !********************************************* saturation calculation ********************************************
        mass1_sum=0d0
        mass2_sum=0d0
        vol1_sum=0d0
        vol2_sum=0d0
        !$omp parallel do reduction (+:mass1_sum,mass2_sum,vol1_sum,vol2_sum)
        do k=n_exclude_inlet+1,nzGlobal-n_exclude_outlet
            mass1_sum = mass1_sum + mass1_0(k)
            mass2_sum = mass2_sum + mass2_0(k)           
            vol1_sum = vol1_sum + vol1_0(k)
            vol2_sum = vol2_sum + vol2_0(k)
            !print*,mass2_0(k),k     
        enddo
        
        saturation =vol1_sum/(vol1_sum+vol2_sum)
        open(unit=14, file='out1.output/saturation.dat' ,status='unknown',position='append')
        write(14,"(I10,5(1x,e14.7))")ntime,saturation,vol1_sum,vol2_sum,mass1_sum,mass2_sum
        close(14)
        
        temp1=0d0
        temp2=0d0
        temp3=0d0
        temp4=0d0
        !$omp parallel do reduction (+:temp1,temp2,temp3,temp4)
        do k=1,nzGlobal
            temp1 = temp1 + mass1_0(k)
            temp2 = temp2 + mass2_0(k)           
            temp3 = temp3 + vol1_0(k)
            temp4 = temp4 + vol2_0(k)
        enddo      

        saturation_full_domain =temp3/(temp3+temp4)   !full domain saturation

        !********************************************* flowrate calculation **********************************************
        fl1_avg = 0d0
        fl2_avg = 0d0
        fl1_avg_whole = 0d0
        fl2_avg_whole = 0d0          
        !$omp parallel do reduction (+:fl1_avg_whole,fl2_avg_whole)
        do k=1,nzGlobal   !average velocity of whole domain
            fl1_avg_whole = fl1_avg_whole + fl1_0(k)
            fl2_avg_whole = fl2_avg_whole + fl2_0(k)
        enddo      
        fl1_avg_whole = fl1_avg_whole/dble(nzGlobal)
        fl2_avg_whole = fl2_avg_whole/dble(nzGlobal)
        fl_avg_whole = fl1_avg_whole + fl2_avg_whole
        
        !$omp parallel do reduction (+:fl1_avg,fl2_avg)
        do k=n_exclude_inlet+1,nzGlobal-n_exclude_outlet   !average velocity of sample domain          
            fl1_avg = fl1_avg + fl1_0(k)
            fl2_avg = fl2_avg + fl2_0(k)
        enddo      
        fl1_avg = fl1_avg/dble(nzGlobal-n_exclude_outlet-n_exclude_inlet)
        fl2_avg = fl2_avg/dble(nzGlobal-n_exclude_outlet-n_exclude_inlet)
        fl_avg = fl1_avg + fl2_avg
 
        !flowrate based on darcy velocity
        temp = fl_avg/ A_xy  ! averaged darcy velocity for the bulk fluid
        ca = temp * la_nu1 / gamma  !based on bulk fluid velocity

        !********************************************* save data ************************************************************
        open(unit=15, file='out1.output/Ca_number.dat' ,status='unknown',position='append')
        write(15,"(I10,4(1x,e14.7))")ntime,ca,umax_global,kinetic_energy(1),kinetic_energy(2)
        close(15)      

        open(unit=13, file='out1.output/flowrate_time.dat' ,status='unknown',position='append')
        write(13,"(I10,6(1x,E14.6))")ntime,fl_avg_whole,fl1_avg_whole,fl2_avg_whole,fl_avg,fl1_avg,fl2_avg
        close(13)

        open(unit=14, file='out1.output/saturation_full_domain.dat' ,status='unknown',position='append')
        write(14,"(I10,5(1x,e14.7))")ntime,saturation_full_domain,temp3,temp4,temp1,temp2
        close(14)         

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
                write(8,"(5(1x,e14.7))")dble(k),fl1_0(k),fl2_0(k),mass1_0(k)/(mass1_0(k)+mass2_0(k)),pre_0(k)/(dble(pore_profile_z(k))+eps)
            enddo
            close(8)
        endif

        deallocate(fl1_0,fl2_0,mass1_0,mass2_0,pre_0)

        ! check simulation status
        if(ieee_is_nan(saturation_full_domain).or.ieee_is_nan(ca))then
            print*, 'Simulation failed due to NAN of "saturation or capillary number"!'
            simulation_end_indicator = 3
        elseif(umax_global > 0.5)then
            print*, 'Simulation failed due to maximum velocity larger than 0.5!'
            simulation_end_indicator = 3
        endif 
        
        if(steady_state_option==3)then   ! steady state based on saturation
            monitor_previous_value = monitor_current_value  ! store previous step
            monitor_current_value  = saturation_full_domain  
            temp = dabs(monitor_current_value-monitor_previous_value)/(dabs(monitor_current_value)+eps) 
            open(unit=15, file='out1.output/steady_monitor_saturation_error.dat' ,status='unknown',position='append')
            write(15,"(I10,(1x,e14.7))")ntime,temp
            close(15)
            if(temp<convergence_criteria.and.ntime>ntime_monitor)then
                print*, 'Simulation converged based on saturation! Relative error =', temp, '; maximum velocity = ', umax_global
                simulation_end_indicator = 1   !successfully end simulation
            endif
        endif
    endif
    call MPI_Bcast(simulation_end_indicator,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)

    return
end subroutine monitor





!===============================================================================================================================================
!---------------------- monitor_multiphase steady flow  ----------------------
!===============================================================================================================================================
!******************************* monitor_multiphase steady flow - based on phase field *************************************
subroutine monitor_multiphase_steady_phasefield
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    USE,INTRINSIC :: IEEE_ARITHMETIC
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k
    integer (kind=1) :: wall_indicator
    real(kind=8) :: tmp1, tmp2, tmp
    real(kind=8) :: umax, d_phi_max, d_phi_max_global, fx,fy,fz

    umax=0d0
    d_phi_max=0d0

    call compute_macro_vars

    !$omp parallel do private(i,j,tmp1,tmp2,wall_indicator) reduction(max:umax,d_phi_max)
    !$acc kernels present(u,v,w,walls,phi)
    !$acc loop reduction(max:umax,d_phi_max)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,k)
                tmp1 = (u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))* (1-wall_indicator)                              
                if(umax<tmp1)then
                    umax=tmp1
                endif
                tmp2 = dabs(phi(i,j,k)-phi_old(i,j,k))* (1-wall_indicator)                              
                if(d_phi_max<tmp2)then
                    d_phi_max=tmp2
                endif               
            enddo
        enddo
    enddo
    !$acc end kernels

    !$omp parallel do private(i,j)
    !$acc kernels present(phi)
    !$acc loop 
    do k=1,nz
        do j=1,ny
            do i=1,nx
                phi_old(i,j,k) = phi(i,j,k)                                  
            enddo
        enddo
    enddo
    !$acc end kernels

    CALL MPI_REDUCE(d_phi_max,d_phi_max_global, 1, MPI_DOUBLE_PRECISION, MPI_max, 0, MPI_COMM_VGRID, ierr)
    CALL MPI_REDUCE(umax,umax_global, 1, MPI_DOUBLE_PRECISION, MPI_max, 0, MPI_COMM_VGRID, ierr)

    if(id.eq.0.and.ntime>ntime_relaxation)then
        if(ieee_is_nan(umax).or.ieee_is_nan(d_phi_max))then
            print*, 'Simulation failed due to NAN!'
            simulation_end_indicator = 3
        else
            umax_global = dsqrt(umax_global)
            open(unit=14, file='out1.output/steady_monitor_max_phi_change.dat' ,status='unknown',position='append')
            write(14,"(I10,2(1x,e14.7))")ntime,d_phi_max_global,umax_global
            close(14) 
            if(umax_global<0.5)then         !maximum velocity smaller than 0.5, otherwise, simulation is considered as failed
                if(d_phi_max_global<convergence_criteria.and.ntime>ntime_monitor)then
                    print*, 'Simulation converged based on maximum local change of phi! Relative error =', d_phi_max_global, '; &
                    & maximum velocity = ', umax_global
                    simulation_end_indicator = 1   !successfully end simulation
                endif
            else
                print*,'Maximum velocity larger than 0.5, simulation failed!!!'
                simulation_end_indicator = 3   !stop simulation due to failure
            endif     
        endif         
    endif

    call MPI_Bcast(simulation_end_indicator,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)

    return
end subroutine monitor_multiphase_steady_phasefield
!******************************* monitor_multiphase steady flow - based on phase field *************************************


!*************************** monitor_multiphase steady flow - based on capillary pressure **********************************
subroutine monitor_multiphase_steady_capillarypressure
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    USE,INTRINSIC :: IEEE_ARITHMETIC
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k,i_w,i_nw,i_w_sum,i_nw_sum
    integer (kind=1) :: wall_indicator
    real(kind=8) :: umax,temp, pre_w, pre_nw, pre_w_sum, pre_nw_sum, pc

    call compute_macro_vars

    !$omp parallel do private(i,j,temp,wall_indicator) reduction(max:umax)
    !$acc kernels present(u,v,w,walls)
    !$acc loop reduction(max:umax)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,k)
                temp = (u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))* (1-wall_indicator)
                if(umax<temp)then
                    umax=temp
                endif
            enddo
        enddo
    enddo
    !$acc end kernels

    pre_w=0d0
    pre_nw=0d0
    i_w=0
    i_nw=0
    !$OMP parallel DO private(i,j) reduction(+:pre_w,pre_nw,i_w,i_nw)
    !$acc kernels  present(phi,rho,walls)
    !$acc loop reduction(+:pre_w,pre_nw,i_w,i_nw) collapse(3)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                if(walls(i,j,k)==0)then
                    if(phi(i,j,k)<-0.99d0)then
                        pre_w = pre_w + rho(i,j,k)
                        i_w = i_w +1
                    endif
                    if(phi(i,j,k)>0.99d0)then
                        pre_nw = pre_nw + rho(i,j,k)
                        i_nw = i_nw +1
                    endif
                endif
            enddo
        enddo
    enddo
    !$acc end kernels

    CALL MPI_REDUCE(pre_w,pre_w_sum, 1, MPI_DOUBLE_PRECISION, MPI_sum, 0, MPI_COMM_VGRID, ierr)
    CALL MPI_REDUCE(pre_nw,pre_nw_sum, 1, MPI_DOUBLE_PRECISION, MPI_sum, 0, MPI_COMM_VGRID, ierr)
    CALL MPI_REDUCE(i_w,i_w_sum, 1, MPI_INTEGER, MPI_sum, 0, MPI_COMM_VGRID, ierr)
    CALL MPI_REDUCE(i_nw,i_nw_sum, 1, MPI_INTEGER, MPI_sum, 0, MPI_COMM_VGRID, ierr)
    CALL MPI_REDUCE(umax,umax_global, 1, MPI_DOUBLE_PRECISION, MPI_max, 0, MPI_COMM_VGRID, ierr)

    if(id.eq.0)then
        umax_global = dsqrt(umax_global)
        pre_w = pre_w_sum/(i_w_sum+eps)
        pre_nw = pre_nw_sum/(i_nw_sum+eps)
        pc=(pre_nw-pre_w)/3d0
        monitor_previous_value = monitor_current_value  ! store previous step
        monitor_current_value  = pc        
        temp = dabs(monitor_current_value-monitor_previous_value)/(dabs(monitor_current_value)+eps)
        open(unit=15, file='out1.output/steady_monitor_capillary_pressure_error.dat' ,status='unknown',position='append')
        write(15,"(I10,6(1x,e14.7))")ntime,temp,pc,pre_w,pre_nw,temp,umax_global
        close(15)

        if(ieee_is_nan(umax).or.ieee_is_nan(pc))then
            print*, 'Simulation failed due to NAN!'
            simulation_end_indicator = 3
        else
            if(umax_global<0.5)then         !maximum velocity smaller than 0.5, otherwise, simulation is considered as failed
                if(temp<convergence_criteria.and.ntime>ntime_monitor)then
                    print*, 'Simulation converged based on change of capillary pressure! Relative error =', temp, '; maximum velocity = ', umax_global
                    simulation_end_indicator = 1   !successfully end simulation
                endif
            else
                print*,'maximum velocity larger than 0.5, simulation failed!!!'
                simulation_end_indicator = 3   !stop simulation due to failure
            endif
        endif
    endif

    call MPI_Bcast(simulation_end_indicator,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)

    return
end subroutine monitor_multiphase_steady_capillarypressure
!******************************* monitor_multiphase steady flow - based on capillary pressure *************************************




!============================================================================================================================================
!---------------------- misc subroutines  ----------------------
!============================================================================================================================================
!*******************************    monitor_breakthrough    *************************************
subroutine monitor_breakthrough
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k,obs_z,itemp,z

    itemp=0
    if(idz==npz-1)then
        obs_z = nz-1   !observation position
        !$OMP PARALLEL DO reduction(+:itemp)
        !$acc kernels present(phi,walls)
        !$acc loop reduction(+:itemp)
        do j=1,ny
            do i=1,nx
                if(walls(i,j,obs_z)==0.and.phi(i,j,obs_z)>0d0)then
                    itemp = itemp + 1
                endif
            enddo
        enddo
        !$acc end kernels
    endif

    call MPI_Reduce(itemp,outlet_phase1_sum,1,MPI_integer,MPI_SUM,0,MPI_COMM_VGRID,ierr)

    if(id==0.and.outlet_phase1_sum>=1)then        !break through detect fluid 1 only
        simulation_end_indicator = 1
        print*,'Breakthrough point reached! Exiting program!'
    endif

    call MPI_Bcast(simulation_end_indicator,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)

    return
end subroutine monitor_breakthrough
!*******************************    monitor_breakthrough    *************************************


!*******************************    calculate saturationm   *************************************
subroutine cal_saturation
    use Misc_module
    use Fluid_singlephase
    use Fluid_multiphase
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k
    integer (kind=1) :: wall_indicator
    real(kind=8) :: ft0,ft1,ft2,ft3,ft4,ft5,ft6,ft7,ft8,ft9,ft10,ft11,ft12,ft13,ft14,ft15,ft16,ft17,ft18,v1,v2

    v1=0d0
    v2=0d0
    
    !$omp parallel do private(i,j,wall_indicator)reduction(+:v1,v2)
    !$acc kernels present(phi,rho,walls)
    !$acc loop reduction(+:v1,v2) collapse(3)
    do k=1,nz      
        do j=1,ny
            do i=1,nx
                wall_indicator = walls(i,j,k)
                v1= v1 + 0.5d0*(1d0+phi(i,j,k))* (1-wall_indicator)
                v2= v2 + 0.5d0*(1d0-phi(i,j,k))* (1-wall_indicator)               
            enddo
        enddo       
    enddo
    !$acc end kernels

    CALL MPI_REDUCE(v1,vol1_sum, 1, MPI_DOUBLE_PRECISION, MPI_sum, 0, MPI_COMM_VGRID, ierr)
    CALL MPI_REDUCE(v2,vol2_sum, 1, MPI_DOUBLE_PRECISION, MPI_sum, 0, MPI_COMM_VGRID, ierr)

    if(id==0)then
        saturation_full_domain = vol1_sum/(vol1_sum+vol2_sum+eps)
    endif

    call MPI_Bcast(saturation_full_domain,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_VGRID,ierr)

    return
end subroutine cal_saturation
!*******************************    calculate saturationm   *************************************
