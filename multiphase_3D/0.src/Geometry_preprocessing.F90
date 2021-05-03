#include "./preprocessor.h"

!=====================================================================================================================================
!---------------------- geometry_preprocessing new ----------------------
! process the geometry before the start of the main iteration: normal directions of the solid surface; extrapolation weights from 
! different directions not recommended for relatively large domain - the processing time before each new or old simulation could be 
! too long (not paralllelized on CPU for the GPU version code); also not optimal for Xeon Phi as the on-package memory is limited
!=====================================================================================================================================
subroutine geometry_preprocessing_new
    use Misc_module
    use Fluid_multiphase, only: theta, ISO8   !contact angle 
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,L,M,N,k,icount,out1,out2,out3,num,icount1,icount2,ncount,iteration,overlap_temp
    real(kind=8) :: tmp
    real(kind=8) :: nwx,nwy,nwz 
    real(kind=8), allocatable, dimension(:,:,:) :: walls_smooth_1, walls_smooth_2   
    integer(kind=1), allocatable, dimension(:,:,:) :: walls_temp1   
    integer :: i_max_iteration, ghost_layers
    integer, dimension(:), parameter ::  iex(0:26)=(/0, 1, -1,  0,  0,  0,  0, 1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1 /)
    integer, dimension(:), parameter ::  iey(0:26)=(/0, 0,  0,  1, -1,  0,  0, 1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1, -1,  1, -1,  1 /)
    integer, dimension(:), parameter ::  iez(0:26)=(/0, 0,  0,  0,  0,  1, -1, 0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1, -1,  1,  1, -1, -1,  1,  1, -1 /)
    real(kind=8), parameter, dimension(:) :: we(0:3) = (/8.0d0/27.0d0, 2.0d0/27.0d0, 1.0d0/54.0d0, 1.0d0/216.0d0/)

    !without contact angle information, used to load geometry data
    type indirect_fluid_boundary_nodes_tmp                   
        integer :: ix,iy,iz
        real(kind=8) :: nwx,nwy,nwz 
    end type indirect_fluid_boundary_nodes_tmp
    type(indirect_solid_boundary_nodes), allocatable, dimension(:) :: solid_boundary_nodes_global
    type(indirect_fluid_boundary_nodes_tmp), allocatable, dimension(:) :: fluid_boundary_nodes_global
    type(indirect_solid_boundary_nodes) :: sdummy
    type(indirect_fluid_boundary_nodes_tmp) :: fdummy
    integer, dimension(:) ::blocklengths(3), types(3)
    integer(KIND=MPI_ADDRESS_KIND)  :: extent,lb ,base ,displacements(3)
    integer :: mpi_solid_bc_type,mpi_fluid_bc_type

    if(id==0)then
        i_max_iteration = 4  !smoothing iterations
        ghost_layers = 6  !used for smoothing iterations
        ghost_layers =  ghost_layers + overlap_phi   ! additional layers used for solid/fluid boundary nodes that outside the domain
        allocate(walls_smooth_1(1-ghost_layers:nxGlobal+ghost_layers,1-ghost_layers:nyGlobal+ghost_layers,1-ghost_layers:nzGlobal+ghost_layers))
        allocate(walls_smooth_2(1-ghost_layers:nxGlobal+ghost_layers,1-ghost_layers:nyGlobal+ghost_layers,1-ghost_layers:nzGlobal+ghost_layers))   
        allocate(walls_temp1(1-ghost_layers:nxGlobal+ghost_layers,1-ghost_layers:nyGlobal+ghost_layers,1-ghost_layers:nzGlobal+ghost_layers))

        !$OMP PARALLEL DO private(i,j) 
        do k=1,nzGlobal
            do j=1,nyGlobal
                do i=1,nxGlobal    
                    walls_temp1(i,j,k) = walls_global(i,j,k)
                enddo
            enddo
        enddo   

        if(kper==0)then  !z direction, not periodic       
            !$OMP PARALLEL DO private(i)
            do j=1,nyGlobal
                do i=1,nxGlobal                      
                    walls_temp1(i,j,1-ghost_layers:0)=  walls_temp1(i,j,1)                  
                    walls_temp1(i,j,nzglobal+1:nzGlobal+ghost_layers)=  walls_temp1(i,j,nzglobal)               
                enddo
            enddo 
        else            !z direction,periodic       
            !$OMP PARALLEL DO private(i)
            do j=1,nyGlobal
                do i=1,nxGlobal                      
                    walls_temp1(i,j,1-ghost_layers:0) =  walls_temp1(i,j,nzglobal-ghost_layers+1:nzglobal)                  
                    walls_temp1(i,j,nzglobal+1:nzGlobal+ghost_layers)=  walls_temp1(i,j,1:ghost_layers)               
                enddo
            enddo 
        endif
        
        if(jper==0)then  !y direction, not periodic       
            !$OMP PARALLEL DO private(i)
            do k=1-ghost_layers,nzGlobal+ghost_layers
                do i=1,nxGlobal                       
                    walls_temp1(i,1-ghost_layers:0,k)=  walls_temp1(i,1,k)                    
                    walls_temp1(i,nyglobal+1:nyGlobal+ghost_layers,k)=  walls_temp1(i,nyGlobal,k)               
                enddo
            enddo 
        else        !y direction,periodic  
            !$OMP PARALLEL DO private(i)   
            do k=1-ghost_layers,nzGlobal+ghost_layers
                do i=1,nxGlobal                       
                    walls_temp1(i,1-ghost_layers:0,k)=  walls_temp1(i,nyglobal-ghost_layers+1:nyglobal,k)                    
                    walls_temp1(i,nyglobal+1:nyGlobal+ghost_layers,k)=  walls_temp1(i,1:ghost_layers,k)               
                enddo
            enddo 
        endif  

        if(iper==0)then  !x direction, not periodic  
            !$OMP PARALLEL DO private(j)
            do k=1-ghost_layers,nzGlobal+ghost_layers
                do j=1-ghost_layers,nyGlobal+ghost_layers                       
                    walls_temp1(1-ghost_layers:0,j,k)=  walls_temp1(1,j,k)                    
                    walls_temp1(nxglobal+1:nxGlobal+ghost_layers,j,k)=  walls_temp1(nxGlobal,j,k)               
                enddo
            enddo
        else          !x direction, periodic  
            !$OMP PARALLEL DO private(j)
            do k=1-ghost_layers,nzGlobal+ghost_layers
                do j=1-ghost_layers,nyGlobal+ghost_layers                       
                    walls_temp1(1-ghost_layers:0,j,k)=  walls_temp1(nxglobal-ghost_layers+1:nxglobal,j,k)                    
                    walls_temp1(nxglobal+1:nxGlobal+ghost_layers,j,k)=  walls_temp1(1:ghost_layers,j,k)               
                enddo
            enddo
        endif

        !$OMP PARALLEL DO private(i,j) 
        do k=1-ghost_layers,nzGlobal+ghost_layers
            do j=1-ghost_layers,nyGlobal+ghost_layers
                do i=1-ghost_layers,nxGlobal+ghost_layers  
                    walls_smooth_1(i,j,k) =  walls_temp1(i,j,k)        
                    walls_smooth_2(i,j,k) =  walls_temp1(i,j,k)                   
                enddo
            enddo
        enddo   

        !extend wall information
        !$OMP PARALLEL DO private(i,j,n) 
        do k=1-ghost_layers+1,nzGlobal+ghost_layers-1
            do j=1-ghost_layers+1,nyGlobal+ghost_layers-1
                do i=1-ghost_layers+1,nxGlobal+ghost_layers-1
                    if(walls_temp1(i,j,k)==1)then
                        do n=1,18
                            if(walls_temp1(i+ex(n),j+ey(n),k+ez(n))<=0)then
                                walls_temp1(i,j,k)=2   !solid boundary nodes
                                exit
                            endif
                        enddo
                    endif   
                    if(walls_temp1(i,j,k)==0)then
                        do n=1,18
                            if(walls_temp1(i+ex(n),j+ey(n),k+ez(n))>=1)then
                                walls_temp1(i,j,k)=-1   !fluid boundary nodes
                                exit
                            endif
                        enddo
                    endif 
                enddo
            enddo
        enddo 

        do iteration=1,i_max_iteration
            !$OMP PARALLEL 
            !$OMP do private(i,j,n,m)   
            do k=1-ghost_layers+1,nzGlobal+ghost_layers-1
                do j=1-ghost_layers+1,nyGlobal+ghost_layers-1
                    do i=1-ghost_layers+1,nxGlobal+ghost_layers-1
                        walls_smooth_2(i,j,k)= 0d0 
                        do n=0,26   
                            m = iex(n)*iex(n) + iey(n)*iey(n) + iez(n)*iez(n)      
                            walls_smooth_2(i,j,k)= walls_smooth_2(i,j,k) + walls_smooth_1(i+iex(n),j+iey(n),k+iez(n))*we(m)
                        enddo
                    enddo
                enddo
            enddo
        
            !$OMP do private(i,j)
            do k=1-ghost_layers+1,nzGlobal+ghost_layers-1
                do j=1-ghost_layers+1,nyGlobal+ghost_layers-1
                    do i=1-ghost_layers+1,nxGlobal+ghost_layers-1             
                        walls_smooth_1(i,j,k)= walls_smooth_2(i,j,k) 
                    enddo
                enddo
            enddo
            !$OMP end PARALLEL  
        enddo

        num_solid_boundary_global = 0 
        num_fluid_boundary_global = 0
        !$omp parallel do private(i,j) reduction (+:num_solid_boundary_global,num_fluid_boundary_global)
        do k=1-overlap_phi,nzGlobal+overlap_phi       !smoothing process includs overlap_phi region used for cnx,cny,cnz calcaulation
            do j=1-overlap_phi,nyGlobal+overlap_phi
                do i=1-overlap_phi,nxGlobal+overlap_phi
                    if(walls_temp1(i,j,k)==2)then
                        num_solid_boundary_global = num_solid_boundary_global + 1
                    endif
                    if(walls_temp1(i,j,k)==-1)then
                        num_fluid_boundary_global = num_fluid_boundary_global + 1
                    endif
                enddo
            enddo
        enddo
        print*, "Total number of solid boundary nodes: ", num_solid_boundary_global
        print*, "Total number of fluid boundary nodes: ", num_fluid_boundary_global
    endif

    call MPI_Bcast(num_solid_boundary_global,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)   
    call MPI_Bcast(num_fluid_boundary_global,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)   
    allocate(solid_boundary_nodes_global(num_solid_boundary_global))  
    allocate(fluid_boundary_nodes_global(num_fluid_boundary_global))  

    if(id==0)then
        icount1 = 0
        icount2 = 0
        do k=1-overlap_phi,nzGlobal+overlap_phi            
            do j=1-overlap_phi,nyGlobal+overlap_phi
                do i=1-overlap_phi,nxGlobal+overlap_phi
                    if(walls_temp1(i,j,k)==2)then
                        icount1 = icount1 + 1
                        solid_boundary_nodes_global(icount1)%ix = i
                        solid_boundary_nodes_global(icount1)%iy = j
                        solid_boundary_nodes_global(icount1)%iz = k
                        ncount = 0
                        solid_boundary_nodes_global(icount1)%la_weight = 0d0
                        do n=1,18
                            if(walls_temp1(i+ex(n),j+ey(n),k+ez(n))<=0)then
                                ncount= ncount + 1
                                solid_boundary_nodes_global(icount1)%la_weight = solid_boundary_nodes_global(icount1)%la_weight + w_equ(n)                  
                                solid_boundary_nodes_global(icount1)%neighbor_list(ncount)=n                    
                            endif
                        enddo
                        solid_boundary_nodes_global(icount1)%i_fluid_num = ncount  !number of fluid nodes connected to the boundary node                             
                    endif
                    if(walls_temp1(i,j,k)==-1)then
                        icount2 = icount2 + 1
                        fluid_boundary_nodes_global(icount2)%ix = i
                        fluid_boundary_nodes_global(icount2)%iy = j
                        fluid_boundary_nodes_global(icount2)%iz = k                        
                    endif
                enddo
            enddo
        enddo  

        !$omp parallel do private (i,j,k,nwx,nwy,nwz,tmp)
        do  num=1,num_fluid_boundary_global
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cx~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            i = fluid_boundary_nodes_global(num)%ix
            j = fluid_boundary_nodes_global(num)%iy
            k = fluid_boundary_nodes_global(num)%iz
                    
            nwx = &
                    
                ISO8(1)*(walls_smooth_2(i+1,j  ,k  ) - walls_smooth_2(i-1,j  ,k  )) + &

                ISO8(2)*(&
                walls_smooth_2(i+1,j+1,k  ) - walls_smooth_2(i-1,j-1,k  ) + &
                walls_smooth_2(i+1,j-1,k  ) - walls_smooth_2(i-1,j+1,k  ) + &
                walls_smooth_2(i+1,j  ,k+1) - walls_smooth_2(i-1,j  ,k-1) + &
                walls_smooth_2(i+1,j  ,k-1) - walls_smooth_2(i-1,j  ,k+1)) + &

                ISO8(3)*(&
                walls_smooth_2(i+1,j+1,k+1) - walls_smooth_2(i-1,j-1,k-1) + &
                walls_smooth_2(i+1,j+1,k-1) - walls_smooth_2(i-1,j-1,k+1) + &
                walls_smooth_2(i+1,j-1,k+1) - walls_smooth_2(i-1,j+1,k-1) + &
                walls_smooth_2(i+1,j-1,k-1) - walls_smooth_2(i-1,j+1,k+1)) + &

                ISO8(4)*(walls_smooth_2(i+2,j  ,k  ) - walls_smooth_2(i-2,j  ,k  )) + &

                ISO8(5)*(&
                walls_smooth_2(i+2,j+1,k  ) - walls_smooth_2(i-2,j-1,k  ) + &
                walls_smooth_2(i+2,j-1,k  ) - walls_smooth_2(i-2,j+1,k  ) + &
                walls_smooth_2(i+2,j  ,k+1) - walls_smooth_2(i-2,j  ,k-1) + &
                walls_smooth_2(i+2,j  ,k-1) - walls_smooth_2(i-2,j  ,k+1) + &
                walls_smooth_2(i+1,j+2,k  ) - walls_smooth_2(i-1,j-2,k  ) + &
                walls_smooth_2(i+1,j-2,k  ) - walls_smooth_2(i-1,j+2,k  ) + &
                walls_smooth_2(i+1,j  ,k+2) - walls_smooth_2(i-1,j  ,k-2) + &
                walls_smooth_2(i+1,j  ,k-2) - walls_smooth_2(i-1,j  ,k+2)) + &

                ISO8(6)*(&
                walls_smooth_2(i+2,j+1,k+1) - walls_smooth_2(i-2,j-1,k-1) + &
                walls_smooth_2(i+2,j+1,k-1) - walls_smooth_2(i-2,j-1,k+1) + &
                walls_smooth_2(i+2,j-1,k+1) - walls_smooth_2(i-2,j+1,k-1) + &
                walls_smooth_2(i+2,j-1,k-1) - walls_smooth_2(i-2,j+1,k+1) + &
                walls_smooth_2(i+1,j+2,k+1) - walls_smooth_2(i-1,j-2,k-1) + &
                walls_smooth_2(i+1,j+2,k-1) - walls_smooth_2(i-1,j-2,k+1) + &
                walls_smooth_2(i+1,j-2,k+1) - walls_smooth_2(i-1,j+2,k-1) + &
                walls_smooth_2(i+1,j-2,k-1) - walls_smooth_2(i-1,j+2,k+1) + &
                walls_smooth_2(i+1,j+1,k+2) - walls_smooth_2(i-1,j-1,k-2) + &
                walls_smooth_2(i+1,j+1,k-2) - walls_smooth_2(i-1,j-1,k+2) + &
                walls_smooth_2(i+1,j-1,k+2) - walls_smooth_2(i-1,j+1,k-2) + &
                walls_smooth_2(i+1,j-1,k-2) - walls_smooth_2(i-1,j+1,k+2)) + &

                ISO8(7)*(&
                walls_smooth_2(i+2,j+2,k  ) - walls_smooth_2(i-2,j-2,k  ) + &
                walls_smooth_2(i+2,j-2,k  ) - walls_smooth_2(i-2,j+2,k  ) + &
                walls_smooth_2(i+2,j  ,k+2) - walls_smooth_2(i-2,j  ,k-2) + &
                walls_smooth_2(i+2,j  ,k-2) - walls_smooth_2(i-2,j  ,k+2))

                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cy~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            nwy = &
                    
                ISO8(1)*(walls_smooth_2(i  ,j+1,k  ) - walls_smooth_2(i  ,j-1,k  )) + &

                ISO8(2)*(&
                walls_smooth_2(i+1,j+1,k  ) - walls_smooth_2(i-1,j-1,k  ) + &
                walls_smooth_2(i-1,j+1,k  ) - walls_smooth_2(i+1,j-1,k  ) + &
                walls_smooth_2(i  ,j+1,k+1) - walls_smooth_2(i  ,j-1,k-1) + &
                walls_smooth_2(i  ,j+1,k-1) - walls_smooth_2(i  ,j-1,k+1)) + &

                ISO8(3)*(&
                walls_smooth_2(i+1,j+1,k+1) - walls_smooth_2(i-1,j-1,k-1) + &
                walls_smooth_2(i+1,j+1,k-1) - walls_smooth_2(i-1,j-1,k+1) + &
                walls_smooth_2(i-1,j+1,k-1) - walls_smooth_2(i+1,j-1,k+1) + &
                walls_smooth_2(i-1,j+1,k+1) - walls_smooth_2(i+1,j-1,k-1)) + &

                ISO8(4)*(walls_smooth_2(i  ,j+2,k  ) - walls_smooth_2(i  ,j-2,k  )) + &

                ISO8(5)*(&
                walls_smooth_2(i+2,j+1,k  ) - walls_smooth_2(i-2,j-1,k  ) + &
                walls_smooth_2(i-2,j+1,k  ) - walls_smooth_2(i+2,j-1,k  ) + &
                walls_smooth_2(i  ,j+2,k+1) - walls_smooth_2(i  ,j-2,k-1) + &
                walls_smooth_2(i  ,j+2,k-1) - walls_smooth_2(i  ,j-2,k+1) + &
                walls_smooth_2(i+1,j+2,k  ) - walls_smooth_2(i-1,j-2,k  ) + &
                walls_smooth_2(i-1,j+2,k  ) - walls_smooth_2(i+1,j-2,k  ) + &
                walls_smooth_2(i  ,j+1,k+2) - walls_smooth_2(i  ,j-1,k-2) + &
                walls_smooth_2(i  ,j+1,k-2) - walls_smooth_2(i  ,j-1,k+2)) + &

                ISO8(6)*(&
                walls_smooth_2(i+2,j+1,k+1) - walls_smooth_2(i-2,j-1,k-1) + &
                walls_smooth_2(i+2,j+1,k-1) - walls_smooth_2(i-2,j-1,k+1) + &
                walls_smooth_2(i-2,j+1,k+1) - walls_smooth_2(i+2,j-1,k-1) + &
                walls_smooth_2(i-2,j+1,k-1) - walls_smooth_2(i+2,j-1,k+1) + &
                walls_smooth_2(i+1,j+2,k+1) - walls_smooth_2(i-1,j-2,k-1) + &
                walls_smooth_2(i+1,j+2,k-1) - walls_smooth_2(i-1,j-2,k+1) + &
                walls_smooth_2(i-1,j+2,k+1) - walls_smooth_2(i+1,j-2,k-1) + &
                walls_smooth_2(i-1,j+2,k-1) - walls_smooth_2(i+1,j-2,k+1) + &
                walls_smooth_2(i+1,j+1,k+2) - walls_smooth_2(i-1,j-1,k-2) + &
                walls_smooth_2(i+1,j+1,k-2) - walls_smooth_2(i-1,j-1,k+2) + &
                walls_smooth_2(i-1,j+1,k+2) - walls_smooth_2(i+1,j-1,k-2) + &
                walls_smooth_2(i-1,j+1,k-2) - walls_smooth_2(i+1,j-1,k+2)) + &

                ISO8(7)*(&
                walls_smooth_2(i+2,j+2,k  ) - walls_smooth_2(i-2,j-2,k  ) + &
                walls_smooth_2(i-2,j+2,k  ) - walls_smooth_2(i+2,j-2,k  ) + &
                walls_smooth_2(i  ,j+2,k+2) - walls_smooth_2(i  ,j-2,k-2) + &
                walls_smooth_2(i  ,j+2,k-2) - walls_smooth_2(i  ,j-2,k+2))
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cz~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            nwz = &
                    
                ISO8(1)*(walls_smooth_2(i  ,j  ,k+1) -walls_smooth_2(i  ,j  ,k-1)) + &

                ISO8(2)*(&
                walls_smooth_2(i  ,j+1,k+1) - walls_smooth_2(i  ,j-1,k-1) + &
                walls_smooth_2(i  ,j-1,k+1) - walls_smooth_2(i  ,j+1,k-1) + &
                walls_smooth_2(i+1,j  ,k+1) - walls_smooth_2(i-1,j  ,k-1) + &
                walls_smooth_2(i-1,j  ,k+1) - walls_smooth_2(i+1,j  ,k-1)) + &

                +ISO8(3)*(&
                walls_smooth_2(i+1,j+1,k+1) - walls_smooth_2(i-1,j-1,k-1) + &
                walls_smooth_2(i+1,j-1,k+1) - walls_smooth_2(i-1,j+1,k-1) + &
                walls_smooth_2(i-1,j+1,k+1) - walls_smooth_2(i+1,j-1,k-1) + &
                walls_smooth_2(i-1,j-1,k+1) - walls_smooth_2(i+1,j+1,k-1)) + &

                ISO8(4)*(walls_smooth_2(i  ,j  ,k+2) - walls_smooth_2(i  ,j  ,k-2)) + &

                ISO8(5)*(&
                walls_smooth_2(i  ,j+1,k+2) - walls_smooth_2(i  ,j-1,k-2) + &
                walls_smooth_2(i  ,j-1,k+2) - walls_smooth_2(i  ,j+1,k-2) + &
                walls_smooth_2(i+2,j  ,k+1) - walls_smooth_2(i-2,j  ,k-1) + &
                walls_smooth_2(i-2,j  ,k+1) - walls_smooth_2(i+2,j  ,k-1) + &
                walls_smooth_2(i  ,j+2,k+1) - walls_smooth_2(i  ,j-2,k-1) + &
                walls_smooth_2(i  ,j-2,k+1) - walls_smooth_2(i  ,j+2,k-1) + &
                walls_smooth_2(i+1,j  ,k+2) - walls_smooth_2(i-1,j  ,k-2) + &
                walls_smooth_2(i-1,j  ,k+2) - walls_smooth_2(i+1,j  ,k-2)) + &

                ISO8(6)*(&
                walls_smooth_2(i+2,j+1,k+1) - walls_smooth_2(i-2,j-1,k-1) + &
                walls_smooth_2(i+2,j-1,k+1) - walls_smooth_2(i-2,j+1,k-1) + &
                walls_smooth_2(i-2,j+1,k+1) - walls_smooth_2(i+2,j-1,k-1) + &
                walls_smooth_2(i-2,j-1,k+1) - walls_smooth_2(i+2,j+1,k-1) + &
                walls_smooth_2(i+1,j+2,k+1) - walls_smooth_2(i-1,j-2,k-1) + &
                walls_smooth_2(i+1,j-2,k+1) - walls_smooth_2(i-1,j+2,k-1) + &
                walls_smooth_2(i-1,j+2,k+1) - walls_smooth_2(i+1,j-2,k-1) + &
                walls_smooth_2(i-1,j-2,k+1) - walls_smooth_2(i+1,j+2,k-1) + &
                walls_smooth_2(i+1,j+1,k+2) - walls_smooth_2(i-1,j-1,k-2) + &
                walls_smooth_2(i+1,j-1,k+2) - walls_smooth_2(i-1,j+1,k-2) + &
                walls_smooth_2(i-1,j+1,k+2) - walls_smooth_2(i+1,j-1,k-2) + &
                walls_smooth_2(i-1,j-1,k+2) - walls_smooth_2(i+1,j+1,k-2)) + &

                ISO8(7)*(&
                walls_smooth_2(i  ,j+2,k+2) - walls_smooth_2(i  ,j-2,k-2) + &
                walls_smooth_2(i  ,j-2,k+2) - walls_smooth_2(i  ,j+2,k-2) + &
                walls_smooth_2(i+2,j  ,k+2) - walls_smooth_2(i-2,j  ,k-2) + &
                walls_smooth_2(i-2,j  ,k+2) - walls_smooth_2(i+2,j  ,k-2))
                        
            tmp = 1d0/(dsqrt (nwx*nwx+nwy*nwy+nwz*nwz) + eps)          !normalized vector
            fluid_boundary_nodes_global(num)%nwx = nwx*tmp
            fluid_boundary_nodes_global(num)%nwy = nwy*tmp
            fluid_boundary_nodes_global(num)%nwz = nwz*tmp
        enddo
        
        deallocate(walls_temp1, walls_smooth_1, walls_smooth_2)
    endif
    
    blocklengths(1) = 4
    blocklengths(2) = 18
    blocklengths(3) = 1
    types(1) = MPI_INTEGER
    types(2) = MPI_INTEGER
    types(3) = MPI_double_precision
    call MPI_GET_ADDRESS(sdummy%ix, displacements(1), ierr) 
    call MPI_GET_ADDRESS(sdummy%neighbor_list(1), displacements(2), ierr)
    call MPI_GET_ADDRESS(sdummy%la_weight, displacements(3), ierr)
    base = displacements(1)
    displacements(1) = displacements(1) - base
    displacements(2) = displacements(2) - base
    displacements(3) = displacements(3) - base
    call MPI_Type_create_struct(3,blocklengths,displacements,types,mpi_solid_bc_type,ierr)
    call MPI_TYPE_COMMIT(mpi_solid_bc_type,ierr)
    ! call MPI_TYPE_get_EXTENT(mpi_solid_bc_type,lb,extent,ierr)
    ! call MPI_Type_size(mpi_solid_bc_type,i,ierr)
    ! print*, extent,'s-type'

    blocklengths(1) = 3
    blocklengths(2) = 3
    types(1) = MPI_INTEGER
    types(2) = MPI_double_precision
    call MPI_GET_ADDRESS(fdummy%ix, displacements(1), ierr) 
    call MPI_GET_ADDRESS(fdummy%nwx, displacements(2), ierr) 
    base = displacements(1)
    displacements(1) = displacements(1) - base
    displacements(2) = displacements(2) - base
    call MPI_Type_create_struct(2,blocklengths,displacements,types,mpi_fluid_bc_type,ierr)
    call MPI_TYPE_COMMIT(mpi_fluid_bc_type,ierr)
    !call MPI_TYPE_get_EXTENT(mpi_fluid_bc_type,lb,extent,ierr)
    !print*, extent,'f-type'

    call MPI_Bcast(solid_boundary_nodes_global,num_solid_boundary_global,mpi_solid_bc_type,0,MPI_COMM_VGRID,ierr)
    call MPI_Bcast(fluid_boundary_nodes_global,num_fluid_boundary_global,mpi_fluid_bc_type,0,MPI_COMM_VGRID,ierr) 

    CALL MPI_CART_COORDS(MPI_COMM_VGRID, id, mpi_dim, mpi_coords, ierr)
    out1=mpi_coords(1)
    out2=mpi_coords(2)
    out3=mpi_coords(3)

    icount=0
    overlap_temp = 3
    !$OMP parallel do private(l,m,n,i,j,k) reduction (+:icount)
    do num=1, num_solid_boundary_global
        l= solid_boundary_nodes_global(num)%ix
        m= solid_boundary_nodes_global(num)%iy
        n= solid_boundary_nodes_global(num)%iz
        i = l-out1*nx
        j = m-out2*ny
        k = n-out3*nz
        !three layers of ghost nodes, used in CSF calculation
        if(i>=1-overlap_temp.and.i<=nx+overlap_temp.and.j>=1-overlap_temp.and.j<=ny+overlap_temp.and.k>=1-overlap_temp.and.k<=nz+overlap_temp)then    
            ! if(i>=1.and.i<=nx.and.j>=1.and.j<=ny.and.k>=1.and.k<=nz)then 
            icount = icount + 1
        endif
    enddo  
    num_solid_boundary = icount
    allocate(solid_boundary_nodes(num_solid_boundary))
    !print*, "total number of solid boundary nodes in local MPI domain",id,"= ", num_solid_boundary
    icount=0
    overlap_temp = 2
    !$OMP parallel do private(l,m,n,i,j,k) reduction (+:icount)
    do num=1, num_fluid_boundary_global
        l= fluid_boundary_nodes_global(num)%ix
        m= fluid_boundary_nodes_global(num)%iy
        n= fluid_boundary_nodes_global(num)%iz
        i = l-out1*nx
        j = m-out2*ny
        k = n-out3*nz 
        if(i>=1-overlap_temp.and.i<=nx+overlap_temp.and.j>=1-overlap_temp.and.j<=ny+overlap_temp.and.k>=1-overlap_temp.and.k<=nz+overlap_temp)then
            icount = icount + 1
        endif
    enddo 
    num_fluid_boundary = icount
    allocate(fluid_boundary_nodes(num_fluid_boundary))
    !print*, "total number of fluid boundary nodes in local MPI domain",id,"= ", num_fluid_boundary

    overlap_temp = 3
    icount=0  
    do num=1, num_solid_boundary_global
        l= solid_boundary_nodes_global(num)%ix
        m= solid_boundary_nodes_global(num)%iy
        n= solid_boundary_nodes_global(num)%iz
        i = l-out1*nx
        j = m-out2*ny
        k = n-out3*nz
        !three layer of ghost nodes, used in CSF calculation
        if(i>=1-overlap_temp.and.i<=nx+overlap_temp.and.j>=1-overlap_temp.and.j<=ny+overlap_temp.and.k>=1-overlap_temp.and.k<=nz+overlap_temp)then   
            icount = icount + 1
            solid_boundary_nodes(icount)%ix = i
            solid_boundary_nodes(icount)%iy = j
            solid_boundary_nodes(icount)%iz = k
            solid_boundary_nodes(icount)%i_fluid_num = solid_boundary_nodes_global(num)%i_fluid_num
            solid_boundary_nodes(icount)%la_weight = solid_boundary_nodes_global(num)%la_weight
            solid_boundary_nodes(icount)%neighbor_list(:)=solid_boundary_nodes_global(num)%neighbor_list(:)
        endif
    enddo  

    icount=0
    overlap_temp = 2  
    do num=1, num_fluid_boundary_global
        l= fluid_boundary_nodes_global(num)%ix
        m= fluid_boundary_nodes_global(num)%iy
        n= fluid_boundary_nodes_global(num)%iz
        i = l-out1*nx
        j = m-out2*ny
        k = n-out3*nz

        if(i>=1-overlap_temp.and.i<=nx+overlap_temp.and.j>=1-overlap_temp.and.j<=ny+overlap_temp.and.k>=1-overlap_temp.and.k<=nz+overlap_temp)then
            icount = icount + 1
            fluid_boundary_nodes(icount)%ix = i
            fluid_boundary_nodes(icount)%iy = j
            fluid_boundary_nodes(icount)%iz = k
            fluid_boundary_nodes(icount)%nwx =  fluid_boundary_nodes_global(num)%nwx
            fluid_boundary_nodes(icount)%nwy =  fluid_boundary_nodes_global(num)%nwy
            fluid_boundary_nodes(icount)%nwz =  fluid_boundary_nodes_global(num)%nwz
            fluid_boundary_nodes(icount)%cos_theta =dcos(theta)         !assign contact angle
        endif
    enddo    

    deallocate(solid_boundary_nodes_global,fluid_boundary_nodes_global)    
    
    return
end subroutine geometry_preprocessing_new



!=====================================================================================================================================
!---------------------- geometry_preprocessing based on existing preprocessed data ----------------------
! load the geometry information from precomputed file. recommended for relatively large domain
! the geometry in the simulation must 100% match the preprocessed geometry data!!!
! Due to the issue with different padding schemes used in different compilers for MPI_Type_create_struct,
! The compiler used to compute the geometry info must be the same with the compiler for simulation
! I.e., Intel - Intel, PGI - PGI, GCC - GCC
!=====================================================================================================================================
subroutine geometry_preprocessing_load
    use Misc_module
    use Fluid_multiphase, only: theta   !contact angle 
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,L,M,N,k,icount,out1,out2,out3,num,overlap_temp
    type indirect_fluid_boundary_nodes_tmp                   !without contact angle information, used to load data
        integer :: ix,iy,iz
        real(kind=8) :: nwx,nwy,nwz 
    end type indirect_fluid_boundary_nodes_tmp
    type(indirect_solid_boundary_nodes), allocatable, dimension(:) :: solid_boundary_nodes_global
    type(indirect_fluid_boundary_nodes_tmp), allocatable, dimension(:) :: fluid_boundary_nodes_global
    type(indirect_solid_boundary_nodes) :: sdummy
    type(indirect_fluid_boundary_nodes_tmp) :: fdummy
    integer, dimension(:) ::blocklengths(3), types(3)
    integer(KIND=MPI_ADDRESS_KIND)  :: extent,lb ,base ,displacements(3)
    integer :: mpi_solid_bc_type,mpi_fluid_bc_type
    character (len=300) :: dummy   !file name
    
    blocklengths(1) = 4
    blocklengths(2) = 18
    blocklengths(3) = 1
    types(1) = MPI_INTEGER
    types(2) = MPI_INTEGER
    types(3) = MPI_DOUBLE_PRECISION

    call MPI_GET_ADDRESS(sdummy%ix, displacements(1), ierr)
    call MPI_GET_ADDRESS(sdummy%neighbor_list(1), displacements(2), ierr)
    call MPI_GET_ADDRESS(sdummy%la_weight, displacements(3), ierr)
    base = displacements(1)
    displacements(1) = displacements(1) - base
    displacements(2) = displacements(2) - base
    displacements(3) = displacements(3) - base
    call MPI_Type_create_struct(3,blocklengths,displacements,types,mpi_solid_bc_type,ierr)
    call MPI_TYPE_COMMIT(mpi_solid_bc_type,ierr)
    !call MPI_TYPE_get_EXTENT(mpi_solid_bc_type,lb,extent,ierr) 
    ! print*, extent,'s-type'

    blocklengths(1) = 3
    blocklengths(2) = 3
    types(1) = MPI_INTEGER
    types(2) = MPI_double_precision
    call MPI_GET_ADDRESS(fdummy%ix, displacements(1), ierr) 
    call MPI_GET_ADDRESS(fdummy%nwx, displacements(2), ierr) 
    base = displacements(1)
    displacements(1) = displacements(1) - base
    displacements(2) = displacements(2) - base
    call MPI_Type_create_struct(2,blocklengths,displacements,types,mpi_fluid_bc_type,ierr)
    call MPI_TYPE_COMMIT(mpi_fluid_bc_type,ierr)  
    ! call MPI_TYPE_get_EXTENT(mpi_fluid_bc_type,lb,extent,ierr)     
    ! print*, extent,'f-type'

    if(id==0)then  
        OPEN(UNIT=11,FILE=trim(adjustl(geo_boundary_file_path)),FORM='unformatted', status='old',access='stream')
        read(11)num_solid_boundary_global,num_fluid_boundary_global
        print*, "total number of solid boundary nodes= ", num_solid_boundary_global
        print*, "total number of fluid boundary nodes= ", num_fluid_boundary_global
    endif   
    call MPI_Bcast(num_solid_boundary_global,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)   
    call MPI_Bcast(num_fluid_boundary_global,1,MPI_INTEGER,0,MPI_COMM_VGRID,ierr)   
    allocate(solid_boundary_nodes_global(num_solid_boundary_global))  
    allocate(fluid_boundary_nodes_global(num_fluid_boundary_global))  
    
    if(id==0)then         
        read(11)solid_boundary_nodes_global
        read(11)fluid_boundary_nodes_global
        close(11)        
        ! print*,'Checking geometry data load status...'
        ! print*,fluid_boundary_nodes_global(111)
        ! print*,'End checking geometry data read correctness'
    endif  

    call MPI_Bcast(solid_boundary_nodes_global,num_solid_boundary_global,mpi_solid_bc_type,0,MPI_COMM_VGRID,ierr) 
    call MPI_Bcast(fluid_boundary_nodes_global,num_fluid_boundary_global,mpi_fluid_bc_type,0,MPI_COMM_VGRID,ierr) 

    CALL MPI_CART_COORDS(MPI_COMM_VGRID, id, mpi_dim, mpi_coords, ierr)
    out1=mpi_coords(1)
    out2=mpi_coords(2)
    out3=mpi_coords(3)

    icount=0
    overlap_temp = 3
    !$OMP parallel do private(l,m,n,i,j,k) reduction (+:icount)
    do num=1, num_solid_boundary_global
        l= solid_boundary_nodes_global(num)%ix
        m= solid_boundary_nodes_global(num)%iy
        n= solid_boundary_nodes_global(num)%iz
        i = l-out1*nx
        j = m-out2*ny
        k = n-out3*nz
        !three layers of ghost nodes, used in CSF calculation
        if(i>=1-overlap_temp.and.i<=nx+overlap_temp.and.j>=1-overlap_temp.and.j<=ny+overlap_temp.and.k>=1-overlap_temp.and.k<=nz+overlap_temp)then    
            icount = icount + 1
        endif
    enddo 
    num_solid_boundary = icount
    allocate(solid_boundary_nodes(num_solid_boundary))
    !print*, "total number of solid boundary nodes in local MPI domain",id,"= ", num_solid_boundary
    icount=0
    overlap_temp = 2
    !$OMP parallel do private(l,m,n,i,j,k) reduction (+:icount)
    do num=1, num_fluid_boundary_global
        l= fluid_boundary_nodes_global(num)%ix
        m= fluid_boundary_nodes_global(num)%iy
        n= fluid_boundary_nodes_global(num)%iz
        i = l-out1*nx
        j = m-out2*ny
        k = n-out3*nz
        if(i>=1-overlap_temp.and.i<=nx+overlap_temp.and.j>=1-overlap_temp.and.j<=ny+overlap_temp.and.k>=1-overlap_temp.and.k<=nz+overlap_temp)then
            icount = icount + 1
        endif
    enddo 
    num_fluid_boundary = icount
    allocate(fluid_boundary_nodes(num_fluid_boundary))
    !print*, "total number of fluid boundary nodes in local MPI domain",id,"= ", num_fluid_boundary

    icount=0
    overlap_temp = 3  
    do num=1, num_solid_boundary_global
        l= solid_boundary_nodes_global(num)%ix
        m= solid_boundary_nodes_global(num)%iy
        n= solid_boundary_nodes_global(num)%iz
        i = l-out1*nx
        j = m-out2*ny
        k = n-out3*nz
        !three layer of ghost nodes, used in CSF calculation
        if(i>=1-overlap_temp.and.i<=nx+overlap_temp.and.j>=1-overlap_temp.and.j<=ny+overlap_temp.and.k>=1-overlap_temp.and.k<=nz+overlap_temp)then   
            icount = icount + 1
            solid_boundary_nodes(icount)%ix = i
            solid_boundary_nodes(icount)%iy = j
            solid_boundary_nodes(icount)%iz = k
            solid_boundary_nodes(icount)%i_fluid_num = solid_boundary_nodes_global(num)%i_fluid_num
            solid_boundary_nodes(icount)%la_weight = solid_boundary_nodes_global(num)%la_weight
            solid_boundary_nodes(icount)%neighbor_list(:)=solid_boundary_nodes_global(num)%neighbor_list(:)
        endif
    enddo  

    icount=0 
    overlap_temp = 2 
    do num=1, num_fluid_boundary_global
        l= fluid_boundary_nodes_global(num)%ix
        m= fluid_boundary_nodes_global(num)%iy
        n= fluid_boundary_nodes_global(num)%iz
        i = l-out1*nx
        j = m-out2*ny
        k = n-out3*nz
        if(i>=1-overlap_temp.and.i<=nx+overlap_temp.and.j>=1-overlap_temp.and.j<=ny+overlap_temp.and.k>=1-overlap_temp.and.k<=nz+overlap_temp)then
            icount = icount + 1
            fluid_boundary_nodes(icount)%ix = i
            fluid_boundary_nodes(icount)%iy = j
            fluid_boundary_nodes(icount)%iz = k
            fluid_boundary_nodes(icount)%nwx =  fluid_boundary_nodes_global(num)%nwx
            fluid_boundary_nodes(icount)%nwy =  fluid_boundary_nodes_global(num)%nwy
            fluid_boundary_nodes(icount)%nwz =  fluid_boundary_nodes_global(num)%nwz
            fluid_boundary_nodes(icount)%cos_theta =dcos(theta)         !assign contact angle
        endif
    enddo  

    deallocate(solid_boundary_nodes_global,fluid_boundary_nodes_global)          
    return
end subroutine geometry_preprocessing_load


!**************************************** check preprocessed data ****************************************
! only used for debugging
subroutine check_geometry_linked_data
    use Misc_module
    use Fluid_multiphase, only: theta   !contact angle 
    use mpi_variable
    IMPLICIT NONE
    include 'mpif.h'
    integer :: i,j,k,m,num
    character (len=30) :: flnm   !file name

    write(flnm,"('solid',i4.4)")id
    open(unit=9+id, file='out3.field_data/'//trim(flnm), status='replace')
    do num=1, num_solid_boundary                        !cover MPI domain with three layer of ghost nodes, to be used in color gradient calculation    
        write(9+id,"(I8,1x,3(I4,1x),e14.7)")num,solid_boundary_nodes(num)%ix,solid_boundary_nodes(num)%iy,solid_boundary_nodes(num)%iz,solid_boundary_nodes(num)%la_weight  
    enddo 
    close(9+id)

    write(flnm,"('fluid',i4.4)")id
    open(unit=9+id, file='out3.field_data/'//trim(flnm), status='replace')
    do num=1, num_fluid_boundary                        !cover MPI domain with three layer of ghost nodes, to be used in color gradient calculation    
        write(9+id,"(I8,3(1x,I4),3(1x,e14.7))")num,fluid_boundary_nodes(num)%ix,fluid_boundary_nodes(num)%iy,fluid_boundary_nodes(num)%iz,&
        & fluid_boundary_nodes(num)%nwx,fluid_boundary_nodes(num)%nwy,fluid_boundary_nodes(num)%nwz  
    enddo 
    close(9+id)

    return
end subroutine check_geometry_linked_data
