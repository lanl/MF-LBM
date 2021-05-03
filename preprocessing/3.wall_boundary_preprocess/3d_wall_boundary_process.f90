module preprocessing
    implicit none
    integer :: nx,ny,nz,nx_sample,ny_sample,nz_sample, num_solid_boundary, num_fluid_boundary
    integer (kind=1), allocatable, dimension(:,:,:) :: walls
    double precision, allocatable, dimension(:,:,:) :: walls_smooth_1, walls_smooth_2
    integer  :: mirror_bc  !mirror boundary condition indicators
    integer  :: wall_data_status  ! 1: existing walldata; 0: creat new wall geometry
    integer  :: i_max_iteration, ghost_layers, overlap_phi
    integer  :: iper, jper, kper  !periodic boundary condition indicators
    !domain_wall_status at the edges of the computation domain: 1-wall; 0-nowall    
    integer :: domain_wall_status_x_min,domain_wall_status_x_max,domain_wall_status_y_min,domain_wall_status_y_max,domain_wall_status_z_min,domain_wall_status_z_max
    
    integer, dimension(:), parameter ::  ex(0:26)=(/0, 1, -1,  0,  0,  0,  0, 1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1 /)
    integer, dimension(:), parameter ::  ey(0:26)=(/0, 0,  0,  1, -1,  0,  0, 1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1, -1,  1, -1,  1 /)
    integer, dimension(:), parameter ::  ez(0:26)=(/0, 0,  0,  0,  0,  1, -1, 0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1, -1,  1,  1, -1, -1,  1,  1, -1 /)
    double precision, parameter, dimension(:) :: w(0:3) = (/8.0d0/27.0d0, 2.0d0/27.0d0, 1.0d0/54.0d0, 1.0d0/216.0d0/)
    double precision, dimension(:), parameter :: w_equ(0:18)=(/1d0/3d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,&
        & 1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)
 
    !gradient calculation weights
    double precision, parameter, dimension(:) :: ISO4(1:2) = (/1.0d0/6.0d0, 1.0d0/12.0d0/)
    double precision, parameter, dimension(:) :: ISO8(1:7) = (/4.0d0/45.0d0, 1.0d0/21.0d0, 2.0d0/105.0d0, 5.0d0/504.0d0, 1.0d0/315.0d0, 1.0d0/630.0d0, 1.0d0/5040.0d0/)
 
    type solid_boundary_nodes
        integer  :: ix,iy,iz,i_fluid_num         !number of fluid nodes connected to the boundary node
        integer, dimension (1:18) :: neighbor_list
        double precision :: la_weight       
    end type solid_boundary_nodes
    
    type fluid_boundary_nodes
        integer  :: ix,iy,iz
        double precision :: snx,sny,snz 
    end type fluid_boundary_nodes
   
    type(solid_boundary_nodes),allocatable,dimension(:) :: boundary_nodes_solid   
    type(fluid_boundary_nodes),allocatable,dimension(:) :: boundary_nodes_fluid   
    
end module preprocessing

!**********************main***************************
program Main
    use preprocessing
    implicit none
    integer :: i,j,k
   
    call initialization

    print*, 'initialization completed!'
   
    call smoothing

    print*, 'smoothing procedure completed!'
    
    call gradient

    print*, 'gradient procedure completed!'
        
    call visual

    print*, 'vtk file with solid wall orientation created!'
	
    call output

    print*, 'Wall boundary nodes info file created!'
	   
    deallocate(boundary_nodes_solid,boundary_nodes_fluid)
    deallocate(walls, walls_smooth_1, walls_smooth_2)
	
end program Main
!**********************main***************************

!*******************************************initialization**************************************************
subroutine initialization
    use preprocessing
    implicit none
    integer :: i,j,k,n,icount1,icount2,num,ncount
    character(len=300) :: wall_file_path, dummy
 
    call system('mkdir ./output')
  
    open(5,file='./input_parameters.txt',status='old')
    read(5,'(A)') dummy
    read(5,'(A)') wall_file_path
    read(5,*)iper,jper,kper
    read(5,*)domain_wall_status_x_min,domain_wall_status_x_max
    read(5,*)domain_wall_status_y_min,domain_wall_status_y_max
    read(5,*)domain_wall_status_z_min,domain_wall_status_z_max
    ! read(5,*)i_max_iteration  !smoothing iterations, default is 4
    ! read(5,*)ghost_layers,overlap_phi   
    close(5) 

    i_max_iteration = 4
    ghost_layers = 6
    overlap_phi = 4

    if(domain_wall_status_x_min==1)print*, 'Apply nonslip boundary condition at x = 1'
    if(domain_wall_status_x_max==1)print*, 'Apply nonslip boundary condition at x = xmax'   
    if(domain_wall_status_y_min==1)print*, 'Apply nonslip boundary condition at y = 1'
    if(domain_wall_status_y_max==1)print*, 'Apply nonslip boundary condition at y = ymax'
    if(domain_wall_status_z_min==1)print*, 'Apply nonslip boundary condition at z = 1'
    if(domain_wall_status_z_max==1)print*, 'Apply nonslip boundary condition at z = zmax' 
    print*,'number of ghost layers for smoothing process =', ghost_layers, 'number of ghost layers for overlaping phi =', overlap_phi    
    if(iper==1)print*,"apply periodic boundary condition along x direction"
    if(jper==1)print*,"apply periodic boundary condition along y direction"
    if(kper==1)print*,"apply periodic boundary condition along z direction"

  
    OPEN(UNIT=9,FILE=trim(adjustl(wall_file_path)),STATUS='OLD',FORM='UNFORMATTED',access='stream')
    read(9)nx_sample,ny_sample,nz_sample    
    print*,'nx_sample = ', nx_sample, ' ny_sample =', ny_sample,' nz_sample =', nz_sample
    nx=nx_sample
    ny=ny_sample
    nz=nz_sample        

    ghost_layers =  ghost_layers + overlap_phi  ! additional layers used for solid/fluid boundary nodes that outside the domain
    allocate(walls(1-ghost_layers:nx+ghost_layers,1-ghost_layers:ny+ghost_layers,1-ghost_layers:nz+ghost_layers))
    allocate(walls_smooth_1(1-ghost_layers:nx+ghost_layers,1-ghost_layers:ny+ghost_layers,1-ghost_layers:nz+ghost_layers),walls_smooth_2(1-ghost_layers:nx+ghost_layers,1-ghost_layers:ny+ghost_layers,1-ghost_layers:nz+ghost_layers))             
    !$OMP PARALLEL 
    !$OMP DO private(i,j) 
    do k=1-ghost_layers,nz+ghost_layers
        do j=1-ghost_layers,ny+ghost_layers
            do i=1-ghost_layers,nx+ghost_layers               
                walls_smooth_1(i,j,k)= 0d0
                walls_smooth_2(i,j,k)= 0d0
                walls(i,j,k)= 0
            enddo
        enddo
    enddo   
    !$OMP END PARALLEL 

    read(9)(((walls(i,j,k),i=1,nx_sample),j=1,ny_sample),k=1,nz_sample)  
    close(9)

    !domain boundary treatment  (order from x to z to account for the channel flow inlet/outlet)      
    if(domain_wall_status_x_min==0)walls(1,:,:)=0     
    if(domain_wall_status_x_min==1)walls(1,:,:)=1 
    if(domain_wall_status_x_max==0)walls(nx,:,:)=0     
    if(domain_wall_status_x_max==1)walls(nx,:,:)=1         
    
    if(domain_wall_status_y_min==0)walls(2:nx-1,1,:)=0     
    if(domain_wall_status_y_min==1)walls(2:nx-1,1,:)=1   
    if(domain_wall_status_y_max==0)walls(2:nx-1,ny,:)=0     
    if(domain_wall_status_y_max==1)walls(2:nx-1,ny,:)=1  
    
    if(domain_wall_status_z_min==0)walls(2:nx-1,2:ny-1,1)=0     
    if(domain_wall_status_z_min==1)walls(2:nx-1,2:ny-1,1)=1   
    if(domain_wall_status_z_max==0)walls(2:nx-1,2:ny-1,nz)=0     
    if(domain_wall_status_z_max==1)walls(2:nx-1,2:ny-1,nz)=1      

    if(kper==0)then  !z direction, not periodic       
        !$OMP PARALLEL DO private(i)
        do j=1,ny
            do i=1,nx                      
                walls(i,j,1-ghost_layers:0)=  walls(i,j,1)                  
                walls(i,j,nz+1:nz+ghost_layers)=  walls(i,j,nz)               
            enddo
        enddo 
    else            !z direction,periodic       
        !$OMP PARALLEL DO private(i)
        do j=1,ny
            do i=1,nx                      
                walls(i,j,1-ghost_layers:0) =  walls(i,j,nz-ghost_layers+1:nz)                  
                walls(i,j,nz+1:nz+ghost_layers)=  walls(i,j,1:ghost_layers)               
            enddo
        enddo 
    endif
        
    if(jper==0)then  !y direction, not periodic       
        !$OMP PARALLEL DO private(i)
        do k=1-ghost_layers,nz+ghost_layers
            do i=1,nx                       
                walls(i,1-ghost_layers:0,k)=  walls(i,1,k)                    
                walls(i,ny+1:ny+ghost_layers,k)=  walls(i,ny,k)               
            enddo
        enddo 
    else        !y direction,periodic  
        !$OMP PARALLEL DO private(i)   
        do k=1-ghost_layers,nz+ghost_layers
            do i=1,nx                       
                walls(i,1-ghost_layers:0,k)=  walls(i,ny-ghost_layers+1:ny,k)                    
                walls(i,ny+1:ny+ghost_layers,k)=  walls(i,1:ghost_layers,k)               
            enddo
        enddo 
    endif  
         
    if(iper==0)then  !x direction, not periodic  
        !$OMP PARALLEL DO private(j)
        do k=1-ghost_layers,nz+ghost_layers
            do j=1-ghost_layers,ny+ghost_layers                       
                walls(1-ghost_layers:0,j,k)=  walls(1,j,k)                    
                walls(nx+1:nx+ghost_layers,j,k)=  walls(nx,j,k)               
            enddo
        enddo
    else          !x direction, periodic  
         !$OMP PARALLEL DO private(j)
        do k=1-ghost_layers,nz+ghost_layers
            do j=1-ghost_layers,ny+ghost_layers                       
                walls(1-ghost_layers:0,j,k)=  walls(nx-ghost_layers+1:nx,j,k)                    
                walls(nx+1:nx+ghost_layers,j,k)=  walls(1:ghost_layers,j,k)               
            enddo
        enddo
    endif
 
    if(mirror_bc==1)then
        do k=1-ghost_layers,nz+ghost_layers
            do i=1-ghost_layers,nx+ghost_layers
                walls(i,1-ghost_layers:0,k)=  walls(i,2:ghost_layers+1,k)
            enddo
        enddo
    endif
     
    !$OMP parallel DO private(i,j) 
    do k=1-ghost_layers,nz+ghost_layers
        do j=1-ghost_layers,ny+ghost_layers
            do i=1-ghost_layers,nx+ghost_layers  
                walls_smooth_1(i,j,k) =  walls(i,j,k)        
                walls_smooth_2(i,j,k) =  walls(i,j,k)                   
            enddo
        enddo
    enddo        
   
     
      !extend wall information
    !$OMP parallel DO private(i,j,n) 
    do k=1-ghost_layers+1,nz+ghost_layers-1
        do j=1-ghost_layers+1,ny+ghost_layers-1
            do i=1-ghost_layers+1,nx+ghost_layers-1
                if(walls(i,j,k)==1)then
                    do n=1,18
                        if(walls(i+ex(n),j+ey(n),k+ez(n))<=0)then
                            walls(i,j,k)=2   !solid boundary nodes
                            exit
                        endif
                    enddo
                endif   
                if(walls(i,j,k)==0)then
                    do n=1,18
                        if(walls(i+ex(n),j+ey(n),k+ez(n))>=1)then
                            walls(i,j,k)=-1   !solid boundary nodes
                            exit
                        endif
                    enddo
                endif 
            enddo
        enddo
    enddo   
  
        
    icount1 = 0 
    icount2 = 0
     !$omp parallel do private(i,j) reduction (+:icount1,icount2)
    do k=1-overlap_phi,nz+overlap_phi             !smoothing process includs overlap_phi region used for cnx,cny,cnz calcaulation
        do j=1-overlap_phi,ny+overlap_phi
            do i=1-overlap_phi,nx+overlap_phi
                if(walls(i,j,k)==2)then
                    icount1 = icount1 + 1
                endif
                if(walls(i,j,k)==-1)then
                    icount2 = icount2 + 1
                endif
            enddo
        enddo
    enddo
    num_solid_boundary = icount1
    num_fluid_boundary = icount2
    print*, 'total number of solid boundary nodes = ', num_solid_boundary
    print*, 'total number of fluid boundary nodes = ', num_fluid_boundary
   
    
    allocate(boundary_nodes_solid(num_solid_boundary))  
    allocate(boundary_nodes_fluid(num_fluid_boundary))
    
    icount1 = 0
    icount2 = 0
    do k=1-overlap_phi,nz+overlap_phi             !smoothing process includs overlap_phi region used for cnx,cny,cnz calcaulation
        do j=1-overlap_phi,ny+overlap_phi
            do i=1-overlap_phi,nx+overlap_phi
                if(walls(i,j,k)==2)then
                    icount1 = icount1 + 1
                    boundary_nodes_solid(icount1)%ix = i
                    boundary_nodes_solid(icount1)%iy = j
                    boundary_nodes_solid(icount1)%iz = k
                    ncount = 0                   
                    boundary_nodes_solid(icount1)%la_weight = 0d0                  
                    
                    do n=1,18
                        if(walls(i+ex(n),j+ey(n),k+ez(n))<=0)then
                            ncount= ncount + 1
                            boundary_nodes_solid(icount1)%la_weight = boundary_nodes_solid(icount1)%la_weight + w_equ(n)      
                            boundary_nodes_solid(icount1)%neighbor_list(ncount)=n                    
                        endif
                    enddo
                    boundary_nodes_solid(icount1)%i_fluid_num = ncount  !number of fluid nodes connected to the boundary node                             
                endif
                if(walls(i,j,k)==-1)then
                    icount2 = icount2 + 1
                    boundary_nodes_fluid(icount2)%ix = i
                    boundary_nodes_fluid(icount2)%iy = j
                    boundary_nodes_fluid(icount2)%iz = k                
                endif
            enddo
        enddo
    enddo
    

 
    return
end subroutine initialization 
!*******************************************initialization**************************************************

!*******************************************Smoothing iterations**************************************************
subroutine smoothing
    use preprocessing
    integer :: i,j,k,m,n
  
    do iteration=1,i_max_iteration
        !$OMP PARALLEL 
        !$OMP do private(i,j,n,m)
        do k=1-ghost_layers+1,nz+ghost_layers-1
            do j=1-ghost_layers+1,ny+ghost_layers-1
                do i=1-ghost_layers+1,nx+ghost_layers-1
                    walls_smooth_2(i,j,k)= 0d0 
                    do n=0,26   
                        m = ex(n)*ex(n) + ey(n)*ey(n) + ez(n)*ez(n)      
                        walls_smooth_2(i,j,k)= walls_smooth_2(i,j,k) + walls_smooth_1(i+ex(n),j+ey(n),k+ez(n))*w(m)
                    enddo
                enddo
            enddo
        enddo
        
         !$OMP do private(i,j)
        do k=1-ghost_layers+1,nz+ghost_layers-1
            do j=1-ghost_layers+1,ny+ghost_layers-1
                do i=1-ghost_layers+1,nx+ghost_layers-1                    
                    walls_smooth_1(i,j,k)= walls_smooth_2(i,j,k) 
                enddo
            enddo
        enddo
       !$OMP end PARALLEL  
    enddo
    
    if(mirror_bc==1)then
        do k=1-ghost_layers,nz+ghost_layers
            do i=1-ghost_layers,nx+ghost_layers                
                walls_smooth_2(i,1-ghost_layers:0,k)=  walls_smooth_2(i,2:ghost_layers+1,k)
            enddo
        enddo 
    endif
      
    return
end subroutine smoothing
!*******************************************Smoothing iterations************************************************** 

!*******************************************gradient calculation**************************************************
subroutine gradient
    use preprocessing
    integer :: i,j,k,m,n
    double precision :: nwx,nwy,nwz,tmp
  
       !$omp parallel do private (i,j,k,nwx,nwy,nwz,tmp)
    do  num=1,num_fluid_boundary
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cx~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        i = boundary_nodes_fluid(num)%ix
        j = boundary_nodes_fluid(num)%iy
        k = boundary_nodes_fluid(num)%iz
                    
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cx~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

            ISO8(3)*(&
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

        tmp = 1d0/(dsqrt (nwx*nwx+nwy*nwy+nwz*nwz) + 1d-16)          !normalized vector
        boundary_nodes_fluid(num)%snx = nwx*tmp 
        boundary_nodes_fluid(num)%sny = nwy*tmp
        boundary_nodes_fluid(num)%snz = nwz*tmp   
        
    enddo
          
    return
end subroutine gradient
!*******************************************gradient calculation************************************************** 

!*******************************************visualization**************************************************


subroutine visual  
    use preprocessing
    IMPLICIT NONE
    integer :: i,j,k,num
    real*4,ALLOCATABLE,DIMENSION(:,:,:) :: nwx,nwy,nwz  
    character :: buffer*80, lf*1, str1*10, str2*10, str3*10
    integer   :: ivtk = 9    
 
    allocate(nwx(1-ghost_layers:nx+ghost_layers,1-ghost_layers:ny+ghost_layers,1-ghost_layers:nz+ghost_layers),&
    & nwy(1-ghost_layers:nx+ghost_layers,1-ghost_layers:ny+ghost_layers,1-ghost_layers:nz+ghost_layers),nwz(1-ghost_layers:nx+ghost_layers,1-ghost_layers:ny+ghost_layers,1-ghost_layers:nz+ghost_layers))
    
      !$omp parallel do
    do k=1-ghost_layers,nz+ghost_layers
        do j=1-ghost_layers,ny+ghost_layers
            do i=1-ghost_layers,nx+ghost_layers                 
                nwx(i,j,k)= 0d0
                nwy(i,j,k)= 0d0
                nwz(i,j,k)= 0d0
            enddo
        enddo
    enddo
    
     !$omp parallel do
    do  num=1,num_fluid_boundary
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cx~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        i = boundary_nodes_fluid(num)%ix
        j = boundary_nodes_fluid(num)%iy
        k = boundary_nodes_fluid(num)%iz
        nwx(i,j,k) = boundary_nodes_fluid(num)%snx
        nwy(i,j,k) = boundary_nodes_fluid(num)%sny
        nwz(i,j,k) = boundary_nodes_fluid(num)%snz
    enddo

    OPEN(UNIT = ivtk, FILE ='./output/visual.vtk', FORM='unformatted',access='stream',status='replace',convert='BIG_ENDIAN')

    lf = char(10) ! line feed character

    buffer = '# vtk DataFile Version 3.0'//lf
    write(ivtk) trim(buffer)

    buffer = 'vtk output'//lf
    write(ivtk) trim(buffer)

    buffer = 'BINARY'//lf
    write(ivtk) trim(buffer)

    buffer = 'DATASET STRUCTURED_POINTS '//lf
    write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')nx+2*ghost_layers
    write(str2(1:10),'(i10)')ny+2*ghost_layers
    write(str3(1:10),'(i10)')nz+2*ghost_layers
    buffer = 'DIMENSIONS '//str1//' '//str2//' '//str3//lf
    write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1-ghost_layers
    write(str2(1:10),'(i10)')1-ghost_layers
    write(str3(1:10),'(i10)')1-ghost_layers
    buffer = 'ORIGIN '//str1//' '//str2//' '//str3//lf
    write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')1
    write(str2(1:10),'(i10)')1
    write(str3(1:10),'(i10)')1
    buffer = 'SPACING '//str1//' '//str2//' '//str3//lf
    write(ivtk) trim(buffer)

    write(str1(1:10),'(i10)')(nx+2*ghost_layers)*(ny+2*ghost_layers)*(nz+2*ghost_layers)
    buffer = 'POINT_DATA '//str1//lf
    write(ivtk) trim(buffer)

    !scalar - walls
    buffer = 'SCALARS walls float'//lf
    write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write(ivtk) trim(buffer)
    write(ivtk)(((real(walls_smooth_2(i,j,k)),i=1-ghost_layers,nx+ghost_layers),j=1-ghost_layers,ny+ghost_layers),k=1-ghost_layers,nz+ghost_layers)

         
    buffer = 'SCALARS velocity_X float'//lf
    write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write(ivtk) trim(buffer)
    write(ivtk)(((nwx(i,j,k),i=1-ghost_layers,nx+ghost_layers),j=1-ghost_layers,ny+ghost_layers),k=1-ghost_layers,nz+ghost_layers)
       
    buffer = 'SCALARS velocity_Y float'//lf
    write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write(ivtk) trim(buffer)
    write(ivtk)(((nwy(i,j,k),i=1-ghost_layers,nx+ghost_layers),j=1-ghost_layers,ny+ghost_layers),k=1-ghost_layers,nz+ghost_layers)

    buffer = 'SCALARS velocity_Z float'//lf
    write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf
    write(ivtk) trim(buffer)
    write(ivtk)(((nwz(i,j,k),i=1-ghost_layers,nx+ghost_layers),j=1-ghost_layers,ny+ghost_layers),k=1-ghost_layers,nz+ghost_layers)

    close(ivtk)
        
    deallocate(nwx,nwy,nwz)

    return
end
!*******************************************visualization**************************************************
	
!*******************************************output**************************************************
subroutine output
    use preprocessing
    integer :: i,j,k,m,n,num
  
    open(11, file='./output/indirect_addressing_walls_boundary.dat', status='replace',form='unformatted',access='stream')
    write(11)num_solid_boundary, num_fluid_boundary
    !    do num=1,num_fluid_boundary
    !        write(11)boundary_nodes_fluid(num)%ix,boundary_nodes_fluid(num)%iy,boundary_nodes_fluid(num)%iz,boundary_nodes_fluid(num)%snx,boundary_nodes_fluid(num)%sny,boundary_nodes_fluid(num)%snz
    !    enddo
    !    do num=1,num_solid_boundary
    !        write(11)boundary_nodes_solid(num)%ix,boundary_nodes_solid(num)%iy,boundary_nodes_solid(num)%iz,boundary_nodes_solid(num)%i_fluid_num
    !        do n=1,18
    !            write(11)boundary_nodes_solid(num)%neighbor_list(n)
    !        enddo
    !        write(11)boundary_nodes_solid(num)%la_weight
    !    enddo
    write(11)boundary_nodes_solid
    write(11)boundary_nodes_fluid
    close(11)   
    
    return
end subroutine output
!*******************************************output************************************************** 
	
