program lbm_2d
  implicit none
  ! Parameters
  real(8), parameter :: PI = 4.0d0*atan(1.0d0)
  integer :: num_threads=16
  integer :: nx,ny
  integer :: ncy, L_front,L_back,L_top,L_bot	
  integer :: nfront,nback,ntop,nbot
  integer, parameter :: q = 9
  integer :: lattice_type   
  integer :: nprocsx, nprocsy
  integer :: iplot, max_iter, isave, irestart, statsbegin,statsend, iplotstats
  integer :: restart_step = 1
  real(8), parameter :: rho0 = 1.0
  integer, allocatable,dimension(:,:) :: isfluid, issolid,isbound
  real(8), allocatable,dimension(:,:) :: x,y,xvar,yvar
  integer :: nbct	!Total number of points on the cylinder boundary
  real*8, allocatable,dimension(:,:) :: boundary
  integer, allocatable,dimension(:) :: bou_i,bou_j,label_bc
  real(8) :: Re,uo,nu,omega, uinlet,vinlet
  real(8), allocatable,dimension(:) :: w,Hxx,Hyy,Hxy,Hxxy,Hxyy,Hxxyy
  integer, allocatable,dimension(:, :) :: e
  real(8), allocatable,dimension(:, :, :) :: f, f_eq, f_tmp,fin,fout
  real(8), allocatable, dimension(:, :) :: rho, p, ux, uy, mxx, myy, mxy,scalar
  real(8), allocatable, dimension(:) :: ut, un, mxxp, myyp, mxyp
  real(8) :: x_c,y_c,r_c,distance, sinth,costh,radius
  real(8), allocatable,dimension(:, :) :: er,er1
  real(8), allocatable,dimension(:) :: theta_c,p_mean_cy,radii,xb,yb,delta_uk,pb,p_fit
  real(8), allocatable,dimension(:) :: x_ref1,y_ref1,x_ref2,y_ref2,x_ref3,y_ref3,x_ref4,y_ref4
  real(8), allocatable,dimension(:) :: ux_ref1,uy_ref1,ux_ref2,uy_ref2,ux_ref3,uy_ref3,ux_ref4,uy_ref4
  real(8), allocatable,dimension(:) :: p_ref1,p_ref2,p_ref3,p_ref4
  integer, allocatable,dimension(:,:) :: i_p2,j_p2,i_p3,j_p3,i_p4,j_p4
  integer :: i, j, k, l, iter, coord_i,coord_j
  real(8) :: cu
  real(8) :: dx, dy, delx
  real(8) :: error1, frob1,max_radii
  real(8) :: L_phy, nu_phy, u_phy, delx_phy, delt_phy
  real(8) :: L_lat, nu_lat, u_lat, delx_lat, delt_lat
  character (len=6) :: output_dir_name = 'output'//char(0)
  character (len=7) :: restart_dir_name = 'restart'//char(0)
  character(len=100) :: filename,filename_bin
  character (len=100) :: bc_type
  logical :: x_periodic, y_periodic, channel_with_cylinder,incomp,channel_with_square,vel_interp
  logical, parameter :: verbose = .true.
  
  	!$ call omp_set_num_threads(num_threads)
  
  	call create_output_dir_if_needed(output_dir_name,6)
  	call create_output_dir_if_needed(restart_dir_name,7)
  	
  	call read_input_file()
  	
  	 !Get domain coordination
  	if (channel_with_cylinder) then
  		nfront = L_front*ncy
   		nback = L_back*ncy
   		ntop = L_top*ncy
   		nbot = L_bot*ncy
   		nx = nfront + ncy + nback
   		ny = ntop + ncy + nbot
   		
   		call allocate_memory()
   		call get_coord()
   	else
   		call allocate_memory()
  	end if
  	call write_grid_plot3d()
   	
   	dx = 1.0d0
   	dy = 1.0d0
   	
	!Flow constants: 
	uinlet  = uo		!Inlet axial velocity
	vinlet  = 0.0d0   	!inlet radial velocity 
	
	! Lattice weights and velocity sets
	w(0:q-1) = (/ 16.0/36.0, 4.0/36.0, 4.0/36.0, 4.0/36.0, 4.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 /)
	e(0:q-1, 1) = (/ 0, 1, 0, -1, 0, 1, -1, -1, 1 /)
	e(0:q-1, 2) = (/ 0, 0, 1, 0, -1, 1, 1, -1, -1 /)
	
	
	nu    = uinlet *(2.0d0*max_radii)/Re;     				!kinematic viscosity 
	omega = 1.0d0/((3.0*nu)+(1.0d0/2.0d0)) 	!relaxation parameter
	
	L_lat = 2.0d0*max_radii
	nu_lat = nu
	u_lat = uo
	delx_lat = 1.0d0
	delt_lat = 1.0d0
	
	L_phy = 1.0d0
	nu_phy = 0.01d0
	u_phy = Re*nu_phy/L_phy
	delx_phy = L_phy/L_lat
	
	delt_phy = (delx_phy**2)*nu_lat/nu_phy
	
	call create_master_p3d_file()
	
	open(unit=10,file='simulation_parameters.dat')
		write(10,*) 'nx,ny', nx,ny
		write(10,*) 'r_cylin', ncy/2
		write(10,*) 'centre_cylin', nfront + (ncy/2),nbot + ncy/2 
		write(10,*) 'Maximum_radius', max_radii
		write(10,*) 'Re:', Re
		write(10,*) 'uo:', uo
		write(10,*) 'nu:', nu
		write(10,*) 'tau:', 1.0d0/omega
		write(10,*) 'dt_phy:', delt_phy
	close(10)
  
	! Second-order Hermite polynomial 
	do k = 0, q-1
		Hxx(k) = e(k, 1)*e(k, 1) - (1.0d0/3.0d0)	!Hxx
		Hyy(k) = e(k, 2)*e(k, 2) - (1.0d0/3.0d0)	!Hyy
		Hxy(k) = e(k, 1)*e(k, 2)					!Hxy
		Hxxy(k) = (e(k, 1)*e(k, 1) - (1.0d0/3.0d0))*e(k, 2)									!Hxxy
		Hxyy(k) = (e(k, 2)*e(k, 2) - (1.0d0/3.0d0))*e(k, 1)									!Hxxy
		Hxxyy(k) = (e(k, 1)*e(k, 1) - (1.0d0/3.0d0))*(e(k, 2)*e(k, 2) - (1.0d0/3.0d0))		!Hxxy
	end do

	 if (irestart == 0) then
	 	! Initialize
	  	rho = rho0
	  	ux = 0.0d0
	  	uy = 0.0d0
	  	call initialize(f, rho, ux, uy)
  		write (filename_bin, 100) 'grid',10000000,'.bin',char(0)
  		call write_function_plot3d(filename_bin)
	 
  	else if (irestart == 1) then 
     	call read_restart_file(iter)
     	restart_step = iter
     	write(*, *) 'Restarting step = ', restart_step
  	end if
	
	
  ! Main loop
  do iter = restart_step+1, max_iter
	
	!Macroscopic variables:
	!$omp parallel do private(i,j) shared( f, e, rho, ux, uy) schedule(guided) 
	do i = 1, nx
		do j = 1, ny
		  rho(i, j) = sum(f(i, j, :))
		  ux(i, j) = sum(f(i, j, :) * e(:, 1)) / rho(i, j)
		  uy(i, j) = sum(f(i, j, :) * e(:, 2)) / rho(i, j)
		end do
    end do
    !$omp end parallel do
	
		
	!Boundary conditions:  MACROSCOPIC (DIRICHLET)
	
	!Right wall: 
	if(.not. x_periodic) then
		rho(nx,1:ny) = rho(nx-1,1:ny)
		ux(nx,1:ny) = ux(nx-1,1:ny)
		uy(nx,1:ny) = uy(nx-1,1:ny)
	end if
	
	!Left wall: INLET
	if(.not. x_periodic) then
		ux(1,1:ny) = uinlet
		uy(1,1:ny) = vinlet
	end if
	
	!Top wall: Wall
	if(.not. y_periodic) then
		ux(1:nx,ny) = 0.0d0
		uy(1:nx,ny) = 0.0d0
	end if
	
	!Bottom wall: wall
	if(.not. y_periodic) then
		ux(1:nx,1) = 0.0d0
		uy(1:nx,1) = 0.0d0
	end if
		
	!Regularized boundary condition: Moments calculations
	!=============================================================================================================================
	
	!-----------------------------------------------EDGES-----------------------------------------------------------------------
	!CYLINDER wall (No-Slip):---------------------------------------------------------------------------------------------------
	
	if(channel_with_cylinder) then
		!velocity interpolation for curved boundary: 
		if(vel_interp) call velocity_interpolation_boundary_nodes()
		
		do l = 1,nbct
			i = bou_i(l)
			j = bou_j(l)
			
			!Original coordinate system:
			call boundary_cases(bou_i(l),bou_j(l),label_bc(l),ux(i,j),uy(i,j))

		end do
	end if
	
	if(.not. x_periodic) then
		do j = 1,ny
			bc_type = "outlet_bc"
			coord_i = nx
			coord_j = j
			if (incomp) call apply_bc_incomp(coord_i,coord_j,bc_type)
			if (.not. incomp) call apply_bc(coord_i,coord_j,bc_type)
		end do
	end if
	
	if(.not. x_periodic) then
		do j = 1,ny
			bc_type = "inlet_bc"
			coord_i = 1
			coord_j = j
			if (incomp) call apply_bc_incomp(coord_i,coord_j,bc_type)
			if (.not. incomp) call apply_bc(coord_i,coord_j,bc_type)
		end do
	end if
	
	if(.not. y_periodic) then
		do i = 1,nx
			bc_type = "topwall_bc"
			coord_i = i
			coord_j = ny
			if (incomp) call apply_bc_incomp(coord_i,coord_j,bc_type)
			if (.not. incomp) call apply_bc(coord_i,coord_j,bc_type)
		end do
	end if
	
	if(.not. y_periodic) then
		do i = 1,nx
			bc_type = "bottomwall_bc"
			coord_i = i
			coord_j = 1
			if (incomp) call apply_bc_incomp(coord_i,coord_j,bc_type)
			if (.not. incomp) call apply_bc(coord_i,coord_j,bc_type)
		end do
	end if
	
		

	!=============================================================================================================================
	
	! do i = 1, nx
      ! do j = 1, ny
        ! rho(i, j) = sum(f(i, j, :))
      ! end do
    ! end do
	
    ! Collision step
	!----------------------------------------------------------------------------------------------------
		!$omp parallel do collapse(2) private(cu, k) shared(e, ux, uy, rho, w, f_eq, fout, f, omega)
    	do i = 1, nx
		  do j = 1, ny
		    do k = 0, q-1
			  cu = e(k, 1) * ux(i, j) + e(k, 2) * uy(i, j)
		      f_eq(i, j, k) = w(k)*rho(i, j)*( 1.0d0 + (3.0d0*cu) + (4.50d0*(cu*cu)) - &
									& 1.50d0*( ux(i, j)*ux(i, j) + uy(i, j)*uy(i, j)) )
		      fout(i, j, k) = f(i, j, k) - omega * (f(i, j, k) - f_eq(i, j, k))
		    end do
		  end do
    	end do
    	!$omp end parallel do
	
	! *** Modifying non-equilibrium moments with post-collisional populations***
	!$omp parallel do collapse(2) shared(rho, mxx,myy,mxy, Hxx,Hyy,Hxy, fout)
	do i = 1, nx
		do j = 1, ny 
			mxx(i, j) =  sum(fout(i, j, :) * Hxx(:)) / rho(i, j)
			myy(i, j) =  sum(fout(i, j, :) * Hyy(:)) / rho(i, j)
			mxy(i, j) =  sum(fout(i, j, :) * Hxy(:)) / rho(i, j)
		end do
    end do
    !$omp end parallel do
	
	
	!Kinetic Projection / Regularization (using modified moments)
	!----------------------------------------------------------------------------------------------------	
		!$omp parallel do collapse(2) private(cu, k) shared(e, ux, uy, rho, w, mxx,myy,mxy,Hxx,Hyy,Hxy, fout)
		do i = 1, nx
		  do j = 1, ny
		    do k = 0, q-1
			  cu = e(k, 1) * ux(i, j) + e(k, 2) * uy(i, j)
		      fout(i, j, k) = w(k)*rho(i, j)*( 1.0d0 + (3.0d0*cu) +  &
								& 4.50d0*mxx(i, j)*Hxx(k) + &
								& 4.50d0*myy(i, j)*Hyy(k) + &
								& 9.00d0*mxy(i, j)*Hxy(k) )
		    end do
		  end do
    	end do
		!$omp end parallel do
	
	
    ! Streaming step
    f_tmp = fout
    do i = 1, nx
      do j = 1, ny
        do k = 0, q-1
          if (i + e(k, 1) > 0 .and. i + e(k, 1) <= nx .and. j + e(k, 2) > 0 .and. j + e(k, 2) <= ny) then
            f(i + e(k, 1), j + e(k, 2), k) = f_tmp(i, j, k)
          end if
        end do
      end do
    end do
	
	if(x_periodic) then
		  do j = 1, ny
			do k = 0, q-1
			  if (1 - e(k, 1) < 1 ) f(1, j, k) = f_tmp(nx, j, k)
			  if (nx - e(k, 1) > nx) f(nx, j, k) = f_tmp(1, j, k)
			end do
		  end do
	end if
	
	if(y_periodic) then
		do i = 1, nx
			do k = 0, q-1
			  if (1 - e(k, 2) < 1 ) f(i, 1, k) = f_tmp(i, ny, k)
			  if (ny - e(k, 2) > ny) f(i, ny, k) = f_tmp(i, 1, k)
			end do
		end do
	end if
	if((iter .ge. statsbegin) .and. (iter .le. statsend)) then
		!Calculating mean quantities
		call calculate_mean(iter)
	end if
		
		frob1=0.0d0
		do i=1,nx
			do j=1,ny
				er(i,j)=abs(ux(i,j)-er1(i,j))
				frob1=frob1+(er(i,j)**2)
			end do
		end do
		error1=sqrt(frob1)
		er1(:,:)=ux(:,:)
	
	print*,iter,error1
	
	if(mod(iter,iplot)==0)then
		write (filename, 100) 'grid',10000000+iter,'.vtk',char(0)
		write (filename_bin, 100) 'grid',10000000+iter,'.bin',char(0)
		100 format (a4,I8,a4,a1)
		if(channel_with_cylinder .or. channel_with_square) then
			do  i = 1, nx
				do  j = 1, ny
					if(issolid(i,j)==1) then
						ux(i,j) = 0.0d0
						uy(i,j) = 0.0d0
					end if
				end do
			end do
		end if
		call write_function_plot3d(filename_bin)
	end if
	
	if (mod(iter, isave) == 0) then 
          call write_restart_file(iter)
    end if
    
    if (iter==40*iplot) then
    	open(unit=105,file='data_center_line.dat')
    		do i = 1, nx
    			write(105,*) x(i,ny/2), rho(i,ny/2),ux(i,ny/2),uy(i,ny/2)
    		end do
    	close(105)
    end if
       
  end do

  
  write (filename, 100) 'grid',10000000+iter,'.vtk',char(0)
  write (filename_bin, 100) 'grid',10000000+iter,'.bin',char(0)
	if(channel_with_cylinder .or. channel_with_square) then
		do  i = 1, nx
			do  j = 1, ny
				if(issolid(i,j)==1) then
					ux(i,j) = 0.0d0
					uy(i,j) = 0.0d0
				end if
			end do
		end do
	end if
  ! call output_vtk(filename)
  call write_function_plot3d(filename_bin)
  

contains

	subroutine calculate_mean(step)
		implicit none
		integer,intent(in) :: step 
		real(8) :: mean_counter
		real(8),dimension(1:nbct) :: pvar,thetavar
		real(8) :: temp_theta,temp_p
		real(8) :: x1,y1,p1,x2,y2,p2,x3,y3,p3,x4,y4,p4
		real(8) :: xd,yd,p_0,p_1
		real(8) :: thet,xs,ys,ps
		integer :: i,j,k
	
		p(:,:) = rho(:,:)/3.0d0
		
		mean_counter = real(step - statsbegin)
		do k = 1,nbct
			i = bou_i(k)
			j = bou_j(k)
			
			!Second reference node: First fluid point
			xd = (x_ref2(k) - x(i_p2(1,k),j_p2(1,k)))/dx
			yd = (y_ref2(k) - y(i_p2(1,k),j_p2(1,k)))/dy
			
			p_0 = (1.0d0 - xd)*p(i_p2(1,k),j_p2(1,k)) + xd*p(i_p2(2,k),j_p2(2,k))
			p_1 = (1.0d0 - xd)*p(i_p2(3,k),j_p2(3,k)) + xd*p(i_p2(4,k),j_p2(4,k))
			
			p_ref2(k) =  (1.0d0 - yd)*p_0 + yd*p_1
			
			!Third reference node: Second fluid point
			xd = (x_ref3(k) - x(i_p3(1,k),j_p3(1,k)))/dx
			yd = (y_ref3(k) - y(i_p3(1,k),j_p3(1,k)))/dy
			
			p_0 = (1.0d0 - xd)*p(i_p3(1,k),j_p3(1,k)) + xd*p(i_p3(2,k),j_p3(2,k))
			p_1 = (1.0d0 - xd)*p(i_p3(3,k),j_p3(3,k)) + xd*p(i_p3(4,k),j_p3(4,k))
			
			p_ref3(k) =  (1.0d0 - yd)*p_0 + yd*p_1
			
			!Fourth reference node: Third fluid point
			xd = (x_ref4(k) - x(i_p4(1,k),j_p4(1,k)))/dx
			yd = (y_ref4(k) - y(i_p4(1,k),j_p4(1,k)))/dy
			
			p_0 = (1.0d0 - xd)*p(i_p4(1,k),j_p4(1,k)) + xd*p(i_p4(2,k),j_p4(2,k))
			p_1 = (1.0d0 - xd)*p(i_p4(3,k),j_p4(3,k)) + xd*p(i_p4(4,k),j_p4(4,k))
			
			p_ref4(k) =  (1.0d0 - yd)*p_0 + yd*p_1
			
			xs = x_ref1(k) - x_c
			ys = y_ref1(k) - y_c
			x2 = x_ref2(k) - x_c
			y2 = y_ref2(k) - y_c
			x3 = x_ref3(k) - x_c
			y3 = y_ref3(k) - y_c
			x4 = x_ref4(k) - x_c
			y4 = y_ref4(k) - y_c
			
			p2 = p_ref2(k)
			p3 = p_ref3(k)
			p4 = p_ref4(k)
			call quadratic_interpolation(x2,y2,p2,x3,y3,p3,x4,y4,p4,xs,ys,ps)
			p_ref1(k) = ps			
			
			p_mean_cy(k) = ((mean_counter*p_mean_cy(k)) + p_ref1(k) )/(mean_counter + 1.0d0)	
  			!-----------------------------------------------------------------------------
  			
		end do
		
    	if(mod(step,iplotstats)==0) then
    		thetavar = theta_c
			pvar = p_mean_cy
			do i = 1, nbct - 1
		        do j = 1, nbct - i
		            if (thetavar(j) > thetavar(j + 1)) then
		            
		                ! Swap theta values
		                temp_theta = thetavar(j)
		                thetavar(j) = thetavar(j + 1)
		                thetavar(j + 1) = temp_theta
		                
		                temp_p = pvar(j)
		                pvar(j) = pvar(j + 1)
		                pvar(j + 1) = temp_p
		            end if
		        end do
		    end do
		 
			call gaussian_smoothing(nbct,thetavar,pvar,p_fit)

			open(unit=101,file="pmean.dat")
				do k = 1, nbct
					if(thetavar(k) .le. PI) then
						write(101,*) abs((thetavar(k)*180.0d0/PI)-180.0d0), pvar(k), p_fit(k), &
								& 2.0d0*(pvar(k) - (1.196512d0/3.0d0))/(1.1965120d0*0.10d0*0.10d0), &
								& 2.0d0*(p_fit(k) - (1.196512d0/3.0d0))/(1.1965120d0*0.10d0*0.10d0)
					end if
				end do
			close(101)
		end if

	end subroutine calculate_mean
	
	subroutine gaussian_smoothing(n_l,theta_l, p_l, p_smooth_l)
    implicit none
    integer, intent(in) :: n_l
    real(8), intent(in) :: theta_l(n_l), p_l(n_l)
    real(8), intent(out) :: p_smooth_l(n_l)
    real(8) :: weight_l, sum_weights_l, sum_values_l
    integer :: i_l, j_l
    real(8) :: sigma_l = 0.1d0

    do i_l = 1, n_l
       sum_weights_l = 0.0d0
       sum_values_l = 0.0d0
       do j_l = 1, n_l
          weight_l = exp(-((theta_l(i_l) - theta_l(j_l))**2) / (2.0d0 * sigma_l**2))
          sum_weights_l = sum_weights_l + weight_l
          sum_values_l = sum_values_l + weight_l * p_l(j_l)
       end do
       p_smooth_l(i_l) = sum_values_l / sum_weights_l
    end do
  end subroutine gaussian_smoothing
	
	subroutine quadratic_interpolation(x1,y1,p1,x2,y2,p2,x3,y3,p3,x_s,y_s,p_s)
		implicit none
		real(8),intent(in) :: x1,y1,p1,x2,y2,p2,x3,y3,p3,x_s,y_s
		real(8),intent(out) :: p_s
		real(8) :: a0,a1,a2,denom
		real(8) :: r1,r2,r3,rs
		
		r1 = sqrt((x1**2) + (y1**2))
		r2 = sqrt((x2**2) + (y2**2))
		r3 = sqrt((x3**2) + (y3**2))
		rs = sqrt((x_s**2) + (y_s**2))	
		
		! Step 1: Compute the coefficients a0, a1, a2 for the quadratic interpolation
		denom = (r1 - r2) * (r1 - r3) * (r2 - r3)
		
		! Calculate coefficients a0, a1, a2
		a0 = (r1*r3*p2*(r3-r1) + (r2**2)*(r3*p1 - r1*p3) + r2*((r1**2)*p3 - (r3**2)*p1))/denom
		a1 = ((r3**2)*(p1-p2) + (r1**2)*(p2-p3) + (r2**2)*(p3-p1) )/denom 
		a2 = (r3*(p2-p1) + r2*(p1-p3) + r1*(p3-p2))/denom
		
		!Step 2: Use the quadratic formula to compute the interpolated pressure at surface point (xs)
		p_s = a0 + a1*rs + a2*(rs**2)
		
	end subroutine quadratic_interpolation
	
	subroutine velocity_interpolation_boundary_nodes()
		implicit none
		real(8) :: term1,term2,term3
		real(8) :: xd,yd,ux_0,uy_0,ux_1,uy_1
		integer :: i,j,k
		
		do k = 1,nbct
			i = bou_i(k)
			j = bou_j(k)
			
			!First reference node: cylinder point
			ux_ref1(k) = 0.0d0
			uy_ref1(k) = 0.0d0
				
			!Second reference node: First fluid point
			xd = (x_ref2(k) - x(i_p2(1,k),j_p2(1,k)))/dx
			yd = (y_ref2(k) - y(i_p2(1,k),j_p2(1,k)))/dy
			
			ux_0 = (1.0d0 - xd)*ux(i_p2(1,k),j_p2(1,k)) + xd*ux(i_p2(2,k),j_p2(2,k))
			uy_0 = (1.0d0 - xd)*uy(i_p2(1,k),j_p2(1,k)) + xd*uy(i_p2(2,k),j_p2(2,k))
			
			ux_1 = (1.0d0 - xd)*ux(i_p2(3,k),j_p2(3,k)) + xd*ux(i_p2(4,k),j_p2(4,k))
			uy_1 = (1.0d0 - xd)*uy(i_p2(3,k),j_p2(3,k)) + xd*uy(i_p2(4,k),j_p2(4,k))
			
			ux_ref2(k) =  (1.0d0 - yd)*ux_0 + yd*ux_1
			uy_ref2(k) =  (1.0d0 - yd)*uy_0 + yd*uy_1
			
			!Third reference node: Second fluid point
			xd = (x_ref3(k) - x(i_p3(1,k),j_p3(1,k)))/dx
			yd = (y_ref3(k) - y(i_p3(1,k),j_p3(1,k)))/dy
			
			ux_0 = (1.0d0 - xd)*ux(i_p3(1,k),j_p3(1,k)) + xd*ux(i_p3(2,k),j_p3(2,k))
			uy_0 = (1.0d0 - xd)*uy(i_p3(1,k),j_p3(1,k)) + xd*uy(i_p3(2,k),j_p3(2,k))
			
			ux_1 = (1.0d0 - xd)*ux(i_p3(3,k),j_p3(3,k)) + xd*ux(i_p3(4,k),j_p3(4,k))
			uy_1 = (1.0d0 - xd)*uy(i_p3(3,k),j_p3(3,k)) + xd*uy(i_p3(4,k),j_p3(4,k))
			
			ux_ref3(k) =  (1.0d0 - yd)*ux_0 + yd*ux_1
			uy_ref3(k) =  (1.0d0 - yd)*uy_0 + yd*uy_1
			
			!Variable interpolation:
			term1 = (2.0d0*delx*delx - (delta_uk(k)**2) + 3.0d0* delta_uk(k)*delx)/(2.0d0*(delx**2))
			term2 =  delta_uk(k)*( delta_uk(k) - 2.0d0*delx)/(delx**2)
			term3 =  delta_uk(k)*( delta_uk(k) - delx)/(2.0d0*(delx**2))
			!Boundary node velocities:
			ux(i,j) = term1*ux_ref1(k) + term2*ux_ref2(k) - term3*ux_ref3(k) 
			uy(i,j) = term1*uy_ref1(k) + term2*uy_ref2(k) - term3*uy_ref3(k) 
			
	
		end do
		
	end subroutine velocity_interpolation_boundary_nodes
		

  subroutine get_coord()
    implicit none
	integer, dimension(0:q-1, 2) :: neigh
	integer :: flag0,flag1,flag2,flag3
	real, parameter :: epsilon = 0.5          ! Tolerance for boundary
    integer :: i, j, k,l, cnt, i1,j1,i2,j2,i3,j3,i_temp,j_temp,ss,ee
    real(8) :: thet,x_s,y_s,thet1,thet2,thet_mid,lambda
    real(8) :: xmid,ymid,epsi,delta_x,delta_y,rr,unit_nx,unit_ny

	epsi = 0.2
    isfluid = 0
	isbound = 0
	issolid = 0
	
	!Cylinder center:
	x_c = real(nfront + (ncy/2) + 1)
	y_c = real(nbot + ncy/2 + 1)
	r_c = real(ncy/2)
	
	do i = 1,nx
		x(i,:) = real(i)
		
	end do
	do j = 1,ny
		y(:,j) = real(j)
	end do
	
	!Identifying the boundary points and solid points
	do  i = 1, nx
		do  j = 1, ny
            distance = sqrt((x(i,j) - x_c)**2 + (y(i,j) - y_c)**2)
            if (abs(distance - r_c) <= epsilon) then
				isbound(i,j) = 1
				
			elseif (abs(distance) <= r_c) then
				issolid(i,j) = 1
				
			else
				isfluid(i,j) = 1
            end if
						
		end do
	end do
	
	do  i = 2, nx-1
		do  j = 2, ny-1
			if(isfluid(i,j)==1) then
				neigh(0,1:2) = (/ i,j /)
				neigh(1,1:2) = (/ i-1,j /)
				neigh(2,1:2) = (/ i,j-1 /)
				neigh(3,1:2) = (/ i+1,j /)
				neigh(4,1:2) = (/ i,j+1 /)
				neigh(5,1:2) = (/ i-1,j-1 /)
				neigh(6,1:2) = (/ i+1,j-1 /)
				neigh(7,1:2) = (/ i+1,j+1 /)
				neigh(8,1:2) = (/ i-1,j+1 /)
				do l = 0,q-1
					if((isfluid(neigh(l,1),neigh(l,2))==0) .and. (isbound(neigh(l,1),neigh(l,2))==0)) then
						isbound(neigh(l,1),neigh(l,2))=1
						issolid(neigh(l,1),neigh(l,2))=0
					end if
				end do
			end if
		end do
	end do
	
	nbct = sum(isbound(:,:))
	call allocate_cylinder_memory()
	

	cnt = 1
	do  i = 1, nx
		do  j = 1, ny
			if(isbound(i,j)==1) then
				bou_i(cnt) = i
				bou_j(cnt) = j
				
				flag0 = 1
				flag1 = 1
				flag2 = 1
				flag3 = 1
				if((isfluid(i-1,j)==0) .and. (isfluid(i-1,j+1)==0) .and. (isfluid(i,j+1)==0))then
					flag0 = 0
				end if
				if((isfluid(i+1,j)==0) .and. (isfluid(i+1,j+1)==0) .and. (isfluid(i,j+1)==0))then
					flag1 = 0
				end if
				if((isfluid(i-1,j)==0) .and. (isfluid(i-1,j-1)==0) .and. (isfluid(i,j-1)==0))then
					flag2 = 0
				end if
				if((isfluid(i+1,j)==0) .and. (isfluid(i+1,j-1)==0) .and. (isfluid(i,j-1)==0))then
					flag3 = 0
				end if
				
				label_bc(cnt) = flag0*1 + flag1*2 + flag2*4 + flag3*8
				
				cnt = cnt+1
			end if
		end do
	end do
	
	do k = 1, nbct
		i = bou_i(k)
		j = bou_j(k)
		radii(k) = sqrt(((x(i,j)-x_c)**2) + ((y(i,j)-y_c)**2) )			
	end do
	max_radii = radii(1)
	do k = 2,nbct
		if(radii(k) .gt. max_radii)then
			max_radii = radii(k)
		end if
	end do

	
	!Finding boundary and two fluid nodes for velocity interpolation
	!-----------------------------------------------------------------------------------------------
	delx = sqrt(2.0d0)
	do k = 1, nbct
			
			rr = sqrt((x(bou_i(k),bou_j(k)) - x_c)**2 + (y(bou_i(k),bou_j(k)) - y_c)**2)

			  if (rr > 1.0e-12) then
				unit_nx = (x(bou_i(k),bou_j(k)) - x_c)/rr
				unit_ny = (y(bou_i(k),bou_j(k)) - y_c)/rr
				xb(k) = x_c + (max_radii*unit_nx)
				yb(k) = y_c + (max_radii*unit_ny)
			  else
				 ! Point is exactly at center
				 unit_nx = 0.0
				 unit_ny = 0.0
				 xb(k) = x_c
				 yb(k) = y_c
			  end if
			  
			  x_ref1(k) = xb(k)
			  y_ref1(k) = yb(k)
			  
			  x_ref2(k) = xb(k) + (delx*unit_nx)
			  y_ref2(k) = yb(k) + (delx*unit_ny)
			  
			  x_ref3(k) = xb(k) + (2.0d0*delx*unit_nx)
			  y_ref3(k) = yb(k) + (2.0d0*delx*unit_ny)
			  
			  x_ref4(k) = xb(k) + (3.0d0*delx*unit_nx)
			  y_ref4(k) = yb(k) + (3.0d0*delx*unit_ny)
			  
			  !delta value for velocity interpolation
			  delta_x = x(bou_i(k),bou_j(k)) - xb(k)
			  delta_y = y(bou_i(k),bou_j(k)) - yb(k)
			  delta_uk(k) = sqrt(delta_x**2 + delta_y**2)
			  
			  !Finding four neighbour points for bi-linear interpolation
			  !   3 --- 4
			  !   |     |
			  !   1 --- 2
			  
			  !Neighbouring fluid nodes for reference point 2
			  i_p2(1,k) = floor(x_ref2(k))
			  i_p2(2,k) = floor(x_ref2(k)) + 1
			  i_p2(3,k) = floor(x_ref2(k))
			  i_p2(4,k) = floor(x_ref2(k)) + 1
			  
			  j_p2(1,k) = floor(y_ref2(k))
			  j_p2(2,k) = floor(y_ref2(k))
			  j_p2(3,k) = floor(y_ref2(k)) + 1
			  j_p2(4,k) = floor(y_ref2(k)) + 1
			  
			  !Neighbouring fluid nodes for reference point 3
			  i_p3(1,k) = floor(x_ref3(k))
			  i_p3(2,k) = floor(x_ref3(k)) + 1
			  i_p3(3,k) = floor(x_ref3(k))
			  i_p3(4,k) = floor(x_ref3(k)) + 1
			  
			  j_p3(1,k) = floor(y_ref3(k))
			  j_p3(2,k) = floor(y_ref3(k))
			  j_p3(3,k) = floor(y_ref3(k)) + 1
			  j_p3(4,k) = floor(y_ref3(k)) + 1
			  
			  !Neighbouring fluid nodes for reference point 4
			  i_p4(1,k) = floor(x_ref4(k))
			  i_p4(2,k) = floor(x_ref4(k)) + 1
			  i_p4(3,k) = floor(x_ref4(k))
			  i_p4(4,k) = floor(x_ref4(k)) + 1
			  
			  j_p4(1,k) = floor(y_ref4(k))
			  j_p4(2,k) = floor(y_ref4(k))
			  j_p4(3,k) = floor(y_ref4(k)) + 1
			  j_p4(4,k) = floor(y_ref4(k)) + 1
			  
			  
  			!theta angles for the boundary (cylinder) points 
  			theta_c(k) = atan(yb(k) - y_c, xb(k) - x_c)
			if(theta_c(k) .lt. 0) theta_c(k) = theta_c(k) + 2.0d0*PI
			
	end do
	!-----------------------------------------------------------------------------------------------
	
	ss = 37
	ee = 37
	open(unit=10, file='normal.dat', status="replace")
		do k = ss, ee
			write(10,*) x(bou_i(k),bou_j(k)),y(bou_i(k),bou_j(k))
			write(10,*) x_ref1(k),y_ref1(k)
			write(10,*) x_ref2(k),y_ref2(k)
			write(10,*) x_ref3(k),y_ref3(k)
			write(10,*) x_ref4(k),y_ref4(k)
		end do
	close(10)
!		open(unit=10, file='neigh2.dat', status="replace")
!		do k = ss, ee
!			write(10,*) x(i_p2(1,k),j_p2(1,k)),y(i_p2(1,k),j_p2(1,k))
!			write(10,*) x(i_p2(2,k),j_p2(2,k)),y(i_p2(2,k),j_p2(2,k))
!			write(10,*) x(i_p2(3,k),j_p2(3,k)),y(i_p2(3,k),j_p2(3,k))
!			write(10,*) x(i_p2(4,k),j_p2(4,k)),y(i_p2(4,k),j_p2(4,k))
!		end do
!	close(10)
!	open(unit=10, file='neigh3.dat', status="replace")
!		do k = ss, ee
!			write(10,*) x(i_p3(1,k),j_p3(1,k)),y(i_p3(1,k),j_p3(1,k))
!			write(10,*) x(i_p3(2,k),j_p3(2,k)),y(i_p3(2,k),j_p3(2,k))
!			write(10,*) x(i_p3(3,k),j_p3(3,k)),y(i_p3(3,k),j_p3(3,k))
!			write(10,*) x(i_p3(4,k),j_p3(4,k)),y(i_p3(4,k),j_p3(4,k))
!		end do
!	close(10)
!	open(unit=10, file='neigh4.dat', status="replace")
!		do k = ss, ee
!			write(10,*) x(i_p4(1,k),j_p4(1,k)),y(i_p4(1,k),j_p4(1,k))
!			write(10,*) x(i_p4(2,k),j_p4(2,k)),y(i_p4(2,k),j_p4(2,k))
!			write(10,*) x(i_p4(3,k),j_p4(3,k)),y(i_p4(3,k),j_p4(3,k))
!			write(10,*) x(i_p4(4,k),j_p4(4,k)),y(i_p4(4,k),j_p4(4,k))
!		end do
!	close(10)
	
	open(unit=10, file='surface2.dat', status="replace")
		do k = 0, nbct-1
			thet = k * (2.0d0 * PI / nbct) 
			x_s = x_c + max_radii*Cos(thet) 
			y_s = y_c + max_radii*Sin(thet)
			write(10,*) x_s,y_s
		end do
	close(10)
	
	
	open(unit=10, file='bound.dat', status="replace")
	do  i = 1, nx
		do  j = 1, ny
			if(isbound(i,j)==1) then
				
					write(10,*) x(i,j), y(i,j)
				
			end if
		end do
	end do	
	close(10)
	
	open(unit=20, file='solid.dat', status="replace")
	do  i = 1, nx
		do  j = 1, ny
			if(issolid(i,j)==1) then
				
					write(20,*) x(i,j), y(i,j)
				
			end if
		end do
	end do	
	close(20)
	
	open(unit=30, file='fluid.dat', status="replace")
	do  i = 1, nx
		do  j = 1, ny
			if(isfluid(i,j)==1) then
				
					write(30,*) x(i,j), y(i,j)
				
			end if
		end do
	end do	
	close(30)
	
  end subroutine get_coord

  subroutine initialize(f, rho, ux, uy)
    real(8), dimension(nx, ny, 0:q-1), intent(out) :: f
    real(8), dimension(nx, ny), intent(out) :: rho, ux, uy
    integer :: i, j, k

    do i = 1, nx
      do j = 1, ny
        rho(i, j) = rho0
        ux(i, j) = 0.0
        uy(i, j) = 0.0
        do k = 0, q-1
          f(i, j, k) = w(k) * rho(i, j)
        end do
      end do
    end do
  end subroutine initialize
  
  subroutine apply_bc(xi,yj,bc_type)
	implicit none
	integer,intent(in) :: xi,yj
	character (len=100), intent(in) :: bc_type
	integer,dimension(0:q-1) :: Is,Os
	real(8) :: rhoI_b,uxI_b,uyI_b,mxxI_b,myyI_b,mxyI_b
	real(8) :: rho_prime_b,Mxx_prime_b,Myy_prime_b,Mxy_prime_b
	
	
	SELECT CASE (bc_type)
		CASE ("outlet_bc")
		  !Right wall (OUTLET):-------------------------------------------------------------------------------------------------------
			rhoI_b = 0.0d0
			uxI_b = 0.0d0
			uyI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 1, 0, 1, 1, 0, 0, 1 /)
			Os = (/ 1, 0, 1, 1, 1, 0, 1, 1, 0 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					uxI_b = uxI_b + f(xi,yj,k)*e(k, 1)
					uyI_b = uyI_b + f(xi,yj,k)*e(k, 2)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			
			uxI_b = uxI_b/rhoI_b
			uyI_b = uyI_b/rhoI_b
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
			
			rho_prime_b = rho(xi,yj)
!			rho_prime_b = (4.0d0*rhoI_b + 3.0d0*rhoI_b*mxxI_b)/(3.0d0 + 3.0d0*ux(xi,yj))
			Mxx_prime_b = (rho_prime_b + 9.0d0*rhoI_b*mxxI_b - 3.0d0*rho_prime_b*ux(xi,yj))/(6.0d0*rho_prime_b)
			Myy_prime_b = 6.0d0*rhoI_b*myyI_b/(5.0d0*rho_prime_b)
			Mxy_prime_b = (6.0d0*rhoI_b*mxyI_b - rho_prime_b*uy(xi,yj))/(3.0d0*rho_prime_b)
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
									& + (3.0d0*uy(xi,yj)*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
			
			
		CASE ("inlet_bc")
			rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 0, 1, 1, 1, 0, 1, 1, 0 /)
			Os = (/ 1, 1, 1, 0, 1, 1, 0, 0, 1 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = (4.0d0*rhoI_b + 3.0d0*rhoI_b*mxxI_b)/(3.0d0 - 3.0d0*ux(xi,yj)) 
			Mxx_prime_b = (rho_prime_b + 9.0d0*rhoI_b*mxxI_b + 3.0d0*rho_prime_b*ux(xi,yj))/(6.0d0*rho_prime_b)
			Myy_prime_b = 6.0d0*rhoI_b*myyI_b/(5.0d0*rho_prime_b)
			Mxy_prime_b = 2.0d0*rhoI_b*mxyI_b/rho_prime_b
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
									& + (3.0d0*uy(xi,yj)*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
			

		CASE ("topwall_bc")
			rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 1, 1, 0, 1, 1, 0, 0 /)
			Os = (/ 1, 1, 0, 1, 1, 0, 0, 1, 1 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = 3.0d0*rhoI_b*( 4.0d0 + 3.0d0*(1.0d0-omega)*myyI_b)/(9.0d0 + omega) 
			Mxx_prime_b = 6.0d0*rhoI_b*mxxI_b/(5.0d0*rho_prime_b)
			Myy_prime_b = (rho_prime_b + 9.0d0*rhoI_b*myyI_b)/(6.0d0*rho_prime_b)
			Mxy_prime_b = 2.0d0*rhoI_b*mxyI_b/rho_prime_b
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
									& + (3.0d0*uy(xi,yj)*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		CASE ("bottomwall_bc")
			rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 0, 1, 1, 0, 0, 1, 1 /)
			Os = (/ 1, 1, 1, 1, 0, 1, 1, 0, 0 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = 3.0d0*rhoI_b*( 4.0d0 + 3.0d0*(1.0d0-omega)*myyI_b)/(9.0d0 + omega) 
			Mxx_prime_b = 6.0d0*rhoI_b*mxxI_b/(5.0d0*rho_prime_b)
			Myy_prime_b = (rho_prime_b + 9.0d0*rhoI_b*myyI_b)/(6.0d0*rho_prime_b)
			Mxy_prime_b = 2.0d0*rhoI_b*mxyI_b/rho_prime_b
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
									& + (3.0d0*uy(xi,yj)*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		CASE DEFAULT
		  print*,"Not a valid boundary type!!!!!!!"
	END SELECT
  
  
  end subroutine apply_bc
  
  subroutine apply_bc_incomp(xi,yj,bc_type)
	implicit none
	integer,intent(in) :: xi,yj
	character (len=100), intent(in) :: bc_type
	integer,dimension(0:q-1) :: Is,Os
	real(8) :: rhoI_b,uxI_b,uyI_b,mxxI_b,myyI_b,mxyI_b
	real(8) :: rho_prime_b,Mxx_prime_b,Myy_prime_b,Mxy_prime_b
	
	
	SELECT CASE (bc_type)
		CASE ("outlet_bc")
		  !Right wall (OUTLET):-------------------------------------------------------------------------------------------------------
			rhoI_b = 0.0d0
			uxI_b = 0.0d0
			uyI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 1, 0, 1, 1, 0, 0, 1 /)
			Os = (/ 1, 0, 1, 1, 1, 0, 1, 1, 0 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					uxI_b = uxI_b + f(xi,yj,k)*e(k, 1)
					uyI_b = uyI_b + f(xi,yj,k)*e(k, 2)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			
			uxI_b = uxI_b/rhoI_b
			uyI_b = uyI_b/rhoI_b
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
			
!			rho_prime_b = 6.0d0*rhoI_b/(5.0d0 + 3.0d0*ux(xi,yj)- 3.0d0*(ux(xi,yj)**2))
			rho_prime_b = rho(xi,yj)
			Mxx_prime_b = ux(xi,yj)**2
			Myy_prime_b = uy(xi,yj)**2
			Mxy_prime_b = (6.0d0*rhoI_b*mxyI_b - rho_prime_b*uy(xi,yj))/(3.0d0*rho_prime_b)
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
									& + (3.0d0*uy(xi,yj)*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
			
			
		CASE ("inlet_bc")
			rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 0, 1, 1, 1, 0, 1, 1, 0 /)
			Os = (/ 1, 1, 1, 0, 1, 1, 0, 0, 1 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = 6.0d0*rhoI_b/(5.0d0 - 3.0d0*ux(xi,yj)- 3.0d0*(ux(xi,yj)**2)) 
			Mxx_prime_b = ux(xi,yj)**2
			Myy_prime_b = 0.0d0
			Mxy_prime_b = 2.0d0*rhoI_b*mxyI_b/rho_prime_b
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
									& + (3.0d0*uy(xi,yj)*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
			

		CASE ("topwall_bc")
			rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 1, 1, 0, 1, 1, 0, 0 /)
			Os = (/ 1, 1, 0, 1, 1, 0, 0, 1, 1 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = 6.0d0*rhoI_b/5.0d0
			Mxx_prime_b = 0.0d0
			Myy_prime_b = 0.0d0
			Mxy_prime_b = 5.0d0*mxyI_b/3.0d0
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
									& + (3.0d0*uy(xi,yj)*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		CASE ("bottomwall_bc")
			rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 0, 1, 1, 0, 0, 1, 1 /)
			Os = (/ 1, 1, 1, 1, 0, 1, 1, 0, 0 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = 6.0d0*rhoI_b/5.0d0
			Mxx_prime_b = 0.0d0
			Myy_prime_b = 0.0d0
			Mxy_prime_b = 5.0d0*mxyI_b/3.0d0
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
									& + (3.0d0*uy(xi,yj)*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		CASE DEFAULT
		  print*,"Not a valid boundary type!!!!!!!"
	END SELECT
  
  
  end subroutine apply_bc_incomp
  
  subroutine boundary_cases(xi,yj,label,uxp,uyp)
	implicit none
	integer,intent(in) :: xi,yj,label
	integer,dimension(0:q-1) :: Is,Os
	real(8),intent(in) :: uxp,uyp
	real(8) :: rhoI_b,mxxI_b,myyI_b,mxyI_b
	real(8) :: rho_prime_b,Mxx_prime_b,Myy_prime_b,Mxy_prime_b
	
	select case (label)
	
		case (0) 
   
		!--------------------------------------------------------------------------------------------------------
		case (1) 
			rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 0, 0, 1, 0, 0, 0, 1 /)
			Os = (/ 1, 0, 1, 1, 0, 0, 1, 0, 0 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
			
			
			rho_prime_b = (12.0d0*rhoI_b*(-3.0d0 - 3.0d0*myyI_b + 3.0d0*mxxI_b*(-1.0d0 + omega) &
							& + 7.0d0*mxyI_b*(-1.0d0 + omega) + 3.0d0*myyI_b*omega)) &
							& / (2.0d0*(-8.0d0 + 7.0d0*uxp - 7.0d0*uyp) + omega*(-9.0d0 + uxp + 15.0d0*(uxp**2) &
							& - uyp + 9.0d0*uxp*uyp + 15.0d0*(uyp**2)))
			
			Mxx_prime_b = -((2.0d0*(-rho_prime_b - 9.0d0*mxxI_b*rhoI_b - 6.0d0*mxyI_b*rhoI_b + 2.0d0*rho_prime_b*uxp &
							& + rho_prime_b*uyp)) / (9.0d0*rho_prime_b))

			Myy_prime_b = (2.0d0*(rho_prime_b + 6.0d0*mxyI_b*rhoI_b + 9.0d0*myyI_b*rhoI_b + rho_prime_b*uxp &
							& + 2.0d0*rho_prime_b*uyp)) / (9.0d0*rho_prime_b)

			Mxy_prime_b = -((-7.0d0*rho_prime_b - 18.0d0*mxxI_b*rhoI_b - 132.0d0*mxyI_b*rhoI_b - 18.0d0*myyI_b*rhoI_b &
							& - 7.0d0*rho_prime_b*uxp + 7.0d0*rho_prime_b*uyp) / (27.0d0*rho_prime_b))
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
									& + (3.0d0*uyp*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		!--------------------------------------------------------------------------------------------------------
		case (2)
		rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 0, 0, 1, 1, 0, 0, 1, 0 /)
			Os = (/ 1, 1, 1, 0, 0, 1, 0, 0, 0 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = (12.0d0*rhoI_b*(-3.0d0 - 3.0d0*myyI_b + 3.0d0*mxxI_b*(-1.0d0 + omega) &
							& - 7.0d0*mxyI_b*(-1.0d0 + omega) + 3.0d0*myyI_b*omega)) &
							& / (-2.0d0*(8.0d0 + 7.0d0*uxp + 7.0d0*uyp) + omega*(-9.0d0 - uxp + 15.0d0*(uxp**2) &
							& - uyp - 9.0d0*uxp*uyp + 15.0d0*(uyp**2)))
			
			Mxx_prime_b = (2.0d0*(rho_prime_b + 9.0d0*mxxI_b*rhoI_b - 6.0d0*mxyI_b*rhoI_b + 2.0d0*rho_prime_b*uxp &
							& - rho_prime_b*uyp)) / (9.0d0*rho_prime_b)

			Myy_prime_b = -((2.0d0*(-rho_prime_b + 6.0d0*mxyI_b*rhoI_b - 9.0d0*myyI_b*rhoI_b + rho_prime_b*uxp &
							& - 2.0d0*rho_prime_b*uyp)) / (9.0d0*rho_prime_b))

			Mxy_prime_b = -((7.0d0*rho_prime_b + 18.0d0*mxxI_b*rhoI_b - 132.0d0*mxyI_b*rhoI_b + 18.0d0*myyI_b*rhoI_b &
							& - 7.0d0*rho_prime_b*uxp - 7.0d0*rho_prime_b*uyp) / (27.0d0*rho_prime_b))
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
									& + (3.0d0*uyp*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		!--------------------------------------------------------------------------------------------------------

		case (3) 
		rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 0, 1, 1, 0, 0, 1, 1 /)
			Os = (/ 1, 1, 1, 1, 0, 1, 1, 0, 0 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = 3.0d0*rhoI_b*( 4.0d0 - 3.0d0*(omega-1.0d0)*myyI_b)&
							& /(9.0d0 + omega + 3.0d0*uyp + 3.0d0*omega*uyp - 6.0d0*omega*(uyp**2)) 
			Mxx_prime_b = 6.0d0*rhoI_b*mxxI_b/(5.0d0*rho_prime_b)
			Myy_prime_b = (rho_prime_b + 9.0d0*rhoI_b*myyI_b + 3.0d0*rho_prime_b*uyp)/(6.0d0*rho_prime_b)
			Mxy_prime_b = (6.0d0*rhoI_b*mxyI_b + rho_prime_b*uxp)/(3.0d0*rho_prime_b)
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
									& + (3.0d0*uyp*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		!--------------------------------------------------------------------------------------------------------
         
		case (4)
			rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 1, 0, 0, 1, 0, 0, 0 /)
			Os = (/ 1, 0, 0, 1, 1, 0, 0, 1, 0 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = (12.0d0*rhoI_b*(-3.0d0 - 3.0d0*myyI_b + 3.0d0*mxxI_b*(-1.0d0 + omega) &
							& - 7.0d0*mxyI_b*(-1.0d0 + omega) + 3.0d0*myyI_b*omega)) &
							& / (2.0d0*(-8.0d0 + 7.0d0*uxp + 7.0d0*uyp) + omega*(-9.0d0 + uxp + 15.0d0*(uxp**2) &
							& + uyp - 9.0d0*uxp*uyp + 15.0d0*(uyp**2)))
			
			Mxx_prime_b = -((2.0d0*(-rho_prime_b - 9.0d0*mxxI_b*rhoI_b + 6.0d0*mxyI_b*rhoI_b + 2.0d0*rho_prime_b*uxp &
							& - rho_prime_b*uyp)) / (9.0d0*rho_prime_b))

			Myy_prime_b = -((2.0d0*(-rho_prime_b + 6.0d0*mxyI_b*rhoI_b - 9.0d0*myyI_b*rhoI_b - rho_prime_b*uxp &
							& + 2.0d0*rho_prime_b*uyp)) / (9.0d0*rho_prime_b))

			Mxy_prime_b = -((7.0d0*rho_prime_b + 18.0d0*mxxI_b*rhoI_b - 132.0d0*mxyI_b*rhoI_b + 18.0d0*myyI_b*rhoI_b &
							& + 7.0d0*rho_prime_b*uxp + 7.0d0*rho_prime_b*uyp) / (27.0d0*rho_prime_b))
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
									& + (3.0d0*uyp*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		!--------------------------------------------------------------------------------------------------------
         

		case (5)
		rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 1, 0, 1, 1, 0, 0, 1 /)
			Os = (/ 1, 0, 1, 1, 1, 0, 1, 1, 0 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = 3.0d0*rhoI_b*( -4.0d0 + 3.0d0*(omega-1.0d0)*mxxI_b)&
							& /(-9.0d0 - omega + 3.0d0*uxp + 3.0d0*omega*uxp + 6.0d0*omega*(uxp**2)) 
			Mxx_prime_b = (rho_prime_b + 9.0d0*rhoI_b*mxxI_b - 3.0d0*rho_prime_b*uxp)/(6.0d0*rho_prime_b)
			Myy_prime_b = 6.0d0*rhoI_b*myyI_b/(5.0d0*rho_prime_b)
			Mxy_prime_b = (6.0d0*rhoI_b*mxyI_b - rho_prime_b*uyp)/(3.0d0*rho_prime_b)
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
									& + (3.0d0*uyp*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		!--------------------------------------------------------------------------------------------------------

		case (6)
	  
		case (7)
		rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 1, 1, 1, 1, 0, 1, 1 /)
			Os = (/ 1, 1, 1, 1, 1, 1, 1, 1, 0 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = (36.0d0*rhoI_b*(-23.0d0 - 3.0d0*myyI_b + 3.0d0*mxxI_b*(-1.0d0 + omega) &
							& - 9.0d0*mxyI_b*(-1.0d0 + omega) + 3.0d0*myyI_b*omega)) &
							& / (6.0d0*(-132.0d0 + 5.0d0*uxp - 5.0d0*uyp) + omega*(-13.0d0 + 69.0d0*(uxp**2) &
							& + 39.0d0*uxp - 207.0d0*uxp*uyp - 39.0d0*uyp + 69.0d0*(uyp**2)))
      							
			Mxx_prime_b = -((-2.0d0*rho_prime_b - 75.0d0*mxxI_b*rhoI_b + 18.0d0*mxyI_b*rhoI_b &
							& - 6.0d0*myyI_b*rhoI_b + 6.0d0*rho_prime_b*uxp - 6.0d0*rho_prime_b*uyp) &
							& / (69.0d0*rho_prime_b))
			Myy_prime_b = -((-2.0d0*rho_prime_b - 75.0d0*myyI_b*rhoI_b + 18.0d0*mxyI_b*rhoI_b &
							& - 6.0d0*mxxI_b*rhoI_b + 6.0d0*rho_prime_b*uxp - 6.0d0*rho_prime_b*uyp) &
							& / (69.0d0*rho_prime_b))
			Mxy_prime_b = -((rho_prime_b + 3.0d0*mxxI_b*rhoI_b - 32.0d0*mxyI_b*rhoI_b + 3.0d0*myyI_b*rhoI_b  &
      						& - 3.0d0*rho_prime_b*uxp + 3.0d0*rho_prime_b*uyp) / (23.0d0*rho_prime_b))
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
									& + (3.0d0*uyp*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		!--------------------------------------------------------------------------------------------------------
	  
		case (8)
		rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 0, 1, 1, 0, 0, 1, 0, 0 /)
			Os = (/ 1, 1, 0, 0, 1, 0, 0, 0, 1 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = (12.0d0*rhoI_b*(-3.0d0 - 3.0d0*myyI_b + 3.0d0*mxxI_b*(-1.0d0 + omega) &
							& + 7.0d0*mxyI_b*(-1.0d0 + omega) + 3.0d0*myyI_b*omega)) &
							& / (-2.0d0*(8.0d0 + 7.0d0*uxp - 7.0d0*uyp) + omega*(-9.0d0 - uxp + 15.0d0*(uxp**2) &
							& + uyp + 9.0d0*uxp*uyp + 15.0d0*(uyp**2)))
			
			Mxx_prime_b = (2.0d0*(rho_prime_b + 9.0d0*mxxI_b*rhoI_b + 6.0d0*mxyI_b*rhoI_b + 2.0d0*rho_prime_b*uxp &
							& + rho_prime_b*uyp)) / (9.0d0*rho_prime_b)

			Myy_prime_b = (2.0d0*(rho_prime_b + 6.0d0*mxyI_b*rhoI_b + 9.0d0*myyI_b*rhoI_b - rho_prime_b*uxp &
							& - 2.0d0*rho_prime_b*uyp)) / (9.0d0*rho_prime_b)

			Mxy_prime_b = -((-7.0d0*rho_prime_b - 18.0d0*mxxI_b*rhoI_b - 132.0d0*mxyI_b*rhoI_b - 18.0d0*myyI_b*rhoI_b &
							& + 7.0d0*rho_prime_b*uxp - 7.0d0*rho_prime_b*uyp) / (27.0d0*rho_prime_b))
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
									& + (3.0d0*uyp*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		!--------------------------------------------------------------------------------------------------------
	  
		case (9)
	  
		case (10)
		rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 0, 1, 1, 1, 0, 1, 1, 0 /)
			Os = (/ 1, 1, 1, 0, 1, 1, 0, 0, 1 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = 3.0d0*rhoI_b*( 4.0d0 - 3.0d0*(omega-1.0d0)*mxxI_b)&
							& /(9.0d0 + omega + 3.0d0*uxp + 3.0d0*omega*uxp - 6.0d0*omega*(uxp**2)) 
			Mxx_prime_b = (rho_prime_b + 9.0d0*rhoI_b*mxxI_b + 3.0d0*rho_prime_b*uxp)/(6.0d0*rho_prime_b)
			Myy_prime_b = 6.0d0*rhoI_b*myyI_b/(5.0d0*rho_prime_b)
			Mxy_prime_b = (6.0d0*rhoI_b*mxyI_b + rho_prime_b*uyp)/(3.0d0*rho_prime_b)
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
									& + (3.0d0*uyp	*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		!--------------------------------------------------------------------------------------------------------
	  
		case (11)
		rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 1, 1, 1, 0, 1, 1, 1 /)
			Os = (/ 1, 1, 1, 1, 1, 1, 1, 0, 1 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = (36.0d0*rhoI_b*(-23.0d0 - 3.0d0*myyI_b + 3.0d0*mxxI_b*(-1.0d0 + omega) &
							& + 9.0d0*mxyI_b*(-1.0d0 + omega) + 3.0d0*myyI_b*omega)) &
							& / (-6.0d0*(132.0d0 + 5.0d0*uxp + 5.0d0*uyp) + omega*(-13.0d0 + 69.0d0*(uxp**2) &
							& - 39.0d0*uxp + 207.0d0*uxp*uyp - 39.0d0*uyp + 69.0d0*(uyp**2)))
      							
			Mxx_prime_b = -((-2.0d0*rho_prime_b - 75.0d0*mxxI_b*rhoI_b - 18.0d0*mxyI_b*rhoI_b &
							& - 6.0d0*myyI_b*rhoI_b - 6.0d0*rho_prime_b*uxp - 6.0d0*rho_prime_b*uyp) &
							& / (69.0d0*rho_prime_b))
			Myy_prime_b = -((-2.0d0*rho_prime_b - 75.0d0*myyI_b*rhoI_b - 18.0d0*mxyI_b*rhoI_b &
							& - 6.0d0*mxxI_b*rhoI_b - 6.0d0*rho_prime_b*uxp - 6.0d0*rho_prime_b*uyp) &
							& / (69.0d0*rho_prime_b))
			Mxy_prime_b = -((-rho_prime_b - 3.0d0*mxxI_b*rhoI_b - 32.0d0*mxyI_b*rhoI_b - 3.0d0*myyI_b*rhoI_b  &
      						& - 3.0d0*rho_prime_b*uxp - 3.0d0*rho_prime_b*uyp) / (23.0d0*rho_prime_b))
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
									& + (3.0d0*uyp*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		!--------------------------------------------------------------------------------------------------------
	  
		case (12)
		rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 1, 1, 0, 1, 1, 0, 0 /)
			Os = (/ 1, 1, 0, 1, 1, 0, 0, 1, 1 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = 3.0d0*rhoI_b*( -4.0d0 + 3.0d0*(omega-1.0d0)*myyI_b)&
							& /(-9.0d0 - omega + 3.0d0*uyp + 3.0d0*omega*uyp + 6.0d0*omega*(uyp**2)) 
			Mxx_prime_b = 6.0d0*rhoI_b*mxxI_b/(5.0d0*rho_prime_b)
			Myy_prime_b = (rho_prime_b + 9.0d0*rhoI_b*myyI_b - 3.0d0*rho_prime_b*uyp)/(6.0d0*rho_prime_b)
			Mxy_prime_b = (6.0d0*rhoI_b*mxyI_b - rho_prime_b*uxp)/(3.0d0*rho_prime_b)
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
									& + (3.0d0*uyp*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		!--------------------------------------------------------------------------------------------------------
	  
		case (13)
		rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 1, 1, 1, 1, 1, 0, 1 /)
			Os = (/ 1, 1, 1, 1, 1, 0, 1, 1, 1 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = (36.0d0*rhoI_b*(-23.0d0 - 3.0d0*myyI_b + 3.0d0*mxxI_b*(-1.0d0 + omega) &
							& + 9.0d0*mxyI_b*(-1.0d0 + omega) + 3.0d0*myyI_b*omega)) &
							& / (6.0d0*(-132.0d0 + 5.0d0*uxp + 5.0d0*uyp) + omega*(-13.0d0 + 69.0d0*(uxp**2) &
							& + 39.0d0*uxp + 207.0d0*uxp*uyp + 39.0d0*uyp + 69.0d0*(uyp**2)))
      							
			Mxx_prime_b = -((-2.0d0*rho_prime_b - 75.0d0*mxxI_b*rhoI_b - 18.0d0*mxyI_b*rhoI_b &
							& - 6.0d0*myyI_b*rhoI_b + 6.0d0*rho_prime_b*uxp + 6.0d0*rho_prime_b*uyp) &
							& / (69.0d0*rho_prime_b))
			Myy_prime_b = -((-2.0d0*rho_prime_b - 75.0d0*myyI_b*rhoI_b - 18.0d0*mxyI_b*rhoI_b &
							& - 6.0d0*mxxI_b*rhoI_b + 6.0d0*rho_prime_b*uxp + 6.0d0*rho_prime_b*uyp) &
							& / (69.0d0*rho_prime_b))
			Mxy_prime_b = -((-rho_prime_b - 3.0d0*mxxI_b*rhoI_b - 32.0d0*mxyI_b*rhoI_b - 3.0d0*myyI_b*rhoI_b  &
      						& + 3.0d0*rho_prime_b*uxp + 3.0d0*rho_prime_b*uyp) / (23.0d0*rho_prime_b))
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
									& + (3.0d0*uyp*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		!--------------------------------------------------------------------------------------------------------
	  
		case (14)
		rhoI_b = 0.0d0
			mxxI_b = 0.0d0
			myyI_b = 0.0d0
			mxyI_b = 0.0d0
			
			Is = (/ 1, 1, 1, 1, 1, 1, 1, 1, 0 /)
			Os = (/ 1, 1, 1, 1, 1, 1, 0, 1, 1 /)
			
			do k = 0, q-1
				if(Is(k)==1) then	
					rhoI_b = rhoI_b + f(xi,yj,k)
					mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
					myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
					mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)
				end if
			end do
			mxxI_b = mxxI_b/rhoI_b
			myyI_b = myyI_b/rhoI_b
			mxyI_b = mxyI_b/rhoI_b
				
			rho_prime_b = (36.0d0*rhoI_b*(-23.0d0 - 3.0d0*myyI_b + 3.0d0*mxxI_b*(-1.0d0 + omega) &
							& - 9.0d0*mxyI_b*(-1.0d0 + omega) + 3.0d0*myyI_b*omega)) &
							& / (-6.0d0*(132.0d0 + 5.0d0*uxp - 5.0d0*uyp) + omega*(-13.0d0 + 69.0d0*(uxp**2) &
							& - 39.0d0*uxp - 207.0d0*uxp*uyp + 39.0d0*uyp + 69.0d0*(uyp**2)))
      							
			Mxx_prime_b = -((-2.0d0*rho_prime_b - 75.0d0*mxxI_b*rhoI_b + 18.0d0*mxyI_b*rhoI_b &
							& - 6.0d0*myyI_b*rhoI_b - 6.0d0*rho_prime_b*uxp + 6.0d0*rho_prime_b*uyp) &
							& / (69.0d0*rho_prime_b))
			Myy_prime_b = -((-2.0d0*rho_prime_b - 75.0d0*myyI_b*rhoI_b + 18.0d0*mxyI_b*rhoI_b &
							& - 6.0d0*mxxI_b*rhoI_b - 6.0d0*rho_prime_b*uxp + 6.0d0*rho_prime_b*uyp) &
							& / (69.0d0*rho_prime_b))
			Mxy_prime_b = -((rho_prime_b + 3.0d0*mxxI_b*rhoI_b - 32.0d0*mxyI_b*rhoI_b + 3.0d0*myyI_b*rhoI_b  &
      						& + 3.0d0*rho_prime_b*uxp - 3.0d0*rho_prime_b*uyp) / (23.0d0*rho_prime_b))
			
			do k = 0, q-1
				f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
									& + (3.0d0*uyp*e(k, 2))  &
									& + (9.0d0*Mxx_prime_b*Hxx(k)/2.0d0) &
									& + (9.0d0*Myy_prime_b*Hyy(k)/2.0d0) &
									& + (9.0d0*Mxy_prime_b*Hxy(k)) )						
			end do
		!--------------------------------------------------------------------------------------------------------
	  
		case (15)
	  
           

      case default
          
         
   end select
	
  end subroutine boundary_cases

   
  subroutine create_output_dir_if_needed (dir_name, length)
    implicit none
    integer, intent(in) :: length
    character (len=length), intent(in) :: dir_name
    character :: delimiter
    logical   :: dir_exists = .false. 

    call get_environment_variable ('DELIMITER', delimiter)
    inquire(file=dir_name,exist=dir_exists)
    if (dir_exists) then 
       write (*, *) 'OUTPUT directory exists...'
    else 
       write (*, *) 'NO OUTPUT directory exists, creating one...'
       call system ('mkdir '//dir_name)
    end if

  end subroutine create_output_dir_if_needed
  
  subroutine write_grid_plot3d()
    implicit none
    integer :: i, j, m 
    integer :: nprocs
    
    ! Compute grid spacing
	  dx = 1.0
	  dy = 1.0

	nprocs = 1
	 do j = 1, ny
     yvar(:,j) = (j - 1) * dy
     do i = 1, nx
        xvar(i,:) = (i - 1) * dx
     end do
  end do

    open(unit = 100, form = 'unformatted', file = output_dir_name//'/'//trim('grid.x'))
    write(100) nprocs
    write(100) (nx, ny, m = 1, nprocs)
    do m = 0, nprocs-1
       write (100) &
            &  (( sngl(xvar(i, j)), i = 1, nx), j = 1, ny), &
            &  (( sngl(yvar(i, j)), i = 1, nx), j = 1, ny)
    end do

    close(100)

  end subroutine write_grid_plot3d
  
  
subroutine write_function_plot3d(filename)
    implicit none

    integer :: i, j, l, m 
    integer :: nblocks,nprocs
	character(len=100), intent(in) :: filename

  
	nprocs = 1
    nblocks = nprocs

    open(unit = 200, form = 'unformatted', file=output_dir_name//'/'//trim(filename))
    write(200) nblocks
    write(200) (nx, ny, 3, l = 0, nprocs-1)
    do m = 0, nprocs-1
       write(200) &
            & (( sngl(rho(i, j)), i = 1, nx), j = 1, ny), &
            & (( sngl(ux(i, j)), i = 1, nx), j = 1, ny), &
            & (( sngl(uy(i, j)), i = 1, nx), j = 1, ny)
    end do

    close(200)

  end subroutine write_function_plot3d
  
  subroutine write_restart_file(soln_time)
    implicit none
    integer, intent(in) :: soln_time 
    integer :: iunit
    character(len=25) :: file_name

    iunit = 10000
    file_name = restart_dir_name//'/'//'restart.dat'
    write(*, *) "Writing restart file at",soln_time,"iteration"

    open(unit = iunit, file = file_name, status = "unknown", form="unformatted")

    write(iunit) soln_time, f, f_eq, ux, uy, rho

    close(iunit)

  end subroutine write_restart_file

  subroutine read_restart_file(soln_time)
    implicit none
    integer, intent(out) :: soln_time 
    integer :: iunit
    character(len=25) :: file_name

    iunit = 10000
    file_name = restart_dir_name//'/'//'restart.dat'
    write(*, *) "Reading restart file from",soln_time,"iteration"

    open(unit = iunit, file = file_name, status = "unknown", form="unformatted")

    read(iunit) soln_time, f, f_eq, ux, uy, rho

    close(iunit)

  end subroutine read_restart_file
  
  subroutine read_input_file()
    implicit none
    integer :: input = 100
    integer :: iread_error = 0
    integer :: i 

    namelist/DomainSize/nx,ny,lattice_type
    namelist/Cylinder/ncy,L_front,L_back,L_top,L_bot
    namelist/Numbers/Re
    namelist/Parallel/nprocsx,nprocsy
    namelist/Controls/uo,iplot,max_iter,isave,irestart,statsbegin,statsend,iplotstats
    namelist/LogicalControls/x_periodic,y_periodic,channel_with_cylinder,channel_with_square,incomp,vel_interp

300 format("Error while reading input.dat file...")
150 if (iread_error .ne. 0) then 
       write(*, 300); stop 
    end if
200 format(a, T5)

    open (unit = input, file = "input.dat", status = "old")
    write(*, *) 'Reading input.dat file ...'

    read(input, DomainSize, iostat=iread_error, err=150)
    write(*, 200) 'DomainSize'
    
    read(input, Cylinder, iostat=iread_error, err=150)
    write(*, 200) 'Cylinder'

    read(input, Numbers, iostat=iread_error, err=150)
    write(*, 200) 'Numbers'

    read(input, Parallel, iostat=iread_error, err=150)
    write(*, 200) 'Parallel'

    read(input, Controls, iostat=iread_error, err=150)
    write(*, 200) 'Controls'
    
    read(input, LogicalControls, iostat=iread_error, err=150)
    write(*, 200) 'LogicalControls'

    close(input)

    if (verbose) then 
       write(*, DomainSize)
       write(*, Cylinder)
       write(*, Numbers)
       write(*, Parallel)
       write(*, Controls)
       write(*, LogicalControls)
    end if

    write(*, '(a40, 5x, g15.5)') 'Reynolds number (Re) is ', Re 

  end subroutine read_input_file
  
  subroutine allocate_memory()
  	implicit none
  	
  	!Allocating variables
  	allocate(isfluid(1:nx,1:ny), issolid(1:nx,1:ny),isbound(1:nx,1:ny))
  	allocate(x(1:nx,1:ny),y(1:nx,1:ny),e(0:q-1, 1:2))
  	allocate(xvar(1:nx,1:ny),yvar(1:nx,1:ny))
  	allocate(w(0:q-1),Hxx(0:q-1),Hyy(0:q-1),Hxy(0:q-1),Hxxy(0:q-1),Hxyy(0:q-1),Hxxyy(0:q-1))
  	allocate(f(1:nx, 1:ny, 0:q-1), f_eq(1:nx, 1:ny, 0:q-1), f_tmp(1:nx, 1:ny, 0:q-1))
  	allocate(fin(1:nx, 1:ny, 0:q-1),fout(1:nx, 1:ny, 0:q-1))
  	allocate(rho(1:nx, 1:ny),p(1:nx, 1:ny), ux(1:nx, 1:ny), uy(1:nx, 1:ny))
  	allocate(mxx(1:nx, 1:ny), myy(1:nx, 1:ny), mxy(1:nx, 1:ny))
  	allocate(er(1:nx, 1:ny),er1(1:nx, 1:ny))
  	allocate(scalar(1:nx, 1:ny))
  
  
  end subroutine allocate_memory
  
  subroutine allocate_cylinder_memory()
  	implicit none
  	
  	!Allocating variables after the get_coord subroutine execution:
  	allocate(bou_i(1:nbct),bou_j(1:nbct))
	allocate(label_bc(1:nbct),theta_c(1:nbct))
	allocate(p_mean_cy(1:nbct))
	allocate(ut(1:nbct),un(1:nbct))
	allocate(radii(1:nbct))
	allocate(xb(1:nbct),yb(1:nbct))
	allocate(delta_uk(1:nbct),pb(1:nbct),p_fit(1:nbct))
	allocate(x_ref1(1:nbct),y_ref1(1:nbct))
	allocate(x_ref2(1:nbct),y_ref2(1:nbct))
	allocate(x_ref3(1:nbct),y_ref3(1:nbct))
	allocate(x_ref4(1:nbct),y_ref4(1:nbct))
	allocate(i_p2(1:4,1:nbct),j_p2(1:4,1:nbct))
	allocate(i_p3(1:4,1:nbct),j_p3(1:4,1:nbct))
	allocate(i_p4(1:4,1:nbct),j_p4(1:4,1:nbct))
	allocate(p_ref1(1:nbct),ux_ref1(1:nbct),uy_ref1(1:nbct))
	allocate(p_ref2(1:nbct),ux_ref2(1:nbct),uy_ref2(1:nbct))
	allocate(p_ref3(1:nbct),ux_ref3(1:nbct),uy_ref3(1:nbct))
	allocate(p_ref4(1:nbct),ux_ref4(1:nbct),uy_ref4(1:nbct))

  end subroutine allocate_cylinder_memory
  
  subroutine create_master_p3d_file()
  	implicit none

  character(len=100) :: filename,filename_bin
  integer(8) :: i,j,istart,iend,istep,iter
  integer(8) :: soln_time
  real(8) :: phy_time
  
  istart = 0
  iend = 10000*iplot

  1005 format (a4,I8,a4)
  
  open(unit=10,file='master.p3d')
  write(10,*)"{"
  write(10,*)
  write(10,*) '"auto-detect-format"', ": true,"
  write(10,*)
   write(10,*)  '"filenames"', ": ["
   write(10,*)
   
  do iter = istart,iend,iplot
	soln_time = iter/iplot
	phy_time = iter*delt_phy
	write (filename_bin, 1005) 'grid',10000000+iter,'.bin'
	write(10,*) '{ "time" :',phy_time,',"xyz" : "grid.x", "function" :' ,' "',trim(filename_bin),'" },'
  end do
  write(10,*)
  write(10,*) "]"
  write(10,*)"}"
  close(10)

  
  end subroutine create_master_p3d_file

end program lbm_2d
