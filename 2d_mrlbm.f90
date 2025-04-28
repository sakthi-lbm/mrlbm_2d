program lbm_2d
    implicit none
    ! Parameters
    real(8), parameter :: PI = 4.0d0*atan(1.0d0)
    integer :: num_threads=8
    integer :: nprocsx, nprocsy
    integer :: iplot, max_iter, isave, irestart, statsbegin,statsend, iplotstats, cycle_period
    integer :: restart_step = 1

    !geometry and grid variables
    integer, allocatable,dimension(:,:) :: isfluid, issolid,isbound
    real(8), allocatable,dimension(:,:) :: x,y,xvar,yvar
    integer :: nbct
    real(8), allocatable,dimension(:,:) :: boundary
    integer, allocatable,dimension(:) :: bou_i,bou_j,label_bc,idx
    integer :: nx,ny
    integer :: ncy, L_front,L_back,L_top,L_bot
    integer :: nfront,nback,ntop,nbot
    real(8) :: x_c,y_c,r_c,distance
    real(8) :: L_phy, nu_phy, u_phy, delx_phy, delt_phy
    real(8) :: L_lat, nu_lat, u_lat, delx_lat, delt_lat
    real(8) :: error1, frob1,max_radii, r_cyl,p_int, ul2norm 
    real(8), allocatable,dimension(:, :) :: er,er1

    !LBM variables
    integer, parameter :: q = 9
    integer :: lattice_type 
    real(8), parameter :: rho0 = 1.0
    real(8) :: Re,uo,nu,omega, uinlet,vinlet
    real(8), allocatable,dimension(:) :: w,Hxx,Hyy,Hxy,Hxxy,Hxyy,Hxxyy
    real(8) :: cu
    real(8) :: dx, dy, delx
    real(8), allocatable,dimension(:, :, :) :: cx_p, cy_p, Hxx_p,Hyy_p,Hxy_p
    integer, allocatable,dimension(:, :) :: e, e_p
    real(8), allocatable,dimension(:, :, :) :: f, f_eq, f_tmp,fin,fout
    real(8), allocatable, dimension(:, :) :: rho, p, ux, uy, mxx, myy, mxy,scalar
    real(8), allocatable, dimension(:, :) :: ux_prev, uy_prev
    real(8), allocatable, dimension(:) ::  cth_cyl, sth_cyl, cths_cyl, sths_cyl
    real(8), allocatable, dimension(:, :) ::  cth_glb, sth_glb, c2th_glb, s2th_glb
    real(8), allocatable, dimension(:) :: uxp_b, uyp_b
    real(8), allocatable, dimension(:,:) :: mxxp, myyp, mxyp

    !Statistics variables
    real(8), allocatable,dimension(:) :: theta_c,p_mean_cy, tw_mean_cy,radii,xb,yb,delta_uk,pb,p_fit, tw_fit
    real(8), allocatable,dimension(:) :: x_ref1,y_ref1,x_ref2,y_ref2,x_ref3,y_ref3,x_ref4,y_ref4
    real(8), allocatable,dimension(:) :: ux_ref1,uy_ref1,ux_ref2,uy_ref2,ux_ref3,uy_ref3,ux_ref4,uy_ref4
    real(8), allocatable,dimension(:) :: var_ref1,var_ref2,var_ref3,var_ref4
    real(8), allocatable,dimension(:) :: mxxp_ref1,mxxp_ref2,mxxp_ref3,mxxp_ref4
    real(8), allocatable,dimension(:) :: myyp_ref1,myyp_ref2,myyp_ref3,myyp_ref4
    integer, allocatable,dimension(:,:) :: i_p2,j_p2,i_p3,j_p3,i_p4,j_p4
    integer, allocatable,dimension(:,:) :: Incs_b, Outs_b
    real(8), allocatable,dimension(:) :: thetas_cy, ps_cy, mxyps_cy, tws_cy
    real(8), allocatable,dimension(:) :: theta_cy, p_cy, mxyp_cy, tw_cy
    real(8), allocatable,dimension(:) :: F_drag, F_lift
    real(8), allocatable,dimension(:) :: p_costheta, p_sintheta, tw_costheta, tw_sintheta
    real(8) :: rho_inf_mean, p_inf_mean, mean_counter, rho_mean_counter, vis_drag, vis_lift, pres_drag, pres_lift
    real(8) :: vis_drag_mean, vis_lift_mean, pres_drag_mean, pres_lift_mean
    real(8) :: F_drag_mean, F_lift_mean, C_drag_mean, C_lift_mean, force_norm, p_norm
    integer :: i, j, k, l, iter, coord_i,coord_j
    integer :: i_probe1, j_probe1, i_probe2, j_probe2, i_probe3, j_probe3, i_probe4, j_probe4

    !logics
    character (len=6) :: output_dir_name = 'output'//char(0)
    character (len=7) :: restart_dir_name = 'restart'//char(0)
    character (len=10) :: probe_dir_name = 'data_probe'//char(0)
    character (len=8) :: geo_dir_name = 'data_geo'//char(0)
    character (len=9) :: mean_dir_name = 'data_mean'//char(0)
    character(len=100) :: filename,filename_bin
    character (len=100) :: bc_type
    logical :: x_periodic, y_periodic, channel_with_cylinder,incomp,channel_with_square
    logical :: vel_interp, mom_interp, rotated_coordinate, post_process
    logical :: pop_collision, mom_collision
    logical, parameter :: verbose = .true.

    pop_collision = .NOT. mom_collision

    !$ call omp_set_num_threads(num_threads)

    call create_output_dir_if_needed(output_dir_name,6)
    call create_output_dir_if_needed(restart_dir_name,7)
    call create_output_dir_if_needed(probe_dir_name,10)
    call create_output_dir_if_needed(geo_dir_name,8)
    call create_output_dir_if_needed(mean_dir_name,9)

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


    !Probes
    !=========================================================================================================
    i_probe1 = (nfront + ncy + 1) + (ncy/2)
    j_probe1 = (nbot + (ncy/2) + 1)
    i_probe2 = (nfront + ncy + 1) + (ncy)
    j_probe2 = (nbot + (ncy/2) + 1)
    i_probe3 = (nfront + ncy + 1) + (3*ncy/2)
    j_probe3 = (nbot + (ncy/2) + 1)
    i_probe4 = (nfront + ncy + 1) + (2*ncy)
    j_probe4 = (nbot + (ncy/2) + 1)

    open(unit=100,file="data_probe/probe_loc.dat")
        write(100,*) "probe 1: ", i_probe1, j_probe1
        write(100,*) "probe 2: ", i_probe2, j_probe2
        write(100,*) "probe 3: ", i_probe3, j_probe3
        write(100,*) "probe 4: ", i_probe4, j_probe4
    close(100)
    !==============================================================================================================

    dx = 1.0d0
    dy = 1.0d0

    !Flow constants: 
    uinlet  = uo         !Inlet axial velocity
    vinlet  = 0.0d0     !inlet radial velocity 

    ! Lattice weights and velocity sets
    w(0:q-1) = (/ 16.0/36.0, 4.0/36.0, 4.0/36.0, 4.0/36.0, 4.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 /)  !w_k
    e(0:q-1, 1) = (/ 0, 1, 0, -1, 0, 1, -1, -1, 1 /)    !cx
    e(0:q-1, 2) = (/ 0, 0, 1, 0, -1, 1, 1, -1, -1 /)    !cy


    nu    = uinlet * (ny-1)/Re
    if (channel_with_cylinder) nu    = uinlet *(2.0d0*r_cyl)/Re;    !kinematic viscosity 
    omega = 1.0d0/((3.0*nu)+(1.0d0/2.0d0)) !relaxation parameter

    L_lat = ny-1
    if (channel_with_cylinder) L_lat = 2.0d0*r_cyl

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
        write(10,*) 'No. of points on diameter', ncy
        write(10,*) 'centre_cylin', nfront + (ncy/2),nbot + ncy/2 
        write(10,*) 'Maximum_radius', r_cyl
        write(10,*) 'No. of points on cylinder:', nbct
        write(10,*) 'Re:', Re
        write(10,*) 'uo:', uo
        write(10,*) 'nu:', nu
        write(10,*) 'tau:', 1.0d0/omega
        write(10,*) 'dt_phy:', delt_phy
    close(10)

    ! Second-order Hermite polynomial 
    do k = 0, q-1
        Hxx(k) = e(k, 1)*e(k, 1) - (1.0d0/3.0d0)    !Hxx
        Hyy(k) = e(k, 2)*e(k, 2) - (1.0d0/3.0d0)    !Hyy
        Hxy(k) = e(k, 1)*e(k, 2)                    !Hxy
        Hxxy(k) = (e(k, 1)*e(k, 1) - (1.0d0/3.0d0))*e(k, 2)                                 !Hxxy
        Hxyy(k) = (e(k, 2)*e(k, 2) - (1.0d0/3.0d0))*e(k, 1)                                 !Hxxy
        Hxxyy(k) = (e(k, 1)*e(k, 1) - (1.0d0/3.0d0))*(e(k, 2)*e(k, 2) - (1.0d0/3.0d0))      !Hxxy
    end do

    ! Lattice Velocities and Second-order Hermite polynomial in rotated coordinates
    do i = 1, nx
        do j = 1, ny
            do k = 0, q-1
                cx_p(i, j, k) = e(k, 1) * cth_glb(i,j) + e(k, 2) * sth_glb(i,j)     !cx'
                cy_p(i, j, k) = -e(k, 1) * sth_glb(i,j) + e(k, 2) * cth_glb(i,j)    !cy'

                Hxx_p(i, j, k) = Hxx(k) * (cth_glb(i,j)**2) + Hyy(k) *(sth_glb(i,j)**2) + Hxy(k) * s2th_glb(i,j)    !Hxx'
                Hyy_p(i, j, k) = Hxx(k) * (sth_glb(i,j)**2) + Hyy(k) *(cth_glb(i,j)**2) - Hxy(k) * s2th_glb(i,j)    !Hyy'
                Hxy_p(i, j, k) = (Hyy(k) - Hxx(k)) * 0.50d0 * s2th_glb(i,j) + Hxy(k) * c2th_glb(i,j)                !Hxy'
            end do
        end do
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

    call finding_incoming_outgoing_pops()

    open(unit=505,file='data_mean/lift_drag.dat')

    ! Main loop
  do iter = restart_step+1, max_iter

    !Macroscopic variables:
    !$omp parallel do private(i,j) shared( f, e, rho, ux, uy) schedule(guided) 
    do i = 1, nx
        do j = 1, ny
            rho(i, j) = sum(f(i, j, :))
            ux(i, j) = sum(f(i, j, :) * e(:, 1)) / rho(i, j)
            uy(i, j) = sum(f(i, j, :) * e(:, 2)) / rho(i, j)
            mxx(i, j) =  sum(f(i, j, :) * Hxx(:)) / rho(i, j)
            myy(i, j) =  sum(f(i, j, :) * Hyy(:)) / rho(i, j)
            mxy(i, j) =  sum(f(i, j, :) * Hxy(:)) / rho(i, j)
        end do
    end do
    !$omp end parallel do

    !Free-density: Mean inlet density
    rho_mean_counter = real(iter - restart_step)
    rho_inf_mean = ((rho_mean_counter*rho_inf_mean) + rho(1,ny/2) )/(rho_mean_counter + 1.0d0)
    p_inf_mean = rho_inf_mean/3.0d0

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

    if(channel_with_cylinder .and. (.not. rotated_coordinate)) then
        !velocity interpolation for curved boundary: 
        if(vel_interp) call velocity_interpolation_boundary_nodes()

        do l = 1,nbct
            i = bou_i(l)
            j = bou_j(l)

            !Original coordinate system:
            call numerical_boundary_cases(bou_i(l),bou_j(l),label_bc(l),ux(i,j),uy(i,j))

        end do
    end if

    if(channel_with_cylinder .and. rotated_coordinate) then
        !velocity interpolation for curved boundary: 
        if(vel_interp) call velocity_interpolation_boundary_nodes()
        if(mom_interp) call moments_interpolation_boundary_nodes()

        do l = 1,nbct
            i = bou_i(l)
            j = bou_j(l)

            !!Rotated coordinate system
            uxp_b(l) = ux(i,j)*cth_glb(i,j) + uy(i,j)*sth_glb(i,j)
            uyp_b(l) = uy(i,j)*cth_glb(i,j) - ux(i,j)*sth_glb(i,j)
            call numerical_boundary_cases_rotation(bou_i(l),bou_j(l),label_bc(l),uxp_b(l),uyp_b(l))
            !call numerical_boundary_cases_rotation_rhoeq(bou_i(l),bou_j(l),label_bc(l),uxp_b(l),uyp_b(l))
            !call numerical_boundary_cases_rotation_weak(bou_i(l),bou_j(l),label_bc(l),uxp_b(l),uyp_b(l))
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
    !Collision - velocity space
    !----------------------------------------------------------------------------------------------------	
    if(pop_collision) then
        ! do i = 1, nx
            ! do j = 1, ny
                ! rho(i, j) = sum(f(i, j, :))
            ! end do
        ! end do

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

        !$omp parallel do collapse(2) shared(rho, mxx,myy,mxy, Hxx,Hyy,Hxy, fout)
        do i = 1, nx
            do j = 1, ny 
                mxx(i, j) =  sum(fout(i, j, :) * Hxx(:)) / rho(i, j)
                myy(i, j) =  sum(fout(i, j, :) * Hyy(:)) / rho(i, j)
                mxy(i, j) =  sum(fout(i, j, :) * Hxy(:)) / rho(i, j)
            end do
        end do
        !$omp end parallel do
    end if


    !----------------------------------------------------------------------------------------------------	
    !Collision - Moments space
    !----------------------------------------------------------------------------------------------------	
    if(mom_collision)then
        !$omp parallel do collapse(2) shared(rho, mxx,myy,mxy, Hxx,Hyy,Hxy, fout)
        do i = 1, nx
            do j = 1, ny 
                mxx(i, j) =  (1.0d0 - omega)*mxx(i, j) + omega*ux(i, j)*ux(i, j)
                myy(i, j) =  (1.0d0 - omega)*myy(i, j) + omega*uy(i, j)*uy(i, j)
                mxy(i, j) =  (1.0d0 - omega)*mxy(i, j) + omega*ux(i, j)*uy(i, j)
            end do
        end do
        !$omp end parallel do
    end if


    !----------------------------------------------------------------------------------------------------	
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

    !Co-efficients Calculation:
    if(channel_with_cylinder) then
        do l = 1,nbct
            i = bou_i(l)
            j = bou_j(l)
            call evaluate_forces(bou_i(l),bou_j(l),label_bc(l), F_drag(l), F_lift(l))
        end do
    end if
    !Coefficients:
    F_drag_mean = sum(F_drag(:))
    F_lift_mean = sum(F_lift(:))
    force_norm = r_cyl*rho_inf_mean*(uo**2)
    C_drag_mean = F_drag_mean/force_norm
    C_lift_mean = F_lift_mean/force_norm

    write(505, *) iter, F_drag_mean, F_lift_mean, C_drag_mean, C_lift_mean


    !----------------------------------------------------------------------------------------------------	
    ! Streaming step
    !----------------------------------------------------------------------------------------------------	
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


    !writing output files and statistics:
    if(post_process) then

        if(iter == statsbegin) open(unit=101,file="data_probe/p_probe.dat")
        if(iter == statsbegin) open(unit=102,file="data_probe/ux_probe.dat")
        if(iter == statsbegin) open(unit=103,file="data_probe/uy_probe.dat")

        write(101,*) iter, rho(i_probe1,j_probe1)/3.0d0, rho(i_probe2,j_probe2)/3.0d0, &
            & rho(i_probe3,j_probe3)/3.0d0, rho(i_probe4,j_probe4)/3.0d0 
        write(102,*) iter, ux(i_probe1,j_probe1)/3.0d0, ux(i_probe2,j_probe2)/3.0d0, &
            & ux(i_probe3,j_probe3)/3.0d0, ux(i_probe4,j_probe4)/3.0d0 
        write(103,*) iter, uy(i_probe1,j_probe1)/3.0d0, uy(i_probe2,j_probe2)/3.0d0, &
            & uy(i_probe3,j_probe3)/3.0d0, uy(i_probe4,j_probe4)/3.0d0

        if(iter == statsend) close(101)	
        if(iter == statsend) close(102)	
        if(iter == statsend) close(103)	


        !Statistics
        !=======================================================================================================

        if((iter .ge. statsbegin) .and. (iter .le. statsend)) then

            call time_averaging_statistics(iter)

            if(mod((iter-statsbegin),cycle_period)==0) then

                call write_statistics_out()

            end if

        end if

        !=======================================================================================================

        !output binary files	
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

    end if      !write output condition

    call calculate_norm(ul2norm)

    if (mod(iter, 100) == 0) then 
        write(*, *) "STEP: ",iter ,"     ", "VEL_NORM: ",  ul2norm 
    end if


    if (mod(iter, isave) == 0) then 
        call write_restart_file(iter)
    end if


 end do         !Main loop end

    close(505)
contains

    subroutine evaluate_forces(xi,yj,label, F_x, F_y)
        implicit none
        integer,intent(in) :: xi,yj,label
        real(8),intent(out) :: F_x, F_y
        integer,dimension(0:q-1) :: Is,Os
        real(8) :: Fx_inc, Fy_inc, Fx_out, Fy_out

        Is(0:q-1) = Incs_b(label,0:q-1)
        Os(0:q-1) = Outs_b(label,0:q-1)

        Fx_inc  = 0.0d0
        Fy_inc  = 0.0d0
        Fx_out  = 0.0d0 
        Fy_out  = 0.0d0

        do k = 0, q-1
            if(Os(k)==1) then
                Fx_out  = Fx_out + fout(xi, yj, k)*e(k, 1)
                Fy_out  = Fy_out + fout(xi, yj, k)*e(k, 2)
            end if

            if(Is(k)==1) then
                Fx_inc  = Fx_inc + f(xi, yj, k)*e(k, 1)
                Fy_inc  = Fy_inc + f(xi, yj, k)*e(k, 2)
            end if
        end do

        F_x = Fx_out + Fx_inc
        F_y = Fy_out + Fy_inc

    end subroutine evaluate_forces

    subroutine write_statistics_out()
        implicit none
        integer :: i, j, k

        do k = 1, nbct
            p_costheta(k)  = p_mean_cy(k) * cths_cyl(k)
            tw_sintheta(k) = tw_mean_cy(k) * sths_cyl(k)
            p_sintheta(k)  = p_mean_cy(k) * sths_cyl(k)
            tw_costheta(k) = tw_mean_cy(k) * cths_cyl(k)
        end do

        call trapezoidal_sub(nbct,thetas_cy,p_costheta, pres_drag)
        call trapezoidal_sub(nbct,thetas_cy,tw_sintheta, vis_drag)
        call trapezoidal_sub(nbct,thetas_cy,p_sintheta, pres_lift)
        call trapezoidal_sub(nbct,thetas_cy,tw_costheta, vis_lift)

        !Force calculation:
        F_drag_mean = r_cyl*(-pres_drag + vis_drag)
        F_lift_mean = r_cyl*(-pres_lift - vis_lift)

        !print*,iter, F_drag_mean, F_lift_mean

        !Coefficients:
        force_norm = r_cyl*rho_inf_mean*(uo**2)
        p_norm = 0.50d0*rho_inf_mean*(uo**2)
        C_drag_mean = F_drag_mean/force_norm
        C_lift_mean = F_lift_mean/force_norm

        open(unit=201,file="data_mean/p_coeff.dat")
            do k = 1, nbct
                if(thetas_cy(k) .le. PI) then
                    write(201,*) abs((thetas_cy(k)*180.0d0/PI)-180.0d0), (p_mean_cy(k) - p_inf_mean)/p_norm, &
                        &  (p_fit(k) - p_inf_mean)/p_norm
                end if
            end do
        close(201)
        open(unit=202,file="data_mean/skin_friction.dat")
            do k = 1, nbct
                if(thetas_cy(k) .le. PI) then
                    write(202,*) abs((thetas_cy(k)*180.0d0/PI)-180.0d0), tw_mean_cy(k)/p_norm, tw_fit(k)/p_norm
                end if
            end do
        close(202)

        open(unit=203,file="data_mean/coefficients.dat")
            write(203,*) "Drag, C_D:", C_drag_mean
            write(203,*) "Lift, C_L:", C_lift_mean
        close(203)

    end subroutine write_statistics_out

    subroutine time_averaging_statistics(step)
        implicit none
        integer, intent(in) :: step
        integer :: i,j,k

        mean_counter = real(step - statsbegin)

        do i = 1, nx
            do j = 1, ny
                rho(i, j) =  sum(f(i, j, :))
                p(i, j) = rho(i, j)/3.0d0
                mxyp(i, j) =  sum(f(i, j, :) * Hxy_p(i, j, :)) / rho(i, j)
            end do
        end do

        call extrapolate_var_on_cylinder(p, p_cy)
        call extrapolate_var_on_cylinder(mxyp, mxyp_cy)

        do k = 1,nbct
            ps_cy(k) = p_cy(idx(k))
            mxyps_cy(k) = mxyp_cy(idx(k))
            cths_cyl(k) = cth_cyl(idx(k))
            sths_cyl(k) = sth_cyl(idx(k))

            !wall shear stress
            tw_cy(k) = (-3.0d0 * nu * omega) * mxyps_cy(k)
            !tw_cy(k) = (-9.0d0 * nu * omega) * ps_cy(k) * mxyps_cy(k)

            !Mean pressure and mean shear stress calculation
            p_mean_cy(k) = ((mean_counter*p_mean_cy(k)) + ps_cy(k) )/(mean_counter + 1.0d0)
            tw_mean_cy(k) = ((mean_counter*tw_mean_cy(k)) + tw_cy(k) )/(mean_counter + 1.0d0)
        end do

        call gaussian_smoothing(nbct,thetas_cy,p_mean_cy,p_fit)
        call gaussian_smoothing(nbct,thetas_cy,tw_mean_cy,tw_fit)

    end subroutine time_averaging_statistics

    subroutine extrapolate_var_on_cylinder(var_int, int_var_cy)
        implicit none
        real(8),dimension(1:nx,1:ny),intent(in) :: var_int
        real(8),dimension(1:nbct),intent(out) :: int_var_cy
        real(8) :: x1,y1,var1,x2,y2,var2,x3,y3,var3,x4,y4,var4
        real(8) :: xd,yd,var_0,var_1
        real(8) :: xs,ys,vars
        integer :: i,j,k

        do k = 1,nbct
            i = bou_i(k)
            j = bou_j(k)

            !Second reference node: First fluid point
            xd = (x_ref2(k) - x(i_p2(1,k),j_p2(1,k)))/dx
            yd = (y_ref2(k) - y(i_p2(1,k),j_p2(1,k)))/dy

            var_0 = (1.0d0 - xd)*var_int(i_p2(1,k),j_p2(1,k)) + xd*var_int(i_p2(2,k),j_p2(2,k))
            var_1 = (1.0d0 - xd)*var_int(i_p2(3,k),j_p2(3,k)) + xd*var_int(i_p2(4,k),j_p2(4,k))

            var_ref2(k) =  (1.0d0 - yd)*var_0 + yd*var_1

            !Third reference node: Second fluid point
            xd = (x_ref3(k) - x(i_p3(1,k),j_p3(1,k)))/dx
            yd = (y_ref3(k) - y(i_p3(1,k),j_p3(1,k)))/dy

            var_0 = (1.0d0 - xd)*var_int(i_p3(1,k),j_p3(1,k)) + xd*var_int(i_p3(2,k),j_p3(2,k))
            var_1 = (1.0d0 - xd)*var_int(i_p3(3,k),j_p3(3,k)) + xd*var_int(i_p3(4,k),j_p3(4,k))

            var_ref3(k) =  (1.0d0 - yd)*var_0 + yd*var_1

            !Fourth reference node: Third fluid point
            xd = (x_ref4(k) - x(i_p4(1,k),j_p4(1,k)))/dx
            yd = (y_ref4(k) - y(i_p4(1,k),j_p4(1,k)))/dy

            var_0 = (1.0d0 - xd)*var_int(i_p4(1,k),j_p4(1,k)) + xd*var_int(i_p4(2,k),j_p4(2,k))
            var_1 = (1.0d0 - xd)*var_int(i_p4(3,k),j_p4(3,k)) + xd*var_int(i_p4(4,k),j_p4(4,k))

            var_ref4(k) =  (1.0d0 - yd)*var_0 + yd*var_1

            xs = x_ref1(k) - x_c
            ys = y_ref1(k) - y_c
            x2 = x_ref2(k) - x_c
            y2 = y_ref2(k) - y_c
            x3 = x_ref3(k) - x_c
            y3 = y_ref3(k) - y_c
            x4 = x_ref4(k) - x_c
            y4 = y_ref4(k) - y_c

            var2 = var_ref2(k)
            var3 = var_ref3(k)
            var4 = var_ref4(k)
            call quadratic_interpolation(x2,y2,var2,x3,y3,var3,x4,y4,var4,xs,ys,vars)
            int_var_cy(k) = vars			
        end do

    end subroutine extrapolate_var_on_cylinder

    subroutine trapezoidal_sub(n_int,x_int, fx_int, integral)
        implicit none
        integer, intent(in) :: n_int
        real(8), intent(in) :: x_int(n_int), fx_int(n_int)
        real(8), intent(out) :: integral
        integer :: i
        real :: dx_int

        integral = 0.0d0
        do i = 1, n_int - 1
            dx_int = x_int(i+1) - x_int(i)
            integral = integral + 0.50d0 * (fx_int(i) + fx_int(i+1)) * dx_int
        end do
    end subroutine trapezoidal_sub

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

    subroutine moments_interpolation_boundary_nodes()
        implicit none
        real(8) :: term1,term2,term3
        real(8) :: xd,yd,mxx_0,myy_0,mxx_1,myy_1
        real(8) :: mxxp1,mxxp2,mxxp3,mxxp4
        real(8) :: myyp1,myyp2,myyp3,myyp4
        integer :: i,j,k

        do k = 1,nbct
            i = bou_i(k)
            j = bou_j(k)

            !First reference node: cylinder point
            mxxp_ref1(k) = 0.0d0
            myyp_ref1(k) = 0.0d0

            !Second reference node: First fluid point
            xd = (x_ref2(k) - x(i_p2(1,k),j_p2(1,k)))/dx
            yd = (y_ref2(k) - y(i_p2(1,k),j_p2(1,k)))/dy

            !Calculation for Mxxp at the fluid nodes for quadratic interpolations
            mxxp1 = mxx(i_p2(1,k),j_p2(1,k))*(cth_glb(i_p2(1,k),j_p2(1,k))**2) &
                & + myy(i_p2(1,k),j_p2(1,k))*(sth_glb(i_p2(1,k),j_p2(1,k))**2) &
                & + mxy(i_p2(1,k),j_p2(1,k))*s2th_glb(i_p2(1,k),j_p2(1,k))
            mxxp2 = mxx(i_p2(2,k),j_p2(2,k))*(cth_glb(i_p2(2,k),j_p2(2,k))**2) &
                & + myy(i_p2(2,k),j_p2(2,k))*(sth_glb(i_p2(2,k),j_p2(2,k))**2) &
                & + mxy(i_p2(2,k),j_p2(2,k))*s2th_glb(i_p2(2,k),j_p2(2,k))
            mxxp3 = mxx(i_p2(3,k),j_p2(3,k))*(cth_glb(i_p2(3,k),j_p2(3,k))**2) &
                & + myy(i_p2(3,k),j_p2(3,k))*(sth_glb(i_p2(3,k),j_p2(3,k))**2) &
                & + mxy(i_p2(3,k),j_p2(3,k))*s2th_glb(i_p2(3,k),j_p2(3,k))
            mxxp4 = mxx(i_p2(4,k),j_p2(4,k))*(cth_glb(i_p2(4,k),j_p2(4,k))**2) &
                & + myy(i_p2(4,k),j_p2(4,k))*(sth_glb(i_p2(4,k),j_p2(4,k))**2) &
                & + mxy(i_p2(4,k),j_p2(4,k))*s2th_glb(i_p2(4,k),j_p2(4,k))

            !Calculation for Myyp at the fluid nodes for quadratic interpolations						
            myyp1 = mxx(i_p2(1,k),j_p2(1,k))*(sth_glb(i_p2(1,k),j_p2(1,k))**2) &
                & + myy(i_p2(1,k),j_p2(1,k))*(cth_glb(i_p2(1,k),j_p2(1,k))**2) &
                & - mxy(i_p2(1,k),j_p2(1,k))*s2th_glb(i_p2(1,k),j_p2(1,k))
            myyp2 = mxx(i_p2(2,k),j_p2(2,k))*(sth_glb(i_p2(2,k),j_p2(2,k))**2) &
                & + myy(i_p2(2,k),j_p2(2,k))*(cth_glb(i_p2(2,k),j_p2(2,k))**2) &
                & - mxy(i_p2(2,k),j_p2(2,k))*s2th_glb(i_p2(2,k),j_p2(2,k))
            myyp3 = mxx(i_p2(3,k),j_p2(3,k))*(sth_glb(i_p2(3,k),j_p2(3,k))**2) &
                & + myy(i_p2(3,k),j_p2(3,k))*(cth_glb(i_p2(3,k),j_p2(3,k))**2) &
                & - mxy(i_p2(3,k),j_p2(3,k))*s2th_glb(i_p2(3,k),j_p2(3,k))
            myyp4 = mxx(i_p2(4,k),j_p2(4,k))*(sth_glb(i_p2(4,k),j_p2(4,k))**2) &
                & + myy(i_p2(4,k),j_p2(4,k))*(cth_glb(i_p2(4,k),j_p2(4,k))**2) &
                & - mxy(i_p2(4,k),j_p2(4,k))*s2th_glb(i_p2(4,k),j_p2(4,k))

            mxx_0 = (1.0d0 - xd)*mxxp1 + xd*mxxp2
            myy_0 = (1.0d0 - xd)*myyp1 + xd*myyp2

            mxx_1 = (1.0d0 - xd)*mxxp3 + xd*mxxp4
            myy_1 = (1.0d0 - xd)*myyp3 + xd*myyp4

            mxxp_ref2(k) =  (1.0d0 - yd)*mxx_0 + yd*mxx_1
            myyp_ref2(k) =  (1.0d0 - yd)*myy_0 + yd*myy_1

            !Third reference node: Second fluid point
            xd = (x_ref3(k) - x(i_p3(1,k),j_p3(1,k)))/dx
            yd = (y_ref3(k) - y(i_p3(1,k),j_p3(1,k)))/dy

            !Calculation for Mxxp at the fluid nodes for quadratic interpolations
            mxxp1 = mxx(i_p3(1,k),j_p3(1,k))*(cth_glb(i_p3(1,k),j_p3(1,k))**2) &
                & + myy(i_p3(1,k),j_p3(1,k))*(sth_glb(i_p3(1,k),j_p3(1,k))**2) &
                & + mxy(i_p3(1,k),j_p3(1,k))*s2th_glb(i_p3(1,k),j_p3(1,k))
            mxxp2 = mxx(i_p3(2,k),j_p3(2,k))*(cth_glb(i_p3(2,k),j_p3(2,k))**2) &
                & + myy(i_p3(2,k),j_p3(2,k))*(sth_glb(i_p3(2,k),j_p3(2,k))**2) &
                & + mxy(i_p3(2,k),j_p3(2,k))*s2th_glb(i_p3(2,k),j_p3(2,k))
            mxxp3 = mxx(i_p3(3,k),j_p3(3,k))*(cth_glb(i_p3(3,k),j_p3(3,k))**2) &
                & + myy(i_p3(3,k),j_p3(3,k))*(sth_glb(i_p3(3,k),j_p3(3,k))**2) &
                & + mxy(i_p3(3,k),j_p3(3,k))*s2th_glb(i_p3(3,k),j_p3(3,k))
            mxxp4 = mxx(i_p3(4,k),j_p3(4,k))*(cth_glb(i_p3(4,k),j_p3(4,k))**2) &
                & + myy(i_p3(4,k),j_p3(4,k))*(sth_glb(i_p3(4,k),j_p3(4,k))**2) &
                & + mxy(i_p3(4,k),j_p3(4,k))*s2th_glb(i_p3(4,k),j_p3(4,k))
            !Calculation for Myyp at the fluid nodes for quadratic interpolations						
            myyp1 = mxx(i_p3(1,k),j_p3(1,k))*(sth_glb(i_p3(1,k),j_p3(1,k))**2) &
                & + myy(i_p3(1,k),j_p3(1,k))*(cth_glb(i_p3(1,k),j_p3(1,k))**2) &
                & - mxy(i_p3(1,k),j_p3(1,k))*s2th_glb(i_p3(1,k),j_p3(1,k))
            myyp2 = mxx(i_p3(2,k),j_p3(2,k))*(sth_glb(i_p3(2,k),j_p3(2,k))**2) &
                & + myy(i_p3(2,k),j_p3(2,k))*(cth_glb(i_p3(2,k),j_p3(2,k))**2) &
                & - mxy(i_p3(2,k),j_p3(2,k))*s2th_glb(i_p3(2,k),j_p3(2,k))
            myyp3 = mxx(i_p3(3,k),j_p3(3,k))*(sth_glb(i_p3(3,k),j_p3(3,k))**2) &
                & + myy(i_p3(3,k),j_p3(3,k))*(cth_glb(i_p3(3,k),j_p3(3,k))**2) &
                & - mxy(i_p3(3,k),j_p3(3,k))*s2th_glb(i_p3(3,k),j_p3(3,k))
            myyp4 = mxx(i_p3(4,k),j_p3(4,k))*(sth_glb(i_p3(4,k),j_p3(4,k))**2) &
                & + myy(i_p3(4,k),j_p3(4,k))*(cth_glb(i_p3(4,k),j_p3(4,k))**2) &
                & - mxy(i_p3(4,k),j_p3(4,k))*s2th_glb(i_p3(4,k),j_p3(4,k))

            mxx_0 = (1.0d0 - xd)*mxxp1 + xd*mxxp2
            myy_0 = (1.0d0 - xd)*myyp1 + xd*myyp2

            mxx_1 = (1.0d0 - xd)*mxxp3 + xd*mxxp4
            myy_1 = (1.0d0 - xd)*myyp3 + xd*myyp4

            mxxp_ref3(k) =  (1.0d0 - yd)*mxx_0 + yd*mxx_1
            myyp_ref3(k) =  (1.0d0 - yd)*myy_0 + yd*myy_1

            !Variable interpolation:
            term1 = (2.0d0*delx*delx - (delta_uk(k)**2) + 3.0d0* delta_uk(k)*delx)/(2.0d0*(delx**2))
            term2 =  delta_uk(k)*( delta_uk(k) - 2.0d0*delx)/(delx**2)
            term3 =  delta_uk(k)*( delta_uk(k) - delx)/(2.0d0*(delx**2))
            !Boundary node moments:
            mxxp(i,j) = term1*mxxp_ref1(k) + term2*mxxp_ref2(k) - term3*mxxp_ref3(k)
            myyp(i,j) = term1*myyp_ref1(k) + term2*myyp_ref2(k) - term3*myyp_ref3(k) 

        end do

    end subroutine moments_interpolation_boundary_nodes


    subroutine get_coord()
        implicit none
        integer, dimension(0:q-1, 2) :: neigh
        integer :: flag0,flag1,flag2,flag3
        real, parameter :: epsilon = 0.5          ! Tolerance for boundary
        integer :: i, j, k,l, cnt, i1,j1,i2,j2,i3,j3,i_temp,j_temp,ss,ee
        real(8) :: thet,x_s,y_s,thet1,thet2,thet_mid,lambda,radi
        real(8) :: xmid,ymid,epsi,delta_x,delta_y,rr,unit_nx,unit_ny,temp

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

        do i = 1, nx
            do j = 1, ny
                radi = sqrt( ((x(i,j) - x_c)**2) + ((y(i,j) - y_c)**2))
                cth_glb(i,j) = (x(i,j) - x_c)/radi
                sth_glb(i,j) = (y(i,j) - y_c)/radi
                s2th_glb(i,j) = 2.0d0*sth_glb(i,j)*cth_glb(i,j)
                c2th_glb(i,j) = (cth_glb(i,j)**2) - (sth_glb(i,j)**2)
            end do
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

        r_cyl = max_radii


        !Finding boundary and two fluid nodes for velocity interpolation
        !-----------------------------------------------------------------------------------------------
        delx = sqrt(2.0d0)
        do k = 1, nbct

            rr = sqrt((x(bou_i(k),bou_j(k)) - x_c)**2 + (y(bou_i(k),bou_j(k)) - y_c)**2)

            if (rr > 1.0e-12) then
                unit_nx = (x(bou_i(k),bou_j(k)) - x_c)/rr
                unit_ny = (y(bou_i(k),bou_j(k)) - y_c)/rr
                xb(k) = x_c + (r_cyl*unit_nx)
                yb(k) = y_c + (r_cyl*unit_ny)
            else
                ! Point is exactly at center
                unit_nx = 0.0
                unit_ny = 0.0
                xb(k) = x_c
                yb(k) = y_c
            end if

            cth_cyl(k) = unit_nx
            sth_cyl(k) = unit_ny

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

        !Sorting theta values:
        ! Initialize index array
        do i = 1, nbct
            idx(i) = i
        end do

        ! Sort indices based on normalized theta
        do i = 1, nbct - 1
            do j = i + 1, nbct
                if (theta_c(idx(i)) > theta_c(idx(j))) then
                    temp = idx(i)
                    idx(i) = idx(j)
                    idx(j) = temp
                end if
            end do
        end do

        do i = 1, nbct
            thetas_cy(i) = theta_c(idx(i))
        end do

        !-----------------------------------------------------------------------------------------------

        ss = 37
        ee = 37
        open(unit=10, file='data_geo/normal.dat', status="replace")
            do k = ss, ee
                write(10,*) x(bou_i(k),bou_j(k)),y(bou_i(k),bou_j(k))
                write(10,*) x_ref1(k),y_ref1(k)
                write(10,*) x_ref2(k),y_ref2(k)
                write(10,*) x_ref3(k),y_ref3(k)
                write(10,*) x_ref4(k),y_ref4(k)
            end do
        close(10)

        open(unit=10, file='data_geo/surface2.dat', status="replace")
            do k = 0, nbct-1
                thet = k * (2.0d0 * PI / nbct) 
                x_s = x_c + r_cyl*Cos(thet) 
                y_s = y_c + r_cyl*Sin(thet)
                write(10,*) x_s,y_s
            end do
        close(10)


        open(unit=10, file='data_geo/bound.dat', status="replace")
        do  i = 1, nx
            do  j = 1, ny
                if(isbound(i,j)==1) then
                    write(10,*) x(i,j), y(i,j)
                end if
            end do
        end do	
        close(10)

        open(unit=20, file='data_geo/data_geosolid.dat', status="replace")
        do  i = 1, nx
            do  j = 1, ny
                if(issolid(i,j)==1) then
                    write(20,*) x(i,j), y(i,j)
                end if
            end do
        end do
        close(20)

        open(unit=30, file='data_geo/fluid.dat', status="replace")
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
            !rho_prime_b = (4.0d0*rhoI_b + 3.0d0*rhoI_b*mxxI_b)/(3.0d0 + 3.0d0*ux(xi,yj))
            Mxx_prime_b = (rho_prime_b + 9.0d0*rhoI_b*mxxI_b - 3.0d0*rho_prime_b*ux(xi,yj))/(6.0d0*rho_prime_b)
            Myy_prime_b = 6.0d0*rhoI_b*myyI_b/(5.0d0*rho_prime_b)
            Mxy_prime_b = (6.0d0*rhoI_b*mxyI_b - rho_prime_b*uy(xi,yj))/(3.0d0*rho_prime_b)

            if(mom_collision)then
                rho(xi, yj) = rho_prime_b
                mxx(xi, yj) = Mxx_prime_b
                myy(xi, yj) = Myy_prime_b
                mxy(xi, yj) = Mxy_prime_b
            end if

            if(pop_collision)then
                do k = 0, q-1
                    f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
                                & + (3.0d0*uy(xi,yj)*e(k, 2))  &
                                & + (4.50d0*Mxx_prime_b*Hxx(k)) &
                                & + (4.50d0*Myy_prime_b*Hyy(k)) &
                                & + (9.0d0*Mxy_prime_b*Hxy(k)) )
                end do
            end if


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

            if(mom_collision)then
                rho(xi, yj) = rho_prime_b
                mxx(xi, yj) = Mxx_prime_b
                myy(xi, yj) = Myy_prime_b
                mxy(xi, yj) = Mxy_prime_b
            end if

            if(pop_collision)then
                do k = 0, q-1
                    f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
                                & + (3.0d0*uy(xi,yj)*e(k, 2))  &
                                & + (4.50d0*Mxx_prime_b*Hxx(k)) &
                                & + (4.50d0*Myy_prime_b*Hyy(k)) &
                                & + (9.0d0*Mxy_prime_b*Hxy(k)) )
                end do
            end if

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

            if(mom_collision)then
                rho(xi, yj) = rho_prime_b
                mxx(xi, yj) = Mxx_prime_b
                myy(xi, yj) = Myy_prime_b
                mxy(xi, yj) = Mxy_prime_b
            end if

            if(pop_collision)then
                do k = 0, q-1
                    f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
                                & + (3.0d0*uy(xi,yj)*e(k, 2))  &
                                & + (4.50d0*Mxx_prime_b*Hxx(k)) &
                                & + (4.50d0*Myy_prime_b*Hyy(k)) &
                                & + (9.0d0*Mxy_prime_b*Hxy(k)) )
                end do
            end if

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

            if(mom_collision)then
                rho(xi, yj) = rho_prime_b
                mxx(xi, yj) = Mxx_prime_b
                myy(xi, yj) = Myy_prime_b
                mxy(xi, yj) = Mxy_prime_b
            end if

            if(pop_collision)then
                do k = 0, q-1
                    f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
                                & + (3.0d0*uy(xi,yj)*e(k, 2))  &
                                & + (4.50d0*Mxx_prime_b*Hxx(k)) &
                                & + (4.50d0*Myy_prime_b*Hyy(k)) &
                                & + (9.0d0*Mxy_prime_b*Hxy(k)) )
                end do
            end if

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

            !rho_prime_b = 6.0d0*rhoI_b/(5.0d0 + 3.0d0*ux(xi,yj)- 3.0d0*(ux(xi,yj)**2))
            rho_prime_b = rho(xi,yj)
            Mxx_prime_b = ux(xi,yj)**2
            Myy_prime_b = uy(xi,yj)**2
            Mxy_prime_b = (6.0d0*rhoI_b*mxyI_b - rho_prime_b*uy(xi,yj))/(3.0d0*rho_prime_b)

            if(mom_collision)then
                rho(xi, yj) = rho_prime_b
                mxx(xi, yj) = Mxx_prime_b
                myy(xi, yj) = Myy_prime_b
                mxy(xi, yj) = Mxy_prime_b
            end if

            if(pop_collision)then
                do k = 0, q-1
                    f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
                                & + (3.0d0*uy(xi,yj)*e(k, 2))  &
                                & + (4.50d0*Mxx_prime_b*Hxx(k)) &
                                & + (4.50d0*Myy_prime_b*Hyy(k)) &
                                & + (9.0d0*Mxy_prime_b*Hxy(k)) )
                end do
            end if
        
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

            if(mom_collision)then
                rho(xi, yj) = rho_prime_b
                mxx(xi, yj) = Mxx_prime_b
                myy(xi, yj) = Myy_prime_b
                mxy(xi, yj) = Mxy_prime_b
            end if
            
            if(pop_collision)then
                do k = 0, q-1
                    f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
                                & + (3.0d0*uy(xi,yj)*e(k, 2))  &
                                & + (4.50d0*Mxx_prime_b*Hxx(k)) &
                                & + (4.50d0*Myy_prime_b*Hyy(k)) &
                                & + (9.0d0*Mxy_prime_b*Hxy(k)) )
                end do
            end if

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

            if(mom_collision)then
                rho(xi, yj) = rho_prime_b
                mxx(xi, yj) = Mxx_prime_b
                myy(xi, yj) = Myy_prime_b
                mxy(xi, yj) = Mxy_prime_b
            end if

            if(pop_collision)then
                do k = 0, q-1
                    f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
                                & + (3.0d0*uy(xi,yj)*e(k, 2))  &
                                & + (4.50d0*Mxx_prime_b*Hxx(k)) &
                                & + (4.50d0*Myy_prime_b*Hyy(k)) &
                                & + (9.0d0*Mxy_prime_b*Hxy(k)) )
                end do
            end if

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

            if(mom_collision)then
                rho(xi, yj) = rho_prime_b
                mxx(xi, yj) = Mxx_prime_b
                myy(xi, yj) = Myy_prime_b
                mxy(xi, yj) = Mxy_prime_b
            end if

            if(pop_collision)then
                do k = 0, q-1
                    f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*ux(xi,yj)*e(k, 1)) &
                                & + (3.0d0*uy(xi,yj)*e(k, 2))  &
                                & + (4.50d0*Mxx_prime_b*Hxx(k)) &
                                & + (4.50d0*Myy_prime_b*Hyy(k)) &
                                & + (9.0d0*Mxy_prime_b*Hxy(k)) )
                end do
            end if

        CASE DEFAULT
            print*,"Not a valid boundary type!!!!!!!"
        END SELECT


    end subroutine apply_bc_incomp

    subroutine numerical_boundary_cases(xi,yj,label,uxp,uyp)
        implicit none
        integer,intent(in) :: xi,yj,label
        real(8),intent(in) :: uxp,uyp
        integer,dimension(0:q-1) :: Is,Os
        real(8),dimension(0:q-1) :: A_i,E_i, B11_i, B22_i, B12_i
        real(8),dimension(1:3,1:3) :: A_coeff
        real(8),dimension(1:3) :: b_coeff
        real(8) :: rhoI_b,mxxI_b,myyI_b,mxyI_b
        real(8) :: rho_prime_b,Mxx_prime_b,Myy_prime_b,Mxy_prime_b
        real(8) :: A_prime, E_prime, G_prime, B11_prime, B22_prime, B12_prime
        real(8) :: D_xx_prime, D_yy_prime, D_xy_prime, J11_prime, J22_prime, J12_prime
        real(8) :: F11_xx_prime, F12_xx_prime, F22_xx_prime, F11_yy_prime, F12_yy_prime, F22_yy_prime
        real(8) :: F11_xy_prime, F12_xy_prime, F22_xy_prime
        real(8) :: J11_xx_star, J22_xx_star, J12_xx_star, J11_yy_star, J22_yy_star, J12_yy_star,&
                    & J11_xy_star, J22_xy_star, J12_xy_star
        real(8) :: L_11_xx, L_22_xx, L_12_xx, L_11_yy, L_22_yy, L_12_yy, L_11_xy, L_22_xy, L_12_xy
        real(8) :: R_xx, R_yy, R_xy, denominator
        integer :: nvar_sys = 3


        Is(0:q-1) = Incs_b(label,0:q-1)
        Os(0:q-1) = Outs_b(label,0:q-1)

        do k = 0, q-1
            A_i(k) = w(k)*( 1.0d0 + 3.0d0*uxp*e(k, 1) + 3.0d0*uyp*e(k, 2) )
            E_i(k) = w(k)*( 1.0d0 + 3.0d0*uxp*e(k, 1) + 3.0d0*uyp*e(k, 2) &
                & + 4.50d0*(uxp**2)*Hxx(k) + 4.50d0*(uyp**2)*Hyy(k) + 9.0d0*(uxp*uyp)*Hxy(k) )
            B11_i(k) = 4.50d0*w(k)*Hxx(k)
            B22_i(k) = 4.50d0*w(k)*Hyy(k)
            B12_i(k) = 4.50d0*w(k)*Hxy(k)
        end do

        !gamma,delta = 1,2,3
        !alpha,beta = x,y,z

        rhoI_b = 0.0d0; mxxI_b = 0.0d0; myyI_b = 0.0d0; mxyI_b = 0.0d0
        A_prime = 0.0d0; E_prime = 0.0d0

        B11_prime = 0.0d0; B22_prime = 0.0d0; B12_prime = 0.0d0
        D_xx_prime = 0.0d0; D_yy_prime = 0.0d0; D_xy_prime = 0.0d0 
        F11_xx_prime = 0.0d0; F12_xx_prime = 0.0d0; F22_xx_prime = 0.0d0
        F11_yy_prime = 0.0d0; F12_yy_prime = 0.0d0; F22_yy_prime = 0.0d0
        F11_xy_prime = 0.0d0; F12_xy_prime = 0.0d0; F22_xy_prime = 0.0d0

        do k = 0, q-1
            if(Os(k)==1) then
                A_prime = A_prime + A_i(k)
                E_prime = E_prime + E_i(k)

                B11_prime = B11_prime + B11_i(k)
                B22_prime = B22_prime + B22_i(k)
                B12_prime = B12_prime + B12_i(k)

            end if

            if(Is(k)==1) then
                rhoI_b = rhoI_b + f(xi,yj,k)
                mxxI_b = mxxI_b + f(xi,yj,k)*Hxx(k)
                myyI_b = myyI_b + f(xi,yj,k)*Hyy(k)
                mxyI_b = mxyI_b + f(xi,yj,k)*Hxy(k)

                D_xx_prime = D_xx_prime + A_i(k)*Hxx(k)
                D_yy_prime = D_yy_prime + A_i(k)*Hyy(k)
                D_xy_prime = D_xy_prime + A_i(k)*Hxy(k)

                F11_xx_prime = F11_xx_prime + B11_i(k)*Hxx(k)
                F22_xx_prime = F22_xx_prime + B22_i(k)*Hxx(k)
                F12_xx_prime = F12_xx_prime + B12_i(k)*Hxx(k)

                F11_yy_prime = F11_yy_prime + B11_i(k)*Hyy(k)
                F22_yy_prime = F22_yy_prime + B22_i(k)*Hyy(k)
                F12_yy_prime = F12_yy_prime + B12_i(k)*Hyy(k)

                F11_xy_prime = F11_xy_prime + B11_i(k)*Hxy(k)
                F22_xy_prime = F22_xy_prime + B22_i(k)*Hxy(k)
                F12_xy_prime = F12_xy_prime + B12_i(k)*Hxy(k)

            end if
        end do

        mxxI_b = mxxI_b/rhoI_b
        myyI_b = myyI_b/rhoI_b
        mxyI_b = mxyI_b/rhoI_b

        G_prime = (1.0d0 - omega)*A_prime + omega*E_prime

        J11_prime = (1.0d0 - omega)*B11_prime
        J22_prime = (1.0d0 - omega)*B22_prime
        J12_prime = (1.0d0 - omega)*B12_prime

        J11_xx_star = J11_prime * mxxI_b
        J22_xx_star = J22_prime * mxxI_b
        J12_xx_star = J12_prime * mxxI_b

        J11_yy_star = J11_prime * myyI_b
        J22_yy_star = J22_prime * myyI_b
        J12_yy_star = J12_prime * myyI_b

        J11_xy_star = J11_prime * mxyI_b
        J22_xy_star = J22_prime * mxyI_b
        J12_xy_star = J12_prime * mxyI_b

        L_11_xx = J11_xx_star - F11_xx_prime
        L_11_yy = J11_yy_star - F11_yy_prime
        L_11_xy = J11_xy_star - F11_xy_prime

        L_22_xx = J22_xx_star - F22_xx_prime
        L_22_yy = J22_yy_star - F22_yy_prime
        L_22_xy = J22_xy_star - F22_xy_prime

        L_12_xx = 2.0d0 * (J12_xx_star - F12_xx_prime)
        L_12_yy = 2.0d0 * (J12_yy_star - F12_yy_prime)
        L_12_xy = 2.0d0 * (J12_xy_star - F12_xy_prime)

        R_xx = D_xx_prime -  G_prime * mxxI_b
        R_yy = D_yy_prime -  G_prime * myyI_b
        R_xy = D_xy_prime -  G_prime * mxyI_b

        A_coeff(1,1:3) = (/ L_11_xx, L_22_xx, L_12_xx /) 
        A_coeff(2,1:3) = (/ L_11_yy, L_22_yy, L_12_yy /)
        A_coeff(3,1:3) = (/ L_11_xy, L_22_xy, L_12_xy /)  

        b_coeff(1:3) = (/ R_xx, R_yy, R_xy /) 


        call solve_system(nvar_sys, A_coeff,b_coeff,Mxx_prime_b,Myy_prime_b,Mxy_prime_b)

        denominator = (1.0d0 - omega)*A_prime + (1.0d0 - omega)*(Mxx_prime_b*B11_prime + Myy_prime_b*B22_prime &
                        & + 2.0d0*Mxy_prime_b*B12_prime) + omega*E_prime
        rho_prime_b = rhoI_b/denominator

        if(mom_collision)then
            rho(xi, yj) = rho_prime_b
            mxx(xi, yj) = Mxx_prime_b
            myy(xi, yj) = Myy_prime_b
            mxy(xi, yj) = Mxy_prime_b
        end if
        
        if(pop_collision)then
            do k = 0, q-1
                f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*e(k, 1)) &
                            & + (3.0d0*uyp*e(k, 2))  &
                            & + (4.50d0*Mxx_prime_b*Hxx(k)) &
                            & + (4.50d0*Myy_prime_b*Hyy(k)) &
                            & + (9.0d0*Mxy_prime_b*Hxy(k)) )
            end do
        end if


    end subroutine numerical_boundary_cases

    subroutine numerical_boundary_cases_rotation(xi,yj,label,uxp,uyp)
        implicit none
        integer,intent(in) :: xi,yj,label
        real(8),intent(in) :: uxp,uyp
        integer,dimension(0:q-1) :: Is,Os
        real(8),dimension(0:q-1) :: A_i,E_i, B11_i, B22_i, B12_i
        real(8),dimension(1:3,1:3) :: A_coeff
        real(8),dimension(1:3) :: b_coeff
        real(8) :: rhoI_b,mxxI_b,myyI_b,mxyI_b
        real(8) :: rho_prime_b,Mxx_prime_b,Myy_prime_b,Mxy_prime_b
        real(8) :: A_prime, E_prime, G_prime, B11_prime, B22_prime, B12_prime
        real(8) :: D_xx_prime, D_yy_prime, D_xy_prime, J11_prime, J22_prime, J12_prime
        real(8) :: F11_xx_prime, F12_xx_prime, F22_xx_prime, F11_yy_prime, F12_yy_prime, F22_yy_prime
        real(8) :: F11_xy_prime, F12_xy_prime, F22_xy_prime
        real(8) :: J11_xx_star, J22_xx_star, J12_xx_star, J11_yy_star, J22_yy_star, J12_yy_star, &
                    & J11_xy_star, J22_xy_star, J12_xy_star
        real(8) :: L_11_xx, L_22_xx, L_12_xx, L_11_yy, L_22_yy, L_12_yy, L_11_xy, L_22_xy, L_12_xy
        real(8) :: R_xx, R_yy, R_xy, denominator
        integer :: nvar_sys = 3


        Is(0:q-1) = Incs_b(label,0:q-1)
        Os(0:q-1) = Outs_b(label,0:q-1)

        do k = 0, q-1
            A_i(k) = w(k)*( 1.0d0 + 3.0d0*uxp*cx_p(xi, yj, k) + 3.0d0*uyp*cy_p(xi, yj, k) )
            E_i(k) = w(k)*( 1.0d0 + 3.0d0*uxp*cx_p(xi, yj, k) + 3.0d0*uyp*cy_p(xi, yj, k) &
                & + 4.50d0*(uxp**2)*Hxx_p(xi, yj, k) + 4.50d0*(uyp**2)*Hyy_p(xi, yj, k) + 9.0d0*(uxp*uyp)*Hxy_p(xi, yj, k) )
            B11_i(k) = 4.50d0*w(k)*Hxx_p(xi, yj, k)
            B22_i(k) = 4.50d0*w(k)*Hyy_p(xi, yj, k)
            B12_i(k) = 4.50d0*w(k)*Hxy_p(xi, yj, k)
        end do

        !gamma,delta = 1,2,3
        !alpha,beta = x,y,z

        rhoI_b = 0.0d0; mxxI_b = 0.0d0; myyI_b = 0.0d0; mxyI_b = 0.0d0
        A_prime = 0.0d0; E_prime = 0.0d0

        B11_prime = 0.0d0; B22_prime = 0.0d0; B12_prime = 0.0d0
        D_xx_prime = 0.0d0; D_yy_prime = 0.0d0; D_xy_prime = 0.0d0 
        F11_xx_prime = 0.0d0; F12_xx_prime = 0.0d0; F22_xx_prime = 0.0d0
        F11_yy_prime = 0.0d0; F12_yy_prime = 0.0d0; F22_yy_prime = 0.0d0
        F11_xy_prime = 0.0d0; F12_xy_prime = 0.0d0; F22_xy_prime = 0.0d0

        do k = 0, q-1
            if(Os(k)==1) then
                A_prime = A_prime + A_i(k)
                E_prime = E_prime + E_i(k)

                B11_prime = B11_prime + B11_i(k)
                B22_prime = B22_prime + B22_i(k)
                B12_prime = B12_prime + B12_i(k)

            end if

            if(Is(k)==1) then
                rhoI_b = rhoI_b + f(xi, yj, k)
                mxxI_b = mxxI_b + f(xi, yj, k)*Hxx_p(xi, yj, k)
                myyI_b = myyI_b + f(xi, yj, k)*Hyy_p(xi, yj, k)
                mxyI_b = mxyI_b + f(xi, yj, k)*Hxy_p(xi, yj, k)

                D_xx_prime = D_xx_prime + A_i(k)*Hxx_p(xi, yj, k)
                D_yy_prime = D_yy_prime + A_i(k)*Hyy_p(xi, yj, k)
                D_xy_prime = D_xy_prime + A_i(k)*Hxy_p(xi, yj, k)

                F11_xx_prime = F11_xx_prime + B11_i(k)*Hxx_p(xi, yj, k)
                F22_xx_prime = F22_xx_prime + B22_i(k)*Hxx_p(xi, yj, k)
                F12_xx_prime = F12_xx_prime + B12_i(k)*Hxx_p(xi, yj, k)

                F11_yy_prime = F11_yy_prime + B11_i(k)*Hyy_p(xi, yj, k)
                F22_yy_prime = F22_yy_prime + B22_i(k)*Hyy_p(xi, yj, k)
                F12_yy_prime = F12_yy_prime + B12_i(k)*Hyy_p(xi, yj, k)

                F11_xy_prime = F11_xy_prime + B11_i(k)*Hxy_p(xi, yj, k)
                F22_xy_prime = F22_xy_prime + B22_i(k)*Hxy_p(xi, yj, k)
                F12_xy_prime = F12_xy_prime + B12_i(k)*Hxy_p(xi, yj, k)

            end if
        end do

        mxxI_b = mxxI_b/rhoI_b
        myyI_b = myyI_b/rhoI_b
        mxyI_b = mxyI_b/rhoI_b



        G_prime = (1.0d0 - omega)*A_prime + omega*E_prime

        J11_prime = (1.0d0 - omega)*B11_prime
        J22_prime = (1.0d0 - omega)*B22_prime
        J12_prime = (1.0d0 - omega)*B12_prime

        J11_xx_star = J11_prime * mxxI_b
        J22_xx_star = J22_prime * mxxI_b
        J12_xx_star = J12_prime * mxxI_b

        J11_yy_star = J11_prime * myyI_b
        J22_yy_star = J22_prime * myyI_b
        J12_yy_star = J12_prime * myyI_b

        J11_xy_star = J11_prime * mxyI_b
        J22_xy_star = J22_prime * mxyI_b
        J12_xy_star = J12_prime * mxyI_b

        L_11_xx = J11_xx_star - F11_xx_prime
        L_11_yy = J11_yy_star - F11_yy_prime
        L_11_xy = J11_xy_star - F11_xy_prime

        L_22_xx = J22_xx_star - F22_xx_prime
        L_22_yy = J22_yy_star - F22_yy_prime
        L_22_xy = J22_xy_star - F22_xy_prime

        L_12_xx = 2.0d0 * (J12_xx_star - F12_xx_prime)
        L_12_yy = 2.0d0 * (J12_yy_star - F12_yy_prime)
        L_12_xy = 2.0d0 * (J12_xy_star - F12_xy_prime)

        R_xx = D_xx_prime -  G_prime * mxxI_b
        R_yy = D_yy_prime -  G_prime * myyI_b
        R_xy = D_xy_prime -  G_prime * mxyI_b

        Mxx_prime_b = mxxp(xi, yj)
        Myy_prime_b = myyp(xi, yj)
        Mxy_prime_b =  (R_xy - Mxx_prime_b * L_11_xy - Myy_prime_b * L_22_xy) / L_12_xy

        denominator = (1.0d0 - omega)*A_prime + (1.0d0 - omega)*(Mxx_prime_b*B11_prime + Myy_prime_b*B22_prime &
            & + 2.0d0*Mxy_prime_b*B12_prime) + omega*E_prime
        rho_prime_b = rhoI_b/denominator
        
        if(mom_collision) then
            rho(xi, yj) = rho_prime_b
            mxx(xi, yj) = Mxx_prime_b * (cth_glb(xi,yj)**2) + Myy_prime_b *(sth_glb(xi,yj)**2) - Mxy_prime_b * s2th_glb(xi,yj)
            myy(xi, yj) = Mxx_prime_b * (sth_glb(xi,yj)**2) + Myy_prime_b *(cth_glb(xi,yj)**2) + Mxy_prime_b * s2th_glb(xi,yj)
            mxy(xi, yj) = (Mxx_prime_b - Myy_prime_b) * 0.50d0 * s2th_glb(xi,yj) + Mxy_prime_b * c2th_glb(xi,yj)
        end if

        if(pop_collision) then
            do k = 0, q-1
                f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*cx_p(xi, yj, k)) &
                            & + (3.0d0*uyp*cy_p(xi, yj, k))  &
                            & + (4.50d0*Mxx_prime_b*Hxx_p(xi, yj, k)) &
                            & + (4.50d0*Myy_prime_b*Hyy_p(xi, yj, k)) &
                            & + (9.0d0*Mxy_prime_b*Hxy_p(xi, yj, k)) )
            end do
        end if


    end subroutine numerical_boundary_cases_rotation
    
    subroutine numerical_boundary_cases_rotation_weak(xi,yj,label,uxp,uyp)
        implicit none
        integer,intent(in) :: xi,yj,label
        real(8),intent(in) :: uxp,uyp
        integer,dimension(0:q-1) :: Is,Os
        real(8),dimension(0:q-1) :: A_i,E_i, B11_i, B22_i, B12_i
        real(8),dimension(1:3,1:3) :: A_coeff
        real(8),dimension(1:3) :: b_coeff
        real(8) :: rhoI_b,mxxI_b,myyI_b,mxyI_b
        real(8) :: rho_prime_b,Mxx_prime_b,Myy_prime_b,Mxy_prime_b
        real(8) :: A_prime, E_prime, G_prime, B11_prime, B22_prime, B12_prime
        real(8) :: D_xx_prime, D_yy_prime, D_xy_prime, J11_prime, J22_prime, J12_prime
        real(8) :: F11_xx_prime, F12_xx_prime, F22_xx_prime, F11_yy_prime, F12_yy_prime, F22_yy_prime
        real(8) :: F11_xy_prime, F12_xy_prime, F22_xy_prime
        real(8) :: J11_xx_star, J22_xx_star, J12_xx_star, J11_yy_star, J22_yy_star, J12_yy_star, &
                    & J11_xy_star, J22_xy_star, J12_xy_star
        real(8) :: L_11_xx, L_22_xx, L_12_xx, L_11_yy, L_22_yy, L_12_yy, L_11_xy, L_22_xy, L_12_xy
        real(8) :: R_xx, R_yy, R_xy, denominator
        integer :: nvar_sys = 3


        Is(0:q-1) = Incs_b(label,0:q-1)
        Os(0:q-1) = Outs_b(label,0:q-1)

        do k = 0, q-1
            A_i(k) = w(k)*( 1.0d0 + 3.0d0*uxp*cx_p(xi, yj, k) + 3.0d0*uyp*cy_p(xi, yj, k) )
            B11_i(k) = 4.50d0*w(k)*Hxx_p(xi, yj, k)
            B22_i(k) = 4.50d0*w(k)*Hyy_p(xi, yj, k)
            B12_i(k) = 4.50d0*w(k)*Hxy_p(xi, yj, k)
        end do

        !gamma,delta = 1,2,3
        !alpha,beta = x,y,z

        rhoI_b = 0.0d0; mxxI_b = 0.0d0; myyI_b = 0.0d0; mxyI_b = 0.0d0
        A_prime = 0.0d0; E_prime = 0.0d0

        B11_prime = 0.0d0; B22_prime = 0.0d0; B12_prime = 0.0d0
        D_xx_prime = 0.0d0; D_yy_prime = 0.0d0; D_xy_prime = 0.0d0 
        F11_xx_prime = 0.0d0; F12_xx_prime = 0.0d0; F22_xx_prime = 0.0d0
        F11_yy_prime = 0.0d0; F12_yy_prime = 0.0d0; F22_yy_prime = 0.0d0
        F11_xy_prime = 0.0d0; F12_xy_prime = 0.0d0; F22_xy_prime = 0.0d0

        do k = 0, q-1

            if(Is(k)==1) then
                rhoI_b = rhoI_b + f(xi, yj, k)
                mxxI_b = mxxI_b + f(xi, yj, k)*Hxx_p(xi, yj, k)
                myyI_b = myyI_b + f(xi, yj, k)*Hyy_p(xi, yj, k)
                mxyI_b = mxyI_b + f(xi, yj, k)*Hxy_p(xi, yj, k)
                
                A_prime = A_prime + A_i(k)

                B11_prime = B11_prime + B11_i(k)
                B22_prime = B22_prime + B22_i(k)
                B12_prime = B12_prime + B12_i(k)

                D_xx_prime = D_xx_prime + A_i(k)*Hxx_p(xi, yj, k)
                D_yy_prime = D_yy_prime + A_i(k)*Hyy_p(xi, yj, k)
                D_xy_prime = D_xy_prime + A_i(k)*Hxy_p(xi, yj, k)

                F11_xx_prime = F11_xx_prime + B11_i(k)*Hxx_p(xi, yj, k)
                F22_xx_prime = F22_xx_prime + B22_i(k)*Hxx_p(xi, yj, k)
                F12_xx_prime = F12_xx_prime + B12_i(k)*Hxx_p(xi, yj, k)

                F11_yy_prime = F11_yy_prime + B11_i(k)*Hyy_p(xi, yj, k)
                F22_yy_prime = F22_yy_prime + B22_i(k)*Hyy_p(xi, yj, k)
                F12_yy_prime = F12_yy_prime + B12_i(k)*Hyy_p(xi, yj, k)

                F11_xy_prime = F11_xy_prime + B11_i(k)*Hxy_p(xi, yj, k)
                F22_xy_prime = F22_xy_prime + B22_i(k)*Hxy_p(xi, yj, k)
                F12_xy_prime = F12_xy_prime + B12_i(k)*Hxy_p(xi, yj, k)

            end if
        end do

        mxxI_b = mxxI_b/rhoI_b
        myyI_b = myyI_b/rhoI_b
        mxyI_b = mxyI_b/rhoI_b

        L_11_xx = B11_prime * mxxI_b - F11_xx_prime
        L_11_yy = B11_prime * myyI_b - F11_yy_prime
        L_11_xy = B11_prime * mxyI_b - F11_xy_prime

        L_22_xx = B22_prime * mxxI_b- F22_xx_prime
        L_22_yy = B22_prime * myyI_b- F22_yy_prime
        L_22_xy = B22_prime * mxyI_b- F22_xy_prime

        L_12_xx = 2.0d0 * (B12_prime * mxxI_b- F12_xx_prime)
        L_12_yy = 2.0d0 * (B12_prime * myyI_b- F12_yy_prime)
        L_12_xy = 2.0d0 * (B12_prime * mxyI_b- F12_xy_prime)

        R_xx = D_xx_prime -  A_prime * mxxI_b
        R_yy = D_yy_prime -  A_prime * myyI_b
        R_xy = D_xy_prime -  A_prime * mxyI_b

        Mxx_prime_b = mxxp(xi, yj)
        Myy_prime_b = myyp(xi, yj)
        Mxy_prime_b =  (R_xy - Mxx_prime_b * L_11_xy - Myy_prime_b * L_22_xy) / L_12_xy

        denominator = A_prime + (Mxx_prime_b*B11_prime + Myy_prime_b*B22_prime  + 2.0d0*Mxy_prime_b*B12_prime)
        rho_prime_b = rhoI_b/denominator
        
        if(mom_collision) then
            rho(xi, yj) = rho_prime_b
            mxx(xi, yj) = Mxx_prime_b * (cth_glb(xi,yj)**2) + Myy_prime_b *(sth_glb(xi,yj)**2) - Mxy_prime_b * s2th_glb(xi,yj)
            myy(xi, yj) = Mxx_prime_b * (sth_glb(xi,yj)**2) + Myy_prime_b *(cth_glb(xi,yj)**2) + Mxy_prime_b * s2th_glb(xi,yj)
            mxy(xi, yj) = (Mxx_prime_b - Myy_prime_b) * 0.50d0 * s2th_glb(xi,yj) + Mxy_prime_b * c2th_glb(xi,yj)
        end if

        if(pop_collision) then
            do k = 0, q-1
                f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*cx_p(xi, yj, k)) &
                            & + (3.0d0*uyp*cy_p(xi, yj, k))  &
                            & + (4.50d0*Mxx_prime_b*Hxx_p(xi, yj, k)) &
                            & + (4.50d0*Myy_prime_b*Hyy_p(xi, yj, k)) &
                            & + (9.0d0*Mxy_prime_b*Hxy_p(xi, yj, k)) )
            end do
        end if


    end subroutine numerical_boundary_cases_rotation_weak
    
    subroutine numerical_boundary_cases_rotation_rhoeq(xi,yj,label,uxp,uyp)
        implicit none
        integer,intent(in) :: xi,yj,label
        real(8),intent(in) :: uxp,uyp
        integer,dimension(0:q-1) :: Is,Os
        real(8),dimension(0:q-1) :: A_i,E_i, B11_i, B22_i, B12_i
        real(8),dimension(1:3,1:3) :: A_coeff
        real(8),dimension(1:3) :: b_coeff
        real(8) :: rhoI_b,mxxI_b,myyI_b,mxyI_b
        real(8) :: rho_prime_b,Mxx_prime_b,Myy_prime_b,Mxy_prime_b
        real(8) :: A_prime, E_prime, G_prime, B11_prime, B22_prime, B12_prime
        real(8) :: D_xx_prime, D_yy_prime, D_xy_prime, J11_prime, J22_prime, J12_prime
        real(8) :: F11_xx_prime, F12_xx_prime, F22_xx_prime, F11_yy_prime, F12_yy_prime, F22_yy_prime
        real(8) :: F11_xy_prime, F12_xy_prime, F22_xy_prime
        real(8) :: J11_xx_star, J22_xx_star, J12_xx_star, J11_yy_star, J22_yy_star, J12_yy_star, &
                    & J11_xy_star, J22_xy_star, J12_xy_star
        real(8) :: L_11_xx, L_22_xx, L_12_xx, L_11_yy, L_22_yy, L_12_yy, L_11_xy, L_22_xy, L_12_xy
        real(8) :: R_xx, R_yy, R_xy, denominator
        integer :: nvar_sys = 3


        Is(0:q-1) = Incs_b(label,0:q-1)
        Os(0:q-1) = Outs_b(label,0:q-1)

        do k = 0, q-1
            A_i(k) = w(k)*( 1.0d0 + 3.0d0*uxp*cx_p(xi, yj, k) + 3.0d0*uyp*cy_p(xi, yj, k) )
            E_i(k) = w(k)*( 1.0d0 + 3.0d0*uxp*cx_p(xi, yj, k) + 3.0d0*uyp*cy_p(xi, yj, k) &
                & + 4.50d0*(uxp**2)*Hxx_p(xi, yj, k) + 4.50d0*(uyp**2)*Hyy_p(xi, yj, k) + 9.0d0*(uxp*uyp)*Hxy_p(xi, yj, k) )
            B11_i(k) = 4.50d0*w(k)*Hxx_p(xi, yj, k)
            B22_i(k) = 4.50d0*w(k)*Hyy_p(xi, yj, k)
            B12_i(k) = 4.50d0*w(k)*Hxy_p(xi, yj, k)
        end do

        !gamma,delta = 1,2,3
        !alpha,beta = x,y,z

        rhoI_b = 0.0d0; mxxI_b = 0.0d0; myyI_b = 0.0d0; mxyI_b = 0.0d0
        E_prime = 0.0d0

        B11_prime = 0.0d0; B22_prime = 0.0d0; B12_prime = 0.0d0
        D_xx_prime = 0.0d0; D_yy_prime = 0.0d0; D_xy_prime = 0.0d0 
        F11_xx_prime = 0.0d0; F12_xx_prime = 0.0d0; F22_xx_prime = 0.0d0
        F11_yy_prime = 0.0d0; F12_yy_prime = 0.0d0; F22_yy_prime = 0.0d0
        F11_xy_prime = 0.0d0; F12_xy_prime = 0.0d0; F22_xy_prime = 0.0d0

        do k = 0, q-1
            if(Is(k)==1) then
                rhoI_b = rhoI_b + f(xi, yj, k)
                mxxI_b = mxxI_b + f(xi, yj, k)*Hxx_p(xi, yj, k)
                myyI_b = myyI_b + f(xi, yj, k)*Hyy_p(xi, yj, k)
                mxyI_b = mxyI_b + f(xi, yj, k)*Hxy_p(xi, yj, k)

                E_prime = E_prime + E_i(k)
                
                D_xx_prime = D_xx_prime + A_i(k)*Hxx_p(xi, yj, k)
                D_yy_prime = D_yy_prime + A_i(k)*Hyy_p(xi, yj, k)
                D_xy_prime = D_xy_prime + A_i(k)*Hxy_p(xi, yj, k)

                F11_xx_prime = F11_xx_prime + B11_i(k)*Hxx_p(xi, yj, k)
                F22_xx_prime = F22_xx_prime + B22_i(k)*Hxx_p(xi, yj, k)
                F12_xx_prime = F12_xx_prime + B12_i(k)*Hxx_p(xi, yj, k)

                F11_yy_prime = F11_yy_prime + B11_i(k)*Hyy_p(xi, yj, k)
                F22_yy_prime = F22_yy_prime + B22_i(k)*Hyy_p(xi, yj, k)
                F12_yy_prime = F12_yy_prime + B12_i(k)*Hyy_p(xi, yj, k)

                F11_xy_prime = F11_xy_prime + B11_i(k)*Hxy_p(xi, yj, k)
                F22_xy_prime = F22_xy_prime + B22_i(k)*Hxy_p(xi, yj, k)
                F12_xy_prime = F12_xy_prime + B12_i(k)*Hxy_p(xi, yj, k)

            end if
        end do

        mxxI_b = mxxI_b/rhoI_b
        myyI_b = myyI_b/rhoI_b
        mxyI_b = mxyI_b/rhoI_b

        L_11_xx = F11_xx_prime
        L_11_yy = F11_yy_prime
        L_11_xy = F11_xy_prime

        L_22_xx = F22_xx_prime
        L_22_yy = F22_yy_prime
        L_22_xy = F22_xy_prime

        L_12_xx = 2.0d0 * F12_xx_prime
        L_12_yy = 2.0d0 * F12_yy_prime
        L_12_xy = 2.0d0 * F12_xy_prime

        R_xx = E_prime * mxxI_b - D_xx_prime
        R_yy = E_prime * myyI_b - D_yy_prime
        R_xy = E_prime * mxyI_b - D_xy_prime

        Mxx_prime_b = mxxp(xi, yj)
        Myy_prime_b = myyp(xi, yj)
        Mxy_prime_b =  (R_xy - Mxx_prime_b * L_11_xy - Myy_prime_b * L_22_xy) / L_12_xy

        rho_prime_b = rhoI_b/E_prime
        
        if(mom_collision) then
            rho(xi, yj) = rho_prime_b
            mxx(xi, yj) = Mxx_prime_b * (cth_glb(xi,yj)**2) + Myy_prime_b *(sth_glb(xi,yj)**2) - Mxy_prime_b * s2th_glb(xi,yj)
            myy(xi, yj) = Mxx_prime_b * (sth_glb(xi,yj)**2) + Myy_prime_b *(cth_glb(xi,yj)**2) + Mxy_prime_b * s2th_glb(xi,yj)
            mxy(xi, yj) = (Mxx_prime_b - Myy_prime_b) * 0.50d0 * s2th_glb(xi,yj) + Mxy_prime_b * c2th_glb(xi,yj)
        end if

        if(pop_collision) then
            do k = 0, q-1
                f(xi,yj,k) = w(k)*rho_prime_b*( 1.0d0 + (3.0d0*uxp*cx_p(xi, yj, k)) &
                            & + (3.0d0*uyp*cy_p(xi, yj, k))  &
                            & + (4.50d0*Mxx_prime_b*Hxx_p(xi, yj, k)) &
                            & + (4.50d0*Myy_prime_b*Hyy_p(xi, yj, k)) &
                            & + (9.0d0*Mxy_prime_b*Hxy_p(xi, yj, k)) )
            end do
        end if


    end subroutine numerical_boundary_cases_rotation_rhoeq

    subroutine solve_system(nvar, A_sys,b_sys,x_11,x_22,x_12)
        implicit none
        integer,intent(in) :: nvar
        real(8),dimension(1:nvar,1:nvar),intent(in) :: A_sys 
        real(8),dimension(1:nvar),intent(in) :: b_sys 
        real(8),intent(out) :: x_11,x_22,x_12

        ! Define the variables
        real :: a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3
        real :: denom

        !SYSTEM OF EQUATIONS
        !	eq1 = a1*mxx + b1*myy + c1*mxy - d1 == 0;
        !	eq2 = a2*mxx + b2*myy + c2*mxy - d2 == 0;
        !	eq3 = a3*mxx + b3*myy + c3*mxy - d3 == 0;

        a1 = A_sys(1,1); b1 = A_sys(1,2); c1 = A_sys(1,3); d1 = b_sys(1); 
        a2 = A_sys(2,1); b2 = A_sys(2,2); c2 = A_sys(2,3); d2 = b_sys(2); 
        a3 = A_sys(3,1); b3 = A_sys(3,2); c3 = A_sys(3,3); d3 = b_sys(3); 

        ! Calculate the denominator (common for all three equations)
        denom = a3*b2*c1 - a2*b3*c1 - a3*b1*c2 + a1*b3*c2 + a2*b1*c3 - a1*b2*c3

        ! Calculate mxx, myy, and mxy
        x_11 = -((-b3*c2*d1 + b2*c3*d1 + b3*c1*d2 - b1*c3*d2 - b2*c1*d3 + b1*c2*d3) / denom)

        x_22 = -((a3*c2*d1 - a2*c3*d1 - a3*c1*d2 + a1*c3*d2 + a2*c1*d3 - a1*c2*d3) / denom)

        x_12 = -((-a3*b2*d1 + a2*b3*d1 + a3*b1*d2 - a1*b3*d2 - a2*b1*d3 + a1*b2*d3) / denom)


    end subroutine solve_system

    subroutine finding_incoming_outgoing_pops()
        implicit none

        Incs_b = 0.0d0
        Outs_b = 0.0d0

        !case-1
        Incs_b(1,0:q-1) = (/ 1, 1, 0, 0, 1, 0, 0, 0, 1 /)
        Outs_b(1,0:q-1) = (/ 1, 0, 1, 1, 0, 0, 1, 0, 0 /)

        !case-2
        Incs_b(2,0:q-1) = (/ 1, 0, 0, 1, 1, 0, 0, 1, 0 /)
        Outs_b(2,0:q-1) = (/ 1, 1, 1, 0, 0, 1, 0, 0, 0 /)

        !case-3
        Incs_b(3,0:q-1) = (/ 1, 1, 0, 1, 1, 0, 0, 1, 1 /)
        Outs_b(3,0:q-1) = (/ 1, 1, 1, 1, 0, 1, 1, 0, 0 /)

        !case-4
        Incs_b(4,0:q-1) = (/ 1, 1, 1, 0, 0, 1, 0, 0, 0 /)
        Outs_b(4,0:q-1) = (/ 1, 0, 0, 1, 1, 0, 0, 1, 0 /)

        !case-5
        Incs_b(5,0:q-1) = (/ 1, 1, 1, 0, 1, 1, 0, 0, 1 /)
        Outs_b(5,0:q-1) = (/ 1, 0, 1, 1, 1, 0, 1, 1, 0 /)

        !case-7
        Incs_b(7,0:q-1) = (/ 1, 1, 1, 1, 1, 1, 0, 1, 1 /)
        Outs_b(7,0:q-1) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 0 /)

        !case-8
        Incs_b(8,0:q-1) = (/ 1, 0, 1, 1, 0, 0, 1, 0, 0 /)
        Outs_b(8,0:q-1) = (/ 1, 1, 0, 0, 1, 0, 0, 0, 1 /)

        !case-10
        Incs_b(10,0:q-1) = (/ 1, 0, 1, 1, 1, 0, 1, 1, 0 /)
        Outs_b(10,0:q-1) = (/ 1, 1, 1, 0, 1, 1, 0, 0, 1 /)

        !case-11
        Incs_b(11,0:q-1) = (/ 1, 1, 1, 1, 1, 0, 1, 1, 1 /)
        Outs_b(11,0:q-1) = (/ 1, 1, 1, 1, 1, 1, 1, 0, 1 /)

        !case-12
        Incs_b(12,0:q-1) = (/ 1, 1, 1, 1, 0, 1, 1, 0, 0 /)
        Outs_b(12,0:q-1) = (/ 1, 1, 0, 1, 1, 0, 0, 1, 1 /)

        !case-13
        Incs_b(13,0:q-1) = (/ 1, 1, 1, 1, 1, 1, 1, 0, 1 /)
        Outs_b(13,0:q-1) = (/ 1, 1, 1, 1, 1, 0, 1, 1, 1 /)

        !case-14
        Incs_b(14,0:q-1) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 0 /)
        Outs_b(14,0:q-1) = (/ 1, 1, 1, 1, 1, 1, 0, 1, 1 /)


    end subroutine finding_incoming_outgoing_pops



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
        real(8),dimension(1:nx, 1:ny) ::  var_mxxp, var_myyp, var_mxyp

        nprocs = 1
        nblocks = nprocs
        do i = 1, nx
            do j = 1, ny 
                var_mxxp(i, j) =  sum(f(i, j, :) * Hxx_p(i, j, :)) / rho(i, j)
                var_myyp(i, j) =  sum(f(i, j, :) * Hyy_p(i, j, :)) / rho(i, j)
                var_mxyp(i, j) =  sum(f(i, j, :) * Hxy_p(i, j, :)) / rho(i, j)
            end do
        end do

        open(unit = 200, form = 'unformatted', file=output_dir_name//'/'//trim(filename))
            write(200) nblocks
            write(200) (nx, ny, 6, l = 0, nprocs-1)
            do m = 0, nprocs-1
                write(200) &
                    & (( sngl(rho(i, j)), i = 1, nx), j = 1, ny), &
                    & (( sngl(ux(i, j)), i = 1, nx), j = 1, ny), &
                    & (( sngl(uy(i, j)), i = 1, nx), j = 1, ny), &
                    & (( sngl(var_mxxp(i, j)), i = 1, nx), j = 1, ny), &
                    & (( sngl(var_myyp(i, j)), i = 1, nx), j = 1, ny), &
                    & (( sngl(var_mxyp(i, j)), i = 1, nx), j = 1, ny)
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
        namelist/Controls/uo,iplot,max_iter,isave,irestart,statsbegin,statsend,iplotstats, cycle_period
        namelist/LogicalControls/post_process, x_periodic,y_periodic,channel_with_cylinder,channel_with_square, &
            & incomp,vel_interp, mom_interp, rotated_coordinate, mom_collision

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
        allocate(ux_prev(1:nx, 1:ny), uy_prev(1:nx, 1:ny))
        allocate(mxx(1:nx, 1:ny), myy(1:nx, 1:ny), mxy(1:nx, 1:ny))
        allocate(mxxp(1:nx, 1:ny), myyp(1:nx, 1:ny), mxyp(1:nx, 1:ny))
        allocate(er(1:nx, 1:ny),er1(1:nx, 1:ny))
        allocate(scalar(1:nx, 1:ny))
        allocate(Incs_b(0:15, 0:q-1),Outs_b(0:15, 0:q-1))
        allocate(cth_glb(1:nx, 1:ny),sth_glb(1:nx, 1:ny))
        allocate(c2th_glb(1:nx, 1:ny),s2th_glb(1:nx, 1:ny))
        allocate(cx_p(1:nx, 1:ny, 0:q-1), cy_p(1:nx, 1:ny, 0:q-1))
        allocate(Hxx_p(1:nx, 1:ny, 0:q-1), Hyy_p(1:nx, 1:ny, 0:q-1), Hxy_p(1:nx, 1:ny, 0:q-1))

    end subroutine allocate_memory

    subroutine allocate_cylinder_memory()
        implicit none

        !Allocating variables after the get_coord subroutine execution:
        allocate(bou_i(1:nbct),bou_j(1:nbct))
        allocate(label_bc(1:nbct),theta_c(1:nbct))
        allocate(p_mean_cy(1:nbct), tw_mean_cy(1:nbct))
        allocate(uxp_b(1:nbct),uyp_b(1:nbct))
        allocate(radii(1:nbct))
        allocate(xb(1:nbct),yb(1:nbct))
        allocate(delta_uk(1:nbct),pb(1:nbct),p_fit(1:nbct),tw_fit(1:nbct))
        allocate(x_ref1(1:nbct),y_ref1(1:nbct))
        allocate(x_ref2(1:nbct),y_ref2(1:nbct))
        allocate(x_ref3(1:nbct),y_ref3(1:nbct))
        allocate(x_ref4(1:nbct),y_ref4(1:nbct))
        allocate(i_p2(1:4,1:nbct),j_p2(1:4,1:nbct))
        allocate(i_p3(1:4,1:nbct),j_p3(1:4,1:nbct))
        allocate(i_p4(1:4,1:nbct),j_p4(1:4,1:nbct))
        allocate(var_ref1(1:nbct),ux_ref1(1:nbct),uy_ref1(1:nbct),mxxp_ref1(1:nbct),myyp_ref1(1:nbct))
        allocate(var_ref2(1:nbct),ux_ref2(1:nbct),uy_ref2(1:nbct),mxxp_ref2(1:nbct),myyp_ref2(1:nbct))
        allocate(var_ref3(1:nbct),ux_ref3(1:nbct),uy_ref3(1:nbct),mxxp_ref3(1:nbct),myyp_ref3(1:nbct))
        allocate(var_ref4(1:nbct),ux_ref4(1:nbct),uy_ref4(1:nbct),mxxp_ref4(1:nbct),myyp_ref4(1:nbct))
        allocate(idx(1:nbct), thetas_cy(1:nbct), ps_cy(1:nbct), mxyps_cy(1:nbct), tws_cy(1:nbct))
        allocate(theta_cy(1:nbct), p_cy(1:nbct), mxyp_cy(1:nbct), tw_cy(1:nbct))
        allocate(p_costheta(1:nbct), p_sintheta(1:nbct), tw_costheta(1:nbct), tw_sintheta(1:nbct))
        allocate(cth_cyl(1:nbct),sth_cyl(1:nbct), cths_cyl(1:nbct),sths_cyl(1:nbct))
        allocate(F_drag(1:nbct), F_lift(1:nbct))

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

        open(unit=10,file='output/master.p3d')
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

    subroutine calculate_norm(l2norm)
        implicit none
        integer :: i, j 
        double precision, intent(out) :: l2norm 
        double precision :: numerator, denominator

        numerator = 0.0d0
        denominator = 0.0d0

        do j = 1, ny
            do i = 1, nx 

            numerator = numerator + sqrt( (ux(i, j) - ux_prev(i, j))**2 + (uy(i, j) - uy_prev(i, j))**2)
            denominator = denominator + sqrt( ux(i, j)**2 + uy(i, j)**2 )

            end do
        end do

        l2norm = numerator/denominator
        ux_prev = ux
        uy_prev = uy 

    end subroutine calculate_norm

end program lbm_2d
