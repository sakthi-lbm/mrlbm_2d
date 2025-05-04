program statistics_calculation
    implicit none
    real(8), parameter :: PI = 4.0d0*atan(1.0d0)
    integer :: unit, ios, count
    character(len=256) :: line
    integer :: nstep, nbct
    integer, parameter :: nmax = 900000
    integer, allocatable, dimension(:) :: time
    real(8), allocatable, dimension(:) :: fd_signal,fl_signal
    real(8), allocatable, dimension(:) :: p1_signal,p2_signal,p3_signal,p4_signal
    real(8), allocatable, dimension(:) :: ux1_signal,ux2_signal,ux3_signal,ux4_signal
    real(8), allocatable, dimension(:) :: uy1_signal,uy2_signal,uy3_signal,uy4_signal
    real(8), allocatable, dimension(:,:) :: ps, ps_fit
    real(8), allocatable, dimension(:) :: theta, ps_avg, ps_fit_avg
    real(8), allocatable, dimension(:) :: Cp, Cp_fit, Cd, Cl
    real(8) :: uo, rho_infty, p_infty, r_cy, D_cy, Cd_avg, Cl_avg
    real(8) :: p_norm, f_norm, Str
    integer :: peak1, peak2, n_cycle
    integer :: period_l, period_d, period1, period2, period3, period4
    integer :: cycle_start_l, cycle_start_d, cycle_start1, cycle_start2, cycle_start3, cycle_start4
    integer :: cycle_end_l, cycle_end_d, cycle_end1, cycle_end2, cycle_end3, cycle_end4
    logical, parameter :: verbose = .true.
    integer :: i, j, k, mean_counter

    n_cycle = 7
    unit=100
    count = 0

    call read_input_file()
    p_infty = rho_infty/3.0d0
    D_cy = 2.0d0*r_cy
    p_norm = 0.50d0 * rho_infty * (uo**2)
    f_norm = 0.50d0 * rho_infty * (uo**2) * D_cy

    !finding the number of lines
    ! open(unit=unit, file='data_probe/test.dat', status='old', action='read')
    open(unit=unit, file='data_mean/forces_time.dat', status='old', action='read')
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            count = count + 1
        end do
    close(unit)
    nstep = count

    call allocate_variables()


    ! Read the data
    
    ! open(unit=unit, file='data_probe/test.dat', status='old')
    !     do i = 1, nstep
    !         read(unit, *) time(i), p1_signal(i), p2_signal(i), p3_signal(i), p4_signal(i)
    !     end do
    ! close(unit)
    ! open(unit=unit, file='data_probe/test.dat', status='old')
    !     do i = 1, nstep
    !         read(unit, *) time(i), ux1_signal(i), ux2_signal(i), ux3_signal(i), ux4_signal(i)
    !     end do
    ! close(unit)
    ! open(unit=unit, file='data_probe/test.dat', status='old')
    !     do i = 1, nstep
    !         read(unit, *) time(i), uy1_signal(i), uy2_signal(i), uy3_signal(i), uy4_signal(i)
    !     end do
    ! close(unit)

    ! !Point-1
    
    ! call find_peaks(nstep, time, uy1_signal, cycle_start1, period1)
    ! call find_peaks(nstep, time, uy2_signal, cycle_start2, period2)
    ! call find_peaks(nstep, time, uy3_signal, cycle_start3, period3)
    ! call find_peaks(nstep, time, uy4_signal, cycle_start4, period4)
    ! print*,cycle_start_l, period_l
    ! print*,cycle_start1, period1
    ! print*,cycle_start2, period2
    ! print*,cycle_start3, period3
    ! print*,cycle_start4, period4

    !Averaging
    !=======================================================================================================
    open(unit=unit, file='data_mean/theta.dat', status='old')
        do i = 1, nbct
            read(unit, *) theta(i)
        end do
    close(unit)
    open(unit=unit, file='data_mean/forces_time.dat', status='old')
        do i = 1, nstep
            read(unit, *) time(i), fd_signal(i), fl_signal(i)
        end do
    close(unit)
    call find_peaks(nstep, time, fl_signal, cycle_start_l, period_l)
    
    cycle_end_l = cycle_start_l + (n_cycle*period_l)
    !Strouhal number
    Str = D_cy/(period_l * uo)

    open(unit=unit, file='data_mean/p_surface_time.dat', status='old', &
            &action='read', form='unformatted')
        do i = 1, nstep
            read(unit) time(i), ( ps(i, k), k = 1, nbct)
        end do
    close(unit)

    open(unit=unit, file='data_mean/p_surface_smooth_time.dat', status='old', &
            & action='read', form='unformatted')
        do i = 1, nstep
            read(unit) time(i), ( ps_fit(i, k), k = 1, nbct)
        end do
    close(unit)

    do i = cycle_start_l, cycle_end_l
        mean_counter = i-cycle_start_l
        do k = 1, nbct
            ps_avg(k) = (mean_counter*ps_avg(k) + ps(i,k))/(mean_counter + 1)
            ps_fit_avg(k) = (mean_counter*ps_fit_avg(k) + ps_fit(i,k))/(mean_counter + 1)

            Cp(k) = (ps_avg(k) - p_infty)/p_norm
            Cp_fit(k) = (ps_fit_avg(k) - p_infty)/p_norm
        end do

        Cd(i) = fd_signal(i)/f_norm
        Cl(i) = fl_signal(i)/f_norm

        Cd_avg = (mean_counter*Cd_avg + Cd(i))/(mean_counter + 1)
        Cl_avg = (mean_counter*Cl_avg + (Cl(i)**2))/(mean_counter + 1)

    end do
    Cl_avg = sqrt(Cl_avg)

    ! open(unit=unit, file='plots/cp.dat')
    !     do  k = 1, nbct
    !         write(unit,*) theta(k), Cp(k)
    !     end do
    ! close(unit)

    ! open(unit=unit, file='plots/cp_fit.dat')
    !     do  k = 1, nbct
    !         write(unit,*) theta(k), Cp_fit(k)
    !     end do
    ! close(unit)

    open(unit=unit, file='plots/cp.dat')
        do  k = 1, nbct
            if(theta(k) .le. PI) then
                write(unit,*) abs((theta(k)*180.0d0/PI)-180.0d0), Cp(k)
            end if
        end do
    close(unit)

    open(unit=unit, file='plots/cp_fit.dat')
    do  k = 1, nbct
        if(theta(k) .le. PI) then
            write(unit,*) abs((theta(k)*180.0d0/PI)-180.0d0), Cp_fit(k)
        end if
    end do
    close(unit)

    open(unit=unit, file='plots/cd_time.dat')
        do  i = cycle_start_l, cycle_end_l
            write(unit,*) real(i-cycle_start_l)/real(period_l), Cd(i)
        end do
    close(unit)
    open(unit=unit, file='plots/cl_time.dat')
        do  i = cycle_start_l, cycle_end_l
            write(unit,*) real(i-cycle_start_l)/real(period_l), Cl(i)
        end do
    close(unit)
    open(unit=unit, file='plots/drag_lift.dat')
        write(unit,*) "Strouhal number (St)  ",":", Str
        write(unit,*) "Drag coefficinet (C_D)  ",":", Cd_avg
        write(unit,*) "Lift coefficinet (C_L ",":", Cl_avg
    close(unit)


contains

subroutine find_peaks(n, step, signal, cycle_start, period)
    implicit none
    integer, intent(in) ::  n
    integer, dimension(1:n), intent(in) :: step
    real(8), dimension(1:n), intent(in) :: signal
    integer, intent(out) :: cycle_start, period
    integer :: i

    ! Find first two local maxima (peaks)
    peak1 = -1
    peak2 = -1
    do i = 2, n-1
        if (signal(i) > signal(i-1) .and. signal(i) > signal(i+1)) then
            if (peak1 == -1) then
                peak1 = i
            elseif (peak2 == -1) then
                peak2 = i
                exit
            end if
        end if
    end do

    if (peak1 == -1 .or. peak2 == -1) then
        print *, 'Error: Could not detect two peaks.'
        stop
    end if

    ! Compute and print cycle information
    cycle_start = peak1
    period = peak2 - peak1

end subroutine find_peaks

subroutine allocate_variables()
    implicit none

    allocate(time(1:nstep))
    allocate(fd_signal(1:nstep), fl_signal(1:nstep))
    allocate(p1_signal(1:nstep), p2_signal(1:nstep), p3_signal(1:nstep), p4_signal(1:nstep))
    allocate(ux1_signal(1:nstep), ux2_signal(1:nstep), ux3_signal(1:nstep), ux4_signal(1:nstep))
    allocate(uy1_signal(1:nstep), uy2_signal(1:nstep), uy3_signal(1:nstep), uy4_signal(1:nstep))
    allocate(ps(1:nstep,1:nbct), ps_fit(1:nstep,1:nbct))
    allocate(theta(1:nbct), ps_avg(1:nbct), ps_fit_avg(1:nbct))
    allocate(Cd(1:nstep), Cl(1:nstep))
    allocate(Cp(1:nbct), Cp_fit(1:nbct))
    

end subroutine allocate_variables

subroutine read_input_file()
    implicit none
    integer :: input = 100
    integer :: iread_error = 0
    integer :: i 

    namelist/parameters/nbct, r_cy, uo, rho_infty

    300 format("Error while reading input.dat file...")
    150 if (iread_error .ne. 0) then 
            write(*, 300); stop 
        end if
    200 format(a, T5)

    open (unit = input, file = "input_stat.dat", status = "old")
        write(*, *) 'Reading input_stat.dat file ...'

        read(input, parameters, iostat=iread_error, err=150)
        write(*, 200) 'parameters'

    close(input)

    ! if (verbose) then 
    !     write(*, parameters)
    ! end if

end subroutine read_input_file

end program statistics_calculation

