program detect_vortex_cycle
    implicit none
    integer, parameter :: nmax = 900000
    real(8), dimension(nmax) :: time
     real(8), dimension(nmax) :: signal
    real(8), dimension(nmax) :: p1_signal,p2_signal,p3_signal,p4_signal
    real(8), dimension(nmax) :: ux1_signal,ux2_signal,ux3_signal,ux4_signal
    real(8), dimension(nmax) :: uy1_signal,uy2_signal,uy3_signal,uy4_signal
    integer :: i, n, peak1, peak2
    real(8) :: period

    ! Read the data
    open(unit=10, file='data_probe/uy_probe.dat', status='old')
    n = 0
    do
        n = n + 1
        read(10, *, end=100) time(n), p1_signal(n), p2_signal(n), p3_signal(n), p4_signal(n)
    end do
	100 continue
    close(10)
    n = n - 1

    ! Find first two local maxima (peaks)
    peak1 = -1
    peak2 = -1
    signal = p4_signal
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
    period = time(peak2) - time(peak1)
    print *, 'Start of vortex cycle: t =', time(peak1)
    print *, 'End   of vortex cycle: t =', time(peak2)
    print *, 'Shedding period      T =', period

end program detect_vortex_cycle

