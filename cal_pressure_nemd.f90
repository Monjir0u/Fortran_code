program calculate_pressure_average
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    character(len=30) :: infile, outfile, dir
    real(dp) :: col2_sum, col4_sum, total_avg
    integer :: count, ios
    real(dp) :: col1, col2, col3, col4
    
    dir = "a0.01"
    infile = trim(dir) // "/pressure.dat"
    outfile = "p/pressure_ave_" // trim(dir) // ".dat"
    
    col2_sum = 0.0_dp
    col4_sum = 0.0_dp
    count = 0
    
    open(10, file=infile, status='old', action='read', iostat=ios)
    if (ios /= 0) then
       print *, "Error opening file: ", infile
       stop
    end if
    
    do
       read(10, *, iostat=ios) col1, col2, col3, col4
       if (ios /= 0) exit
       col2_sum = col2_sum + col2
       col4_sum = col4_sum + col4
       count = count + 1
    end do
    close(10)
    
    if (count > 0) then
       total_avg = (col2_sum + col4_sum) / (2 * count)
    else
       print *, "No data read from file."
       stop
    end if
    
    open(20, file=outfile, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
       print *, "Error opening file: ", outfile
       stop
    end if
    
    write(20, '(1X,E14.7)') total_avg
    close(20)
    
    print *, "Averages written to", outfile
  end program calculate_pressure_average
  