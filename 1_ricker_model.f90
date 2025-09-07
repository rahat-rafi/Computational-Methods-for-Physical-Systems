program ricker_model
    implicit none

    ! all the variables necessary are here
    integer, parameter :: num_generations = 50, num_r = 4
    real, parameter :: initial_x = 1.1, small_threshold = 0.2   ! when fractional populatiion is 20%
    real :: r_values(num_r), x_current(num_r), r_prime(num_r), growth_factor
    integer :: generation, r_index, file_unit, iostat
    
    ! given values in the question
    r_values = [1.9, 2.9, 3.3, 3.6]
    
    ! how population is denoted in the equation
    ! x = p / K, so p = x * K
    x_current = initial_x
    
     ! calculating r' = log(r) as specified
    do r_index = 1, num_r
        r_prime(r_index) = log(r_values(r_index))
    end do
    
    ! necessary output file for plotting
    open(newunit=file_unit, file='1_Ricker_model.txt', status='replace', action='write', iostat=iostat)
    if (iostat /= 0) then
        error stop "Error opening output file!"
    end if
    
    ! header for the output file
    write(file_unit, '(A10,4A15,4A15,A15)') 'Generation', 'x_r=1.9', 'x_r=2.9', 'x_r=3.3', 'x_r=3.6'
    
    ! for generation 1
    write(file_unit, '(I10,4F15.6,4F15.6,F15.6)') 1, x_current
    
    ! for generations 2 to 50
    do generation = 2, num_generations
        do r_index = 1, num_r
             ! Use r' = log(r) only when x is small, otherwise use original r
            if (x_current(r_index) < small_threshold) then
                growth_factor = r_prime(r_index)
            else
                growth_factor = r_values(r_index)
            end if
            
            x_current(r_index) = x_current(r_index) * exp(growth_factor * (1.0 - x_current(r_index)))
        end do
        
        write(file_unit, '(I10,4F15.6,4F15.6,F15.6)') generation, x_current
        
    end do
    
    close(file_unit)
    print *, "Simulation complete. Data written to 1_Ricker_model.txt"

end program ricker_model
