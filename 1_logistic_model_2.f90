program logistic_decreasing_k
    implicit none

    ! all the variables necessary are here
    integer, parameter :: num_generations = 50, num_r = 4
    real, parameter :: initial_x = 0.1, k0 = 1.0, k_decrease = 0.01
    real :: r_values(num_r), x_current(num_r), p_current(num_r)
    real :: k_current
    integer :: generation, r_index, file_unit, iostat
    
    ! given values in the question
    r_values = [1.9, 2.9, 3.3, 3.6]
    k_current = k0
    
    ! how population is denoted in the equation
    ! x = p / K, so p = x * K
    x_current = initial_x
    p_current = x_current * k_current
    
    ! necessary output file for plotting
    open(newunit=file_unit, file='1_logistic_model_2.txt', &
         status='replace', action='write', iostat=iostat)
    if (iostat /= 0) then
        error stop "Error opening output file!"
    end if
    
    ! header for the output file
    write(file_unit, '(A10,4A15,4A15,A15)') 'Generation', 'x_r=1.9', 'x_r=2.9', &
    	'x_r=3.3', 'x_r=3.6', 'p_r=1.9', 'p_r=2.9', 'p_r=3.3', 'p_r=3.6', 'K'
    
    ! for generation 1
    write(file_unit, '(I10,4F15.6,4F15.6,F15.6)') 1, x_current, p_current, k_current
    
    ! for generations 2 to 50
    do generation = 2, num_generations
        k_current = k_current - k_decrease
        
        do r_index = 1, num_r
            x_current(r_index) = r_values(r_index) * x_current(r_index) * (1.0 - x_current(r_index))
            
            p_current(r_index) = x_current(r_index) * k_current
        end do
        
        write(file_unit, '(I10,4F15.6,4F15.6,F15.6)') generation, x_current, p_current, k_current
        
    end do
    
    close(file_unit)
    print *, "Simulation complete. Data written to 1_logistic_model_2.txt"
    print *, "Final carrying capacity K =", k_current

end program logistic_decreasing_k
