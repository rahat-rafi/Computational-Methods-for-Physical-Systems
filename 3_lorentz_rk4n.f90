program lorentz_rk4
    implicit none
    
    integer, parameter :: dp = selected_real_kind(15, 300)
    integer, parameter :: n = 6  ! 3 positions + 3 velocities
    
    ! All the variables needed
    real(dp) :: y(n), err(n), t, dt, t_end
    real(dp) :: r0(3), v0(3), exact_x, exact_y, exact_z
    real(dp) :: error_x, error_y, error_z, max_error
    integer :: i, n_steps, file_unit
    character(len=50) :: filename
    
    ! I created an array for the values of dt
    real(dp) :: dt_list(4) = [1.0_dp, 0.1_dp, 1.0e-4_dp, 1.0e-5_dp]
    integer :: dt_index
    
    ! Initial conditions for r & v
    r0 = [-1.0_dp, 0.0_dp, 0.0_dp]  ! Initial position
    v0 = [0.0_dp, 1.0_dp, 1.0_dp]   ! Initial velocity
    
    t_end = 100.0_dp
    
    do dt_index = 1, size(dt_list)
        dt = dt_list(dt_index)
        
        write(filename, '(A, F0.5, A)') '3_lorentz_dt_', dt, '.dat'
        
        ! Resetting initial conditions from the question
        y(1:3) = r0
        y(4:6) = v0
        
        n_steps = nint(t_end / dt)
        
        open(newunit=file_unit, file=trim(filename), status='replace')
        write(file_unit, '(A)') 't x y z max_error error_x'
        
        ! I am writing the initial value here
        t = 0.0_dp
        max_error = 0
        error_x = 0
        write(file_unit, '(11ES15.6)') t, y(1), y(2), y(3), max_error, error_x
        do i = 1, n_steps
            call rk4n(dt, y, n, err)
            t = t + dt
            
            ! Calculating the exact solution for the differential equation
            exact_x = -cos(t)
            exact_y = sin(t)
            exact_z = t
            
            ! Error calculation
            error_x = y(1) - exact_x
            error_y = y(2) - exact_y
            error_z = y(3) - exact_z
            max_error = maxval([abs(error_x), abs(error_y), abs(error_z)])
            
            write(file_unit, '(11ES15.6)') t, y(1), y(2), y(3), max_error, error_x
        end do
        
        close(file_unit)
    end do
    
    do dt_index = 1, size(dt_list)
        dt = dt_list(dt_index)
        write(filename, '(A, F0.5, A)') '3_lorentz_dt_', dt, '.dat'
    end do

contains

! ===== RK4N Subroutine =====
subroutine rk4n(dx, y, n, err)
    real(dp), intent(in) :: dx
    integer, intent(in) :: n
    real(dp), intent(inout) :: y(n)
    real(dp), intent(out) :: err(n)
    
    real(dp) :: y1(n), y2(n), y3(n)
    real(dp) :: k1(n), k2(n), k3(n), k4(n)
    real(dp) :: rksum
    integer :: ivar
    
    call derivs(y, k1, n)
    k1 = k1 * dx
    y1 = y + 0.5_dp * k1
    
    call derivs(y1, k2, n)
    k2 = k2 * dx
    y2 = y + 0.5_dp * k2
    
    call derivs(y2, k3, n)
    k3 = k3 * dx
    y3 = y + k3
    
    call derivs(y3, k4, n)
    k4 = k4 * dx
    
    ! Update and compute error
    do ivar = 1, n
        rksum = (1.0_dp/6.0_dp) * (k1(ivar) + 2.0_dp*k2(ivar) + 2.0_dp*k3(ivar) + k4(ivar))
        err(ivar) = abs(rksum - k2(ivar))
        y(ivar) = y(ivar) + rksum
    end do
    
end subroutine rk4n

! ===== Derivatives =====
subroutine derivs(y, dydt, n)
    integer, intent(in) :: n
    real(dp), intent(in) :: y(n)
    real(dp), intent(out) :: dydt(n)
    
    dydt(1) = y(4)  
    dydt(2) = y(5)   
    dydt(3) = y(6)  
    
    dydt(4) = y(5)   
    dydt(5) = -y(4)  
    dydt(6) = 0.0_dp 
    
end subroutine derivs

end program lorentz_rk4
