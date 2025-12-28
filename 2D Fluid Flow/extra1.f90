module relax_cfd
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)
contains

  ! psi_ij = (1/5)[ SC + (1/4) SS ] + (1/5) psi_ij + (1/40) SCO
  subroutine relax_psi_sweep(psi, vort, nx, ny, wrel, ixst, ixfn, iyst, iyfn, dpmx)
    integer, intent(in)    :: nx, ny
    integer, intent(inout) :: ixst, ixfn, iyst, iyfn
    real(dp), intent(inout):: psi(nx,ny)
    real(dp), intent(in)   :: vort(nx,ny), wrel
    real(dp), intent(out)  :: dpmx
    integer :: i, j, incx, incy
    real(dp) :: oldv, updt, delta, SC, SS, SCO

    dpmx = 0.0_dp
    incx = merge( 1, -1, ixfn>=ixst )
    incy = merge( 1, -1, iyfn>=iyst )

    do j = iyst, iyfn, incy
      do i = ixst, ixfn, incx
        oldv = psi(i,j)

        SC  = psi(i+1,j)   + psi(i-1,j)   + psi(i,j+1)   + psi(i,j-1)
        SS  = psi(i+1,j+1) + psi(i-1,j+1) + psi(i-1,j-1) + psi(i+1,j-1)
        SCO = vort(i+1,j)  + vort(i-1,j)  + vort(i,j+1)  + vort(i,j-1)

        updt     = (SC + 0.25_dp*SS)/5.0_dp + vort(i,j)/5.0_dp + SCO/40.0_dp
        psi(i,j) = (1.0_dp - wrel)*oldv + wrel*updt

        delta = abs(psi(i,j) - oldv)
        if (delta > dpmx) dpmx = delta
      end do
    end do
  end subroutine relax_psi_sweep

  ! del^2 omega = R * (psi_y omega_x − psi_x omega_y)
  ! omega_ij = (1/5)[ SC + (1/4) SS ] - (3R/10) * F
  ! F = 0.25 * [ (psi_{i,j+1}-psi_{i,j-1})(omega_{i+1,j}-omega_{i-1,j})
  !             - (psi_{i+1,j}-psi_{i-1,j})(omega_{i,j+1}-omega_{i,j-1}) ]
  subroutine relax_vort_sweep(psi, vort, nx, ny, wrel, R, ixst, ixfn, iyst, iyfn, domx)
    integer, intent(in)    :: nx, ny
    integer, intent(inout) :: ixst, ixfn, iyst, iyfn
    real(dp), intent(in)   :: psi(nx,ny), wrel, R
    real(dp), intent(inout):: vort(nx,ny)
    real(dp), intent(out)  :: domx
    integer :: i, j, incx, incy
    real(dp) :: oldv, updt, delta, SCO, SSO, Fterm

    domx = 0.0_dp
    incx = merge( 1, -1, ixfn>=ixst )
    incy = merge( 1, -1, iyfn>=iyst )

    do j = iyst, iyfn, incy
      do i = ixst, ixfn, incx
        oldv = vort(i,j)

        SCO = vort(i+1,j)   + vort(i-1,j)   + vort(i,j+1)   + vort(i,j-1)
        SSO = vort(i+1,j+1) + vort(i-1,j+1) + vort(i-1,j-1) + vort(i+1,j-1)

        Fterm = 0.25_dp * ( (psi(i  ,j+1)-psi(i  ,j-1))*(vort(i+1,j)-vort(i-1,j))  &
                          - (psi(i+1,j)-psi(i-1,j))*(vort(i  ,j+1)-vort(i  ,j-1)) )

        updt      = (SCO + 0.25_dp*SSO)/5.0_dp - 3.0_dp*(R*Fterm)/10.0_dp
        vort(i,j) = (1.0_dp - wrel)*oldv + wrel*updt

        delta = abs(vort(i,j) - oldv)
        if (delta > domx) domx = delta
      end do
    end do
  end subroutine relax_vort_sweep

end module relax_cfd
!======================================================================
program rect_obstruction_wstudy
  use relax_cfd
  implicit none

  !---------------- geometry parameters (Case 1) ----------------------
  integer, parameter :: n1 = 400   ! upstream length
  integer, parameter :: n2 = 400   ! downstream length
  integer, parameter :: nw = 50    ! block width
  integer, parameter :: nt = 50    ! block height from bottom
  integer, parameter :: nh = 400   ! height above block

  integer, parameter :: nx = n1 + nw + n2 + 1   ! = 851  (x: 1..851)
  integer, parameter :: ny = nt + nh + 1        ! = 451  (y: 1..451)

  real(dp), parameter :: R         = 1.0_dp
  integer,  parameter :: niter     = 200000
  integer,  parameter :: iter_print = 5000
  integer,  parameter :: iter_switch = 20000

  real(dp), parameter :: w_basic   = 0.10_dp   ! basic w
  real(dp), parameter :: w_high    = 1.00_dp   ! higher w used after iter_switch

  ! subset spacing for final field output
  integer,  parameter :: stepx = 10, stepy = 10

  ! Running both simulations: baseline and switched-w
  call run_simulation(.false., 'extra1_basic_field_200K_step10.dat')
  call run_simulation(.true.,  'extra1_wswitch_field_200K_step10.dat')

contains

  !====================================================================
  subroutine run_simulation(do_switch, outfile)
    logical,      intent(in) :: do_switch
    character(*), intent(in) :: outfile

    ! field & mask arrays
    real(dp), allocatable :: psi(:,:), vort(:,:)
    logical,  allocatable :: is_solid(:,:)

    ! indices and geometry
    integer :: i, j, iter
    integer :: ix_block_left, ix_block_right, jy_block_top

    ! relaxation diagnostics
    real(dp) :: dpmx, domx, dtmp, wrel

    ! rectangles for fluid regions (I, II, III)
    integer :: ixA1, ixA2, ixB1, ixB2, ixC1, ixC2
    integer :: jyA1, jyA2, jyB1, jyB2, jyC1, jyC2

    ! for velocity output
    real(dp) :: u, v, x, y
    integer, parameter :: uout = 20

    !--------------------------- allocating arrays -----------------------------
    allocate(psi(nx,ny), vort(nx,ny), is_solid(nx,ny))

    !------------------------ defining obstacle indices ------------------------
    ix_block_left  = n1 + 1          ! left edge of block  (ix = 401)
    ix_block_right = n1 + nw + 1     ! right edge of block (ix = 451)
    jy_block_top   = nt + 1          ! top of block        (iy = 51)

    !-------------------------- initial uniform flow -------------------------
    ! Starting from unimpeded uniform flow everywhere: psi = y, vort = 0
    do j = 1, ny
       do i = 1, nx
          psi(i,j)     = real(j-1, dp)
          vort(i,j)    = 0.0_dp
          is_solid(i,j)= .false.
       end do
    end do

    !-------------------------- marking rectangular obstacle -------------------
    do j = 1, jy_block_top
       do i = ix_block_left, ix_block_right
          is_solid(i,j) = .true.
          psi(i,j)      = 0.0_dp          ! no-slip: psi constant at centerline value
          vort(i,j)     = 0.0_dp
       end do
    end do

    !---------------------- applying Case 1 BCs once ---------------------------
    call apply_bc_case1(psi, vort, nx, ny, is_solid, n1, n2, nw, nt, nh)

    !------------------------- defining fluid rectangles -----------------------
    ! Region I: left of the obstruction
    ixA1 = 2
    ixA2 = n1
    jyA1 = 2
    jyA2 = nt + 1

    ! Region II: right of the obstruction
    ixB1 = n1 + nw + 2
    ixB2 = n1 + n2 + nw
    jyB1 = 2
    jyB2 = nt + 1

    ! Region III: above the obstruction
    ixC1 = 2
    ixC2 = n1 + n2 + nw
    jyC1 = nt + 2
    jyC2 = nt + nh            ! = ny - 1

    !---------------------- announcing run type --------------------------------
    if (.not. do_switch) then
       write(*,'(/,"===== BASELINE RUN: constant w =",f6.2," =====")') w_basic
    else
       write(*,'(/,"===== SWITCHED-w RUN: w=",f6.2," then ",f6.2," after iter ",i8," =====")') &
            w_basic, w_high, iter_switch
    end if
    call flush(6)

    !---------------------------- main iteration loop ------------------------
    do iter = 1, niter

       ! choosing wrel depending on whether we switch or not
       if (.not. do_switch) then
          wrel = w_basic
       else
          if (iter <= iter_switch) then
             wrel = w_basic
          else
             wrel = w_high
          end if
       end if

       ! 1) Applying outer BCs and enforce psi in solid region
       call apply_bc_case1(psi, vort, nx, ny, is_solid, n1, n2, nw, nt, nh)

       !---------------- psi relaxation (forward + backward) ----------------
       dpmx = 0.0_dp

       ! forward sweeps
       call relax_psi_sweep(psi, vort, nx, ny, wrel, ixA1, ixA2, jyA1, jyA2, dtmp)
       dpmx = max(dpmx, dtmp)

       call relax_psi_sweep(psi, vort, nx, ny, wrel, ixB1, ixB2, jyB1, jyB2, dtmp)
       dpmx = max(dpmx, dtmp)

       call relax_psi_sweep(psi, vort, nx, ny, wrel, ixC1, ixC2, jyC1, jyC2, dtmp)
       dpmx = max(dpmx, dtmp)

       ! backward sweeps
       call relax_psi_sweep(psi, vort, nx, ny, wrel, ixC2, ixC1, jyC1, jyC2, dtmp)
       dpmx = max(dpmx, dtmp)

       call relax_psi_sweep(psi, vort, nx, ny, wrel, ixB2, ixB1, jyB1, jyB2, dtmp)
       dpmx = max(dpmx, dtmp)

       call relax_psi_sweep(psi, vort, nx, ny, wrel, ixA2, ixA1, jyA1, jyA2, dtmp)
       dpmx = max(dpmx, dtmp)

       ! Re-enforcing psi BCs after psi relaxation
       call apply_bc_case1(psi, vort, nx, ny, is_solid, n1, n2, nw, nt, nh)

       ! Applying vorticity BCs at object walls
       call vortbc_rect(psi, vort, nx, ny, n1, nw, nt)

       !---------------- vorticity relaxation (forward + backward) ----------
       domx = 0.0_dp

       ! forward sweeps
       call relax_vort_sweep(psi, vort, nx, ny, wrel, R, ixA1, ixA2, jyA1, jyA2, dtmp)
       domx = max(domx, dtmp)

       call relax_vort_sweep(psi, vort, nx, ny, wrel, R, ixB1, ixB2, jyB1, jyB2, dtmp)
       domx = max(domx, dtmp)

       call relax_vort_sweep(psi, vort, nx, ny, wrel, R, ixC1, ixC2, jyC1, jyC2, dtmp)
       domx = max(domx, dtmp)

       ! backward sweeps
       call relax_vort_sweep(psi, vort, nx, ny, wrel, R, ixC2, ixC1, jyC1, jyC2, dtmp)
       domx = max(domx, dtmp)

       call relax_vort_sweep(psi, vort, nx, ny, wrel, R, ixB2, ixB1, jyB1, jyB2, dtmp)
       domx = max(domx, dtmp)

       call relax_vort_sweep(psi, vort, nx, ny, wrel, R, ixA2, ixA1, jyA1, jyA2, dtmp)
       domx = max(domx, dtmp)

       ! Enforcing BCs again for vorticity & psi at boundaries and object
       call apply_bc_case1(psi, vort, nx, ny, is_solid, n1, n2, nw, nt, nh)
       call vortbc_rect(psi, vort, nx, ny, n1, nw, nt)

       ! Progress display
       if (mod(iter, iter_print) == 0) then
          write(*,'("iter=",i8,"  wrel=",f6.2,"  dpmx=",1pe12.5,"  domx=",1pe12.5)') &
               iter, wrel, dpmx, domx
          call flush(6)
       end if

    end do   ! end iter loop

    write(*,'("Run finished for ",a,", final dpmx=",1pe12.5," domx=",1pe12.5)') trim(outfile), dpmx, domx
    call flush(6)

    !----------------------- output final field (x,y,u,v,psi,vort) ----------
    open(unit=uout, file=outfile, status='replace')
    write(uout,'("# x  y  u  v  psi  vort")')

    do j = 2, ny-1, stepy
       y = real(j-1, dp)
       do i = 2, nx-1, stepx
          ! skip obstacle interior
          if (.not. is_solid(i,j)) then
             x = real(i-1, dp)
             u = 0.5_dp * ( psi(i, j+1) - psi(i, j-1) )   ! u = dψ/dy
             v = -0.5_dp * ( psi(i+1, j) - psi(i-1, j) )  ! v = -dψ/dx
             write(uout,'(6f14.6)') x, y, u, v, psi(i,j), vort(i,j)
          end if
       end do
       write(uout,*)
    end do

    call flush(uout)
    close(uout)

    deallocate(psi, vort, is_solid)

  end subroutine run_simulation

  !====================================================================
  ! Case 1 boundary conditions
  subroutine apply_bc_case1(psi, vort, nx, ny, is_solid, n1, n2, nw, nt, nh)
    integer,  intent(in)    :: nx, ny, n1, n2, nw, nt, nh
    real(dp), intent(inout) :: psi(nx,ny), vort(nx,ny)
    logical,  intent(in)    :: is_solid(nx,ny)
    integer :: i, j

    ! Upstream edge (i=1): unimpeded uniform flow
    do j = 1, ny
       psi(1 ,j)  = real(j-1, dp)
       vort(1 ,j) = 0.0_dp
    end do

    ! Downstream edge (i=nx): unimpeded uniform flow
    do j = 1, ny
       psi(nx,j)  = real(j-1, dp)
       vort(nx,j) = 0.0_dp
    end do

    ! Top edge (j=ny): uniform flow
    do i = 1, nx
       psi(i,ny)  = real(ny-1, dp)
       vort(i,ny) = 0.0_dp
    end do

    ! Bottom left edge (centerline, to left of block)
    do i = 1, n1+1
       psi(i,1)  = 0.0_dp
       vort(i,1) = 0.0_dp
    end do

    ! Bottom right edge (centerline, to right of block)
    do i = n1+nw+1, nx
       psi(i,1)  = 0.0_dp
       vort(i,1) = 0.0_dp
    end do

    ! Ensuring psi is fixed inside the solid block
    do j = 1, ny
       do i = 1, nx
          if (is_solid(i,j)) then
             psi(i,j)  = 0.0_dp
          end if
       end do
    end do
  end subroutine apply_bc_case1

  !====================================================================
  ! Vorticity BCs at object walls (like vortbc in your F77 code)
  subroutine vortbc_rect(psi, vort, nx, ny, n1, nw, nt)
    integer,  intent(in)    :: nx, ny, n1, nw, nt
    real(dp), intent(inout) :: psi(nx,ny), vort(nx,ny)
    integer :: i, j

    ! Top side of the object: iy = nt+1, x from n1+1 to n1+nw+1
    j = nt + 1
    do i = n1+1, n1+nw+1
       psi(i,j)   = 0.0_dp
       vort(i,j)  = -2.0_dp * psi(i,j+1)
    end do

    ! Left side of the object: ix = n1+1, y from 1 to nt+1
    i = n1 + 1
    do j = 1, nt+1
       psi(i,j)   = 0.0_dp
       vort(i,j)  = -2.0_dp * psi(i-1,j)
    end do

    ! Right side of the object: ix = n1+nw+1, y from 1 to nt+1
    i = n1 + nw + 1
    do j = 1, nt+1
       psi(i,j)   = 0.0_dp
       vort(i,j)  = -2.0_dp * psi(i+1,j)
    end do

  end subroutine vortbc_rect

end program rect_obstruction_wstudy

