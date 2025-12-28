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

        SC  = psi(i+1,j)   + psi(i-1,j)   + psi(i,  j+1) + psi(i,  j-1)
        SS  = psi(i+1,j+1) + psi(i-1,j+1) + psi(i-1,j-1) + psi(i+1,j-1)
        SCO = vort(i+1,j)  + vort(i-1,j)  + vort(i,  j+1) + vort(i,  j-1)

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

        SCO = vort(i+1,j)   + vort(i-1,j)   + vort(i,  j+1) + vort(i,  j-1)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program rect_obstruction_case2
  use relax_cfd
  implicit none

  !---------------- geometry parameters (matching prof’s definitions) ---------
  integer, parameter :: n1 = 400   ! upstream length
  integer, parameter :: n2 = 400   ! downstream length
  integer, parameter :: nw = 50    ! block width
  integer, parameter :: nt = 100    ! block height from bottom
  integer, parameter :: nh = 400   ! height above block

  integer, parameter :: nx = n1 + nw + n2 + 1   ! = 851  (x: 1..851)
  integer, parameter :: ny = nt + nh + 1        ! = 451  (y: 1..451)

  real(dp), parameter :: R    = 1.0_dp
  real(dp), parameter :: wrel = 0.1_dp
  integer,  parameter :: niter = 200000       ! for case 2

  ! "manageable subset" spacing for final velocity field
  integer,  parameter :: stepx = 10, stepy = 10

  !----------------------- field & mask arrays ---------------------------------
  real(dp), allocatable :: psi(:,:), vort(:,:)
  logical,  allocatable :: is_solid(:,:)

  ! indices and geometry
  integer :: i, j, iter
  integer :: ix_block_left, ix_block_right, jy_block_top

  ! relaxation diagnostics
  real(dp) :: dpmx, domx, dtmp

  ! rectangles for fluid regions (I, II, III)
  integer :: ixA1, ixA2, ixB1, ixB2, ixC1, ixC2
  integer :: jyA1, jyA2, jyB1, jyB2, jyC1, jyC2

  ! for velocity / psi output
  real(dp) :: u, v, x, y
  integer  :: ix_snap

  integer, parameter :: uout  = 20   ! final field output unit
  integer, parameter :: usnap = 30   ! psi(ix=390) snapshot output unit

  !--------------------------- allocate arrays ---------------------------------
  allocate(psi(nx,ny), vort(nx,ny), is_solid(nx,ny))

  !------------------------ define obstacle indices ----------------------------
  ix_block_left  = n1 + 1          ! left edge of block  (ix = 401)
  ix_block_right = n1 + nw + 1     ! right edge of block (ix = 451)
  jy_block_top   = nt + 1          ! top of block        (iy = 51)

  !-------------------------- initial uniform flow -----------------------------
  ! Starting from unimpeded uniform flow everywhere: psi = y, vort = 0
  do j = 1, ny
     do i = 1, nx
        psi(i,j)     = real(j-1, dp)
        vort(i,j)    = 0.0_dp
        is_solid(i,j)= .false.
     end do
  end do

  !-------------------------- mark rectangular obstacle -----------------------
  do j = 1, jy_block_top
     do i = ix_block_left, ix_block_right
        is_solid(i,j) = .true.
        psi(i,j)      = 0.0_dp          ! no-slip: psi constant at centerline value
        vort(i,j)     = 0.0_dp
     end do
  end do

  !---------------------- applying Case 2 BCs once --------------------------------
  call apply_bc_case2(psi, vort, nx, ny, is_solid, n1, n2, nw, nt, nh)

  !------------------------- define fluid rectangles ---------------------------
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

  !------------------- open snapshot file for psi(ix=390) ---------------------
  ix_snap = 390
  open(unit=usnap, file='case2_psi_ix390_snap_h100.dat', status='replace', action='write')
  write(usnap,'("# iter  iy   y    u(ix390)   v(ix390)   psi(ix390)   vort(ix390)")')
  call flush(usnap)

  !---------------------------- main iteration loop ---------------------------
  do iter = 1, niter

     ! 1) Applying Case 2 BCs and enforce psi in solid region
     call apply_bc_case2(psi, vort, nx, ny, is_solid, n1, n2, nw, nt, nh)

     !---------------- psi relaxation (forward + backward) --------------------
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
     call apply_bc_case2(psi, vort, nx, ny, is_solid, n1, n2, nw, nt, nh)

     ! Apply vorticity BCs at object walls
     call vortbc_rect(psi, vort, nx, ny, n1, nw, nt)

     !---------------- vorticity relaxation (forward + backward) --------------
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

     ! Enforcing BCs for vorticity and psi again
     call apply_bc_case2(psi, vort, nx, ny, is_solid, n1, n2, nw, nt, nh)
     call vortbc_rect(psi, vort, nx, ny, n1, nw, nt)

     ! Simple progress display
     if (mod(iter, 5000) == 0) then
        write(*,'("iter=",i7,"  dpmx=",1pe12.5,"  domx=",1pe12.5)') iter, dpmx, domx
        call flush(6)
     end if

     !------------- psi/u/v/vort at ix=390 snapshots every 20000 iterations ---
     if (mod(iter, 20000) == 0) then
        do j = 2, ny-1
           y = real(j-1, dp)
           u = 0.5_dp * ( psi(ix_snap, j+1) - psi(ix_snap, j-1) )
           v = -0.5_dp * ( psi(ix_snap+1, j) - psi(ix_snap-1, j) )
           write(usnap,'(i9,1x,i6,1x,f10.4,4(1x,f14.6))') iter, j, y, u, v,          &
                psi(ix_snap,j), vort(ix_snap,j)
        end do
        write(usnap,*)
        call flush(usnap)
     end if

  end do

  write(*,*) 'Finished Case 2: final dpmx, domx = ', dpmx, domx
  call flush(6)

  !----------------------- output final field (x,y,u,v,psi,vort) --------------
  open(unit=uout, file='case2_field_200K_step10_h100.dat', status='replace', action='write')
  write(uout,'("# x  y  u  v  psi  vort")')

  do j = 2, ny-1, stepy
     y = real(j-1, dp)
     do i = 2, nx-1, stepx
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
  close(usnap)

  deallocate(psi, vort, is_solid)

contains

  !---------------------- Case 2 boundary conditions --------------------------
  !
  ! Upstream edge (i=1):
  !   - zero vorticity
  !   - no x-derivative of psi  => psi(1,j) = psi(2,j)
  !
  ! Downstream edge (i=nx):
  !   - no x-derivative for both psi and vort:
  !       psi(nx,j)  = psi(nx-1,j)
  !       vort(nx,j) = vort(nx-1,j)
  !
  ! Top edge (j=ny):
  !   - dψ/dy = 1  => psi(i,ny) = psi(i,ny-1) + 1
  !   - zero vorticity
  !
  ! Bottom edges: same centerline condition as case 1 (psi=0, vort=0).
  !
  subroutine apply_bc_case2(psi, vort, nx, ny, is_solid, n1, n2, nw, nt, nh)
    integer,  intent(in)    :: nx, ny, n1, n2, nw, nt, nh
    real(dp), intent(inout) :: psi(nx,ny), vort(nx,ny)
    logical,  intent(in)    :: is_solid(nx,ny)
    integer :: i, j

    ! Upstream edge (i=1): zero vorticity + no x-derivative of psi
    do j = 1, ny
       psi(1,j)  = psi(2,j)        ! dψ/dx = 0
       vort(1,j) = 0.0_dp          ! zero vorticity
    end do

    ! Downstream edge (i=nx): no x-derivative for both psi and vort
    do j = 1, ny
       psi(nx,j)  = psi(nx-1,j)    ! dψ/dx = 0
       vort(nx,j) = vort(nx-1,j)   ! dω/dx = 0
    end do

    ! Top edge (j=ny): dψ/dy = 1, zero vorticity
    do i = 1, nx
       psi(i,ny)  = psi(i,ny-1) + 1.0_dp   ! Δy = 1 => dψ/dy = 1
       vort(i,ny) = 0.0_dp                 ! zero vorticity
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

    ! Ensure psi is fixed inside the solid block
    do j = 1, ny
       do i = 1, nx
          if (is_solid(i,j)) then
             psi(i,j) = 0.0_dp
          end if
       end do
    end do
  end subroutine apply_bc_case2

  !-------------------- vorticity BCs at object walls (like vortbc) ----------
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

end program rect_obstruction_case2

