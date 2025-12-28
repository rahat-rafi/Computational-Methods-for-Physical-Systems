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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program rect_obstruction_Re_sweep
  use relax_cfd
  implicit none

  !---------------- geometry parameters (same as Case 1) ----------------------
  integer, parameter :: n1 = 400   ! upstream length
  integer, parameter :: n2 = 400   ! downstream length
  integer, parameter :: nw = 50    ! block width
  integer, parameter :: nt = 50    ! block height from bottom
  integer, parameter :: nh = 400   ! height above block

  integer, parameter :: nx = n1 + nw + n2 + 1   ! = 851  (x: 1..851)
  integer, parameter :: ny = nt + nh + 1        ! = 451  (y: 1..451)

  real(dp), parameter :: wrel = 0.1_dp
  integer,  parameter :: niter = 200000        ! iterations per R

  ! "manageable subset" spacing for output field
  integer,  parameter :: stepx = 10, stepy = 10

  ! Reynolds numbers to investigate
  integer, parameter :: nR = 4
  integer, dimension(nR), parameter :: R_int_list = (/ 1, 3, 5, 7 /)
  real(dp) :: R

  !----------------------- field & mask arrays ---------------------------------
  real(dp), allocatable :: psi(:,:), vort(:,:)
  logical,  allocatable :: is_solid(:,:)

  ! indices and geometry
  integer :: i, j, iter, kR
  integer :: ix_block_left, ix_block_right, jy_block_top

  ! relaxation diagnostics
  real(dp) :: dpmx, domx, dtmp

  ! rectangles for fluid regions (I, II, III)
  integer :: ixA1, ixA2, ixB1, ixB2, ixC1, ixC2
  integer :: jyA1, jyA2, jyB1, jyB2, jyC1, jyC2

  ! for velocity output
  real(dp) :: u, v, x, y

  integer, parameter :: uout = 20  ! output unit
  character(len=80) :: fname

  !--------------------------- allocating arrays ---------------------------------
  allocate(psi(nx,ny), vort(nx,ny), is_solid(nx,ny))

  !------------------------ defining obstacle indices ----------------------------
  ix_block_left  = n1 + 1          ! left edge of block  (ix = 401)
  ix_block_right = n1 + nw + 1     ! right edge of block (ix = 451)
  jy_block_top   = nt + 1          ! top of block        (iy = 51)

  !-------------------- geometry: solid block mask (fixed) --------------------
  do j = 1, ny
     do i = 1, nx
        is_solid(i,j) = .false.
     end do
  end do

  do j = 1, jy_block_top
     do i = ix_block_left, ix_block_right
        is_solid(i,j) = .true.
     end do
  end do

  !------------------------- defining fluid rectangles (fixed) -------------------
  ! Region I: left of the obstruction (Region I)
  ixA1 = 2
  ixA2 = n1
  jyA1 = 2
  jyA2 = nt + 1

  ! Region II: right of the obstruction (Region II)
  ixB1 = n1 + nw + 2
  ixB2 = n1 + n2 + nw
  jyB1 = 2
  jyB2 = nt + 1

  ! Region III: above the obstruction (Region III)
  ixC1 = 2
  ixC2 = n1 + n2 + nw
  jyC1 = nt + 2
  jyC2 = nt + nh            ! = ny - 1

  !====================================================================
  !  Looping over Reynolds numbers: R = 1, 3, 5, 7
  !====================================================================
  do kR = 1, nR
     R = real(R_int_list(kR), dp)

     if (kR == 1) then
        !-------------------------- initial uniform flow ----------------------
        do j = 1, ny
           do i = 1, nx
              psi(i,j)  = real(j-1, dp)    ! uniform flow: psi ~ y
              vort(i,j) = 0.0_dp
           end do
        end do

        ! inside solid: psi = 0
        do j = 1, jy_block_top
           do i = ix_block_left, ix_block_right
              psi(i,j)  = 0.0_dp
              vort(i,j) = 0.0_dp
           end do
        end do
     else
        ! For kR > 1, we KEEP psi,vort from previous R as starting field.
        ! Only boundaries will be reset inside apply_bc_case1 / vortbc_rect.
        write(*,*) 'Starting new R from previous solution, R = ', R
     end if

     ! Applying Case 1 boundary conditions once before iterations
     call apply_bc_case1(psi, vort, nx, ny, is_solid, n1, n2, nw, nt, nh)
     call vortbc_rect(psi, vort, nx, ny, n1, nw, nt)

     write(*,*) '========================================================='
     write(*,'("Beginning iterations for R =", f6.2)') R
     write(*,*) '========================================================='

     !---------------------------- main iteration loop -----------------------
     do iter = 1, niter

        ! 1) Applying outer unimpeded BCs and enforce psi in solid
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

        ! Enforcing BCs for vorticity at top/inlet/outlet and object again
        call apply_bc_case1(psi, vort, nx, ny, is_solid, n1, n2, nw, nt, nh)
        call vortbc_rect(psi, vort, nx, ny, n1, nw, nt)

        ! Progress display
        if (mod(iter, 5000) == 0) then
           write(*,'("R=",f5.2,"  iter=",i7,"  dpmx=",1pe12.5,"  domx=",1pe12.5)') &
                R, iter, dpmx, domx
           call flush(6)
        end if

     end do

     write(*,'("Finished R=",f6.2,"  final dpmx=",1pe12.5,"  domx=",1pe12.5)') R, dpmx, domx

     !----------------------- output final field for this R -------------------
     write(fname,'("extra_field_R",i0,"_200K_step10.dat")') R_int_list(kR)
     open(unit=uout, file=fname, status='replace', action='write')
     write(uout,'("# R =",f6.2)') R
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

     write(*,'("Wrote field for R=",f6.2," to ",a)') R, trim(fname)
     call flush(6)

  end do  ! end loop over R

  deallocate(psi, vort, is_solid)

contains

  !---------------------- Case 1 boundary conditions --------------------------
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
             ! vorticity inside solid is not used; we leave it as-is or 0
          end if
       end do
    end do
  end subroutine apply_bc_case1

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

end program rect_obstruction_Re_sweep

