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
  subroutine relax_vort_sweep(psi, vort, nx, ny, wrel, R, ixst, ixfn, iyst, iyfn, dvmx)
    integer, intent(in)    :: nx, ny
    integer, intent(inout) :: ixst, ixfn, iyst, iyfn
    real(dp), intent(in)   :: psi(nx,ny), wrel, R
    real(dp), intent(inout):: vort(nx,ny)
    real(dp), intent(out)  :: dvmx
    integer :: i, j, incx, incy
    real(dp) :: oldv, updt, delta, SCO, SSO, Fterm

    dvmx = 0.0_dp

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
        if (delta > dvmx) dvmx = delta
      end do
    end do
  end subroutine relax_vort_sweep

end module relax_cfd
!===============================================================================
program rect_obstruction_case1
  use relax_cfd
  implicit none

  !---------------- main parameters & IO ---------------------------------------
  integer :: nmax, ndec
  integer :: nx, ny, n1, n2, nw, nt, nh
  real(dp) :: w, r, eps
  character(len=70) :: outflow, outpsi

  ! arrays
  real(dp), allocatable :: psi(:,:), vort(:,:)

  ! iteration control
  integer :: iter
  real(dp) :: dpmx, dp_loc, dv_loc

  ! region indices
  integer :: ixst, ixfn, iyst, iyfn
  integer :: ix, iy

  ! velocity & coordinates for output
  real(dp) :: u, v, x, y

  !---------------------------------------------------------------------------
  ! Input: same structure as your professor's code
  !   1) nmax, ndec
  !   2) nx, ny, n1, n2, nw, nt, nh
  !   3) w, r, eps
  !   4) outflow, outpsi
  !---------------------------------------------------------------------------
  read(*,*) nmax, ndec
  read(*,*) nx, ny, n1, n2, nw, nt, nh
  read(*,*) w, r, eps
  read(*,*) outflow, outpsi

  allocate(psi(nx,ny), vort(nx,ny))

  open(25, file='flow.txt', status='replace', action='write')
  open(35, file='psi.txt',  status='replace', action='write')

  !----------------- initial condition (same logic as initial) ----------------
  call initial_case1(psi, vort, nx, ny, n1, nw, nt)

  !----------------- main iteration loop --------------------------------------
  do iter = 1, nmax
     dpmx = 0.0_dp

     !=================== Forward relaxation ==========================
     ! Region I: Left of the obstruction
     ixst = 2
     ixfn = n1
     iyst = 2
     iyfn = nt+1
     call relax_psi_sweep (psi, vort, nx, ny, w, ixst, ixfn, iyst, iyfn, dp_loc)
     call vortbc_case1(psi, vort, nx, ny, n1, nw, nt)
     call relax_vort_sweep(psi, vort, nx, ny, w, r, ixst, ixfn, iyst, iyfn, dv_loc)
     if (dp_loc > dpmx) dpmx = dp_loc

     ! Region II: Right of the obstruction
     ixst = n1+nw+2
     ixfn = n1+n2+nw
     iyst = 2
     iyfn = nt+1
     call relax_psi_sweep (psi, vort, nx, ny, w, ixst, ixfn, iyst, iyfn, dp_loc)
     call vortbc_case1(psi, vort, nx, ny, n1, nw, nt)
     call relax_vort_sweep(psi, vort, nx, ny, w, r, ixst, ixfn, iyst, iyfn, dv_loc)
     if (dp_loc > dpmx) dpmx = dp_loc

     ! Region III: Above the obstruction
     ixst = 2
     ixfn = n1+n2+nw
     iyst = nt+2
     iyfn = nt+nh
     call relax_psi_sweep (psi, vort, nx, ny, w, ixst, ixfn, iyst, iyfn, dp_loc)
     call vortbc_case1(psi, vort, nx, ny, n1, nw, nt)
     call relax_vort_sweep(psi, vort, nx, ny, w, r, ixst, ixfn, iyst, iyfn, dv_loc)
     if (dp_loc > dpmx) dpmx = dp_loc

     !=================== Backward relaxation =========================
     ! Region III: above, backward
     ixst = n1+n2+nw
     ixfn = 2
     iyst = nt+nh
     iyfn = nt+2
     call relax_psi_sweep (psi, vort, nx, ny, w, ixst, ixfn, iyst, iyfn, dp_loc)
     call vortbc_case1(psi, vort, nx, ny, n1, nw, nt)
     call relax_vort_sweep(psi, vort, nx, ny, w, r, ixst, ixfn, iyst, iyfn, dv_loc)
     if (dp_loc > dpmx) dpmx = dp_loc

     ! Region II: right, backward
     ixst = n1+n2+nw
     ixfn = n1+nw+2
     iyst = nt+1
     iyfn = 2
     call relax_psi_sweep (psi, vort, nx, ny, w, ixst, ixfn, iyst, iyfn, dp_loc)
     call vortbc_case1(psi, vort, nx, ny, n1, nw, nt)
     call relax_vort_sweep(psi, vort, nx, ny, w, r, ixst, ixfn, iyst, iyfn, dv_loc)
     if (dp_loc > dpmx) dpmx = dp_loc

     ! Region I: left, backward
     ixst = n1
     ixfn = 2
     iyst = nt+1
     iyfn = 2
     call relax_psi_sweep (psi, vort, nx, ny, w, ixst, ixfn, iyst, iyfn, dp_loc)
     call vortbc_case1(psi, vort, nx, ny, n1, nw, nt)
     call relax_vort_sweep(psi, vort, nx, ny, w, r, ixst, ixfn, iyst, iyfn, dv_loc)
     if (dp_loc > dpmx) dpmx = dp_loc

     !=================== convergence / progress ======================
     if (dpmx < eps) then
        write(*,*) 'Converged or stopped at iter = ', iter, ' dpmx = ', dpmx
        call flush(6)
        exit
     end if

     if (mod(iter, 1000) == 0) then
        write(*,'("iter=",i7,"  dpmx=",1pe12.5)') iter, dpmx
        call flush(6)
     end if

     !================= write line at ix=390 every ndec =================
     if (mod(iter, ndec) == 0) then
        ix = 390
        if (ix >= 2 .and. ix <= nx-1) then
           write(35,'("# iter=",i8)') iter
           do iy = 2, nt+nh
              if (iy >= 2 .and. iy <= ny-1) then
                 u = 0.5_dp*(psi(ix,iy+1) - psi(ix,iy-1))
                 v = -0.5_dp*(psi(ix+1,iy) - psi(ix-1,iy))
                 write(35,'(2I8,4F14.6)') iter, iy, u, v, psi(ix,iy), vort(ix,iy)
              end if
           end do
           call flush(35)
        end if
     end if

  end do
  !===================== end of iteration loop ================================

  write(*,*) 'Finished Case 1, final dpmx = ', dpmx
  call flush(6)

  !---------------------- final coarse output to outflow ----------------------
  write(25,'("# x   y   u   v   psi   vort")')

  do ix = 5, nx, 10
     do iy = 5, ny, 10
        if (.not.((ix > n1) .and. (ix < (n1+nw+1)) .and. (iy < (nt+1)))) then
           if (ix >= 2 .and. ix <= nx-1 .and. iy >= 2 .and. iy <= ny-1) then
              x = real(ix-1, dp)
              y = real(iy-1, dp)
              u = 0.5_dp * (psi(ix,iy+1) - psi(ix,iy-1))
              v = -0.5_dp * (psi(ix+1,iy) - psi(ix-1,iy))
              write(25,'(6F14.6)') x, y, u, v, psi(ix,iy), vort(ix,iy)
           end if
        end if
     end do
  end do
  call flush(25)

  close(25)
  close(35)
  deallocate(psi, vort)

contains
  !========================= initial (Case 1) ================================
  subroutine initial_case1(psi, vort, nx, ny, n1, nw, nt)
    integer,  intent(in)    :: nx, ny, n1, nw, nt
    real(dp), intent(inout) :: psi(nx,ny), vort(nx,ny)
    integer :: ix, iy
    real(dp) :: d2x, d2y, frac
    integer :: ixmn, ixmx, iymn, iymx

    ! upstream edge (i=1)
    ix = 1
    do iy = 1, ny
       psi(ix,iy)  = 1.0_dp*iy - 1.0_dp
       vort(ix,iy) = 0.0_dp
    end do

    ! downstream edge (i=nx)
    ix = nx
    do iy = 1, ny
       psi(ix,iy)  = 1.0_dp*iy - 1.0_dp
       vort(ix,iy) = 0.0_dp
    end do

    ! Top edge (j=ny)
    iy = ny
    do ix = 1, nx
       psi(ix,iy)  = 1.0_dp*iy - 1.0_dp
       vort(ix,iy) = 0.0_dp
    end do

    ! Bottom left edge
    iy = 1
    do ix = 1, n1+1
       psi(ix,iy)  = 0.0_dp
       vort(ix,iy) = 0.0_dp
    end do

    ! Bottom right edge
    iy = 1
    do ix = n1+nw+1, nx
       psi(ix,iy)  = 0.0_dp
       vort(ix,iy) = 0.0_dp
    end do

    ! Top edge of obstruction
    iy = nt+1
    do ix = n1+1, n1+nw+1
       psi(ix,iy) = 0.0_dp
    end do

    ! Left edge of obstruction
    ix = n1+1
    do iy = 1, nt+1
       psi(ix,iy) = 0.0_dp
    end do

    ! Right edge of obstruction
    ix = n1+nw+1
    do iy = 1, nt+1
       psi(ix,iy) = 0.0_dp
    end do

    ! Region I: Psi interior
    ixmn = 1
    ixmx = n1+1
    do ix = 2, n1
       do iy = 2, nt+1
          frac = real(ix-ixmn,dp) / real(ixmx-ixmn,dp)
          psi(ix,iy) = psi(ixmn,iy) + frac*(psi(ixmx,iy) - psi(ixmn,iy))
       end do
    end do

    ! Region II: Psi interior
    ixmn = n1+nw+1
    ixmx = nx
    do ix = ixmn+1, nx-1
       do iy = 2, nt+1
          frac = real(ix-ixmn,dp) / real(ixmx-ixmn,dp)
          psi(ix,iy) = psi(ixmn,iy) + frac*(psi(ixmx,iy) - psi(ixmn,iy))
       end do
    end do

    ! Region III: Psi interior
    iymn = nt+1
    iymx = ny
    do ix = 2, nx-1
       do iy = nt+2, ny-1
          frac = real(iy-iymn,dp) / real(iymx-iymn,dp)
          psi(ix,iy) = psi(ix,iymn) + frac*(psi(ix,iymx) - psi(ix,iymn))
       end do
    end do

    ! Vorticity on object walls
    ! left edge of object
    ix = n1+1
    do iy = 2, nt+1
       vort(ix,iy) = -2.0_dp*(psi(ix-1,iy) - psi(ix,iy))
    end do

    ! right edge of object
    ix = n1+nw+1
    do iy = 2, nt+1
       vort(ix,iy) = -2.0_dp*(psi(ix+1,iy) - psi(ix,iy))
    end do

    ! top edge of object
    iy = nt+1
    do ix = n1+1, n1+nw+1
       vort(ix,iy) = -2.0_dp*(psi(ix,iy+1) - psi(ix,iy))
    end do

    ! Vorticity interior - Region I
    do ix = 2, n1
       do iy = 2, nt+1
          d2x = psi(ix+1,iy) - 2.0_dp*psi(ix,iy) + psi(ix-1,iy)
          d2y = psi(ix,iy+1) - 2.0_dp*psi(ix,iy) + psi(ix,iy-1)
          vort(ix,iy) = -d2x - d2y
       end do
    end do

    ! Region II
    do ix = n1+nw+2, nx-1
       do iy = 2, nt+1
          d2x = psi(ix+1,iy) - 2.0_dp*psi(ix,iy) + psi(ix-1,iy)
          d2y = psi(ix,iy+1) - 2.0_dp*psi(ix,iy) + psi(ix,iy-1)
          vort(ix,iy) = -d2x - d2y
       end do
    end do

    ! Region III
    do ix = 2, nx-1
       do iy = nt+2, ny-1
          d2x = psi(ix+1,iy) - 2.0_dp*psi(ix,iy) + psi(ix-1,iy)
          d2y = psi(ix,iy+1) - 2.0_dp*psi(ix,iy) + psi(ix,iy-1)
          vort(ix,iy) = -d2x - d2y
       end do
    end do
  end subroutine initial_case1

  !========================= vorticity BC on object ===========================
  subroutine vortbc_case1(psi, vort, nx, ny, n1, nw, nt)
    integer,  intent(in)    :: nx, ny, n1, nw, nt
    real(dp), intent(inout) :: psi(nx,ny), vort(nx,ny)
    integer :: ix, iy

    ! Top side of the object
    iy = nt+1
    do ix = n1+1, n1+nw+1
       psi(ix,iy)  = 0.0_dp
       vort(ix,iy) = -2.0_dp*psi(ix,iy+1)
    end do

    ! Left side of the object
    ix = n1+1
    do iy = 1, nt+1
       psi(ix,iy)  = 0.0_dp
       vort(ix,iy) = -2.0_dp*psi(ix-1,iy)
    end do

    ! Right side of the object
    ix = n1+nw+1
    do iy = 1, nt+1
       psi(ix,iy)  = 0.0_dp
       vort(ix,iy) = -2.0_dp*psi(ix+1,iy)
    end do

  end subroutine vortbc_case1

end program rect_obstruction_case1

