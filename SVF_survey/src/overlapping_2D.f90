module procedures
  implicit none
contains

  subroutine linked_list(ngrid, rgrid, gridmin, gridmax, ll, lirst, pos)
    implicit none
    integer :: i, ng, ipx, ipy
    integer, intent(in) :: ngrid
    real, intent(in) :: gridmin, gridmax, rgrid
    real, dimension(:,:), intent(in) :: pos
    integer, dimension(:,:), intent(out) :: lirst
    integer, dimension(:), intent(out) :: ll

    ng = size(pos, dim=2)
    lirst = 0
    ll = 0
    do i = 1, ng
      ipx = int((pos(1, i) - gridmin) / rgrid + 1.)
      ipy = int((pos(2, i) - gridmin) / rgrid + 1.)
      if(ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid) then
         lirst(ipx, ipy) = i
      end if
    end do

    do i = 1, ng
      ipx = int((pos(1, i) - gridmin) / rgrid + 1.)
      ipy = int((pos(2, i) - gridmin) / rgrid + 1.)

      if (ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid) then
        ll(lirst(ipx, ipy)) = i
        lirst(ipx, ipy) = i
      endif
    end do

  end subroutine linked_list

  character(len=20) function str(k)
    implicit none
    integer, intent(in) :: k
    write(str, *) k
    str = adjustl(str)
  end function str

end module procedures


program overlapping_2D
  use procedures
  implicit none

  integer :: i, j, count, nv, nr
  integer :: ix, iy, ipx, ipy, ndif
  integer, parameter :: ngrid = 100

  integer, dimension(ngrid, ngrid) :: lirst
  integer, dimension(:), allocatable :: ll

  real :: overlap, px, py, rvoid, ng, nden
  real :: disx, disy, dis, rgrid, gridmin, gridmax

  integer, dimension(:), allocatable :: ng_arr, mark

  real, dimension(:,:), allocatable :: pos_arr
  real, dimension(:), allocatable :: rvoid_arr, nden_arr

  character(len=500) :: input_voids, output_voids
  character(len=10) :: overlap_char, gridmin_char, gridmax_char
  character(len=1)  :: creturn = achar(13)

  logical :: cond

  if (iargc() .ne. 5) then
    write(*,*) 'Some arguments are missing.'
    write(*,*) '1) input_voids'
    write(*,*) '2) output_voids'
    write(*,*) '3) overlap'
    write(*,*) '4) gridmin'
    write(*,*) '5) gridmax'
    write(*,*) ''
    stop
  end if

  call getarg(1, input_voids)
  call getarg(2, output_voids)
  call getarg(3, overlap_char)
  call getarg(4, gridmin_char)
  call getarg(5, gridmax_char)

  read(overlap_char, *) overlap
  read(gridmin_char, *) gridmin
  read(gridmax_char, *) gridmax

  write(*,*) '-----------------------'
  write(*,*) 'Running overlapping_2D.exe'
  write(*,*) 'Input parameters:'
  write(*,*) ''
  write(*,*) 'input_voids: ', trim(input_voids)
  write(*,*) 'output_voids: ', trim(output_voids)
  write(*,*) 'overlap: ', overlap_char
  write(*,*) 'gridmin: ', gridmin_char, ' Mpc'
  write(*,*) 'gridmax: ', gridmax_char, ' Mpc'
  write(*,*) ''

  overlap = 1 - overlap

  nv = 0
  open(10, file=input_voids, status='old')
  do
    read(10, *, end=10)
    nv = nv + 1
  end do
  10 rewind(10)
  write(*, *) 'Number of voids: ', nv

  allocate(pos_arr(2, nv), rvoid_arr(nv), ng_arr(nv), nden_arr(nv), mark(nv))

  do i = 1, nv
    read(10, *) px, py, rvoid, ng, nden
    pos_arr(1, i) = px
    pos_arr(2, i) = py
    rvoid_arr(i) = rvoid
    nden_arr(i) = nden
    ng_arr(i) = ng
  end do
  close(10)

  mark = 0
  rgrid = (gridmax - gridmin) / ngrid

  write(*,*) 'min_rvoid:', maxval(rvoid_arr)
  write(*,*) 'max_rvoid:', minval(rvoid_arr)
  write(*,*) 'rgrid: ', rgrid, ' Mpc'
  write(*,*) ''

  allocate(ll(nv))
  call linked_list(ngrid, rgrid, gridmin, gridmax, ll, lirst, pos_arr)

  do i = 1, nv - 1
    if (mod(i, int(1e3)) .eq. 1) then
      write(*,*) 'Void ', i, ' out of ', nv
    end if

    if (mark(i) .ne. 1) then

      px = pos_arr(1, i)
      py = pos_arr(2, i)

      ndif = int(overlap * (rvoid_arr(i) + rvoid_arr(i + 1)) / rgrid + 1)

      ipx = int((px - gridmin) / rgrid + 1)
      ipy = int((py - gridmin) / rgrid + 1)

      do ix = ipx - ndif, ipx + ndif
        do iy = ipy - ndif, ipy + ndif

          j = lirst(ix, iy)
          if (j .ne. 0) then
            do
              j = ll(j)
              disx = pos_arr(1, j) - px
              disy = pos_arr(2, j) - py

              dis = sqrt(disx ** 2 + disy ** 2)

              cond = rvoid_arr(j) .le. rvoid_arr(i) .and. dis .lt.&
              & overlap * (rvoid_arr(i) + rvoid_arr(j))

              if (cond .and. i .ne. j) mark(j) = 1

              if (j .eq. lirst(ix, iy)) exit

            end do
          end if
        end do
      end do
    end if
  end do

  count = 0
  open(11, file=output_voids, status='unknown')

  do i = 1, nv
    if(mark(i) .eq. 0) then
      count = count + 1
      write(11, '(3F10.3, 1I10, 1F10.3)') &
      & pos_arr(1,i), pos_arr(2,i), rvoid_arr(i), ng_arr(i), nden_arr(i)
    endif
  end do

  write(*,*) ''
  write(*,*) 'Final number of voids: ', count

end program overlapping_2D
