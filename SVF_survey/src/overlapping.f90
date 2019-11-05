module procedures
  implicit none
contains

  subroutine linked_list(ngrid, rgrid, gridmin, gridmax, ll, lirst, pos)
    implicit none
    integer :: i, ng, ipx, ipy, ipz
    integer, intent(in) :: ngrid
    real*8, intent(in) :: gridmin, gridmax, rgrid
    real*8, dimension(:,:), intent(in) :: pos
    integer, dimension(:,:,:), intent(out) :: lirst
    integer, dimension(:), intent(out) :: ll

    ng = size(pos, dim=2)
    lirst = 0
    ll = 0
    do i = 1, ng
      ipx = int((pos(1, i) - gridmin) / rgrid + 1.)
      ipy = int((pos(2, i) - gridmin) / rgrid + 1.)
      ipz = int((pos(3, i) - gridmin) / rgrid + 1.)
      if(ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.&
      ipz.gt.0.and.ipz.le.ngrid) lirst(ipx, ipy, ipz) = i
    end do

    do i = 1, ng
      ipx = int((pos(1, i) - gridmin) / rgrid + 1.)
      ipy = int((pos(2, i) - gridmin) / rgrid + 1.)
      ipz = int((pos(3, i) - gridmin) / rgrid + 1.)

      if (ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.ipz&
      &.gt.0.and.ipz.le.ngrid) then
        ll(lirst(ipx, ipy, ipz)) = i
        lirst(ipx, ipy, ipz) = i
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


program overlapping
  use procedures
  implicit none

  integer :: i, j, count, nv, nr
  integer :: ix, iy, iz, ipx, ipy, ipz, ndif
  integer, parameter :: ngrid = 100

  integer, dimension(ngrid, ngrid, ngrid) :: lirst
  integer, dimension(:), allocatable :: ll

  real*8 :: overlap, px, py, pz, rvoid, ng, nden
  real*8 :: disx, disy, disz, dis, rgrid, gridmin, gridmax

  integer, dimension(:), allocatable :: ng_arr, mark

  real*8, dimension(:,:), allocatable :: pos_arr
  real*8, dimension(:), allocatable :: rvoid_arr, nden_arr

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
  write(*,*) 'Running overlapping.exe'
  write(*,*) 'Input parameters:'
  write(*,*) ''
  write(*,*) 'input_voids: ', trim(input_voids)
  write(*,*) 'output_voids: ', trim(output_voids)
  write(*,*) 'overlap: ', overlap_char, ' Mpc/h'
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

  allocate(pos_arr(3, nv), rvoid_arr(nv), ng_arr(nv),&
  & nden_arr(nv), mark(nv))


  do i = 1, nv
    read(10, *) px, py, pz, rvoid, ng, nden
    pos_arr(1, i) = px
    pos_arr(2, i) = py
    pos_arr(3, i) = pz
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
      pz = pos_arr(3, i)

      ndif = int(overlap * (rvoid_arr(i) + rvoid_arr(i + 1)) / rgrid + 1)

      ipx = int((px - gridmin) / rgrid + 1)
      ipy = int((py - gridmin) / rgrid + 1)
      ipz = int((pz - gridmin) / rgrid + 1)

      do ix = ipx - ndif, ipx + ndif
        do iy = ipy - ndif, ipy + ndif
          do iz = ipz - ndif, ipz + ndif

            j = lirst(ix, iy, iz)
            if (j .ne. 0) then
              do
                j = ll(j)
                disx = pos_arr(1, j) - px
                disy = pos_arr(2, j) - py
                disz = pos_arr(3, j) - pz

                dis = sqrt(disx ** 2 + disy ** 2 + disz ** 2)

                cond = rvoid_arr(j) .le. rvoid_arr(i) .and. dis .lt.&
                & overlap * (rvoid_arr(i) + rvoid_arr(j))

                if (cond .and. i .ne. j) mark(j) = 1

                if (j .eq. lirst(ix, iy, iz)) exit

              end do
            end if
          end do
        end do
      end do
    end if
  end do

  count = 0
  open(11, file=output_voids, status='unknown')

  do i = 1, nv
    if(mark(i) .eq. 0) then
      count = count + 1
      write(11, '(4F10.3, 1I10, 2F10.3)') &
      & pos_arr(1,i), pos_arr(2,i), pos_arr(3,i), rvoid_arr(i), ng_arr(i),&
      & nden_arr(i)
    endif
  end do

  write(*,*) ''
  write(*,*) 'Final number of voids: ', count

end program overlapping
