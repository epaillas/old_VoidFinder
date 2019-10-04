program overlapping_2D
implicit none

integer :: i, j, jj, counter, ix, iy, ix2, iy2
integer :: nv, nr, ngrid, ndif
integer :: indx, indy, ipx, ipy

integer, dimension(:, :), allocatable :: lirst
integer, dimension(:), allocatable :: ll

real :: overlap, boxsize, rgrid
real :: xread, yread, xv, yv, rread, ngread, ndenread
real :: disx, disy, dis

integer, dimension(:), allocatable :: ngv, mark

real, dimension(:,:), allocatable :: pos_voids, pos_rand
real, dimension(:), allocatable :: rv, nden

character(len=500) :: input_voids, output_voids
character(len=10) :: overlap_char, ngrid_char, box_char
character(len=1)  :: creturn = achar(13)

logical :: condition

if (iargc() .ne. 5) then
  write(*,*) 'Some arguments are missing.'
  write(*,*) '1) input_voids'
  write(*,*) '2) output_voids'
  write(*,*) '3) boxsize'
  write(*,*) '4) overlap'
  write(*,*) '5) ngrid'
  write(*,*) ''
  stop
end if

call getarg(1, input_voids)
call getarg(2, output_voids)
call getarg(3, box_char)
call getarg(4, overlap_char)
call getarg(5, ngrid_char)

read(box_char, *) boxsize
read(overlap_char, *) overlap
read(ngrid_char, *) ngrid

write(*,*) '-----------------------'
write(*,*) 'Running overlapping_2D.exe'
write(*,*) 'Input parameters:'
write(*,*) ''
write(*,*) 'input_voids: ', trim(input_voids)
write(*,*) 'output_voids: ', trim(output_voids)
write(*,*) 'boxsize: ', box_char, ' Mpc/h'
write(*,*) 'overlap: ', overlap_char, ' Mpc/h'
write(*,*) 'ngrid: ', ngrid_char

overlap = 1 - overlap

nv = 0
open(10, file=input_voids, status='old')
do
  read(10, *, end=10)
  nv = nv + 1
end do
10 rewind(10)
write(*, *) 'n_voids: ', nv

allocate(pos_voids(2, nv), rv(nv), ngv(nv), mark(nv), nden(nv))

do i = 1, nv
  read(10, *) xread, yread, rread, ngread, ndenread
  pos_voids(1, i) = xread
  pos_voids(2, i) = yread
  rv(i) = rread
  nden(i) = ndenread
  ngv(i) = ngread
end do
close(10)

mark = 0
rgrid = boxsize / ngrid
write(*, *) 'min_rvoid: ', minval(rv)
write(*, *) 'max_rvoid: ', maxval(rv)
write(*,*) 'rgrid: ', rgrid, ' Mpc'
write(*,*) ''

allocate(ll(nv))
allocate(lirst(ngrid, ngrid))

lirst = 0
ll = 0

do i = 1, nv
  indx = int(pos_voids(1, i) / rgrid + 1.)
  indy = int(pos_voids(2, i) / rgrid + 1.)

  if (indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid) then
    lirst(indx, indy) = i
  end if
end do

do i = 1, nv
  indx = int(pos_voids(1, i) / rgrid + 1.)
  indy = int(pos_voids(2, i) / rgrid + 1.)

  if (indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid) then
    ll(lirst(indx, indy)) = i
    lirst(indx, indy) = i
  end if
end do

do i = 1, nv - 1
  if (mod(i, int(1e3)) .eq. 1) then
    write(*, 101, advance='no' ) creturn , i , nv
    101 format(a, ' Void ', i8 ,' out of  ', i8)
  end if

  if (mark(i) .ne. 1) then

    xv = pos_voids(1, i)
    yv = pos_voids(2, i)

    ndif = int(overlap * (rv(i) + rv(i + 1)) / rgrid + 1)

    ipx = int(xv / rgrid + 1)
    ipy = int(yv / rgrid + 1)

    do ix = ipx - ndif, ipx + ndif
      do iy = ipy - ndif, ipy + ndif

        ix2 = ix
        iy2 = iy

        if (ix2 .gt. ngrid) ix2 = ix2 - ngrid
        if (ix2 .lt. 1) ix2 = ix2 + ngrid
        if (iy2 .gt. ngrid) iy2 = iy2 - ngrid
        if (iy2 .lt. 1) iy2 = iy2 + ngrid

        j = lirst(ix2, iy2)
        if (j .ne. 0) then
          do
            j = ll(j)
            disx = pos_voids(1, j) - xv
            disy = pos_voids(2, j) - yv

            if (disx .lt. -boxsize/2) disx = disx + boxsize
            if (disx .gt. boxsize/2) disx = disx - boxsize
            if (disy .lt. -boxsize/2) disy = disy + boxsize
            if (disy .gt. boxsize/2) disy = disy - boxsize

            dis = sqrt(disx ** 2 + disy ** 2)

            condition = rv(j) .le. rv(i) .and. dis .lt. overlap * (rv(i) + rv(j))
            if (condition .and. i .ne. j) mark(j) = 1

            if (j .eq. lirst(ix2, iy2)) exit

          end do
        end if
      end do
    end do
  end if
end do

counter = 0
open(12, file=output_voids, status='unknown')

do i = 1, nv
  if(mark(i) .eq. 0) then
    counter = counter + 1
    write(12, '(3F10.3, 1I10, 1F10.3)') &
    & pos_voids(1,i), pos_voids(2,i), rv(i), ngv(i), nden(i)
  endif
end do

write(*,*) ''
write(*,*) 'Final number of voids: ', counter

end program overlapping_2D
