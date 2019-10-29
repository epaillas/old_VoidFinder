PROGRAM void_den_profile
implicit none

integer, parameter:: dp=kind(0.d0) ! Double precision

real(dp) :: thresh, gridmin, gridmax, rgrid
real(dp) :: posx, posy, posz, disx, disy, disz, dis
real(dp) :: minposx, minposy, minposz
real(dp) :: xvc, yvc, zvc, rv, rwidth, rmax, rmin
real(dp) :: tmp1, tmp2

integer*4 :: ng, nr, nv, nc, nrbin, rind, ng_loc, nr_loc
integer*4 :: id, ierr, process_num, iargc
integer*4 :: filenumber
integer*4 :: i, j, k, ii, ix, iy, iz
integer*4 :: indx, indy, indz
integer*4 :: ipx, ipy, ipz, ndif
integer*4 :: cell_count, ncell
integer*4 :: ngrid

integer*4, dimension(:, :, :), allocatable :: lirst_data, nlirst_data
integer*4, dimension(:, :, :), allocatable :: lirst_rand, nlirst_rand
integer*4, dimension(:), allocatable :: ll_data, ll_rand

real(dp), allocatable, dimension(:,:)  :: pos_data, pos_rand
real(dp), dimension(:), allocatable :: rbin_data, rbin_rand, rbin, rbin_edges, nden

character(20), external :: str
character(len=500) :: input_tracers, input_randoms, input_centres, output_den
character(len=10) :: rmax_char, rmin_char, nrbin_char, gridmin_char, gridmax_char
character(len=100) :: fmt_char
character(len=1)  :: creturn = achar(13)

id = 0

if (iargc() .ne. 9) then
  if (id == 0) write(*,*) 'Some arguments are missing.'
  if (id == 0) write(*,*) '1) input_data'
  if (id == 0) write(*,*) '2) input_randoms'
  if (id == 0) write(*,*) '3) input_centres'
  if (id == 0) write(*,*) '4) output_den'
  if (id == 0) write(*,*) '5) rmin'
  if (id == 0) write(*,*) '6) rmax'
  if (id == 0) write(*,*) '7) nrbin'
  if (id == 0) write(*,*) '8) gridmin'
  if (id == 0) write(*,*) '9) gridmax'
  if (id == 0) write(*,*) ''
  stop
end if

call getarg(1, input_tracers)
call getarg(2, input_randoms)
call getarg(3, input_centres)
call getarg(4, output_den)
call getarg(5, rmin_char)
call getarg(6, rmax_char)
call getarg(7, nrbin_char)
call getarg(8, gridmin_char)
call getarg(9, gridmax_char)

read(gridmin_char, *) gridmin
read(gridmax_char, *) gridmax
read(rmin_char, *) rmin
read(rmax_char, *) rmax
read(nrbin_char, *) nrbin

if (id == 0) write(*,*) '-----------------------'
if (id == 0) write(*,*) 'Running step0.exe'
if (id == 0) write(*,*) 'Input parameters:'
if (id == 0) write(*,*) ''
if (id == 0) write(*, *) 'mpi_processes: ', process_num
if (id == 0) write(*, *) 'input_tracers: ', trim(input_tracers)
if (id == 0) write(*, *) 'input_randoms: ', trim(input_randoms)
if (id == 0) write(*, *) 'input_centres: ', trim(input_centres)
if (id == 0) write(*, *) 'output_den: ', trim(output_den)
if (id == 0) write(*, *) 'gridmin: ', trim(gridmin_char), ' Mpc'
if (id == 0) write(*, *) 'gridmax: ', trim(gridmax_char), ' Mpc'
if (id == 0) write(*, *) 'rmin: ', trim(rmin_char), ' Mpc'
if (id == 0) write(*, *) 'rmax: ', trim(rmax_char), ' Mpc'
if (id == 0) write(*, *) 'nrbin: ', trim(nrbin_char)
if (id == 0) write(*,*) ''

open(10, file=input_tracers, status='old', form='unformatted')
read(10) ng
allocate(pos_data(3, ng))
read(10) pos_data
close(10)
if (id == 0) write(*,*) 'ntracers: ', ng

open(11, file=input_randoms, status='old', form='unformatted')
read(11) nr
allocate(pos_rand(3, nr))
read(11) pos_rand
close(11)
if (id == 0) write(*,*) 'nrandoms: ', nr

open(12, file=input_centres, status='old')
nc = 0
do
  read(12, *, end=12)
  nc = nc + 1
end do
12 rewind(12)
if (id == 0) write(*,*) 'Number of centres: ', nc

! Memory allocation
allocate(pos_data(3, ng))
allocate(pos_rand(3, nr))

! Read data catalogue positions
do i = 1, ng
  read(10, *) posx, posy, posz
  pos_data(1, i) = posx
  pos_data(2, i) = posy
  pos_data(3, i) = posz
end do
close(10)

! Read rand catalogue positions
do i = 1, nr
  read(11, *) posx, posy, posz
  pos_rand(1, i) = posx
  pos_rand(2, i) = posy
  pos_rand(3, i) = posz
end do
close(11)

! Construct data catalogue linked list
ngrid = 120
rgrid = (gridmax - gridmin) / ngrid

if (id == 0) write(*,*) 'rgrid: ', rgrid
if (id == 0) write(*,*) ''

allocate(ll_data(ng))
allocate(lirst_data(ngrid,ngrid,ngrid))
allocate(nlirst_data(ngrid,ngrid,ngrid))

lirst_data = 0
ll_data = 0

do i = 1, ng
  indx = int((pos_data(1, i) - gridmin) / rgrid + 1.)
  indy = int((pos_data(2, i) - gridmin) / rgrid + 1.)
  indz = int((pos_data(3, i) - gridmin) / rgrid + 1.)

  if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
  indz.gt.0.and.indz.le.ngrid) then
    lirst_data(indx, indy, indz) = i
    nlirst_data(indx, indy, indz) = nlirst_data(indx, indy, indz) + 1
  end if
end do

do i = 1, ng
  indx = int((pos_data(1, i) - gridmin) / rgrid + 1.)
  indy = int((pos_data(2, i) - gridmin) / rgrid + 1.)
  indz = int((pos_data(3, i) - gridmin) / rgrid + 1.)

  if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
  indz.gt.0.and.indz.le.ngrid) then
    ll_data(lirst_data(indx, indy, indz)) = i
    lirst_data(indx, indy, indz) = i
  end if
end do

! Construct random catalogue linked list
allocate(ll_rand(nr))
allocate(lirst_rand(ngrid,ngrid,ngrid))
allocate(nlirst_rand(ngrid,ngrid,ngrid))

lirst_rand = 0
ll_rand = 0

do i = 1, nr
  indx = int((pos_rand(1, i) - gridmin) / rgrid + 1.)
  indy = int((pos_rand(2, i) - gridmin) / rgrid + 1.)
  indz = int((pos_rand(3, i) - gridmin) / rgrid + 1.)

  if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
  indz.gt.0.and.indz.le.ngrid) then
    lirst_rand(indx, indy, indz) = i
    nlirst_rand(indx, indy, indz) =  nlirst_rand(indx, indy, indz) + 1
  end if
end do

do i = 1, nr
  indx = int((pos_rand(1, i) - gridmin) / rgrid + 1.)
  indy = int((pos_rand(2, i) - gridmin) / rgrid + 1.)
  indz = int((pos_rand(3, i) - gridmin) / rgrid + 1.)

  if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
  indz.gt.0.and.indz.le.ngrid) then
    ll_rand(lirst_rand(indx, indy, indz)) = i
    lirst_rand(indx, indy, indz) = i
  end if
end do

rwidth = rmax / nrbin
allocate(rbin_data(nrbin))
allocate(rbin_rand(nrbin))
allocate(nden(nrbin))
allocate(rbin(nrbin))
allocate(rbin_edges(nrbin+1))

rwidth = (rmax - rmin) / nrbin
do i = 1, nrbin + 1
  rbin_edges(i) = rmin+(i-1)*(rmax-rwidth)/(nrbin-1.)
end do
do i = 1, nrbin
  rbin(i) = rbin_edges(i+1)-rwidth/2.
end do

open(13, file=output_den, status='unknown')
fmt_char = trim(nrbin_char) // trim('f10.3)')
write(13, fmt='(' // fmt_char) (rbin(i), i = 1, nrbin)

cell_count = 0

! Loop over prospective void centres
do i = 1, nc
  read(12,*) xvc, yvc, zvc, rv

  cell_count = cell_count + 1

  if (mod(cell_count, int(1e2)) .eq. 1) then
    write(*, 101, advance='no') creturn , cell_count , nc
    101 format(a , 'Centre ', i10 ,' out of', i10)
  end if

  rbin_data = 0
  rbin_rand = 0

  ipx = int((xvc - gridmin) / rgrid + 1.)
  ipy = int((yvc - gridmin) / rgrid + 1.)
  ipz = int((zvc - gridmin) / rgrid + 1.)

  ndif = int(rmax * rv / rgrid + 1)

  do ix = ipx - ndif, ipx + ndif, 1
    do iy = ipy - ndif, ipy + ndif, 1
      do iz = ipz - ndif, ipz + ndif, 1

        ncell = sqrt(real((ix - ipx) ** 2 + (iy - ipy)** 2&
        & + (iz - ipz) ** 2))

        if(ncell .gt. ndif + 1) cycle

        ! Data catalogue
        ii = lirst_data(ix, iy, iz)
        if (ii .ne. 0) then
          do
            ii = ll_data(ii)

            disx = pos_data(1, ii) - xvc
            disy = pos_data(2, ii) - yvc
            disz = pos_data(3, ii) - zvc

            dis = sqrt(disx ** 2 + disy ** 2 + disz ** 2)
            dis = dis / rv

            if (dis .lt. rmax) then
              rind = int(dis / rwidth + 1)
              rbin_data(rind) = rbin_data(rind) + 1
            end if

            if (ii .eq. lirst_data(ix, iy, iz)) exit
          end do
        end if

        ! Random catalogue
        ii = lirst_rand(ix, iy, iz)
        if (ii .ne. 0) then
          do
            ii = ll_rand(ii)

            disx = pos_rand(1, ii) - xvc
            disy = pos_rand(2, ii) - yvc
            disz = pos_rand(3, ii) - zvc

            dis = sqrt(disx ** 2 + disy ** 2 + disz ** 2)
            dis = dis / rv

            if (dis .lt. rmax) then
              rind = int(dis / rwidth + 1)
              rbin_rand(rind) = rbin_rand(rind) + 1
            end if

            if (ii .eq. lirst_rand(ix, iy, iz)) exit
          end do
        end if
      end do
    end do
  end do

  ! Normalize random counts
  rbin_rand = rbin_rand  / (nr * 1./ng)

  ! Calculate density profile
  nden = rbin_data / rbin_rand
  write(13, fmt='(f10.3,' // fmt_char) rv, (nden(j), j = 1, nrbin)

end do

! Deallocate memory and end the program

close(13)
deallocate(pos_data)
deallocate(pos_rand)

if (id == 0) print*, ''
if (id == 0) write(*,*) '-----------------------'

end PROGRAM void_den_profile

character(len=20) function str(k)
  implicit none
  ! Convert an integer to string.
  integer, intent(in) :: k
  write (str, *) k
  str = adjustl(str)

end function str
