PROGRAM void_den_profile
implicit none

integer, parameter:: dp=kind(0.d0) ! Double precision

real(dp) :: gridmin, gridmax, rgrid, ncell
real(dp) :: posx, posy, posz, disx, disy, disz
real(dp) :: comx, comy, comz, pi, sigma
real(dp) :: xvc, yvc, zvc, rv, min_rvoid, max_rvoid
real(dp) :: swidth, piwidth, smax, smin, pimax, pimin

integer*4 :: ng, nr, nc, nsbin, npibin, sind, piind
integer*4 :: id, iargc
integer*4 :: i, j, ii, ix, iy, iz
integer*4 :: indx, indy, indz
integer*4 :: ipx, ipy, ipz, ndif
integer*4 :: cell_count
integer*4 :: ngrid

integer*4, dimension(:, :, :), allocatable :: lirst_data, nlirst_data
integer*4, dimension(:, :, :), allocatable :: lirst_rand, nlirst_rand
integer*4, dimension(:), allocatable :: ll_data, ll_rand

real(dp), dimension(3) :: com, r
real(dp), allocatable, dimension(:,:)  :: pos_data, pos_rand
real(dp), dimension(:), allocatable :: sbin, sbin_edges, pibin, pibin_edges
real(dp), dimension(:,:), allocatable :: VG, VR, xi

character(20), external :: str
character(len=500) :: input_tracers, input_randoms, input_centres, output_den
character(len=10) :: smax_char, smin_char, nsbin_char, gridmin_char, gridmax_char
character(len=100) :: min_rvoid_char, max_rvoid_char
character(len=1)  :: creturn = achar(13)

id = 0

if (iargc() .ne. 11) then
  if (id == 0) write(*,*) 'Some arguments are missing.'
  if (id == 0) write(*,*) '1) input_data'
  if (id == 0) write(*,*) '2) input_randoms'
  if (id == 0) write(*,*) '3) input_centres'
  if (id == 0) write(*,*) '4) output_den'
  if (id == 0) write(*,*) '5) smin'
  if (id == 0) write(*,*) '6) smax'
  if (id == 0) write(*,*) '7) nsbin'
  if (id == 0) write(*,*) '8) gridmin'
  if (id == 0) write(*,*) '9) gridmax'
  if (id == 0) write(*,*) '10) min_rvoid'
  if (id == 0) write(*,*) '11) max_rvoid'
  if (id == 0) write(*,*) ''
  stop
end if

call getarg(1, input_tracers)
call getarg(2, input_randoms)
call getarg(3, input_centres)
call getarg(4, output_den)
call getarg(5, smin_char)
call getarg(6, smax_char)
call getarg(7, nsbin_char)
call getarg(8, gridmin_char)
call getarg(9, gridmax_char)
call getarg(10, min_rvoid_char)
call getarg(11, max_rvoid_char)

read(gridmin_char, *) gridmin
read(gridmax_char, *) gridmax
read(smin_char, *) smin
read(smax_char, *) smax
read(nsbin_char, *) nsbin
read(min_rvoid_char, *) min_rvoid
read(max_rvoid_char, *) max_rvoid

if (id == 0) write(*,*) '-----------------------'
if (id == 0) write(*,*) 'Running void_gal_CCF.exe'
if (id == 0) write(*,*) 'Input parameters:'
if (id == 0) write(*,*) ''
if (id == 0) write(*, *) 'input_tracers: ', trim(input_tracers)
if (id == 0) write(*, *) 'input_randoms: ', trim(input_randoms)
if (id == 0) write(*, *) 'input_centres: ', trim(input_centres)
if (id == 0) write(*, *) 'output_den: ', trim(output_den)
if (id == 0) write(*, *) 'gridmin: ', trim(gridmin_char), ' Mpc'
if (id == 0) write(*, *) 'gridmax: ', trim(gridmax_char), ' Mpc'
if (id == 0) write(*, *) 'smin: ', trim(smin_char), ' Mpc'
if (id == 0) write(*, *) 'smax: ', trim(smax_char), ' Mpc'
if (id == 0) write(*, *) 'nsbin: ', trim(nsbin_char)
if (id == 0) write(*, *) 'min_rvoid: ', trim(min_rvoid_char), ' Mpc'
if (id == 0) write(*, *) 'max_rvoid: ', trim(max_rvoid_char), ' Mpc'
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

! Construct data catalogue linked list
ngrid = 100
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

npibin = nsbin
allocate(sbin(nsbin))
allocate(sbin_edges(nsbin+1))
allocate(pibin(npibin))
allocate(pibin_edges(npibin+1))
allocate(VG(nsbin, npibin))
allocate(VR(nsbin, npibin))
allocate(xi(nsbin, npibin))

swidth = (smax - smin) / nsbin
do i = 1, nsbin + 1
  sbin_edges(i) = smin+(i-1)*(smax-swidth)/(nsbin-1.)
end do
do i = 1, nsbin
  sbin(i) = sbin_edges(i+1)-swidth/2.
end do

pimin = smin
pimax = smax

piwidth = (pimax - pimin) / npibin
do i = 1, npibin + 1
  pibin_edges(i) = pimin+(i-1)*(pimax-piwidth)/(npibin-1.)
end do
do i = 1, npibin
  pibin(i) = pibin_edges(i+1)-piwidth/2.
end do

VG = 0
VR = 0
cell_count = 0

do i = 1, nc
  read(12,*) xvc, yvc, zvc, rv

  cell_count = cell_count + 1

  if (mod(cell_count, int(1e2)) .eq. 1) then
    write(*, 101, advance='no') creturn , cell_count , nc
    101 format(a , 'Centre ', i10 ,' out of', i10)
  end if

  if (rv .lt. min_rvoid .or. rv .gt. max_rvoid) cycle

  ipx = int((xvc - gridmin) / rgrid + 1.)
  ipy = int((yvc - gridmin) / rgrid + 1.)
  ipz = int((zvc - gridmin) / rgrid + 1.)

  !ndif = int(smax * rv / rgrid + 1)
  ndif = int(smax / rgrid + 1)

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

            comx = (xvc + pos_data(1, ii)) / 2
            comy = (yvc + pos_data(2, ii)) / 2
            comz = (zvc + pos_data(3, ii)) / 2

            r = (/ disx, disy, disz /)
            com = (/ comx, comy, comz /)

            pi = abs(dot_product(r, com)) / norm2(com)
            sigma = sqrt(norm2(r)**2 - pi**2)

            !pi = pi / rv
            !sigma = sigma / rv

            if (sigma .lt. smax .and. pi .lt. smax) then
              sind = int(sigma / swidth + 1)
              piind = int(pi / piwidth + 1)
              VG(sind, piind) = VG(sind, piind) + 1
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

            comx = (xvc + pos_rand(1, ii)) / 2
            comy = (yvc + pos_rand(2, ii)) / 2
            comz = (zvc + pos_rand(3, ii)) / 2

            r = (/ disx, disy, disz /)
            com = (/ comx, comy, comz /)

            pi = abs(dot_product(r, com)) / norm2(com)
            sigma = sqrt(norm2(r)**2 - pi**2)

            !pi = pi / rv
            !sigma = sigma / rv

            if (sigma .lt. smax .and. pi .lt. smax) then
              sind = int(sigma / swidth + 1)
              piind = int(pi / piwidth + 1)
              VR(sind, piind) = VR(sind, piind) + 1
            end if

            if (ii .eq. lirst_rand(ix, iy, iz)) exit
          end do
        end if
      end do
    end do
  end do
end do

xi = (nr * 1./ng) * (VG * 1./VR) - 1

open(13, file=output_den, status='unknown')
do i = 1, nsbin
  do j = 1, npibin
    write(13, fmt='(7f10.5)') sbin(i), sbin_edges(i), sbin_edges(i + 1),&
    & pibin(j), pibin_edges(j), pibin_edges(j + 1), xi(i, j)
  end do
end do

close(13)
deallocate(pos_data)
deallocate(pos_rand)

end PROGRAM void_den_profile

character(len=20) function str(k)
  implicit none
  ! Convert an integer to string.
  integer, intent(in) :: k
  write (str, *) k
  str = adjustl(str)

end function str
