program vg_ccf_monopole
implicit none

real*8 :: boxsize, box2
real*8 :: xvc, yvc, zvc, rvoid
real*8 :: gx, gy, gz
real*8 :: disx, disy, disz, rescaled_distance
real*8 :: dis
real*8 :: dbin, dmin, dmax
real*8 :: rhomed
real*8 :: rcell
real*8 :: pi = 4.*atan(1.)
real*8 :: rmin, rmax

real*8, allocatable, dimension (: , :) :: pos
real*8, dimension(:), allocatable :: weight, w
real*8, allocatable, dimension (: , :) :: DD, RR, dprofile
real*8, allocatable, dimension(:) :: bins, bin_edges
real*8, allocatable, dimension(:) :: r

integer*8 :: np
integer*8 :: nvc
integer*8 :: nbins
integer*8 :: ii, i, j, k, ix, iy, iz, ipx, ipy, ipz, ix2, iy2, iz2
integer*8 :: ndif
integer*8 :: ngridll, ncell
integer*8 :: indx, indy, indz
integer*8 :: bin_index
integer*8 :: nperx, npery, nperz
integer*8 :: neigh, counter

integer*8, dimension(:, :, :), allocatable :: nincell, ngrid, lirst, nlirst
integer*8, dimension(:), allocatable :: ll 

character(len=600) :: input_tracers, input_voids, output_file
character(len=10) :: box_char
character(len=10) :: nbins_char
character(len=10) :: dmin_char
character(len=10) :: dmax_char
character(len=100) :: fmt_char

! Get files
if (iargc() .ne. 7) then
  write(*,*) 'Wrong number of arguments.'
  write(*,*) 'Usage: ./vg_ccf_monopole.exe tracers_file voids_file output_file boxsize&
  & nbins dmin dmax'
  write(*,*) 'Exiting...'
  stop
end if

call getarg(1, input_tracers)
call getarg(2, input_voids)
call getarg(3, output_file)
call getarg(4, box_char)
call getarg(5, nbins_char)
call getarg(6, dmin_char)
call getarg(7, dmax_char)

read(box_char, *) boxsize
read(nbins_char, *) nbins
read(dmin_char, *) dmin
read(dmax_char, *) dmax

box2 = boxsize / 2.

write(*,*) 'Computing the 3D void galaxy density profile with the following.&
& parameters:'
write(*,*) ''
write(*,*) 'Tracers file: ', trim(input_tracers)
write(*,*) ''
write(*,*) 'Voids file: ', trim(input_voids)
write(*,*) ''
write(*,*) 'Output file: ', trim(output_file)
write(*,*) ''
write(*,*) 'Size of the simulation box: ', trim(box_char // ' Mpc/h')
write(*,*) 'Radial range: ', trim(dmin_char // dmax_char // ' void radii')
write(*,*) 'Number of bins: ', trim(nbins_char)
write(*,*) ''

! Count the number of tracers
open(10, file=input_tracers, status='old')
np = 0
do
  read(10, *, end=10)
  np = np + 1
end do
10 rewind(10)
write(*,*) 'Number of tracers: ', np

! Count the number of voids
nvc = 0 ! Total number of voids
open(11, file=input_voids, status='old')
do
  read(11, *, end=11)
  nvc = nvc + 1
end do
11 rewind(11)
write(*,*) 'Number of voids: ', nvc


! Memory allocation
allocate(bin_edges(nbins + 1))
allocate(bins(nbins))
allocate(DD(nvc, nbins))
allocate(RR(nvc, nbins))
allocate(r(nvc))

! Construct bins for the density profile
dbin = (dmax - dmin) / nbins

do i = 1, nbins + 1
  bin_edges(i) = dmin + (i - 1) * (dmax - dbin) / (nbins - 1.)
end do

do i = 1, nbins
  bins(i) = bin_edges(i + 1) - dbin / 2.
end do

! Store tracers positions in arrays
allocate(pos(3, np))
allocate(weight(np))

do i = 1, np
  read(10, *) gx, gy, gz
  pos(1, i) = gx
  pos(2, i) = gy
  pos(3, i) = gz
  weight(i) = 1.
end do
close(10)


! Mean density inside the box
!rhomed = sum(weight) / (boxsize ** 3)
rhomed = 1e7 / boxsize**3
print*, rhomed

! Construct linked list for tracers
write(*,*) ''
write(*,*) 'Constructing linked list...'
ngridll = 60
allocate(lirst(ngridll, ngridll, ngridll))
allocate(nlirst(ngridll, ngridll, ngridll))
allocate(ll(np))
rcell = (boxsize) / real(ngridll)

lirst = 0
ll = 0

do i = 1, np
  indx = int((pos(1, i)) / rcell + 1.)
  indy = int((pos(2, i)) / rcell + 1.)
  indz = int((pos(3, i)) / rcell + 1.)

  if(indx.gt.0.and.indx.le.ngridll.and.indy.gt.0.and.indy.le.ngridll.and.&
  indz.gt.0.and.indz.le.ngridll)lirst(indx,indy,indz)=i

  if(indx.gt.0.and.indx.le.ngridll.and.indy.gt.0.and.indy.le.ngridll.and.&
  indz.gt.0.and.indz.le.ngridll)nlirst(indx,indy,indz) = &
  nlirst(indx, indy, indz) + 1
end do

do i = 1, np
  indx = int((pos(1, i))/ rcell + 1.)
  indy = int((pos(2, i))/ rcell + 1.)
  indz = int((pos(3, i))/ rcell + 1.)
  if(indx.gt.0.and.indx.le.ngridll.and.indy.gt.0.and.indy.le.ngridll.and.&
  &indz.gt.0.and.indz.le.ngridll) then
    ll(lirst(indx,indy,indz)) = i
    lirst(indx,indy,indz) = i
  endif
end do

write(*,*) 'Linked list successfully constructed'
write(*,*) ''
write(*,*) 'Starting loop over voids...'

! Compute the density profile for each void
do i = 1, nvc ! For each void
  read (11, *) xvc, yvc, zvc, rvoid ! Read its data

  r(i) = rvoid

  if (mod(i, int(1e3)) .eq. 1) then
    write(*,*) 'Center', i, 'of', nvc
  end if

  ipx = int((xvc) / rcell + 1.)
  ipy = int((yvc) / rcell + 1.)
  ipz = int((zvc) / rcell + 1.)

  ndif = int((dmax * rvoid / rcell + 1.))

  neigh = 0

  do ix = ipx - ndif, ipx + ndif
    do iy = ipy - ndif, ipy + ndif
      do iz = ipz - ndif, ipz + ndif

        ix2 = ix
        iy2 = iy
        iz2 = iz

        if (ix2 .gt. ngridll) ix2 = ix2 - ngridll
        if (ix2 .lt. 1) ix2 = ix2 + ngridll
        if (iy2 .gt. ngridll) iy2 = iy2 - ngridll
        if (iy2 .lt. 1) iy2 = iy2 + ngridll
        if (iz2 .gt. ngridll) iz2 = iz2 - ngridll
        if (iz2 .lt. 1) iz2 = iz2 + ngridll

        ii = lirst(ix2,iy2,iz2)
        if(ii.ne.0) then
          do
            ii = ll(ii)
            neigh = neigh + 1
            disx = pos(1, ii) - xvc
            disy = pos(2, ii) - yvc
            disz = pos(3, ii) - zvc

            ! Periodic boundary conditions
            if (disx .lt. -box2) disx = disx + boxsize
            if (disx .gt. box2) disx = disx - boxsize
            if (disy .lt. -box2) disy = disy + boxsize
            if (disy .gt. box2) disy = disy - boxsize
            if (disz .lt. -box2) disz = disz + boxsize
            if (disz .gt. box2) disz = disz - boxsize

            dis = sqrt(disx ** 2 + disy ** 2 + disz ** 2)

            rescaled_distance = dis / rvoid

            if (rescaled_distance .le. dmax ) then

              bin_index = int((rescaled_distance - dmin) / dbin + 1)
              DD(i, bin_index) = DD(i, bin_index) + weight(ii)

            end if

            if(ii.eq.lirst(ix2,iy2,iz2)) exit

          end do
        end if
      end do
    end do
  end do

  ! Volume of a spherical shell at r
  do k = 1, nbins
    RR(i, k) =  4./3 * pi * ((rvoid * (bins(k) + dbin/2.)) ** 3 - &
    & (rvoid * (bins(k) - dbin/2.)) ** 3)
  end do

end do ! End loop over voids

write(*,*) ''
write(*,*) 'Calculation finished. Writing output...'

! Density profiles in units of the mean density
allocate(dprofile(nvc, nbins))
dprofile = (DD / RR) * 1. / rhomed

open(12, file=output_file, status='unknown')
fmt_char = trim(nbins_char // 'f10.3)')
write(12, fmt='(' // fmt_char) (bins(i), i = 1, nbins)

do i = 1, nvc
    write(12, fmt='(f10.3,' // fmt_char) r(i), (dprofile(i, j), j = 1, nbins)
end do

close(10)
close(11)
close(12)
close(13)

end program vg_ccf_monopole
