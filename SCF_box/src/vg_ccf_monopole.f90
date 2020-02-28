program vg_ccf_monopole
  implicit none
  
  integer, parameter:: dp=kind(0.d0)
  
  real(dp) :: rgrid, boxsize, vol, rhomed
  real(dp) :: posx, posy, posz, disx, disy, disz, dis
  real(dp) :: xvc, yvc, zvc, rv, min_rv, max_rv, median_rv
  real(dp) :: rwidth, rmax, rmin
  real(dp) :: pi = 4.*atan(1.)
  
  integer*4 :: ng, nc, nrbin, rind
  integer*4 :: id, iargc
  integer*4 :: i, j, ii, jj, ix, iy, iz, ix2, iy2, iz2
  integer*4 :: indx, indy, indz, nrows, ncols
  integer*4 :: ipx, ipy, ipz, ndif
  integer*4 :: ngrid
  
  integer*4, dimension(:, :, :), allocatable :: lirst, nlirst
  integer*4, dimension(:), allocatable :: ll
  
  real(dp), dimension(3) :: r
  real(dp), allocatable, dimension(:,:)  :: pos_data
  real(dp), dimension(:), allocatable :: rbin, rbin_edges
  real(dp), dimension(:), allocatable :: VG, VR, xi
  
  character(20), external :: str
  character(len=500) :: input_tracers, input_centres, output_den
  character(len=10) :: rmax_char, rmin_char, nrbin_char, ngrid_char
  character(len=100) :: min_rv_char, max_rv_char, box_char
  character(len=1)  :: creturn = achar(13)
  
  if (iargc() .ne. 10) then
      write(*,*) 'Some arguments are missing.'
      write(*,*) '1) input_data'
      write(*,*) '2) input_centres'
      write(*,*) '3) output_den'
      write(*,*) '4) boxsize'
      write(*,*) '5) rmin'
      write(*,*) '6) rmax'
      write(*,*) '7) nrbin'
      write(*,*) '8) min_rv'
      write(*,*) '9) max_rv'
      write(*,*) '10) ngrid'
      write(*,*) ''
      stop
    end if
    
  call getarg(1, input_tracers)
  call getarg(2, input_centres)
  call getarg(3, output_den)
  call getarg(4, box_char)
  call getarg(5, rmin_char)
  call getarg(6, rmax_char)
  call getarg(7, nrbin_char)
  call getarg(8, min_rv_char)
  call getarg(9, max_rv_char)
  call getarg(10, ngrid_char)
  
  read(box_char, *) boxsize
  read(rmin_char, *) rmin
  read(rmax_char, *) rmax
  read(nrbin_char, *) nrbin
  read(min_rv_char, *) min_rv
  read(max_rv_char, *) max_rv
  read(ngrid_char, *) ngrid
  
  write(*,*) '-----------------------'
  write(*,*) 'Running vg_ccf_monopole.exe'
  write(*,*) 'input parameters:'
  write(*,*) ''
  write(*, *) 'input_tracers: ', trim(input_tracers)
  write(*, *) 'input_centres: ', trim(input_centres)
  write(*, *) 'boxsize: ', trim(box_char)
  write(*, *) 'output_den: ', trim(output_den)
  write(*, *) 'rmin: ', trim(rmin_char), ' Mpc'
  write(*, *) 'rmax: ', trim(rmax_char), ' Mpc'
  write(*, *) 'nrbin: ', trim(nrbin_char)
  write(*, *) 'min_rv: ', trim(min_rv_char), ' Mpc'
  write(*, *) 'max_rv: ', trim(max_rv_char), ' Mpc'
  write(*, *) 'ngrid: ', trim(ngrid_char)
  write(*,*) ''

  open(10, file=input_tracers, status='old', form='unformatted')
  read(10) nrows
  read(10) ncols
  allocate(pos_data(ncols, nrows))
  read(10) pos_data
  close(10)
  ng = nrows
  if (id == 0) write(*,*) 'ntracers: ', ng

  nc = 0
  open(11, file=input_centres, status='old')
  do
    read(11, *, end=11)
    nc = nc + 1
  end do
  11 rewind(11)
  write(*,*) 'Number of voids: ', nc
  
  allocate(rbin(nrbin))
  allocate(rbin_edges(nrbin+1))
  allocate(VG(nrbin))
  allocate(VR(nrbin))
  allocate(xi(nrbin))
  
  rwidth = (rmax - rmin) / nrbin
  do i = 1, nrbin + 1
    rbin_edges(i) = rmin+(i-1)*rwidth
  end do
  do i = 1, nrbin
    rbin(i) = rbin_edges(i+1)-rwidth/2.
  end do
  
  ! Mean density inside the box
  rhomed = ng / (boxsize ** 3)
  
  ! Construct linked list for tracers
  write(*,*) ''
  write(*,*) 'Constructing linked list...'
  allocate(lirst(ngrid, ngrid, ngrid))
  allocate(nlirst(ngrid, ngrid, ngrid))
  allocate(ll(ng))
  rgrid = (boxsize) / real(ngrid)
  
  lirst = 0
  ll = 0
  
  do i = 1, ng
    indx = int((pos_data(1, i)) / rgrid + 1.)
    indy = int((pos_data(2, i)) / rgrid + 1.)
    indz = int((pos_data(3, i)) / rgrid + 1.)
  
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
    indz.gt.0.and.indz.le.ngrid)lirst(indx,indy,indz)=i
  
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
    indz.gt.0.and.indz.le.ngrid)nlirst(indx,indy,indz) = &
    nlirst(indx, indy, indz) + 1
  end do
  
  do i = 1, ng
    indx = int((pos_data(1, i))/ rgrid + 1.)
    indy = int((pos_data(2, i))/ rgrid + 1.)
    indz = int((pos_data(3, i))/ rgrid + 1.)
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
    &indz.gt.0.and.indz.le.ngrid) then
      ll(lirst(indx,indy,indz)) = i
      lirst(indx,indy,indz) = i
    endif
  end do
  
  write(*,*) 'Linked list successfully constructed'
  write(*,*) ''
  write(*,*) 'Starting loop over voids...'
  
  VG = 0
  VR = 0
  
  do i = 1, nc
    read (11, *) xvc, yvc, zvc, rv
  
    if (rv .lt. min_rv .or. rv .gt. max_rv) cycle
  
    ipx = int((xvc) / rgrid + 1.)
    ipy = int((yvc) / rgrid + 1.)
    ipz = int((zvc) / rgrid + 1.)
  
    !ndif = int((rmax * rv / rgrid + 1.))
    ndif = int(rmax / rgrid + 1.)
  
    do ix = ipx - ndif, ipx + ndif
      do iy = ipy - ndif, ipy + ndif
        do iz = ipz - ndif, ipz + ndif
  
          ix2 = ix
          iy2 = iy
          iz2 = iz
  
          if (ix2 .gt. ngrid) ix2 = ix2 - ngrid
          if (ix2 .lt. 1) ix2 = ix2 + ngrid
          if (iy2 .gt. ngrid) iy2 = iy2 - ngrid
          if (iy2 .lt. 1) iy2 = iy2 + ngrid
          if (iz2 .gt. ngrid) iz2 = iz2 - ngrid
          if (iz2 .lt. 1) iz2 = iz2 + ngrid
  
          ii = lirst(ix2,iy2,iz2)
          if(ii.ne.0) then
            do
              ii = ll(ii)
              disx = pos_data(1, ii) - xvc
              disy = pos_data(2, ii) - yvc
              disz = pos_data(3, ii) - zvc
  
              if (disx .lt. -boxsize/2) disx = disx + boxsize
              if (disx .gt. boxsize/2) disx = disx - boxsize
              if (disy .lt. -boxsize/2) disy = disy + boxsize
              if (disy .gt. boxsize/2) disy = disy - boxsize
              if (disz .lt. -boxsize/2) disz = disz + boxsize
              if (disz .gt. boxsize/2) disz = disz - boxsize
  
              r = (/ disx, disy, disz /)
              dis = norm2(r)
  
              if (dis .lt. rmax) then
                rind = int((dis - rmin) / rwidth + 1)
                VG(rind) = VG(rind) + 1
              end if
  
              if(ii.eq.lirst(ix2,iy2,iz2)) exit
  
            end do
          end if
        end do
      end do
    end do
  
    do ii = 1, nrbin
      vol = 4./3 * pi * ((rbin(ii) + rwidth/2.) ** 3 - &
      & (rbin(ii) - rwidth/2.) ** 3)

      VR(ii) = VR(ii) + rhomed * vol
    end do
  
  end do
  
  write(*,*) ''
  write(*,*) 'Calculation finished. Writing output...'
  
  xi = (VG * 1./VR) - 1
  
  open(12, file=output_den, status='unknown')
  do i = 1, nrbin
    write(12, fmt='(4f10.5)') rbin(i), rbin_edges(i), rbin_edges(i + 1), xi(i)
  end do
  
  stop
  
  end program vg_ccf_monopole
  