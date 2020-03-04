module procedures
  implicit none
contains

  subroutine linked_list(ngrid, rgrid, ll, lirst, pos_data)
    implicit none
    integer :: i, ng, ipx, ipy, ipz
    real*8, intent(in) :: rgrid
    integer, intent(in) :: ngrid
    real*8, dimension(:,:), intent(in) :: pos_data
    integer, dimension(:,:,:), intent(out) :: lirst
    integer, dimension(:), intent(out) :: ll

    ng = size(pos_data, dim=2)
    lirst = 0
    ll = 0
    do i = 1, ng
      ipx = int(pos_data(1, i) / rgrid + 1.)
      ipy = int(pos_data(2, i) / rgrid + 1.)
      ipz = int(pos_data(3, i) / rgrid + 1.)
      if(ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.&
      ipz.gt.0.and.ipz.le.ngrid) lirst(ipx, ipy, ipz) = i
    end do

    do i = 1, ng
      ipx = int(pos_data(1, i) / rgrid + 1.)
      ipy = int(pos_data(2, i) / rgrid + 1.)
      ipz = int(pos_data(3, i) / rgrid + 1.)

      if (ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.ipz&
      &.gt.0.and.ipz.le.ngrid) then
        ll(lirst(ipx, ipy, ipz)) = i
        lirst(ipx, ipy, ipz) = i
      endif
    end do

  end subroutine linked_list

  character(len=20) function str(k)
    implicit none
    integer*4, intent(in) :: k
    write(str, *) k
    str = adjustl(str)
  end function str


end module procedures


PROGRAM recentering
  use mpi
  use procedures
  implicit none

  real*8 :: boxsize, delta, rgrid, rho_mean, nden
  real*8 :: px, py, pz, disx, disy, disz, dis
  real*8 :: pxr, pyr, pzr, rvoidr
  real*8 :: rvoid, rwidth, rvoidmax
  real*8 :: rnd, rnd_phi, rnd_theta, rnd_rvoid
  real*8 :: rnd_px, rnd_py, rnd_pz, rnd_ng, rnd_nden
  real*8 :: pi = 4.*atan(1.)

  integer :: ng, nc, nv, rind, stuck, nrows, ncols
  integer :: id, ierr, process_num, iargc, filenumber
  integer :: i, j, k, ii, ix, iy, iz, ix2, iy2, iz2
  integer :: ipx, ipy, ipz, ndif, ngrid
  integer, parameter ::  nrbin = 1000, nrc = 128

  integer, dimension(:,:,:), allocatable :: lirst
  integer, dimension(:), allocatable :: ll

  real*8, allocatable, dimension(:,:)  :: pos_data
  real*8, dimension(nrbin) :: rbin, cum_rbin

  character(len=500) :: input_tracers, input_centres, output_voids
  character(len=10) :: box_char, rvoidmax_char, delta_char, ngrid_char
  character(len=1)  :: creturn = achar(13)

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, process_num, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, id, ierr)


  if (iargc() .ne. 7) then
    if (id == 0) write(*,*) 'recentering.exe: some parameters are missing.'
    if (id == 0) write(*,*) ''
    if (id == 0) write(*,*) '1) input_tracers'
    if (id == 0) write(*,*) '2) input_centres'
    if (id == 0) write(*,*) '3) output_voids'
    if (id == 0) write(*,*) '4) boxsize'
    if (id == 0) write(*,*) '5) density_threshold'
    if (id == 0) write(*,*) '6) rvoidmax'
    if (id == 0) write(*,*) '7) ngrid'
    stop
  end if

  call getarg(1, input_tracers)
  call getarg(2, input_centres)
  call getarg(3, output_voids)
  call getarg(4, box_char)
  call getarg(5, delta_char)
  call getarg(6, rvoidmax_char)
  call getarg(7, ngrid_char)

  read(box_char, *) boxsize
  read(rvoidmax_char, *) rvoidmax
  read(delta_char, *) delta
  read(ngrid_char, *) ngrid


  if (id == 0) write(*,*) '-----------------------'
  if (id == 0) write(*,*) 'Running recentering.exe'
  if (id == 0) write(*,*) 'Input parameters:'
  if (id == 0) write(*,*) ''
  if (id == 0) write(*,*) 'mpi_processes: ', process_num
  if (id == 0) write(*,*) 'input_tracers: ', trim(input_tracers)
  if (id == 0) write(*,*) 'input_centres: ', trim(input_centres)
  if (id == 0) write(*,*) 'output_voids: ', trim(output_voids)
  if (id == 0) write(*,*) 'boxsize: ', trim(box_char), ' Mpc/h'
  if (id == 0) write(*,*) 'rvoidmax: ', trim(rvoidmax_char), ' Mpc/h'
  if (id == 0) write(*,*) 'density_threshold: ', trim(delta_char), ' * rho_mean'
  if (id == 0) write(*,*) 'random_centres: ', nrc
  if (id == 0) write(*,*) 'ngrid: ', ngrid, ' Mpc/h'
  if (id == 0) write(*,*) ''

  if (process_num .gt. 1) then
    output_voids = trim(output_voids) // '.' // trim(str(id))
  end if

  open(10, file=input_tracers, status='old', form='unformatted')
  read(10) nrows
  read(10) ncols
  allocate(pos_data(ncols, nrows))
  read(10) pos_data
  close(10)
  ng = nrows
  if (id == 0) write(*,*) 'ntracers: ', ng

  open(11, file=input_centres, status='old')
  nc = 0
  do
    read(11, *, end=11)
    nc = nc + 1
  end do
  11 rewind(11)
  if (id == 0) write(*,*) 'n_centres: ', nc
  if (id == 0) write(*,*) ''

  
  rho_mean = ng * 1./(boxsize ** 3)
  rgrid = boxsize / ngrid
  ndif = int(rvoidmax / rgrid + 1.)
  rwidth = rvoidmax / nrbin

  allocate(ll(ng))
  allocate(lirst(ngrid, ngrid, ngrid))
  call linked_list(ngrid, rgrid, ll, lirst, pos_data)

  filenumber = id + 20
  open(filenumber, file=output_voids, status='unknown')
  nv = 0
  do i = 1, nc
    read(11,*) px, py, pz, rvoid, ng, nden

    pxr = px
    pyr = py
    pzr = pz
    rvoidr = rvoid

    ! if (id == 0 .and. mod(i, int(1e2)) .eq. 1) then
    !   write(*, 101, advance='no' ) creturn , i , nc
    !   101 format(a , 'Centre ', i10 ,' out of', i10)
    ! end if

    if(mod(i, process_num) .eq. id) then

      stuck = 0
      do j = 1, nrc
        rbin = 0
        cum_rbin = 0

        call random_number(rnd)
        rnd_phi = rnd * 2 * pi
        call random_number(rnd)
        rnd_theta = acos(rnd * 2 - 1)

        rnd_px = rvoid/4 * sin(rnd_theta) * cos(rnd_phi) + px
        rnd_py = rvoid/4 * sin(rnd_theta) * sin(rnd_phi) + py
        rnd_pz = rvoid/4 * cos(rnd_theta) + pz

        if (sqrt((rnd_px-pxr)**2 + (rnd_py-pyr)**2 + (rnd_pz-pzr)**2 )&
        & .gt. rvoidr) cycle

        if (rnd_px .lt. 0) rnd_px = rnd_px + boxsize
        if (rnd_px .gt. boxsize) rnd_px = rnd_px - boxsize
        if (rnd_py .lt. 0) rnd_py = rnd_py + boxsize
        if (rnd_py .gt. boxsize) rnd_py = rnd_py - boxsize
        if (rnd_pz .lt. 0) rnd_pz = rnd_pz + boxsize
        if (rnd_pz .gt. boxsize) rnd_pz = rnd_pz - boxsize

        ipx = int(rnd_px / rgrid + 1.)
        ipy = int(rnd_py / rgrid + 1.)
        ipz = int(rnd_pz / rgrid + 1.)

        do ix = ipx - ndif, ipx + ndif, 1
          do iy = ipy - ndif, ipy + ndif, 1
            do iz = ipz - ndif, ipz + ndif, 1

              if (sqrt(real((ix-ipx)**2 +(iy-ipy)**2+(iz-ipz)**2)).gt.ndif+1) cycle

              ix2 = ix
              iy2 = iy
              iz2 = iz

              if (ix2 .gt. ngrid) ix2 = ix2 - ngrid
              if (ix2 .lt. 1) ix2 = ix2 + ngrid
              if (iy2 .gt. ngrid) iy2 = iy2 - ngrid
              if (iy2 .lt. 1) iy2 = iy2 + ngrid
              if (iz2 .gt. ngrid) iz2 = iz2 - ngrid
              if (iz2 .lt. 1) iz2 = iz2 + ngrid

              ii = lirst(ix2, iy2, iz2)
              if (ii .ne. 0) then
                do
                  ii = ll(ii)

                  disx = pos_data(1, ii) - rnd_px
                  disy = pos_data(2, ii) - rnd_py
                  disz = pos_data(3, ii) - rnd_pz

                  if (disx .lt. -boxsize/2) disx = disx + boxsize
                  if (disx .gt. boxsize/2) disx = disx - boxsize
                  if (disy .lt. -boxsize/2) disy = disy + boxsize
                  if (disy .gt. boxsize/2) disy = disy - boxsize
                  if (disz .lt. -boxsize/2) disz = disz + boxsize
                  if (disz .gt. boxsize/2) disz = disz - boxsize

                  if (sqrt(disx ** 2 + disy ** 2 + disz ** 2)&
                  & .lt. rvoidmax) then
                    dis = sqrt(disx ** 2 + disy ** 2 + disz ** 2)
                    rind = int(dis / rwidth + 1)
                    rbin(rind) = rbin(rind) + 1
                  end if

                  if (ii .eq. lirst(ix2, iy2, iz2)) exit
                end do
              end if

            end do
          end do
        end do

        stuck = stuck + 1

        cum_rbin(1) = rbin(1)
        do ii = 2, nrbin
          cum_rbin(ii) =  cum_rbin(ii - 1) + rbin(ii)
        end do

        do ii = nrbin, 1, -1
          rnd_rvoid = rwidth * ii
          rnd_ng = int(cum_rbin(ii))
          rnd_nden = cum_rbin(ii) / (4.e0/3.e0 * pi * rnd_rvoid ** 3)
          if (rnd_nden .lt. delta * rho_mean .and. rnd_rvoid .gt. rvoid&
          & .and. rnd_ng .gt. 0) then
            rvoid = rnd_rvoid
            px = rnd_px
            py = rnd_py
            pz = rnd_pz
            ng = rnd_ng
            nden = rnd_nden / rho_mean
            stuck = 0
            exit
          end if
        end do

        if (stuck .gt. 64) exit

      end do

      write(filenumber, '(4F10.3, 1I10, 1F10.3)') px, py, pz, rvoid, ng, nden

    end if
  end do

  close(filenumber)
  deallocate(pos_data)

  if (id == 0) write(*,*) ''
  if (id == 0) write(*,*) ''

  call MPI_Finalize(ierr)

end PROGRAM recentering
