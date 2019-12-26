module procedures
  implicit none
contains

  subroutine linked_list(ngrid, rgrid, ll, lirst, pos_data)
    implicit none
    integer*4 :: i, ng, ipx, ipy
    integer*4, intent(in) :: ngrid
    real*8, intent(in) :: rgrid
    real*8, dimension(:,:), intent(in) :: pos_data
    integer*4, dimension(:,:), intent(out) :: lirst
    integer*4, dimension(:), intent(out) :: ll

    ng = size(pos_data, dim=2)
    lirst = 0
    ll = 0
    do i = 1, ng
      ipx = int(pos_data(1, i) / rgrid + 1.)
      ipy = int(pos_data(2, i) / rgrid + 1.)
      if(ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid) then
        lirst(ipx, ipy) = i
      end if
    end do

    do i = 1, ng
      ipx = int(pos_data(1, i) / rgrid + 1.)
      ipy = int(pos_data(2, i) / rgrid + 1.)

      if (ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid) then
        ll(lirst(ipx, ipy)) = i
        lirst(ipx, ipy) = i
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


PROGRAM grow_spheres
  use mpi
  use procedures
  implicit none

  real*8 :: boxsize, delta, rgrid, rho_mean, nden
  real*8 :: px, py, disx, disy, dis
  real*8 :: rvoid, rwidth, rvoidmax
  real*8 :: pi = 4.*atan(1.)

  integer*4 :: ng, nc, nv, rind, nrows, ncols
  integer*4 :: id, ierr, process_num, iargc, filenumber
  integer*4 :: i, j, k, ii, ix, iy, ix2, iy2
  integer*4 :: ipx, ipy, ndif, ngrid
  integer*4, parameter :: nrbin = 1000

  real*8, allocatable, dimension(:,:)  :: field, centres
  real*8, dimension(nrbin) :: rbin, cum_rbin, counter

  character(len=500) :: input_field, input_centres, output_voids
  character(len=10) :: box_char, rvoidmax_char, delta_char, ngrid_char
  character(len=1)  :: creturn = achar(13)

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, process_num, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, id, ierr)


  if (iargc() .ne. 6) then
    if (id == 0) write(*,*) 'grow_spheres.exe: some parameters are missing.'
    if (id == 0) write(*,*) ''
    if (id == 0) write(*,*) '1) input_field'
    if (id == 0) write(*,*) '2) input_centres'
    if (id == 0) write(*,*) '3) output_voids'
    if (id == 0) write(*,*) '4) boxsize'
    if (id == 0) write(*,*) '5) density_threshold'
    if (id == 0) write(*,*) '6) rvoidmax'
    stop
  end if

  call getarg(1, input_field)
  call getarg(2, input_centres)
  call getarg(3, output_voids)
  call getarg(4, box_char)
  call getarg(5, delta_char)
  call getarg(6, rvoidmax_char)

  read(box_char, *) boxsize
  read(rvoidmax_char, *) rvoidmax
  read(delta_char, *) delta


  if (id == 0) write(*,*) '-----------------------'
  if (id == 0) write(*,*) 'Running grow_spheres_DF.exe'
  if (id == 0) write(*,*) 'Input parameters:'
  if (id == 0) write(*,*) ''
  if (id == 0) write(*,*) 'mpi_processes: ', process_num
  if (id == 0) write(*,*) 'input_field: ', trim(input_field)
  if (id == 0) write(*,*) 'input_centres: ', trim(input_centres)
  if (id == 0) write(*,*) 'output_voids: ', trim(output_voids)
  if (id == 0) write(*,*) 'boxsize: ', trim(box_char), ' Mpc/h'
  if (id == 0) write(*,*) 'rvoidmax: ', trim(rvoidmax_char), ' Mpc/h'
  if (id == 0) write(*,*) 'density_threshold: ', trim(delta_char), ' * rho_mean'

  if (process_num .gt. 1) then
    output_voids = trim(output_voids) // '.' // trim(str(id))
  end if

  open(10, file=input_field, status='old', form='unformatted')
  read(10) nrows
  read(10) ncols
  allocate(field(ncols, nrows))
  read(10) field
  ngrid = nrows
  field = field + 1

  open(11, file=input_centres, status='old', form='unformatted')
  read(11) nrows
  ncols = 2 ! need to rerun centres
  allocate(centres(ncols, nrows))
  read(11) centres
  nc = nrows

  if (id == 0) write(*,*) 'ncentres: ', nc
  if (id == 0) write(*,*) 'xmin, xmax: ', minval(centres(1,:)), maxval(centres(1,:))

  rho_mean = sum(field) / size(field)
  rgrid = boxsize / ngrid
  ndif = int(rvoidmax / rgrid + 1)
  rwidth = rvoidmax / nrbin

  if (id == 0) write(*,*) 'rho_mean: ', rho_mean
  if (id == 0) write(*,*) 'rgrid: ', rgrid, ' Mpc/h'
  if (id == 0) write(*,*) 'rwidth: ', rwidth, ' Mpc/h'
  if (id == 0) write(*,*) 'ndif: ', ndif
  if (id == 0) write(*,*) ''

  filenumber = id + 20
  open(filenumber, file=output_voids, status='unknown')
  nv = 0
  do i = 1, nc

    px = centres(1, i)
    py = centres(2, i)
    counter = 0


    if(mod(i, process_num) .eq. id) then

        rbin = 0
        cum_rbin = 0

        ipx = int(px / rgrid + 1)
        ipy = int(py / rgrid + 1)

        do ix = ipx - ndif, ipx + ndif, 1
          do iy = ipy - ndif, ipy + ndif, 1

            if (sqrt(real((ix - ipx) ** 2 + (iy - ipy)** 2))&
            & .gt. ndif + 1) cycle

            ix2 = ix
            iy2 = iy

            if (ix2 .gt. ngrid) ix2 = ix2 - ngrid
            if (ix2 .lt. 1) ix2 = ix2 + ngrid
            if (iy2 .gt. ngrid) iy2 = iy2 - ngrid
            if (iy2 .lt. 1) iy2 = iy2 + ngrid

            disx = (ix2 - 0.5) * rgrid - px
            disy = (iy2 - 0.5) * rgrid - py

            if (disx .lt. -boxsize/2) disx = disx + boxsize
            if (disx .gt. boxsize/2) disx = disx - boxsize
            if (disy .lt. -boxsize/2) disy = disy + boxsize
            if (disy .gt. boxsize/2) disy = disy - boxsize
            
            dis = sqrt(disx ** 2 + disy ** 2)

            if (dis .lt. rvoidmax) then
              rind = int(dis / rwidth + 1)
              rbin(rind) = rbin(rind) + field(ix2, iy2)
              counter(rind) = counter(rind) + 1
            end if

          end do
        end do

        cum_rbin(1) = rbin(1)
        do ii = 2, nrbin
          cum_rbin(ii) =  cum_rbin(ii - 1) + rbin(ii)
        end do

        do ii = nrbin, 1, -1
          rvoid = rwidth * ii
          nden = cum_rbin(ii) / sum(counter(1:ii))
          if (nden .lt. delta * rho_mean) then
            nv = nv + 1
            write(filenumber, '(4F10.3)') &
            px, py, rvoid, nden / rho_mean
            exit
          end if
        end do
    end if
  end do

  close(filenumber)

  if (id == 0) write(*,*) ''
  if (id == 0) write(*,*) ''
  if (id == 0) write(*,*) nv, ' voids found by rank 0.'

  call MPI_Finalize(ierr)

end PROGRAM grow_spheres
