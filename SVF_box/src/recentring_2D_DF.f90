module procedures
  implicit none
contains

subroutine hpsort_eps_epw (n, ra, ind, eps)
  implicit none  
  integer, intent(in)   :: n  
  real*8, intent(in)  :: eps
  integer :: ind (n)  
  real*8 :: ra (n)
  integer :: i, ir, j, l, iind  
  real*8 :: rra  

  IF (ind (1) .eq.0) then  
     DO i = 1, n  
        ind (i) = i  
     ENDDO
  ENDIF
  IF (n.lt.2) return  
  l = n / 2 + 1  
  ir = n  

  sorting: do 
    IF ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
    ELSE  
       rra  = ra (ir)  
       iind = ind (ir)  
       ra (ir) = ra (1)  
       ind (ir) = ind (1)  
       ir = ir - 1  
       IF ( ir .eq. 1 ) then  
          ra (1)  = rra  
          ind (1) = iind  
          exit sorting  
       ENDIF
    ENDIF
    i = l  
    j = l + l  
    DO while ( j .le. ir )  
       IF ( j .lt. ir ) then  
          IF ( hslt( ra (j),  ra (j + 1) ) ) then  
             j = j + 1  
          ENDIF
       ENDIF
       IF ( hslt( rra, ra (j) ) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       ELSE
          j = ir + 1  
       ENDIF
    ENDDO
    ra (i) = rra  
    ind (i) = iind  

  END DO sorting    
contains 

  logical function hslt( a, b )
    REAL*8 :: a, b
    IF( abs(a-b) <  eps ) then
      hslt = .false.
    ELSE
      hslt = ( a < b )
    end if
  end function hslt

  !
end subroutine hpsort_eps_epw

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
  real*8 :: px, py, disx, disy, dis
  real*8 :: pxr, pyr, rvoidr, eps, quant
  real*8 :: rvoid, rwidth, rvoidmax
  real*8 :: rnd, rnd_phi, rnd_theta, rnd_rvoid
  real*8 :: rnd_px, rnd_py, rnd_ng, rnd_nden
  real*8 :: pi = 4.*atan(1.)

  integer :: ng, nc, nv, rind, stuck, nrows, ncols
  integer :: id, ierr, process_num, iargc, filenumber
  integer :: i, j, k, ii, ix, iy, ix2, iy2
  integer :: ipx, ipy, ndif, ngrid, ncells
  integer, parameter ::  nrbin = 1000, nrc = 128

  integer*4, dimension(:), allocatable :: ind

  real*8, dimension(:,:), allocatable :: field
  real*8, dimension(:), allocatable :: flat_field
  real*8, dimension(nrbin) :: rbin, cum_rbin, counter

  character(len=500) :: input_field, input_centres, output_voids
  character(len=10) :: box_char, rvoidmax_char, delta_char, ngrid_char
  character(len=1)  :: creturn = achar(13)

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, process_num, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, id, ierr)


  if (iargc() .ne. 6) then
    if (id == 0) write(*,*) 'recentering.exe: some parameters are missing.'
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
  if (id == 0) write(*,*) 'Running recentering.exe'
  if (id == 0) write(*,*) 'Input parameters:'
  if (id == 0) write(*,*) ''
  if (id == 0) write(*,*) 'mpi_processes: ', process_num
  if (id == 0) write(*,*) 'input_field: ', trim(input_field)
  if (id == 0) write(*,*) 'input_centres: ', trim(input_centres)
  if (id == 0) write(*,*) 'output_voids: ', trim(output_voids)
  if (id == 0) write(*,*) 'boxsize: ', trim(box_char), ' Mpc/h'
  if (id == 0) write(*,*) 'rvoidmax: ', trim(rvoidmax_char), ' Mpc/h'
  if (id == 0) write(*,*) 'density_threshold: ', trim(delta_char), ' * rho_mean'
  if (id == 0) write(*,*) 'random_centres: ', nrc
  if (id == 0) write(*,*) ''

  if (process_num .gt. 1) then
    output_voids = trim(output_voids) // '.' // trim(str(id))
  end if

  open(10, file=input_field, status='old', form='unformatted')
  read(10) nrows
  read(10) ncols
  allocate(field(ncols, nrows))
  read(10) field
  close(10)
  ngrid = nrows
  field = field + 1

  ! find 0.25 quantile of field
  allocate(flat_field(nrows*ncols))
  allocate(ind(nrows*ncols))
  eps = 1e-10
  flat_field = pack(field, .true.)
  call hpsort_eps_epw (nrows * ncols, flat_field, ind, eps)
  quant = flat_field((nrows*ncols)/10)

  open(11, file=input_centres, status='old')
  nc = 0
  do
    read(11, *, end=11)
    nc = nc + 1
  end do
  11 rewind(11)
  if (id == 0) write(*,*) 'quantile 0.1: ', quant
  if (id == 0) write(*,*) 'n_centres: ', nc
  if (id == 0) write(*,*) ''

  
  rho_mean = sum(field) / size(field)
  rgrid = boxsize / ngrid
  ndif = int(rvoidmax / rgrid + 1.)
  rwidth = rvoidmax / nrbin

  filenumber = id + 20
  open(filenumber, file=output_voids, status='unknown')
  nv = 0
  do i = 1, nc
    read(11,*) px, py, rvoid, nden

    pxr = px
    pyr = py
    rvoidr = rvoid

    if(mod(i, process_num) .eq. id) then

      print*, i, mod(i, process_num)

      stuck = 0
      do j = 1, nrc
        rbin = 0
        cum_rbin = 0
        counter = 0

        call random_number(rnd)
        rnd_phi = rnd * 2 * pi

        rnd_px = rvoid/4 * cos(rnd_phi) + px
        rnd_py = rvoid/4 * sin(rnd_phi) + py

        if (sqrt((rnd_px-pxr)**2 + (rnd_py-pyr)**2)&
        & .gt. rvoidr) cycle

        if (rnd_px .lt. 0) rnd_px = rnd_px + boxsize
        if (rnd_px .gt. boxsize) rnd_px = rnd_px - boxsize
        if (rnd_py .lt. 0) rnd_py = rnd_py + boxsize
        if (rnd_py .gt. boxsize) rnd_py = rnd_py - boxsize

        ipx = int(rnd_px / rgrid + 1.)
        ipy = int(rnd_py / rgrid + 1.)

        do ix = ipx - ndif, ipx + ndif, 1
          do iy = ipy - ndif, ipy + ndif, 1

            if (sqrt(real((ix-ipx)**2 +(iy-ipy)**2)).gt.ndif+1) cycle

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

        stuck = stuck + 1

        cum_rbin(1) = rbin(1)
        do ii = 2, nrbin
          cum_rbin(ii) =  cum_rbin(ii - 1) + rbin(ii)
        end do

        do ii = nrbin, 1, -1
          rnd_rvoid = rwidth * ii
          ncells = sum(counter(1:ii))
          rnd_nden = cum_rbin(ii) / ncells
          if (rnd_nden .lt. -0.01 .and. rnd_rvoid .gt. rvoid &
          .and. ncells .ge. 1) then
            rvoid = rnd_rvoid
            px = rnd_px
            py = rnd_py
            nden = rnd_nden / rho_mean
            stuck = 0
            exit
          end if
        end do

        if (stuck .gt. 64) exit

      end do

      write(filenumber, '(4F10.3)') px, py, rvoid, nden

    end if
  end do

  close(filenumber)

  if (id == 0) write(*,*) ''
  if (id == 0) write(*,*) ''

  call MPI_Finalize(ierr)

end PROGRAM recentering
