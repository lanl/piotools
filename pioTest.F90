program test
  use pio_interface
  use iso_c_binding
  implicit none
!  use pio_interface
  character(len=1024) :: fname
  real(c_double) :: t0
  integer(c_int64_t) :: ncell, iStart, iCount, i, j
  integer :: ID, ndim, idim, nMat, imat
  real(c_double), pointer, dimension(:) :: vcell, center1, center2
  type(pio_2d_t) :: fvol
  
  
  if (command_argument_count() /= 1) then
     call get_command_argument(0,fname)
     stop 'Usage: ' // trim(fname) // ' <pio_file>'
  end if
  call get_command_argument(1,fname)
  write(*,*) 'fname = ', trim(fname)

  write(*,*) "______________________"
  ID = 1
  write(*,*) "FROM FORTRAN"
  call pio_init(1, fname, 0)
  if(pio_exists(id,'vcell', 0)) then
     t0 = pio_now()
     vcell => pio_get_d(id, 'vcell', 0)
     t0 = pio_now() - t0
  end if
  ncell = pio_nCell(id)
  ndim = pio_nDim(id)
  iStart = ncell/2+1
  iCount = 10
  iCount = min(10000, ncell - iStart);
  
  !   center1 => pio_get_range_d(ID, "cell_center", 1, iStart, iCount)
  !   center2 => pio_get_range_d(ID, "cell_center", 2, iStart, iCount)
  center1 => pio_get_d(ID, "cell_center", 1)
  center2 => pio_get_d(ID, "cell_center", 2)
  do i=1,min(10,iCount)
     write(*,*) i + iStart-1, ":", center1(i), " ", center2(i)
  end do

  nMat = pio_nmat(id)
  t0 = pio_now()
  fvol = pio_get_range_matvar(ID, "chunk_vol", iStart, iCount)
  t0 = pio_now() - t0

  write(*,*) "_____________Volume chunks 1-10: dt=", t0
  nMat = size(fvol%data, 1)
  do i = 1,min(10,iCount)
     !write(*,"(i8,' ',2(G20.12,1x), 2(f20.8,1x))") i+iStart-1, &
     !(fvol%data(j)%p(i), j=1,nMat), center1(i+iStart-1), center2(i+iStart-1)
     write(*,*) i+iStart-1, (fvol%data(j)%p(i), j=1,nMat), center1(i+iStart-1), center2(i+iStart-1)
  end do

  call pio_release(fvol)
  call pio_release(center1);
  call pio_release(center2);
  call pio_release(id);
end program test

