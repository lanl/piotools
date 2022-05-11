! ========================================================================================
!  (C) (or copyright) 2022. Triad National Security, LLC. All rights reserved.
! 
!  This program was produced under U.S. Government contract 89233218CNA000001
!  for Los Alamos National Laboratory (LANL), which is operated by Triad
!  National Security, LLC for the U.S. Department of Energy/National Nuclear
!  Security Administration. All rights in the program are reserved by Triad
!  National Security, LLC, and the U.S. Department of Energy/National Nuclear
!  Security Administration. The Government is granted for itself and others
!  acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license
!  in this material to reproduce, prepare derivative works, distribute copies to
!  the public, perform publicly and display publicly, and to permit others to do
!  so.
! ========================================================================================

!
! F90 functions defined in pioInterface()
! Version: 220509: initial definition

#define C0(x) (trim(x)//char(0))

module pio_interface
  use iso_c_binding

  implicit none

  public
  interface
     subroutine pio_init_c(ID, fname, verbose) BIND(C, name="pio_init")
       use iso_c_binding
       implicit none
       integer(c_int), VALUE, intent(in) :: ID
       character(c_char), intent(in) :: fname(1)
       integer(c_int), VALUE, intent(in) :: verbose
     end subroutine pio_init_c

     subroutine pio_release(ID) BIND(C, name="pio_release")
       use iso_c_binding
       implicit none
       integer(c_int), VALUE, intent(in) :: id
     end subroutine pio_release

     subroutine pio_release_d(fptr) BIND(C, name="pio_release_d")
       use iso_c_binding
       implicit none
       integer(c_int), VALUE, intent(in) :: id
       double(c_double), pointer, dimension(:), intent(inout) :: fptr
     end subroutine pio_release_d

     subroutine pio_release_i64(fptr) BIND(C, name="pio_release_i64")
       use iso_c_binding
       implicit none
       integer(c_int), VALUE, intent(in) :: id
       integer(c_int64_t), pointer, dimension(:), intent(inout) :: fptr
     end subroutine pio_release_i64

     integer(c_int64_t) function pio_ncell(ID) BIND(C, name="pio_nCell")
       use iso_c_binding
       implicit none
       integer(c_int), VALUE, intent(in) :: id
     end function pio_nCell

     integer(c_int) function pio_ndim(ID) BIND(C, name="pio_nDim")
       use iso_c_binding
       implicit none
       integer(c_int), VALUE, intent(in) :: id
     end function pio_ndim

     integer(c_int) function pio_nmat(ID) BIND(C, name="pio_nMat")
       use iso_c_binding
       implicit none
       integer(c_int), VALUE, intent(in) :: id
     end function pio_nmat

     logical(c_bool) function pio_exists_c(ID, var) BIND(C, name="pio_exists")
       use iso_c_binding
       implicit none
       integer(c_int), VALUE, intent(in) :: id
       character(c_char), intent(in) :: var(1)
     end function pio_exists_c

     integer(c_int64_t) function pio_width_c(ID, var) BIND(C, name="pio_width")
       use iso_c_binding
       implicit none
       integer(c_int), VALUE, intent(in) :: id
       character(c_char), intent(in) :: var(1)
     end function pio_width_c

     integer(c_int) function pio_length_c(ID, var) BIND(C, name="pio_length")
       use iso_c_binding
       implicit none
       integer(c_int), VALUE, intent(in) :: id
       character(c_char), intent(in) :: var(1)
     end function pio_length_c

     function pio_get_d_c(ID, var, index) BIND(C, name="pio_get_d")
       use iso_c_binding
       implicit none
       type(c_ptr) :: pio_get_d_c
       integer(c_int), VALUE, intent(in) :: id
       character(c_char), intent(in) :: var(1)
       integer(c_int), VALUE, intent(in) :: index
     end function pio_get_d_c
       
     function pio_get_range_d_c(ID, var, index, iStart, nCount) BIND(C, name="pio_get_d")
       use iso_c_binding
       implicit none
       type(c_ptr) :: pio_get_range_d_c
       integer(c_int), VALUE, intent(in) :: id
       character(c_char), intent(in) :: var(1)
       integer(c_int), VALUE, intent(in) :: index
       integer(c_int64_t), VALUE, intent(in) :: iStart
       integer(c_int64_t), VALUE, intent(in) :: nCount
     end function pio_get_range_d_c

     function pio_get_i64_c(ID, var, index) BIND(C, name="pio_get_d")
       use iso_c_binding
       implicit none
       type(c_ptr) :: pio_get_i64_c
       integer(c_int), VALUE, intent(in) :: id
       character(c_char), intent(in) :: var(1)
       integer(c_int), VALUE, intent(in) :: index
     end function pio_get_i64_c
       
     function pio_get_range_i64_c(ID, var, index, iStart, nCount) BIND(C, name="pio_get_d")
       use iso_c_binding
       implicit none
       type(c_ptr) :: pio_get_range_i64_c
       integer(c_int), VALUE, intent(in) :: id
       character(c_char), intent(in) :: var(1)
       integer(c_int), VALUE, intent(in) :: index
       integer(c_int64_t), VALUE, intent(in) :: iStart
       integer(c_int64_t), VALUE, intent(in) :: nCount
     end function pio_get_range_i64_c
       
    end interface

contains


  subroutine pio_init(ID, fname, verbose)
    use iso_c_binding
    implicit none
    integer,  intent(in) :: ID
    character*(*), intent(in) :: fname
    integer, intent(in) :: verbose
    call pio_init_c(ID, C0(fname), verbose)
  end subroutine pio_init

  integer(c_int64_t) function pio_width(ID, var)
    use iso_c_binding
    implicit none
    integer(c_int), VALUE, intent(in) :: id
    character*(*), intent(in) :: var
    pio_width = pio_width_c(ID, C0(var))
  end function pio_width

  integer(c_int64_t)  function pio_length(ID, var)
    use iso_c_binding
    implicit none
    integer(c_int), VALUE, intent(in) :: id
    character*(*), intent(in) :: var
    pio_length = pio_length_c(ID, C0(var))
  end function pio_length

  function pio_center(ID, index) result(values) 
    use iso_c_binding
    implicit none
    integer(c_int64_t), pointer, dimension (:) :: values
    interface
       function pio_center_c(ID, index) BIND(C, name="pio_center")
         use iso_c_binding
         implicit none
         type(c_ptr) :: pio_center_c
         integer(c_int), VALUE, intent(in) :: id
         integer(c_int), VALUE, intent(in) :: index
       end function pio_center_c
    end interface
    integer, intent(in) :: ID
    integer, intent(in) :: index
    integer(c_int64_t) :: n
    type(c_ptr) :: cptr
    n = pio_nCell(ID)
    cptr = pio_center_c(ID, index)
    call c_f_pointer(cptr, values, [n])
  end function pio_center

  function pio_daughter(ID) result(values) 
    use iso_c_binding
    implicit none
    integer(c_int64_t), pointer, dimension (:) :: values
    interface
       function pio_daughter_c(ID) BIND(C, name="pio_daughter")
         use iso_c_binding
         implicit none
         type(c_ptr) :: pio_daughter_c
         integer(c_int), VALUE, intent(in) :: id
       end function pio_daughter_c
    end interface
    integer, intent(in) :: ID
    integer(c_int64_t) :: n
    type(c_ptr) :: cptr
    n = pio_nCell(ID)
    cptr = pio_daughter_c(ID)
    call c_f_pointer(cptr, values, [n])
  end function pio_daughter

  logical(c_bool) function pio_exists(ID, var)
    implicit none
    integer, intent(in) :: id
    character*(*), intent(in) :: var
    pio_exists = pio_exists_c(id, C0(var))
  end function pio_exists

  function pio_get_d(ID, var, index) result(values)
    use iso_c_binding
    implicit none
    real(c_double), pointer, dimension (:) :: values
    integer, intent(in) :: id
    character*(*), intent(in) :: var
    integer, intent(in) :: index
    integer*8 :: n
    type(c_ptr) :: cptr
    n = pio_length(ID, C0(var))
    cptr = pio_get_d_c(ID, C0(var), index)
    call c_f_pointer(cptr, values, [n])
  end function pio_get_d

  function pio_get_range_d(ID, var, index, iStart, nCount) result(values)
    implicit none
    real(c_double), pointer, dimension (:) :: values
    integer, intent(in) :: id
    character*(*), intent(in) :: var
    integer, intent(in) :: index
    integer(c_int64_t), VALUE, intent(in) :: iStart
    integer(c_int64_t), VALUE, intent(in) :: nCount
    type(c_ptr) :: cptr    
    cptr = pio_get_range_d_c(ID, C0(var), index, iStart, nCount)
    call c_f_pointer(cptr, values, [nCount])
  end function pio_get_range_d

  function pio_get_i64(ID, var, index) result(values)
    use iso_c_binding
    implicit none
    integer(c_int64_t), pointer, dimension (:) :: values
    integer, intent(in) :: id
    character*(*), intent(in) :: var
    integer, intent(in) :: index
    integer*8 :: n
    type(c_ptr) :: cptr
    n = pio_length(ID, C0(var))
    cptr = pio_get_i64_c(ID, C0(var), index)
    call c_f_pointer(cptr, values, [n])
  end function pio_get_i64

  function pio_get_range_i64(ID, var, index, iStart, nCount) result(values)
    implicit none
    integer(c_int64_t), pointer, dimension (:) :: values
    integer, intent(in) :: id
    character*(*), intent(in) :: var
    integer, intent(in) :: index
    integer(c_int64_t), VALUE, intent(in) :: iStart
    integer(c_int64_t), VALUE, intent(in) :: nCount
    type(c_ptr) :: cptr
    cptr = pio_get_range_i64_c(ID, C0(var), index, iStart, nCount)
    call c_f_pointer(cptr, values, [nCount])
  end function pio_get_range_i64

end module pio_interface
