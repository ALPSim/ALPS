!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ALPS Project: Algorithms and Libraries for Physics Simulations
!
! ALPS Libraries
!
! Copyright (C) 2011 by Synge Todo <wistaria@comp-phys.org>
!
! This software is part of the ALPS libraries, published under the ALPS
! Library License; you can use, redistribute it and/or modify it under
! the terms of the license, either version 1 or (at your option) any later
! version.
! 
! You should have received a copy of the ALPS Library License along with
! the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
! available from http://alps.comp-phys.org/.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
! FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
! SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
! FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
! ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
! DEALINGS IN THE SOFTWARE.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! subroutine alps_run
subroutine alps_run()
  return
end subroutine alps_run

! subroutine alps_init
subroutine alps_init(caller)
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer :: caller(2)
  real*8 :: dval
  integer :: ival
  character*32 :: cval

  write(0, *) ' ##### alps_init() ##### '
  call alps_get_parameter(dval, 'X', ALPS_DOUBLE_PRECISION, caller)
  write(0, *) '  parameter X     = ',dval

  call alps_get_parameter(ival, 'Y', ALPS_INT, caller)
  write(0, *) '  parameter Y     = ',ival

  call alps_get_parameter(cval, "WORLD", ALPS_CHAR, caller)
  write(0, *) '  parameter WORLD = ',cval

  ival = 0
  call alps_parameter_defined(ival, 'Z', caller);
  write(0, *) '  defined parameter Z = ',ival, ''

  return
end subroutine alps_init

! alps_init_observables
subroutine alps_init_observables(caller)
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer caller(2)
  ! write(*,*) 'alps_init_observables()'
  return
end subroutine alps_init_observables

! alps_progerss
subroutine alps_progress(prgrs, caller)
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer caller(2)
  real*8 prgrs
  ! write(*,*) 'alps_init_observables()'
  prgrs = 1.0
end subroutine alps_progress

! alps_is_thermalized
subroutine alps_is_thermalized(thrmlz, caller)
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer caller(2)
  integer thrmlz
  ! write(*, *) 'alps_is_thermalized()'
  thrmlz = 1
  return
end subroutine alps_is_thermalized

! alps_save
subroutine alps_save(caller)
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer caller(2)
  integer ival(4)
  real rval(4)
  character*64 cval(4)
  integer i
  ! write(*, *) 'alps_save()'

  do i = 1, 4
     ival(i) = i * 10
     rval(i) = i * i * 0.1
  end do
  cval(1) = "one"
  cval(2) = "TWO"
  cval(3) = "Three"
  cval(4) = "four"

  call alps_dump(ival, 4, ALPS_INT, caller)
  call alps_dump(rval, 4, ALPS_REAL, caller)
  call alps_dump(cval, 4, ALPS_CHAR, caller)

  return
end subroutine alps_save

! alps_load
subroutine alps_load(caller)
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer caller(2)
  !   write(*, *) 'alps_load()'
  return
end subroutine alps_load

! alps_finalize
subroutine alps_finalize(caller)
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer caller(2)
  ! write(*, *) 'alps_finalize()'
  return
end subroutine alps_finalize
