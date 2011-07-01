! オリジナルコード：宮下精二著「熱・統計力学」(培風館 1993年) p.263

module ising_mod
  implicit none
  real, parameter :: V0 = .465661288D-9

  integer, allocatable, dimension(:) :: IP, IM
  integer, allocatable, dimension(:,:) :: IS
  real*8, allocatable, dimension(:) :: P
  integer :: K, MCS, INT, L, IX
  real :: TEMP
  !$omp threadprivate (K, MCS, INT, TEMP, IP, IM, P, IS, IX, L)
end module ising_mod

! subroutine alps_init
subroutine alps_init(caller)
  use ising_mod
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer :: caller(2)
  integer :: i, j
  real*8 :: W

  call alps_get_parameter(TEMP, "TEMPERATURE", ALPS_REAL, caller)
  call alps_get_parameter(L, "L", ALPS_INT, caller)
  call alps_get_parameter(MCS, "MCS", ALPS_INT, caller)
  call alps_get_parameter(INT, "INT", ALPS_INT, caller)
  call alps_get_parameter(IX, "WORKER_SEED", ALPS_INT, caller)

  allocate( IP(L) )
  allocate( IM(L) )
  allocate( P(-4:4) )
  allocate( IS(L, L) )

  K = 0

  do i = -4, 4
     W = exp(float(i)/TEMP)
     P(i) = W / (W + 1/W)
  end do

  do i = 1, L
     IP(i) = i + 1
     IM(i) = i - 1
  end do

  do i = 1, L
     do j = 1, L
        IS(i, j) = 1
     end do
  end do

  IP(L) = 1
  IM(1) = L

  return
end subroutine alps_init

! subroutine alps_run
subroutine alps_run(caller)
  use ising_mod
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer :: caller(2)
  integer :: i, j, M
  real*8 :: EN, MG

  do i = 1, L
     do j = 1, L
        M = IS(IP(i), j) + IS(i, IP(j)) + IS(IM(i), j) + IS(i, IM(j))
        IS(i, j) = -1

        IX = IAND(IX * 5 * 11, 2147483647) ! this should be replaced with a better RNG
        if(P(M).gt.V0*IX) IS(i, j) = 1
     end do
  end do

  EN = 0.0D0
  MG = 0.0D0
  do i = 1, L
     do j = 1, L
        EN = EN + IS(i, j) * (IS(IP(i), j) + IS(i, IP(j)))
        MG = MG + IS(i, j)
     end do
  end do

 call alps_accumulate_observable(EN, 1, ALPS_DOUBLE_PRECISION, "Energy", caller)
 call alps_accumulate_observable(MG, 1, ALPS_DOUBLE_PRECISION, "Magnetization", caller)
  K = K + 1

  return
end subroutine alps_run

! alps_init_observables
subroutine alps_init_observables(caller)
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer :: caller(2)

  call alps_init_observable(1, ALPS_REAL, "Energy", caller)
  call alps_init_observable(1, ALPS_REAL, "Magnetization", caller)

  return
end subroutine alps_init_observables

! alps_progerss
subroutine alps_progress(prgrs, caller)
  use ising_mod
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer :: caller(2)
  real*8 :: prgrs

  prgrs = K / (INT + MCS)

end subroutine alps_progress

! alps_is_thermalized
subroutine alps_is_thermalized(thrmlz, caller)
  use ising_mod
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer :: caller(2)
  integer :: thrmlz

  if(K >= INT) then
     thrmlz = 1
  else
     thrmlz = 0
  end if

  return
end subroutine alps_is_thermalized

! alps_save
subroutine alps_save(caller)
  use ising_mod
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer caller(2)

  call alps_dump(K, 1, ALPS_INT, caller)
  call alps_dump(IX, 1, ALPS_INT, caller)
  call alps_dump(IS, L * L, ALPS_INT, caller)

  return
end subroutine alps_save

! alps_load
subroutine alps_load(caller)
  use ising_mod
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer :: caller(2)

  call alps_restore(K, 1, ALPS_INT, caller)
  call alps_restore(IX, 1, ALPS_INT, caller)
  call alps_restore(IS, L * L, ALPS_INT, caller)

  return
end subroutine alps_load

! alps_finalize
subroutine alps_finalize(caller)
  use ising_mod
  implicit none
  include "alps/fortran/alps_fortran.h"
  integer :: caller(2)

  deallocate(IP)
  deallocate(IM)
  deallocate(P)
  deallocate(IS)

  return
end subroutine alps_finalize

