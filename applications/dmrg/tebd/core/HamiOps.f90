!    Copyright (C) 2010  M. L. Wall and L. D. Carr, Colorado School of Mines
!    This file is part of OpenSourceTEBD.
!
!    This software is part of the ALPS libraries, published under the ALPS
!    Library License; you can use, redistribute it and/or modify it under
!    the terms of the license, either version 1 or (at your option) any later
!    version.
!
!    You should have received a copy of the ALPS Library License along with
!    the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
!    available from http://alps.comp-phys.org/.
!
!    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!    FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
!    SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
!    FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
!    ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!    DEALINGS IN THE SOFTWARE.
MODULE HamiOps
!
! Purpose: Module which contains operators to define various Hamiltonians, 
!procedures to generate Hamiltonians,  procedures to generate fock spaces
!and initial states with and without internal degrees of freedom, and vector 
!coupling coefficients for OpenSourceTEBD v3.0
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!       8/18/10  M. L. Wall	alpha release
!
USE GlobalData	
USE LinearOps

IMPLICIT NONE

TYPE secOrdProp
	TYPE(matrix), POINTER :: U1(:)
END TYPE secOrdProp

TYPE fourOrdProp
	TYPE(matrix), POINTER :: U1(:)
	TYPE(matrix), POINTER :: U2(:)
END TYPE fourOrdProp
	
TYPE Hamiltonian
	TYPE(matrix), POINTER :: H(:)
	TYPE(DecomposedMPO), POINTER :: HMPO(:)
END TYPE Hamiltonian 

TYPE spinParams
	!Spin model is generally H=(j_{xy}/2)*(s_+s_-+h.c.)+j_z s_zs_z-h s_z-gam s_x +d s_z^2+k \vec{S}_i\cdot \vec{S}_j
	!Global parameters
	REAL(KIND=rKIND) :: Jz !s_zs_z heisenberg coupling
	REAL(KIND=rKIND) :: Jxy ! 1/2(s_+s_-+h.c.) heisenberg coupling
	REAL(KIND=rKIND) :: h !magnetic field along z
	REAL(KIND=rKIND) :: gam !magnetic field along x
	REAL(KIND=rKIND) :: d !single-ion anisotropy
	REAL(KIND=rKIND) :: k !biquadratic
	!position dependent parameters-not currently supported
!	REAL(KIND=rKIND) :: J0
!	REAL(KIND=rKIND) :: Jz0
!	REAL(KIND=rKIND) :: Jxy0
!	REAL(KIND=rKIND) :: J1
!	REAL(KIND=rKIND) :: Jp
!	REAL(KIND=rKIND) :: Jzp
!	REAL(KIND=rKIND) :: Jxyp
!	REAL(KIND=rKIND) :: Jz1
!	REAL(KIND=rKIND) :: Jxy1
END TYPE spinParams

TYPE bosonHubbardParams
	REAL(KIND=rKIND) :: mu
	REAL(KIND=rKIND) :: t
	REAL(KIND=rKIND) :: V
	REAL(KIND=rKIND) :: U
	!position dependent parameters-not currently supported
!	REAL(KIND=rKIND) :: tp
!	REAL(KIND=rKIND) :: Vp
!	REAL(KIND=rKIND) :: t0
!	REAL(KIND=rKIND) :: t1
!	REAL(KIND=rKIND) :: V0
!	REAL(KIND=rKIND) :: V1
END TYPE bosonHubbardParams

TYPE hardcorebosonParams
	REAL(KIND=rKIND) :: mu
	REAL(KIND=rKIND) :: t
	REAL(KIND=rKIND) :: V
	!position dependent parameters-not currently supported
!	REAL(KIND=rKIND) :: tp
!	REAL(KIND=rKIND) :: Vp
!	REAL(KIND=rKIND) :: t0
!	REAL(KIND=rKIND) :: t1
!	REAL(KIND=rKIND) :: V0
!	REAL(KIND=rKIND) :: V1
END TYPE hardcorebosonParams

TYPE fermionHubbardParams
	REAL(KIND=rKIND) :: mu
	REAL(KIND=rKIND) :: t
	REAL(KIND=rKIND) :: V
	REAL(KIND=rKIND) :: U
	!position dependent parameters-not currently supported
!	REAL(KIND=rKIND) :: tp
!	REAL(KIND=rKIND) :: Vp
!	REAL(KIND=rKIND) :: t0
!	REAL(KIND=rKIND) :: t1
!	REAL(KIND=rKIND) :: V0
!	REAL(KIND=rKIND) :: V1
END TYPE fermionHubbardParams

TYPE spinlessfermionsParams
	REAL(KIND=rKIND) :: mu
	REAL(KIND=rKIND) :: t
	REAL(KIND=rKIND) :: V
	!position dependent parameters-not currently supported
!	REAL(KIND=rKIND) :: tp
!	REAL(KIND=rKIND) :: Vp
!	REAL(KIND=rKIND) :: t0
!	REAL(KIND=rKIND) :: t1
!	REAL(KIND=rKIND) :: V0
!	REAL(KIND=rKIND) :: V1
END TYPE spinlessfermionsParams

TYPE HamiParams
	TYPE(spinParams) :: sp
	TYPE(bosonhubbardparams) :: bp
	TYPE(hardcorebosonparams) :: hcbp
	TYPE(fermionhubbardparams) :: fhp
	TYPE(spinlessfermionsParams) :: sfp
END TYPE HamiParams

TYPE HamiParamslist
	TYPE(spinParams), POINTER :: sp(:)
	TYPE(bosonhubbardparams), POINTER :: bp(:)
	TYPE(hardcorebosonparams), POINTER :: hcbp(:)
	TYPE(fermionhubbardparams), POINTER :: fhp(:)
	TYPE(spinlessfermionsParams), POINTER :: sfp(:)
END TYPE HamiParamslist

TYPE rtpData
	!indexed by numquenches only
	REAL(KIND=rKind), ALLOCATABLE :: tau(:) !timescale of :th quench
	INTEGER, ALLOCATABLE :: nsteps(:) !number of timesteps for the :th quench
	INTEGER, ALLOCATABLE :: stepsForStore(:) !number of timesteps between each measurement
	INTEGER, ALLOCATABLE :: nparams(:) !number of Hami params being quenched

	!indexed by numquenches then nparams
	Type(vector), POINTER :: pow(:) !power of :th quench
	Type(Charvec), POINTER :: gquench(:) !define which parameter the :th quench is changing
	Type(vector), POINTER :: gi(:) !initial value of :th quench
	Type(vector), POINTER :: gf(:) !final value of :th quench
!	REAL(KIND=rKind), ALLOCATABLE :: tau(:) !timescale of :th quench
!	REAL(KIND=rKind), ALLOCATABLE :: pow(:) !power of :th quench
!	CHARACTER(len=10), ALLOCATABLE :: gquench(:) !define which parameter the :th quench is changing
!	REAL(KIND=rKind), ALLOCATABLE :: gi(:) !initial value of :th quench
!	REAL(KIND=rKind), ALLOCATABLE :: gf(:) !final value of :th quench
!	INTEGER, ALLOCATABLE :: nsteps(:) !number of timesteps for the :th quench
!	INTEGER, ALLOCATABLE :: stepsForStore(:) !number of timesteps between each measurement
END TYPE rtpData

!Spin system parameters
REAL(KIND=rKind) :: spin !Spin of particles
INTEGER :: spinSize !Shorthand for 2*spin+1

! *** General operators ***
TYPE(matrix) :: one_op !Unit operator
TYPE(matrix) :: DiagOp !Diagonal operator (e.g. Sz, total number)

! *** Spinless operators ***
TYPE(matrix) :: a_op !Destruction operator
TYPE(matrix) :: t_op !Tunneling operator
TYPE(matrix) :: n_op !Tunneling operator

! *** Fermionic phase operator ***
TYPE(matrix) :: fermiPhase_op

! *** Internal degrees of freedom variables ***
TYPE(vectorInt) :: Conserv !Conserv holds an Abelian conserved quantity

! *** Internal degrees of freedom operators ***
!The list index is the internal degrees of freedom index
TYPE(matrix), POINTER :: a_opS(:) !Destruction operator.  Index is idof component
TYPE(matrix) :: Sx 
TYPE(matrix) :: Splus 
TYPE(matrix) :: Sminus 
TYPE(matrix) :: Sz 
TYPE(matrix) :: exchange_xy 
TYPE(matrix) :: biquadratic



! Interfaces for generic procedures
INTERFACE AllocateProp
	MODULE PROCEDURE AllocateProp_s, AllocateProp_f
END INTERFACE AllocateProp

INTERFACE DeAllocateProp
	MODULE PROCEDURE DeAllocateProp_s, DeAllocateProp_f
END INTERFACE DeAllocateProp

INTERFACE ConstructPropagators
	MODULE PROCEDURE ConstructPropagators_s, ConstructPropagators_f
END INTERFACE ConstructPropagators

INTERFACE ASSIGNMENT (=)
	MODULE PROCEDURE spinParams_assign, bosonHubbardParams_assign, hardcorebosonParams_assign,&
			 fermionHubbardParams_assign, spinlessfermionsParams_assign
END INTERFACE

INTERFACE DestroySplitOperator
	MODULE PROCEDURE DestroySplitOperator_s, DestroySplitOperator_v
END INTERFACE DestroySplitOperator


CONTAINS


SUBROUTINE spinParams_assign(spinP, spinPP)
!
!Purpose: overload the assignment opeartor to assign a set of spin parameters to another set
!
IMPLICIT NONE
TYPE(spinParams), INTENT(OUT) :: spinP
TYPE(spinParams), INTENT(IN) :: spinPP

spinP%Jz=spinPP%Jz
spinP%Jxy=spinPP%Jxy
spinP%h=spinPP%h
spinP%gam=spinPP%gam
spinP%k=spinPP%k
spinP%d=spinPP%d

END SUBROUTINE spinParams_assign

SUBROUTINE bosonhubbardParams_assign(bosonP, bosonPP)
!
!Purpose: overload the assignment opeartor to assign a set of boson hubbard parameters to another set
!
IMPLICIT NONE
TYPE(bosonHubbardParams), INTENT(OUT) :: bosonP
TYPE(bosonHubbardParams), INTENT(IN) :: bosonPP

bosonP%t=bosonPP%t
bosonP%u=bosonPP%U
bosonP%V=bosonPP%V
bosonP%mu=bosonPP%mu

END SUBROUTINE bosonhubbardParams_assign

SUBROUTINE hardcorebosonParams_assign(bosonP, bosonPP)
!
!Purpose: overload the assignment opeartor to assign a set of hardcore boson parameters to another set
!
IMPLICIT NONE
TYPE(hardcorebosonParams), INTENT(OUT) :: bosonP
TYPE(hardcorebosonParams), INTENT(IN) :: bosonPP

bosonP%t=bosonPP%t
bosonP%V=bosonPP%V
bosonP%mu=bosonPP%mu

END SUBROUTINE hardcorebosonParams_assign

SUBROUTINE fermionhubbardParams_assign(fermionP, fermionPP)
!
!Purpose: overload the assignment opeartor to assign a set of fermion hubbard parameters to another set
!
IMPLICIT NONE
TYPE(fermionHubbardParams),INTENT(OUT) :: fermionP
TYPE(fermionHubbardParams), INTENT(IN) :: fermionPP

fermionP%t=fermionPP%t
fermionP%u=fermionPP%U
fermionP%V=fermionPP%V
fermionP%mu=fermionPP%mu

END SUBROUTINE fermionhubbardParams_assign

SUBROUTINE spinLessfermionsParams_assign(fermionP, fermionPP)
!
!Purpose: overload the assignment opeartor to assign a set of spinless fermions parameters to another set
!
IMPLICIT NONE
TYPE(spinLessfermionsParams), INTENT(OUT) :: fermionP
TYPE(spinLessfermionsParams),  INTENT(IN) :: fermionPP

fermionP%t=fermionPP%t
fermionP%V=fermionPP%V
fermionP%mu=fermionPP%mu

END SUBROUTINE spinLessfermionsParams_assign

SUBROUTINE DefineModel(HamiType)
!
!Purpose: Allocate model operators
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) :: HamiType

SELECT CASE(HamiType)
	CASE ('spin')
		CALL SetupSpin()
	CASE ('boson Hubbard')
		CALL SetupBosonHubbard()
	CASE ('hardcore boson')
		CALL SetupHardCoreBoson()
	CASE ('fermion Hubbard')
		CALL SetupFermionHubbard()
	CASE ('spinless fermions')
		CALL SetupSpinlessFermions()
	CASE DEFAULT
		PRINT *, "Hamiltonian type not recognized!"
		PRINT *, "Use 'spin', 'boson Hubbard', 'hardcore boson', 'fermion Hubbard', or 'spinless fermions'."
		STOP 
END SELECT

END SUBROUTINE DefineModel

SUBROUTINE DestroyModel(HamiType)
!
!Purpose: Deallocate the model operators
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) :: HamiType

SELECT CASE(HamiType)
	CASE ('spin')
		CALL DestroySpin()
	CASE ('boson Hubbard')
		CALL DestroyBosonHubbard()
	CASE ('hardcore boson')
		CALL DestroyHardCoreBoson()
	CASE ('fermion Hubbard')
		CALL DestroyFermionHubbard()
	CASE ('spinless fermions')
		CALL DestroySpinlessFermions()
	CASE DEFAULT
		PRINT *, "Hamiltonian type not recognized!"
		PRINT *, "Use 'spin', 'boson Hubbard', 'hardcore boson', 'fermion Hubbard', or 'spinless fermions'."
		STOP 
END SELECT

END SUBROUTINE DestroyModel


SUBROUTINE SetupSpin()
!
!Purpose: Define the Spin model
!
IMPLICIT NONE
INTEGER :: i
TYPE(spinParams) :: params
NAMELIST /spinp/ spin!, params%J, params%Jz, params%Jxy, params%h, params%gam, params%d, params%k

!get spin parameters from input
OPEN(138,file=nmlName)
READ(138,spinp)
CLOSE(138)

spinSize=FLOOR(2.0*spin+1.0)
localSize=spinSize
Nmax=spinSize-1

!Build spin operators
ALLOCATE(one_op%m(localSize, localSize))
ALLOCATE(Diagop%m(localSize, localSize))
ALLOCATE(Splus%m(localSize, localSize))
ALLOCATE(Sminus%m(localSize, localSize))
ALLOCATE(Sx%m(localSize, localSize))
ALLOCATE(Sz%m(localSize, localSize))
ALLOCATE(exchange_xy%m(localSize*localSize, localSize*localSize))
ALLOCATE(biquadratic%m(localSize*localSize, localSize*localSize))

one_op%m=0.0_rKind
DO i=1,localSize,1
	one_op%m(i,i)=1.0_rKind
END DO

Splus%m=0.0_rKind
Sminus%m=0.0_rKind
Sz%m=0.0_rKind

DO i=1,spinSize
	Sz%m(i,i)=spin-(i-1)*1.0_rKind
	IF(i.ne.spinSize) THEN
		Splus%m(i,i+1)=SQRT((spin+(spin-(i-1)*1.0_rKind))*(spin+1.0_rKind-(spin-(i-1)*1.0_rKind)))
	END IF
END DO
Sminus%m=TRANSPOSE(Splus%m)
Sx%m=0.5_rKind*(Splus%m+Sminus%m)

exchange_xy%m=0.5_rKind*(TensorProd(sPlus%m,sMinus%m)+TRANSPOSE(TensorProd(sPlus%m,sMinus%m)))
biquadratic%m=MATMUL(exchange_xy%m+TensorProd(Sz%m,Sz%m),exchange_xy%m+TensorProd(Sz%m,Sz%m))
Diagop%m=Sz%m

END SUBROUTINE SetupSpin

SUBROUTINE DestroySpin()
!
!Purpose: Destroy the Spin model
!
IMPLICIT NONE
INTEGER :: i
TYPE(spinParams) :: params

!Destroy spin operators
DEALLOCATE(one_op%m)
DEALLOCATE(Diagop%m)
DEALLOCATE(Splus%m)
DEALLOCATE(Sminus%m)
DEALLOCATE(Sx%m)
DEALLOCATE(Sz%m)
DEALLOCATE(exchange_xy%m)
DEALLOCATE(biquadratic%m)

END SUBROUTINE DestroySpin

SUBROUTINE SetupBosonHubbard()
!
!Purpose: Define the Boson hubbard model
!
IMPLICIT NONE
INTEGER :: i
TYPE(bosonhubbardParams) :: params
NAMELIST /bosonp/ Nmax!, params%J, params%Jz, params%Jxy, params%h, params%gam, params%d, params%k

!get boson parameters from input
OPEN(138,file=nmlName)
READ(138,bosonp)
CLOSE(138)

localSize=Nmax+1

!Build boson operators
ALLOCATE(one_op%m(localSize, localSize))
ALLOCATE(Diagop%m(localSize, localSize))
ALLOCATE(a_op%m(localSize, localSize))
ALLOCATE(n_op%m(localSize, localSize))
ALLOCATE(t_op%m(localSize*localSize, localSize*localSize))

one_op%m=0.0_rKind
a_op%m=0.0_rKind
DO i=1,localSize,1
	one_op%m(i,i)=1.0_rKind
	IF(i.ne.localSize) THEN
		a_op%m(i,i+1)=REAL(SQRT(i*1.0_rKind),KIND=rKind)
	END IF
END DO
!Define the tunneling operator	
t_op%m = tensorProd(TRANSPOSE(a_op%m),a_op%m) + TRANSPOSE(tensorProd(TRANSPOSE(a_op%m),a_op%m))
!Define the number operator
n_op%m=MATMUL(TRANSPOSE(a_op%m),a_op%m)
Diagop%m=n_op%m

END SUBROUTINE SetupBosonHubbard

SUBROUTINE DestroyBosonHubbard()
!
!Purpose: Destroy the Boson hubbard model
!
IMPLICIT NONE
INTEGER :: i

!Destroy boson operators
DEALLOCATE(one_op%m)
DEALLOCATE(Diagop%m)
DEALLOCATE(a_op%m)
DEALLOCATE(n_op%m)
DEALLOCATE(t_op%m)

END SUBROUTINE DestroyBosonHubbard

SUBROUTINE SetuphardcoreBoson()
!
!Purpose: Define the hardcore boson model
!
IMPLICIT NONE
INTEGER :: i
TYPE(hardcorebosonParams) :: params

Nmax=1
localSize=Nmax+1

!Build boson operators
ALLOCATE(one_op%m(localSize, localSize))
ALLOCATE(Diagop%m(localSize, localSize))
ALLOCATE(a_op%m(localSize, localSize))
ALLOCATE(n_op%m(localSize, localSize))
ALLOCATE(t_op%m(localSize*localSize, localSize*localSize))

one_op%m=0.0_rKind
a_op%m=0.0_rKind
DO i=1,localSize,1
	one_op%m(i,i)=1.0_rKind
	IF(i.ne.localSize) THEN
		a_op%m(i,i+1)=REAL(SQRT(i*1.0_rKind),KIND=rKind)
	END IF
END DO
!Define the tunneling operator	
t_op%m = tensorProd(TRANSPOSE(a_op%m),a_op%m) + TRANSPOSE(tensorProd(TRANSPOSE(a_op%m),a_op%m))
!Define the number operator
n_op%m=MATMUL(TRANSPOSE(a_op%m),a_op%m)
Diagop%m=n_op%m

END SUBROUTINE SetuphardcoreBoson

SUBROUTINE DestroyhardcoreBoson()
!
!Purpose: Destroy the hardcore boson model
!
IMPLICIT NONE
INTEGER :: i

DEALLOCATE(one_op%m)
DEALLOCATE(Diagop%m)
DEALLOCATE(a_op%m)
DEALLOCATE(n_op%m)
DEALLOCATE(t_op%m)

END SUBROUTINE DestroyhardcoreBoson

SUBROUTINE SetupFermionHubbard()
!
!Purpose: Define the fermion hubbard model
!
IMPLICIT NONE
INTEGER :: i, j, fermiswitch, k,l
TYPE(fermionhubbardParams) :: params
TYPE(matrix) :: stateList,newState
REAL(KIND=rKind) :: fPhase

Nmax=2
spin=0.5_rKind
spinSize=2
localSize=4
idof=.true.

ALLOCATE(a_opS(spinSize))
DO i=1,spinSize
	ALLOCATE(a_opS(i)%m(localSize, localSize))
END DO
ALLOCATE(Splus%m(localSize, localSize))
ALLOCATE(Sminus%m(localSize, localSize))
ALLOCATE(Sz%m(localSize, localSize))
ALLOCATE(one_op%m(localSize, localSize))
ALLOCATE(Diagop%m(localSize, localSize))
ALLOCATE(n_op%m(localSize, localSize))
ALLOCATE(fermiPhase_op%m(localSize, localSize))
ALLOCATE(t_op%m(localSize*localSize,localSize*localSize))

ALLOCATE(stateList%m(localSize, spinSize))
ALLOCATE(newState%m(1, spinSize))
ALLOCATE(Conserv%vi(localSize))
fermiSwitch=1
CALL onsiteStateListIdof(stateList%m, spinSize, fermiSwitch)

one_op%m=0.0_rKind
DO i=1,localSize,1
	one_op%m(i,i)=1.0_rKind
END DO

!Define the destruction operator with the fermi phase
DO k=1,spinSize,1
	a_opS(k)%m = 0.0_rKind
	DO j = 1, localSize
		newState%m(1, :) = stateList%m(j, :)
		newState%m(1, k) = newState%m(1, k) - 1.0_rKind
		DO i = 1, localSize
			fPhase=0
			DO l=1,k-1
				fPhase=fPhase+FLOOR(REAL(stateList%m(j, l)))
			END DO		
			a_opS(k)%m(i, j) = ((-1.0_rKind)**fPhase)*kronDelta(stateList%m(i, :), &
							newState%m(1, :), spinSize)
		END DO			
	END DO
END DO

!Define spin operators
Splus%m=0.5_8*MATMUL(TRANSPOSE(a_opS(1)%m),a_opS(2)%m)
Sminus%m=TRANSPOSE(Splus%m)
Sz%m=0.5_rKind*(MATMUL(TRANSPOSE(a_opS(1)%m),a_opS(1)%m)-MATMUL(TRANSPOSE(a_opS(2)%m),a_opS(2)%m))

n_op%m=0.0_rKind
DO k=1,spinSize
	n_op%m=n_op%m+MATMUL(TRANSPOSE(a_opS(k)%m),a_opS(k)%m)
END DO
Diagop%m=n_op%m

!Construct the fermi phase operator to define phase for n-site observables
fermiPhase_op%m=0.0_rKind
DO i=1,localSize
	fermiPhase_op%m(i,i)=(-1.0_rKIND)**(REAL(n_op%m(i,i)))
END DO

!Store the on-site quantum numbers
Conserv%vi=0
DO i=1,localSize,1
	Conserv%vi(i)=FLOOR(REAL(n_op%m(i,i)))
END DO

!Define the tunneling operator with fermi phase
t_op%m=0.0
DO k=1,spinSize,1
	t_op%m=t_op%m+MATMUL(TRANSPOSE(TensorProd(a_opS(k)%m,one_op%m)),TensorProd(fermiPhase_op%m,a_opS(k)%m))&
						+MATMUL(TRANSPOSE(TensorProd(fermiPhase_op%m,a_opS(k)%m)),TensorProd(a_opS(k)%m,one_op%m))
END DO

DEALLOCATE(stateList%m)
DEALLOCATE(newState%m)

END SUBROUTINE SetupFermionHubbard

SUBROUTINE DestroyFermionHubbard()
!
!Purpose: Destroy the fermion hubbard model
!
IMPLICIT NONE
INTEGER :: i, j, fermiswitch, k

DO i=1,spinSize
	DEALLOCATE(a_opS(i)%m)
END DO
DEALLOCATE(a_opS)

DEALLOCATE(Splus%m)
DEALLOCATE(Sminus%m)
DEALLOCATE(Sz%m)
DEALLOCATE(one_op%m)
DEALLOCATE(Diagop%m)
DEALLOCATE(n_op%m)
DEALLOCATE(fermiPhase_op%m)
DEALLOCATE(t_op%m)
DEALLOCATE(Conserv%vi)


END SUBROUTINE DestroyFermionHubbard


SUBROUTINE SetupspinlessFermions()
!
!Purpose: Define the spinless fermions model
!
IMPLICIT NONE
INTEGER :: i
TYPE(spinlessfermionsParams) :: params

Nmax=1
localSize=Nmax+1

!Build fermi operators
ALLOCATE(one_op%m(localSize, localSize))
ALLOCATE(Diagop%m(localSize, localSize))
ALLOCATE(a_op%m(localSize, localSize))
ALLOCATE(n_op%m(localSize, localSize))
ALLOCATE(fermiPhase_op%m(localSize, localSize))
ALLOCATE(t_op%m(localSize*localSize, localSize*localSize))

one_op%m=0.0_rKind
a_op%m=0.0_rKind
DO i=1,localSize,1
	one_op%m(i,i)=1.0_rKind
	IF(i.ne.localSize) THEN
		a_op%m(i,i+1)=REAL(SQRT(i*1.0_rKind),KIND=rKind)
	END IF
END DO

!Define the number operator
n_op%m=MATMUL(TRANSPOSE(a_op%m),a_op%m)
Diagop%m=n_op%m

!Construct the fermi phase operator to define phase for n-site observables
fermiPhase_op%m=0.0_rKind
DO i=1,localSize
	fermiPhase_op%m(i,i)=(-1.0_rKIND)**(REAL(n_op%m(i,i)))
END DO

!Define the tunneling operator	
t_op%m=MATMUL(TRANSPOSE(TensorProd(a_op%m,one_op%m)),TensorProd(fermiPhase_op%m,a_op%m))&
					+MATMUL(TRANSPOSE(TensorProd(fermiPhase_op%m,a_op%m)),TensorProd(a_op%m,one_op%m))

END SUBROUTINE SetupspinlessFermions

SUBROUTINE DestroyspinlessFermions()
!
!Purpose: Destroy the spinless fermions model
!
IMPLICIT NONE
INTEGER :: i

!destroy fermi operators
DEALLOCATE(one_op%m)
DEALLOCATE(Diagop%m)
DEALLOCATE(a_op%m)
DEALLOCATE(n_op%m)
DEALLOCATE(fermiPhase_op%m)
DEALLOCATE(t_op%m)

END SUBROUTINE DestroyspinlessFermions

SUBROUTINE SetHamiltonian(H,HamiType,Hparams)
!
!Purpose: Setup a Hamiltonian based on parameters in a namelist file-used for itp and initial rtp
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(HamiParams) :: Hparams
CHARACTER(len=*), INTENT(IN) :: HamiType

SELECT CASE(HamiType)
	CASE ('spin')
		CALL ReadSpinHami(Hparams%sp)
		CALL SetSpinHami(H,Hparams%sp)
	CASE ('boson Hubbard')
		CALL ReadBosonHubbardHami(Hparams%bp)
		CALL SetBosonHubbardHami(H,Hparams%bp)
	CASE ('hardcore boson')
		CALL ReadHardcoreBosonHami(Hparams%hcbp)
		CALL SetHardCoreBosonHami(H,Hparams%hcbp)
	CASE ('fermion Hubbard')
		CALL ReadFermionHubbardHami(Hparams%fhp)
		CALL SetFermionHubbardHami(H,Hparams%fhp)
	CASE ('spinless fermions')
		CALL ReadspinlessFermionsHami(Hparams%sfp)
		CALL SetSpinlessFermionsHami(H,Hparams%sfp)
	CASE DEFAULT
		PRINT *, "Hamiltonian type not recognized!"
		PRINT *, "Use 'spin', 'boson Hubbard', 'hardcore boson', 'fermion Hubbard', or 'spinless fermions'."
		STOP 
END SELECT

END SUBROUTINE SetHamiltonian

SUBROUTINE ReadSpinHami(spinP)
!
!Purpose: read spin parameters in from NameList
!
IMPLICIT NONE
TYPE(spinParams) :: spinP
NAMELIST /sp/ spinP

OPEN(138,file=nmlName)
READ(138,sp)
CLOSE(138)

END SUBROUTINE ReadSpinHami

SUBROUTINE ReadBosonHubbardHami(bosonp)
!
!Purpose: read boson hubbard parameters in from NameList
!
IMPLICIT NONE
TYPE(bosonHubbardParams) :: bosonp
NAMELIST /bp/ bosonp

OPEN(138,file=nmlName)
READ(138,bp)
CLOSE(138)

END SUBROUTINE ReadBosonHubbardHami

SUBROUTINE ReadHardcoreBosonHami(hcbosonp)
!
!Purpose: read hardcoreboson parameters in from NameList
!
IMPLICIT NONE
TYPE(hardcorebosonParams) :: hcbosonp
NAMELIST /hcbp/ hcbosonp

OPEN(138,file=nmlName)
READ(138,hcbp)
CLOSE(138)

END SUBROUTINE ReadHardcoreBosonHami

SUBROUTINE ReadfermionHubbardHami(fermiP)
!
!Purpose: read boson hubbard parameters in from NameList
!
IMPLICIT NONE
TYPE(fermionHubbardParams) :: fermiP
NAMELIST /fp/ fermiP

OPEN(138,file=nmlName)
READ(138,fp)
CLOSE(138)

END SUBROUTINE ReadfermionHubbardHami

SUBROUTINE ReadspinlessfermionsHami(sfermiP)
!
!Purpose: read boson hubbard parameters in from NameList
!
IMPLICIT NONE
TYPE(spinlessfermionsParams) :: sfermiP
NAMELIST /sfp/ sfermiP

OPEN(138,file=nmlName)
READ(138,sfp)
CLOSE(138)

END SUBROUTINE ReadspinlessfermionsHami

SUBROUTINE SetSpinHami(H,spinP)
!
!Purpose: Setup the spin hamiltonian
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(spinParams) :: spinP
INTEGER :: i
TYPE(matrix) :: szSq

ALLOCATE(szSq%m(localSize,localSize))
szSq%m=MATMUL(Sz%m, Sz%m)

ALLOCATE(H(systemSize-1))
DO i=1,systemSize-1
	ALLOCATE(H(i)%m(localSize*localSize,localSize*localSize))
	H(i)%m=spinP%Jz*TensorProd(Sz%m, Sz%m)+spinP%Jxy*exchange_xy%m+spinp%K*biquadratic%m
	H(i)%m=H(i)%m-spinP%h*HamiOneSite(Sz)-spinP%gam*HamiOneSite(Sx)+spinP%D*HamiOneSIte(Szsq)
	IF(i==1) THEN
		H(i)%m=H(i)%m-spinP%h*HamiLeft(Sz)-spinP%gam*HamiLeft(Sx)+spinP%D*HamiLeft(Szsq)
	ELSE IF(i==systemSize-1) THEN
		H(i)%m=H(i)%m-spinP%h*HamiRight(Sz)-spinP%gam*HamiRight(Sx)+spinP%D*HamiRight(Szsq)
	END IF
END DO

DEALLOCATE(szSq%m)

END SUBROUTINE SetSpinHami

SUBROUTINE SetBosonHubbardHami(H,bosonP)
!
!Purpose: Setup the boson hubbard hamiltonian
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(bosonHubbardParams) :: bosonP
INTEGER :: i
TYPE(matrix) :: nnmo

ALLOCATE(nnmo%m(locaLSize,localSize))
nnmo%m=0.5_rKind*(MATMUL(n_op%m,n_op%m)-MATMUL(n_op%m,one_op%m))

ALLOCATE(H(systemSize-1))
DO i=1,systemSize-1
	ALLOCATE(H(i)%m(localSize*localSize,localSize*localSize))
	H(i)%m=bosonP%V*TensorProd(n_op%m, n_op%m)-bosonP%t*t_op%m
	H(i)%m=H(i)%m-bosonP%mu*HamiOneSite(n_op)+bosonP%U*HamiOneSite(nnmo)
	IF(i==1) THEN
		H(i)%m=H(i)%m-bosonP%mu*HamiLeft(n_op)+bosonP%U*HamiLeft(nnmo)
	ELSE IF(i==systemSize-1) THEN
		H(i)%m=H(i)%m-bosonP%mu*HamiRight(n_op)+bosonP%U*HamiRight(nnmo)
	END IF
END DO

DEALLOCATE(nnmo%m)

END SUBROUTINE SetBosonHubbardHami

SUBROUTINE SetHardcoreBosonHami(H,hcbosonP)
!
!Purpose: Setup the hardcore boson hamiltonian
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(hardcorebosonParams) :: hcbosonP
INTEGER :: i

ALLOCATE(H(systemSize-1))
DO i=1,systemSize-1
	ALLOCATE(H(i)%m(localSize*localSize,localSize*localSize))
	H(i)%m=hcbosonP%V*TensorProd(n_op%m, n_op%m)-hcbosonP%t*t_op%m
	H(i)%m=H(i)%m-hcbosonP%mu*HamiOneSite(n_op)
	IF(i==1) THEN
		H(i)%m=H(i)%m-hcbosonP%mu*HamiLeft(n_op)
	ELSE IF(i==systemSize-1) THEN
		H(i)%m=H(i)%m-hcbosonP%mu*HamiRight(n_op)
	END IF
END DO

END SUBROUTINE SetHardcoreBosonHami


SUBROUTINE SetfermionHubbardHami(H,fermiP)
!
!Purpose: Setup the fermion hubbard hamiltonian
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(fermionHubbardParams) :: fermiP
INTEGER :: i
TYPE(matrix) :: nupndown

ALLOCATE(nupndown%m(locaLSize,localSize))
nupndown%m=MATMUL(MATMUL(TRANSPOSE(a_opS(1)%m),a_opS(1)%m),MATMUL(TRANSPOSE(a_opS(2)%m),a_opS(2)%m))

ALLOCATE(H(systemSize-1))
DO i=1,systemSize-1
	ALLOCATE(H(i)%m(localSize*localSize,localSize*localSize))
	H(i)%m=fermiP%V*TensorProd(n_op%m, n_op%m)-fermiP%t*t_op%m
	H(i)%m=H(i)%m-fermiP%mu*HamiOneSite(n_op)+fermiP%U*HamiOneSite(nupndown)
	IF(i==1) THEN
		H(i)%m=H(i)%m-fermiP%mu*HamiLeft(n_op)+fermiP%U*HamiLeft(nupndown)
	ELSE IF(i==systemSize-1) THEN
		H(i)%m=H(i)%m-fermiP%mu*HamiRight(n_op)+fermiP%U*HamiRight(nupndown)
	END IF
END DO

DEALLOCATE(nupndown%m)

END SUBROUTINE SetfermionHubbardHami

SUBROUTINE SetspinlessfermionsHami(H,sfermiP)
!
!Purpose: Setup the spinless fermions hamiltonian
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(spinlessfermionsParams) :: sfermiP
INTEGER :: i

ALLOCATE(H(systemSize-1))
DO i=1,systemSize-1
	ALLOCATE(H(i)%m(localSize*localSize,localSize*localSize))
	H(i)%m=sfermiP%V*TensorProd(n_op%m, n_op%m)-sfermiP%t*t_op%m
	H(i)%m=H(i)%m-sfermiP%mu*HamiOneSite(n_op)
	IF(i==1) THEN
		H(i)%m=H(i)%m-sfermiP%mu*HamiLeft(n_op)
	ELSE IF(i==systemSize-1) THEN
		H(i)%m=H(i)%m-sfermiP%mu*HamiRight(n_op)
	END IF
END DO

END SUBROUTINE SetspinlessfermionsHami

SUBROUTINE SetHamiltonianList(H,HamiType,Hparams,HparamsList,rtpD)
!
!Purpose: Setup list of rtp Hamiltonian parameters based on initial Hami parameters Hparams and rtp data
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
CHARACTER(len=*) :: HamiType
TYPE(HamiParams) :: Hparams
TYPE(HamiParamsList) :: HparamsList
TYPE(rtpData) :: rtpD
REAL(KIND=rKIND) :: step
INTEGER :: i, j,p, counter, stride



SELECT CASE(HamiType)
	CASE ('spin')
		!Allocate total number of time steps worth of parameter sets
		ALLOCATE(HparamsList%sp(SUM(rtpD%nsteps)+1))
		HparamsList%sp(1)=Hparams%sp
		stride=0
		DO i=1, numQuenches
			!Initialize
			DO j=1,rtpD%nsteps(i)
				HparamsList%sp(stride+j+1)=HparamsList%sp(stride+j)
			END DO

			step=rtpD%tau(i)/((rtpD%nsteps(i))*1.0_rKind)
			!Loop over different params
			DO p=1,rtpD%nparams(i)
				!catch errors			
				IF(rtpD%pow(i)%v(p).lt.0.0_8) THEN
					STOP "Negative quench powers not supported!"
				END IF
				!Cycle if power close to 0
				IF(rtpD%pow(i)%v(p).le.10.0_8**(-10)) THEN
					CYCLE
				END IF

			SELECT CASE(rtpD%gquench(i)%v(p))
				CASE('Jz')
					DO j=1,rtpD%nsteps(i)
						HparamsList%sp(stride+j+1)%Jz=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('Jxy')
					DO j=1,rtpD%nsteps(i)
						HparamsList%sp(stride+j+1)%Jxy=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('H')
					DO j=1,rtpD%nsteps(i)
						HparamsList%sp(stride+j+1)%h=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('Gamma')
					DO j=1,rtpD%nsteps(i)
						HparamsList%sp(stride+j+1)%gam=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('D')
					DO j=1,rtpD%nsteps(i)
						HparamsList%sp(stride+j+1)%d=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('K')
					DO j=1,rtpD%nsteps(i)
						HparamsList%sp(stride+j+1)%k=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE DEFAULT
					PRINT *, "quench parameter not recognized!"
					PRINT *, "use Jz, Jxy, H, Gamma, D, or K for spin models!"
					STOP
			END SELECT
			END DO
			stride=stride+rtpD%nsteps(i)
		END DO	
	CASE ('boson Hubbard')
		!Allocate total number of time steps worth of parameter sets
		stride=0
		ALLOCATE(HparamsList%bp(SUM(rtpD%nsteps)+1))
		HparamsList%bp(1)=Hparams%bp
		DO i=1, numQuenches
			!Initialize
			DO j=1,rtpD%nsteps(i)
				HparamsList%bp(stride+j+1)=HparamsList%bp(stride+j)
			END DO

			step=rtpD%tau(i)/((rtpD%nsteps(i))*1.0_rKind)
			!Loop over different params
			DO p=1,rtpD%nparams(i)
				!catch errors			
				IF(rtpD%pow(i)%v(p).lt.0.0_8) THEN
					STOP "Negative quench powers not supported!"
				END IF
				!Cycle if power close to 0
				IF(rtpD%pow(i)%v(p).le.10.0_8**(-10)) THEN
					CYCLE
				END IF
			SELECT CASE(rtpD%gquench(i)%v(p))
				CASE('mu')
					DO j=1,rtpD%nsteps(i)
						HparamsList%bp(stride+j+1)%mu=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('t')
					DO j=1,rtpD%nsteps(i)
						HparamsList%bp(stride+j+1)%t=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('V')
					DO j=1,rtpD%nsteps(i)
						HparamsList%bp(stride+j+1)%V=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('U')
					DO j=1,rtpD%nsteps(i)
						HparamsList%bp(stride+j+1)%U=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE DEFAULT
					PRINT *, "quench parameter not recognized!"
					PRINT *, "use mu, t, V, or U for boson hubbard models!"
					STOP
			END SELECT
			END DO
			stride=stride+rtpD%nsteps(i)
		END DO

	CASE ('hardcore boson')
		!Allocate total number of time steps worth of parameter sets
		ALLOCATE(HparamsList%hcbp(SUM(rtpD%nsteps)+1))
		stride=0
		HparamsList%hcbp(1)=Hparams%hcbp
		DO i=1, numQuenches
			!Initialize
			DO j=1,rtpD%nsteps(i)
				HparamsList%hcbp(stride+j+1)=HparamsList%hcbp(stride+j)
			END DO

			step=rtpD%tau(i)/((rtpD%nsteps(i))*1.0_rKind)
			!Loop over different params
			DO p=1,rtpD%nparams(i)
				!catch errors			
				IF(rtpD%pow(i)%v(p).lt.0.0_8) THEN
					STOP "Negative quench powers not supported!"
				END IF
				!Cycle if power close to 0
				IF(rtpD%pow(i)%v(p).le.10.0_8**(-10)) THEN
					CYCLE
				END IF

			SELECT CASE(rtpD%gquench(i)%v(p))
				CASE('mu')
					DO j=1,rtpD%nsteps(i)
						HparamsList%hcbp(stride+j+1)%mu=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('t')
					DO j=1,rtpD%nsteps(i)
						HparamsList%hcbp(stride+j+1)%t=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('V')
					DO j=1,rtpD%nsteps(i)
						HparamsList%hcbp(stride+j+1)%V=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE DEFAULT
					PRINT *, "quench parameter not recognized!"
					PRINT *, "use mu, t, or V for hardcore boson models!"
					STOP
			END SELECT
			END DO
!			DO j=1,rtpD%nsteps(i)+1
!					print *, 'I P j quench params',i,p, HparamsList%hcbp(j)%mu, HparamsList%hcbp(j)%t, HparamsList%hcbp(j)%V
!			END DO
			stride=stride+rtpD%nsteps(i)
		END DO	


	CASE ('fermion Hubbard')
		!Allocate total number of time steps worth of parameter sets
		ALLOCATE(HparamsList%fhp(SUM(rtpD%nsteps)+1))
		HparamsList%fhp(1)=Hparams%fhp
		stride=0
		DO i=1, numQuenches
			!Initialize
			DO j=1,rtpD%nsteps(i)
				HparamsList%fhp(stride+j+1)=HparamsList%fhp(stride+j)
			END DO

			step=rtpD%tau(i)/((rtpD%nsteps(i))*1.0_rKind)
			!Loop over different params
			DO p=1,rtpD%nparams(i)
				!catch errors			
				IF(rtpD%pow(i)%v(p).lt.0.0_8) THEN
					STOP "Negative quench powers not supported!"
				END IF
				!Cycle if power close to 0
				IF(rtpD%pow(i)%v(p).le.10.0_8**(-10)) THEN
					CYCLE
				END IF

			SELECT CASE(rtpD%gquench(i)%v(p))
				CASE('mu')
					DO j=1,rtpD%nsteps(i)
						HparamsList%fhp(stride+j+1)%mu=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('t')
					DO j=1,rtpD%nsteps(i)
						HparamsList%fhp(stride+j+1)%t=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('V')
					DO j=1,rtpD%nsteps(i)
						HparamsList%fhp(stride+j+1)%V=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('U')
					DO j=1,rtpD%nsteps(i)
						HparamsList%fhp(stride+j+1)%U=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE DEFAULT
					PRINT *, "quench parameter not recognized!"
					PRINT *, "use mu, t, V, or U for fermion hubbard models!"
					STOP
			END SELECT
			END DO
			stride=stride+rtpD%nsteps(i)
		END DO

	CASE ('spinless fermions')
		!Allocate total number of time steps worth of parameter sets
		ALLOCATE(HparamsList%sfp(SUM(rtpD%nsteps)+1))
		stride=0
		HparamsList%sfp(1)=Hparams%sfp
		DO i=1, numQuenches
			!Initialize
			DO j=1,rtpD%nsteps(i)
				HparamsList%sfp(stride+j+1)=HparamsList%sfp(j)
			END DO

			step=rtpD%tau(i)/((rtpD%nsteps(i))*1.0_rKind)
			!Loop over different params
			DO p=1,rtpD%nparams(i)
				!catch errors			
				IF(rtpD%pow(i)%v(p).lt.0.0_8) THEN
					STOP "Negative quench powers not supported!"
				END IF
				!Cycle if power close to 0
				IF(rtpD%pow(i)%v(p).le.10.0_8**(-10)) THEN
					CYCLE
				END IF


			SELECT CASE(rtpD%gquench(i)%v(p))
				CASE('mu')
					DO j=1,rtpD%nsteps(i)
						HparamsList%sfp(stride+j+1)%mu=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('t')
					DO j=1,rtpD%nsteps(i)
						HparamsList%sfp(stride+j+1)%t=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE('V')
					DO j=1,rtpD%nsteps(i)
						HparamsList%sfp(stride+j+1)%V=rtpd%gi(i)%v(p)+((j*step/rtpD%tau(i))**rtpD%pow(i)%v(p))*(rtpd%gf(i)%v(p)-rtpd%gi(i)%v(p))
					END DO
				CASE DEFAULT
					PRINT *, "quench parameter not recognized!"
					PRINT *, "use mu, t, or V for spinless fermion models!"
					STOP
			END SELECT
			END DO
			stride=stride+rtpD%nsteps(i)
		END DO	
	CASE DEFAULT
		PRINT *, "Hamiltonian type not recognized!"
		PRINT *, "Use 'spin', 'boson Hubbard', 'hardcore boson', 'fermion Hubbard', or 'spinless fermions'."
		STOP 
END SELECT


END SUBROUTINE SetHamiltonianList


SUBROUTINE UpdateHamiltonian(H,HamiType,HparamsList,i)
!
!Purpose: Setup a Hamiltonian based on parameters in a namelist file-used for itp and initial rtp
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(HamiParamsList) :: HparamsList
INTEGER, INTENT(IN) :: i
CHARACTER(len=*), INTENT(IN) :: HamiType

SELECT CASE(HamiType)
	CASE ('spin')
		CALL SetSpinHami(H,HparamsList%sp(i))
	CASE ('boson Hubbard')
		CALL SetBosonHubbardHami(H,HparamsList%bp(i))
	CASE ('hardcore boson')
		CALL SetHardCoreBosonHami(H,HparamsList%hcbp(i))
	CASE ('fermion Hubbard')
		CALL SetFermionHubbardHami(H,HparamsList%fhp(i))
	CASE ('spinless fermions')
		CALL SetSpinlessFermionsHami(H,HparamsList%sfp(i))
	CASE DEFAULT
		PRINT *, "Hamiltonian type not recognized!"
		PRINT *, "Use 'spin', 'boson Hubbard', 'hardcore boson', 'fermion Hubbard', or 'spinless fermions'."
		STOP 
END SELECT

END SUBROUTINE UpdateHamiltonian

SUBROUTINE AllocateOps(Ops,numops,opsize)
!
!Purpose: Allocate a numops length list of opsizeXopsize matrices, name it Ops
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Ops(:)
INTEGER, INTENT(IN) :: numops,opsize
INTEGER :: i

ALLOCATE(Ops(numops), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate Ops'
END IF 
DO i=1,numops
	ALLOCATE(Ops(i)%m(opsize,opsize), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate Ops'
	END IF 
END DO

END SUBROUTINE AllocateOps
	
SUBROUTINE DeallocateOps(Ops,numops)
!
!Purpose: Deallocate a numops length list of opsizeXopsize matrices
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Ops(:)
INTEGER, INTENT(IN) :: numops
INTEGER :: i
DO i=1,numops
	DEALLOCATE(Ops(i)%m, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate Ops'
	END IF 
END DO
DEALLOCATE(Ops, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate Ops'
END IF 

END SUBROUTINE DeallocateOps


SUBROUTINE SplitOperator(Op12, OpD, trunc)
!
!Purpose: Express a two-site operator as a sum of one site operators, i.e.
!Op12%m([j1,j2],[i1,i2]) -> \sum_{d} OpD%l(d)%m(i1,j1)*OpD%r(d)%m(i2,j2),
!using the singular value decompositon.  
!The OPTIONAL argument trunc specifies a truncation of the bond dimension for
!discarded weight=trunc
!    
IMPLICIT NONE   
COMPLEX(KIND=rKind), INTENT(INOUT) :: Op12(:,:)  
TYPE(matrix) ::  temp, U, V
TYPE(DecomposedMPO) :: OpD
INTEGER :: i,j,k,l,d, counter, p
TYPE(vector) :: S
REAL(KIND=rKind) :: truncInner, norm, summer
REAL(KIND=rKind), OPTIONAL :: trunc

IF(PRESENT(trunc)) THEN
	truncInner=trunc
ELSE
	truncInner=10.0_rKind**(-15)
END IF

d=INT(SQRT(REAL(SIZE(Op12,1))))
ALLOCATE(temp%m(d*d,d*d))
!Op12([i j], [k l])-> Op12([i k], [j l])
DO i=1,d
	DO j=1,d
		DO k=1,d
			DO l=1,d
				temp%m((i-1)*d+k,(j-1)*d+l)=Op12((i-1)*d+j,(k-1)*d+l)
			END DO		
		END DO
	END DO
END DO

!Op12_{ij}->\sum Op1_{ik} S_k Op2_{kj}
CALL svd2(U,S,V,temp%m)
DEALLOCATE(temp%m)

norm=0.0_rKind
DO i=1,SIZE(S%v)
	norm=norm+S%v(i)*S%v(i)
END DO

norm=SQRT(norm)

summer=0.0_rKind
counter=0
DO i=1,SIZE(S%v)
	summer=summer+S%v(i)*S%v(i)
	S%v(i)=SQRT(S%v(i))
	counter=counter+1
	IF(1.0_rKind-SQRT(summer)/norm.lt.truncInner) EXIT
END DO

ALLOCATE(OpD%l(counter),OpD%r(counter))

DO k=1,counter
	ALLOCATE(OpD%l(k)%m(d,d),OpD%r(k)%m(d,d))
	DO i=1,d
		DO j=1,d	
			OpD%l(k)%m(i,j)=U%m((i-1)*d+j,k)*S%v(k)
			OpD%r(k)%m(i,j)=V%m(k,(i-1)*d+j)*S%v(k)
		END DO
	END DO
END DO

DEALLOCATE(U%m, V%m, S%v)

END SUBROUTINE SplitOperator

SUBROUTINE SplitOperatorList(Op12List, OpDList, trunc)
!
!Purpose: Express a two-site operator list as a sum of one site operators, i.e.
!Op12(L)%m([j1,j2],[i1,i2]) -> \sum_{d} OpD(L)%l(d)%m(i1,j1)*OpD(L)%r(d)%m(i2,j2),
!using the singular value decompositon.  
!    
IMPLICIT NONE     
TYPE(matrix), POINTER :: Op12List(:)
TYPE(DecomposedMPO), POINTER :: OpDList(:)
REAL(KIND=rKind), OPTIONAL :: trunc
INTEGER :: numList, i

numList=SIZE(Op12List)
ALLOCATE(OpDList(numList))
DO i=1,numList
	IF(PRESENT(trunc)) THEN
		CALL SplitOperator(Op12List(i)%m, OpDList(i), trunc)
	ELSE
		CALL SplitOperator(Op12List(i)%m, OpDList(i))
	END IF
END DO

END SUBROUTINE SplitOperatorList

SUBROUTINE DestroySplitOperator_s(MPO)
!
!Purpose:Deallocate a Decomposed Operator
!
TYPE(DecomposedMPO) :: MPO
INTEGER :: i

DO i=1,SIZE(MPO%l)
	DEALLOCATE(MPO%l(i)%m)
	DEALLOCATE(MPO%r(i)%m)
END DO

DEALLOCATE(MPO%l, MPO%r)

END SUBROUTINE DestroySplitOperator_s

SUBROUTINE DestroySplitOperator_v(MPO)
!
!Purpose:Deallocate a Decomposed Operator
!
TYPE(DecomposedMPO), POINTER :: MPO(:)
INTEGER :: i, j

DO j=1,SIZE(MPO)
	DO i=1,SIZE(MPO(j)%l)
		DEALLOCATE(MPO(j)%l(i)%m)
		DEALLOCATE(MPO(j)%r(i)%m)
	END DO
	DEALLOCATE(MPO(j)%l, MPO(j)%r)
END DO
DEALLOCATE(MPO)

END SUBROUTINE DestroySplitOperator_v

FUNCTION HamiOneSite(Op)
!
!Purpose: Return the contribution of the one-site operator Op to a Hamiltonian in TEBD form
!
IMPLICIT NONE
TYPE(matrix), INTENT(IN) :: Op
COMPLEX(KIND=rKind) :: HamiOneSite(localSize*localSize,localSize*localSize)

HamiOneSite=0.5_rKind*(TensorProd(Op%m,one_op%m)+TensorProd(one_op%m,Op%m))

END FUNCTION HamiOneSite


FUNCTION HamiLeft(Op)
!
!Purpose: Return the contribution of the one-site operator Op to a Hamiltonian in TEBD form
!
IMPLICIT NONE
TYPE(matrix), INTENT(IN) :: Op
COMPLEX(KIND=rKind) :: HamiLeft(localSize*localSize,localSize*localSize)

HamiLeft=0.5_rKind*TensorProd(Op%m,one_op%m)

END FUNCTION HamiLeft

FUNCTION HamiRight(Op)
!
!Purpose: Return the contribution of the one-site operator Op to a Hamiltonian in TEBD form
!
IMPLICIT NONE
TYPE(matrix), INTENT(IN) :: Op
COMPLEX(KIND=rKind) :: HamiRight(localSize*localSize,localSize*localSize)

HamiRight=0.5_rKind*TensorProd(one_op%m,Op%m)

END FUNCTION HamiRight


SUBROUTINE AllocateProp_s(U)
!
!Purpose: Allocate the Propagator based on boundary conditions
! and order of trotter decomposition
!
IMPLICIT NONE
TYPE(secOrdProp) :: U
INTEGER :: i

ALLOCATE(U%U1(systemSize-1))
DO i=1,systemSize-1
	ALLOCATE(U%U1(i)%m(localSize*localSize,localSize*localSize))
END DO

END SUBROUTINE AllocateProp_s

SUBROUTINE AllocateProp_f(U)
!
!Purpose: Allocate the Propagator based on boundary conditions
! and order of trotter decomposition
!
IMPLICIT NONE
TYPE(fourOrdProp) :: U
INTEGER :: i

ALLOCATE(U%U1(systemSize-1))
ALLOCATE(U%U2(systemSize-1))
DO i=1,systemSize-1
	ALLOCATE(U%U1(i)%m(localSize*localSize,localSize*localSize))
	ALLOCATE(U%U2(i)%m(localSize*localSize,localSize*localSize))
END DO

END SUBROUTINE AllocateProp_f

SUBROUTINE DeAllocateProp_s(U)
!
!Purpose: Allocate the Propagator based on boundary conditions
! and order of trotter decomposition
!
IMPLICIT NONE
TYPE(secOrdProp) :: U
INTEGER :: i

DO i=1,systemSize-1
	DEALLOCATE(U%U1(i)%m)
END DO
DEALLOCATE(U%U1)

END SUBROUTINE DeAllocateProp_s

SUBROUTINE DeAllocateProp_f(U)
!
!Purpose: Allocate the Propagator based on boundary conditions
! and order of trotter decomposition
!
IMPLICIT NONE
TYPE(fourOrdProp) :: U
INTEGER :: i

DO i=1,systemSize-1
	DEALLOCATE(U%U1(i)%m)
	DEALLOCATE(U%U2(i)%m)
END DO

DEALLOCATE(U%U1)
DEALLOCATE(U%U2)

END SUBROUTINE DeAllocateProp_f

SUBROUTINE ConstructPropagators_s(H, U, dt)
!
!Purpose: Construct the Trotter-Suzuki propagator U from the Hamiltonian H
!Using a routine adaptive to different boundary conditions and trotter schemes
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(secOrdProp) :: U
COMPLEX(KIND=rKind), INTENT(IN) :: dt
INTEGER :: i

DO i=1,systemSize-1
	CALL Matrix_Exponential(H(i)%m, U%U1(i)%m, dt, localSize*localSize)
END DO

END SUBROUTINE ConstructPropagators_s

SUBROUTINE ConstructPropagators_f(H, U, dt)
!
!Purpose: Construct the Trotter-Suzuki propagator U from the Hamiltonian H
!Using a routine adaptive to different boundary conditions and trotter schemes
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(fourOrdProp) :: U
COMPLEX(KIND=rKind), INTENT(IN) :: dt
INTEGER :: i

DO i=1,systemSize-1
	CALL Matrix_Exponential(H(i)%m, U%U1(i)%m, dt, localSize*localSize)
	CALL Matrix_Exponential(H(i)%m, U%U2(i)%m, -2.0_rKind*dt, localSize*localSize)
END DO

END SUBROUTINE ConstructPropagators_f

SUBROUTINE onsiteStateListIdof(list, idofSize, fermiSwitch)
!
!Purpose: Outer routine of a subroutine that returns the list of on-site states for  
!an idofSize-size internal Hilbert space truncated at maxfilling particles per site.
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: list(:, :)
INTEGER, INTENT(INOUT) :: idofSize
INTEGER i,j,n(idofSize),counter,k,l, dum3, dum4
INTEGER, INTENT(INOUT), OPTIONAL :: fermiSwitch

list = 0.0_rKind
n=0.0_rKind
counter=0
IF(PRESENT(fermiSwitch)) THEN
!Loop over total number of particles nmax, beginning with nMax=0	
	DO k=1,Nmax+1,1	
		n=0.0_rKind
		n(1)=k-1
		!FORTRAN calls by reference, and so all
		!objects passed into a recursive function
		!need to be variables-NOT indicies
		dum3=k-1
		dum4=1
		CALL onsiteIdofInner(list,dum3,counter,dum4,n, idofSize, fermiSwitch)
	END DO
ELSE
!Loop over total number of particles nmax, beginning with nMax=0	
	DO k=1,Nmax+1,1	
		n=0.0_rKind
		n(1)=k-1
		!FORTRAN calls by reference, and so all
		!objects passed into a recursive function
		!need to be variables-NOT indicies
		dum3=k-1
		dum4=1
		CALL onsiteIdofInner(list,dum3,counter,dum4,n, idofSize)
	END DO
END IF
END SUBROUTINE onsiteStateListIdof

RECURSIVE SUBROUTINE onsiteIdofInner(list, nmax,  counter, m,n, idofSize, fermiSwitch)
!
!Purpose: Inner routine of a subroutine that returns the list of on-site states for  
!an idofSize-size internal Hilbert space truncated at maxfilling particles per site.
!
!See manual for more detail
!
IMPLICIT NONE
!All passed variables in a recursive function HAVE to be INOUTs
INTEGER, INTENT(INOUT) :: nmax
INTEGER, INTENT(INOUT) :: m
INTEGER, INTENT(INOUT) :: counter
INTEGER, INTENT(INOUT) :: idofSize
INTEGER, INTENT(INOUT) :: n(:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: list(:, :)
INTEGER, INTENT(INOUT), OPTIONAL :: fermiSwitch
INTEGER i,j,ntemp(idofSize),k,l,dum1, dum2

!!!!!!!!!!!!!!!!Fermi Statistics!!!!!!!!!!!!!!
IF(PRESENT(fermiSwitch)) THEN
	!If the state is stuff,0,0,...,0 just count it
	IF (nmax==0) THEN
		IF(MAXVAL(n).le.1) THEN
			counter=counter+1
			list(counter, :)=n(:)
		END IF
		!If the state is stuff,1,0,...,0 then
		!move the 1 over sequentially, counting
		!each time
	ELSE IF (nmax==1) THEN
		IF(MAXVAL(n).le.1) THEN
			DO i=m,idofSize,1
				!define a temporary array
				ntemp=n
				!move a particle out of the ith state into
				!the mth state
				ntemp(m)=ntemp(m)-1
				ntemp(i)=ntemp(i)+1
				counter=counter+1
				list(counter, :)=ntemp(:)
			END DO
			!If the state is stuff,k,0,...,0 then
			!perform the algorithm on the subspace beginning
			!with k-see manual for more detail
		END IF
	ELSE
		!Loop over putting all k particles in the ith component
		DO i=m,idofSize,1
			!Zero out the subspace array	
			n(m:idofSize)=0
			n(i)=nmax
			!Loop over recursive calls
			DO l=1,nmax,1
				ntemp=n
				ntemp(i)=ntemp(i)+1-l
				ntemp(i+1)=ntemp(i+1)+l-1
				!FORTRAN calls by reference, and so all
				!things passed into a recursive function
				!need to be variables-NOT indices
				dum1=l-1
				dum2=i+1
				CALL onsiteIdofInner(list,dum1,counter,dum2,ntemp, idofSize, fermiSwitch)
			END DO
		END DO
	END IF
!!!!!!!!!!!!!!!!Bose Statistics!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE
	!If the state is stuff,0,0,...,0 just count it
	IF (nmax==0) THEN
		counter=counter+1
		list(counter, :)=n(:)
		!If the state is stuff,1,0,...,0 then
		!move the 1 over sequentially, counting
		!each time
	ELSE IF (nmax==1) THEN
		DO i=m,idofSize,1
			!define a temporary array
			ntemp=n
			!move a particle out of the ith state into
			!the mth state
			ntemp(m)=ntemp(m)-1
			ntemp(i)=ntemp(i)+1
			counter=counter+1
			list(counter, :)=ntemp(:)
		END DO
		!If the state is stuff,k,0,...,0 then
		!perform the algorithm on the subspace beginning
		!with k-see manual for more detail
	ELSE
		!Loop over putting all k particles in the ith component
		DO i=m,idofSize,1
			!Zero out the subspace array	
			n(m:idofSize)=0
			n(i)=nmax
			!Loop over recursive calls
			DO l=1,nmax,1
				ntemp=n
				ntemp(i)=ntemp(i)+1-l
				ntemp(i+1)=ntemp(i+1)+l-1
				!FORTRAN calls by reference, and so all
				!things passed into a recursive function
				!need to be variables-NOT indices
				dum1=l-1
				dum2=i+1
				CALL onsiteIdofInner(list,dum1,counter,dum2,ntemp, idofSize)
			END DO
		END DO
	END IF
END IF
END SUBROUTINE onsiteIdofInner

END MODULE HamiOps
