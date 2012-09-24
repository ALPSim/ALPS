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
MODULE StateOps
!
! Purpose: Module to define and manipulate states in the vidal representation
!
! Record of Revisions
!	Date		Programmer	Description of change
!	====		==========	=====================
!       8/18/10  M. L. Wall	alpha release
!
USE GlobalData
USE LinearOps
USE HamiOps

IMPLICIT NONE


INTERFACE AllocateGamLam
	MODULE PROCEDURE AllocateGamLam_t, AllocateGamLam_m
END INTERFACE AllocateGamLam

INTERFACE DeAllocateGamLam
	MODULE PROCEDURE DeAllocateGamLam_t, DeAllocateGamLam_m
END INTERFACE DeAllocateGamLam

INTERFACE CopyGamLam
	MODULE PROCEDURE CopyGamLam_t, CopyGamLam_m
END INTERFACE CopyGamLam

INTERFACE InitialSetNC
MODULE PROCEDURE InitialSetNC_t, InitialSetNC_m
END INTERFACE InitialSetNC

INTERFACE ProductStateMPD
MODULE PROCEDURE ProductStateMPD_t, ProductStateMPD_m
END INTERFACE ProductStateMPD

CONTAINS

SUBROUTINE AllocateGamLam_t(Gammas, Lambdas, chi)
!
!Purpose: Allocate gammas and Lambdas based on a value of chi.  Gammas are tensors
! as is appropriate for non-nonumber conserving systems or number conserving systems
! with internal degrees of freedom
!
IMPLICIT NONE  
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER, INTENT(IN) :: chi
INTEGER :: i, d1, d2, dl, dr, counter, comp

!Gammas live on sites-there are systemSize of them
ALLOCATE(Gammas(systemSize))
!Lambdas live on links, the extra 2 assist in computing two-site observables
ALLOCATE(Lambdas(systemSize+1) , STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate lambdas'
END IF 

dr=MIN(localSize,chi)
ALLOCATE(Gammas(1)%t(1,localSize,dr))
Gammas(1)%t=0.0_rKind
ALLOCATE(Lambdas(1)%v(1))
Lambdas(1)%v=0.0_rKind
ALLOCATE(Lambdas(2)%v(dr))
Lambdas(2)%v=0.0_rKind
ALLOCATE(Gammas(systemSize)%t(dr,localSize,1))
Gammas(systemSize)%t=0.0_rKind
ALLOCATE(Lambdas(systemSize+1)%v(1))
Lambdas(systemSize+1)%v=0.0_rKind

counter=-1
!Find the greatest power of d that is less than the maximal chi
DO i=0,systemSize
	IF(localSize**i.gt.chi) THEN
		EXIT
	ELSE
		counter=counter+1
	END IF
END DO
!Work with the exponents instead to avoid very large integers
DO i=2,systemSize-1
	comp=MIN(i-1,systemSize-i+1)
	IF(comp.le.counter) THEN
		dl=localSize**comp
	ELSE
		dl=chi
	END IF
	comp=MIN(i,systemSize-i)
	IF(comp.le.counter) THEN
		dr=localSize**comp
	ELSE
		dr=chi
	END IF
	ALLOCATE(Gammas(i)%t(dl,localSize,dr), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate gammas'
	END IF 
	Gammas(i)%t=0.0_8
	ALLOCATE(Lambdas(i+1)%v(dr), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate lambdas'
	END IF 
	Lambdas(i+1)%v=0.0_rKind
END DO

END SUBROUTINE AllocateGamLam_t

SUBROUTINE AllocateGamLam_m(Gammas, Lambdas, chi)
!
!Purpose: Allocate gammas and Lambdas based on a value of chi.  Gammas are matrices
! (no on-site index), as is appropriate for number conserving systems without internal degrees of freedom
!
IMPLICIT NONE  
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER, INTENT(IN) :: chi
INTEGER :: i, d1, d2, dl, dr, counter, comp

!Gammas live on sites-there are systemSize of them
ALLOCATE(Gammas(systemSize))
!Lambdas live on links, the extra 2 assist in computing two-site observables
ALLOCATE(Lambdas(systemSize+1) , STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate lambdas'
END IF 

dr=MIN(localSize,chi)
ALLOCATE(Gammas(1)%m(1,dr))
Gammas(1)%m=0.0_rKind
ALLOCATE(Lambdas(1)%v(1))
Lambdas(1)%v=0.0_rKind
ALLOCATE(Lambdas(2)%v(dr))
Lambdas(2)%v=0.0_rKind
ALLOCATE(Gammas(systemSize)%m(dr,1))
Gammas(systemSize)%m=0.0_rKind
ALLOCATE(Lambdas(systemSize+1)%v(1))
Lambdas(systemSize+1)%v=0.0_rKind

counter=-1
!Find the greatest power of d that is less than the maximal chi
DO i=0,systemSize
	IF(localSize**i.gt.chi) THEN
		EXIT
	ELSE
		counter=counter+1
	END IF
END DO
!Work with the exponents instead to avoid very large integers
DO i=2,systemSize-1
	comp=MIN(i-1,systemSize-i+1)
	IF(comp.le.counter) THEN
		dl=localSize**comp
	ELSE
		dl=chi
	END IF
	comp=MIN(i,systemSize-i)
	IF(comp.le.counter) THEN
		dr=localSize**comp
	ELSE
		dr=chi
	END IF
	ALLOCATE(Gammas(i)%m(dl,dr), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate gammas'
	END IF 
	Gammas(i)%m=0.0_8
	ALLOCATE(Lambdas(i+1)%v(dr), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate lambdas'
	END IF 
	Lambdas(i+1)%v=0.0_rKind
END DO

END SUBROUTINE AllocateGamLam_m

SUBROUTINE CopyGamLam_t(GammasCopy, LambdasCopy, GammasOrig, LambdasOrig)
!
!Purpose: Copy Gammas and Lambdas from Orig to Copy
!
IMPLICIT NONE
TYPE(tensor), POINTER :: GammasCopy(:), GammasOrig(:)
TYPE(vector), POINTER :: LambdasCopy(:), LambdasOrig(:)
INTEGER :: i, alpha, j, beta, chimin, chiminr
 
DO i=1,(systemSize+1)
	chiMin=MIN(SIZE(LambdasCopy(i)%v),SIZE(LambdasOrig(i)%v))
	LambdasCopy(i)%v(1:chiMin)=LambdasOrig(i)%v(1:chiMin)
END DO
	
DO i=1,systemSize
	chiMin=MIN(SIZE(GammasCopy(i)%t,1),SIZE(GammasOrig(i)%t,1))
	chiMinr=MIN(SIZE(GammasCopy(i)%t,3),SIZE(GammasOrig(i)%t,3))
	GammasCopy(i)%t(1:chiMin,:,1:chiMinr)=GammasOrig(i)%t(1:chiMin,:,1:chiMinr)
END DO
END SUBROUTINE CopyGamLam_t

SUBROUTINE CopyGamLam_m(GammasCopy, LambdasCopy, GammasOrig, LambdasOrig)
!
!Purpose: Copy Gammas and Lambdas from Orig to Copy-number conserving version
!
IMPLICIT NONE
TYPE(matrix), POINTER :: GammasCopy(:), GammasOrig(:)
TYPE(vector), POINTER :: LambdasCopy(:), LambdasOrig(:)
INTEGER :: i, alpha, j, beta, chimin, chiminr
 
DO i=1,(systemSize+1)
	chiMin=MIN(SIZE(LambdasCopy(i)%v),SIZE(LambdasOrig(i)%v))
	LambdasCopy(i)%v(1:chiMin)=LambdasOrig(i)%v(1:chiMin)
END DO
	
DO i=1,systemSize
	chiMin=MIN(SIZE(GammasCopy(i)%m,1),SIZE(GammasOrig(i)%m,1))
	chiMinr=MIN(SIZE(GammasCopy(i)%m,2),SIZE(GammasOrig(i)%m,2))
	GammasCopy(i)%m(1:chiMin,1:chiMinr)=GammasOrig(i)%m(1:chiMin,1:chiMinr)
END DO
END SUBROUTINE CopyGamLam_m

SUBROUTINE DeallocateGamLam_t(Gammas, Lambdas)
!
!Purpose: Deallocate gammas and Lambdas
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)	
INTEGER :: i
!Deallocate each site/link object
DO i=1,SIZE(Gammas)
	DEALLOCATE(Gammas(i)%t, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate Gammas'
	END IF 
	DEALLOCATE(Lambdas(i)%v, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate Lambdas'
	END IF 
END DO
!Deallocate the list of objects
DEALLOCATE(Lambdas(SIZE(Lambdas))%v, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate Lambdas'
END IF 
DEALLOCATE(Gammas, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate Gammas'
END IF 
DEALLOCATE(Lambdas, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate Lambdas'
END IF 

END SUBROUTINE DeallocateGamLam_t

SUBROUTINE DeallocateGamLam_m(Gammas, Lambdas)
!
!Purpose: Deallocate gammas and Lambdas
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)	
INTEGER :: i
!Deallocate each site/link object
DO i=1,SIZE(Gammas)
	DEALLOCATE(Gammas(i)%m, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate Gammas'
	END IF 
	DEALLOCATE(Lambdas(i)%v, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate Lambdas'
	END IF 
END DO
!Deallocate the list of objects
DEALLOCATE(Lambdas(SIZE(Lambdas))%v, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate Lambdas'
END IF 
DEALLOCATE(Gammas, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate Gammas'
END IF 
DEALLOCATE(Lambdas, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate Lambdas'
END IF 

END SUBROUTINE DeallocateGamLam_m

SUBROUTINE AllocateLabel(LabelLeft, LabelRight, chi)
!
!Purpose: Allocate labels needed for conserving a single Abelian symmetry
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: chi
INTEGER :: i, d1, d2, dl, dr, counter, comp

!LableLeft(l)%vi(alpha) gives the conserved quantity associated with the alpha^th 
!left Schmidt vector for a bipartite splitting at link l, likewise for LabelRight
ALLOCATE(LabelLeft(systemSize+1), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate labelLeft'
END IF 
ALLOCATE(LabelRight(systemSize+1), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate labelRight'
END IF 


counter=-1
DO i=0,systemSize
	IF(localSize**i.gt.chi) THEN
		EXIT
	ELSE
		counter=counter+1
	END IF
END DO

ALLOCATE(LabelLeft(1)%vi(1),LabelRight(1)%vi(1))
LabelLeft(1)%vi=10000
LabelRight(1)%vi=10000
DO i=1,systemSize
	comp=MIN(i,systemSize-i)
	IF(comp.le.counter) THEN
		dr=localSize**comp
	ELSE
		dr=chi
	END IF
	ALLOCATE(LabelLeft(i+1)%vi(dr), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate LL'
	END IF 
	LabelLeft(i+1)%vi=10000
	ALLOCATE(LabelRight(i+1)%vi(dr), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate LL'
	END IF 
	LabelRight(i+1)%vi=10000
END DO

END SUBROUTINE AllocateLabel	

SUBROUTINE CopyLabel(LabLCopy, LabRCopy, LabLOrig, LabROrig)
!
!Purpose: Copy Labels from Orig to Copy
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabLCopy(:), LabRCopy(:), LabLOrig(:), LabROrig(:)
INTEGER :: i, alpha, chimin

DO i=1,(systemSize+1)
	LabLCopy(i)%vi=10000
	LabRCopy(i)%vi=10000
	chimin=MIN(SIZE(LabLCopy(i)%vi),SIZE(LabLOrig(i)%vi))
	LabLCopy(i)%vi(1:chiMin)=LabLOrig(i)%vi(1:chiMin)
	LabRCopy(i)%vi(1:chiMin)=LabROrig(i)%vi(1:chiMin)
END DO

END SUBROUTINE CopyLabel

SUBROUTINE DeallocateLabel(LabelLeft, LabelRight)
!
!Purpose: Deallocate labels needed for conserving a single Abelian symmetry
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: i

DO i=1, (systemSize+1)
	DEALLOCATE(LabelLeft(i)%vi, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate LabelLeft'
	END IF 
	DEALLOCATE(LabelRight(i)%vi, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate LabelRight'
	END IF 
END DO
DEALLOCATE(LabelLeft, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate LabelLeft'
END IF 
DEALLOCATE(LabelRight, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate LabelRight'
END IF 

END SUBROUTINE DeallocateLabel


SUBROUTINE ProductStateMPD_t(Gammas, Lambdas, carray)
!
!Purpose: Construct the Vidal decomposition of a product state whose coefficients
!		  are stored in carray
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER :: i, j, m, d, chi
COMPLEX(KIND=rKind), INTENT(IN) :: carray(:, :)

m = SIZE(carray, 2)
d = SIZE(carray, 1)
!Check to be sure that M and d read in agree with the global variables
IF (m /= systemSize) THEN
	PRINT *, 'systemSize parameter conflicts with input data in ProductStateMPD' 
END IF	
IF (d /= localSize) THEN
	PRINT *, 'localSize parameter conflicts with input data in ProductStateMPD'
END IF
	
DO i = 1, systemSize
	Gammas(i)%t = CMPLX(0.0, KIND=rKind)
	Lambdas(i)%v = 0.0_rKind
	Lambdas(i)%v(1) = 1.0_rKind ! Assign the first component of each lambda the value 1, as this is a product state.
	DO j = 1, SIZE(Gammas(1)%t,2)
		Gammas(i)%t(1, j, 1) = carray(j, i) ! Assign the alpha=1 values of Gammas to be the coefficients of each on-site state.
	END DO
END DO
Lambdas(systemSize+1)%v = 0.0_rKind
Lambdas(systemSize+1)%v(1) = 1.0_rKind

END SUBROUTINE ProductStateMPD_t

SUBROUTINE ProductStateMPD_m(Gammas,LabelLeft,LabelRight, Lambdas, carray)
!
!Purpose: Construct the Vidal decomposition of a product state whose coefficients
!		  are stored in carray
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: i, j, m, d, chi
COMPLEX(KIND=rKind), INTENT(IN) :: carray(:, :)

m = SIZE(carray, 2)
d = SIZE(carray, 1)
!Check to be sure that M and d read in agree with the global variables
IF (m /= systemSize) THEN
	PRINT *, 'systemSize parameter conflicts with input data in ProductStateMPD' 
END IF	
IF (d /= localSize) THEN
	PRINT *, 'localSize parameter conflicts with input data in ProductStateMPD'
END IF

LabelLeft(1)%vi(1)=0	
DO i = 1, systemSize
	Gammas(i)%m = CMPLX(0.0, KIND=rKind)
	Lambdas(i)%v = 0.0_rKind
	Lambdas(i)%v(1) = 1.0_rKind ! Assign the first component of each lambda the value 1, as this is a product state.
	DO j = 1, localSize
		IF(carray(j,i).ne.0.0_rKind) THEN
			LabelLeft(i+1)%vi(1)=LabelLeft(i)%vi(1)+j-1
			Gammas(i)%m(1, 1) = carray(j, i) ! Assign the alpha=1 values of Gammas to be the coefficients of each on-site state.
		END IF
	END DO
END DO

DO i=1,systemSize+1
	LabelRight(i)%vi(1)=totQ-LabelLeft(i)%vi(1)
END DO

Lambdas(systemSize+1)%v = 0.0_rKind
Lambdas(systemSize+1)%v(1) = 1.0_rKind

END SUBROUTINE ProductStateMPD_m


SUBROUTINE ProductStateLabels(LabelLeft, LabelRight, carray,intDegFree)
!
!Purpose: Construct the lists of number conserving labels of a product state whose coefficients
!		  are stored in carray
!
!WARNING: The input state should be an eigenstate of total number.  This routine will not check for this!
!
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: i, j, m, d, chi
COMPLEX(KIND=rKind), INTENT(IN) :: carray(:, :)
INTEGER, INTENT(IN), OPTIONAL :: intDegFree

m = SIZE(carray, 2)
d = SIZE(carray, 1)

!Check to be sure that M and d read in agree with the global variables
IF (m /= systemSize) THEN
	PRINT *, 'systemSize parameter conflicts with input data in ProductStateMPD' 
END IF	
IF (d /= localSize) THEN
	PRINT *, 'localSize parameter conflicts with input data in ProductStateMPD'
END IF

!If internal degrees of freedom are present, find the occupation of each component
IF(PRESENT(intDegFree)) THEN
	LabelLeft(1)%vi(1)=0.0_rKind
	DO i=1,systemSize
		DO j=1,localSize
			IF(ABS(carray(j,i)).ne.0.0_rKind) THEN
				LabelLeft(i+1)%vi(1)=LabelLeft(i)%vi(1)+Conserv%vi(j)
				EXIT
			END IF
		END DO
	END DO
	ELSE
	LabelLeft(1)%vi(1)=0.0_rKind
	DO i=1,systemSize
		DO j=1,localSize
			IF(ABS(carray(j,i)).ne.0.0_rKind) THEN
				LabelLeft(i+1)%vi(1)=LabelLeft(i)%vi(1)+j-1
			END IF
		END DO
	END DO
END IF

!Construct LabelRight from LabelLeft
DO i=1,(systemSize+1),1
	LabelRight(i)%vi(1)=totQ-LabelLeft(i)%vi(1)
END DO
	
END SUBROUTINE ProductStateLabels


SUBROUTINE AllStates(Gammas, Lambdas)
!
!Purpose: Creates an initial state that is a product of local states 
! which contain all possible states in the same amount.  Used as initial ITP
! state for number non-conserving code.
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)	
INTEGER :: i, j

DO i=1,systemSize
	Gammas(i)%t=CMPLX(0.0, KIND=rKind)
	Lambdas(i)%v=0.0_rKind
	Lambdas(i)%v(1)=1.0_rKind
	!Each state is weighted equally by normalization by 1/sqrt(d)
	Gammas(i)%t(1,:,1)=CMPLX((1.0_rKind)/SQRT(SIZE(Gammas(i)%t,2)*1.0_rKind),KIND=rKind)
END DO
Lambdas(systemSize+1)%v=0.0_rKind
Lambdas(systemSize+1)%v(1)=1.0_rKind
END SUBROUTINE AllStates

SUBROUTINE InitialSetNC_t(Gammas, Lambdas, LabelLeft, LabelRight, intDegFree, domain, rank)
!
!Purpose: Creates an initial state consistent with number conservation
!
!   The algorithm places the particles in a "wedding cake" structure
!   with the greatest number of particles in the center of the cake ("tops")
!   A lesser number surrounding this peak ("center"), and no particles in the gap between
!   the cake and the edge ("hole").  This mimics the ground state in a harmonic trap.
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!   See the manual for details/examples concerning the algorithm
!
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN), OPTIONAL :: intDegFree, rank
INTEGER, INTENT(IN), OPTIONAL :: domain(2)
INTEGER :: minc, maxc
REAL(KIND=rKind) ::  sumsq1, sumsq2
INTEGER :: i, j, k, l, hole, tops, center,nj, locSysSize, counter
REAL(KIND=rKind), ALLOCATABLE :: rthing1(:),rthing2(:)

IF(PRESENT(rank)) THEN
	CALL Seed_Init(rank)
ELSE
	CALL Seed_Init()
END IF

IF(totQ==0) THEN
	DO i=1,(systemSize+1),1
		LabelLeft(i)%vi(1)=0
		LabelRight(i)%vi(1)=0
		Lambdas(i)%v(1)=1.0_rKind
	END DO
	DO i=1,systemSize
		Gammas(i)%t(1,1,1)=1.0_rKind
	END DO
ELSE
	IF(PRESENT(domain)) THEN
		locSysSize=domain(2)-domain(1)+1
	ELSE
		locSysSize=systemSize
	END IF

	IF(Nmax*locSysSize.lt.totQ) THEN
		PRINT *, "!!!!WARNING: totQ < Nmax*locSysSize!!!!"
		totQ=Nmax*locSysSize
	END IF
	! If the number of sites with the least number of particles (hole) is even
	If((MOD(locSysSize-MOD(totQ,locSysSize),2)==0).OR. &
	! or there is no hole (integer filling)
	(MOD(totQ,locSysSize)==0)) THEN
	!!If we have unit filling there is no hole
		If(MOD(totQ,locSysSize)==0) THEN
			hole=0
			!for non-unit filling the hole is
		ELSE
			hole=(locSysSize-MOD(totQ,locSysSize))/2
		END IF
		! Else if the number of holes is odd
	ELSE
		hole=(locSysSize-MOD(totQ,locSysSize)+1)/2
	END IF
	! Number of sites that have the greatest number of particles
	! number in the "top of the wedding cake"
	tops=MOD(totQ-1,locSysSize)+1

	! Number of sites that have the lesser number of particles
	! number in the "center of the wedding cake"
	!Floor ensures that when N=L we fill uniformly
	center=FLOOR(REAL(totQ,KIND=8)/REAL(locSysSize,KIND=8)-10.0**(-8))+1
	! LabelLeft in 0<=link<=hole
	DO i=1,(hole+1),1
		IF(i==1) THEN
			!There are never any particles to the left of the system
			LabelLeft(i)%vi(1)=0
		ELSE
			!Count the cumulative number to the left
			LabelLeft(i)%vi(1)=LabelLeft(i-1)%vi(1)+center-1
		END IF
	END DO
	! LabelLeft in hole+1<=link<=hole+top
	DO i=(hole+2),(hole+tops+1),1
		LabelLeft(i)%vi(1)=LabelLeft(i-1)%vi(1)+center
	END DO
	! LabelLeft in hole+top+1<=link<=locSysSize+1
	DO i=(hole+tops+2),(locSysSize+1),1
		LabelLeft(i)%vi(1)=LabelLeft(i-1)%vi(1)+center-1
	END DO
	! Construct LabelRight from LabelLeft
	DO i=1,(locSysSize+1),1
		LabelRight(i)%vi(1)=totQ-LabelLeft(i)%vi(1)
	END DO
	! Construct Lambdas
	DO i=1,(locSysSize+1),1
		Lambdas(i)%v(1)=1.0_rKind
	END DO
	!Internal degree(s) of freedom present
	!Find the number of states with number=center-1 and center particles
	!Store these in minc and maxc
	minc=0.0_rKind
	maxc=0.0_rKind
	DO j=1,localSize,1
		nj=Conserv%vi(j)
		IF(nj==center-1) THEN
			minc=minc+1
		ELSE IF(nj==center) THEN
			maxc=maxc+1
		END IF
	END DO
	ALLOCATE(rthing1(minc),rthing2(maxc))
	! Construct Gammas
	! Gammas in 1<=site<=hole
	DO i=1,hole,1
		DO l=1,minc
			CALL RANDOM_NUMBER(rthing1(l)) 
		END DO
		rthing1=rthing1/SQRT(DOT_PRODUCT(rthing1,rthing1))
		!Sum over internal degree of freedom, weighting each randomly
		counter=1
		DO j=1,localSize,1
			nj=Conserv%vi(j)
			IF(nj==center-1) THEN
				Gammas(i)%t(1,j,1)=rthing1(counter)
				counter=counter+1
			END IF
		END DO
	END DO
	! Gammas in hole+1<=site<=hole+top
	DO i=hole+1,(hole+tops),1
		DO l=1,maxc
			CALL RANDOM_NUMBER(rthing2(l)) 
		END DO
		rthing2=rthing2/SQRT(DOT_PRODUCT(rthing2,rthing2))
		!Sum over internal degree of freedom, weighting each randomly
		counter=1
		DO j=1,localSize,1
			nj=Conserv%vi(j)
			IF(nj==center) THEN
				Gammas(i)%t(1,j,1)=rthing2(counter)
				counter=counter+1
			END IF
		END DO
	END DO
	! Gammas in hole+top+1<=site<=systemSize
	DO i=(hole+tops+1),locSysSize,1
		DO l=1,minc
			CALL RANDOM_NUMBER(rthing1(l)) 
		END DO
		rthing1=rthing1/SQRT(DOT_PRODUCT(rthing1,rthing1))
		!Sum over internal degree of freedom, weighting each randomly
		counter=1
		DO j=1,localSize,1
			nj=Conserv%vi(j)
			IF(nj==center-1) THEN
				Gammas(i)%t(1,j,1)=rthing1(counter)
				counter=counter+1
			END IF
		END DO
	END DO
	DEALLOCATE(rthing1, rthing2)
	DO i=locSysSize+1,(systemSize+1),1
		LabelLeft(i)%vi(1)=totQ
		LabelRight(i)%vi(1)=0
		Lambdas(i)%v(1)=1.0_rKind
	END DO
	DO i=locSysSize+1,systemSize
		Gammas(i)%t(1,1,1)=1.0_rKind
	END DO
END IF

END SUBROUTINE InitialSetNC_t

SUBROUTINE InitialSetNC_m(Gammas, Lambdas, LabelLeft, LabelRight, domain, rank)
!
!Purpose: Creates an initial state consistent with number conservation
!
!   The algorithm places the particles in a "wedding cake" structure
!   with the greatest number of particles in the center of the cake ("tops")
!   A lesser number surrounding this peak ("center"), and no particles in the gap between
!   the cake and the edge ("hole").  This mimics the ground state in a harmonic trap.
!
!   See the manual for details/examples concerning the algorithm
!
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN), OPTIONAL ::  rank
INTEGER, INTENT(IN), OPTIONAL :: domain(2)
INTEGER :: minc, maxc
REAL(KIND=rKind) ::  sumsq1, sumsq2
INTEGER :: i, j, k, l, hole, tops, center,nj, locSysSize, counter
REAL(KIND=rKind), ALLOCATABLE :: rthing1(:),rthing2(:)

IF(PRESENT(rank)) THEN
	CALL Seed_Init(rank)
ELSE
	CALL Seed_Init()
END IF

IF(totQ==0) THEN
	DO i=1,(systemSize+1),1
		LabelLeft(i)%vi(1)=0
		LabelRight(i)%vi(1)=0
		Lambdas(i)%v(1)=1.0_rKind
	END DO
	DO i=1,systemSize
		Gammas(i)%m(1,1)=1.0_rKind
	END DO
ELSE
	IF(PRESENT(domain)) THEN
		locSysSize=domain(2)-domain(1)+1
	ELSE
		locSysSize=systemSize
	END IF

	IF(Nmax*locSysSize.lt.totQ) THEN
		PRINT *, "!!!!WARNING: totQ < Nmax*locSysSize!!!!"
		totQ=Nmax*locSysSize
	END IF
	! If the number of sites with the least number of particles (hole) is even
	If((MOD(locSysSize-MOD(totQ,locSysSize),2)==0).OR. &
	! or there is no hole (integer filling)
	(MOD(totQ,locSysSize)==0)) THEN
	!!If we have unit filling there is no hole
		If(MOD(totQ,locSysSize)==0) THEN
			hole=0
			!for non-unit filling the hole is
		ELSE
			hole=(locSysSize-MOD(totQ,locSysSize))/2
		END IF
		! Else if the number of holes is odd
	ELSE
		hole=(locSysSize-MOD(totQ,locSysSize)+1)/2
	END IF
	! Number of sites that have the greatest number of particles
	! number in the "top of the wedding cake"
	tops=MOD(totQ-1,locSysSize)+1

	! Number of sites that have the lesser number of particles
	! number in the "center of the wedding cake"
	!Floor ensures that when N=L we fill uniformly
	center=FLOOR(REAL(totQ,KIND=8)/REAL(locSysSize,KIND=8)-10.0**(-8))+1
	! LabelLeft in 0<=link<=hole
	DO i=1,(hole+1),1
		IF(i==1) THEN
			!There are never any particles to the left of the system
			LabelLeft(i)%vi(1)=0
		ELSE
			!Count the cumulative number to the left
			LabelLeft(i)%vi(1)=LabelLeft(i-1)%vi(1)+center-1
		END IF
	END DO
	! LabelLeft in hole+1<=link<=hole+top
	DO i=(hole+2),(hole+tops+1),1
		LabelLeft(i)%vi(1)=LabelLeft(i-1)%vi(1)+center
	END DO
	! LabelLeft in hole+top+1<=link<=locSysSize+1
	DO i=(hole+tops+2),(locSysSize+1),1
		LabelLeft(i)%vi(1)=LabelLeft(i-1)%vi(1)+center-1
	END DO
	! Construct LabelRight from LabelLeft
	DO i=1,(locSysSize+1),1
		LabelRight(i)%vi(1)=totQ-LabelLeft(i)%vi(1)
	END DO
	! Construct Lambdas
	DO i=1,(locSysSize+1),1
		Lambdas(i)%v(1)=1.0_rKind
	END DO
	! Construct Gammas
	! Gammas in 1<=site<=hole
	DO i=1,hole,1
		Gammas(i)%m(1,1)=CMPLX(1.0,KIND=rKind)
	END DO
	! Gammas in hole+1<=site<=hole+top
	DO i=hole+1,(hole+tops),1
		Gammas(i)%m(1,1)=CMPLX(1.0,KIND=rKind)
	END DO
	! Gammas in hole+top+1<=site<=systemSize
	DO i=(hole+tops+1),locSysSize,1
		Gammas(i)%m(1,1)=CMPLX(1.0,KIND=rKind)
	END DO

	DO i=locSysSize+1,(systemSize+1),1
		LabelLeft(i)%vi(1)=totQ
		LabelRight(i)%vi(1)=0
		Lambdas(i)%v(1)=1.0_rKind
	END DO
	DO i=locSysSize+1,systemSize
		Gammas(i)%m(1,1)=1.0_rKind
	END DO
END IF

END SUBROUTINE InitialSetNC_m


END MODULE StateOps
