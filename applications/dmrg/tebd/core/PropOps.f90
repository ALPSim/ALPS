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
MODULE PropOps
!
! Purpose: Module to propagate TEBD form wavefunctions
! for OpenSourceTEBD v3.0
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!       8/18/10  M. L. Wall	alpha release
!
USE GlobalData
USE LinearOps
USE HamiOps
USE StateOps
USE LocalOps
USE ObOps
IMPLICIT NONE

INTERFACE LeftSweep
	MODULE PROCEDURE LeftSweep_I, LeftSweep_u
END INTERFACE LeftSweep

INTERFACE RightSweep
	MODULE PROCEDURE RightSweep_I, RightSweep_u
END INTERFACE RightSweep

INTERFACE TrotterStep
	MODULE PROCEDURE TrotterStep_t, TrotterStep_f
END INTERFACE TrotterStep

INTERFACE TrotterStepNC
MODULE PROCEDURE TrotterStepNC_tt, TrotterStepNC_ft, TrotterStepNC_tm, TrotterStepNC_fm
END INTERFACE TrotterStepNC

INTERFACE LeftSweepNC
	MODULE PROCEDURE LeftSweepNC_It, LeftSweepNC_ut,LeftSweepNC_Im, LeftSweepNC_um
END INTERFACE LeftSweepNC

INTERFACE RightSweepNC
	MODULE PROCEDURE RightSweepNC_It, RightSweepNC_ut,RightSweepNC_Im, RightSweepNC_um
END INTERFACE RightSweepNC

INTERFACE CanonicalIzeNC
	MODULE PROCEDURE CanonicalizeNC_t, CanonicalizeNC_m
END INTERFACE CanonicalizeNC

INTERFACE ImagTimePropNC
MODULE PROCEDURE ImagTimePropNC_t, ImagTimePropNC_m
END INTERFACE ImagTimePropNC
	
CONTAINS

SUBROUTINE LeftSweep_I(Gammas,Lambdas,localTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(OUT) :: localTruncerr
INTEGER :: i
REAL(KIND=rKind) :: truncSum
TYPE(matrix) :: eyeMat

ALLOCATE(eyeMat%m(localSize*localSize,localSize*localSize))
eyeMat%m=0.0_rKind
DO i=1,SIZE(eyeMat%m,1)
	eyeMat%m(i,i)=1.0_rKind
END DO
truncSum=0.0_rKind

i=systemSize-1
CALL TwoSiteOpEdge(i,eyeMat%m,Gammas,Lambdas,localTruncerr)
truncSum=truncSum+localTruncerr
DO i=systemSize-2,1,(-1)
	CALL TwoSiteOpR(i,eyeMat%m,Gammas,Lambdas,localTruncerr)
    truncSum=truncSum+localTruncerr
END DO

DEALLOCATE(eyeMat%m)
localTruncerr=truncSum

END SUBROUTINE LeftSweep_I

SUBROUTINE LeftSweep_u(U,Gammas,Lambdas,localTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(OUT) :: localTruncerr
REAL(KIND=rKind) :: truncSum
INTEGER :: i
TYPE(matrix), POINTER :: U(:)

truncSum=0.0_rKind

i=systemSize-1
CALL TwoSiteOpEdge(i,U(i)%m,Gammas,Lambdas,localTruncerr)
truncSum=truncSum+localTruncerr
DO i=systemSize-2,1,(-1)
	CALL TwoSiteOpR(i,U(i)%m,Gammas,Lambdas,localTruncerr)
    truncSum=truncSum+localTruncerr
END DO
localTruncerr=truncSum


END SUBROUTINE LeftSweep_u

SUBROUTINE RightSweep_I(Gammas,Lambdas,localTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(OUT) :: localTruncerr
REAL(KIND=rKind) :: truncSum
INTEGER :: i
TYPE(matrix) :: eyeMat

ALLOCATE(eyeMat%m(localSize*localSize,localSize*localSize))
eyeMat%m=0.0_rKind
DO i=1,SIZE(eyeMat%m,1)
	eyeMat%m(i,i)=1.0_rKind
END DO

truncSum=0.0_rKind

i=1
CALL TwoSiteOpEdge(i,eyeMat%m,Gammas,Lambdas,localTruncerr)
truncSum=truncSum+localTruncerr

DO i=2,systemSize-1
	CALL TwoSiteOpL(i,eyeMat%m,Gammas,Lambdas,localTruncerr)
    truncSum=truncSum+localTruncerr
END DO

DEALLOCATE(eyeMat%m)
localTruncerr=truncSum

END SUBROUTINE RightSweep_I

SUBROUTINE RightSweep_u(U,Gammas,Lambdas,localTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(OUT) :: localTruncerr
REAL(KIND=rKind) :: truncSum
INTEGER :: i
TYPE(matrix), POINTER :: U(:)

truncSum=0.0_rKind

i=1
CALL TwoSiteOpEdge(i,U(i)%m,Gammas,Lambdas,localTruncerr)
truncSum=truncSum+localTruncerr

DO i=2,systemSize-1
	CALL TwoSiteOpL(i,U(i)%m,Gammas,Lambdas,localTruncerr)
    truncSum=truncSum+localTruncerr
END DO

localTruncerr=truncSum

END SUBROUTINE RightSweep_u

SUBROUTINE TrotterStep_t(Udt, Gammas, Lambdas, totalTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(secOrdProp) :: Udt
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(INOUT) ::totalTruncerr 
REAL(KIND=rKind) :: localTruncerr

CALL Rightsweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr

END SUBROUTINE TrotterStep_t

SUBROUTINE TrotterStep_f(Udt, Gammas, Lambdas, totalTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(fourOrdProp) :: Udt
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(INOUT) ::totalTruncerr 
REAL(KIND=rKind) :: localTruncerr

CALL Rightsweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL Rightsweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Udt%U2,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL Rightsweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL Rightsweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL Rightsweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL Rightsweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL Rightsweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL Rightsweep(Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL Rightsweep(Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL Rightsweep(Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL Rightsweep(Udt%U2,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL Rightsweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweep(Udt%U1,Gammas,Lambdas,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr

END SUBROUTINE TrotterStep_f

SUBROUTINE Canonicalize(Gammas,Lambdas)
!
!Purpose: Make all Bipartite splittings canonical. 
!Used to reorthogonalize the Schmidt basis after an ITP timestep.
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind) :: trunc

CALL Rightsweep(Gammas,Lambdas, trunc)
CALL Leftsweep(Gammas,Lambdas, trunc)

END SUBROUTINE Canonicalize

SUBROUTINE ImagTimeProp(H, GammasInner, LambdasInner)
!
!Purpose: Imaginary time propagaton algorithm-find the ground state
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(DecomposedMPO), POINTER :: MPO(:)
TYPE(tensor), POINTER :: GammasOuter(:), GammasInner(:)
TYPE(vector), POINTER :: LambdasOuter(:), LambdasInner(:)
TYPE(secOrdProp) :: U2
TYPE(fourOrdProp) :: U4
COMPLEX(KIND=rKind) :: eye, dt
COMPLEX(KIND=rKIND) :: numberInner
REAL(KIND=rKind) :: LrgstLamsBfr(systemSize+1), LrgstLamsAft(systemSize+1), GapBfrAft(systemSize+1) 
REAL(KIND=rKind) :: cenLamBfr, cenLamAft, testTol, convCr, totalTruncerr, energy
INTEGER, DIMENSION(1) :: lrgstPoint
INTEGER :: i, j, l,  chiInc, l1, numIter

!!! Define the trotter time steps for ITP. 
eye=CMPLX(0.0,1.0,KIND=rKind)
IF(trotterOrder==2) THEN
	CALL AllocateProp(U2)
ELSE
	CALL AllocateProp(U4)
END IF

numIter=SIZE(chivals)
CALL SplitOperatorList(H, MPO, 10.0_rKind**(-10))

DO i=1,numIter
	chiLimit=chivals(i)
	dt=-eye*dtITPvals(i)
	convCr=convCriterion(i)
	IF(i.ne.1) THEN
		stepsForjudge=FLOOR(stepsForJudge*ABS(dt/dtITPvals(i-1)))
	END IF
	IF(trotterOrder==2) THEN
		CALL ConstructPropagators(H, U2, dt)
	ELSE
		CALL ConstructPropagators(H, U4, dt)
	END IF

	IF(print_switch) THEN
		PRINT *, 'iteration',i,'chi',chiLimit,'dt',ABS(dt)
	END IF
	!Store the largest lambda of each splitting initially
	DO l=1,systemSize+1
		LrgstLamsBfr(l)=LambdasInner(l)%v(1)
	END DO
	!!! We iterate the time step until the iteration time reaches 'maxITPsteps' 
	!!! or until the convergence criterion is satisfied.
	DO j=1,maxITPsteps
		!!! This 'if' statement is for judging the convergence.
		!!! We judge the convergence when the number of iterations j is a multiple of 'stepsForJudge'.
		IF(MOD(j,stepsForJudge)==0) THEN
			CALL Canonicalize(GammasInner,LambdasInner)
			CALL TotalOneSite(numberInner, DiagOp%m,GammasInner, LambdasInner)
			CALL TotalEnergy(energy,MPO, GammasInner, LambdasInner)
			!Store the largest lambda of each splitting after stepsForJudge time steps
			DO l=1,systemSize+1
				LrgstLamsAft(l)=LambdasInner(l)%v(1)
				!Find the percent differences of each lambda
				GapBfrAft(l)=ABS((LrgstLamsBfr(l)-LrgstLamsAft(l))/LrgstLamsBfr(l))
			END DO
			!The lambda with the largest percent difference after stepsForJudge time steps determines convergence
			testTol=MAXVAL(GapBfrAft)
			!Store the location of the lambda with the largest percent difference after stepsForJudge time steps 
			lrgstPoint=MAXLOC(GapBfrAft)
			IF(print_switch) THEN
				PRINT *, 'ITP step j', j, 'lambda with largest percent difference', LambdasInner(lrgstPoint(1))%v(1), &
						'found at position', lrgstPoint(1)
				PRINT *, 'Percent difference', testTol,'convergence Criterion', convCr
				PRINT *, 'Diagonal expectation', REAL(numberInner), 'Energy',energy
			END IF
			IF(testTol < convCr) EXIT
			!Store the new largest lambda of each splitting
				DO l=1,systemSize+1
					LrgstLamsBfr(l)=LrgstLamsAft(l)
				END DO
		END IF

		!Time step
		IF(trotterOrder==2) THEN
			CALL TrotterStep(U2, GammasInner, LambdasInner, totalTruncerr)
		ELSE
			CALL TrotterStep(U4, GammasInner, LambdasInner, totalTruncerr)
		END IF
	END DO
END DO
IF(trotterOrder==2) THEN
	CALL DeallocateProp(U2) 
ELSE
	CALL DeallocateProp(U4) 
END IF

CALL DestroySplitOperator(MPO)


END SUBROUTINE ImagTimeProp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Number conserving method starts !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE LeftSweepNC_It(Gammas,Lambdas, LabelLeft, LabelRight,localTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: localTruncerr
REAL(KIND=rKind) :: truncSum
INTEGER :: i
TYPE(matrix) :: eyeMat

ALLOCATE(eyeMat%m(localSize*localSize,localSize*localSize))
eyeMat%m=0.0_rKind
DO i=1,SIZE(eyeMat%m,1)
	eyeMat%m(i,i)=1.0_rKind
END DO

truncSum=0.0_rKind
i=systemSize-1
	CALL TwoSiteOpNCEdge(i,eyeMat%m,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr

DO i=systemSize-2,1,(-1)
	CALL TwoSiteOpNCR(i,eyeMat%m,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr
END DO
localTruncerr=truncSum
DEALLOCATE(eyeMat%m)

END SUBROUTINE LeftSweepNC_It

SUBROUTINE LeftSweepNC_Im(Gammas,Lambdas, LabelLeft, LabelRight,localTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: localTruncerr
REAL(KIND=rKind) :: truncSum
INTEGER :: i
TYPE(matrix) :: eyeMat

ALLOCATE(eyeMat%m(localSize*localSize,localSize*localSize))
eyeMat%m=0.0_rKind
DO i=1,SIZE(eyeMat%m,1)
	eyeMat%m(i,i)=1.0_rKind
END DO

truncSum=0.0_rKind
i=systemSize-1
	CALL TwoSiteOpNCEdge(i,eyeMat%m,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr

DO i=systemSize-2,1,(-1)
	CALL TwoSiteOpNCR(i,eyeMat%m,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr
END DO
localTruncerr=truncSum
DEALLOCATE(eyeMat%m)

END SUBROUTINE LeftSweepNC_Im

SUBROUTINE LeftSweepNC_ut(U,Gammas,Lambdas, LabelLeft, LabelRight, localTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: localTruncerr
REAL(KIND=rKind) :: truncSum
INTEGER :: i
TYPE(matrix), POINTER :: U(:)

truncSum=0.0_rKind
i=systemSize-1
	CALL TwoSiteOpNCEdge(i,U(i)%m,Gammas,Lambdas, LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr

DO i=systemSize-2,1,(-1)
	CALL TwoSiteOpNCR(i,U(i)%m,Gammas,Lambdas, LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr
END DO
localTruncerr=truncSum
END SUBROUTINE LeftSweepNC_ut

SUBROUTINE LeftSweepNC_um(U,Gammas,Lambdas, LabelLeft, LabelRight, localTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: localTruncerr
REAL(KIND=rKind) :: truncSum
INTEGER :: i
TYPE(matrix), POINTER :: U(:)

truncSum=0.0_rKind
i=systemSize-1
	CALL TwoSiteOpNCEdge(i,U(i)%m,Gammas,Lambdas, LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr

DO i=systemSize-2,1,(-1)
	CALL TwoSiteOpNCR(i,U(i)%m,Gammas,Lambdas, LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr
END DO
localTruncerr=truncSum
END SUBROUTINE LeftSweepNC_um

SUBROUTINE RightSweepNC_It(Gammas,Lambdas, LabelLeft, LabelRight,localTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: localTruncerr
REAL(KIND=rKind) :: truncSum
INTEGER :: i
TYPE(matrix) :: eyeMat

ALLOCATE(eyeMat%m(localSize*localSize,localSize*localSize))
eyeMat%m=0.0_rKind
DO i=1,SIZE(eyeMat%m,1)
	eyeMat%m(i,i)=1.0_rKind
END DO

truncSum=0.0_rKind

i=1
	CALL TwoSiteOpNCEdge(i,eyeMat%m,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr

DO i=2,systemSize-1
	CALL TwoSiteOpNCL(i,eyeMat%m,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr
END DO

localTruncerr=truncSum
DEALLOCATE(eyeMat%m)

END SUBROUTINE RightSweepNC_It

SUBROUTINE RightSweepNC_Im(Gammas,Lambdas, LabelLeft, LabelRight,localTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: localTruncerr
REAL(KIND=rKind) :: truncSum
INTEGER :: i
TYPE(matrix) :: eyeMat

ALLOCATE(eyeMat%m(localSize*localSize,localSize*localSize))
eyeMat%m=0.0_rKind
DO i=1,SIZE(eyeMat%m,1)
	eyeMat%m(i,i)=1.0_rKind
END DO

truncSum=0.0_rKind
i=1
	CALL TwoSiteOpNCEdge(i,eyeMat%m,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr
DO i=2,systemSize-1
	CALL TwoSiteOpNCL(i,eyeMat%m,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr
END DO
localTruncerr=truncSum
DEALLOCATE(eyeMat%m)

END SUBROUTINE RightSweepNC_Im

SUBROUTINE RightSweepNC_ut(U,Gammas,Lambdas, LabelLeft, LabelRight, localTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: localTruncerr
REAL(KIND=rKind) :: truncSum
INTEGER :: i
TYPE(matrix), POINTER :: U(:)

truncSum=0.0_rKind

i=1
	CALL TwoSiteOpNCEdge(i,U(i)%m,Gammas,Lambdas, LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr

DO i=2,systemSize-1
	CALL TwoSiteOpNCL(i,U(i)%m,Gammas,Lambdas, LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr
END DO
localTruncerr=truncSum
END SUBROUTINE RightSweepNC_ut

SUBROUTINE RightSweepNC_um(U,Gammas,Lambdas, LabelLeft, LabelRight, localTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: localTruncerr
REAL(KIND=rKind) :: truncSum
INTEGER :: i
TYPE(matrix), POINTER :: U(:)

truncSum=0.0_rKind

i=1
	CALL TwoSiteOpNCEdge(i,U(i)%m,Gammas,Lambdas, LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr

DO i=2,systemSize-1
	CALL TwoSiteOpNCL(i,U(i)%m,Gammas,Lambdas, LabelLeft, LabelRight,localTruncerr)
    truncSum=truncSum+localTruncerr
END DO
localTruncerr=truncSum
END SUBROUTINE RightSweepNC_um


SUBROUTINE TrotterStepNC_tt(Udt, Gammas, Lambdas,LabelLeft, LabelRight, totalTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(secOrdProp) :: Udt
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(INOUT) ::totalTruncerr 
REAL(KIND=rKind) :: localTruncerr

CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr

END SUBROUTINE TrotterStepNC_tt

SUBROUTINE TrotterStepNC_tm(Udt, Gammas, Lambdas,LabelLeft, LabelRight, totalTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(secOrdProp) :: Udt
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(INOUT) ::totalTruncerr 
REAL(KIND=rKind) :: localTruncerr

CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr

END SUBROUTINE TrotterStepNC_tm

SUBROUTINE TrotterStepNC_ft(Udt, Gammas, Lambdas,LabelLeft, LabelRight, totalTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(fourOrdProp) :: Udt
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(INOUT) ::totalTruncerr 
REAL(KIND=rKind) :: localTruncerr

CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U2,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U2,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr

END SUBROUTINE TrotterStepNC_ft

SUBROUTINE TrotterStepNC_fm(Udt, Gammas, Lambdas,LabelLeft, LabelRight, totalTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction 
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(fourOrdProp) :: Udt
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(INOUT) ::totalTruncerr 
REAL(KIND=rKind) :: localTruncerr



CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U2,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U2,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL RightsweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr
CALL LeftSweepNC(Udt%U1,Gammas,Lambdas,LabelLeft, LabelRight,localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr

END SUBROUTINE TrotterStepNC_fm


SUBROUTINE CanonicalizeNC_t(Gammas,Lambdas,LabelLeft, LabelRight)
!
!Purpose: Make all Bipartite splittings canonical. 
!Used to reorthogonalize the Schmidt basis after an ITP timestep.
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind) :: trunc

CALL RightsweepNC(Gammas,Lambdas,LabelLeft, LabelRight, trunc)
CALL LeftsweepNC(Gammas,Lambdas,LabelLeft, LabelRight, trunc)

END SUBROUTINE CanonicalizeNC_t

SUBROUTINE CanonicalizeNC_m(Gammas,Lambdas,LabelLeft, LabelRight)
!
!Purpose: Make all Bipartite splittings canonical. 
!Used to reorthogonalize the Schmidt basis after an ITP timestep.
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind) :: trunc

CALL RightsweepNC(Gammas,Lambdas,LabelLeft, LabelRight, trunc)
CALL LeftsweepNC(Gammas,Lambdas,LabelLeft, LabelRight, trunc)

END SUBROUTINE CanonicalizeNC_m


SUBROUTINE ImagTimePropNC_t(H, GammasInner, LambdasInner,LabelLeftInner, LabelRightInner)
!
!Purpose: Imaginary time propagaton algorithm-find the ground state
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(DecomposedMPO), POINTER :: MPO(:)
TYPE(tensor), POINTER :: GammasOuter(:), GammasInner(:)
TYPE(vector), POINTER :: LambdasOuter(:), LambdasInner(:)
TYPE(vectorInt), POINTER :: LabelLeftOuter(:), LabelRightOuter(:)
TYPE(vectorInt), POINTER :: LabelLeftInner(:), LabelRightInner(:)
TYPE(secOrdProp) :: U2
TYPE(fourOrdProp) :: U4
COMPLEX(KIND=rKIND) :: numberInner
COMPLEX(KIND=rKind) :: eye, dt
REAL(KIND=rKind) :: LrgstLamsBfr(systemSize+1), LrgstLamsAft(systemSize+1), GapBfrAft(systemSize+1) 
REAL(KIND=rKind) ::  cenLamBfr, cenLamAft, testTol, convCr, totalTruncerr, energy
INTEGER, DIMENSION(1) :: lrgstPoint
INTEGER :: i, j, l,  chiInc, l1, numIter

!!! Define the trotter time steps for ITP. 
eye=CMPLX(0.0,1.0,KIND=rKind)
IF(trotterOrder==2) THEN
	CALL AllocateProp(U2)
ELSE
	CALL AllocateProp(U4)
END IF

numIter=SIZE(chivals)
CALL SplitOperatorList(H, MPO, 10.0_rKind**(-10))

DO i=1,numIter
	chiLimit=chivals(i)
	dt=-eye*dtITPvals(i)
	convCr=convCriterion(i)
	IF(i.ne.1) THEN
		stepsForjudge=FLOOR(stepsForJudge*ABS(dt/dtITPvals(i-1)))
	END IF
	IF(trotterOrder==2) THEN
		CALL ConstructPropagators(H, U2, dt)
	ELSE
		CALL ConstructPropagators(H, U4, dt)
	END IF

	IF(print_switch) THEN
		PRINT *, 'iteration',i,'chi',chiLimit,'dt',ABS(dt)
	END IF
	!Store the largest lambda of each splitting initially
	DO l=1,systemSize+1
		LrgstLamsBfr(l)=LambdasInner(l)%v(1)
	END DO
	!!! We iterate the time step until the iteration time reaches 'maxITPsteps' 
	!!! or until the convergence criterion is satisfied.
	DO j=1,maxITPsteps
		!!! This 'if' statement is for judging the convergence.
		!!! We judge the convergence when the number of iterations j is a multiple of 'stepsForJudge'.
		IF(MOD(j,stepsForJudge)==0) THEN
			CALL CanonicalizeNC(GammasInner,LambdasInner,LabelLeftInner,LabelRightInner)
			CALL TotalOneSite(numberInner, DiagOp%m,GammasInner, LambdasInner)
			CALL TotalEnergy(energy,MPO, GammasInner, LambdasInner)
			!Store the largest lambda of each splitting after stepsForJudge time steps
			DO l=1,systemSize+1
				LrgstLamsAft(l)=LambdasInner(l)%v(1)
				!Find the percent differences of each lambda
				GapBfrAft(l)=ABS((LrgstLamsBfr(l)-LrgstLamsAft(l))/LrgstLamsBfr(l))
			END DO
			!The lambda with the largest percent difference after stepsForJudge time steps determines convergence
			testTol=MAXVAL(GapBfrAft)
			!Store the location of the lambda with the largest percent difference after stepsForJudge time steps 
			lrgstPoint=MAXLOC(GapBfrAft)
			IF(print_switch) THEN
				PRINT *, 'ITP step j', j, 'lambda with largest percent difference', LambdasInner(lrgstPoint(1))%v(1), &
						'found at position', lrgstPoint(1)
				PRINT *, 'Percent difference', testTol,'convergence Criterion', convCr
				PRINT *, 'Diagonal expectation', REAL(numberInner), 'Energy',energy
			END IF
			IF(testTol < convCr) EXIT
			!Store the new largest lambda of each splitting
				DO l=1,systemSize+1
					LrgstLamsBfr(l)=LrgstLamsAft(l)
				END DO
		END IF

		!Time step
		IF(trotterOrder==2) THEN
			CALL TrotterStepNC(U2, GammasInner, LambdasInner,LabelLeftInner,LabelRightInner, totalTruncerr)
		ELSE
			CALL TrotterStepNC(U4, GammasInner, LambdasInner,LabelLeftInner,LabelRightInner, totalTruncerr)
		END IF
	END DO

END DO
IF(trotterOrder==2) THEN
	CALL DeallocateProp(U2) 
ELSE
	CALL DeallocateProp(U4) 
END IF

CALL DestroySplitOperator(MPO)


END SUBROUTINE ImagTimePropNC_t

SUBROUTINE ImagTimePropNC_m(H, GammasInner, LambdasInner,LabelLeftInner, LabelRightInner)
!
!Purpose: Imaginary time propagaton algorithm-find the ground state
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(DecomposedMPO), POINTER :: MPO(:)
TYPE(matrix), POINTER :: GammasOuter(:), GammasInner(:)
TYPE(vector), POINTER :: LambdasOuter(:), LambdasInner(:)
TYPE(vectorInt), POINTER :: LabelLeftOuter(:), LabelRightOuter(:)
TYPE(vectorInt), POINTER :: LabelLeftInner(:), LabelRightInner(:)
TYPE(secOrdProp) :: U2
TYPE(fourOrdProp) :: U4
COMPLEX(KIND=rKind) :: numberInner
COMPLEX(KIND=rKind) :: eye, dt
REAL(KIND=rKind) :: LrgstLamsBfr(systemSize+1), LrgstLamsAft(systemSize+1), GapBfrAft(systemSize+1) 
REAL(KIND=rKind) ::  cenLamBfr, cenLamAft, testTol, convCr, totalTruncerr, energy
INTEGER, DIMENSION(1) :: lrgstPoint
INTEGER :: i, j, l,  chiInc, l1, numIter

!!! Define the trotter time steps for ITP. 
eye=CMPLX(0.0,1.0,KIND=rKind)
IF(trotterOrder==2) THEN
	CALL AllocateProp(U2)
ELSE
	CALL AllocateProp(U4)
END IF

numIter=SIZE(chivals)
CALL SplitOperatorList(H, MPO, 10.0_rKind**(-10))

DO i=1,numIter
	chiLimit=chivals(i)
	dt=-eye*dtITPvals(i)
	convCr=convCriterion(i)
	IF(i.ne.1) THEN
		stepsForjudge=FLOOR(stepsForJudge*ABS(dt/dtITPvals(i-1)))
	END IF
	IF(trotterOrder==2) THEN
		CALL ConstructPropagators(H, U2, dt)
	ELSE
		CALL ConstructPropagators(H, U4, dt)
	END IF

	IF(print_switch) THEN
		PRINT *, 'iteration',i,'chi',chiLimit,'dt',ABS(dt)
	END IF
	!Store the largest lambda of each splitting initially
	DO l=1,systemSize+1
		LrgstLamsBfr(l)=LambdasInner(l)%v(1)
	END DO
	!!! We iterate the time step until the iteration time reaches 'maxITPsteps' 
	!!! or until the convergence criterion is satisfied.
	DO j=1,maxITPsteps
		!!! This 'if' statement is for judging the convergence.
		!!! We judge the convergence when the number of iterations j is a multiple of 'stepsForJudge'.
		IF(MOD(j,stepsForJudge)==0) THEN
			CALL CanonicalizeNC(GammasInner,LambdasInner,LabelLeftInner,LabelRightInner)
			CALL TotalOneSite(numberInner, DiagOp%m,GammasInner, LambdasInner,LabelLeftInner)
			CALL TotalEnergy(energy,MPO, GammasInner, LambdasInner,LabelLeftInner)
			!Store the largest lambda of each splitting after stepsForJudge time steps
			DO l=1,systemSize+1
				LrgstLamsAft(l)=LambdasInner(l)%v(1)
				!Find the percent differences of each lambda
				GapBfrAft(l)=ABS((LrgstLamsBfr(l)-LrgstLamsAft(l))/LrgstLamsBfr(l))
			END DO
			!The lambda with the largest percent difference after stepsForJudge time steps determines convergence
			testTol=MAXVAL(GapBfrAft)
			!Store the location of the lambda with the largest percent difference after stepsForJudge time steps 
			lrgstPoint=MAXLOC(GapBfrAft)
			IF(print_switch) THEN
				PRINT *, 'ITP step j', j, 'lambda with largest percent difference', LambdasInner(lrgstPoint(1))%v(1), &
						'found at position', lrgstPoint(1)
				PRINT *, 'Percent difference', testTol,'convergence Criterion', convCr
				PRINT *, 'Diagonal expectation', REAL(numberInner), 'Energy',energy
			END IF
			IF(testTol < convCr) EXIT
			!Store the new largest lambda of each splitting
				DO l=1,systemSize+1
					LrgstLamsBfr(l)=LrgstLamsAft(l)
				END DO
		END IF

		!Time step
		IF(trotterOrder==2) THEN
			CALL TrotterStepNC(U2, GammasInner, LambdasInner,LabelLeftInner,LabelRightInner, totalTruncerr)
		ELSE
			CALL TrotterStepNC(U4, GammasInner, LambdasInner,LabelLeftInner,LabelRightInner, totalTruncerr)
		END IF
	END DO

END DO
IF(trotterOrder==2) THEN
	CALL DeallocateProp(U2) 
ELSE
	CALL DeallocateProp(U4) 
END IF

CALL DestroySplitOperator(MPO)


END SUBROUTINE ImagTimePropNC_m


END MODULE PropOps
