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
MODULE ObOps
!
! Purpose: Module to compute observables
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
IMPLICIT NONE

!Type for local measures
TYPE mlocal
	COMPLEX(KIND=rKIND), POINTER :: Op(:,:)
	COMPLEX(KIND=rKIND), POINTER :: value(:)
END TYPE mlocal

!Type for site-averaged measures
TYPE mavg
	COMPLEX(KIND=rKIND), POINTER :: Op(:,:)
	COMPLEX(KIND=rKIND) :: value
END TYPE mavg

!Type for correlation functions
TYPE mcorr
	COMPLEX(KIND=rKIND), POINTER :: Op(:,:)
	TYPE(DecomposedMPO) :: MPO
	COMPLEX(KIND=rKIND), POINTER :: value(:,:)
END TYPE mcorr

!Type for correlation functions with fermi phases
TYPE mcorrf
	COMPLEX(KIND=rKIND), POINTER :: Op(:,:)
	TYPE(DecomposedMPO) :: MPO
	COMPLEX(KIND=rKIND), POINTER :: value(:,:)
END TYPE mcorrf

!Type for entropies
TYPE entropy
	REAL(KIND=rKIND) :: qme
	REAL(KIND=rKIND), POINTER :: vN(:)
	REAL(KIND=rKIND), POINTER :: chain(:)
	REAL(KIND=rKIND), POINTER :: tbvN(:,:)
END TYPE entropy

!Type for overlaps
TYPE moverlap
	COMPLEX(KIND=rKind) :: value
	LOGICAL :: m
	TYPE(tensor), POINTER :: Gammas(:)
	TYPE(vector), POINTER :: Lambdas(:)
	TYPE(matrix), POINTER :: Gammasm(:)
	TYPE(vectorInt), POINTER :: LL(:)
	TYPE(vectorInt), POINTER :: LR(:)
END TYPE moverlap

!Type for all collected measures
TYPE measure
	LOGICAL :: enPres
	TYPE(DecomposedMPO), POINTER :: HMPO(:)
	REAL(KIND=rKind) :: en
	LOGICAL :: localpres
	TYPE(mlocal), POINTER :: local(:)
	LOGICAL :: avgpres
	TYPE(mavg), POINTER :: avg(:)
	INTEGER :: entpres
	TYPE(entropy) :: ent
	LOGICAL :: corrpres
	TYPE(mcorr), POINTER :: corr(:)
	LOGICAL :: fermicorrpres
	TYPE(mcorrf), POINTER :: fermicorr(:)
	LOGICAL :: overlapPres
	TYPE(moverlap), POINTER :: overlap(:)
END TYPE measure


!*** INTERFACES


INTERFACE OneSiteExpVal
	MODULE PROCEDURE OneSiteExpVal_c, OneSiteExpVal_mc
END INTERFACE  OneSiteExpVal

INTERFACE TwoSiteExpValG
MODULE PROCEDURE TwoSiteExpValG_ct, TwoSiteExpValG_cm
END INTERFACE  TwoSiteExpValG


INTERFACE FormSingleSiteRho
	MODULE PROCEDURE FormSingleSiteRho_t, FormSingleSiteRho_to, FormSingleSiteRho_m, FormSingleSiteRho_mo
END INTERFACE FormSingleSiteRho

INTERFACE SingleSiteDensityMatrix
	MODULE PROCEDURE SingleSiteDensityMatrix_t, SingleSiteDensityMatrix_m
END INTERFACE SingleSiteDensityMatrix

INTERFACE TotalOneSite
	MODULE PROCEDURE  TotalOneSite_ct, TotalOneSite_cm
END INTERFACE TotalOneSite
	
INTERFACE GKernel
	MODULE PROCEDURE GKernel_t, GKernel_m
END INTERFACE GKernel

INTERFACE GNext
	MODULE PROCEDURE GNext_t, GNext_m, GNext_mt
END INTERFACE GNext

INTERFACE GContraction
	MODULE PROCEDURE GContraction_t,GContraction_to, GContraction_m,GContraction_mo
END INTERFACE GContraction

INTERFACE TotalEnergy
	MODULE PROCEDURE TotalEnergy_t,TotalEnergy_m
END INTERFACE TotalEnergy

INTERFACE LocalEnergy
	MODULE PROCEDURE LocalEnergy_t, LocalEnergy_m
END INTERFACE LocalEnergy

INTERFACE LocalEntropyDist
	MODULE PROCEDURE LocalEntropyDist_t, LocalEntropyDist_m
END INTERFACE LocalEntropyDist

INTERFACE MeyerQmeasure
	MODULE PROCEDURE MeyerQmeasure_t, MeyerQmeasure_m
END INTERFACE MeyerQmeasure

INTERFACE EvaluateMeasures
	MODULE PROCEDURE EvaluateMeasures_t, EvaluateMeasures_m
END INTERFACE EvaluateMeasures

INTERFACE InnerProduct
	MODULE PROCEDURE InnerProduct_t, InnerProduct_m
END INTERFACE InnerProduct

INTERFACE AllocateOverlap
	MODULE PROCEDURE AllocateOverlap_t,AllocateOverlap_tl, AllocateOverlap_m
END INTERFACE AllocateOverlap


CONTAINS


SUBROUTINE FormSingleSiteRho_to(rho1,  Gamma1, Lambda1)
!
!Purpose: Calculate the single site density matrix on the site where Gamma resides, store in rho1
!
IMPLICIT NONE
COMPLEX(KIND=rKind) :: rho1(:,:)
REAL(KIND=rKind) ::  Lambda1(:)
COMPLEX(KIND=rKind) :: Gamma1(:,:,:)
INTEGER :: chi0, chi1, i, j, alpha, beta

chi0 = SIZE(Gamma1,1)
chi1 = SIZE(Gamma1,3)
rho1=0.0_rKind
DO alpha=1,chi0
	DO beta=1,chi1
		IF(Lambda1(beta).gt.1E-15) THEN
		DO i=1,SIZE(Gamma1,2)
			rho1(i,:) = rho1(i,:)+Gamma1(alpha,i,beta)*Lambda1(beta) &
							*Lambda1(beta)*CONJG(Gamma1(alpha,:,beta))
		END DO
		END IF
	END DO
END DO
END SUBROUTINE FormSingleSiteRho_to

SUBROUTINE FormSingleSiteRho_t(rho1, Lambda0,Gamma1, Lambda1)
!
!Purpose: Calculate the single site density matrix on the site where Gamma resides, store in rho1
!
IMPLICIT NONE
COMPLEX(KIND=rKind) :: rho1(:,:)
REAL(KIND=rKind) :: Lambda0(:), Lambda1(:)
COMPLEX(KIND=rKind) :: Gamma1(:,:,:)
INTEGER :: chi0, chi1, i, j, alpha, beta

chi0 = SIZE(Gamma1,1)
chi1 = SIZE(Gamma1,3)
rho1=0.0_rKind
DO alpha=1,chi0
	IF(Lambda0(alpha).gt.1E-15) THEN
	DO beta=1,chi1
		IF(Lambda1(beta).gt.1E-15) THEN
		DO i=1,SIZE(Gamma1,2)
			rho1(i,:) = rho1(i,:)+Lambda0(alpha)*Gamma1(alpha,i,beta) &
							*CONJG(Gamma1(alpha,:,beta))*Lambda0(alpha)
		END DO
		END IF
	END DO
	END IF
END DO
END SUBROUTINE FormSingleSiteRho_t

SUBROUTINE FormSingleSiteRho_m(rho1, Lambda0, Gamma1, Lambda1, LabelLeft0,LabelLeft1)
!
!Purpose: Calculate the single site density matrix on the site where Gamma resides, store in rho1
!
IMPLICIT NONE
COMPLEX(KIND=rKind) :: rho1(:,:)
REAL(KIND=rKind) :: Lambda0(:), Lambda1(:)
COMPLEX(KIND=rKind) :: Gamma1(:,:)
INTEGER :: LabelLeft0(:), LabelLeft1(:)
INTEGER :: chi0, chi1, i, j, alpha, beta, loci

chi0 = SIZE(Gamma1,1)
chi1 = SIZE(Gamma1,2)
rho1=0.0_rKind
DO alpha=1,chi0
	IF(LabelLeft0(alpha).lt.10000) THEN
	DO beta=1,chi1
		IF(LabelLeft1(beta).lt.10000) THEN
		loci=LabelLeft1(beta)-LabelLeft0(alpha)+1
		IF((1.le.loci).and.(localSize.ge.loci)) THEN
		rho1(loci,loci) = rho1(loci,loci)+Lambda0(alpha)*Gamma1(alpha,beta) &
							*CONJG(Gamma1(alpha,beta))*Lambda0(alpha)
		END IF
		END IF
	END DO
	END IF
END DO
END SUBROUTINE FormSingleSiteRho_m

SUBROUTINE FormSingleSiteRho_mo(rho1,  Gamma1, Lambda1, LabelLeft0,LabelLeft1)
!
!Purpose: Calculate the single site density matrix on the site where Gamma resides, store in rho1
!
IMPLICIT NONE
COMPLEX(KIND=rKind) :: rho1(:,:)
REAL(KIND=rKind) ::  Lambda1(:)
COMPLEX(KIND=rKind) :: Gamma1(:,:)
INTEGER :: LabelLeft0(:), LabelLeft1(:)
INTEGER :: chi0, chi1, i, j, alpha, beta, loci

chi0 = SIZE(Gamma1,1)
chi1 = SIZE(Gamma1,2)
rho1=0.0_rKind
DO alpha=1,chi0
	IF(LabelLeft0(alpha).lt.10000) THEN
	DO beta=1,chi1
		IF(LabelLeft1(beta).lt.10000) THEN
		loci=LabelLeft1(beta)-LabelLeft0(alpha)+1
		IF((1.le.loci).and.(localSize.ge.loci)) THEN
		rho1(loci,loci) = rho1(loci,loci)+Gamma1(alpha,beta)*Lambda1(beta) &
							*Lambda1(beta)*CONJG(Gamma1(alpha,beta))
		END IF
		END IF
	END DO
	END IF
END DO
END SUBROUTINE FormSingleSiteRho_mo
		
SUBROUTINE SingleSiteDensityMatrix_t(rho,Gammas,Lambdas)
!
!Purpose: Calculate the single site density matrix on every site, store in rho
!
IMPLICIT NONE
TYPE(matrix), INTENT(OUT) :: rho(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER :: i

DO i=1,systemSize
	IF(i==1) THEN
		CALL FormSingleSiteRho(rho(i)%m, Gammas(i)%t, Lambdas(i+1)%v)
	ELSE
		CALL FormSingleSiteRho(rho(i)%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
	END IF
END DO
END SUBROUTINE SingleSiteDensityMatrix_t

SUBROUTINE SingleSiteDensityMatrix_m(rho,Gammas,Lambdas, LabelLeft)
!
!Purpose: Calculate the single site density matrix on every site, store in rho
!
IMPLICIT NONE
TYPE(matrix), INTENT(OUT) :: rho(:)
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
INTEGER :: i

DO i=1,systemSize
	IF(i==1) THEN
		CALL FormSingleSiteRho(rho(i)%m,  Gammas(i)%m, Lambdas(i+1)%v,LabelLeft(i)%vi, LabelLeft(i+1)%vi)
	ELSE
		CALL FormSingleSiteRho(rho(i)%m, Lambdas(i)%v, Gammas(i)%m, Lambdas(i+1)%v,LabelLeft(i)%vi, LabelLeft(i+1)%vi)
	END IF
END DO
END SUBROUTINE SingleSiteDensityMatrix_m


!!!!!!!!!!!!!!!!!!!!!BEGIN CONTENTS OF INTERFACE OneSiteExpVal!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OneSiteExpVal_c(expList,Op, Gammas, Lambdas)
!
!Purpose: Calculate expectation of the complex operator Op on every site
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: expList(systemSize)
COMPLEX(KIND=rKIND), INTENT(IN) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate rho in OneSiteExpVal_c'
END IF
		
DO i=1,systemSize,1
	rho%m=0.0_rKind
	expList(i)=0.0_rKind
	IF(i==1) THEN
		CALL FormSingleSiteRho(rho%m, Gammas(i)%t, Lambdas(i+1)%v)	
	ELSE
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
	END IF
	expList(i)=TraceMatmul(Op,rho%m)
END DO
				
DEALLOCATE(rho%m, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate rho in OneSiteExpVal_c'
END IF
END SUBROUTINE OneSiteExpVal_c


SUBROUTINE OneSiteExpVal_mc(expList,Op, Gammas, Lambdas, LabelLeft)
!
!Purpose: Calculate expectation of the complex operator Op on every site
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: expList(systemSize)
COMPLEX(KIND=rKIND), INTENT(IN) :: Op(:,:)
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
TYPE(matrix) :: rho
INTEGER :: i,j

ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate rho in OneSiteExpVal_r'
END IF
		
DO i=1,systemSize,1
	rho%m=0.0_rKind
	expList(i)=0.0_rKind
	IF(i==1) THEN
		CALL FormSingleSiteRho(rho%m, Gammas(i)%m, Lambdas(i+1)%v, LabelLeft(i)%vi, LabelLeft(i+1)%vi)	
	ELSE
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%m, Lambdas(i+1)%v, LabelLeft(i)%vi, LabelLeft(i+1)%vi)	
	END IF
	expList(i)=REAL(TraceMatmul(Op,rho%m),KIND=rKind)
END DO
				
DEALLOCATE(rho%m, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate rho in OneSiteExpVal_r'
END IF

END SUBROUTINE OneSiteExpVal_mc
!!!!!!!!!!!!!!!!!!!!!END CONTENTS OF INTERFACE OneSiteExpVal!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TotalOneSite_ct(total,Op, Gammas, Lambdas)
!
!Purpose: Calculate the total number value of Op across all lattice sites: complex version
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: total
COMPLEX(KIND=rKind) ,INTENT(IN) :: Op(:,:)
COMPLEX(KIND=rKind) :: expList(systemSize)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
COMPLEX(KIND=rKind) :: dumNum
INTEGER :: i

CALL OneSIteExpVal(expList,Op, Gammas, Lambdas)
total=SUM(expList)

END SUBROUTINE TotalOneSite_ct


SUBROUTINE TotalOneSite_cm(total,Op, Gammas, Lambdas,LabelLeft)
!
!Purpose: Calculate the total number value of Op across all lattice sites: complex version
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: total
COMPLEX(KIND=rKind) ,INTENT(IN) :: Op(:,:)
COMPLEX(KIND=rKind) :: expList(systemSize)
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
TYPE(matrix) :: rho
COMPLEX(KIND=rKind) :: dumNum
INTEGER :: i

CALL OneSiteExpVal_mc(expList,Op, Gammas, Lambdas, LabelLeft)
total=SUM(expList)


END SUBROUTINE TotalOneSite_cm

SUBROUTINE Qdepletion(depletion, rho, population, CW)
!
!Purpose: Calculate the Quantum depletion
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT), OPTIONAL :: CW(:)
COMPLEX(KIND=rKind), INTENT(IN) :: rho(:,:)
REAL(KIND=rKind), INTENT(IN) :: population
REAL(KIND=rKind), INTENT(OUT) :: depletion
REAL(KIND=rKind), ALLOCATABLE :: S(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: U(:,:), VT(:,:), work(:,:), rwork(:,:), rhoCopy(:,:)
INTEGER :: workSize, info, l1, l2
CHARACTER(1) :: jobu, jobvt

IF(population==0.0) THEN
depletion=0.0_rKind
	IF(PRESENT(CW)) THEN
		DO l1=1,systemSize
			CW(l1)=0.0_rKind
		END DO
	END IF
ELSE

!Allocate the variables needed to perform an SVD on rho		
	jobu='A'
	jobvt='A'
	workSize=5*systemSize
	ALLOCATE(S(systemSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate S in Qdepletion'
			END IF
	ALLOCATE(U(systemSize,systemSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate U in Qdepletion'
			END IF
	ALLOCATE(VT(systemSize,systemSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate VT in Qdepletion'
			END IF
	ALLOCATE(rhoCopy(systemSize,systemSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rhoCopy in Qdepletion'
			END IF
	ALLOCATE(work(workSize,workSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate work in Qdepletion'
			END IF
	ALLOCATE(rwork(workSize,workSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rwork in Qdepletion'
			END IF
	DO l1=1,systemSize
		DO l2=1,systemSize
			rhoCopy(l1,l2)=rho(l1,l2)
		END DO
	END DO
	!Perform an SVD on rho
	CALL ZGESVD(jobu, jobvt, systemSize, systemSize, rhoCopy, systemSize, S, U, systemSize, VT, systemSize, &
				work, workSize, rwork, info)


	depletion=1.0_rKind-S(1)/population

IF(PRESENT(CW)) THEN
	DO l1=1,systemSize
		CW(l1)=SQRT(S(1))*VT(1,l1)
	END DO
END IF

	!Deallocate the unnecessary variables
	DEALLOCATE(U,VT,work,rwork,S, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SVD variables in Qdepletion'
			END IF
	DEALLOCATE(rhoCopy, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rhocopy in Qdepletion'
			END IF
END IF
END SUBROUTINE Qdepletion

REAL(KIND=rKind) FUNCTION MeyerQmeasure_t(Gammas, Lambdas)
!
!Purpose: Calculate the Meyer Q-measure (average local impurity)
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
COMPLEX(KIND=rKind) :: dumPur
INTEGER :: i
REAL(KIND=rKind) :: totalPurity

ALLOCATE(rho%m(localSize, localSize), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate rho in MeyerQMeasure'
END IF
totalPurity = 0.0_rKind
! Calculate total purity, i.e., sum of tr(rho**2).
DO i = 1, systemSize
	IF(i==1) THEN
		CALL FormSingleSiteRho(rho%m, Gammas(i)%t, Lambdas(i+1)%v)
	ELSE
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
	END IF
	dumPur=TraceMatmul(rho%m, rho%m)
	totalPurity = totalPurity + REAL(dumPur, KIND=rKind)
END DO
MeyerQmeasure_t = localSize*1.0_rKind/(localSize*1.0_rKind-1.0_rKind) * (1.0_rKind - totalPurity/systemSize*1.0_rKind) ! Calculate average impurity.
DEALLOCATE(rho%m, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate rho in MeyerQMeasure'
END IF
END FUNCTION MeyerQmeasure_t

REAL(KIND=rKind) FUNCTION MeyerQmeasure_m(Gammas, Lambdas,LabelLeft)
!
!Purpose: Calculate the Meyer Q-measure (average local impurity)
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
TYPE(matrix) :: rho
COMPLEX(KIND=rKind) :: dumPur
INTEGER :: i
REAL(KIND=rKind) :: totalPurity

ALLOCATE(rho%m(localSize, localSize), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate rho in MeyerQMeasure'
END IF
totalPurity = 0.0_rKind
! Calculate total purity, i.e., sum of tr(rho**2).
DO i = 1, systemSize
	IF(i==1) THEN
		CALL FormSingleSiteRho(rho%m, Gammas(i)%m, Lambdas(i+1)%v,LabelLeft(i)%vi,LabelLeft(i+1)%vi)
	ELSE
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%m, Lambdas(i+1)%v,LabelLeft(i)%vi,LabelLeft(i+1)%vi)
	END IF
	dumPur=TraceMatmul(rho%m, rho%m)
	totalPurity = totalPurity + REAL(dumPur, KIND=rKind)
END DO
MeyerQmeasure_m = localSize*1.0_rKind/(localSize*1.0_rKind-1.0_rKind) * (1.0_rKind - totalPurity/systemSize*1.0_rKind) ! Calculate average impurity.
DEALLOCATE(rho%m, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate rho in MeyerQMeasure'
END IF
END FUNCTION MeyerQmeasure_m

REAL(KIND=rKind) FUNCTION ChainEntropy(link,Lambdas)
!
!Purpose: Calculate the entropy of entanglement of the chain to the left of link
!with the chain to the right on link (in the MPS approximation)
!
INTEGER, INTENT(IN) :: link
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind) :: temp
INTEGER :: chi, i

chi=SIZE(Lambdas(link)%v)
temp = 0.0_rKind
DO i=1,chi
	IF(Lambdas(link)%v(i).ne.0.0_rKind) THEN
		temp=temp-Lambdas(link)%v(i)*Lambdas(link)%v(i)*LOG(ABS(Lambdas(link)%v(i)*Lambdas(link)%v(i)))
	END IF
END DO
ChainEntropy=temp

END FUNCTION ChainEntropy

SUBROUTINE LocalEntropyDist_t(entDist, Gammas, Lambdas,tsalliSq)
!
!Purpose: Calculate the local entropy at each lattice site.
!If the OPTIONAL argument tsalliSq is present, compute the Tsallis entropy
!S_q=(1/(1-q))(Tr( rho^q) -1), else compute the von Neumann entropy
!S=-Tr( rho Log_d rho)
!
!See manual for more detail
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: entDist(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
REAL(KIND=rKind), INTENT(IN), OPTIONAL :: tsallisQ
INTEGER :: i, j, workSize, info
REAL(KIND=rKind) :: evals(localSize), temp
COMPLEX(KIND=rKind), ALLOCATABLE :: workArr(:), rworkArr(:)


!Allocate single-site density matrix
ALLOCATE(rho%m(localSize, localSize), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate rho in LocalEntropyDist'
END IF
!Allocate workspace for ZHEEV
workSize = 2*localSize-1
ALLOCATE(workarr(workSize), rworkarr(3*localSize-2), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate workarr in LocalEntropyDist'
END IF

DO i = 1, systemSize
	temp = 0.0_rKind
!Form single-site density matrix
	CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
!Diagonalize density matrix
	CALL ZHEEV('N', 'U', localSize, rho%m, localSize, evals, workArr, workSize, rworkArr, info) ! Diagonalize single site density matrix.

IF(PRESENT(tsalliSq)) THEN
	!Compute the Tsallis entropy
	temp=-(1.0_rKind/(1.0_rKind-tsallisQ))
	DO j=1,localSize
		IF(evals(j).ne.0.0_rKind) THEN
			temp=temp+(1.0_rKind/(1.0_rKind-tsalliSq))*(ABS(evals(j))**tsalliSq)
		ELSE
			temp=temp+0.0_rKind
		END IF
	END DO
	entDist(i) = temp
ELSE
!Compute vonNeumann entropy
	DO j = 1, localSize
		IF(evals(j).ne.0.0_rKind) THEN
		        temp = temp + evals(j)*LOG(ABS(evals(j)))/LOG(1.0_rKind*localSize) ! Compute tr(rho*log_d(rho)).
		ELSE
			temp=temp+0.0_rKind
		END IF

	END DO
	entDist(i) = -1.0_rKind*temp ! Set value of entropy on-site.
END IF

END DO

DEALLOCATE(rho%m, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate rho in LocalEntropyDist'
END IF
DEALLOCATE(workarr, rworkarr, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate variables in LocalEntropyDist'
END IF
END SUBROUTINE LocalEntropyDist_t

SUBROUTINE LocalEntropyDist_m(entDist, Gammas, Lambdas,LabelLeft,tsalliSq)
!
!Purpose: Calculate the local entropy at each lattice site.
!If the OPTIONAL argument tsalliSq is present, compute the Tsallis entropy
!S_q=(1/(1-q))(Tr( rho^q) -1), else compute the von Neumann entropy
!S=-Tr( rho Log_d rho)
!
!See manual for more detail
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: entDist(:)
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
TYPE(matrix) :: rho
REAL(KIND=rKind), INTENT(IN), OPTIONAL :: tsallisQ
INTEGER :: i, j, workSize, info
REAL(KIND=rKind) :: evals(localSize), temp
COMPLEX(KIND=rKind), ALLOCATABLE :: workArr(:), rworkArr(:)


!Allocate single-site density matrix
ALLOCATE(rho%m(localSize, localSize), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate rho in LocalEntropyDist'
END IF
!Allocate workspace for ZHEEV
workSize = 2*localSize-1
ALLOCATE(workarr(workSize), rworkarr(3*localSize-2), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate workarr in LocalEntropyDist'
END IF

DO i = 1, systemSize
	temp = 0.0_rKind
!Form single-site density matrix
	CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%m, Lambdas(i+1)%v, LabelLeft(i)%vi, LabelLeft(i+1)%vi)
!Diagonalize density matrix
	CALL ZHEEV('N', 'U', localSize, rho%m, localSize, evals, workArr, workSize, rworkArr, info) ! Diagonalize single site density matrix.

IF(PRESENT(tsalliSq)) THEN
	!Compute the Tsallis entropy
	temp=-(1.0_rKind/(1.0_rKind-tsallisQ))
	DO j=1,localSize
		IF(evals(j).ne.0.0_rKind) THEN
			temp=temp+(1.0_rKind/(1.0_rKind-tsalliSq))*(ABS(evals(j))**tsalliSq)
		ELSE
			temp=temp+0.0_rKind
		END IF
	END DO
	entDist(i) = temp
ELSE
!Compute vonNeumann entropy
	DO j = 1, localSize
		IF(evals(j).ne.0.0_rKind) THEN
		        temp = temp + evals(j)*LOG(ABS(evals(j)))/LOG(1.0_rKind*localSize) ! Compute tr(rho*log_d(rho)).
		ELSE
			temp=temp+0.0_rKind
		END IF

	END DO
	entDist(i) = -1.0_rKind*temp ! Set value of entropy on-site.
END IF

END DO

DEALLOCATE(rho%m, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate rho in LocalEntropyDist'
END IF
DEALLOCATE(workarr, rworkarr, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate variables in LocalEntropyDist'
END IF
END SUBROUTINE LocalEntropyDist_m

SUBROUTINE GKernel_t(gee,Gamma,GammaP)
!
!Purpose: First step in calculating the two-site observable
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix) :: gee
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:,:),GammaP(:,:,:)
INTEGER :: alpha, beta, betapr,i

ALLOCATE(gee%m(SIZE(Gamma,1),SIZE(Gammap,1)))
gee%m=0.0_rKind
!$OMP PARALLEL DO PRIVATE(beta, betapr,i) schedule(dynamic)
DO alpha=1,SIZE(Gamma,1)
	DO beta=1,SIZE(GammaP,1)
		DO betapr=1,SIZE(Gamma,3)
			DO i=1,SIZE(Gamma,2),1
			gee%m(alpha,beta)=gee%m(alpha,beta)+CONJG(GammaP(beta,i,betapr))*Gamma(alpha,i,betapr)
			END DO
		END DO
	END DO
END DO
!$OMP END PARALLEL DO
END SUBROUTINE GKernel_t

SUBROUTINE GKernel_m(gee,Gamma,GammaP, LabelLeft, LabelLeftp)
!
!Purpose: First step in calculating the two-site observable
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix) :: gee
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:),GammaP(:,:,:)
INTEGER, INTENT(IN) :: LabelLeft(:), LabelLeftp(:)
INTEGER :: alpha, beta, betapr,i, loci

ALLOCATE(gee%m(SIZE(Gamma,1),SIZE(Gammap,1)))
gee%m=0.0_rKind
!$OMP PARALLEL DO PRIVATE(beta,betapr,loci) schedule(dynamic)
DO alpha=1,SIZE(Gamma,1)
	IF(LabelLeft(alpha).lt.10000) THEN
	DO beta=1,SIZE(GammaP,1)
		IF(LabelLeft(beta).lt.10000) THEN
		DO betapr=1,SIZE(Gamma,2)
			IF(LabelLeftp(betapr).lt.10000) THEN
			loci=LabelLeftp(betapr)-LabelLeft(alpha)+1
				IF((1.le.loci).and.(loci.le.localSize)) THEN
					gee%m(alpha,beta)=gee%m(alpha,beta)+CONJG(GammaP(beta,loci,betapr))*Gamma(alpha,betapr)
				END IF
			END IF
		END DO
		END IF
	END DO
	END IF
END DO
!$OMP END PARALLEL DO
END SUBROUTINE GKernel_m


SUBROUTINE GNext_t(gee,Gamma,GammaP)
!
!Purpose:  Recursive step in calculating the two-site observable
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix) :: gee
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:,:),GammaP(:,:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: temp(:,:,:)
INTEGER :: alpha, beta, alphapr,betapr,i

ALLOCATE(temp(SIZE(Gamma,1),SIZE(Gamma,2),SIZE(gee%m,2)))
temp=0.0_rKind
!$OMP PARALLEL DO PRIVATE(alphapr,i) schedule(dynamic)
DO alpha=1,SIZE(Gamma,1)
	DO alphapr=1,SIZE(Gamma,3)
		DO i=1,SIZE(Gamma,2)
			temp(alpha,i,:)=temp(alpha,i,:)+Gamma(alpha,i,alphapr)*gee%m(alphapr,:)
		END DO
	END DO
END DO 
!$OMP END PARALLEL DO

DEALLOCATE(gee%m)
ALLOCATE(gee%m(SIZE(Gamma,1),SIZE(GammaP,1)))

gee%m=0.0_rKind
!$OMP PARALLEL DO PRIVATE(i, beta, betapr) schedule(dynamic)
DO beta=1,SIZE(GammaP,1)
	DO betapr=1,SIZE(GammaP,3)
		DO i=1,SIZE(Gammap,2)
			gee%m(:,beta)=gee%m(:,beta)+CONJG(GammaP(beta,i,betapr))*temp(:,i,betapr)
		END DO
	END DO
END DO
!$OMP END PARALLEL DO

DEALLOCATE(temp)

END SUBROUTINE GNext_t

SUBROUTINE GNext_mt(gee,Gamma,GammaP, LabelLeft, LabelLeftp)
!
!Purpose:  Recursive step in calculating the two-site observable
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix) :: gee
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:),GammaP(:,:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: temp(:,:,:)
INTEGER, INTENT(IN) :: LabelLeft(:), LabelLeftp(:)
INTEGER :: alpha, beta, alphapr,betapr,i, loci


ALLOCATE(temp(SIZE(Gamma,1),SIZE(GammaP,2),SIZE(gee%m,2)))
temp=0.0_rKind
!$OMP PARALLEL DO PRIVATE(alphapr,loci) schedule(dynamic)
DO alpha=1,SIZE(Gamma,1)
	IF(LabelLeft(alpha).lt.10000) THEN
	DO alphapr=1,SIZE(Gamma,2)
			IF(LabelLeftp(alphapr).lt.10000) THEN
			loci=LabelLeftp(alphapr)-LabelLeft(alpha)+1
			IF((1.le.loci).and.(loci.le.localSize)) THEN
				temp(alpha,loci,:)=temp(alpha,loci,:)+Gamma(alpha,alphapr)*gee%m(alphapr,:)
			END IF
			END IF
	END DO
	END IF
END DO 
!$OMP END PARALLEL DO

DEALLOCATE(gee%m)
ALLOCATE(gee%m(SIZE(Gamma,1),SIZE(GammaP,1)))

gee%m=0.0_rKind
!$OMP PARALLEL DO PRIVATE(i, beta, betapr) schedule(dynamic)
DO beta=1,SIZE(GammaP,1)
	IF(LabelLeft(beta).lt.10000) THEN
	DO betapr=1,SIZE(GammaP,3)
		IF(LabelLeftp(betapr).lt.10000) THEN
		DO i=1,SIZE(Gammap,2)
			gee%m(:,beta)=gee%m(:,beta)+CONJG(GammaP(beta,i,betapr))*temp(:,i,betapr)
		END DO
		END IF
	END DO
	END IF
END DO
!$OMP END PARALLEL DO
DEALLOCATE(temp)

END SUBROUTINE GNext_mt

SUBROUTINE GNext_m(gee,Gamma,GammaP, LabelLeft, LabelLeftp)
!
!Purpose:  Recursive step in calculating the two-site observable
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix) :: gee
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:),GammaP(:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: temp(:,:,:)
INTEGER, INTENT(IN) :: LabelLeft(:), LabelLeftp(:)
INTEGER :: alpha, beta, alphapr,betapr,i, loci


ALLOCATE(temp(SIZE(Gamma,1),localSize,SIZE(gee%m,2)))
temp=0.0_rKind
!$OMP PARALLEL DO PRIVATE(alphapr,loci) schedule(dynamic)
DO alpha=1,SIZE(Gamma,1)
	IF(LabelLeft(alpha).lt.10000) THEN
	DO alphapr=1,SIZE(Gamma,2)
			IF(LabelLeftp(alphapr).lt.10000) THEN
			loci=LabelLeftp(alphapr)-LabelLeft(alpha)+1
			IF((1.le.loci).and.(loci.le.localSize)) THEN
				temp(alpha,loci,:)=temp(alpha,loci,:)+Gamma(alpha,alphapr)*gee%m(alphapr,:)
			END IF
			END IF
	END DO
	END IF
END DO 
!$OMP END PARALLEL DO
DEALLOCATE(gee%m)
ALLOCATE(gee%m(SIZE(Gamma,1),SIZE(GammaP,1)))

gee%m=0.0_rKind
!$OMP PARALLEL DO PRIVATE(loci, beta, betapr) schedule(dynamic)
DO beta=1,SIZE(GammaP,1)
	IF(LabelLeft(beta).lt.10000) THEN
	DO betapr=1,SIZE(GammaP,2)
		IF(LabelLeftp(betapr).lt.10000) THEN
			loci=LabelLeftp(betapr)-LabelLeft(beta)+1
			IF((1.le.loci).and.(loci.le.localSize)) THEN			
			gee%m(:,beta)=gee%m(:,beta)+CONJG(GammaP(beta,betapr))*temp(:,loci,betapr)
			END IF
		END IF
	END DO
	END IF
END DO
!$OMP END PARALLEL DO
DEALLOCATE(temp)

END SUBROUTINE GNext_m


SUBROUTINE GContraction_t(obsv,gee,Gamma,GammaP,Lambda1,Lambda2)
!
!Purpose: Final step in calculating the two-site observable
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind) :: obsv
TYPE(matrix) :: gee
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:,:),GammaP(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda1(:),Lambda2(:)
INTEGER :: alpha, beta, betapr,i

obsv=0.0_rKind
!$OMP PARALLEL DO REDUCTION(+:obsv) PRIVATE( beta, betapr,i) schedule(dynamic)
DO alpha=1,SIZE(Gamma,3)
	DO beta=1,SIZE(Gamma,3)
		DO betapr=1,SIZE(Gamma,1)
			DO i=1,SIZE(Gamma,2),1
			obsv=obsv+(Lambda1(betapr)**2)*CONJG(GammaP(betapr,i,beta))*Gamma(betapr,i,alpha)*gee%m(alpha,beta)
			END DO
		END DO
	END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE GContraction_t

SUBROUTINE GContraction_to(obsv,gee,Gamma,GammaP,Lambda2)
!
!Purpose: Final step in calculating the two-site observable
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind) :: obsv
TYPE(matrix) :: gee
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:,:),GammaP(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda2(:)
INTEGER :: alpha, beta, betapr,i

obsv=0.0_rKind
!$OMP PARALLEL DO REDUCTION(+:obsv) PRIVATE( beta, betapr,i) schedule(dynamic)
DO alpha=1,SIZE(Gamma,3)
	DO beta=1,SIZE(Gamma,3)
		DO betapr=1,SIZE(Gamma,1)
			DO i=1,SIZE(Gamma,2),1
			obsv=obsv+CONJG(GammaP(betapr,i,beta))*Gamma(betapr,i,alpha)*&
			Lambda2(alpha)*Lambda2(beta)*gee%m(alpha,beta)
			END DO
		END DO
	END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE GContraction_to

SUBROUTINE GContraction_m(obsv,gee,Gamma,GammaP,Lambda1,Lambda2,LabelLeft1,LabelLeft2)
!
!Purpose: Final step in calculating the two-site observable
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind) :: obsv
TYPE(matrix) :: gee
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:),GammaP(:,:,:)
INTEGER, INTENT(IN) :: LabelLeft1(:), LabelLeft2(:)
REAL(KIND=rKind), INTENT(IN) :: Lambda1(:),Lambda2(:)
INTEGER :: alpha, beta, betapr,i, loci

obsv=0.0_rKind
!$OMP PARALLEL DO REDUCTION(+:obsv) PRIVATE( beta, betapr,loci) schedule(dynamic)
DO alpha=1,SIZE(Gamma,2)
	IF(LabelLeft2(alpha).lt.10000) THEN
	DO beta=1,SIZE(Gamma,2)
		IF(LabelLeft2(beta).lt.10000) THEN
		DO betapr=1,SIZE(Gamma,1)
			IF(LabelLeft1(betapr).lt.10000) THEN
				loci=LabelLeft2(alpha)-LabelLeft1(betapr)+1
				IF((1.le.loci).and.(loci.le.localSize)) THEN
					obsv=obsv+(Lambda1(betapr)**2)*CONJG(GammaP(betapr,loci,beta))*Gamma(betapr,alpha)*gee%m(alpha,beta)
				END IF
			END IF
		END DO
		END IF
	END DO
	END IF
END DO
!$OMP END PARALLEL DO

END SUBROUTINE GContraction_m

SUBROUTINE GContraction_mo(obsv,gee,Gamma,GammaP,Lambda2,LabelLeft1,LabelLeft2)
!
!Purpose: Final step in calculating the two-site observable
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind) :: obsv
TYPE(matrix) :: gee
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:),GammaP(:,:,:)
INTEGER, INTENT(IN) :: LabelLeft1(:), LabelLeft2(:)
REAL(KIND=rKind), INTENT(IN) :: Lambda2(:)
INTEGER :: alpha, beta, betapr,i, loci

obsv=0.0_rKind
!$OMP PARALLEL DO REDUCTION(+:obsv) PRIVATE( beta, betapr,loci) schedule(dynamic)
DO alpha=1,SIZE(Gamma,2)
	IF(LabelLeft2(alpha).lt.10000) THEN
	DO beta=1,SIZE(Gamma,2)
		IF(LabelLeft2(beta).lt.10000) THEN
		DO betapr=1,SIZE(Gamma,1)
			IF(LabelLeft1(betapr).lt.10000) THEN
				loci=LabelLeft2(alpha)-LabelLeft1(betapr)+1
				IF((1.le.loci).and.(loci.le.localSize)) THEN
					obsv=obsv+CONJG(GammaP(betapr,loci,beta))*Gamma(betapr,alpha)*&
					Lambda2(alpha)*Lambda2(beta)*gee%m(alpha,beta)
				END IF
			END IF
		END DO
		END IF
	END DO
	END IF
END DO
!$OMP END PARALLEL DO

END SUBROUTINE GContraction_mo

!!!!!!!!!!!!!!!!!!! BEGIN CONTENTS OF INTERFACE TwoSiteExpValG !!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TwoSiteExpValG_ct(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
!
!Purpose: Calculate the expectation value of the complex two-site operator Op1XOp2 at every pair of sites
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)		
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: gee
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2, i,alpha, beta, betapr
INTEGER, INTENT(IN), OPTIONAL :: phaseStat

DO l1=1,systemSize
	!On-diagonal elements use Op1
	IF(l1==1) THEN
		CALL FormSingleSiteRho(rho1,  Gammas(l1)%t, Lambdas(l1+1)%v)
	ELSE
		CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
	END IF
	observable(l1,l1)=TraceMatmul(MATMUL(Op1,Op2),rho1)
END DO				
		
DO l2=systemSize,2,-1
	ALLOCATE(GammaP(SIZE(Gammas(l2)%t,1),SIZE(Gammas(l2)%t,2),SIZE(Gammas(l2)%t,3)))		
	GammaP = Gammas(l2)%t
	!Compute the initial G 
	CALL OneSiteOp(Op2,GammaP)
	CALL GKernel(gee,Gammas(l2)%t,GammaP)
	DEALLOCATE(GammaP)
	DO l1=(l2-1),1,(-1)
		ALLOCATE(GammaP(SIZE(Gammas(l1)%t,1),SIZE(Gammas(l1)%t,2),SIZE(Gammas(l1)%t,3)))		
		GammaP = Gammas(l1)%t					
		!Fermi Phase for final
		IF(PRESENT(phaseStat)) THEN
			CALL OneSiteOp(fermiPhase_op%m,GammaP)
		END IF
		!Compute final G
		CALL OneSiteOp(Op1,GammaP)
		IF(l1==1) THEN
			CALL GContraction(observable(l1,l2),gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)
		ELSE
			CALL GContraction(observable(l1,l2),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
		END IF
		observable(l2,l1)=CONJG(observable(l1,l2))
		GammaP = Gammas(l1)%t
		!Fermi Phase for next
		IF(PRESENT(phaseStat)) THEN
			CALL OneSiteOp(fermiPhase_op%m,GammaP)
		END IF
		!Next G
		CALL GNext(gee,Gammas(l1)%t,GammaP)
		DEALLOCATE(GammaP)
	END DO
	DEALLOCATE(gee%m)
END DO

END SUBROUTINE TwoSiteExpValG_ct

SUBROUTINE TwoSiteExpValG_cm(observable, Op1, Op2, Gammas, Lambdas,LabelLeft, phaseStat)
!
!Purpose: Calculate the expectation value of the real two-site operator Op1XOp2 at every pair of sites
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)		
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: gee
TYPE(tensor) :: GammaP
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2, i,alpha, beta, betapr
INTEGER, INTENT(IN), OPTIONAL :: phaseStat


DO l1=1,systemSize
	!On-diagonal elements use Op1
	IF(l1==1) THEN
		CALL FormSingleSiteRho(rho1,  Gammas(l1)%m, Lambdas(l1+1)%v,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
	ELSE
		CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%m, Lambdas(l1+1)%v,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
	END IF
	observable(l1,l1)=TraceMatmul(MATMUL(Op1,Op2),rho1)
END DO				
		
DO l2=systemSize,2,-1
	!Compute the initial G 
	CALL OneSiteOp(Op2,Gammas(l2)%m,GammaP,LabelLeft(l2)%vi,LabelLeft(l2+1)%vi)
	CALL GKernel(gee,Gammas(l2)%m,GammaP%t,LabelLeft(l2)%vi,LabelLeft(l2+1)%vi)
	DEALLOCATE(GammaP%t)
	DO l1=(l2-1),1,(-1)
		!Fermi Phase for final
		IF(PRESENT(phaseStat)) THEN
			CALL OneSiteOp(fermiPhase_op%m,Gammas(l1)%m,GammaP,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
			!Compute final G
			CALL OneSiteOp(Op1,GammaP%t)
		IF(l1==1) THEN
			CALL GContraction(observable(l1,l2),gee,Gammas(l1)%m,GammaP%t,Lambdas(l1+1)%v,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
		ELSE
			CALL GContraction(observable(l1,l2),gee,Gammas(l1)%m,GammaP%t,Lambdas(l1)%v,Lambdas(l1+1)%v,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
		END IF
			observable(l2,l1)=CONJG(observable(l1,l2))
			DEALLOCATE(GammaP%t)
			CALL OneSiteOp(fermiPhase_op%m,Gammas(l1)%m,GammaP,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
			!Next G
			CALL GNext(gee,Gammas(l1)%m,GammaP%t,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
			DEALLOCATE(GammaP%t)
		ELSE

			!Compute final G
			CALL OneSiteOp(Op1,Gammas(l1)%m,GammaP,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
		IF(l1==1) THEN
			CALL GContraction(observable(l1,l2),gee,Gammas(l1)%m,GammaP%t,Lambdas(l1+1)%v,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
		ELSE
			CALL GContraction(observable(l1,l2),gee,Gammas(l1)%m,GammaP%t,Lambdas(l1)%v,Lambdas(l1+1)%v,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
		END IF
			DEALLOCATE(GammaP%t)
			observable(l2,l1)=CONJG(observable(l1,l2))
			!Next G
			CALL GNext(gee,Gammas(l1)%m,Gammas(l1)%m,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
		END IF
	END DO
	DEALLOCATE(gee%m)
END DO

END SUBROUTINE TwoSiteExpValG_cm

!!!!!!!!!!!!!!!!!!! END CONTENTS OF INTERFACE TwoSiteExpValG !!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE LocalEnergy_t(energy, mpo, Gammas, Lambdas)
!
!Purpose: Calculate the energy associated with each lattice bond.
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: energy(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(DecomposedMPO), POINTER :: MPO(:)
TYPE(tensor) :: GammaP
TYPE(matrix) :: gee
COMPLEX(KIND=rKind) :: rho2(localSize*localSize, localSize*localSize)
COMPLEX(KIND=rKind) :: dumEn
INTEGER :: i, l1, l2, numOps	
energy = 0.0_rKind

DO l2=systemSize,2,-1
	numops=SIZE(MPO(l2-1)%l)
	DO i=1,numOps
	!Compute the initial G 
	ALLOCATE(GammaP%t(SIZE(Gammas(l2)%t,1),SIZE(Gammas(l2)%t,2),SIZE(Gammas(l2)%t,3)))
	GammaP%t=Gammas(l2)%t
	CALL OneSiteOp(MPO(l2-1)%r(i)%m,GammaP%t)
	CALL GKernel(gee,Gammas(l2)%t,GammaP%t)
	DEALLOCATE(GammaP%t)
	l1=(l2-1)
	!Compute final G
	ALLOCATE(GammaP%t(SIZE(Gammas(l1)%t,1),SIZE(Gammas(l1)%t,2),SIZE(Gammas(l1)%t,3)))
	GammaP%t=Gammas(l1)%t
	CALL OneSiteOp(MPO(l2-1)%l(i)%m,GammaP%t)
	IF(l1==1) THEN
		CALL GContraction(dumEn,gee,Gammas(l1)%t,GammaP%t,Lambdas(l1+1)%v)
	ELSE
		CALL GContraction(dumEn,gee,Gammas(l1)%t,GammaP%t,Lambdas(l1)%v,Lambdas(l1+1)%v)
	END IF
	DEALLOCATE(GammaP%t)
	energy(l1) = energy(l1)+REAL(dumEn, KIND=rKind)
	DEALLOCATE(gee%m)
	END DO
END DO

END SUBROUTINE LocalEnergy_t


SUBROUTINE LocalEnergy_m(energy, MPO, Gammas, Lambdas, LabelLeft)
!
!Purpose: Calculate the energy associated with each lattice bond.
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: energy(:)
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(DecomposedMPO), POINTER :: MPO(:)
TYPE(tensor) :: GammaP
TYPE(matrix) :: gee
COMPLEX(KIND=rKind) :: rho2(localSize*localSize, localSize*localSize)
COMPLEX(KIND=rKind) :: dumEn
INTEGER :: i, l1, l2, numOps	
energy = 0.0_rKind


DO l2=systemSize,2,-1
	numops=SIZE(MPO(l2-1)%l)
	DO i=1,numOps
	!Compute the initial G 
	CALL OneSiteOp(MPO(l2-1)%r(i)%m,Gammas(l2)%m,GammaP,LabelLeft(l2)%vi,LabelLeft(l2+1)%vi)
	CALL GKernel(gee,Gammas(l2)%m,GammaP%t,LabelLeft(l2)%vi,LabelLeft(l2+1)%vi)
	DEALLOCATE(GammaP%t)
	l1=(l2-1)
	!Compute final G
	CALL OneSiteOp(MPO(l2-1)%l(i)%m,Gammas(l1)%m,GammaP,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
	IF(l1==1) THEN
		CALL GContraction(dumEn,gee,Gammas(l1)%m,GammaP%t,Lambdas(l1+1)%v,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
	ELSE
		CALL GContraction(dumEn,gee,Gammas(l1)%m,GammaP%t,Lambdas(l1)%v,Lambdas(l1+1)%v,LabelLeft(l1)%vi,LabelLeft(l1+1)%vi)
	END IF
	DEALLOCATE(GammaP%t)
	energy(l1) = energy(l1)+REAL(dumEn, KIND=rKind)

	DEALLOCATE(gee%m)
	END DO
END DO


END SUBROUTINE LocalEnergy_m

SUBROUTINE TotalEnergy_t(energy, MPO, Gammas, Lambdas)
!
!Purpose: Calculate the energy eigenvalue associated with the Hamiltonian H
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: energy
REAL(KIND=rKind), ALLOCATABLE :: locEnergy(:)
TYPE(DecomposedMPO), POINTER :: mpo(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)

ALLOCATE(locEnergy(systemSize-1))
CALL LocalEnergy(locEnergy, mpo, Gammas, Lambdas)
energy=SUM(locEnergy)
DEALLOCATE(locEnergy)

END SUBROUTINE TotalEnergy_t


SUBROUTINE TotalEnergy_m(energy, MPO, Gammas, Lambdas,LabelLeft)
!
!Purpose: Calculate the energy eigenvalue associated with the Hamiltonian H
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: energy
REAL(KIND=rKind), ALLOCATABLE :: locEnergy(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
TYPE(DecomposedMPO), POINTER :: mpo(:)
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)

ALLOCATE(locEnergy(systemSize-1))
CALL LocalEnergy(locEnergy, mpo, Gammas, Lambdas,LabelLeft)
energy=SUM(locEnergy)
DEALLOCATE(locEnergy)

END SUBROUTINE TotalEnergy_m

COMPLEX(KIND=rKind) FUNCTION InnerProduct_t(GammasL, LambdasL, GammasR, LambdasR)
!
!Purpose: Calculate the inner product of the wavefunction by the Gammas and Lambdas.
!This routine can also be used to compute the Fidelity or Loschmidt echo if GammasL/LambdasL
!is some initial state and GammasR/LambdasR some time evolved final state
!
IMPLICIT NONE
TYPE(tensor), POINTER :: GammasL(:), GammasR(:)
TYPE(vector), POINTER :: LambdasL(:), LambdasR(:)
COMPLEX(KIND=rKind) :: temp
COMPLEX(KIND=rKind), ALLOCATABLE :: normKernel(:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaTemp(:,:,:)
INTEGER :: alpha,beta,gamma,eta,i,n,chi,chil1,chil2,chir1,chir2, locSyssize, loclocSize
locSyssize=SIZE(GammasL)


chir1=SIZE(GammasL(1)%t,3)
chir2=SIZE(GammasR(1)%t,3)
loclocSize=SIZE(GammasR(1)%t,2)
ALLOCATE(normKernel(chir1,chir2))
DO alpha=1,chir1
	DO beta=1,chir2
		normKernel(alpha,beta)=CMPLX(0.0,KIND=rKind)
		DO i=1,loclocSize
			normKernel(alpha,beta)=normKernel(alpha,beta)+LambdasL(2)%v(alpha)*CONJG(GammasL(1)%t(1,i,alpha)) &
						*GammasR(1)%t(1,i,beta)*LambdasR(2)%v(beta)
		END DO
	END DO
END DO
DO n=2,(locSyssize-1)
	chil1=SIZE(GammasL(n)%t,1)
	chil2=SIZE(GammasR(n)%t,1)
	chir1=SIZE(GammasL(n)%t,3)
	chir2=SIZE(GammasR(n)%t,3)
	loclocSize=SIZE(GammasR(n)%t,2)
	ALLOCATE(GammaTemp(chil1,loclocSize,chir2))
	GammaTemp=0.0_rKind

	DO beta=1,chir2
		DO i=1,loclocSize
			DO eta=1,chil2
				GammaTemp(:,i,beta) = GammaTemp(:,i,beta) &
			    + normKernel(:,eta)*GammasR(n)%t(eta,i,beta)
			END DO
		END DO
	END DO
	DEALLOCATE(normKernel)
	ALLOCATE(normKernel(chir1,chir2))
	normKernel=0.0_rKind
	DO alpha=1,chir1
		DO beta=1,chir2
			DO gamma=1,chil1
				DO i=1,SIZE(GammasL(n)%t,2)
					normKernel(alpha,beta) = normKernel(alpha,beta) &
						   + CONJG(GammasL(n)%t(gamma,i,alpha)) &
						   * GammaTemp(gamma,i,beta)
				END DO
			END DO
		END DO
	END DO
	DEALLOCATE(GammaTemp)
END DO
temp=CMPLX(0.0,KIND=rKind);

DO alpha=1,SIZE(GammasL(locsysSize)%t,1)
	DO beta=1,SIZE(GammasR(locsysSize)%t,1)
		DO i=1,SIZE(GammasL(locSyssize)%t,2)
			temp = temp+CONJG(GammasL(locSyssize)%t(alpha,i,1))*normKernel(alpha,beta) &
				*GammasR(locSyssize)%t(beta,i,1)
		END DO
	END DO
END DO

DEALLOCATE(normKernel)
InnerProduct_t=temp

END FUNCTION InnerProduct_t

COMPLEX(KIND=rKind) FUNCTION InnerProduct_m(GammasL, LambdasL,LabelL, GammasR, LambdasR,LabelR)
!
!Purpose: Calculate the inner product of the wavefunction by the Gammas and Lambdas.
!This routine can also be used to compute the Fidelity or Loschmidt echo if GammasL/LambdasL
!is some initial state and GammasR/LambdasR some time evolved final state
!
IMPLICIT NONE
TYPE(matrix), POINTER :: GammasL(:), GammasR(:)
TYPE(vectorInt), POINTER :: LabelL(:), LabelR(:)
TYPE(vector), POINTER :: LambdasL(:), LambdasR(:)
COMPLEX(KIND=rKind) :: temp
COMPLEX(KIND=rKind), ALLOCATABLE :: normKernal(:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaTemp(:,:,:)
INTEGER :: alpha,beta,gamma,eta,i,n,chi,chil1,chil2,chir1,chir2, locSyssize, loclocSize, loci, locj
locSyssize=SIZE(GammasL)


chir1=SIZE(GammasL(1)%m,2)
chir2=SIZE(GammasR(1)%m,2)
loclocSize=localSize
ALLOCATE(normKernal(chir1,chir2))
normKernal=0.0_rKind
DO alpha=1,chir1
	IF(LabelL(2)%vi(alpha).lt.10000) THEN
	DO beta=1,chir2
		IF(LabelR(2)%vi(beta).lt.10000) THEN
			IF(LabelL(2)%vi(alpha)==LabelR(2)%vi(beta)) THEN
				normKernal(alpha,beta)=normKernal(alpha,beta)+LambdasL(2)%v(alpha)*CONJG(GammasL(1)%m(1,alpha)) &
							*GammasR(1)%m(1,beta)*LambdasR(2)%v(beta)
			END IF
		END IF
	END DO
	END IF
END DO
DO n=2,(locSyssize-1)
	chil1=SIZE(GammasL(n)%m,1)
	chil2=SIZE(GammasR(n)%m,1)
	chir1=SIZE(GammasL(n)%m,2)
	chir2=SIZE(GammasR(n)%m,2)
	ALLOCATE(GammaTemp(chil1,localSize,chir2))
	GammaTemp=0.0_rKind
	DO beta=1,chir2
		IF(LabelR(n+1)%vi(beta).lt.10000) THEN
		DO eta=1,chil2
			IF(LabelR(n)%vi(eta).lt.10000) THEN
			loci=LabelR(n+1)%vi(beta)-LabelR(n)%vi(eta)+1
			IF((1.le.loci).and.(localSize.ge.loci)) THEN
			GammaTemp(:,loci,beta) = GammaTemp(:,loci,beta) &
		    + normKernal(:,eta)*GammasR(n)%m(eta,beta)
			END IF
			END IF
		END DO
		END IF
	END DO
	DEALLOCATE(normKernal)
	ALLOCATE(normKernal(chir1,chir2))
	normKernal=0.0_rKind
	DO alpha=1,chir1
		IF(LabelL(n+1)%vi(alpha).lt.10000) THEN
		DO beta=1,chir2
			IF(LabelR(n+1)%vi(beta).lt.10000) THEN
			DO gamma=1,chil1
				IF(LabelL(n)%vi(gamma).lt.10000) THEN
					loci=LabelL(n+1)%vi(alpha)-LabelL(n)%vi(gamma)+1
					IF((1.le.loci).and.(localSize.ge.loci)) THEN
					normKernal(alpha,beta) = normKernal(alpha,beta) &
						   + CONJG(GammasL(n)%m(gamma,alpha)) &
						   * GammaTemp(gamma,loci,beta)
					END IF
				END IF
			END DO
			END IF
		END DO
		END IF
	END DO
	DEALLOCATE(GammaTemp)
END DO
temp=CMPLX(0.0,KIND=rKind);

DO alpha=1,SIZE(GammasL(locsysSize)%m,1)
	IF(LabelL(locsysSize)%vi(alpha).lt.10000) THEN
	DO beta=1,SIZE(GammasR(locsysSize)%m,1)
		IF(LabelR(locsysSize)%vi(beta).lt.10000) THEN
			IF(LabelL(locsysSize)%vi(alpha)==LabelR(locsysSize)%vi(beta)) THEN
		temp = temp+CONJG(GammasL(locSyssize)%m(alpha,1))*normKernal(alpha,beta) &
			*GammasR(locSyssize)%m(beta,1)
			END IF
		END IF
	END DO
	END IF
END DO

DEALLOCATE(normKernal)
InnerProduct_m=temp

END FUNCTION InnerProduct_m

SUBROUTINE SetupMeasures(Measures,HamiType)
!
!Purpose:Setup the measure derived type for the TEBD ALPS interface
!
TYPE(measure) :: Measures
CHARACTER(len=*) :: HamiType
INTEGER :: i


SELECT CASE(HamiType)
	CASE ('spin')
		!Sz, Sx,Sz^2, Sx^2, SzSz, SxSx
		CALL AllocateMeasures(Measures,.TRUE.,4, 4, 2, 0, 1, 1)
		Measures%local(1)%op=Sz%m
		Measures%local(2)%op=MATMUL(Sz%m,Sz%m)
		Measures%local(3)%op=Sx%m
		Measures%local(4)%op=MATMUL(Sx%m,Sx%m)
		Measures%avg(1)%op=Sz%m
		Measures%avg(2)%op=MATMUL(Sz%m,Sz%m)
		Measures%avg(3)%op=Sx%m
		Measures%avg(4)%op=MATMUL(Sx%m,Sx%m)
		Measures%Corr(1)%op=TensorProd(Sz%m,Sz%m)
		Measures%Corr(2)%op=TensorProd(Sx%m,Sx%m)
		DO i=1,2
			CALL SplitOperator(Measures%Corr(i)%op,Measures%Corr(i)%mpo)
		END DO
	CASE ('boson Hubbard')
		!n, n^2, nn, a^{\dagger} a
		CALL AllocateMeasures(Measures,.TRUE.,2, 2, 2, 0, 1, 1)
		Measures%local(1)%op=n_op%m
		Measures%local(2)%op=MATMUL(n_op%m,n_op%m)
		Measures%avg(1)%op=n_op%m
		Measures%avg(2)%op=MATMUL(n_op%m,n_op%m)
		Measures%Corr(1)%op=TensorProd(n_op%m,n_op%m)
		Measures%Corr(2)%op=TensorProd(TRANSPOSE(a_op%m),a_op%m)
		DO i=1,2
			CALL SplitOperator(Measures%Corr(i)%op,Measures%Corr(i)%mpo)
		END DO
	CASE ('hardcore boson')
		!n, n^2, nn, a^{\dagger} a
		CALL AllocateMeasures(Measures,.TRUE.,2, 2, 2, 0, 1, 1)
		Measures%local(1)%op=n_op%m
		Measures%local(2)%op=MATMUL(n_op%m,n_op%m)
		Measures%avg(1)%op=n_op%m
		Measures%avg(2)%op=MATMUL(n_op%m,n_op%m)
		Measures%Corr(1)%op=TensorProd(n_op%m,n_op%m)
		Measures%Corr(2)%op=TensorProd(TRANSPOSE(a_op%m),a_op%m)
		DO i=1,2
			CALL SplitOperator(Measures%Corr(i)%op,Measures%Corr(i)%mpo)
		END DO
	CASE ('fermion Hubbard')
		!n, n^2,Sz, Sz^2, nn,SzSz, a^{\dagger} a
		CALL AllocateMeasures(Measures,.TRUE.,4, 4, 2, 2, 1, 1)
		Measures%local(1)%op=n_op%m
		Measures%local(2)%op=MATMUL(n_op%m,n_op%m)
		Measures%local(3)%op=Sz%m
		Measures%local(4)%op=MATMUL(Sz%m,Sz%m)
		Measures%avg(1)%op=n_op%m
		Measures%avg(2)%op=MATMUL(n_op%m,n_op%m)
		Measures%avg(3)%op=Sz%m
		Measures%avg(4)%op=MATMUL(Sz%m,Sz%m)
		Measures%Corr(1)%op=TensorProd(n_op%m,n_op%m)
		Measures%Corr(2)%op=TensorProd(Sz%m,Sz%m)
		Measures%fermicorr(1)%op=TensorProd(TRANSPOSE(a_opS(1)%m),a_opS(1)%m)
		Measures%fermicorr(2)%op=TensorProd(TRANSPOSE(a_opS(2)%m),a_opS(2)%m)
		DO i=1,2
			CALL SplitOperator(Measures%Corr(i)%op,Measures%Corr(i)%mpo)
		END DO
		DO i=1,2
			CALL SplitOperator(Measures%fermiCorr(i)%op,Measures%fermiCorr(i)%mpo)
		END DO
	CASE ('spinless fermions')
		!n, n^2, nn, a^{\dagger} a
		CALL AllocateMeasures(Measures,.TRUE.,2, 2, 1, 1, 1, 1)
		Measures%local(1)%op=n_op%m
		Measures%local(2)%op=MATMUL(n_op%m,n_op%m)
		Measures%avg(1)%op=n_op%m
		Measures%avg(2)%op=MATMUL(n_op%m,n_op%m)
		Measures%Corr(1)%op=TensorProd(n_op%m,n_op%m)
		Measures%fermiCorr(1)%op=TensorProd(TRANSPOSE(a_op%m),a_op%m)
		DO i=1,1
			CALL SplitOperator(Measures%Corr(i)%op,Measures%Corr(i)%mpo)
		END DO
		DO i=1,1
			CALL SplitOperator(Measures%fermiCorr(i)%op,Measures%fermiCorr(i)%mpo)
		END DO
	CASE DEFAULT
		PRINT *, "Hamiltonian type not recognized!"
		PRINT *, "Use 'spin', 'boson Hubbard', 'hardcore boson', 'fermion Hubbard', or 'spinless fermions'."
		STOP 
END SELECT


END SUBROUTINE SetupMeasures

SUBROUTINE AllocateMeasures(Measures,en,numLocal, numAvg, numCorr, numFermiCorr, numEnt, numOverlaps)
!
!Purpose: Allocate the measure derived type to hold num* measures of * kind
!
IMPLICIT NONE
TYPE(measure) :: Measures
LOGICAL, INTENT(IN) :: en
INTEGER, INTENT(IN) :: numLocal, numAvg, numCorr, numFermiCorr, numEnt, numOverlaps
INTEGER :: i

Measures%enPres=en

IF(numlocal.gt.0) THEN
	ALLOCATE(Measures%local(numLocal))
	DO i=1,numLocal
		ALLOCATE(Measures%local(i)%Op(localSize,localSize))
		ALLOCATE(Measures%local(i)%value(systemSize))
	END DO
	Measures%localpres=.TRUE.
ELSE
	Measures%localpres=.FALSE.
END IF

IF(numavg.gt.0) THEN
	ALLOCATE(Measures%avg(numAvg))
	DO i=1,numAvg
		ALLOCATE(Measures%avg(i)%Op(localSize,localSize))
	END DO
	Measures%avgpres=.TRUE.
ELSE
	Measures%avgpres=.FALSE.
END IF

SELECT CASE (numEnt)
CASE (1)
	Measures%entpres=1
	ALLOCATE(Measures%ent%vN(systemSize))
	ALLOCATE(Measures%ent%chain(systemSize+1))
CASE DEFAULT
	Measures%entpres=0
END SELECT

IF(numCorr.gt.0) THEN
	ALLOCATE(Measures%corr(numCorr))
	DO i=1,numCorr
		ALLOCATE(Measures%corr(i)%Op(localSize*localSize,localSize*localSize))	
		ALLOCATE(Measures%corr(i)%value(systemSize,systemSize))
	END DO
	Measures%corrpres=.TRUE.
ELSE
	Measures%corrpres=.FALSE.
END IF

IF(numfermicorr.gt.0) THEN
	ALLOCATE(Measures%fermicorr(numFermiCorr))
	DO i=1,numFermiCorr
		ALLOCATE(Measures%fermicorr(i)%Op(localSize*localSize,localSize*localSize))	
		ALLOCATE(Measures%fermicorr(i)%value(systemSize,systemSize))
	END DO
	Measures%fermicorrpres=.TRUE.
ELSE
	Measures%fermicorrpres=.FALSE.
END IF

IF(numOverlaps.gt.0) THEN
	ALLOCATE(Measures%overlap(numOverlaps))
	Measures%overlapPres=.TRUE.
ELSE
	Measures%overlapPres=.FALSE.
END IF

END SUBROUTINE AllocateMeasures


SUBROUTINE AllocateOverlap_t(Measures, overlapnum, Gammas, Lambdas)
!
!Purpose: Allocate the overlapnum^th overlap MPS
!
IMPLICIT NONE
TYPE(measure) :: Measures
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER :: overlapNum
INTEGER :: chiMax, i

chiMax=1
DO i=2,systemSize+1
	chiMax=MAX(chiMax,SIZE(Lambdas(i)%v,1))
END DO

Measures%overlap(overlapNum)%m=.FALSE.
CALL AllocateGamLam(Measures%overlap(overlapNum)%Gammas,Measures%overlap(overlapNum)%Lambdas,chiMax)
CALL CopyGamLam(Measures%overlap(overlapNum)%Gammas,Measures%overlap(overlapNum)%Lambdas,Gammas,Lambdas)

END SUBROUTINE AllocateOverlap_t

SUBROUTINE AllocateOverlap_tl(Measures, overlapnum, Gammas, Lambdas,LabelLeft, LabelRight)
!
!Purpose: Allocate the overlapnum^th overlap MPS
!
IMPLICIT NONE
TYPE(measure) :: Measures
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: overlapNum
INTEGER :: chiMax, i

chiMax=1
DO i=2,systemSize+1
	chiMax=MAX(chiMax,SIZE(Lambdas(i)%v,1))
END DO

Measures%overlap(overlapNum)%m=.FALSE.
CALL AllocateGamLam(Measures%overlap(overlapNum)%Gammas,Measures%overlap(overlapNum)%Lambdas,chiMax)
CALL CopyGamLam(Measures%overlap(overlapNum)%Gammas,Measures%overlap(overlapNum)%Lambdas,Gammas,Lambdas)
CALL AllocateLabel(Measures%overlap(overlapNum)%LL, Measures%overlap(overlapNum)%LR, chiMax)
CALL CopyLabel(Measures%overlap(overlapNum)%LL, Measures%overlap(overlapNum)%LR, LabelLeft, LabelRight)

END SUBROUTINE AllocateOverlap_tl

SUBROUTINE AllocateOverlap_m(Measures, overlapnum, Gammas, Lambdas,LabelLeft, LabelRight)
!
!Purpose: Allocate the overlapnum^th overlap MPS
!
IMPLICIT NONE
TYPE(measure) :: Measures
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: overlapNum
INTEGER :: chiMax, i

chiMax=1
DO i=2,systemSize+1
	chiMax=MAX(chiMax,SIZE(Lambdas(i)%v,1))
END DO

Measures%overlap(overlapNum)%m=.TRUE.
CALL AllocateGamLam(Measures%overlap(overlapNum)%Gammasm,Measures%overlap(overlapNum)%Lambdas,chiMax)
CALL CopyGamLam(Measures%overlap(overlapNum)%Gammasm,Measures%overlap(overlapNum)%Lambdas,Gammas,Lambdas)
CALL AllocateLabel(Measures%overlap(overlapNum)%LL, Measures%overlap(overlapNum)%LR, chiMax)
CALL CopyLabel(Measures%overlap(overlapNum)%LL, Measures%overlap(overlapNum)%LR, LabelLeft, LabelRight)

END SUBROUTINE AllocateOverlap_m

SUBROUTINE DEAllocateOverlap(Measures)
!
!Purpose: Allocate the overlapnum^th overlap MPS
!
IMPLICIT NONE
TYPE(measure) :: Measures
INTEGER :: chiMax, i
DO i=1, SIZE(Measures%overlap)
	IF(Measures%overlap(i)%m) THEN
		CALL DeallocateGamLam(Measures%overlap(i)%Gammasm,Measures%overlap(i)%Lambdas)
		CALL DEallocateLabel(Measures%overlap(i)%LL, Measures%overlap(i)%LR)
	ELSE
		CALL DeAllocateGamLam(Measures%overlap(i)%Gammas,Measures%overlap(i)%Lambdas)
	END IF
END DO

END SUBROUTINE DeAllocateOverlap

SUBROUTINE DeallocateMeasures(Measures)
!
!Purpose: Deallocate the measure derived type
!
IMPLICIT NONE
TYPE(measure) :: Measures
INTEGER :: i

IF(Measures%localpres) THEN
	DO i=1,SIZE(Measures%local)
		DEALLOCATE(Measures%local(i)%Op)
		DEALLOCATE(Measures%local(i)%value)
	END DO
	DEALLOCATE(Measures%local)
END IF


IF(Measures%avgpres) THEN
	DO i=1,SIZE(Measures%avg)
		DEALLOCATE(Measures%avg(i)%Op)
	END DO
	DEALLOCATE(Measures%avg)
END IF

SELECT CASE (Measures%entpres)
CASE (1)
	DEALLOCATE(Measures%ent%vN)
	DEALLOCATE(Measures%ent%chain)
CASE DEFAULT
END SELECT

IF(Measures%corrpres) THEN
	DO i=1,SIZE(Measures%corr)
		DEALLOCATE(Measures%corr(i)%value)	
		DEALLOCATE(Measures%corr(i)%Op)
		CALL DestroySplitOperator(Measures%corr(i)%mpo)
	END DO
	DEALLOCATE(Measures%corr)
END IF

IF(Measures%fermicorrpres) THEN
	DO i=1,SIZE(Measures%fermicorr)
		DEALLOCATE(Measures%fermicorr(i)%value)	
		DEALLOCATE(Measures%fermicorr(i)%Op)
		CALL DestroySplitOperator(Measures%fermicorr(i)%mpo)
	END DO
	DEALLOCATE(Measures%fermicorr)
END IF

IF(Measures%overlapPres) THEN
	CALL DeallocateOverlap(Measures)
END IF

END SUBROUTINE DeallocateMeasures

SUBROUTINE EvaluateMeasures_t(Measures, Gammas, Lambdas)
!
!Purpose: Evaluate the measures stored in the measure derived type
!
IMPLICIT NONE
TYPE(measure) :: Measures
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
TYPE(fourtensor) :: Theta, Thetafermi
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
REAL(KIND=rKIND) :: temp
REAL(KIND=rKIND), ALLOCATABLE :: evals(:)
INTEGER :: i,j, l1,l2, numO
INTEGER :: workSize, info, numOps
COMPLEX(KIND=rKind), ALLOCATABLE :: workArr(:), rworkArr(:)
TYPE(DecomposedMPO) :: MPO

!!!!!!!!!!!!Local Observables!!!!!!!!!!!!!!
ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate rho in EvaluateMeasures'
END IF

!Zero out average measures
IF(Measures%avgpres) THEN
	DO i=1,SIZE(Measures%avg)
		Measures%avg(i)%value=0.0_rKind
	END DO
END IF
!Initialize single-site entropies
IF(Measures%entpres.ge.1) THEN
	Measures%ent%qme=localSize*1.0_rKind/(localSize*1.0_rKind-1.0_rKind)
!Allocate workspace for ZHEEV
	workSize = 2*localSize-1
	ALLOCATE(workarr(workSize), rworkarr(3*localSize-2),evals(localSize), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate workarr in LocalEntropyDist'
	END IF
END IF


!Construct the single-site density matrix at every site
DO i=1,systemSize,1
	rho%m=0.0_rKind
	IF(i==1) THEN
		CALL FormSingleSiteRho(rho%m,  Gammas(i)%t, Lambdas(i+1)%v)
	ELSE
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
	END IF

!local measures
IF(Measures%localpres) THEN
	DO numO=1,SIZE(Measures%local)
		Measures%local(numO)%value(i)=0.0_rKind
		Measures%local(numO)%value(i)=TraceMatmul(Measures%local(numO)%Op,rho%m)
	END DO
END IF

!average measures
IF(Measures%avgpres) THEN
	DO numO=1,SIZE(Measures%avg)
		Measures%avg(numO)%value=Measures%avg(numO)%value+TraceMatmul(Measures%avg(numO)%Op,rho%m)/(systemSize*1.0_rKind)
	END DO

END IF


!One-body entropies
IF(Measures%entpres.ge.1) THEN
	Measures%ent%qme=Measures%ent%qme-TraceMatmul(rho%m, rho%m)*localSize/(1.0_rKind*systemSize*(localSize*1.0_rKind-1.0_rKind))
	CALL ZHEEV('N', 'U', localSize, rho%m, localSize, evals, workArr, workSize, rworkArr, info)
	Measures%ent%vN(i)=0.0_rKind
	DO j = 1, localSize
		IF(evals(j).ne.0.0_rKind) THEN
        	Measures%ent%vN(i) = Measures%ent%vN(i)-evals(j)*LOG(ABS(evals(j)))/LOG(1.0_rKind*localSize) ! Compute tr(rho*log_d(rho)).
		ELSE
			Measures%ent%vN(i)=Measures%ent%vN(i)+0.0_rKind
		END IF
	END DO
END IF

END DO				

!Chain entropy
IF(Measures%entpres.ge.1) THEN
	DO i=1,systemSize+1
		Measures%ent%chain(i)=ChainEntropy(i,Lambdas)
	END DO
END IF
		
DEALLOCATE(rho%m, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate rho in EvaluateMeasures'
END IF

IF(Measures%entpres.ge.1) THEN
	DEALLOCATE(workarr, rworkarr,evals, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate variables in LocalEntropyDist'
			END IF
END IF

IF(Measures%enPres) THEN
	!Compute the energy
	CALL TotalEnergy(Measures%en, Measures%HMPO, Gammas, Lambdas)
END IF


IF(measures%corrpres) THEN
	ALLOCATE(rho%m(systemSize,systemSize))
	DO i=1,SIZE(Measures%corr)
		Measures%corr(i)%value=0.0_rKind
		DO j=1,SIZE(Measures%corr(i)%MPO%l)
			CALL TwoSiteExpValG(rho%m, Measures%corr(i)%MPO%l(j)%m, Measures%corr(i)%MPO%r(j)%m, Gammas, Lambdas)
			Measures%corr(i)%value=Measures%corr(i)%value+rho%m		
		END DO
	END DO
	DEALLOCATE(rho%m)
END IF
IF(measures%fermicorrpres) THEN
	ALLOCATE(rho%m(systemSize,systemSize))
	DO i=1,SIZE(Measures%fermicorr)
		Measures%fermicorr(i)%value=0.0_rKind
		DO j=1,SIZE(Measures%fermicorr(i)%MPO%l)
			CALL TwoSiteExpValG(rho%m, Measures%fermicorr(i)%MPO%l(j)%m, Measures%fermicorr(i)%MPO%r(j)%m, Gammas, Lambdas, phaseStat=1)
			Measures%fermicorr(i)%value=Measures%fermicorr(i)%value+rho%m
		END DO
	END DO
	DEALLOCATE(rho%m)
END IF

!Measure overlap
IF(Measures%overlapPres) THEN
	DO i=1,SIZE(Measures%overlap)
		Measures%overlap(i)%value=InnerProduct(Gammas,Lambdas, Measures%overlap(i)%Gammas,Measures%overlap(i)%Lambdas)
	END DO
END IF

END SUBROUTINE EvaluateMeasures_t

SUBROUTINE EvaluateMeasures_m(Measures, Gammas, Lambdas,LabelLeft)
!
!Purpose: Evaluate the measures stored in the measure derived type
!
IMPLICIT NONE
TYPE(measure) :: Measures
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
TYPE(matrix) :: rho
TYPE(fourtensor) :: Theta, Thetafermi
TYPE(matrix) :: GammaP
REAL(KIND=rKIND) :: temp
REAL(KIND=rKIND), ALLOCATABLE :: evals(:)
INTEGER :: i,j, l1,l2, numO
INTEGER :: workSize, info, numOps
COMPLEX(KIND=rKind), ALLOCATABLE :: workArr(:), rworkArr(:)
TYPE(DecomposedMPO) :: MPO

!!!!!!!!!!!!Local Observables!!!!!!!!!!!!!!
ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate rho in EvaluateMeasures'
END IF

!Zero out average measures
IF(Measures%avgpres) THEN
	DO i=1,SIZE(Measures%avg)
		Measures%avg(i)%value=0.0_rKind
	END DO
END IF
!Initialize single-site entropies
IF(Measures%entpres.ge.1) THEN
	Measures%ent%qme=localSize*1.0_rKind/(localSize*1.0_rKind-1.0_rKind)
!Allocate workspace for ZHEEV
	workSize = 2*localSize-1
	ALLOCATE(workarr(workSize), rworkarr(3*localSize-2),evals(localSize), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate workarr in LocalEntropyDist'
	END IF
END IF


!Construct the single-site density matrix at every site
DO i=1,systemSize,1
	rho%m=0.0_rKind
	IF(i==1) THEN
		CALL FormSingleSiteRho(rho%m,  Gammas(i)%m, Lambdas(i+1)%v, LabelLeft(i)%vi,LabelLeft(i+1)%vi)
	ELSE
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%m, Lambdas(i+1)%v, LabelLeft(i)%vi,LabelLeft(i+1)%vi)
	END IF

!local measures
IF(Measures%localpres) THEN
	DO numO=1,SIZE(Measures%local)
		Measures%local(numO)%value(i)=0.0_rKind
		Measures%local(numO)%value(i)=TraceMatmul(Measures%local(numO)%Op,rho%m)
	END DO
END IF

!average measures
IF(Measures%avgpres) THEN
	DO numO=1,SIZE(Measures%avg)
		Measures%avg(numO)%value=Measures%avg(numO)%value+TraceMatmul(Measures%avg(numO)%Op,rho%m)/(systemSize*1.0_rKind)
	END DO

END IF


!One-body entropies
IF(Measures%entpres.ge.1) THEN
	Measures%ent%qme=Measures%ent%qme-TraceMatmul(rho%m, rho%m)*localSize/(1.0_rKind*systemSize*(localSize*1.0_rKind-1.0_rKind))
	CALL ZHEEV('N', 'U', localSize, rho%m, localSize, evals, workArr, workSize, rworkArr, info)
	Measures%ent%vN(i)=0.0_rKind
	DO j = 1, localSize
		IF(evals(j).ne.0.0_rKind) THEN
        	Measures%ent%vN(i) = Measures%ent%vN(i)-evals(j)*LOG(ABS(evals(j)))/LOG(1.0_rKind*localSize) ! Compute tr(rho*log_d(rho)).
		ELSE
			Measures%ent%vN(i)=Measures%ent%vN(i)+0.0_rKind
		END IF
	END DO
END IF


END DO				

!Chain entropy
IF(Measures%entpres.ge.1) THEN
	DO i=1,systemSize+1
		Measures%ent%chain(i)=ChainEntropy(i,Lambdas)
	END DO
END IF
		
DEALLOCATE(rho%m, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate rho in EvaluateMeasures'
END IF

IF(Measures%entpres.ge.1) THEN
	DEALLOCATE(workarr, rworkarr,evals, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate variables in LocalEntropyDist'
			END IF
END IF

IF(Measures%enPres) THEN
	!Compute the energy
	CALL TotalEnergy(Measures%en, Measures%HMPO, Gammas, Lambdas, LabelLeft)
END IF



IF(measures%corrpres) THEN
	ALLOCATE(rho%m(systemSize,systemSize))
	DO i=1,SIZE(Measures%corr)
		Measures%corr(i)%value=0.0_rKind
		DO j=1,SIZE(Measures%corr(i)%MPO%l)
			CALL TwoSiteExpValG(rho%m, Measures%corr(i)%MPO%l(j)%m, Measures%corr(i)%MPO%r(j)%m, Gammas, Lambdas, LabelLeft)
			Measures%corr(i)%value=Measures%corr(i)%value+rho%m
		END DO
	END DO
	DEALLOCATE(rho%m)
	DO numO=1,SIZE(Measures%local)
	END DO
END IF
IF(measures%fermicorrpres) THEN
	ALLOCATE(rho%m(systemSize,systemSize))
	DO i=1,SIZE(Measures%fermicorr)
		Measures%fermicorr(i)%value=0.0_rKind
		DO j=1,SIZE(Measures%fermicorr(i)%MPO%l)
			CALL TwoSiteExpValG(rho%m, Measures%fermicorr(i)%MPO%l(j)%m, &
		Measures%fermicorr(i)%MPO%r(j)%m, Gammas, Lambdas, LabelLeft, phaseStat=1)
			Measures%fermicorr(i)%value=Measures%fermicorr(i)%value+rho%m
		END DO
	END DO
	DEALLOCATE(rho%m)
END IF

IF(Measures%overlapPres) THEN
	DO i=1,SIZE(Measures%overlap)
		Measures%overlap(i)%value=InnerProduct(Gammas,Lambdas,LabelLeft, &
		Measures%overlap(i)%Gammasm,Measures%overlap(i)%Lambdas,Measures%overlap(i)%LL)
	END DO
END IF


END SUBROUTINE EvaluateMeasures_m


END MODULE ObOps
