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
MODULE LocalOps
!
! Purpose: Module to perform local operations
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
USE omp_lib

IMPLICIT NONE


INTERFACE OneSiteOp
	MODULE PROCEDURE OneSiteOp_rt,OneSiteOp_ct,OneSiteOp_rtm,OneSiteOp_ctm,OneSiteOp_rm,OneSiteOp_cm
END INTERFACE OneSiteOp


INTERFACE TwoSiteOpNCR
	MODULE PROCEDURE TwoSiteOpNCR_t, TwoSiteOpNCR_m
END INTERFACE TwoSiteOpNCR

INTERFACE TwoSiteOpNCL
	MODULE PROCEDURE TwoSiteOpNCL_t, TwoSiteOpNCL_m
END INTERFACE TwoSiteOpNCL

INTERFACE TwoSiteOpNCEdge
	MODULE PROCEDURE TwoSiteOpNCEdge_t, TwoSiteOpNCEdge_m
END INTERFACE TwoSiteOpNCEdge

INTERFACE FormThetaNCL
	MODULE PROCEDURE FormThetaNCL_t, FormThetaNCL_m
END INTERFACE FormThetaNCL

INTERFACE FormThetaNCR
	MODULE PROCEDURE FormThetaNCR_t, FormThetaNCR_m
END INTERFACE FormThetaNCR

INTERFACE FormThetaNCEdge
	MODULE PROCEDURE FormThetaNCEdge_t, FormThetaNCEdge_m
END INTERFACE FormThetaNCEdge

INTERFACE FormGamma1NC
	MODULE PROCEDURE FormGamma1NC_t,FormGamma1NC_m
END INTERFACE FormGamma1NC

INTERFACE FormGamma2NC
	MODULE PROCEDURE FormGamma2NC_t,FormGamma2NC_m
END INTERFACE FormGamma2NC




! These variables are set for using the subroutine "ZGESVD" in LAPACK, which performs an SVD on a general matrix.
CHARACTER(1) :: jobu_SVD, jobvt_SVD
INTEGER :: matrixSizeSM_SVD, workSizeSM_SVD, matrixSizeLG_SVD, workSizeLG_SVD, &
matrixSize_SVD, workSize_SVD, info_SVD, matrixSizeL_SVD, matrixSizeT_SVD

REAL(KIND=rKind), ALLOCATABLE :: rworkSM_SVD(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: workSM_SVD(:)
REAL(KIND=rKind), ALLOCATABLE :: rworkLG_SVD(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: workLG_SVD(:)
REAL(KIND=rKind), ALLOCATABLE :: rwork_SVD(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: work_SVD(:)	

LOGICAL :: tRenorm=.FALSE., lrenorm=.true.
	
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Number Non-conserving method starts !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OneSiteOp_rt(Op1,Gamma)
!
!Purpose: Perform the real one-site operation Op1 on the Gamma specified
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: Op1(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma(:,:,:)
INTEGER alpha
	
DO alpha = 1,SIZE(Gamma,1)
	Gamma(alpha,:,:) = MATMUL(Op1,Gamma(alpha,:,:));
END DO		
END SUBROUTINE OneSiteOp_rt

SUBROUTINE OneSiteOp_ct(Op1,Gamma)
!
!Purpose: Perform the complex one-site operation Op1 on the Gamma specified
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma(:,:,:)
INTEGER alpha

DO alpha = 1,SIZE(Gamma,1)
	Gamma(alpha,:,:) = MATMUL(Op1,Gamma(alpha,:,:));
END DO		
END SUBROUTINE OneSiteOp_ct

SUBROUTINE FormThetaL(Theta,  Gamma1, Lambda, Gamma2)
!
!Purpose: Form Theta as defined in the Manual
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma1(:,:,:), Gamma2(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda(:)
COMPLEX(KIND=rKind), INTENT(OUT) :: Theta(:,:,:,:)
INTEGER :: chi0,chi1,chi2,alpha,beta,gamma,i,j,k,l
chi0 = SIZE(Gamma1,1)
chi1 = SIZE(Gamma1,3)
chi2 = SIZE(Gamma2,3)
Theta=0.0_rKind

!$OMP PARALLEL DO PRIVATE(gamma,beta,i,j) schedule(dynamic)		
DO alpha=1,chi0
		DO i=1,SIZE(Gamma1,2)
			DO gamma=1,chi1
		Theta(alpha,i,:,:)=Theta(alpha,i,:,:) & 
		+ Lambda(alpha)*Gamma1(alpha,i,gamma)*Gamma2(gamma,:,:)
			END DO 
		END DO 
END DO 
!$OMP END PARALLEL DO

END SUBROUTINE FormThetaL

SUBROUTINE FormThetaR(Theta,   Gamma1, Lambda, Gamma2)
!
!Purpose: Form Theta as defined in the Manual
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma1(:,:,:), Gamma2(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda(:)
COMPLEX(KIND=rKind), INTENT(OUT) :: Theta(:,:,:,:)
INTEGER :: chi0,chi1,chi2,alpha,beta,gamma,i,j,k,l
chi0 = SIZE(Gamma1,1)
chi1 = SIZE(Gamma1,3)
chi2 = SIZE(Gamma2,3)
Theta=0.0_rKind

!$OMP PARALLEL DO PRIVATE(gamma,beta,i,j) schedule(dynamic)		
DO alpha=1,chi2
		DO i=1,SIZE(Gamma2,2)
			DO gamma=1,chi1
		Theta(:,:,i,alpha)=Theta(:,:,i,alpha) & 
		+ Gamma1(:,:,gamma)*Gamma2(gamma,i,alpha)*Lambda(alpha)
			END DO 
		END DO 
END DO 
!$OMP END PARALLEL DO

END SUBROUTINE FormThetaR

SUBROUTINE FormThetaEdge(Theta,  Gamma1, Lambda, Gamma2)
!
!Purpose: Form Theta as defined in the Manual
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma1(:,:,:), Gamma2(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda(:)
COMPLEX(KIND=rKind), INTENT(OUT) :: Theta(:,:,:,:)
INTEGER :: chi0,chi1,chi2,alpha,beta,gamma,i,j,k,l
chi0 = SIZE(Gamma1,1)
chi1 = SIZE(Gamma1,3)
chi2 = SIZE(Gamma2,3)
Theta=0.0_rKind

!$OMP PARALLEL DO PRIVATE(gamma,beta,i,j) schedule(dynamic)		
DO alpha=1,chi2
		DO i=1,SIZE(Gamma2,2)
			DO gamma=1,chi1
		Theta(:,:,i,alpha)=Theta(:,:,i,alpha) & 
		+ Gamma1(:,:,gamma)*Lambda(gamma)*Gamma2(gamma,i,alpha)
			END DO 
		END DO 
END DO 
!$OMP END PARALLEL DO

END SUBROUTINE FormThetaEdge

SUBROUTINE ReshapeTheta(Theta,ThetaRS)
!
!Purpose: Reshape the chiXdXdXchi 4-tensor Theta into a (chi d)X(chi d) matrix ThetaRS and (optionally) renormalize such that
!\sum_{aijb}|Theta_{ab}^{ij}|^2=1
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind), INTENT(OUT) :: ThetaRS(SIZE(Theta,1)*SIZE(Theta,2),SIZE(Theta,3)*SIZE(Theta,4))
COMPLEX(KIND=rKind) :: normSQcmplx
REAL(KIND=rKind) :: norm
INTEGER :: chi0,chi1,alpha,beta,gamma,i,j,k,l
		
chi0 = SIZE(Theta,1)
chi1 = SIZE(Theta,4)

!$OMP PARALLEL DO schedule(dynamic)
DO alpha=1,chi0
	DO beta=1,chi1
		DO i=1,SIZE(Theta,2)
			DO j=1,SIZE(Theta,3)
				ThetaRS((i-1)*chi0+alpha,(j-1)*chi1+beta)=Theta(alpha,i,j,beta)
			END DO 
		END DO 
	END DO 
END DO
!$OMP END PARALLEL DO

END SUBROUTINE ReshapeTheta

SUBROUTINE ThetaOperation(Op2,Theta)
!
!Purpose: Operate the complex two-site operation Op2 on Theta
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN)  :: Op2(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind) :: ThetaTemp(SIZE(Theta,1),SIZE(Theta,2),SIZE(Theta,3),SIZE(Theta,4))
INTEGER :: chi0,chi2,alpha,beta,gamma,i,j,k,l
		
chi0 = SIZE(Theta,1)
chi2 = SIZE(Theta,4)
ThetaTemp=Theta

Theta=0.0_rKind
!$OMP PARALLEL DO PRIVATE(i,j,k,l) schedule(dynamic)
DO alpha=1,chi0				
DO i=1,SIZE(Theta,2)
	DO j=1,SIZE(Theta,3)
		DO k=1,SIZE(Theta,2)
			DO l=1,SIZE(Theta,3)
				Theta(alpha,i,j,:)=Theta(alpha,i,j,:) & 
						+ Op2((i-1)*SIZE(Theta,3)+j,(k-1)*SIZE(Theta,3)+l) &
						* ThetaTemp(alpha,k,l,:)
			END DO 
		END DO 
	END DO 
END DO 
END DO
!$OMP END PARALLEL DO
	
END SUBROUTINE ThetaOperation

SUBROUTINE SVDTruncation(link, MatrixIn, S, U, V)
!
!Purpose: Perform an SVD on the reshaped Theta, keep only the largest chi singular values
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind), INTENT(INOUT) :: MatrixIn(:,:)       
TYPE(vector) :: S
TYPE(matrix) :: U, V
INTEGER :: i, j, dchi1, dchi2, M, N, minn

M=SIZE(MatrixIn,1)
N=SIZE(MatrixIn,2)

minn=MINVAL((/M,N/))

ALLOCATE(U%m(M,minn), V%m(minn,N), S%v(minn))

ALLOCATE(workLG_SVD(5*MAXVAL((/M,N/))))
ALLOCATE(rworkLG_SVD(5*minn))

!Call the LAPACK routine ZGESVD, which performs a SVD on a general matrix
	CALL ZGESVD('S', 'S', M, N, MatrixIn, M, S%v, U%m, & 
				M, V%m, minn, workLG_SVD, 5*MAXVAL((/M,N/)), rworkLG_SVD, info_SVD)

DEALLOCATE(workLG_SVD,rworkLG_SVD)

END SUBROUTINE SVDTruncation


SUBROUTINE FormLambda(Lambda, truncerr, S)
!
!Purpose: Form Lambda from the singular values
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: truncerr
TYPE(vector) :: Lambda
REAL(KIND=rKind), INTENT(IN) :: S(:)
REAL(KIND=rKind) :: norm, norm2
INTEGER :: alpha, counter

!Renormalize the kept singular values such that the sum of the squares is 1
norm=0.0_rKind
norm2=DOT_PRODUCT(S,S)
counter=0
DO alpha=1, SIZE(S)
	IF((ABS(1.0_rKind-norm/norm2).ge.truncLimit).and.(counter.lt.chiLimit)) THEN
		norm=norm + S(alpha) ** 2
		counter=counter+1
	ELSE
		norm=norm
	END IF
END DO

!Schmidt error		
truncerr =ABS( 1.0_rKind - norm/norm2)
!Redefine Lambda	
DEALLOCATE(Lambda%v)
ALLOCATE(Lambda%v(counter))	
DO alpha=1, counter
	IF(ABS(S(alpha))>10.0_rKind**(-15)) THEN
		Lambda%v(alpha)=S(alpha)/SQRT(norm)
	ELSE
		Lambda%v(alpha)=0.0_rKind
	END IF
END DO		
END SUBROUTINE FormLambda

SUBROUTINE FormGamma1(Gamma1,U,chi0,chi1)
!
!Purpose: Form Gamma to the left of the link from the SVD.
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma1(:,:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: U(:,:)
INTEGER, INTENT(IN) :: chi0, chi1
INTEGER :: alpha, beta, i

Gamma1=0.0_rKind
!$OMP PARALLEL DO PRIVATE(i) schedule(dynamic)
DO alpha=1,chi0
	DO i=1,SIZE(Gamma1,2)
		Gamma1(alpha,i,:)=U((i-1)*chi0+alpha,1:chi1)
	END DO 
END DO		
!$OMP END PARALLEL DO			

END SUBROUTINE FormGamma1
	
SUBROUTINE FormGamma2(Gamma2,V,chi1,chi2)
!
!Purpose: Form Gamma to the right of the link from the SVD.
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma2(:,:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: V(:,:)
INTEGER, INTENT(IN) :: chi1, chi2
INTEGER :: alpha, beta, i
		
Gamma2=0.0_rKind
!$OMP PARALLEL DO PRIVATE(i) schedule(dynamic)
DO beta=1,chi2
	DO i=1,SIZE(Gamma2,2)
		Gamma2(:,i,beta)=V(1:chi1,(i-1)*chi2+beta)
	END DO 
END DO	
!$OMP END PARALLEL DO
END SUBROUTINE FormGamma2

SUBROUTINE TwoSiteOpL(link,Op2,Gammas,Lambdas,truncerr)
!
!Purpose: Perform the complex two-site operation Op2 on the sites neighboring link
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
COMPLEX(KIND=rKind) :: ThetaRS(SIZE(Gammas(link)%t,1)*SIZE(Gammas(link)%t,2), &
SIZE(Gammas(link+1)%t,2)*SIZE(Gammas(link+1)%t,3))
TYPE(matrix) :: U, V
TYPE(vector) :: S
INTEGER :: chi0, chi1, chi2, ld

chi0=Size(Gammas(link)%t,1)
chi1=Size(Gammas(link)%t,3)
chi2=Size(Gammas(link+1)%t,3)

CALL FormThetaL(Theta,Gammas(link)%t,Lambdas(link)%v,Gammas(link+1)%t)
CALL ThetaOperation(Op2,Theta)
CALL ReshapeTheta(Theta,ThetaRS)
CALL SVDTruncation(link,ThetaRS, S, U, V)

CALL FormLambda(Lambdas(link+1), truncerr, S%v)
ld=SIZE(Gammas(link)%t,2)
chi1=SIZE(Lambdas(link+1)%v)
DEALLOCATE(Gammas(link)%t)
ALLOCATE(Gammas(link)%t(chi0,ld,chi1))

CALL FormGamma1(Gammas(link)%t,U%m,chi0,chi1)
ld=SIZE(Gammas(link+1)%t,2)
DEALLOCATE(Gammas(link+1)%t)
ALLOCATE(Gammas(link+1)%t(chi1,ld,chi2))

CALL FormGamma2(Gammas(link+1)%t,V%m,chi1,chi2)
DEALLOCATE(U%m,V%m,S%v)

END SUBROUTINE TwoSiteOpL

SUBROUTINE TwoSiteOpR(link,Op2,Gammas,Lambdas,truncerr)
!
!Purpose: Perform the complex two-site operation Op2 on the sites neighboring link
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
COMPLEX(KIND=rKind) :: ThetaRS(SIZE(Gammas(link)%t,1)*SIZE(Gammas(link)%t,2), &
SIZE(Gammas(link+1)%t,2)*SIZE(Gammas(link+1)%t,3))
TYPE(matrix) :: U, V
TYPE(vector) :: S
INTEGER :: chi0, chi1, chi2, ld

chi0=Size(Gammas(link)%t,1)
chi1=Size(Gammas(link)%t,3)
chi2=Size(Gammas(link+1)%t,3)

CALL FormThetaR(Theta,Gammas(link)%t,Lambdas(link+2)%v,Gammas(link+1)%t)
CALL ThetaOperation(Op2,Theta)
CALL ReshapeTheta(Theta,ThetaRS)
CALL SVDTruncation(link,ThetaRS, S, U, V)

CALL FormLambda(Lambdas(link+1), truncerr, S%v)
ld=SIZE(Gammas(link)%t,2)
chi1=SIZE(Lambdas(link+1)%v)
DEALLOCATE(Gammas(link)%t)
ALLOCATE(Gammas(link)%t(chi0,ld,chi1))

CALL FormGamma1(Gammas(link)%t,U%m,chi0,chi1)
ld=SIZE(Gammas(link+1)%t,2)
DEALLOCATE(Gammas(link+1)%t)
ALLOCATE(Gammas(link+1)%t(chi1,ld,chi2))

CALL FormGamma2(Gammas(link+1)%t,V%m,chi1,chi2)
DEALLOCATE(U%m,V%m,S%v)

END SUBROUTINE TwoSiteOpR

SUBROUTINE TwoSiteOpEdge(link,Op2,Gammas,Lambdas,truncerr)
!
!Purpose: Perform the complex two-site operation Op2 on the sites neighboring link
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
COMPLEX(KIND=rKind) :: ThetaRS(SIZE(Gammas(link)%t,1)*SIZE(Gammas(link)%t,2), &
SIZE(Gammas(link+1)%t,2)*SIZE(Gammas(link+1)%t,3))
TYPE(matrix) :: U, V
TYPE(vector) :: S
INTEGER :: chi0, chi1, chi2, ld

chi0=Size(Gammas(link)%t,1)
chi1=Size(Gammas(link)%t,3)
chi2=Size(Gammas(link+1)%t,3)

CALL FormThetaEdge(Theta,Gammas(link)%t,Lambdas(link+1)%v,Gammas(link+1)%t)
CALL ThetaOperation(Op2,Theta)
CALL ReshapeTheta(Theta,ThetaRS)
CALL SVDTruncation(link,ThetaRS, S, U, V)

CALL FormLambda(Lambdas(link+1), truncerr, S%v)
ld=SIZE(Gammas(link)%t,2)
chi1=SIZE(Lambdas(link+1)%v)
DEALLOCATE(Gammas(link)%t)
ALLOCATE(Gammas(link)%t(chi0,ld,chi1))

CALL FormGamma1(Gammas(link)%t,U%m,chi0,chi1)
ld=SIZE(Gammas(link+1)%t,2)
DEALLOCATE(Gammas(link+1)%t)
ALLOCATE(Gammas(link+1)%t(chi1,ld,chi2))

CALL FormGamma2(Gammas(link+1)%t,V%m,chi1,chi2)
DEALLOCATE(U%m,V%m,S%v)

END SUBROUTINE TwoSiteOpEdge
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Number conserving method starts !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE AllocateIndexLR(indL, indR, BlockSize)
!
!Purpose: Allocate indices for the Block diagonal Theta
!		  Used in Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(matrixInt), POINTER :: indL(:), indR(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: k, NumOfBlocks

NumOfBlocks=SIZE(BlockSize,1)
ALLOCATE(indL(NumOfBlocks), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate indexLeft'
END IF 
ALLOCATE(indR(NumOfBlocks), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate indexRight'
END IF 
DO k=1,NumOfBlocks
	ALLOCATE(indL(k)%mi(BlockSize(k,2),2), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate indexLeft'
	END IF 
	ALLOCATE(indR(k)%mi(BlockSize(k,3),2), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate indexRight'
	END IF 
	indL(k)%mi=0
	indR(k)%mi=0
END DO
END SUBROUTINE AllocateIndexLR

SUBROUTINE DeallocateIndexLR(indL, indR)
!
!Purpose: Deallocate indicies for the Block diagonal Theta
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(matrixInt), POINTER :: indL(:), indR(:)
INTEGER :: k, NumOfBlocks

NumOfBlocks=SIZE(indL,1)
DO k=1,NumOfBlocks
	DEALLOCATE(indL(k)%mi, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate indL'
	END IF 
	DEALLOCATE(indR(k)%mi, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate indR'
	END IF 
END DO
DEALLOCATE(indL, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate indL'
END IF 
DEALLOCATE(indR, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate indR'
END IF 

END SUBROUTINE DeallocateIndexLR

SUBROUTINE AllocateBlockTheta(BlockTheta, BlockSize)
!
!Purpose: Allocate Block diagonal Theta based on the number and size of each block
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(matrix), POINTER :: BlockTheta(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: k, NumOfBlocks

NumOfBlocks=SIZE(BlockSize,1)
ALLOCATE(BlockTheta(NumOfBlocks))
DO k=1,NumOfBlocks
	!If the block has both left and right Schmidt vectors (i.e. is of nonzero dimension) 
	!then allocate to full size
	IF((BlockSize(k,2).ne.0).AND.(BlockSize(k,3).ne.0)) THEN
		ALLOCATE(BlockTheta(k)%m(BlockSize(k,2),BlockSize(k,3)), STAT=statInt)
		IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate BlockTheta'
		END IF 
	ELSE
		!ELSE allocate a 1 by 1 block
		ALLOCATE(BlockTheta(k)%m(1,1), STAT=statInt)
		IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate BlockTheta'
		END IF 
	END IF
	BlockTheta(k)%m=CMPLX(0.0,KIND=rKind)
END DO
END SUBROUTINE AllocateBlockTheta

SUBROUTINE DeallocateBlockTheta(BlockTheta)
!
!Purpose: Deallocate Block diagonal Theta based on the number and size of each block
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(matrix), POINTER :: BlockTheta(:)
INTEGER :: k, NumOfBlocks

NumOfBlocks=SIZE(BlockTheta,1)
DO k=1,NumOfBlocks
	DEALLOCATE(BlockTheta(k)%m, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate BlockTheta'
	END IF 
END DO
DEALLOCATE(BlockTheta, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate BlockTheta'
END IF 

END SUBROUTINE DeallocateBlockTheta

SUBROUTINE AllocateUSV(US, SS, VS, BlockSize)
!
!Purpose: Allocate stuff for the singular value decomposition LAPACK routine
!		  Specific to Block diagonal Theta 
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(matrix), POINTER :: US(:), VS(:)
TYPE(vector), POINTER :: SS(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: k, NumOfBlocks

NumOfBlocks=SIZE(BlockSize,1)
!Allocate a list of SVD variables
!The index corresponds to a "block" of the block diagonal theta
ALLOCATE(US(NumOfBlocks))
ALLOCATE(SS(NumOfBlocks))
ALLOCATE(VS(NumOfBlocks))
DO k=1,NumOfBlocks
	!Allocate U, S, and V only if the block is greater than 1X1 in size
	IF((BlockSize(k,2).ne.0).AND.(BlockSize(k,3).ne.0)) THEN
		ALLOCATE(US(k)%m(BlockSize(k,2),BlockSize(k,2)), STAT=statInt)
		IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate US'
		END IF 
		ALLOCATE(SS(k)%v(MIN(BlockSize(k,2),BlockSize(k,3))), STAT=statInt)
		IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate SS'
		END IF 
		ALLOCATE(VS(k)%m(BlockSize(k,3),BlockSize(k,3)), STAT=statInt)
		IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate VS'
		END IF 
	ELSE
		ALLOCATE(US(k)%m(1,1), STAT=statInt)
		IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate US'
		END IF 
		ALLOCATE(SS(k)%v(1), STAT=statInt)
		IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate SS'
		END IF 
		ALLOCATE(VS(k)%m(1,1), STAT=statInt)
		IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate VS'
		END IF 
	END IF
	US(k)%m=CMPLX(0.0,KIND=rKind)
	SS(k)%v=0.0_rKind
	VS(k)%m=CMPLX(0.0,KIND=rKind)
END DO

END SUBROUTINE AllocateUSV

SUBROUTINE DeallocateUSV(US, SS, VS)
!
!Purpose: Allocate stuff for the singular value decomposition LAPACK routine
!		  Specific to Block diagonal Theta 
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(matrix), POINTER :: US(:), VS(:)
TYPE(vector), POINTER :: SS(:)
INTEGER :: k, NumOfBlocks

NumOfBlocks=SIZE(SS,1)
DO k=1,NumOfBlocks
	DEALLOCATE(US(k)%m, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate US'
	END IF 
	DEALLOCATE(SS(k)%v, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate SS'
	END IF 
	DEALLOCATE(VS(k)%m, STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate VS'
	END IF 
END DO
DEALLOCATE(US, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate US'
END IF 
DEALLOCATE(SS, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate SS'
END IF 
DEALLOCATE(VS, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate VS'
END IF 

END SUBROUTINE DeallocateUSV
	
SUBROUTINE AllocateSSflat(ssflat, BlockSize)
!
!Purpose: Allocate vector to hold the singular values from all blocks
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(vector) :: ssflat
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: k, NumOfBlocks, SizeOfSS

NumOfBlocks=SIZE(BlockSize,1)
SizeOfSS=0
DO k=1, NumOfBlocks
	!There are min(m,n) singular values of an mXn matrix
	!Count the singular values from each block
	SizeOfSS=SizeOfSS+MIN(BlockSize(k,2),BlockSize(k,3))
END DO
ALLOCATE(ssflat%v(SizeOfSS), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate ssflat'
END IF 
ssflat%v=0.0_rKind

END SUBROUTINE AllocateSSflat
	

SUBROUTINE OneSiteOp_rtm(Op1,Gamma,Gammat, LabelLeft, LabelLeftp)
!
!Purpose: Perform the real one-site operation Op1 on the Gamma specified.
! Because the new MPS tensor Op1*Gamma may not transform irreducibly (or
! in the same irrep) we store it in its own tensor.  This routine is
! used primarily in GKernalC and GContractionC
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: Op1(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma(:,:)
TYPE(tensor) :: Gammat
INTEGER, INTENT(INOUT) :: LabelLeft(:), LabelLeftp(:)
INTEGER :: alpha, i,beta, loci

ALLOCATE(Gammat%t(SIZE(Gamma,1),localSize,SIZE(Gamma,2)))

Gammat%t=0.0_rKind

DO alpha = 1,SIZE(Gamma,1)
	IF(LabelLeft(alpha).lt.10000) THEN
	DO beta=1,SIZE(Gamma,2)
		IF(LabelLeftp(beta).lt.10000) THEN
			loci=LabelLeftp(beta)-LabelLeft(alpha)+1
			IF((1.le.loci).and.(loci.le.localSize)) THEN
			DO i=1,localSize
				Gammat%t(alpha,i,beta)=Gamma(alpha,beta)*Op1(i,loci)
			END DO
			END IF
		ELSE
			EXIT
		END IF
	END DO
	ELSE
		EXIT
	END IF

END DO		
END SUBROUTINE OneSiteOp_rtm

SUBROUTINE OneSiteOp_ctm(Op1,Gamma,Gammat, LabelLeft, LabelLeftp)
!
!Purpose: Perform the complex one-site operation Op1 on the Gamma specified
! Because the new MPS tensor Op1*Gamma may not transform irreducibly (or
! in the same irrep) we store it in its own tensor.  This routine is
! used primarily in GKernalC and GContractionC
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma(:,:)
TYPE(tensor) :: Gammat
INTEGER, INTENT(INOUT) :: LabelLeft(:), LabelLeftp(:)
INTEGER :: alpha, i,beta, loci

ALLOCATE(Gammat%t(SIZE(Gamma,1),localSize,SIZE(Gamma,2)))

Gammat%t=0.0_rKind

DO alpha = 1,SIZE(Gamma,1)
	IF(LabelLeft(alpha).lt.10000) THEN
	DO beta=1,SIZE(Gamma,2)
		IF(LabelLeftp(beta).lt.10000) THEN
			loci=LabelLeftp(beta)-LabelLeft(alpha)+1
			IF((1.le.loci).and.(loci.le.localSize)) THEN
				DO i=1,localSize
					Gammat%t(alpha,i,beta)=Gamma(alpha,beta)*Op1(i,loci)
				END DO
			END IF
		ELSE
			EXIT
		END IF
	END DO
	ELSE
		EXIT
	END IF

END DO	
	
END SUBROUTINE OneSiteOp_ctm

SUBROUTINE OneSiteOp_rm(Op1,Gamma, LabelLeft, LabelLeftp)
!
!Purpose: Perform the real one-site operation Op1 on the Gamma specified.
! This routine assumes that the operator is diagonal in the on-site number
! eigenbasis so that the Gamma matrices can remain as matrices
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: Op1(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma(:,:)
INTEGER, INTENT(INOUT) :: LabelLeft(:), LabelLeftp(:)
INTEGER :: alpha, i,beta, loci


DO alpha = 1,SIZE(Gamma,1)
	IF(LabelLeft(alpha).lt.10000) THEN
	DO beta=1,SIZE(Gamma,2)
		IF(LabelLeftp(beta).lt.10000) THEN
			loci=LabelLeftp(beta)-LabelLeft(alpha)+1
			IF((1.le.loci).and.(loci.le.localSize)) THEN
				Gamma(alpha,beta)=Gamma(alpha,beta)*Op1(loci,loci)
			END IF
		ELSE
			EXIT
		END IF
	END DO
	ELSE
		EXIT
	END IF

END DO		
END SUBROUTINE OneSiteOp_rm

SUBROUTINE OneSiteOp_cm(Op1,Gamma, LabelLeft, LabelLeftp)
!
!Purpose: Perform the real one-site operation Op1 on the Gamma specified.
! This routine assumes that the operator is diagonal in the on-site number
! eigenbasis so that the Gamma matrices can remain as matrices
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma(:,:)
INTEGER, INTENT(INOUT) :: LabelLeft(:), LabelLeftp(:)
INTEGER :: alpha, i,beta, loci


DO alpha = 1,SIZE(Gamma,1)
	IF(LabelLeft(alpha).lt.10000) THEN
	DO beta=1,SIZE(Gamma,2)
		IF(LabelLeftp(beta).lt.10000) THEN
			loci=LabelLeftp(beta)-LabelLeft(alpha)+1
			IF((1.le.loci).and.(loci.le.localSize)) THEN
				Gamma(alpha,beta)=Gamma(alpha,beta)*Op1(loci,loci)
			END IF
		ELSE
			EXIT
		END IF
	END DO
	ELSE
		EXIT
	END IF

END DO	
END SUBROUTINE OneSiteOp_cm


SUBROUTINE SVDInitNC(k, BlockSize)
!
!Purpose: Allocate variables for SVD
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: k
INTEGER, INTENT(IN) :: BlockSize(:,:)

IF((BlockSize(k,2).ne.0).and.(BlockSize(k,3).ne.0)) THEN
	jobu_SVD='A'
	jobvt_SVD='A'
	matrixSizeL_SVD=BlockSize(k,2)
	matrixSizeT_SVD=BlockSize(k,3)
	workSize_SVD=5*MAX(matrixSizeL_SVD,matrixSizeT_SVD)
	ALLOCATE(work_SVD(workSize_SVD), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate SVDNC variables'
	END IF
	ALLOCATE(rwork_SVD(5*MIN(matrixSizeL_SVD,matrixSizeT_SVD)), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate SVDNC variables'
	END IF
ELSE
	jobu_SVD='A'
	jobvt_SVD='A'
	matrixSizeL_SVD=1
	matrixSizeT_SVD=1
	workSize_SVD=5*MAX(matrixSizeL_SVD,matrixSizeT_SVD)
	ALLOCATE(work_SVD(workSize_SVD), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate SVDNC variables'
	END IF
	ALLOCATE(rwork_SVD(5*MIN(matrixSizeL_SVD,matrixSizeT_SVD)), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate SVDNC variables'
	END IF
END IF
END SUBROUTINE SVDInitNC

SUBROUTINE SVDFinishNC()
!
!Purpose: Deallocate variables for SVD
!
IMPLICIT NONE

DEALLOCATE(work_SVD, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate SVDNC variables'
END IF
DEALLOCATE(rwork_SVD, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate SVDNC variables'
END IF
END SUBROUTINE SVDFinishNC

SUBROUTINE FormThetaNCEdge_t(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
!
!Purpose: Form Theta consistent with number conservation
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind), INTENT(OUT) :: Theta(:,:,:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: chi, i, j, alpha, beta, gamma, lftN1, lftN2, rgtN1, rgtN2, ni, nj
	
		
chi=SIZE(Lambdas(link+1)%v,1)
! Initialize Theta
Theta=CMPLX(0.0,KIND=rKind)
			
!If internal degrees of freedom are present, we must sum over them with
!a constraint that number is conserved		
!$OMP PARALLEL DO PRIVATE(gamma,beta,i,j,ni,nj,lftN1,lftN2,rgtN1,rgtN2) schedule(dynamic)
DO alpha=1,SIZE(Gammas(link)%t,1)
	IF(LabelLeft(link)%vi(alpha).lt.1000) THEN
	DO gamma=1,SIZE(Gammas(link+1)%t,3)
		IF(LabelLeft(link+2)%vi(gamma).lt.1000) THEN
		DO beta=1,SIZE(Gammas(link)%t,3)
!We still have to sum over internal degrees of freedom
			DO i=1,SIZE(Gammas(link)%t,2),1
!number on left site
			ni=Conserv%vi(i)
			DO j=1,SIZE(Gammas(link+1)%t,2),1
!number on right site
			nj=Conserv%vi(j)
! Number of particles on the left side of 'link'.
					lftN1=LabelLeft(link)%vi(alpha)
! Number of particles on the left side of 'link+1'.
					lftN2=LabelLeft(link+2)%vi(gamma)
! So that lambda is not equal to zero
!Label*=1000 means that it has been left uninitialized (incoming state does not have enough entanglement)
						IF((MAX(lftN1,lftN2,rgtN1,rgtN2)<1000) &
! So that lftN1+n_i=lftN2 and rgtN1+n_j=rgtN2
						.AND.(lftN1+ni+nj==lftN2)) THEN

						Theta(alpha,i,j,gamma) = &
						Theta(alpha,i,j,gamma) + &
						Gammas(link)%t(alpha,i,beta) * &
						Lambdas(link+1)%v(beta) * Gammas(link+1)%t(beta,j,gamma)
						END IF
			END DO
			END DO
		END DO
		END IF
	END DO
	END IF
END DO	
!$OMP END PARALLEL DO	

END SUBROUTINE FormThetaNCEdge_t

SUBROUTINE FormThetaNCEdge_m(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
!
!Purpose: Form Theta consistent with number conservation
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind), INTENT(OUT) :: Theta(:,:,:,:)
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: chi, i, j, alpha, beta, gamma, lftN1, lftN2, rgtN1, rgtN2, ni, nj
	
		
chi=SIZE(Lambdas(link+1)%v,1)
! Initialize Theta
Theta=CMPLX(0.0,KIND=rKind)
			
!$OMP PARALLEL DO PRIVATE(gamma,beta,lftN1,lftN2,rgtN1,rgtN2) schedule(dynamic)
DO alpha=1,SIZE(Gammas(link)%m,1)
	IF(LabelLeft(link)%vi(alpha).lt.10000) THEN
	DO gamma=1,SIZE(Gammas(link+1)%m,2)
		IF(LabelLeft(link+2)%vi(gamma).lt.10000) THEN
		DO beta=1,SIZE(Gammas(link)%m,2)
			IF(LabelLeft(link+1)%vi(beta).lt.10000) THEN
! Number of particles on the left side of 'link'.
			lftN1=LabelLeft(link)%vi(alpha)
! Number of particles on the left side of 'link+1'.
			lftN2=LabelLeft(link+1)%vi(beta)
! Number of particles on the right side of 'link+2'.
			rgtN1=LabelRight(link+2)%vi(gamma)
! Number of particles on the right side of 'link+1'.
			rgtN2=LabelRight(link+1)%vi(beta)
! So that lambda is not equal to zero
			IF((lftN2-lftN1+1<=localSize).AND.(lftN2-lftN1+1>=1) &
! 1<=j<=d
			.AND.(rgtN2-rgtN1+1<=localSize).AND.(rgtN2-rgtN1+1>=1)) THEN
! So that lftN1+i-1=lftN2 and rgtN1+j-1=rgtN2
				Theta(alpha,lftN2-lftN1+1,rgtN2-rgtN1+1,gamma) = &
				Theta(alpha,lftN2-lftN1+1,rgtN2-rgtN1+1,gamma) + &
				Gammas(link)%m(alpha,beta) * &
				Lambdas(link+1)%v(beta) * Gammas(link+1)%m(beta,gamma)
			END IF
			END IF
		END DO
		END IF
	END DO
	END IF
END DO
!$OMP END PARALLEL DO
		
END SUBROUTINE FormThetaNCEdge_m

SUBROUTINE FormThetaNCR_t(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
!
!Purpose: Form Theta consistent with number conservation
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind), INTENT(OUT) :: Theta(:,:,:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: chi, i, j, alpha, beta, gamma, lftN1, lftN2, rgtN1, rgtN2, ni, nj
	
		
chi=SIZE(Lambdas(link+1)%v,1)
! Initialize Theta
Theta=CMPLX(0.0,KIND=rKind)
			
!If internal degrees of freedom are present, we must sum over them with
!a constraint that number is conserved		
!$OMP PARALLEL DO PRIVATE(gamma,beta,i,j,ni,nj,lftN1,lftN2,rgtN1,rgtN2) schedule(dynamic)
DO alpha=1,SIZE(Gammas(link)%t,1)
	IF(LabelLeft(link)%vi(alpha).lt.1000) THEN
	DO gamma=1,SIZE(Gammas(link+1)%t,3)
		IF(LabelLeft(link+2)%vi(gamma).lt.1000) THEN
		DO beta=1,SIZE(Gammas(link)%t,3)
!We still have to sum over internal degrees of freedom
			DO i=1,SIZE(Gammas(link)%t,2),1
!number on left site
			ni=Conserv%vi(i)
			DO j=1,SIZE(Gammas(link+1)%t,2),1
!number on right site
			nj=Conserv%vi(j)
! Number of particles on the left side of 'link'.
					lftN1=LabelLeft(link)%vi(alpha)
! Number of particles on the left side of 'link+1'.
					lftN2=LabelLeft(link+2)%vi(gamma)
! So that lambda is not equal to zero
!Label*=1000 means that it has been left uninitialized (incoming state does not have enough entanglement)
						IF((MAX(lftN1,lftN2,rgtN1,rgtN2)<1000) &
! So that lftN1+n_i=lftN2 and rgtN1+n_j=rgtN2
						.AND.(lftN1+ni+nj==lftN2)) THEN

						Theta(alpha,i,j,gamma) = &
						Theta(alpha,i,j,gamma) + &
						Gammas(link)%t(alpha,i,beta)  &
						 * Gammas(link+1)%t(beta,j,gamma)*Lambdas(link+2)%v(gamma)
						END IF
			END DO
			END DO
		END DO
		END IF
	END DO
	END IF
END DO	
!$OMP END PARALLEL DO	

END SUBROUTINE FormThetaNCR_t

SUBROUTINE FormThetaNCR_m(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
!
!Purpose: Form Theta consistent with number conservation
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind), INTENT(OUT) :: Theta(:,:,:,:)
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: chi, i, j, alpha, beta, gamma, lftN1, lftN2, rgtN1, rgtN2, ni, nj
	
		
chi=SIZE(Lambdas(link+1)%v,1)
! Initialize Theta
Theta=CMPLX(0.0,KIND=rKind)
			
!$OMP PARALLEL DO PRIVATE(gamma,beta,lftN1,lftN2,rgtN1,rgtN2) schedule(dynamic)
DO alpha=1,SIZE(Gammas(link)%m,1)
	IF(LabelLeft(link)%vi(alpha).lt.10000) THEN
	DO gamma=1,SIZE(Gammas(link+1)%m,2)
		IF(LabelLeft(link+2)%vi(gamma).lt.10000) THEN
		DO beta=1,SIZE(Gammas(link)%m,2)
			IF(LabelLeft(link+1)%vi(beta).lt.10000) THEN
! Number of particles on the left side of 'link'.
			lftN1=LabelLeft(link)%vi(alpha)
! Number of particles on the left side of 'link+1'.
			lftN2=LabelLeft(link+1)%vi(beta)
! Number of particles on the right side of 'link+2'.
			rgtN1=LabelRight(link+2)%vi(gamma)
! Number of particles on the right side of 'link+1'.
			rgtN2=LabelRight(link+1)%vi(beta)
! So that lambda is not equal to zero
			IF((lftN2-lftN1+1<=localSize).AND.(lftN2-lftN1+1>=1) &
! 1<=j<=d
			.AND.(rgtN2-rgtN1+1<=localSize).AND.(rgtN2-rgtN1+1>=1)) THEN
! So that lftN1+i-1=lftN2 and rgtN1+j-1=rgtN2
				Theta(alpha,lftN2-lftN1+1,rgtN2-rgtN1+1,gamma) = &
				Theta(alpha,lftN2-lftN1+1,rgtN2-rgtN1+1,gamma) + &
				Gammas(link)%m(alpha,beta)  &
				 * Gammas(link+1)%m(beta,gamma)*Lambdas(link+2)%v(gamma)
			END IF
			END IF
		END DO
		END IF
	END DO
	END IF
END DO
!$OMP END PARALLEL DO
		
END SUBROUTINE FormThetaNCR_m

SUBROUTINE FormThetaNCL_t(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
!
!Purpose: Form Theta consistent with number conservation
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind), INTENT(OUT) :: Theta(:,:,:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: chi, i, j, alpha, beta, gamma, lftN1, lftN2, rgtN1, rgtN2, ni, nj
	
		
chi=SIZE(Lambdas(link+1)%v,1)
! Initialize Theta
Theta=CMPLX(0.0,KIND=rKind)
			
!If internal degrees of freedom are present, we must sum over them with
!a constraint that number is conserved		
!$OMP PARALLEL DO PRIVATE(gamma,beta,i,j,ni,nj,lftN1,lftN2,rgtN1,rgtN2) schedule(dynamic)
DO alpha=1,SIZE(Gammas(link)%t,1)
	IF(LabelLeft(link)%vi(alpha).lt.1000) THEN
	DO gamma=1,SIZE(Gammas(link+1)%t,3)
		IF(LabelLeft(link+2)%vi(gamma).lt.1000) THEN
		DO beta=1,SIZE(Gammas(link)%t,3)
!We still have to sum over internal degrees of freedom
			DO i=1,SIZE(Gammas(link)%t,2),1
!number on left site
			ni=Conserv%vi(i)
			DO j=1,SIZE(Gammas(link+1)%t,2),1
!number on right site
			nj=Conserv%vi(j)
! Number of particles on the left side of 'link'.
					lftN1=LabelLeft(link)%vi(alpha)
! Number of particles on the left side of 'link+1'.
					lftN2=LabelLeft(link+2)%vi(gamma)
! So that lambda is not equal to zero
!Label*=1000 means that it has been left uninitialized (incoming state does not have enough entanglement)
						IF((MAX(lftN1,lftN2,rgtN1,rgtN2)<1000) &
! So that lftN1+n_i=lftN2 and rgtN1+n_j=rgtN2
						.AND.(lftN1+ni+nj==lftN2)) THEN

						Theta(alpha,i,j,gamma) = &
						Theta(alpha,i,j,gamma) + &
						Lambdas(link)%v(alpha)*Gammas(link)%t(alpha,i,beta)  &
						 * Gammas(link+1)%t(beta,j,gamma)
						END IF
			END DO
			END DO
		END DO
		END IF
	END DO
	END IF
END DO	
!$OMP END PARALLEL DO	

END SUBROUTINE FormThetaNCL_t

SUBROUTINE FormThetaNCL_m(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
!
!Purpose: Form Theta consistent with number conservation
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind), INTENT(OUT) :: Theta(:,:,:,:)
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: chi, i, j, alpha, beta, gamma, lftN1, lftN2, rgtN1, rgtN2, ni, nj
	
		
chi=SIZE(Lambdas(link+1)%v,1)
! Initialize Theta
Theta=CMPLX(0.0,KIND=rKind)
			
!$OMP PARALLEL DO PRIVATE(gamma,beta,lftN1,lftN2,rgtN1,rgtN2) schedule(dynamic)
DO alpha=1,SIZE(Gammas(link)%m,1)
	IF(LabelLeft(link)%vi(alpha).lt.10000) THEN
	DO gamma=1,SIZE(Gammas(link+1)%m,2)
		IF(LabelLeft(link+2)%vi(gamma).lt.10000) THEN
		DO beta=1,SIZE(Gammas(link)%m,2)
			IF(LabelLeft(link+1)%vi(beta).lt.10000) THEN
! Number of particles on the left side of 'link'.
			lftN1=LabelLeft(link)%vi(alpha)
! Number of particles on the left side of 'link+1'.
			lftN2=LabelLeft(link+1)%vi(beta)
! Number of particles on the right side of 'link+2'.
			rgtN1=LabelRight(link+2)%vi(gamma)
! Number of particles on the right side of 'link+1'.
			rgtN2=LabelRight(link+1)%vi(beta)
! So that lambda is not equal to zero
			IF((lftN2-lftN1+1<=localSize).AND.(lftN2-lftN1+1>=1) &
! 1<=j<=d
			.AND.(rgtN2-rgtN1+1<=localSize).AND.(rgtN2-rgtN1+1>=1)) THEN
! So that lftN1+i-1=lftN2 and rgtN1+j-1=rgtN2
				Theta(alpha,lftN2-lftN1+1,rgtN2-rgtN1+1,gamma) = &
				Theta(alpha,lftN2-lftN1+1,rgtN2-rgtN1+1,gamma) + &
				Lambdas(link)%v(alpha)*Gammas(link)%m(alpha,beta)  &
				 * Gammas(link+1)%m(beta,gamma)
			END IF
			END IF
		END DO
		END IF
	END DO
	END IF
END DO
!$OMP END PARALLEL DO
		
END SUBROUTINE FormThetaNCL_m

SUBROUTINE ThetaOperationNC(link, Op2, Theta, LabelLeft, LabelRight, intDegFree)
!
!Purpose:Operate on Theta with the real two-site operation Op2 consistent with number conservation
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind) :: ThetaTemp(SIZE(Theta,1),SIZE(Theta,2),SIZE(Theta,3),SIZE(Theta,4))
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: chi, i,j,l, k, alpha, gamma, lftN, rgtN,nk,nj,nl,ni
		
chi=SIZE(Theta,1)
ThetaTemp=CMPLX(0.0,KIND=rKind)
		
!Internal degrees of freedom present	
IF(PRESENT(intDegFree)) THEN
!$OMP PARALLEL PRIVATE(gamma,k,i,ni,nk,j,l,nj,nl,lftN,rgtN)
!$OMP DO schedule(dynamic)
	DO alpha=1,SIZE(Theta,1)
		IF(LabelLeft(link)%vi(alpha).lt.1000) THEN	
		DO gamma=1,SIZE(Theta,4)
			IF(LabelRight(link+2)%vi(gamma).lt.1000) THEN
			DO k=1,SIZE(Theta,2)
				nk=Conserv%vi(k)
				DO j=1, SIZE(Theta,3)
				nj=Conserv%vi(j)
! Number of particles on the left side of 'link'.
					lftN=LabelLeft(link)%vi(alpha)
! Number of particles on the right side of 'link+2'.
					rgtN=LabelRight(link+2)%vi(gamma)
! So that corresponding lambda is not equal to zero.
					IF((MAX(lftN,rgtN)<1000).AND. &
					(nj==totQ-lftN-rgtN-nk)) THEN
						ThetaTemp(alpha,k,j,gamma)=CMPLX(0.0,KIND=rKind)
						DO i=1,SIZE(Theta,2)
						ni=Conserv%vi(i)
						DO l=1,SIZE(Theta,3)
						nl=Conserv%vi(l)
! 0<=n_i<=maxfilling
							IF(nl==totQ-lftN-rgtN-ni) THEN
								ThetaTemp(alpha,k,j,gamma) = &
								ThetaTemp(alpha,k,j,gamma) + &
								Op2((k-1)*SIZE(Theta,2)+j, (i-1)*SIZE(Theta,3)+l) * &
								Theta(alpha,i,l,gamma)
							END IF
						END DO
						END DO
					END IF
				END DO
			END DO
			END IF
		END DO
		END IF
	END DO
!$OMP END DO
!$OMP END PARALLEL
!Internal degrees of freedom absent	
ELSE	
!$OMP PARALLEL PRIVATE(gamma,k,i,lftN,rgtN)
!$OMP DO schedule(dynamic)		
		DO alpha=1,SIZE(Theta,1)
			IF(LabelLeft(link)%vi(alpha).lt.1000) THEN	
			DO gamma=1,SIZE(Theta,4)
				IF(LabelRight(link+2)%vi(gamma).lt.1000) THEN
				DO k=1,SIZE(Theta,2)
! Number of particles on the left side of 'link'.
					lftN=LabelLeft(link)%vi(alpha)
! Number of particles on the right side of 'link+2'.
					rgtN=LabelRight(link+2)%vi(gamma)
! So that corresponding lambda is not equal to zero.
					IF((MAX(lftN,rgtN)<1000).AND. &
! 0<=k<=maxfilling
					(totQ-lftN-rgtN-k+1>=0).AND.(totQ-lftN-rgtN-k+1<=Nmax)) THEN
						ThetaTemp(alpha,k,totQ-lftN-rgtN-k+2,gamma)=CMPLX(0.0,KIND=rKind)
						DO i=1,SIZE(Theta,2)
! 0<=i<=maxfilling
							IF((totQ-lftN-rgtN-i+1>=0).AND.(totQ-lftN-rgtN-i+1<=Nmax)) THEN
! Due to number conservation, l and j have to be equal to totQ-lftN-rgtN-k+1 and totQ-lftN-rgtN-i+1.
								ThetaTemp(alpha,k,totQ-lftN-rgtN-k+2,gamma) = &
								ThetaTemp(alpha,k,totQ-lftN-rgtN-k+2,gamma) + &
								Op2((k-1)*SIZE(Theta,2)+totQ-lftN-rgtN-k+2, (i-1)*SIZE(Theta,3)+totQ-lftN-rgtN-i+2) * &
								Theta(alpha,i,totQ-lftN-rgtN-i+2,gamma)
							END IF
						END DO
					END IF
				END DO
				END IF
			END DO
			END IF
		END DO
!$OMP END DO
!$OMP END PARALLEL	
END IF
Theta=ThetaTemp

END SUBROUTINE ThetaOperationNC
	
SUBROUTINE minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
!
!Purpose:Find the minimum and maximum values of the number on the left=N_L(\alpha)+N_S(i) and on the right=N_R(\alpha )+N_S(j)
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: link
INTEGER, INTENT(OUT) :: minNL, maxNL, minNR, maxNR
INTEGER :: chi, alpha, gamma
INTEGER :: Label0(MAXVAL((/SIZE(LabelLeft(link)%vi,1),SIZE(LabelRight(link+2)%vi,1)/)))
		
chi=SIZE(LabelLeft(link)%vi,1)
!as few as 0 particles on site i
minNL=MINVAL(LabelLeft(link)%vi)
Label0=0
DO alpha=1,SIZE(LabelLeft(link)%vi)
	!Label=10000 means the state is notinitialized (not enough entanglement to need it)
	IF(LabelLeft(link)%vi(alpha)==10000) THEN
		Label0(alpha)=0
	ELSE
		Label0(alpha)=LabelLeft(link)%vi(alpha)
	END IF
END DO
!up to maxfilling particles allowed on site i
maxNL=MAXVAL(Label0)+Nmax
!as few as 0 particles on site j
minNR=MINVAL(LabelRight(link+2)%vi)
Label0=0
DO gamma=1,SIZE(LabelRight(link+2)%vi)
	!Label=10000 means the state is notinitialized (not enough entanglement to need it)
	IF(LabelRight(link+2)%vi(gamma)==10000) THEN
		Label0(gamma)=0
	ELSE
		Label0(gamma)=LabelRight(link+2)%vi(gamma)
	END IF
END DO
!up to maxfilling particles allowed on site j
maxNR=MAXVAL(Label0)+Nmax
		
!If the onsite maximum number puts us over the total number, make the total number our bound
IF(maxNL>totQ-minNR) THEN
	maxNL=totQ-minNR
ELSE
	minNR=totQ-maxNL
END IF
		
IF(maxNR>totQ-minNL) THEN
	maxNR=totQ-minNL
ELSE
	minNL=totQ-maxNR
END IF
END SUBROUTINE minmaxNLR
	
SUBROUTINE SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight, intDegFree)
!
!Purpose:The reshaped Theta will be block diagonal due to number conservation, 
!with each block corresponding to a fixed number on the left.  This subroutine finds how many
!blocks exist and their sizes
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(OUT) :: BlockSize(:,:)
INTEGER, INTENT(IN) :: link, minNL, maxNL, minNR, maxNR
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: chi, alpha, gamma, i, j, k, NLeft, NRight, ni,nj

chi=SIZE(LabelLeft(link)%vi,1)
BlockSize=0
!least number on the left corresponds to index 1
!Blocksize(i,1)=(number on the left)_i
DO k=minNL,maxNL,1
	BlockSize(k-minNL+1,1)=k
END DO

!Internal degrees of freedom present	
IF(PRESENT(intDegFree)) THEN
!Count the size of each block (of fixed number on the left) on the left		
	DO alpha=1,SIZE(LabelLeft(link)%vi)
		If(LabelLeft(link)%vi(alpha)<1000) THEN
			DO i=1,localSize
				ni=Conserv%vi(i)
				NLeft=LabelLeft(link)%vi(alpha)+ni
				IF((minNL<=NLeft).AND.(NLeft<=maxNL)) THEN
					BlockSize(Nleft-minNL+1,2)=BlockSize(Nleft-minNL+1,2)+1
				END IF
			END DO
		END IF					
	END DO

!Count the size of each block (of fixed number on the left) on the right		
	DO gamma=1,SIZE(LabelRight(link+2)%vi)
		IF(LabelRight(link+2)%vi(gamma)<1000) THEN
			DO j=1,localSize
				nj=Conserv%vi(j)

				NRight=LabelRight(link+2)%vi(gamma)+nj
				IF((minNR<=NRight).AND.(NRight<=maxNR)) THEN
					BlockSize(totQ-NRight-minNL+1,3)=BlockSize(totQ-NRight-minNL+1,3)+1
				END IF
			END DO
		END IF
	END DO

!Internal degrees of freedom absent	
ELSE

!Count the size of each block (of fixed number on the left) on the left				
	DO alpha=1,SIZE(LabelLeft(link)%vi)
		If(LabelLeft(link)%vi(alpha)<1000) THEN
			DO i=1,localSize
				NLeft=LabelLeft(link)%vi(alpha)+i-1
				IF((minNL<=NLeft).AND.(NLeft<=maxNL)) THEN
					BlockSize(Nleft-minNL+1,2)=BlockSize(Nleft-minNL+1,2)+1
				END IF
			END DO
		END IF					
	END DO
!Count the size of each block (of fixed number on the left) on the right		
	DO gamma=1,SIZE(LabelRight(link+2)%vi)
		IF(LabelRight(link+2)%vi(gamma)<1000) THEN
		DO j=1,localSize
				NRight=LabelRight(link+2)%vi(gamma)+j-1
				IF((minNR<=NRight).AND.(NRight<=maxNR)) THEN
					BlockSize(totQ-NRight-minNL+1,3)=BlockSize(totQ-NRight-minNL+1,3)+1
				END IF
			END DO
		END IF
	END DO
		
END IF
END SUBROUTINE SizeOfBlocks
	
	
SUBROUTINE IndexLeft(link, indL, minNL, maxNL, LabelLeft, intDegFree)
!
!Purpose: Find the on-site indices i and left schmidt indices alpha that 
!correspond to the fixed number on the left indexed by indind
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:)
TYPE(matrixInt), POINTER :: indL(:)
INTEGER, INTENT(IN) :: link, minNL, maxNL
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: chi, i, alpha, indind,ni
INTEGER :: turnL(maxNL-minNL+1)

chi=SIZE(LabelLeft(link)%vi,1)
turnL=1

!Internal degrees of freedom present	
IF(PRESENT(intDegFree)) THEN
	DO alpha=1,SIZE(LabelLeft(link)%vi)
		IF(LabelLeft(link)%vi(alpha)<1000) THEN
		DO i=1,localSize
			ni=Conserv%vi(i)
!Check that the number on the left is an allowed value
			IF((minNL<=LabelLeft(link)%vi(alpha)+ni) &
			.AND.(LabelLeft(link)%vi(alpha)+ni<=maxNL)) THEN
!This index is the number
				indind=LabelLeft(link)%vi(alpha)+ni-minNL+1
!Index of the schmidt vector
				indL(indind)%mi(turnL(indind),1)=alpha
!Index of the on-site state
				indL(indind)%mi(turnL(indind),2)=i
!Next allowed state index
				turnL(indind)=turnL(indind)+1
			END IF
		END DO
		END IF
	END DO
!Internal degrees of freedom absent
ELSE	

	DO alpha=1,SIZE(LabelLeft(link)%vi)
		IF(LabelLeft(link)%vi(alpha)<1000) THEN
		DO i=1,localSize
			IF((minNL<=LabelLeft(link)%vi(alpha)+i-1) &
				.AND.(LabelLeft(link)%vi(alpha)+i-1<=maxNL)) THEN
			indind=LabelLeft(link)%vi(alpha)+i-1-minNL+1
			indL(indind)%mi(turnL(indind),1)=alpha
			indL(indind)%mi(turnL(indind),2)=i
			turnL(indind)=turnL(indind)+1
			END IF
		END DO
		END IF
	END DO
END IF


END SUBROUTINE IndexLeft
	
SUBROUTINE IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight,intDegFree)
!
!Purpose: Find the on-site indices j and right schmidt gamma alpha that 
!correspond to the fixed number on the left indexed by indind
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelRight(:)
TYPE(matrixInt), POINTER :: indR(:)
INTEGER, INTENT(IN) :: link, minNL, maxNL, minNR, maxNR
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: chi, j, gamma, indind,nj
INTEGER :: turnR(maxNL-minNL+1)
	chi=SIZE(LabelRight(link+2)%vi,1)
	turnR=1

!Internal degrees of freedom present	
IF(PRESENT(intDegFree)) THEN
	DO gamma=1,SIZE(LabelRight(link+2)%vi)
		IF(LabelRight(link+2)%vi(gamma)<1000) THEN
			DO j=1,localSize
			nj=Conserv%vi(j)
!Same as above, now for the right
				IF((minNR<=LabelRight(link+2)%vi(gamma)+nj) &
					.AND.(LabelRight(link+2)%vi(gamma)+nj<=maxNR)) THEN
				indind=totQ-(LabelRight(link+2)%vi(gamma)+nj)-minNL+1
				indR(indind)%mi(turnR(indind),1)=gamma
				indR(indind)%mi(turnR(indind),2)=j
				turnR(indind)=turnR(indind)+1
				END IF
			END DO
		END IF
	END DO
!Internal degrees of freedom absent
ELSE

	DO gamma=1,SIZE(LabelRight(link+2)%vi)
		IF(LabelRight(link+2)%vi(gamma)<1000) THEN
		DO j=1,localSize
			IF((minNR<=LabelRight(link+2)%vi(gamma)+j-1) &
				.AND.(LabelRight(link+2)%vi(gamma)+j-1<=maxNR)) THEN
				!totQ-number on right
			indind=totQ-(LabelRight(link+2)%vi(gamma)+j-1)-minNL+1
			!indind is block number, turnR(indind) is the position in the block
			indR(indind)%mi(turnR(indind),1)=gamma
			indR(indind)%mi(turnR(indind),2)=j
			turnR(indind)=turnR(indind)+1
			END IF
		END DO
		END IF
	END DO
END IF


END SUBROUTINE IndexRight
	
SUBROUTINE FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
!
!Purpose: Form the NumBlocks Blocks that make up the full theta
!
IMPLICIT NONE
TYPE(matrix), POINTER :: BlockTheta(:)
TYPE(matrixInt), POINTER :: indL(:), indR(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Theta(:,:,:,:)
INTEGER :: chi, NumOfBlocks, k, falp, fgamm, lgt, trv

chi=SIZE(Theta,1)
NumOfBlocks=SIZE(BlockTheta,1)

!$OMP PARALLEL DO PRIVATE(lgt,trv,falp,fgamm) schedule(dynamic)		
DO k=1,NumOfBlocks
	IF((BlockSize(k,2).ne.0).AND.(BlockSize(k,3).ne.0)) THEN
		lgt=BlockSize(k,2)
		trv=BlockSize(k,3)
		DO falp=1,lgt
			DO fgamm=1,trv
				BlockTheta(k)%m(falp,fgamm)= &
				Theta(indL(k)%mi(falp,1),indL(k)%mi(falp,2),indR(k)%mi(fgamm,2),indR(k)%mi(fgamm,1))
			END DO
		END DO
	ELSE
	BlockTheta(k)%m=0.0_rKind
	END IF
END DO
!$OMP END PARALLEL DO

END SUBROUTINE FormBlockTheta


SUBROUTINE SVDNC(US, SS, VS, BlockTheta, BlockSize)
!
!Purpose: Perform an SVD on each one of the Blocks from FormBlockTheta
!
IMPLICIT NONE
TYPE(matrix), POINTER :: US(:), VS(:), BlockTheta(:)
TYPE(vector), POINTER :: SS(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: k, NumOfBlocks

NumOfBlocks=SIZE(BlockTheta,1)

!$OMP PARALLEL DO 	
DO k=1,NumOfBlocks
	!Allocate SVD variables for the block
	CALL ZGESVD_Wrapper(US(k)%m,SS(k)%v,VS(k)%m,BlockTheta(k)%m)
END DO
!$OMP END PARALLEL DO	
END SUBROUTINE SVDNC
	
SUBROUTINE FlattenSS(SS, ssfl, BlockSize)
!
!Purpose: Compile Singular values from all blocks into one vector ssfl
!
IMPLICIT NONE
TYPE(vector), POINTER :: SS(:)
REAL(KIND=rKind), INTENT(OUT) :: ssfl(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: k, fbeta, NumOfBlocks, SizeOfSS, turn

NumOfBlocks=SIZE(BlockSize,1)
SizeOfSS=SIZE(ssfl,1)
turn=1
DO k=1,NumOfBlocks
	!There are min(m,n) singular values of an mXn matrix
	DO fbeta=1,MIN(BlockSize(k,2),BlockSize(k,3))
		ssfl(turn)=SS(k)%v(fbeta)
		turn=turn+1
	END DO
END DO
  
END SUBROUTINE FlattenSS
	
SUBROUTINE Ordering(RealArray, order)
!
!Purpose: Create an array of the indices of the singular values from greatest to least
! Should eventualy be replaced by a more efficient sorting algorithm
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: RealArray(:)
INTEGER, INTENT(OUT) :: order(:)
INTEGER :: SizeOfSS, beta, MaxElement
REAL(KIND=rKind) :: DeformedArray(SIZE(RealArray,1))
INTEGER :: MINLOC_array(1)

order=0
SizeOfSS=SIZE(RealArray,1)
MaxElement=MAXVAL(RealArray)
DeformedArray=RealArray
!Work backwards
DO beta=SizeOfSS,1,-1
	!Find position of minimum value
	MINLOC_array=MINLOC(DeformedArray)
	!Place that in the next available position
	order(beta)=MINLOC_array(1)
	!Remove that element from consideration
	DeformedArray(MINLOC_array(1))=MaxElement+10.0
END DO


END SUBROUTINE Ordering
	
SUBROUTINE JudgePosition(Position, order, BlockSize)
!
!Purpose: Find the "inverse map" to ordering above i.e. find the block index and
! the position within the block for each value in the new singular value ordering
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(OUT) :: Position(:,:)
INTEGER, INTENT(IN) :: order(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: ordertemp(SIZE(order,1))
INTEGER :: k, beta, NumOfBlocks, sz

Position=0
sz=SIZE(order,1)
NumOfBlocks=SIZE(BlockSize,1)
ordertemp=order
DO beta=1,sz
	DO k=1,NumOfBlocks
		IF(ordertemp(beta)<=0) EXIT
		ordertemp(beta)=ordertemp(beta)-MIN(BlockSize(k,2),BlockSize(k,3))
		Position(beta,1)=k
		Position(beta,2)=ordertemp(beta)+MIN(BlockSize(k,2),BlockSize(k,3))
	END DO
END DO
END SUBROUTINE JudgePosition
	
SUBROUTINE FormLambdaNC(Lambda1, truncerr, ssfl, order)
!
!Purpose: Form Lambdas from the ordered singular values
!
IMPLICIT NONE
TYPE(vector) :: Lambda1
REAL(KIND=rKind), INTENT(OUT) :: truncerr
REAL(KIND=rKind), INTENT(IN) :: ssfl(:)
INTEGER, INTENT(IN) :: order(:)
REAL(KIND=rKind), ALLOCATABLE :: tempLam(:)
INTEGER :: chi, beta, sz, counter
REAL(KIND=rKind) ::  norm, norm2

!PRINT *, 'Lambda bf', Lambda1%v
norm2=DOT_PRODUCT(ssfl,ssfl)
norm=0.0_8

counter=0
sz=SIZE(order,1)
ALLOCATE(tempLam(sz))
DO beta=1,sz
	IF((ABS(1.0_rKind-norm/norm2).ge.truncLimit).and.(counter.lt.chiLimit)) THEN
		norm=norm+ssfl(order(beta))*ssfl(order(beta))
		tempLam(beta)=ssfl(order(beta))
		counter=counter+1
	ELSE
		tempLam(beta)=0.0_rKind
	END IF
END DO

DEALLOCATE(Lambda1%v)
ALLOCATE(Lambda1%v(counter))
truncerr = ABS(1.0_rKind - norm/norm2)
DO beta=1,counter
	Lambda1%v(beta)=tempLam(beta)/SQRT(norm)
END DO


DEALLOCATE(tempLam)

END SUBROUTINE FormLambdaNC

SUBROUTINE UpdateLabelLeft(link, chi, LabelLeft, minNL, Position, ssfl, order)
!
!Purpose: Update LabelLeft based on the new signular value ordering
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:)
INTEGER, INTENT(IN) :: link, minNL
INTEGER, INTENT(IN) :: Position(:,:)
REAL(KIND=rKind), INTENT(IN) :: ssfl(:)
INTEGER, INTENT(IN) :: order(:)
INTEGER :: order2(SIZE(order,1))
INTEGER :: chi, beta, sz, l1

sz=SIZE(order,1)
IF(chi.ne.SIZE(LabelLeft(link+1)%vi,1)) THEN
	DEALLOCATE(LabelLeft(link+1)%vi)
	ALLOCATE(LabelLeft(link+1)%vi(chi))
END IF



DO l1=1,sz
	order2(order(l1))=l1
END DO
DO beta=1,chi
	IF(beta<=sz) THEN
		!Truncate Schmidt index at chi
		IF(order2(order(beta))<=chi) THEN
			!Keep only the singular values greater than 10**(-15)
			IF(ssfl(order(beta))>10.0_rKind**(-15)) THEN
				LabelLeft(link+1)%vi(beta)=minNL+Position(beta,1)-1
			ELSE
				!If the state is simply unused, make LabelLeft 10000
				LabelLeft(link+1)%vi(beta)=10000
			END IF
		END IF
	ELSE
		!If the state is simply unused, make LabelLeft 10000
		LabelLeft(link+1)%vi(beta)=10000
	END IF

END DO	
END SUBROUTINE UpdateLabelLeft
	
SUBROUTINE UpdateLabelRight(link, LabelLeft, LabelRight)
!
!Purpose: Update LabelRight based on the new signular value ordering
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: link
INTEGER :: chi, beta

chi=SIZE(LabelLeft(link+1)%vi,1)
IF(chi.ne.SIZE(LabelRight(link+1)%vi,1)) THEN
	DEALLOCATE(LabelRight(link+1)%vi)
	ALLOCATE(LabelRight(link+1)%vi(chi))
END IF
DO beta=1,chi
	IF(LabelLeft(link+1)%vi(beta)<1000) THEN
		LabelRight(link+1)%vi(beta)=totQ-LabelLeft(link+1)%vi(beta)
	ELSE
		!If the state is simply unused, make LabelRight 10000
		LabelRight(link+1)%vi(beta)=10000
	END IF
END DO
END SUBROUTINE UpdateLabelRight


SUBROUTINE FormGamma1NC_t( Gamma1, US, indL, order, Position, BlockSize)
!
!Purpose: Form Left Gamma from the left SVD matrix
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma1(:,:,:)
TYPE(matrix), POINTER :: US(:)
TYPE(matrixInt), POINTER :: indL(:)
INTEGER, INTENT(IN) :: order(:)
INTEGER, INTENT(IN) :: Position(:,:), BlockSize(:,:)
INTEGER :: sz, chi, falp, beta, bn, psn,i,j,k
sz=SIZE(order,1)

Gamma1=CMPLX(0.0,KIND=rKind)
!$OMP PARALLEL DO PRIVATE(bn, psn, falp) schedule(dynamic)
DO beta=1,SIZE(Gamma1,3)
	IF(beta<=sz) THEN
		bn=Position(beta,1)
		psn=Position(beta,2)
		DO falp=1,BlockSize(bn,2)
			Gamma1(indL(bn)%mi(falp,1),indL(bn)%mi(falp,2),beta)= US(bn)%m(falp,psn)
		END DO
	END IF
END DO
!$OMP END PARALLEL DO
END SUBROUTINE FormGamma1NC_t

SUBROUTINE FormGamma1NC_m( Gamma1, US, indL, order, Position, BlockSize)
!
!Purpose: Form Left Gamma from the left SVD matrix
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma1(:,:)
TYPE(matrix), POINTER :: US(:)
TYPE(matrixInt), POINTER :: indL(:)
INTEGER, INTENT(IN) :: order(:)
INTEGER, INTENT(IN) :: Position(:,:), BlockSize(:,:)
INTEGER :: sz, chi, falp, beta, bn, psn, i,j,k
sz=SIZE(order,1)

Gamma1=CMPLX(0.0,KIND=rKind)
!$OMP PARALLEL DO PRIVATE(bn, psn, falp) schedule(dynamic)
DO beta=1,SIZE(Gamma1,2)
	IF(beta<=sz) THEN
		bn=Position(beta,1)
		psn=Position(beta,2)
		DO falp=1,BlockSize(bn,2)
			Gamma1(indL(bn)%mi(falp,1),beta)= US(bn)%m(falp,psn)
		END DO
	END IF
END DO
!$OMP END PARALLEL DO
END SUBROUTINE FormGamma1NC_m

SUBROUTINE FormGamma2NC_t(Gamma2, VS, indR, order, Position, BlockSize)
!
!Purpose: Form right Gamma from the right SVD matrix
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma2(:,:,:)
TYPE(matrix), POINTER :: VS(:)
TYPE(matrixInt), POINTER :: indR(:)
INTEGER, INTENT(IN) :: order(:)
INTEGER, INTENT(IN) :: Position(:,:), BlockSize(:,:)
INTEGER :: sz, chi, beta, fgamm, bn, psn,i,j,k
sz=SIZE(order,1)

Gamma2=CMPLX(0.0,KIND=rKind)
!$OMP PARALLEL DO PRIVATE(bn, psn, fgamm) schedule(dynamic)
DO beta=1,SIZE(Gamma2,1)
	IF(beta<=sz) THEN
		bn=Position(beta,1)
		psn=Position(beta,2)
		DO fgamm=1,BlockSize(bn,3)
			Gamma2(beta,indR(bn)%mi(fgamm,2),indR(bn)%mi(fgamm,1))= VS(bn)%m(psn,fgamm)
		END DO
	END IF
END DO
!$OMP END PARALLEL DO
END SUBROUTINE FormGamma2NC_t

SUBROUTINE FormGamma2NC_m(Gamma2,  VS, indR, order, Position, BlockSize)
!
!Purpose: Form right Gamma from the right SVD matrix
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma2(:,:)
TYPE(matrix), POINTER :: VS(:)
TYPE(matrixInt), POINTER :: indR(:)
INTEGER, INTENT(IN) :: order(:)
INTEGER, INTENT(IN) :: Position(:,:), BlockSize(:,:)
INTEGER :: sz, chi, beta, fgamm, bn, psn,i,j
sz=SIZE(order,1)


Gamma2=CMPLX(0.0,KIND=rKind)
!$OMP PARALLEL DO PRIVATE(bn, psn, fgamm) schedule(dynamic)
DO beta=1,SIZE(Gamma2,1)
	IF(beta<=sz) THEN
		bn=Position(beta,1)
		psn=Position(beta,2)
		DO fgamm=1,BlockSize(bn,3)
			Gamma2(beta,indR(bn)%mi(fgamm,1))= VS(bn)%m(psn,fgamm)
		END DO
	END IF
END DO
!$OMP END PARALLEL DO

END SUBROUTINE FormGamma2NC_m


SUBROUTINE TwoSiteOpNCEdge_t(link, Op2, Gammas, Lambdas, LabelLeft, LabelRight, truncerr)
!
!Purpose: Perform the two-site operation Op2 on the sites divided by link
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
TYPE(vector), POINTER :: SS(:)
TYPE(vector) :: ssflat
TYPE(vectorInt) :: ord
TYPE(matrix), POINTER :: BlockTheta(:), US(:), VS(:)
TYPE(matrixInt), POINTER :: indL(:), indR(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
					SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
INTEGER, DIMENSION(:,:), ALLOCATABLE :: BlockSize, Position
INTEGER :: minNL, maxNL, minNR, maxNR, dl, bd

CALL FormThetaNCEdge(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
CALL ThetaOperationNC(link, Op2, Theta, LabelLeft, LabelRight,1)
CALL minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
ALLOCATE(BlockSize(maxNL-minNL+1,3), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate BlockSize in TwoSiteOpNC'
END IF
CALL SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight,1)
CALL AllocateIndexLR(indL, indR , BlockSize)
CALL IndexLeft(link, indL, minNL, maxNL, LabelLeft,1)
CALL IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight,1)
CALL AllocateBlockTheta(BlockTheta, BlockSize)
CALL FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
CALL AllocateUSV(US, SS, VS, BlockSize)
CALL SVDNC(US, SS, VS, BlockTheta, BlockSize)
CALL AllocateSSflat(ssflat, BlockSize)
CALL FlattenSS(SS, ssflat%v, BlockSize)
ALLOCATE(ord%vi(SIZE(ssflat%v,1)), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate ord in TwoSiteOpNC'
END IF
CALL Ordering(ssflat%v, ord%vi)
ALLOCATE(Position(SIZE(ssflat%v,1),2), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate Position in TwoSiteOpNC'
END IF
CALL JudgePosition(Position, ord%vi, BlockSize)
CALL FormLambdaNC(Lambdas(link+1), truncerr, ssflat%v, ord%vi)
CALL UpdateLabelLeft(link, SIZE(Lambdas(link+1)%v,1),LabelLeft, minNL, Position, ssflat%v, ord%vi)
CALL UpdateLabelRight(link, LabelLeft, LabelRight)
dl=SIZE(Gammas(link)%t,2)
bd=SIZE(Gammas(link)%t,1)
DEALLOCATE(Gammas(link)%t)
ALLOCATE(Gammas(link)%t(bd,dl,SIZE(Lambdas(link+1)%v,1)))
Gammas(link)%t=0.0_rKind
CALL FormGamma1NC( Gammas(link)%t, US, indL, ord%vi, Position, BlockSize)
dl=SIZE(Gammas(link+1)%t,2)
bd=SIZE(Gammas(link+1)%t,3)
DEALLOCATE(Gammas(link+1)%t)
ALLOCATE(Gammas(link+1)%t(SIZE(Lambdas(link+1)%v,1),dl,bd))
Gammas(link+1)%t=0.0_rKind
CALL FormGamma2NC(Gammas(link+1)%t, VS, indR, ord%vi, Position, BlockSize)


DEALLOCATE(BlockSize, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate BlockSize in TwoSiteOpNC'
END IF
CALL DeallocateIndexLR(indL,indR)
CALL DeallocateBlockTheta(BlockTheta)
CALL DeallocateUSV(US, SS, VS)
DEALLOCATE(ssflat%v, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate ssflat in TwoSiteOpNC'
END IF
DEALLOCATE(ord%vi, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate ord in TwoSiteOpNC'
END IF
DEALLOCATE(Position, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate position in TwoSiteOpNC'
END IF

END SUBROUTINE TwoSiteOpNCEdge_t

SUBROUTINE TwoSiteOpNCEdge_m(link, Op2, Gammas, Lambdas, LabelLeft, LabelRight, truncerr)
!
!Purpose: Perform the two-site operation Op2 on the sites divided by link
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
TYPE(vector), POINTER :: SS(:)
TYPE(vector) :: ssflat
TYPE(vectorInt) :: ord
TYPE(matrix), POINTER :: BlockTheta(:), US(:), VS(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%m,1),localSize, &
					localSize,SIZE(Gammas(link+1)%m,2))
TYPE(matrixInt), POINTER :: indL(:), indR(:)
INTEGER, DIMENSION(:,:), ALLOCATABLE :: BlockSize, Position
INTEGER :: minNL, maxNL, minNR, maxNR, bd


	CALL FormThetaNCEdge(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
	CALL ThetaOperationNC(link, Op2, Theta, LabelLeft, LabelRight)
	CALL minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
	ALLOCATE(BlockSize(maxNL-minNL+1,3), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate BlockSize in TwoSiteOpNC'
	END IF
	CALL SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight)
	CALL AllocateIndexLR(indL, indR , BlockSize)
	CALL IndexLeft(link, indL, minNL, maxNL, LabelLeft)
	CALL IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight)
	CALL AllocateBlockTheta(BlockTheta, BlockSize)
	CALL FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
	CALL AllocateUSV(US, SS, VS, BlockSize)
	CALL SVDNC(US, SS, VS, BlockTheta, BlockSize)
	CALL AllocateSSflat(ssflat, BlockSize)
	CALL FlattenSS(SS, ssflat%v, BlockSize)
	ALLOCATE(ord%vi(SIZE(ssflat%v,1)), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate ord in TwoSiteOpNC'
	END IF
	CALL Ordering(ssflat%v, ord%vi)
	ALLOCATE(Position(SIZE(ssflat%v,1),2), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate Position in TwoSiteOpNC'
	END IF
	CALL JudgePosition(Position, ord%vi, BlockSize)
	CALL FormLambdaNC(Lambdas(link+1), truncerr, ssflat%v, ord%vi)
	CALL UpdateLabelLeft(link, SIZE(Lambdas(link+1)%v,1),LabelLeft, minNL, Position, ssflat%v, ord%vi)
	CALL UpdateLabelRight(link, LabelLeft, LabelRight)
	bd=SIZE(Gammas(link)%m,1)
	DEALLOCATE(Gammas(link)%m)
	ALLOCATE(Gammas(link)%m(bd,SIZE(Lambdas(link+1)%v,1)))
	Gammas(link)%m=0.0_rKind
	CALL FormGamma1NC( Gammas(link)%m, US, indL, ord%vi, Position, BlockSize)
	bd=SIZE(Gammas(link+1)%m,2)
	DEALLOCATE(Gammas(link+1)%m)
	ALLOCATE(Gammas(link+1)%m(SIZE(Lambdas(link+1)%v,1),bd))
	Gammas(link+1)%m=0.0_rKind
	CALL FormGamma2NC(Gammas(link+1)%m, VS, indR, ord%vi, Position, BlockSize)

		
DEALLOCATE(BlockSize, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate BlockSize in TwoSiteOpNC'
END IF
CALL DeallocateIndexLR(indL,indR)
CALL DeallocateBlockTheta(BlockTheta)
CALL DeallocateUSV(US, SS, VS)
DEALLOCATE(ssflat%v, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate ssflat in TwoSiteOpNC'
END IF
DEALLOCATE(ord%vi, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate ord in TwoSiteOpNC'
END IF
DEALLOCATE(Position, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate position in TwoSiteOpNC'
END IF

END SUBROUTINE TwoSiteOpNCEdge_m

SUBROUTINE TwoSiteOpNCR_t(link, Op2, Gammas, Lambdas, LabelLeft, LabelRight, truncerr)
!
!Purpose: Perform the two-site operation Op2 on the sites divided by link
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
TYPE(vector), POINTER :: SS(:)
TYPE(vector) :: ssflat
TYPE(vectorInt) :: ord
TYPE(matrix), POINTER :: BlockTheta(:), US(:), VS(:)
TYPE(matrixInt), POINTER :: indL(:), indR(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
					SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
INTEGER, DIMENSION(:,:), ALLOCATABLE :: BlockSize, Position
INTEGER :: minNL, maxNL, minNR, maxNR, dl, bd

CALL FormThetaNCR(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
CALL ThetaOperationNC(link, Op2, Theta, LabelLeft, LabelRight,1)
CALL minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
ALLOCATE(BlockSize(maxNL-minNL+1,3), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate BlockSize in TwoSiteOpNC'
END IF
CALL SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight,1)
CALL AllocateIndexLR(indL, indR , BlockSize)
CALL IndexLeft(link, indL, minNL, maxNL, LabelLeft,1)
CALL IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight,1)
CALL AllocateBlockTheta(BlockTheta, BlockSize)
CALL FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
CALL AllocateUSV(US, SS, VS, BlockSize)
CALL SVDNC(US, SS, VS, BlockTheta, BlockSize)
CALL AllocateSSflat(ssflat, BlockSize)
CALL FlattenSS(SS, ssflat%v, BlockSize)
ALLOCATE(ord%vi(SIZE(ssflat%v,1)), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate ord in TwoSiteOpNC'
END IF
CALL Ordering(ssflat%v, ord%vi)
ALLOCATE(Position(SIZE(ssflat%v,1),2), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate Position in TwoSiteOpNC'
END IF
CALL JudgePosition(Position, ord%vi, BlockSize)
CALL FormLambdaNC(Lambdas(link+1), truncerr, ssflat%v, ord%vi)
CALL UpdateLabelLeft(link, SIZE(Lambdas(link+1)%v,1),LabelLeft, minNL, Position, ssflat%v, ord%vi)
CALL UpdateLabelRight(link, LabelLeft, LabelRight)
dl=SIZE(Gammas(link)%t,2)
bd=SIZE(Gammas(link)%t,1)
DEALLOCATE(Gammas(link)%t)
ALLOCATE(Gammas(link)%t(bd,dl,SIZE(Lambdas(link+1)%v,1)))
Gammas(link)%t=0.0_rKind
CALL FormGamma1NC( Gammas(link)%t, US, indL, ord%vi, Position, BlockSize)
dl=SIZE(Gammas(link+1)%t,2)
bd=SIZE(Gammas(link+1)%t,3)
DEALLOCATE(Gammas(link+1)%t)
ALLOCATE(Gammas(link+1)%t(SIZE(Lambdas(link+1)%v,1),dl,bd))
Gammas(link+1)%t=0.0_rKind
CALL FormGamma2NC(Gammas(link+1)%t, VS, indR, ord%vi, Position, BlockSize)


DEALLOCATE(BlockSize, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate BlockSize in TwoSiteOpNC'
END IF
CALL DeallocateIndexLR(indL,indR)
CALL DeallocateBlockTheta(BlockTheta)
CALL DeallocateUSV(US, SS, VS)
DEALLOCATE(ssflat%v, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate ssflat in TwoSiteOpNC'
END IF
DEALLOCATE(ord%vi, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate ord in TwoSiteOpNC'
END IF
DEALLOCATE(Position, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate position in TwoSiteOpNC'
END IF

END SUBROUTINE TwoSiteOpNCR_t

SUBROUTINE TwoSiteOpNCR_m(link, Op2, Gammas, Lambdas, LabelLeft, LabelRight, truncerr)
!
!Purpose: Perform the two-site operation Op2 on the sites divided by link
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
TYPE(vector), POINTER :: SS(:)
TYPE(vector) :: ssflat
TYPE(vectorInt) :: ord
TYPE(matrix), POINTER :: BlockTheta(:), US(:), VS(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%m,1),localSize, &
					localSize,SIZE(Gammas(link+1)%m,2))
TYPE(matrixInt), POINTER :: indL(:), indR(:)
INTEGER, DIMENSION(:,:), ALLOCATABLE :: BlockSize, Position
INTEGER :: minNL, maxNL, minNR, maxNR, bd


	CALL FormThetaNCR(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
	CALL ThetaOperationNC(link, Op2, Theta, LabelLeft, LabelRight)
	CALL minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
	ALLOCATE(BlockSize(maxNL-minNL+1,3), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate BlockSize in TwoSiteOpNC'
	END IF
	CALL SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight)
	CALL AllocateIndexLR(indL, indR , BlockSize)
	CALL IndexLeft(link, indL, minNL, maxNL, LabelLeft)
	CALL IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight)
	CALL AllocateBlockTheta(BlockTheta, BlockSize)
	CALL FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
	CALL AllocateUSV(US, SS, VS, BlockSize)
	CALL SVDNC(US, SS, VS, BlockTheta, BlockSize)
	CALL AllocateSSflat(ssflat, BlockSize)
	CALL FlattenSS(SS, ssflat%v, BlockSize)
	ALLOCATE(ord%vi(SIZE(ssflat%v,1)), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate ord in TwoSiteOpNC'
	END IF
	CALL Ordering(ssflat%v, ord%vi)
	ALLOCATE(Position(SIZE(ssflat%v,1),2), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate Position in TwoSiteOpNC'
	END IF
	CALL JudgePosition(Position, ord%vi, BlockSize)
	CALL FormLambdaNC(Lambdas(link+1), truncerr, ssflat%v, ord%vi)
	CALL UpdateLabelLeft(link, SIZE(Lambdas(link+1)%v,1),LabelLeft, minNL, Position, ssflat%v, ord%vi)
	CALL UpdateLabelRight(link, LabelLeft, LabelRight)
	bd=SIZE(Gammas(link)%m,1)
	DEALLOCATE(Gammas(link)%m)
	ALLOCATE(Gammas(link)%m(bd,SIZE(Lambdas(link+1)%v,1)))
	Gammas(link)%m=0.0_rKind
	CALL FormGamma1NC( Gammas(link)%m, US, indL, ord%vi, Position, BlockSize)
	bd=SIZE(Gammas(link+1)%m,2)
	DEALLOCATE(Gammas(link+1)%m)
	ALLOCATE(Gammas(link+1)%m(SIZE(Lambdas(link+1)%v,1),bd))
	Gammas(link+1)%m=0.0_rKind
	CALL FormGamma2NC(Gammas(link+1)%m, VS, indR, ord%vi, Position, BlockSize)

		
DEALLOCATE(BlockSize, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate BlockSize in TwoSiteOpNC'
END IF
CALL DeallocateIndexLR(indL,indR)
CALL DeallocateBlockTheta(BlockTheta)
CALL DeallocateUSV(US, SS, VS)
DEALLOCATE(ssflat%v, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate ssflat in TwoSiteOpNC'
END IF
DEALLOCATE(ord%vi, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate ord in TwoSiteOpNC'
END IF
DEALLOCATE(Position, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate position in TwoSiteOpNC'
END IF

END SUBROUTINE TwoSiteOpNCR_m

SUBROUTINE TwoSiteOpNCL_t(link, Op2, Gammas, Lambdas, LabelLeft, LabelRight, truncerr)
!
!Purpose: Perform the two-site operation Op2 on the sites divided by link
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
TYPE(vector), POINTER :: SS(:)
TYPE(vector) :: ssflat
TYPE(vectorInt) :: ord
TYPE(matrix), POINTER :: BlockTheta(:), US(:), VS(:)
TYPE(matrixInt), POINTER :: indL(:), indR(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
					SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
INTEGER, DIMENSION(:,:), ALLOCATABLE :: BlockSize, Position
INTEGER :: minNL, maxNL, minNR, maxNR, dl, bd

CALL FormThetaNCL(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
CALL ThetaOperationNC(link, Op2, Theta, LabelLeft, LabelRight,1)
CALL minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
ALLOCATE(BlockSize(maxNL-minNL+1,3), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate BlockSize in TwoSiteOpNC'
END IF
CALL SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight,1)
CALL AllocateIndexLR(indL, indR , BlockSize)
CALL IndexLeft(link, indL, minNL, maxNL, LabelLeft,1)
CALL IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight,1)
CALL AllocateBlockTheta(BlockTheta, BlockSize)
CALL FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
CALL AllocateUSV(US, SS, VS, BlockSize)
CALL SVDNC(US, SS, VS, BlockTheta, BlockSize)
CALL AllocateSSflat(ssflat, BlockSize)
CALL FlattenSS(SS, ssflat%v, BlockSize)
ALLOCATE(ord%vi(SIZE(ssflat%v,1)), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate ord in TwoSiteOpNC'
END IF
CALL Ordering(ssflat%v, ord%vi)
ALLOCATE(Position(SIZE(ssflat%v,1),2), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate Position in TwoSiteOpNC'
END IF
CALL JudgePosition(Position, ord%vi, BlockSize)
CALL FormLambdaNC(Lambdas(link+1), truncerr, ssflat%v, ord%vi)
CALL UpdateLabelLeft(link, SIZE(Lambdas(link+1)%v,1),LabelLeft, minNL, Position, ssflat%v, ord%vi)
CALL UpdateLabelRight(link, LabelLeft, LabelRight)
dl=SIZE(Gammas(link)%t,2)
bd=SIZE(Gammas(link)%t,1)
DEALLOCATE(Gammas(link)%t)
ALLOCATE(Gammas(link)%t(bd,dl,SIZE(Lambdas(link+1)%v,1)))
Gammas(link)%t=0.0_rKind
CALL FormGamma1NC( Gammas(link)%t, US, indL, ord%vi, Position, BlockSize)
dl=SIZE(Gammas(link+1)%t,2)
bd=SIZE(Gammas(link+1)%t,3)
DEALLOCATE(Gammas(link+1)%t)
ALLOCATE(Gammas(link+1)%t(SIZE(Lambdas(link+1)%v,1),dl,bd))
Gammas(link+1)%t=0.0_rKind
CALL FormGamma2NC(Gammas(link+1)%t, VS, indR, ord%vi, Position, BlockSize)


DEALLOCATE(BlockSize, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate BlockSize in TwoSiteOpNC'
END IF
CALL DeallocateIndexLR(indL,indR)
CALL DeallocateBlockTheta(BlockTheta)
CALL DeallocateUSV(US, SS, VS)
DEALLOCATE(ssflat%v, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate ssflat in TwoSiteOpNC'
END IF
DEALLOCATE(ord%vi, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate ord in TwoSiteOpNC'
END IF
DEALLOCATE(Position, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate position in TwoSiteOpNC'
END IF

END SUBROUTINE TwoSiteOpNCL_t

SUBROUTINE TwoSiteOpNCL_m(link, Op2, Gammas, Lambdas, LabelLeft, LabelRight, truncerr)
!
!Purpose: Perform the two-site operation Op2 on the sites divided by link
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(matrix), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
TYPE(vector), POINTER :: SS(:)
TYPE(vector) :: ssflat
TYPE(vectorInt) :: ord
TYPE(matrix), POINTER :: BlockTheta(:), US(:), VS(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%m,1),localSize, &
					localSize,SIZE(Gammas(link+1)%m,2))
TYPE(matrixInt), POINTER :: indL(:), indR(:)
INTEGER, DIMENSION(:,:), ALLOCATABLE :: BlockSize, Position
INTEGER :: minNL, maxNL, minNR, maxNR, bd


	CALL FormThetaNCL(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
	CALL ThetaOperationNC(link, Op2, Theta, LabelLeft, LabelRight)
	CALL minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
	ALLOCATE(BlockSize(maxNL-minNL+1,3), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate BlockSize in TwoSiteOpNC'
	END IF
	CALL SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight)
	CALL AllocateIndexLR(indL, indR , BlockSize)
	CALL IndexLeft(link, indL, minNL, maxNL, LabelLeft)
	CALL IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight)
	CALL AllocateBlockTheta(BlockTheta, BlockSize)
	CALL FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
	CALL AllocateUSV(US, SS, VS, BlockSize)
	CALL SVDNC(US, SS, VS, BlockTheta, BlockSize)
	CALL AllocateSSflat(ssflat, BlockSize)
	CALL FlattenSS(SS, ssflat%v, BlockSize)
	ALLOCATE(ord%vi(SIZE(ssflat%v,1)), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate ord in TwoSiteOpNC'
	END IF
	CALL Ordering(ssflat%v, ord%vi)
	ALLOCATE(Position(SIZE(ssflat%v,1),2), STAT=statInt)
	IF(statInt.ne.0) THEN
		PRINT *, 'Failed to allocate Position in TwoSiteOpNC'
	END IF
	CALL JudgePosition(Position, ord%vi, BlockSize)
	CALL FormLambdaNC(Lambdas(link+1), truncerr, ssflat%v, ord%vi)
	CALL UpdateLabelLeft(link, SIZE(Lambdas(link+1)%v,1),LabelLeft, minNL, Position, ssflat%v, ord%vi)
	CALL UpdateLabelRight(link, LabelLeft, LabelRight)
	bd=SIZE(Gammas(link)%m,1)
	DEALLOCATE(Gammas(link)%m)
	ALLOCATE(Gammas(link)%m(bd,SIZE(Lambdas(link+1)%v,1)))
	Gammas(link)%m=0.0_rKind
	CALL FormGamma1NC( Gammas(link)%m, US, indL, ord%vi, Position, BlockSize)
	bd=SIZE(Gammas(link+1)%m,2)
	DEALLOCATE(Gammas(link+1)%m)
	ALLOCATE(Gammas(link+1)%m(SIZE(Lambdas(link+1)%v,1),bd))
	Gammas(link+1)%m=0.0_rKind
	CALL FormGamma2NC(Gammas(link+1)%m, VS, indR, ord%vi, Position, BlockSize)

		
DEALLOCATE(BlockSize, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate BlockSize in TwoSiteOpNC'
END IF
CALL DeallocateIndexLR(indL,indR)
CALL DeallocateBlockTheta(BlockTheta)
CALL DeallocateUSV(US, SS, VS)
DEALLOCATE(ssflat%v, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate ssflat in TwoSiteOpNC'
END IF
DEALLOCATE(ord%vi, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate ord in TwoSiteOpNC'
END IF
DEALLOCATE(Position, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate position in TwoSiteOpNC'
END IF

END SUBROUTINE TwoSiteOpNCL_m


END MODULE LocalOps
