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
MODULE LinearOps
!
! Purpose: Module Containing derived types/allocations/matrix manipulations
! for OpenSourceTEBD v3.0
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!       8/18/10  M. L. Wall	alpha release
!
USE GlobalData
IMPLICIT NONE

! *** DERIVED TYPES ***

TYPE vector
	REAL(KIND=rKind), POINTER :: v(:)
END TYPE vector

TYPE vectorComplex
	COMPLEX(KIND=rKind), POINTER :: vc(:)
END TYPE vectorComplex

TYPE Charvec
	CHARACTER(len=10), POINTER :: v(:)
END TYPE Charvec

TYPE matrix
	COMPLEX(KIND=rKind), POINTER :: m(:,:)
END TYPE matrix

TYPE tensor
	COMPLEX(KIND=rKind), POINTER :: t(:,:,:)
END TYPE tensor
	
TYPE vectorInt
	INTEGER, POINTER :: vi(:)
END TYPE vectorInt
	
TYPE matrixInt
	INTEGER, POINTER :: mi(:,:)
END TYPE matrixInt
	
TYPE fourtensor
	COMPLEX(KIND=rKind), POINTER :: ft(:,:,:,:)
END TYPE fourtensor

TYPE DecomposedMPO
	TYPE(matrix), POINTER :: l(:)
	TYPE(matrix), POINTER :: r(:)
END TYPE DecomposedMPO

!*** INTERFACES
INTERFACE GEMM
	MODULE PROCEDURE zgemm_wrapper, dgemm_wrapper
END INTERFACE GEMM

!Matrix exponential of real/complex matrices
INTERFACE matrix_exponential
MODULE PROCEDURE matrix_exponential_r,&
				 matrix_exponential_c
END INTERFACE  matrix_exponential

INTERFACE tensorProd
MODULE PROCEDURE tensorProd_r,&
				 tensorProd_c,&
				 tensorProd_cv,&
				 tensorProd_rc,&
				 tensorProd_cr
END INTERFACE  tensorProd

!Trace of A*B AB=rr/cc/rc/cr
INTERFACE TraceMatmul
MODULE PROCEDURE  TraceMatmul_rf,&
				 TraceMatmul_cf,&
				 TraceMatmul_rcf,&
				 TraceMatmul_crf
END INTERFACE  TraceMatmul

INTERFACE ZGESVD_Wrapper
	MODULE PROCEDURE ZGESVD_Wrapperu, ZGESVD_Wrappera
END INTERFACE ZGESVD_Wrapper

!kronecker delta of real/complex vectors
INTERFACE kronDelta
MODULE PROCEDURE kronDelta_r,&
				 kronDelta_c
END INTERFACE  kronDelta

CONTAINS

!!!!!!!!!!!!!!!BEGIN CONTENTS OF INTERFACE matrix_exponential!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE zgemm_wrapper(C,A,B,alphaIn,betaIn)
!
!Purpose: Multiply A and B together and store in C.  alphaIn is an optional scaling parameter
! and betaIn specifies that beta*C is to be added to A*B
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: C(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: A(:,:), B(:,:)
COMPLEX(KIND=rKind), INTENT(IN), OPTIONAL :: alphaIn, betaIn
COMPLEX(KIND=rKind) :: alpha, beta

IF(PRESENT(alphaIn)) THEN
	alpha=alphaIn
ELSE
	alpha=1.0_rKind
END IF

IF(PRESENT(betaIn)) THEN
	beta=betaIn
ELSE
	beta=0.0_rKind
END IF

CALL ZGEMM('N','N',SIZE(A,1),SIZE(B,2),SIZE(A,2),alpha, A,SIZE(A,1),B,SIZE(B,1), beta,C,SIZE(A,1))


END SUBROUTINE zgemm_wrapper

SUBROUTINE dgemm_wrapper(C,A,B,alphaIn,betaIn)
!
!Purpose: Multiply A and B together and store in C.  alphaIn is an optional scaling parameter
! and betaIn specifies that beta*C is to be added to A*B
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(INOUT) :: C(:,:)
REAL(KIND=rKind), INTENT(IN) :: A(:,:), B(:,:)
REAL(KIND=rKind), INTENT(IN), OPTIONAL :: alphaIn, betaIn
REAL(KIND=rKind) :: alpha, beta

IF(PRESENT(alphaIn)) THEN
	alpha=alphaIn
ELSE
	alpha=1.0_rKind
END IF

IF(PRESENT(betaIn)) THEN
	beta=betaIn
ELSE
	beta=0.0_rKind
END IF

CALL DGEMM('N','N',SIZE(A,1),SIZE(B,2),SIZE(A,2),alpha, A,SIZE(A,1),B,SIZE(B,1), beta,C,SIZE(A,1))


END SUBROUTINE dgemm_wrapper

SUBROUTINE ZGESVD_Wrapperu(U,S,V,A)
!
!Purpose: Perform an economy SVD on A using ZGESVD
!    
IMPLICIT NONE     
COMPLEX(KIND=rKind), INTENT(INOUT) :: A(:,:), U(:,:), V(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: S(:)
COMPLEX(KIND=rKind) :: work(5*MAX(SIZE(A,1),SIZE(A,2)))
REAL(KIND=rKind) :: rwork(5*MIN(SIZE(A,1),SIZE(A,2)))
INTEGER :: info


CALL ZGESVD('S','A',SIZE(A,1), SIZE(A,2), A, SIZE(A,1),S, U, SIZE(A,1), V, SIZE(A,2),&
		work, 5*MAX(SIZE(A,1),SIZE(A,2)), rwork, info)

END SUBROUTINE ZGESVD_Wrapperu

SUBROUTINE ZGESVD_Wrappera(U,S,V,A)
!
!Purpose: Perform an SVD on A using ZGESVD (return all columns of U and V)
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: A(:,:)       
TYPE(vector) :: S
TYPE(matrix) :: U, V
COMPLEX(KIND=rKind), ALLOCATABLE :: work(:)
REAL(KIND=rKind), ALLOCATABLE :: rwork(:)
INTEGER :: info
INTEGER :: M, N

M=SIZE(A,1)
N=SIZE(A,2)

ALLOCATE(U%m(M,M), V%m(N,N), S%v(MINVAL((/M,N/))))

ALLOCATE(work(5*MAXVAL((/M,N/))))
ALLOCATE(rwork(5*MINVAL((/M,N/))))

!Call the LAPACK routine ZGESVD, which performs a SVD on a general matrix
	CALL ZGESVD('A', 'A', M, N, A, M, S%v, U%m, & 
				M, V%m, N, work, 5*MAXVAL((/M,N/)), rwork, info)

DEALLOCATE(work,rwork)


END SUBROUTINE ZGESVD_Wrappera


SUBROUTINE EconomyZGESVD_Wrapper(U,S,V,A)
!
!Purpose: Perform an economy SVD, allocate U, S, and V to their proper sizes
!
IMPLICIT NONE
COMPLEX(KIND=8), INTENT(INOUT) :: A(:,:)       
TYPE(vector) :: S
TYPE(matrix) :: U, V
COMPLEX(KIND=8), ALLOCATABLE :: work(:)
REAL(KIND=8), ALLOCATABLE :: rwork(:)
INTEGER :: info
INTEGER :: M, N

M=SIZE(A,1)
N=SIZE(A,2)

ALLOCATE(U%m(M,MINVAL((/M,N/))), V%m(MINVAL((/M,N/)),N), S%v(MINVAL((/M,N/))))

ALLOCATE(work(5*MAXVAL((/M,N/))))
ALLOCATE(rwork(5*MINVAL((/M,N/))))

!Call the LAPACK routine ZGESVD, which performs a SVD on a general matrix
	CALL ZGESVD('S', 'S', M, N, A, M, S%v, U%m, & 
				M, V%m, MIN(M,N), work, 5*MAXVAL((/M,N/)), rwork, info)

DEALLOCATE(work,rwork)


END SUBROUTINE EconomyZGESVD_Wrapper

SUBROUTINE svd2(U,S,V,A)
!
!Purpose: Perform an SVD on A or A^T, whichever is more efficient
!    
IMPLICIT NONE     
COMPLEX(KIND=rKind), INTENT(INOUT) :: A(:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: B(:,:)
TYPE(matrix) :: U, V, VT
TYPE(vector) :: S
INTEGER :: m, n

m=SIZE(A,1)
n=SIZE(A,2)

IF(m.ge.n) THEN
ALLOCATE(U%m(m,n),S%v(n),V%m(n,n))
	CALL ZGESVD_Wrapper(U%m,S%v,V%m,A)
ELSE
ALLOCATE(U%m(m,m),S%v(m),VT%m(n,m))
ALLOCATE(B(n,m))
B=TRANSPOSE(A)
	CALL ZGESVD_Wrapper(VT%m,S%v,U%m, B)
	DEALLOCATE(B)
	ALLOCATE(V%m(m,n))
	V%m=TRANSPOSE(VT%m)
	U%m=TRANSPOSE(U%m)
DEALLOCATE(VT%m)
END IF	

END SUBROUTINE svd2


SUBROUTINE matrix_exponential_r(A,Exp_A,tau,n)
!
!Purpose: If matrix_exponential is called with argument types A=real, Exp_A=real, tau=Real
!         Then compute Exp_A=EXP(-tau*A)
!
!Based on routines by schneider, b. i.(nsf)
!
IMPLICIT NONE
INTEGER                                :: n, i, j, k, info
REAL(KIND=rKind), DIMENSION(:,:)                 :: A
REAL(KIND=rKind), DIMENSION(:,:)                 :: Exp_A
REAL(KIND=rKind), DIMENSION(:),   ALLOCATABLE    :: Eigen_Values
REAL(KIND=rKind), DIMENSION(:,:), ALLOCATABLE    :: Eigen_Vectors
REAL(KIND=rKind), DIMENSION(:,:), ALLOCATABLE    :: Temp
REAL(KIND=rKind), DIMENSION(:),   ALLOCATABLE    :: Scratch
REAL(KIND=rKind)                                 :: tau
REAL(KIND=rKind)				:: expeig
CHARACTER (LEN=80)                     :: title
! Allocate storage for diagonalization routine.
ALLOCATE ( Eigen_Values(n), Eigen_Vectors(n,n), Scratch(10*n), Temp(n,n) , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate matrix_exp variables'
			END IF 
Eigen_Vectors = A
!Call LAPACK routine to diagonalize double precision real symmetric matrix
CALL DSYEV('v','l',n,Eigen_Vectors,n,Eigen_Values,           &
              Scratch,10*n,info)
			  
! Form the matrix with exponentials of the eigenvalues on the diagonal
! Then similarity transform back into the original basis
DO i=1,n
  	DO j=1,n
	Exp_A(i,j)=0.0_rKind
		DO k=1,n
  	  	expeig=exp(-tau*Eigen_Values(k))
        Exp_A(i,j) = Exp_A(i,j) + Eigen_Vectors(i,k)*expeig*Eigen_Vectors(j,k)
		END DO
    END DO
END DO

! Deallocate the unneeded storage
DEALLOCATE ( Eigen_Values, Eigen_Vectors, Scratch, Temp  , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate matrix_exp variables'
			END IF   
END SUBROUTINE matrix_exponential_r


SUBROUTINE matrix_exponential_c(A,Exp_A,t,n)
!
!Purpose: If matrix_exponential is called with argument types A=complex, Exp_A=complex, t=complex
!         Then compute Exp_A=EXP(-i*t*A)
!
!Based on routines by schneider, b. i.(nsf)
!
IMPLICIT NONE
INTEGER                                    :: n, i, j, k, info
COMPLEX(KIND=rKind), DIMENSION(:,:)                 :: A
COMPLEX(KIND=rKind), DIMENSION(:,:)                 :: Exp_A
COMPLEX(KIND=rKind), DIMENSION(:,:), ALLOCATABLE    :: Eigen_Vectors
REAL(KIND=rKind),     DIMENSION(:),   ALLOCATABLE    :: Eigen_Values
COMPLEX(KIND=rKind), DIMENSION(:),   ALLOCATABLE    :: Workv
COMPLEX(KIND=rKind), DIMENSION(:,:), ALLOCATABLE    :: Temp
REAL(KIND=rKind),     DIMENSION(:),   ALLOCATABLE    :: Rworkv
COMPLEX(KIND=rKind)                                 :: t
CHARACTER (LEN=80)                         :: title
COMPLEX(KIND=rKind)                                 :: eye=(0.0_rKind,1.0_rKind)
COMPLEX(KIND=rKind)								 :: expeig
! Allocate some storage for diagonalization routine.
ALLOCATE ( Eigen_Values(n), Eigen_Vectors(n,n), Workv(10*n), Rworkv(10*n), Temp(n,n)  , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate matrix_exp variables'
			END IF   
Eigen_Vectors = A
!Call LAPACK routine to diagonalize double precision hermitian matrix
CALL ZHEEV('v','l',n,Eigen_Vectors,n,Eigen_Values,              &
              Workv,10*n,Rworkv,info)

! Form the matrix with exponentials of the eigenvalues on the diagonal
! Then similarity transform back into the original basis
DO i=1,n
	DO j=1,n
	Exp_A(i,j)=CMPLX(0.0,KIND=rKind)
		DO k=1,n
  	  	expeig=exp(-eye*t*Eigen_Values(k))
        Exp_A(i,j) = Exp_A(i,j) + Eigen_Vectors(i,k)*expeig*CONJG(Eigen_Vectors(j,k))
		END DO
    END DO
END DO
  ! Deallocate the unneeded storage
DEALLOCATE ( Eigen_Values, Eigen_Vectors, Workv, Rworkv, Temp  , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate matrix_exp variables'
			END IF     
END SUBROUTINE matrix_exponential_c


!!!!!!!!!!!!!!!END CONTENTS OF INTERFACE matrix_exponential!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!BEGIN CONTENTS OF INTERFACE tensorProd!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION tensorProd_r(A,B)
!
!Purpose: Return the tensor product of real matrices A and B
!    
IMPLICIT NONE     
REAL(KIND=8), INTENT(IN) :: A(:,:), B(:,:)
REAL(KIND=8) :: tensorProd_r(SIZE(A,1)*SIZE(B,1),SIZE(A,2)*SIZE(B,2))
INTEGER i,j,k,l, dA1, dA2, dB1, dB2
dA1 = SIZE(A,1)
dA2 = SIZE(A,2)
dB1 = SIZE(B,1)
dB2 = SIZE(B,2)
	DO i=1,dA1
		DO j=1,dB1
			DO k=1,dA2
				DO l=1,dB2
					tensorProd_r((i-1)*dB1+j,(k-1)*dB2+l)=A(i,k)*B(j,l)
				END DO
			END DO
		END DO
	END DO
END FUNCTION tensorProd_r

FUNCTION tensorProd_c(A,B)
!
!Purpose: Return the tensor product of complex matrices A and B
!  
IMPLICIT NONE
COMPLEX(KIND=8), INTENT(IN) :: A(:,:), B(:,:)
COMPLEX(KIND=8) :: tensorProd_c(SIZE(A,1)*SIZE(B,1),SIZE(A,2)*SIZE(B,2))
INTEGER i,j,k,l, dA1, dA2, dB1, dB2
dA1 = SIZE(A,1)
dA2 = SIZE(A,2)
dB1 = SIZE(B,1)
dB2 = SIZE(B,2)
	DO i=1,dA1
		DO j=1,dB1
			DO k=1,dA2
				DO l=1,dB2
					tensorProd_c((i-1)*dB1+j,(k-1)*dB2+l)=A(i,k)*B(j,l)
				END DO
			END DO
		END DO
	END DO
END FUNCTION tensorProd_c

FUNCTION tensorProd_cv(A,B)
!
!Purpose: Return the tensor product of complex matrices A and B
!  
IMPLICIT NONE
COMPLEX(KIND=8), INTENT(IN) :: A(:), B(:)
COMPLEX(KIND=8) :: tensorProd_cv(SIZE(A),SIZE(B))
INTEGER i,j, dA1, dB1
dA1 = SIZE(A)
dB1 = SIZE(B)
	DO i=1,dA1
		DO j=1,dB1
			tensorProd_cv(i,j)=A(i)*B(j)
		END DO
	END DO
END FUNCTION tensorProd_cv

FUNCTION tensorProd_rc(A,B)
!
!Purpose: Return the tensor product of real matrix A and complex matrix B
!  
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: A(:,:)
COMPLEX(KIND=8), INTENT(IN) :: B(:,:)
COMPLEX(KIND=8) :: tensorProd_rc(SIZE(A,1)*SIZE(B,1),SIZE(A,2)*SIZE(B,2))
INTEGER i,j,k,l, dA1, dA2, dB1, dB2
dA1 = SIZE(A,1)
dA2 = SIZE(A,2)
dB1 = SIZE(B,1)
dB2 = SIZE(B,2)
	DO i=1,dA1
		DO j=1,dB1
			DO k=1,dA2
				DO l=1,dB2
					tensorProd_rc((i-1)*dB1+j,(k-1)*dB2+l)=A(i,k)*B(j,l)
				END DO
			END DO
		END DO
	END DO
END FUNCTION tensorProd_rc

FUNCTION tensorProd_cr(A,B)
!
!Purpose: Return the tensor product of complex matrix A and real matrix B
!  
IMPLICIT NONE
COMPLEX(KIND=8), INTENT(IN) :: A(:,:)
REAL(KIND=8), INTENT(IN) :: B(:,:)
COMPLEX(KIND=8) :: tensorProd_cr(SIZE(A,1)*SIZE(B,1),SIZE(A,2)*SIZE(B,2))
INTEGER i,j,k,l, dA1, dA2, dB1, dB2
dA1 = SIZE(A,1)
dA2 = SIZE(A,2)
dB1 = SIZE(B,1)
dB2 = SIZE(B,2)
	DO i=1,dA1
		DO j=1,dB1
			DO k=1,dA2
				DO l=1,dB2
					tensorProd_cr((i-1)*dB1+j,(k-1)*dB2+l)=A(i,k)*B(j,l)
				END DO
			END DO
		END DO
	END DO
END FUNCTION tensorProd_cr

!!!!!!!!!!!!!!!END CONTENTS OF INTERFACE tensorProd!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!! BEGIN CONTENTS OF INTERFACE kronDelta !!!!!!!!!!!!!!!!!!!!!!!!

INTEGER FUNCTION kronDelta_r(vec1, vec2, dim)
!
!Purpose: Kronecker delta function defined for two real vectors vec1 and vec2.
!
IMPLICIT NONE	
REAL(KIND=rKind), INTENT(IN) :: vec1(:), vec2(:)
INTEGER, INTENT(IN) :: dim
INTEGER :: dim1, dim2, i, j
INTEGER :: booles
	dim1 = SIZE(vec1)
	dim2 = SIZE(vec2)
		IF (dim1 /= dim .OR. dim2 /= dim) THEN
			STOP "Dimensions of input vectors in function kronDelta must be the same."
		END IF
		DO i = 1, dim
			IF (vec1(i) == vec2(i)) THEN
				booles = 1
			ELSE
				booles = 0
				EXIT
			END IF
		END DO
	kronDelta_r = booles
END FUNCTION kronDelta_r	

INTEGER FUNCTION kronDelta_c(vec1, vec2, dim)
!
!Purpose: Kronecker delta function defined for two complex vectors vec1 and vec2.
!
IMPLICIT NONE	
COMPLEX(KIND=rKind), INTENT(IN) :: vec1(:), vec2(:)
INTEGER, INTENT(IN) :: dim
INTEGER :: dim1, dim2, i, j
INTEGER :: booles
	dim1 = SIZE(vec1)
	dim2 = SIZE(vec2)
		IF (dim1 /= dim .OR. dim2 /= dim) THEN
			STOP "Dimensions of input vectors in function kronDelta must be the same."
		END IF
		DO i = 1, dim
			IF (vec1(i) == vec2(i)) THEN
				booles = 1
			ELSE
				booles = 0
				EXIT
			END IF
		END DO
	kronDelta_c = booles
END FUNCTION kronDelta_c	
!!!!!!!!!!!!!!!!!!! END CONTENTS OF INTERFACE kronDelta !!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!! BEGIN CONTENTS OF INTERFACE TraceMatmul !!!!!!!!!!!!!!!!!!!!!!!!

REAL(KIND=8) FUNCTION TraceMatmul_rf(A,B)
!
!Purpose: Function to calculate the trace of real matrices A*B
!
!See manual for more detail
!
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: A(:,:), B(:,:)
REAL(KIND=8), ALLOCATABLE :: C(:,:)
INTEGER i

ALLOCATE(C(SIZE(A,1),SIZE(B,2)))
C=0.0_8
CALL GEMM(C,A,B)
TraceMatmul_rf=0.0_8
DO i=1,SIZE(A,1)
	TraceMatmul_rf=TraceMatmul_rf+C(i,i)
END DO
DEALLOCATE(C)

END FUNCTION TraceMatmul_rf

COMPLEX(KIND=8) FUNCTION TraceMatmul_cf(A,B)
!
!Purpose: Function to calculate the trace of complex matrices A*B
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=8), INTENT(IN) :: A(:,:), B(:,:)
COMPLEX(KIND=8), ALLOCATABLE :: C(:,:)
INTEGER i

ALLOCATE(C(SIZE(A,1),SIZE(B,2)))
C=0.0_8
CALL GEMM(C,A,B)
TraceMatmul_cf=0.0_8
DO i=1,SIZE(A,1)
	TraceMatmul_cf=TraceMatmul_cf+C(i,i)
END DO
DEALLOCATE(C)

END FUNCTION TraceMatmul_cf

COMPLEX(KIND=8) FUNCTION TraceMatmul_rcf(A,B)
!
!Purpose: Function to calculate the trace of real/complex matrices A*B
!
!See manual for more detail
!
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: A(:,:)
COMPLEX(KIND=8), INTENT(IN) :: B(:,:)
COMPLEX(KIND=8), ALLOCATABLE :: C(:,:)
INTEGER i

ALLOCATE(C(SIZE(A,1),SIZE(B,2)))
C=0.0_8
CALL GEMM(C,CMPLX(1.0_8,0.0_8)*A,B)
TraceMatmul_rcf=0.0_8
DO i=1,SIZE(A,1)
	TraceMatmul_rcf=TraceMatmul_rcf+C(i,i)
END DO
DEALLOCATE(C)

END FUNCTION TraceMatmul_rcf

COMPLEX(KIND=8) FUNCTION TraceMatmul_crf(A,B)
!
!Purpose: Function to calculate the trace of complex/real matrices A*B
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=8), INTENT(IN) :: A(:,:)
REAL(KIND=8), INTENT(IN) :: B(:,:)
COMPLEX(KIND=8), ALLOCATABLE :: C(:,:)
INTEGER i

ALLOCATE(C(SIZE(A,1),SIZE(B,2)))
C=0.0_8
CALL GEMM(C,A,CMPLX(1.0_8,0.0_8)*B)
TraceMatmul_crf=0.0_8
DO i=1,SIZE(A,1)
	TraceMatmul_crf=TraceMatmul_crf+C(i,i)
END DO
DEALLOCATE(C)

END FUNCTION TraceMatmul_crf

!!!!!!!!!!!!!!!!!!! END CONTENTS OF INTERFACE TraceMatmul !!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Seed_Init(rank)
!
!Purpose:Randomly seed the random number generator using the intrinsic system clock
! OPTIONAL argument rank takes in the rank of a given processor to ensure no correlation
! of the random numbers on different processors
!
IMPLICIT NONE
INTEGER,INTENT(IN), OPTIONAL :: rank
INTEGER :: seedsize !Size of seed-returned by intrinsic procedure
INTEGER :: clock !Clock time
INTEGER, DIMENSION(:), ALLOCATABLE :: seed !Seed passed to intrinsic procedure
INTEGER :: i !dummy integer

CALL RANDOM_SEED(size=seedsize) !Ask for the size of the seed
ALLOCATE(seed(seedsize)) !Allocate the seed using the given seed size

CALL SYSTEM_CLOCK(COUNT=clock) !Get the clock time

IF(PRESENT(rank)) THEN
	seed=clock+(rank+1)*17*(/(i-1,i=1,seedsize)/) !Generate the seed using the clock and rank
ELSE
	seed=clock+17*(/(i-1,i=1,seedsize)/) !Generate the seed using the clock
END IF

CALL RANDOM_SEED(PUT=seed) !Seed the random number generator using the given seed
END SUBROUTINE Seed_Init

REAL(KIND=rKind) FUNCTION Rand_Num(min,max)
!
!Purpose:Generate a random number in the range [min,max)
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: min, max !Bounds
REAL(KIND=rKind) :: dumrand !Dummy real

CALL RANDOM_NUMBER(dumrand) !Get a random number in the range [0,1)

Rand_Num=(max-min)*dumrand+min !Put the number in the proper range
END FUNCTION Rand_Num

REAL(KIND=rKind) FUNCTION Rand_Numn()
!
!Purpose:Generate a normally distributed random number using the Box-Muller transformation
!
IMPLICIT NONE

REAL(KIND=rKind) :: dumrand, dumrandp

CALL RANDOM_NUMBER(dumrand) !Get a random number in the range [0,1)
CALL RANDOM_NUMBER(dumrandp) !Get another random number in the range [0,1)
!Box-Mueller transform the two uniform random numbers (pick only one of the transformed pair)
Rand_Numn=Sqrt(-2.0_rKind*Log(dumrand))*cos(2.0_rKind*3.1415926535897_rKind*dumrandp)

END  FUNCTION Rand_Numn

INTEGER FUNCTION Factorial(n)
!
!Purpose: Return the factorial of an integer n
!  
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
INTEGER :: k
Factorial = 1

DO k = 2, n
	Factorial = Factorial * k
END DO

END FUNCTION Factorial


REAL(KIND=rKind) FUNCTION BinomialCoef(n,m)
!
!Purpose: Return the Binomial Coefficient _nC_m
!
IMPLICIT NONE  
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(IN) :: m
INTEGER :: k
BinomialCoef=1.0_rKind

DO k=1,m,1
	BinomialCoef=BinomialCoef*(n-k+1)*1.0_rKind/(k*1.0_rKind)		
END DO

END FUNCTION BinomialCoef



END MODULE LinearOps
