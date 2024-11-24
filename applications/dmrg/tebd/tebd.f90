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
PROGRAM TEBDInterface
!
! Purpose: Main program interface between Open Source TEBD and the ALPS libraries
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!       8/18/10  M. L. Wall	alpha release
!       9/13/10  M. L. Wall	Changed quantum number i/o
!
USE GlobalData
USE LinearOps
USE HamiOps
USE StateOps
USE LocalOps
USE ObOps
USE PropOps
USE Hdf5Interface
USE omp_lib

IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:) !List of gamma tensors
TYPE(matrix), POINTER :: Gammasm(:) !List of gamma tensors
TYPE(vector), POINTER :: Lambdas(:) !List of Lambda vectors
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:) !Lists of number conserving vectors
TYPE(matrix), POINTER :: H(:) !Hamiltonian
TYPE(matrix) :: rho
TYPE(fourordProp) :: Urtp
TYPE(secordProp) :: Urtp2
TYPE(spinParams) :: sp
TYPE(HamiParams) :: Hparams
TYPE(HamiParamsList) :: HparamsList
TYPE(rtpData) :: rtpD
COMPLEX(KIND=rKind), ALLOCATABLE :: carray(:,:)
TYPE(measure) :: Measures !Measures derived type
REAL(KIND=rKIND) :: tick, tock, energy, time, localTruncerr,rtpTruncLimit, totalTruncerr !Timing Variables
INTEGER :: i,j,k,l,p,alpha,beta, numProcs, counter,rtpChiLimit, inCounter !Dummy integers
INTEGER :: overlapnum
CHARACTER(16) :: iWstring, mstring, anothermString
CHARACTER(64) :: fmtName
!Read in input parameters
NAMELIST /SystemSettings/ systemSize,HamiType, initialState, &
rtp, qSwitch,qtype, totQ, numThr,chiLimit, truncLimit, simId, print_switch
NAMELIST /ITPsettings/ numITP, itpFileName
NAMELIST /RTPsettings/ numQuenches, rtpFileName
NAMELIST /spinpar/ sp


PRINT *, " Time-Evolving Block Decimation Program"
PRINT *, "  available from http://alps.comp-phys.org/"
PRINT *, "  copyright(c) 2010 by Michael L. Wall <mwall@mines.edu>"
PRINT *, "                       Lincoln D. Carr <lcarr@mines.edu>"
PRINT *, "  For details see the publication: "
PRINT *, "  B. Bauer et al., J. Stat. Mech. (2011) P05001."
PRINT *, ''
PRINT *, " Based on the ALPS libraries version 2.0"!ALPS_VERSION,
PRINT *, "  available from http://alps.comp-phys.org/"
PRINT *, "  copyright (c) 1994-2011 by the ALPS collaboration."
PRINT *, "  Consult the web page for license details."
PRINT *, "  For details see the publications:"
PRINT *, "  A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007)."
PRINT *, "  B. Bauer et al., J. Stat. Mech. (2011) P05001."


CALL CPU_TIME(tick)
overlapnum=1
!Get system Settings
CALL GET_COMMAND_ARGUMENT(1,nmlName)
OPEN(138,file=nmlName)
READ(138,SystemSettings)
CLOSE(138)

rtpChiLimit=chiLimit
rtpTruncLimit=truncLimit

#if defined(_OPENMP)
!Set number of openmp threads
numProcs=omp_get_num_procs()
IF(numThr.le.numProcs) THEN
	CALL omp_set_num_threads(numThr)
ELSE
	PRINT *, 'Warning: Number of requested threads exceeds number of processors!'
	PRINT *, 'Setting number of threads to',numProcs
	numThr=numProcs
	CALL omp_set_num_threads(numThr)
END IF
#else
numProcs=1
numThr=1
#endif

!Output SystemSettings to hdf5
CALL WriteSystemSettings()

!Setup model/basis operators
CALL DefineModel(HamiType)

!Define Initial State
IF(initialState=='ground') THEN
	IF(print_switch) THEN
		PRINT *, 'Beginning ground state calculation!'
	END IF
	itp=.true.
	!get number of sweeps and name of itp input file
	OPEN(138,file=nmlName)
	READ(138,ITPsettings)
	CLOSE(138)
	!Read in ITP convergence params
	ALLOCATE(chivals(numITP), dtITPvals(numITP), convCriterion(numITP))
	WRITE(mString,'(I4)') numItp
	OPEN(139, file=itpFileName)
	fmtName='('//TRIM(ADJUSTL(mString))//'I16)'
	READ(139,fmtname) (chivals(j),j=1,numITP)
	fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'
	READ(139,fmtname) (dtITPvals(j),j=1,numITP)
	READ(139,fmtname) (convCriterion(j),j=1,numITP)
	CLOSE(139)

	!Read in Hamiltonian parameters
	CALL SetHamiltonian(H,HamiType, Hparams)
	!Write to Parameter group of output file
	CALL WriteITPSettings(HamiType, Hparams)

	!Quantum number conserving initial State
	IF(qSwitch) THEN
		IF(idof) THEN
			CALL AllocateGamLam(Gammas,Lambdas,chivals(1))
			CALL AllocateLabel(LabelLeft,LabelRight,chivals(1))
			CALL InitialSetNC(Gammas,Lambdas,LabelLeft, LabelRight, intdegfree=1)
			CALL ImagTimePropNC(H,Gammas,Lambdas, LabelLeft, LabelRight)
		ELSE
			CALL AllocateGamLam(Gammasm,Lambdas,chivals(1))
			CALL AllocateLabel(LabelLeft,LabelRight,chivals(1))
			CALL InitialSetNC(Gammasm,Lambdas,LabelLeft, LabelRight)
			CALL ImagTimePropNC(H,Gammasm,Lambdas, LabelLeft, LabelRight)
		END IF
	!No quantum numbers conserved
	ELSE
		CALL AllocateGamLam(Gammas,Lambdas,chivals(1))
		CALL AllStates(Gammas,Lambdas)
		CALL ImagTimeProp(H,Gammas,Lambdas)
	END IF
	IF(print_switch) THEN
		PRINT *, 'Ground state calculation finished!'
	END IF
ELSE IF(initialState=='kink') THEN
	qSwitch=.true.
	totQ=FLOOR(0.5_rKind*systemSize)
	CALL AllocateGamLam(Gammasm, Lambdas, 1)
	CALL AllocateLabel(LabelLeft,LabelRight,1)
	!Allocate a matrix to imprint the initial state
	ALLOCATE(carray(localSize,systemSize))
	!Define the kink state uuuuuu...ud...dddddd
	DO i=1,systemSize,1
		carray(:,i)=0.0_rKind
		IF(i.le.FLOOR(0.5_rKind*systemSize)) THEN
		carray(1,i)=1.0_rKind
		ELSE
		carray(localSize,i)=1.0_rKind
		END IF
	END DO
	CALL ProductStateMPD(Gammasm,LabelLeft,LabelRight, Lambdas, carray)
ELSE
	PRINT *, 'InitialState not recognized!'
	PRINT *, "Use 'ground' or 'kink'!"
	STOP
END IF 


!Setup RTP
!Get # of quenches
OPEN(138,file=nmlName)
READ(138,RTPsettings)
CLOSE(138)

IF(numQuenches.gt.0) THEN
	ALLOCATE(rtpD%tau(numQuenches), rtpD%pow(numQuenches), rtpD%gi(numQuenches), &
	rtpD%gf(numQuenches),rtpD%gquench(numQuenches),rtpD%nsteps(numQuenches),rtpD%stepsForStore(numQuenches))
	ALLOCATE(rtpD%nparams(numQuenches))
	!Read in RTP params
	WRITE(mString,'(I4)') numQuenches
	OPEN(137, file=rtpFileName)
	!Read in those things indexed by numQuenches
	fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'
	READ(137,fmtname) (rtpD%tau(j),j=1,numQuenches)
	fmtName='('//TRIM(ADJUSTL(mString))//'I16)'
	READ(137,fmtname) (rtpD%nsteps(j),j=1,numQuenches)
	READ(137,fmtname) (rtpD%stepsForStore(j),j=1,numQuenches)
	READ(137,fmtname) (rtpD%nparams(j),j=1,numQuenches)
	DO j=1,numQuenches
		ALLOCATE(rtpD%gquench(j)%v(rtpD%nparams(j)),rtpD%pow(j)%v(rtpD%nparams(j)),&
		rtpD%gi(j)%v(rtpD%nparams(j)),rtpD%gf(j)%v(rtpD%nparams(j)))
	END DO
	DO j=1,numQuenches
		WRITE(anothermString,'(I4)') rtpD%nparams(j)
		fmtName='('//TRIM(ADJUSTL(anothermString))//'E30.15)'
		READ(137,fmtname) (rtpD%pow(j)%v(p),p=1,rtpD%nparams(j))
	END DO
	DO j=1,numQuenches
		WRITE(anothermString,'(I4)') rtpD%nparams(j)
		fmtName='('//TRIM(ADJUSTL(anothermString))//'E30.15)'
		READ(137,fmtname) (rtpD%gi(j)%v(p),p=1,rtpD%nparams(j))
	END DO
	DO j=1,numQuenches
		WRITE(anothermString,'(I4)') rtpD%nparams(j)
		fmtName='('//TRIM(ADJUSTL(anothermString))//'E30.15)'
		READ(137,fmtname) (rtpD%gf(j)%v(p),p=1,rtpD%nparams(j))
	END DO
	DO j=1,numQuenches
		WRITE(anothermString,'(I4)') rtpD%nparams(j)
		fmtName='('//TRIM(ADJUSTL(anothermString))//'A10)'
		READ(137,fmtname) (rtpD%gquench(j)%v(p),p=1,rtpD%nparams(j))
	END DO
	CLOSE(137)

	!Get Initial Hamiltonian
	CALL SetHamiltonian(H,HamiType, Hparams)

	!Setup quenches
	CALL SetHamiltonianList(H,HamiType, Hparams,HparamsList, rtpD)
	!Write RTP parameters to parameter group of output file
	CALL WriteRTPSettings(HamiType,Hparams,rtpD)

	!Setup measurements
	CALL SetupMeasures(Measures,HamiType)
	IF(qswitch) THEN
		IF(idof) THEN
			CALL AllocateOverlap(Measures, overlapnum, Gammas, Lambdas) 
		ELSE
			CALL AllocateOverlap(Measures, overlapnum, Gammasm, Lambdas,LabelLeft, LabelRight) 
		END IF
	ELSE
		CALL AllocateOverlap(Measures, overlapnum, Gammas, Lambdas) 
	END IF

END IF


chiLimit=rtpChiLimit
truncLimit=rtpTruncLimit

!Allocate Propagator
CALL AllocateProp(Urtp2) !second order
time=0.0_rKind
totalTruncerr=0.0_rKind
counter=1
inCounter=1
DO i=1,numQuenches
	IF(i==1) THEN
		IF(print_switch) THEN
			PRINT *, 'Beginning real time propagation!'
		END IF
	END IF
	!Setup dt
	dtRTP=rtpD%tau(i)/(2.0_8*rtpD%nsteps(i)*1.0_rKind)
	IF(print_switch) THEN
		PRINT *, 'dtrtp', dtrtp
	END IF
	DO j=1,rtpD%nsteps(i)
		!Quench Hamiltonian and form propagator
		CALL UpdateHamiltonian(H,HamiType,HparamsList,counter)
		CALL ConstructPropagators(H,Urtp2,dtrtp)
		!Propagate forward one time step
		IF(qswitch) THEN
			!NC step
			IF(idof) THEN
				CALL TrotterStepNC(Urtp2, Gammas, Lambdas,LabelLeft,LabelRight, localTruncerr)
			ELSE
				CALL TrotterStepNC(Urtp2, Gammasm, Lambdas,LabelLeft,LabelRight, localTruncerr)
			END IF			
		ELSE
			CALL TrotterStep(Urtp2, Gammas, Lambdas, localTruncerr)
		END IF
		totalTruncerr=totalTruncerr+localTruncerr
		!Measurement step
		IF(MOD(counter,rtpD%stepsForStore(i))==0) THEN
			CALL SplitOperatorList(H,Measures%HMPO)
			!Bring into fully canonical form
			IF(qswitch) THEN
				IF(idof) THEN
					CALL CanonicalizeNC(Gammas,Lambdas,LabelLeft, LabelRight)
				ELSE
					CALL CanonicalizeNC(Gammasm,Lambdas,LabelLeft, LabelRight)
				END IF			
			ELSE
				CALL Canonicalize(Gammas,Lambdas)
			END IF			
			!Calculate measurements
			IF(qswitch.and..not.idof) THEN
				CALL EvaluateMeasures(Measures, Gammasm, Lambdas,LabelLeft)
			ELSE
				CALL EvaluateMeasures(Measures, Gammas, Lambdas)
			END IF
			!Output measurements
			CALL WriteMeasures(Measures,HamiType,inCounter, time,totalTruncerr, HparamsList,counter)
			inCounter=inCounter+1
			IF(print_switch) THEN
				PRINT *, 'time & truncation error',time, totalTruncerr
			END IF
			CALL DestroySplitOperator(Measures%HMPO)
		END IF
		time=time+2.0_8*ABS(dtRtp)
		counter=counter+1
	END DO
	IF(i==numQuenches) THEN
		IF(print_switch) THEN
			PRINT *, 'Finished real time propagation!'
		END IF
	END IF
END DO
CALL DeallocateProp(Urtp2)


!Clean up
IF(qswitch) THEN
	IF(idof) THEN
		CALL DeallocateGamLam(Gammas,Lambdas)
	ELSE
		CALL DeallocateGamLam(Gammasm,Lambdas)
	END IF
	CALL DeallocateLabel(LabelLeft,LabelRight)
ELSE
	CALL DeallocateGamLam(Gammas,Lambdas)
END IF

CALL DeallocateOps(H,systemSize-1)
CALL DestroyModel(Hamitype)
CALL DeallocateMeasures(Measures)
CALL FinalizeHdf5()

CALL CPU_TIME(tock)
IF(print_switch) THEN
	PRINT *, 'elapsed time', tock-tick, 'seconds.  Total truncation error is', totaltruncErr
END IF

!Clean up files
INQUIRE(file=rtpFileName,EXIST=fileExist)
IF(fileExist) THEN
	OPEN(137,file=rtpFileName)
	CLOSE(137, STATUS='DELETE')
END IF
INQUIRE(file=nmlName,EXIST=fileExist)
IF(fileExist) THEN
	OPEN(138,file=nmlName)
	CLOSE(138, STATUS='DELETE')
END IF
INQUIRE(file=itpFileName,EXIST=fileExist)
IF(fileExist) THEN
	OPEN(139,file=itpFileName)
	CLOSE(139, STATUS='DELETE')
END IF



END PROGRAM TEBDInterface
