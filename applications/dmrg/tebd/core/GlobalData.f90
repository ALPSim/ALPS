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
MODULE GlobalData
!
! Purpose: Module to declare global variables for OpenSourceTEBD v3.0
!
! Record of Revisions
!	Date		Programmer	Description of change
!	====		==========	=====================
!       8/18/10  M. L. Wall	alpha release
!
IMPLICIT NONE

!Select the precision and range of declared real variables
INTEGER, PARAMETER :: precis=15
INTEGER, PARAMETER :: range=30
INTEGER, PARAMETER :: rKind=SELECTED_REAL_KIND(p=precis,r=range)

! *** GLOBAL INPUT PARAMETERS
INTEGER :: systemSize !Number of lattice sites
INTEGER :: Nmax !Maximum number of particles allowed per lattice site
INTEGER :: localSize !On-site hilbert space dimension
INTEGER :: trotterOrder=2 !Order of the trotter decomposition, can be 2 or 4 (hardwired to 2 for ALPS)
INTEGER :: numITP !Number of ITP iterations
INTEGER, ALLOCATABLE :: chivals(:) !Entanglement cutoffs for ITP sweeps
REAL(KIND=rKind), ALLOCATABLE :: dtITPvals(:) !Time steps for ITP sweeps
INTEGER :: stepsForJudge=100 !Number of time steps before ITP convergence is checked
INTEGER :: maxITPsteps=2000 !greatest allowed number of ITP steps
REAL(KIND=rKind), ALLOCATABLE ::  convCriterion(:) !Convergence criteria for ITP sweeps 
REAL(KIND=rKIND) :: truncLimit=10.0_rKind**(-12)
COMPLEX(KIND=rKind) :: dtRTP=0.1_rKind !Time step for RTP
INTEGER ::  totalStep=10000 !Total number of RTP steps
INTEGER :: stepsForStore=10 !Number of RTP steps between each output write
INTEGER :: numThr !number of openmp threads
INTEGER :: simId !unique simulation ID

! *** I/O Data
CHARACTER(32) :: itpDir='ITPDATA/' !Directory where ITP data is stored
CHARACTER(32) :: rtpDir='RTPDATA/' !Directory where ITP data is stored
CHARACTER(132) :: HamiType !Name of desired Hamiltonian
CHARACTER(132) :: qType !Type of quantum number conserved
CHARACTER(132) :: initialState !Name of desired Hamiltonian
CHARACTER(132) :: nmlName !NameList input name
CHARACTER(132) :: itpFileName !NameList input name
CHARACTER(132) :: rtpFileName !NameList input name
CHARACTER(132) :: outputName !File for output

! *** Switches
LOGICAL :: itp=.FALSE. !Specifies whether ITP is desired
LOGICAL :: fileExist !flag for checking existence of files
LOGICAL :: rtp=.FALSE. !Specifies whether RTP is desired
LOGICAL :: print_switch=.TRUE. !Toggle printing to screen
LOGICAL :: qSwitch=.FALSE. !Toggle number conservation
LOGICAL :: idof=.FALSE. !Carry around on-site dimension for qc code
LOGICAL :: readITP=.FALSE. !Toggle read in of Initial state (Phase Diagram code)
LOGICAL :: writeITP=.FALSE. !Toggle output of final state (Phase diagram Code)
LOGICAL :: readRTP=.FALSE. !Toggle read in of Initial state (Phase Diagram code)
LOGICAL :: writeRTP=.FALSE. !Toggle output of final state (Phase diagram Code)
CHARACTER :: ITPopenKind='B' !'B' means output the MPS in binary, 'S' means output the MPS in scientific notation
CHARACTER :: RTPopenKind='B' !'B' means output the MPS in binary, 'S' means output the MPS in scientific notation

! ** Status checks
INTEGER :: statInt !Status integer to ensure allocation/deallocation happens properly
INTEGER :: fileStatus !Status integer to ensure file open happens properly

! *** Symmetry conservation parameters
INTEGER :: totQ=4 !Total conserved quantity Q

!Global RTP data
INTEGER :: numQuenches !Number of quenches
INTEGER :: chiLimit


END MODULE GlobalData
