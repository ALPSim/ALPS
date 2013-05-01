C Copyright (c) 1993 by Seiji Miyashita
C (宮下精二著「熱・統計力学」(培風館 1993年) p.263)
C MONTE CARLO SIMULATION : THERMAL BATH METHOD
C   2D ISING MODEL (SQUARE LATTICE)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IS(20,20),IP(20),IM(20),P(-4:4),A(4)
C PARAMETERS
      DATA TEMP/2.5/, L/10/, MCS/1000/, INT/1000/
      DATA IX/1234567/, V0/.465661288D-9/
C TABLES
      DO 10 I=-4,4
      W=EXP(FLOAT(I)/TEMP)
 10   P(I)=W/(W+1/W)
      DO 11 I=1,L
      IP(I)=I+1
 11   IM(I)=I-1
      IP(L)=1
      IM(1)=L
C INITIAL CONFIGURATION
      DO 20 I=1,L
      DO 20 J=1,L
 20   IS(I,J)=1
C ACCUMULATION DATA RESET
      DO 21 I=1,4
 21   A(I)=0.0
C SIMULATION
      DO 30 K=1,MCS+INT
      KIJ=0
      DO 31 I=1,L
      DO 31 J=1,L
      M=IS(IP(I),J)+IS(I,IP(J))+IS(IM(I),J)+IS(I,IM(J))
      KIJ=KIJ+1
      IS(I,J)=-1
      IX=IAND(IX*5*11,2147483647)
      IF(P(M).GT.V0*IX) IS(I,J)=1
 31   CONTINUE
C DATA
      IF(K.LE.INT) GOTO 30
      EN=0
      MG=0
      DO 40 I=1,L
      DO 40 J=1,L
      EN=EN+IS(I,J)*(IS(IP(I),J)+IS(I,IP(J)))
 40   MG=MG+IS(I,J)
      A(1)=A(1)+EN
      A(2)=A(2)+EN**2
      A(3)=A(3)+MG
      A(4)=A(4)+MG**2
 30   CONTINUE
C STATISTICS
      DO 50 I=1,4
 50   A(I)=A(I)/MCS
      C=(A(2)-A(1)**2)/L**2/TEMP**2
      X=(A(4)-A(3)**2)/L**2/TEMP
      ENG=A(1)/L**2
      AMG=A(3)/L**2
      WRITE(6,100) TEMP,L,ENG,C,AMG,X
 100  FORMAT(' TEMP=',F10.5,' SIZE=',I5,
     * /' ENG =',F10.5,' C   =',F10.5,
     * /' MAG =',F10.5,' X   =',F10.5)
      END
