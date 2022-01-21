	SUBROUTINE INIT

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
	
	INCLUDE 'mpif.h'
	
	NITER=0
	NITER_TIME=0

	RTIME =0.0D0

    N=10000
	
	IF(RANK.EQ.0) THEN
	
    OPEN(UNIT=22,FILE='FILENAME.DAT',STATUS='UNKNOWN')
        DO I=1,N
        WRITE(22,22) I,I
22      FORMAT (TR1,I6,'_CENT.PLT',TR1,I6,'_NODE.PLT')
        ENDDO
    CLOSE(22)
    
    OPEN(UNIT=22,FILE='FILENAME.DAT',STATUS='OLD')
    OPEN(UNIT=23,FILE='FILENAME_NEW.DAT',STATUS='UNKNOWN')
        DO I=1,N
        READ(22,*) FILENAME,FILENAME2
        WRITE(23,23) IS,FILENAME,IS,FILENAME2
23      FORMAT (TR1,I1,A,TR1,I1,A)
        ENDDO
    CLOSE(22)
    CLOSE(23)
    
    ENDIF
    
	IF(RANK.EQ.0)	OPEN(UNIT=22,FILE='FILENAME_NEW.DAT',STATUS='OLD')

!------------------------------------------------------------------------------
!   INITIALIZING 
	U(IS:IE,JS:JE,KS:KE)   = 1.0D0
	V(IS:IE,JS:JE,KS:KE)   = 0.0D0
    W(IS:IE,JS:JE,KS:KE)   = 0.0D0
	
	T(IS:IE,JS:JE,KS:KE)   = 0.0D0
    
	P(IS:IE,JS:JE,KS:KE)   = 0.0D0
	PLI(IS:IE,JS:JE,KS:KE) = 0.0D0
    
    NLX(IS:IE,JS:JE,KS:KE) = 0.0D0
	NLY(IS:IE,JS:JE,KS:KE) = 0.0D0
	NLZ(IS:IE,JS:JE,KS:KE) = 0.0D0
	
	ADV(IS:IE,JS:JE,KS:KE) = 0.0D0

!   RESTART ROUTINE
	IF(IRESTART.EQ.0)	THEN
	IF(RANK.EQ.0)		THEN
		PRINT *,'Program start...'
		OPEN(UNIT=11,FILE='NU_CD_CL.PLT',STATUS='UNKNOWN')
	    WRITE(11,*) 'variables=TIME,NU,CD,CL,CPB,CD_P,CD_V,CL_P,CL_V,cfl_max'
	ENDIF
	ENDIF
	
    IF(IRESTART.EQ.1) THEN
	IF(RANK.EQ.0)		THEN
	    PRINT *,'Reading restart files...'
	    OPEN(UNIT=11,FILE='NU_CD_CL.PLT',ACCESS='APPEND')
	    OPEN(UNIT=16,FILE='RESTART.DAT',STATUS='OLD')
		READ(16,*) RE,PR,DTIME
		READ(16,*) RTIME,NITER_TIME
		DO K=KS,KE
		DO J=JS,JE
		DO I=IS,IE
			READ(16,*) NLX(I,J,K),NLY(I,J,K),NLZ(I,J,K),ADV(I,J,K)
			READ(16,*) U(I,J,K),V(I,J,K),W(I,J,K),P(I,J,K),T(I,J,K)
		ENDDO
		ENDDO
		ENDDO
		CLOSE(16)
	ENDIF
	ENDIF
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!   ANGULAR VECTOR (RADIAN)
	DO K=KS,KE
	DO J=JS,JE
	DO I=IS,IE
		RADIUS(I,J,K) = DSQRT(X(I,J,K)**2+Y(I,J,K)**2)
	ENDDO
    ENDDO
	ENDDO
	
	DO K=KS,KE
	DO J=JS,JE
	DO I=IS,IE
      THETA(I,J,K) = DATAN(Y(I,J,K)/X(I,J,K))
 	  IF     (THETA(I,J,K).GE.0.0D0.AND.DASIN(Y(I,J,K)/RADIUS(I,J,K)).GT.0.0D0) THEN
		THETA(I,J,K) = THETA(I,J,K)
 	  
	  ELSEIF (THETA(I,J,K).LT.0.0D0.AND.DASIN(Y(I,J,K)/RADIUS(I,J,K)).GT.0.0D0) THEN
		THETA(I,J,K) = PI+THETA(I,J,K)

 	  ELSEIF (THETA(I,J,K).GE.0.0D0.AND.DASIN(Y(I,J,K)/RADIUS(I,J,K)).LE.0.0D0) THEN
		THETA(I,J,K) = PI+THETA(I,J,K)

	  ELSEIF (THETA(I,J,K).LT.0.0D0.AND.DASIN(Y(I,J,K)/RADIUS(I,J,K)).LT.0.0D0) THEN
		THETA(I,J,K) = 2.0D0*PI+THETA(I,J,K)
	  ENDIF
	ENDDO
	ENDDO
    ENDDO
!-------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------
!   OPEN(UNIT=34,FILE='THETA.PLT',STATUS='UNKNOWN')
!   WRITE(34,*) 'title=""'
!	WRITE(34,*) 'variables=x,y,THETA'
!	WRITE(34,*) 'zone t="',RTIME,'",i=',NI+1,',j=',NJ+1
!	DO K=SZ,NK
!	DO J=SY,NJ
!	DO I=SX,NI
!		WRITE(34,*) X(I,J,K),Y(I,J,K),THETA(I,J,K)
!	ENDDO
!	ENDDO
!	ENDDO
!	CLOSE(34)
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!   SURFACE AREA ON THE CELL'S CENTER
    DUM = 0.D0	!SUM_AREA = 0.D0
    
    DO K=KS1,KE1
    DO J=JS1,JE1
    I=IS1
    L1 =DSQRT((XG(I,J+1,K)-XG(I,J,K))**2+(YG(I,J+1,K)-YG(I,J,K))**2+(ZG(I,J+1,K)-ZG(I,J,K))**2)
	L2 =DSQRT((XG(I,J,K+1)-XG(I,J,K))**2+(YG(I,J,K+1)-YG(I,J,K))**2+(ZG(I,J,K+1)-ZG(I,J,K))**2)
	L3 =DSQRT((XG(I,J+1,K)-XG(I,J,K+1))**2+(YG(I,J+1,K)-YG(I,J,K+1))**2+(ZG(I,J+1,K)-ZG(I,J,K+1))**2)
	TH1=DACOS((L3**2-L2**2-L1**2)/(-2.D0*L1*L2))
	AREA(I,J,K) = 0.5D0*L1*L2*DSIN(TH1)
    
    L1 =DSQRT((XG(I,J+1,K+1)-XG(I,J+1,K))**2+(YG(I,J+1,K+1)-YG(I,J+1,K))**2+(ZG(I,J+1,K+1)-ZG(I,J+1,K))**2)
	L2 =DSQRT((XG(I,J+1,K+1)-XG(I,J,K+1))**2+(YG(I,J+1,K+1)-YG(I,J,K+1))**2+(ZG(I,J+1,K+1)-ZG(I,J,K+1))**2)
!	L3 =DSQRT((XG(I,J+1,K)-XG(I,J,K+1))**2+(YG(I,J+1,K)-YG(I,J,K+1))**2+(ZG(I,J+1,K)-ZG(I,J,K+1))**2)
	TH1=DACOS((L3**2-L2**2-L1**2)/(-2.D0*L1*L2))
	AREA(I,J,K) = AREA(I,J,K) +0.5D0*L1*L2*DSIN(TH1)
	
	DUM = DUM +	AREA(I,J,K)	!SUM_AREA = SUM_AREA + AREA(I,J,K)
    
    ENDDO
    ENDDO
WRITE(*,*) RANK,'INI_ING'
	CALL MPI_ALLREDUCE(DUM,SUM_AREA,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ERR)

	IF(RANK.EQ.0) THEN
    WRITE (*,2) NI,NJ,NK,NUM_GRID,SUM_AREA
    WRITE (*,3) RE,PR,DTIME
2   FORMAT(TR1'NI =',I4,TR3'NJ =',I4,TR3'NK =',I4,TR4'Cells =',I7,TR3'Surface Area =',F10.6)
3   FORMAT(TR1'Re =',F10.2,TR4'Pr =',F7.2,TR4'dtime =',F10.4,TR12)
	ENDIF
WRITE(*,*) RANK,'INI_END'
!--------------------------------------------------------------------------------------
	END SUBROUTINE INIT