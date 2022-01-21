	SUBROUTINE INIT

    USE COMDAT_SHARED
	
	IMPLICIT NONE
	SAVE

	NITER=0
	NITER_TIME=0

	RTIME =0.0D0

    PI=4.0D0*DATAN(1.0D0)
    
    N=10000
    
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
        WRITE(23,23) SX,FILENAME,SX,FILENAME2
23      FORMAT (TR1,I1,A,TR1,I1,A)
        ENDDO
    CLOSE(22)
    CLOSE(23)
        OPEN(UNIT=22,FILE='FILENAME_NEW.DAT',STATUS='OLD')

!------------------------------------------------
!   INITIALIZING 
	U(:,:,:)   = 1.0D0
	U_N(:,:,:) = 1.0D0

	V(:,:,:)   = 0.0D0
	V_N(:,:,:) = 0.0D0

	T(:,:,:)   = 0.0D0
	T_N(:,:,:) = 0.0D0
    
	P(:,:,:) = 0.0D0
	PLI(:,:,:) = 0.0D0
    
    NLX(:,:,:) = 0.0D0
	NLY(:,:,:) = 0.0D0
!------------------------------------------------
!   BOUNDARY CONDITION
    CALL BCU
    CALL BCV
    CALL BCP
    CALL BCT
!------------------------------------------------------------------------------
!   RESTART ROUTINE
	IF(IRESTART.EQ.0) THEN
		PRINT *,'Program start...'
		OPEN(UNIT=11,FILE='NU_CD_CL.PLT',STATUS='UNKNOWN')
	    WRITE(11,*) 'variables=TIME,NU,CD,CL,CPB,CD_P,CD_V,CL_P,CL_V,cfl_max'
	ENDIF
	
    IF(IRESTART.EQ.1) THEN
	    PRINT *,'Reading restart files...'
	    OPEN(UNIT=11,FILE='NU_CD_CL.PLT',ACCESS='APPEND')
	    OPEN(UNIT=16,FILE='RESTART.DAT',STATUS='OLD')
		READ(16,*) RE,PR,DTIME,VEL_ROT
		READ(16,*) RTIME,NITER_TIME
		DO K=2,NKM1  
		DO J=sy1-1,ey1+1
		DO I=sx1-1,ex1+1
			READ(16,*) NLX(I,J,K),NLY(I,J,K),ADV(I,J,K)
			READ(16,*) U(I,J,K),V(I,J,K),P(I,J,K),T(I,J,K)
		ENDDO
		ENDDO
		ENDDO
    CLOSE(16)
	ENDIF
!------------------------------------------------------------------------------
	U_N(:,:,:)=U(:,:,:)
	V_N(:,:,:)=V(:,:,:)
    T_N(:,:,:)=T(:,:,:)
!---------------------------------------------------------------------------------
!   ANGULAR VECTOR
!   CELL NODE ANGULAR VECTOR IN RADIAN
    DO K=2,NKM1
	DO I=1,NI
	DO J=1,NJ
		RADIUSG(I,J,K) = (XG(I,J,K)**2+YG(I,J,K)**2)**0.5
	ENDDO
    ENDDO
	ENDDO
    DO K=2,NKM1
	DO I=1,NI
	DO J=2,NJ-1
      THETAG(I,J,K) = DATAN(YG(I,J,K)/XG(I,J,K)) ! CELL NODE ANGULAR VECTOR IN RADIAN 
 	  IF     (THETAG(I,J,K).GE.0.0D0.AND.DASIN(YG(I,J,K)/RADIUSG(I,J,K)).GT.0.0D0) THEN
		THETAG(I,J,K) = THETAG(I,J,K)
 	  
	  ELSEIF (THETAG(I,J,K).LT.0.0D0.AND.DASIN(YG(I,J,K)/RADIUSG(I,J,K)).GT.0.0D0) THEN
		THETAG(I,J,K) = PI+THETAG(I,J,K)

 	  ELSEIF (THETAG(I,J,K).GE.0.0D0.AND.DASIN(YG(I,J,K)/RADIUSG(I,J,K)).LE.0.0D0) THEN
		THETAG(I,J,K) = PI+THETAG(I,J,K)

	  ELSEIF (THETAG(I,J,K).LT.0.0D0.AND.DASIN(YG(I,J,K)/RADIUSG(I,J,K)).LT.0.0D0) THEN
		THETAG(I,J,K) = 2.0D0*PI+THETAG(I,J,K)
	  ENDIF
	ENDDO
	ENDDO
    ENDDO
    DO K=2,NKM1
	DO I=1,NI
	J=1
	THETAG(I,J,K) = 0.D0
	J=NJ
	THETAG(I,J,K) = 2.D0*PI
	ENDDO
	ENDDO

!   CELL CENTER ANGULAR VECTOR IN RADIAN
    DO K=2,NKM1
	DO I=0,NI
	DO J=0,NJ
		RADIUS(I,J,K) = (X(I,J,K)**2+Y(I,J,K)**2)**0.5
	ENDDO
    ENDDO
	ENDDO
    DO K=2,NKM1
	DO I=0,NI
	DO J=0,NJ
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
!	DO K=2,NKM1
!	DO J=0,NJ
!	DO I=0,NI
!		WRITE(34,*) X(I,J,K),Y(I,J,K),THETA(I,J,K)
!	ENDDO
!	ENDDO
!	ENDDO
!	CLOSE(34)	
!   OPEN(UNIT=34,FILE='THETAG.PLT',STATUS='UNKNOWN')
!   WRITE(34,*) 'title=""'
!	WRITE(34,*) 'variables=xG,yG,THETAG'
!	WRITE(34,*) 'zone t="',RTIME,'",i=',NI,',j=',NJ
!	DO K=2,NKM1
!	DO J=1,NJ
!	DO I=1,NI
!		WRITE(34,*) XG(I,J,K),YG(I,J,K),THETAG(I,J,K)
!	ENDDO
!	ENDDO
!	ENDDO
!	CLOSE(34)
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!   SURFACE AREA ON THE CELL'S CENTER
    SUM_AREA = 0.D0
    
    DO K=2,NKM1
    DO J=1,NJ-1
    I=1
    AREA(I,J,K) = ((XG(I,J+1,K)-XG(I,J,K))**2.D0+(YG(I,J+1,K)-YG(I,J,K))**2.D0)**0.5D0
    SUM_AREA = SUM_AREA + AREA(I,J,K)
    ENDDO
    ENDDO

    PE = RE*PR
    
    WRITE (*,2) NI,NJ,NUM_GRID,SUM_AREA
    WRITE (*,3) RE,PR,VEL_ROT,DTIME
2   FORMAT(TR1'NI =',I4,TR4'NJ =',I4,TR6'Cell number =',I7,TR4'Surface Area =',F10.8)    
3   FORMAT(TR1'Re =',F10.2,TR4'Pr =',F7.2,TR4'Alpha =',F6.2,TR4'dtime =',F10.4,TR12)
    PAUSE
!--------------------------------------------------------------------------------------
	END SUBROUTINE INIT