	SUBROUTINE INIT

    USE COMDAT_SHARED
	
	IMPLICIT NONE
	SAVE

	NITER=0
	NITER_TIME=0

	RTIME =0.0D0

    PI=4.0D0*DATAN(1.0D0)
    
!------------------------------------------------
!   INITIALIZING 
	U(:,:,:)   = 1.0D0
	U_N(:,:,:) = 1.0D0

	V(:,:,:)   = 0.0D0
	V_N(:,:,:) = 0.0D0

	T(:,:,:)   = 0.0D0
	T_N(:,:,:) = 0.0D0
    
	P(:,:,:) = 0.0D0
    
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
	    CLOSE(11)
	ENDIF
	
    IF(IRESTART.EQ.1) THEN
	    PRINT *,'Reading restart files...'
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
		
		RTIME_MAX = RTIME + 15.D0
		
	ENDIF
!------------------------------------------------------------------------------
	U_N(:,:,:)=U(:,:,:)
	V_N(:,:,:)=V(:,:,:)
    T_N(:,:,:)=T(:,:,:)
!---------------------------------------------------------------------------------
!   ANGULAR VECTOR (RADIAN)
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
!	WRITE(34,*) 'variables=x,y,THETA,J,Q11,Q12,Q21,Q22'
!	WRITE(34,*) 'zone t="',RTIME,'",i=',NI+1,',j=',NJ+1
!	DO K=2,NKM1
!	DO J=0,NJ
!	DO I=0,NI
!		WRITE(34,*) X(I,J,K),Y(I,J,K),THETA(I,J,K),DJAC(I,J,K),Q11(I,J,K),Q12(I,J,K),Q21(I,J,K),Q22(I,J,K)
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