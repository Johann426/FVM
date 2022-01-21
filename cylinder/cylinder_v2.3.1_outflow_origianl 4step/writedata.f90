    SUBROUTINE WRITEDATA

	USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
    
    IF(MOD(NITER_TIME,NITER_W).EQ.0) THEN

	NUM_WRITE = NUM_WRITE + 1
	
	CONTOUR = 'OOOO_OOOO.PLT'
	
	WRITE(CONTOUR(1:4),'(I4.4)') INT4(NUM_WRITE)
	
!   INTERPOLATION FOR NODAL RESULTS
	K=KMON
	DO J=1,NJ
	DO I=1,NI
		PG(I,J,K) = 0.25D0*(P(I,J,K)+P(I,J-1,K)+P(I-1,J,K)+P(I-1,J-1,K))
		UG(I,J,K) = 0.25D0*(U(I,J,K)+U(I,J-1,K)+U(I-1,J,K)+U(I-1,J-1,K))
		VG(I,J,K) = 0.25D0*(V(I,J,K)+V(I,J-1,K)+V(I-1,J,K)+V(I-1,J-1,K))
		TG(I,J,K) = 0.25D0*(T(I,J,K)+T(I,J-1,K)+T(I-1,J,K)+T(I-1,J-1,K))
	ENDDO
	ENDDO
! BRANCH CUT
	DO I=1,NI
	J=1
		PG(I,J,K) = 0.5D0 * ( PG(I,J,K) + PG(I,NJ,K) )
		UG(I,J,K) = 0.5D0 * ( UG(I,J,K) + UG(I,NJ,K) )
		VG(I,J,K) = 0.5D0 * ( VG(I,J,K) + VG(I,NJ,K) )
        TG(I,J,K) = 0.5D0 * ( TG(I,J,K) + TG(I,NJ,K) )
        
		PG(I,NJ,K) = PG(I,J,K)
		UG(I,NJ,K) = UG(I,J,K)
		VG(I,NJ,K) = VG(I,J,K)
		TG(I,NJ,K) = TG(I,J,K)
	ENDDO

!------------------------------------------------------------------------
!   WRITING LOOP
	CALL CALC_VOR ! VORTICITY POST-PROCESSING
!   NODAL RESULTS
!	WRITE(CONTOUR(6:9),'(A4)') 'NODE'
!	OPEN(UNIT=1000,FILE=CONTOUR,STATUS='UNKNOWN')
!	WRITE(1000,*) 'TITLE="','RE',RE,'PR',PR,'"'
!	WRITE(1000,*) 'variables=x,y,u,v,P,Vorticity,T'
!    WRITE(1000,24) NI,NJ,RTIME
!24  FORMAT(TR1'ZONE T="WHOLE", I=',I6,', j=',I6,', SOLUTIONTIME=',F10.4)
!	K=KMON
!	DO J=1,NJ
!	DO I=1,NI
!		WRITE(1000,*) XG(I,J,K),YG(I,J,K),UG(I,J,K),VG(I,J,K),PG(I,J,K),TG(I,J,K)
!	ENDDO
!	ENDDO
!	CLOSE(1000)
!   CELL CENTER RESULTS
	WRITE(CONTOUR(6:9),'(A4)') 'CENT'
	OPEN(UNIT=100000,FILE=CONTOUR,STATUS='UNKNOWN')
	WRITE(100000,*) 'TITLE="','RE',RE,'PR',PR,'"'
	WRITE(100000,*) 'variables=x,y,u,v,P,Vorticity,T'
    WRITE(100000,25) NI+1,NJ+1,RTIME
25  FORMAT(TR1'ZONE T="WHOLE", I=',I6,', j=',I6,', SOLUTIONTIME=',F10.4)
	K=KMON
	DO J=0,NJ
	DO I=0,NI
		WRITE(100000,*) X(I,J,K),Y(I,J,K),U(I,J,K),V(I,J,K),P(I,J,K),VOR(I,J,K),T(I,J,K)
	ENDDO
	ENDDO
	CLOSE(100000)

    ENDIF
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!   RESTART FILE
    IF(MOD(NITER_TIME,NITER_R).EQ.0) THEN
    OPEN(UNIT=16,FILE='RESTART.DAT',STATUS='UNKNOWN')
        WRITE(16,*) RE,PR,DTIME,VEL_ROT
        WRITE(16,*) RTIME,NITER_TIME
        DO K=2,NKM1
	    DO J=SY,EY
	    DO I=SX,EX
		    WRITE(16,*) NLX(I,J,K),NLY(I,J,K),ADV(I,J,K)
	    	WRITE(16,*) U(I,J,K),V(I,J,K),P(I,J,K),T(I,J,K)
    	ENDDO
	    ENDDO
	    ENDDO
    CLOSE(16)
    ENDIF
!------------------------------------------------------------------------
	END SUBROUTINE WRITEDATA