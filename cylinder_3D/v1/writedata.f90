    SUBROUTINE WRITEDATA

	USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
    
    IF(MOD(NITER_TIME,NITER_W).EQ.0) THEN
	READ(22,*) FILENAME,FILENAME2
	MT = NITER_TIME / NITER_W

!   INTERPOLATION FOR NODAL RESULTS
	DO K=1,NK
	DO J=1,NJ
	DO I=1,NI
		PG(I,J,K) = 0.125D0*(P(I,J,K)+P(I,J-1,K)+P(I-1,J,K)+P(I-1,J-1,K)    &
                        +P(I,J,K-1)+P(I,J-1,K-1)+P(I-1,J,K-1)+P(I-1,J-1,K-1))
		UG(I,J,K) = 0.125D0*(U(I,J,K)+U(I,J-1,K)+U(I-1,J,K)+U(I-1,J-1,K)    &
		                +U(I,J,K-1)+U(I,J-1,K-1)+U(I-1,J,K-1)+U(I-1,J-1,K-1))
		VG(I,J,K) = 0.125D0*(V(I,J,K)+V(I,J-1,K)+V(I-1,J,K)+V(I-1,J-1,K)    &
		                +V(I,J,K-1)+V(I,J-1,K-1)+V(I-1,J,K-1)+V(I-1,J-1,K-1))
		WG(I,J,K) = 0.125D0*(W(I,J,K)+W(I,J-1,K)+W(I-1,J,K)+W(I-1,J-1,K)    &
		                +W(I,J,K-1)+W(I,J-1,K-1)+W(I-1,J,K-1)+W(I-1,J-1,K-1))
		TG(I,J,K) = 0.125D0*(T(I,J,K)+T(I,J-1,K)+T(I-1,J,K)+T(I-1,J-1,K)    &
		                +T(I,J,K-1)+T(I,J-1,K-1)+T(I-1,J,K-1)+T(I-1,J-1,K-1))
	ENDDO
	ENDDO
	ENDDO

!	BRANCH CUT
	DO K=1,NK
	DO I=1,NI
	J=1
		PG(I,J,K) = 0.5D0 * ( PG(I,J,K) + PG(I,NJ,K) )
		UG(I,J,K) = 0.5D0 * ( UG(I,J,K) + UG(I,NJ,K) )
		VG(I,J,K) = 0.5D0 * ( VG(I,J,K) + VG(I,NJ,K) )
		WG(I,J,K) = 0.5D0 * ( WG(I,J,K) + WG(I,NJ,K) )
        TG(I,J,K) = 0.5D0 * ( TG(I,J,K) + TG(I,NJ,K) )
        
		PG(I,NJ,K) = PG(I,J,K)
		UG(I,NJ,K) = UG(I,J,K)
		VG(I,NJ,K) = VG(I,J,K)
		WG(I,NJ,K) = WG(I,J,K)
		TG(I,NJ,K) = TG(I,J,K)
	ENDDO
	ENDDO
	
!------------------------------------------------------------------------
!   WRITING LOOP
	CALL CALC_VOR ! VORTICITY POST-PROCESSING
!   NODAL RESULTS
	OPEN(UNIT=1000+MT,FILE=FILENAME2,STATUS='UNKNOWN')
	WRITE(1000+MT,*) 'TITLE="','RE',RE,'PR',PR,'"'
	WRITE(1000+MT,*) 'variables=x,y,z,u,v,w,P,T'
    WRITE(1000+MT,24) NI,NJ,NK,RTIME
24  FORMAT(TR1'ZONE T="WHOLE", I=',I6,', j=',I6,', k=',I6,', SOLUTIONTIME=',F10.4)
	DO K=1,NK !K=KMON
	DO J=1,NJ
	DO I=1,NI
		WRITE(1000+MT,*) XG(I,J,K),YG(I,J,K),ZG(I,J,K),UG(I,J,K),VG(I,J,K),WG(I,J,K),PG(I,J,K),TG(I,J,K)
	ENDDO
	ENDDO
	ENDDO
	CLOSE(1000+MT)
!   CELL CENTER RESULTS
	OPEN(UNIT=100000+MT,FILE=FILENAME,STATUS='UNKNOWN')
	WRITE(100000+MT,*) 'TITLE="','RE',RE,'PR',PR,'"'
	WRITE(100000+MT,*) 'variables=x,y,z,u,v,w,P,Vor_Z,T'
    WRITE(100000+MT,25) NI+1,NJ+1,NK+1,RTIME
25  FORMAT(TR1'ZONE T="WHOLE", I=',I6,', j=',I6,', k=',I6,', SOLUTIONTIME=',F10.4)
	DO K=0,NK !K=KMON
	DO J=0,NJ
	DO I=0,NI
		WRITE(100000+MT,*) X(I,J,K),Y(I,J,K),Z(I,J,K),U(I,J,K),V(I,J,K),W(I,J,K),P(I,J,K),VOR(I,J,K),T(I,J,K)
	ENDDO
	ENDDO
	ENDDO
	CLOSE(100000+MT)

    ENDIF
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!   RESTART FILE
    IF(RTIME.GE.RTIME_MAX) THEN
    OPEN(UNIT=16,FILE='RESTART.DAT',STATUS='UNKNOWN')
        WRITE(16,*) RE,PR,DTIME,VEL_ROT
        WRITE(16,*) RTIME,NITER_TIME
        DO K=SZ,EZ
	    DO J=SY,EY
	    DO I=SX,EX
		    WRITE(16,*) NLX(I,J,K),NLY(I,J,K),NLZ(I,J,K),ADV(I,J,K)
	    	WRITE(16,*) U(I,J,K),V(I,J,K),W(I,J,K),P(I,J,K),T(I,J,K)
    	ENDDO
	    ENDDO
	    ENDDO
    CLOSE(16)
    ENDIF
!------------------------------------------------------------------------

	END SUBROUTINE WRITEDATA