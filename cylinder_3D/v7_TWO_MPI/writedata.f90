    SUBROUTINE WRITEDATA

	USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
    
    N=MOD(NITER_TIME,NITER_W)
    
    IF(N.EQ.0) THEN
	
	NUM_WRITE = NUM_WRITE + 1
	
	CALL CALC_VOR ! VORTICITY POST-PROCESSING
	
!	CALL MPI_X(VOR_Z)
	
	CALL MPI_X(VOR_Y)
	
	CALL MPI_X(VOR_X)
	
!   INTERPOLATION FOR NODAL RESULTS
	DO K=KS1,KE
	DO J=JS1,JE
	DO I=IS1,IE
		PG(I,J,K) = 0.125D0*(P(I,J,K)+P(I,J-1,K)+P(I-1,J,K)+P(I-1,J-1,K)		&
                        +P(I,J,K-1)+P(I,J-1,K-1)+P(I-1,J,K-1)+P(I-1,J-1,K-1))
		UG(I,J,K) = 0.125D0*(U(I,J,K)+U(I,J-1,K)+U(I-1,J,K)+U(I-1,J-1,K)		&
		                +U(I,J,K-1)+U(I,J-1,K-1)+U(I-1,J,K-1)+U(I-1,J-1,K-1))
		VG(I,J,K) = 0.125D0*(V(I,J,K)+V(I,J-1,K)+V(I-1,J,K)+V(I-1,J-1,K)		&
		                +V(I,J,K-1)+V(I,J-1,K-1)+V(I-1,J,K-1)+V(I-1,J-1,K-1))
		WG(I,J,K) = 0.125D0*(W(I,J,K)+W(I,J-1,K)+W(I-1,J,K)+W(I-1,J-1,K)		&
		                +W(I,J,K-1)+W(I,J-1,K-1)+W(I-1,J,K-1)+W(I-1,J-1,K-1))
		TG(I,J,K) = 0.125D0*(T(I,J,K)+T(I,J-1,K)+T(I-1,J,K)+T(I-1,J-1,K)		&
		                +T(I,J,K-1)+T(I,J-1,K-1)+T(I-1,J,K-1)+T(I-1,J-1,K-1))

	ENDDO
	ENDDO
	ENDDO

!	BRANCH CUT
	DO K=KS1,KE
	DO I=IS1,IE
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
	
	CONTOUR = 'OOOO_OO.PLT'
	
	WRITE(CONTOUR(1:4),'(I4.4)') INT4(NUM_WRITE)
	
	WRITE(CONTOUR(6:7),'(I2.2)') INT4(RANK)

!	!   NODAL RESULTS
!	OPEN(UNIT=1000+RANK,FILE=CONTOUR,STATUS='UNKNOWN')
!	WRITE(1000+RANK,*) 'TITLE="','RE',RE,'PR',PR,'"'
!	WRITE(1000+RANK,*) 'variables=x,y,z,u,v,w,P,T,vorX,vorY,vorZ'
!    WRITE(1000+RANK,24) NI,NJ,KE-KS,RTIME
!24  FORMAT(TR1'ZONE T="WHOLE", I=',I6,', j=',I6,', k=',I6,', SOLUTIONTIME=',F10.4)
!	DO K=KS1,KE
!    DO J=JS1,JE
!    DO I=IS1,IE
!		WRITE(1000+RANK,*) XG(I,J,K),YG(I,J,K),ZG(I,J,K),UG(I,J,K),VG(I,J,K),WG(I,J,K),PG(I,J,K),TG(I,J,K),VOR_X(I,J,K),VOR_Y(I,J,K),VOR_Z(I,J,K)
!	ENDDO
!	ENDDO
!	ENDDO
!	CLOSE(1000+RANK)
	
	!   CELL CENTER RESULTS
	OPEN(UNIT=100000+RANK,FILE=CONTOUR,STATUS='UNKNOWN')
	WRITE(100000+RANK,*) 'TITLE="','RE',RE,'PR',PR,'"'
	WRITE(100000+RANK,*) 'variables=x,y,z,u,v,w,P,T,vorX,vorY,vorZ'
	WRITE(100000+RANK,25) NI+1,NJ+1,KE-KS+1,RTIME
25	FORMAT(TR1'ZONE T="WHOLE", I=',I6,', J=',I6,', K=',I6,', SOLUTIONTIME=',F10.4)
	DO K=KS,KE
	DO J=JS,JE
	DO I=IS,IE
		WRITE(100000+RANK,*) X(I,J,K),Y(I,J,K),Z(I,J,K),U(I,J,K),V(I,J,K),W(I,J,K),P(I,J,K),T(I,J,K),VOR_X(I,J,K),VOR_Y(I,J,K),VOR_Z(I,J,K)
	ENDDO
	ENDDO
	ENDDO
	CLOSE(100000+RANK)
	
	ENDIF ! IF(MOD(NITER_TIME,NITER_W).EQ.0) THEN
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!   RESTART FILE
    IF(RTIME.GE.RTIME_MAX) THEN
    
    FILENAME2 = 'RESTART_OO.DAT'
	WRITE(FILENAME2(9:10),'(I2.2)') INT4(RANK)
	OPEN(UNIT=16+RANK,FILE=FILENAME2,STATUS='UNKNOWN')
        WRITE(16+RANK,*) RE,PR,DTIME	!,VEL_ROT
        WRITE(16+RANK,*) RTIME,NITER_TIME
		DO K=KS,KE
		DO J=JS,JE
		DO I=IS,IE
		    WRITE(16+RANK,*) NLX(I,J,K),NLY(I,J,K),NLZ(I,J,K),ADV(I,J,K)
	    	WRITE(16+RANK,*) U(I,J,K),V(I,J,K),W(I,J,K),P(I,J,K),T(I,J,K)
    	ENDDO
	    ENDDO
	    ENDDO
    CLOSE(16+RANK)
    	
    ENDIF	! IF(RTIME.GE.RTIME_MAX) THEN
!------------------------------------------------------------------------

	END SUBROUTINE WRITEDATA