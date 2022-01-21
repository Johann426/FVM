	SUBROUTINE BCV

	USE COMDAT_SHARED
	
	IMPLICIT NONE
	SAVE
!-----------------------------------------
!   CYLINDER SURFACE
!-----------------------------------------
	DO K=2,NKM1
	DO J=0,NJ
	I=0
		V_F = VEL_ROT * COS(THETA(I,J,K))
    AP(SX1,J,K) = AP(SX1,J,K) +AW(SX1,J,K)
    SU(SX1,J,K) = SU(SX1,J,K) +2.D0*AW(SX1,J,K)*V_F
    AW(SX1,J,K) = 0.0D0
	! CELL-CENTER BOUNDARY VALUE
	V(I,J,K)   = -V(I+1,J,K) + 2.0D0*V_F
	ENDDO
	ENDDO
!-------------------------------------------
!   FAR FIELD INLET & OUTLET REGION
!-------------------------------------------
	DO K=2,NKM1
	DO J=0,NJ
	I=NI
	    IF (THETA(I,J,K).GE.0.5D0*PI .AND. THETA(I,J,K).LE.1.5D0*PI) THEN
	    V_F = VIN
	    ELSE
	    V_F = 0.5D0*(V_N(I,J,K)+V_N(I-1,J,K)) - DTIME*Uin       &
				*(0.5D0*(XIX(I,J,K)+XIX(I-1,J,K))*(V_N(I,J,K)-V_N(I-1,J,K))			&
				+ 0.5D0*(ETX(I,J,K)+ETX(I-1,J,K))*0.25D0*(V_N(I,J+1,K)+V_N(I-1,J+1,K)-V_N(I,J-1,K)-V_N(I-1,J-1,K)))
	    ENDIF
    AP(NIM1,J,K) = AP(NIM1,J,K) +AE(NIM1,J,K)
    SU(NIM1,J,K) = SU(NIM1,J,K) +2.D0*AE(NIM1,J,K)*V_F
    AE(NIM1,J,K) = 0.D0
    ! CELL-CENTER BOUNDARY VALUE
	V(I,J,K) = -V(I-1,J,K) + 2.0D0*V_F
    ENDDO
	ENDDO
!-------------------------------------------
!   BRANCH CUT REGION
!-------------------------------------------
	DO K=2,NKM1
	DO I=0,NI
    J=0
	    V(I,J,K) = V(I,NJM1,K)
	J=NJ
		V(I,J,K) = V(I,1,K)
	ENDDO
	ENDDO
	
	END SUBROUTINE BCV