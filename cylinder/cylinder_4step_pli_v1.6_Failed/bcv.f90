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
		V_F(SX1,J,K) = VEL_ROT * COS(THETA(I,J,K))

    AP(SX1,J,K) = AP(SX1,J,K) +AW(SX1,J,K)
    SU(SX1,J,K) = SU(SX1,J,K) +2.D0*AW(SX1,J,K)*V_F(SX1,J,K)
    AW(SX1,J,K) = 0.0D0
	! CELL-CENTER BOUNDARY VALUE
	V(I,J,K)   = -V(I+1,J,K) + 2.0D0*V_F(SX1,J,K)
	ENDDO
	ENDDO
!-------------------------------------------
!   FAR FIELD INLET & OUTLET REGION
!-------------------------------------------
	DO K=2,NKM1
	DO J=0,NJ
	I=NI
	    IF (THETA(I,J,K).GE.0.5D0*PI .AND. THETA(I,J,K).LE.1.5D0*PI) THEN
	    V_F(NIM1,J,K) = VIN

	    ELSE
	    V_F(NIM1,J,K) = 0.5D0*(V_N(I,J,K)+V_N(I-1,J,K)) - DTIME*Uin       &
                    *(0.5D0*(V_N(I,J,K)+V_N(I-1,J,K))-V_N(I-1,J,K)) &
                    /(0.5D0*(  X(I,J,K)+  X(I-1,J,K)) - X(I-1,J,K))	
	    ENDIF
    AP(NIM1,J,K) = AP(NIM1,J,K) +AE(NIM1,J,K)
    SU(NIM1,J,K) = SU(NIM1,J,K) +2.D0*AE(NIM1,J,K)*V_F(NIM1,J,K)
    AE(NIM1,J,K) = 0.D0
    ! CELL-CENTER BOUNDARY VALUE
	V(I,J,K) = -V(I-1,J,K) + 2.0D0*V_F(NIM1,J,K)
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