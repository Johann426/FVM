	SUBROUTINE CONTRA_VEL

	USE COMDAT_SHARED 
      
	IMPLICIT NONE
    SAVE
!--------------------------------------------------------------------------------------------------------
!   CONTRA-VARIENT VELOCITY
    U1 = XIX*U + XIY*V !+ XIZ*W
    U2 = ETX*U + ETY*V !+ ETZ*W
!   WC = ZTX*U + ZTY*V !+ ZTZ*W

!   VOLUME FLUX : CONTRA-VARIENT VELOCITY MULTIPLIED BY J^-1
    UE(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0*(DJAC(SX1:EX1,SY1:EY1,2:NKM1)*U1(SX1:EX1,SY1:EY1,2:NKM1)   &
                                +DJAC(SX1+1:EX1+1,SY1:EY1,2:NKM1)*U1(SX1+1:EX1+1,SY1:EY1,2:NKM1))
    UW(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0*(DJAC(SX1:EX1,SY1:EY1,2:NKM1)*U1(SX1:EX1,SY1:EY1,2:NKM1)   &
                                +DJAC(SX1-1:EX1-1,SY1:EY1,2:NKM1)*U1(SX1-1:EX1-1,SY1:EY1,2:NKM1))
    VN(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0*(DJAC(SX1:EX1,SY1:EY1,2:NKM1)*U2(SX1:EX1,SY1:EY1,2:NKM1)   &
                                +DJAC(SX1:EX1,SY1+1:EY1+1,2:NKM1)*U2(SX1:EX1,SY1+1:EY1+1,2:NKM1))
    VS(SX1:EX1,SY1:EY1,2:NKM1) = 0.5D0*(DJAC(SX1:EX1,SY1:EY1,2:NKM1)*U2(SX1:EX1,SY1:EY1,2:NKM1)   &
                                +DJAC(SX1:EX1,SY1-1:EY1-1,2:NKM1)*U2(SX1:EX1,SY1-1:EY1-1,2:NKM1))
    
!--------------------------------------------------------------------------------------------------------
!   CELL-FAVE VOLUME FLUX
    IF(FLAG) THEN
    AP(SX,SY1:EY1,2:NKM1) = AP(SX1,SY1:EY1,2:NKM1)
    AP(EX,SY1:EY1,2:NKM1) = AP(EX1,SY1:EY1,2:NKM1)
    AP(SX1:EX1,SY,2:NKM1) = AP(SX1:EX1,EY1,2:NKM1)  !AP(SX1:EX1,SY,2:NKM1) = AP(SX1:EX1,SY1,2:NKM1)
    AP(SX1:EX1,EY,2:NKM1) = AP(SX1:EX1,SY1,2:NKM1)  !AP(SX1:EX1,EY,2:NKM1) = AP(SX1:EX1,EY1,2:NKM1)

!   VOLUME FLUX NORMAL TO XI-CONSTANT SURFACE
    DO K=2,NKM1
    DO J=SY1,EY1
    DO I=SX1,EX1-1
    UE(I,J,K) = UE(I,J,K) +0.25D0*(					&
    +2.D0*RE*DJAC(I+1,J,K)*Q11(I+1,J,K)/AP(I+1,J,K)	&
    *(P(I+2,J,K)-2.D0*P(I+1,J,K)+P(I,J,K))			&
    -2.D0*RE*DJAC(I,J,K)*Q11(I,J,K)/AP(I,J,K)		&
    *(P(I+1,J,K)-2.D0*P(I,J,K)+P(I-1,J,K))			&
									)
    ENDDO
    ENDDO
    ENDDO

    DO K=2,NKM1
    DO J=SY1,EY1
    DO I=SX1+1,EX1
    UW(I,J,K) = UW(I,J,K) +0.25D0*(		            &
    +2.D0*RE*DJAC(I,J,K)*Q11(I,J,K)/AP(I,J,K)		&
    *(P(I+1,J,K)-2.D0*P(I,J,K)+P(I-1,J,K))			&
    -2.D0*RE*DJAC(I-1,J,K)*Q11(I-1,J,K)/AP(I-1,J,K) &
    *(P(I,J,K)-2.D0*P(I-1,J,K)+P(I-2,J,K))			&
									)
    ENDDO
    ENDDO
    ENDDO

!   VOLUME FLUX NORMAL TO ET-CONSTANT SURFACE
    DO K=2,NKM1
    DO J=SY1,EY1-1
    DO I=SX1,EX1
    VN(I,J,K) = VN(I,J,K) +0.25D0*(					&
    +2.D0*RE*DJAC(I,J+1,K)*Q22(I,J+1,K)/AP(I,J+1,K) &
    *(P(I,J+2,K)-2.D0*P(I,J+1,K)+P(I,J,K))			&
    -2.D0*RE*DJAC(I,J,K)*Q22(I,J,K)/AP(I,J,K)		&
    *(P(I,J+1,K)-2.D0*P(I,J,K)+P(I,J-1,K))			&
									)
    ENDDO
    VN(I,EY1,K) = VN(I,EY1,K) +0.25D0*(						&
    +2.D0*RE*DJAC(I,EY1+1,K)*Q22(I,EY1+1,K)/AP(I,EY1+1,K)	&
    *(P(I,SY1+1,K)-2.D0*P(I,EY1+1,K)+P(I,EY1,K))			&
    -2.D0*RE*DJAC(I,EY1,K)*Q22(I,EY1,K)/AP(I,EY1,K)         &
    *(P(I,EY1+1,K)-2.D0*P(I,EY1,K)+P(I,EY1-1,K))			&
										)
    ENDDO
    ENDDO

    DO K=2,NKM1
    DO J=SY1+1,EY1
    DO I=SX1,EX1
    VS(I,J,K) = VS(I,J,K) +0.25D0*(					&
    +2.D0*RE*DJAC(I,J,K)*Q22(I,J,K)/AP(I,J,K)		&
    *(P(I,J+1,K)-2.D0*P(I,J,K)+P(I,J-1,K))			&
    -2.D0*RE*DJAC(I,J-1,K)*Q22(I,J-1,K)/AP(I,J-1,K) &
    *(P(I,J,K)-2.D0*P(I,J-1,K)+P(I,J-2,K))			&
									)
    ENDDO
    VS(I,SY1,K) = VS(I,SY1,K) +0.25D0*(						&
    +2.D0*RE*DJAC(I,SY1,K)*Q22(I,SY1,K)/AP(I,SY1,K)			&
    *(P(I,SY1+1,K)-2.D0*P(I,SY1,K)+P(I,SY1-1,K))			&
    -2.D0*RE*DJAC(I,SY1-1,K)*Q22(I,SY1-1,K)/AP(I,SY1-1,K)	&
    *(P(I,SY1,K)-2.D0*P(I,SY1-1,K)+P(I,EY1-1,K))			&
										)
    ENDDO
    ENDDO
    ENDIF
!--------------------------------------------------------------------------------------------------------
!   BOUNDARY CONDITION FOR VOLUME FLUX
    UE(SX,SY:EY,2:NKM1) = 0.D0
    UW(SX1,SY:EY,2:NKM1) = 0.D0
!--------------------------------------------------------------------------------------------------------
	END SUBROUTINE CONTRA_VEL