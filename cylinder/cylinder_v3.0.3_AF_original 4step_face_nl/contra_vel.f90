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

!   VOLUME FLUX NORMAL TO XI-CONSTANT SURFACE
    DO K=2,NKM1
    DO J=SY1,EY1
    DO I=SX1,EX1-1
    UE(I,J,K) = UE(I,J,K) +0.25D0               &
    *(                                              &
    +2.D0*RE*XIX(I+1,J,K)*DJAC(I+1,J,K)/AP(I+1,J,K)                                                                     &
    *(DJAC(I+2,J,K)*XIX(I+2,J,K)*P(I+2,J,K)-2.D0*DJAC(I+1,J,K)*XIX(I+1,J,K)*P(I+1,J,K)+DJAC(I,J,K)*XIX(I,J,K)*P(I,J,K)) &
    -2.D0*RE*XIX(I,J,K)*DJAC(I,J,K)/AP(I,J,K)                                                                           &
    *(DJAC(I+1,J,K)*XIX(I+1,J,K)*P(I+1,J,K)-2.D0*DJAC(I,J,K)*XIX(I,J,K)*P(I,J,K)+DJAC(I-1,J,K)*XIX(I-1,J,K)*P(I-1,J,K)) &
    +2.D0*RE*XIY(I+1,J,K)*DJAC(I+1,J,K)/AP(I+1,J,K)                                                                     &
    *(DJAC(I+2,J,K)*XIY(I+2,J,K)*P(I+2,J,K)-2.D0*DJAC(I+1,J,K)*XIY(I+1,J,K)*P(I+1,J,K)+DJAC(I,J,K)*XIY(I,J,K)*P(I,J,K)) &
    -2.D0*RE*XIY(I,J,K)*DJAC(I,J,K)/AP(I,J,K)                                                                           &
    *(DJAC(I+1,J,K)*XIY(I+1,J,K)*P(I+1,J,K)-2.D0*DJAC(I,J,K)*XIY(I,J,K)*P(I,J,K)+DJAC(I-1,J,K)*XIY(I-1,J,K)*P(I-1,J,K)) &
    )
    ENDDO
    ENDDO
    ENDDO

    DO K=2,NKM1
    DO J=SY1,EY1
    DO I=SX1+1,EX1
    UW(I,J,K) = UW(I,J,K) +0.25D0               &
    *(                                              &
    +2.D0*RE*XIX(I,J,K)*DJAC(I,J,K)/AP(I,J,K)                                                                           &
    *(DJAC(I+1,J,K)*XIX(I+1,J,K)*P(I+1,J,K)-2.D0*DJAC(I,J,K)*XIX(I,J,K)*P(I,J,K)+DJAC(I-1,J,K)*XIX(I-1,J,K)*P(I-1,J,K)) &
    -2.D0*RE*XIX(I-1,J,K)*DJAC(I-1,J,K)/AP(I-1,J,K)                                                                     &
    *(DJAC(I,J,K)*XIX(I,J,K)*P(I,J,K)-2.D0*DJAC(I-1,J,K)*XIX(I-1,J,K)*P(I-1,J,K)+DJAC(I-2,J,K)*XIX(I-2,J,K)*P(I-2,J,K)) &
    +2.D0*RE*XIY(I,J,K)*DJAC(I,J,K)/AP(I,J,K)                                                                           &
    *(DJAC(I+1,J,K)*XIY(I+1,J,K)*P(I+1,J,K)-2.D0*DJAC(I,J,K)*XIY(I,J,K)*P(I,J,K)+DJAC(I-1,J,K)*XIY(I-1,J,K)*P(I-1,J,K)) &
    -2.D0*RE*XIY(I-1,J,K)*DJAC(I-1,J,K)/AP(I-1,J,K)                                                                     &
    *(DJAC(I,J,K)*XIY(I,J,K)*P(I,J,K)-2.D0*DJAC(I-1,J,K)*XIY(I-1,J,K)*P(I-1,J,K)+DJAC(I-2,J,K)*XIY(I-2,J,K)*P(I-2,J,K)) &
    )
    ENDDO
    ENDDO
    ENDDO

!   VOLUME FLUX NORMAL TO ET-CONSTANT SURFACE
    DO K=2,NKM1
    DO I=SX1,EX1
    DO J=SY1,EY1-1
    VN(I,J,K) = VN(I,J,K) +0.25D0               &
    *(                                              &
    +2.D0*RE*ETX(I,J+1,K)*DJAC(I,J+1,K)/AP(I,J+1,K)                                                                     &
    *(DJAC(I,J+2,K)*ETX(I,J+2,K)*P(I,J+2,K)-2.D0*DJAC(I,J+1,K)*ETX(I,J+1,K)*P(I,J+1,K)+DJAC(I,J,K)*ETX(I,J,K)*P(I,J,K)) &
    -2.D0*RE*ETX(I,J,K)*DJAC(I,J,K)/AP(I,J,K)                                                                           &
    *(DJAC(I,J+1,K)*ETX(I,J+1,K)*P(I,J+1,K)-2.D0*DJAC(I,J,K)*ETX(I,J,K)*P(I,J,K)+DJAC(I,J-1,K)*ETX(I,J-1,K)*P(I,J-1,K)) &
    +2.D0*RE*ETY(I,J+1,K)*DJAC(I,J+1,K)/AP(I,J+1,K)                                                                     &
    *(DJAC(I,J+2,K)*ETY(I,J+2,K)*P(I,J+2,K)-2.D0*DJAC(I,J+1,K)*ETY(I,J+1,K)*P(I,J+1,K)+DJAC(I,J,K)*ETY(I,J,K)*P(I,J,K)) &
    -2.D0*RE*ETY(I,J,K)*DJAC(I,J,K)/AP(I,J,K)                                                                           &
    *(DJAC(I,J+1,K)*ETY(I,J+1,K)*P(I,J+1,K)-2.D0*DJAC(I,J,K)*ETY(I,J,K)*P(I,J,K)+DJAC(I,J-1,K)*ETY(I,J-1,K)*P(I,J-1,K)) &
    )
    ENDDO
    VN(I,EY1,K) = VN(I,EY1,K) +0.25D0           &
    *(                                              &
    +2.D0*RE*ETX(I,EY1+1,K)*DJAC(I,EY1+1,K)/AP(I,SY1,K)                                                                     &
    *(DJAC(I,SY1+1,K)*ETX(I,SY1+1,K)*P(I,SY1+1,K)-2.D0*DJAC(I,EY1+1,K)*ETX(I,EY1+1,K)*P(I,EY1+1,K)+DJAC(I,EY1,K)*ETX(I,EY1,K)*P(I,EY1,K)) &
    -2.D0*RE*ETX(I,EY1,K)*DJAC(I,EY1,K)/AP(I,EY1,K)                                                                           &
    *(DJAC(I,EY1+1,K)*ETX(I,EY1+1,K)*P(I,EY1+1,K)-2.D0*DJAC(I,EY1,K)*ETX(I,EY1,K)*P(I,EY1,K)+DJAC(I,EY1-1,K)*ETX(I,EY1-1,K)*P(I,EY1-1,K)) &
    +2.D0*RE*ETY(I,EY1+1,K)*DJAC(I,EY1+1,K)/AP(I,SY1,K)                                                                     &
    *(DJAC(I,SY1+1,K)*ETY(I,SY1+1,K)*P(I,SY1+1,K)-2.D0*DJAC(I,EY1+1,K)*ETY(I,EY1+1,K)*P(I,EY1+1,K)+DJAC(I,EY1,K)*ETY(I,EY1,K)*P(I,EY1,K)) &
    -2.D0*RE*ETY(I,EY1,K)*DJAC(I,EY1,K)/AP(I,EY1,K)                                                                           &
    *(DJAC(I,EY1+1,K)*ETY(I,EY1+1,K)*P(I,EY1+1,K)-2.D0*DJAC(I,EY1,K)*ETY(I,EY1,K)*P(I,EY1,K)+DJAC(I,EY1-1,K)*ETY(I,EY1-1,K)*P(I,EY1-1,K)) &
    )
    ENDDO
    ENDDO

    DO K=2,NKM1
    DO I=SX1,EX1
    DO J=SY1+1,EY1
    VS(I,J,K) = VS(I,J,K) +0.25D0               &
    *(                                              &
    +2.D0*RE*ETX(I,J,K)*DJAC(I,J,K)/AP(I,J,K)                                                                           &
    *(DJAC(I,J+1,K)*ETX(I,J+1,K)*P(I,J+1,K)-2.D0*DJAC(I,J,K)*ETX(I,J,K)*P(I,J,K)+DJAC(I,J-1,K)*ETX(I,J-1,K)*P(I,J-1,K)) &
    -2.D0*RE*ETX(I,J-1,K)*DJAC(I,J-1,K)/AP(I,J-1,K)                                                                     &
    *(DJAC(I,J,K)*ETX(I,J,K)*P(I,J,K)-2.D0*DJAC(I,J-1,K)*ETX(I,J-1,K)*P(I,J-1,K)+DJAC(I,J-2,K)*ETX(I,J-2,K)*P(I,J-2,K)) &
    +2.D0*RE*ETY(I,J,K)*DJAC(I,J,K)/AP(I,J,K)                                                                           &
    *(DJAC(I,J+1,K)*ETY(I,J+1,K)*P(I,J+1,K)-2.D0*DJAC(I,J,K)*ETY(I,J,K)*P(I,J,K)+DJAC(I,J-1,K)*ETY(I,J-1,K)*P(I,J-1,K)) &
    -2.D0*RE*ETY(I,J-1,K)*DJAC(I,J-1,K)/AP(I,J-1,K)                                                                     &
    *(DJAC(I,J,K)*ETY(I,J,K)*P(I,J,K)-2.D0*DJAC(I,J-1,K)*ETY(I,J-1,K)*P(I,J-1,K)+DJAC(I,J-2,K)*ETY(I,J-2,K)*P(I,J-2,K)) &
    )
    ENDDO
    VS(I,SY1,K) = VS(I,SY1,K) +0.25D0           &
    *(                                              &
    +2.D0*RE*ETX(I,SY1,K)*DJAC(I,SY1,K)/AP(I,SY1,K)                                                                           &
    *(DJAC(I,SY1+1,K)*ETX(I,SY1+1,K)*P(I,SY1+1,K)-2.D0*DJAC(I,SY1,K)*ETX(I,SY1,K)*P(I,SY1,K)+DJAC(I,SY1-1,K)*ETX(I,SY1-1,K)*P(I,SY1-1,K)) &
    -2.D0*RE*ETX(I,SY1-1,K)*DJAC(I,SY1-1,K)/AP(I,EY1,K)                                                                     &
    *(DJAC(I,SY1,K)*ETX(I,SY1,K)*P(I,SY1,K)-2.D0*DJAC(I,SY1-1,K)*ETX(I,SY1-1,K)*P(I,SY1-1,K)+DJAC(I,EY1-1,K)*ETX(I,EY1-1,K)*P(I,EY1-1,K)) &
    +2.D0*RE*ETY(I,SY1,K)*DJAC(I,SY1,K)/AP(I,SY1,K)                                                                           &
    *(DJAC(I,SY1+1,K)*ETY(I,SY1+1,K)*P(I,SY1+1,K)-2.D0*DJAC(I,SY1,K)*ETY(I,SY1,K)*P(I,SY1,K)+DJAC(I,SY1-1,K)*ETY(I,SY1-1,K)*P(I,SY1-1,K)) &
    -2.D0*RE*ETY(I,SY1-1,K)*DJAC(I,SY1-1,K)/AP(I,EY1,K)                                                                     &
    *(DJAC(I,SY1,K)*ETY(I,SY1,K)*P(I,SY1,K)-2.D0*DJAC(I,SY1-1,K)*ETY(I,SY1-1,K)*P(I,SY1-1,K)+DJAC(I,EY1-1,K)*ETY(I,EY1-1,K)*P(I,EY1-1,K)) &
    )
    ENDDO
    ENDDO
    ENDIF
!--------------------------------------------------------------------------------------------------------
!   BOUNDARY CONDITION FOR VOLUME FLUX
    UW(SX1,SY1:EY1,2:NKM1) = 0.D0
!--------------------------------------------------------------------------------------------------------
	END SUBROUTINE CONTRA_VEL