	SUBROUTINE CONTRA_VEL

	USE COMDAT_SHARED 
      
	IMPLICIT NONE
    SAVE
!--------------------------------------------------------------------------------------------------------
!   CONTRA-VARIENT VELOCITY
    U1 = XIX*U + XIY*V + XIZ*W
    U2 = ETX*U + ETY*V + ETZ*W
    U3 = ZTX*U + ZTY*V + ZTZ*W

!   VOLUME FLUX : CONTRA-VARIENT VELOCITY MULTIPLIED BY J^-1
    UE(SX1:EX1,SY1:EY1,SZ1:EZ1) =0.5D0*	&
		(DJAC(SX1:EX1,SY1:EY1,SZ1:EZ1)*U1(SX1:EX1,SY1:EY1,SZ1:EZ1)+DJAC(SX1+1:EX1+1,SY1:EY1,SZ1:EZ1)*U1(SX1+1:EX1+1,SY1:EY1,SZ1:EZ1))
    
    UW(SX1:EX1,SY1:EY1,SZ1:EZ1) =0.5D0*	&
		(DJAC(SX1:EX1,SY1:EY1,SZ1:EZ1)*U1(SX1:EX1,SY1:EY1,SZ1:EZ1)+DJAC(SX1-1:EX1-1,SY1:EY1,SZ1:EZ1)*U1(SX1-1:EX1-1,SY1:EY1,SZ1:EZ1))
    
    VN(SX1:EX1,SY1:EY1,SZ1:EZ1) =0.5D0*	&
		(DJAC(SX1:EX1,SY1:EY1,SZ1:EZ1)*U2(SX1:EX1,SY1:EY1,SZ1:EZ1)+DJAC(SX1:EX1,SY1+1:EY1+1,SZ1:EZ1)*U2(SX1:EX1,SY1+1:EY1+1,SZ1:EZ1))
    
    VS(SX1:EX1,SY1:EY1,SZ1:EZ1) =0.5D0*	&
		(DJAC(SX1:EX1,SY1:EY1,SZ1:EZ1)*U2(SX1:EX1,SY1:EY1,SZ1:EZ1)+DJAC(SX1:EX1,SY1-1:EY1-1,SZ1:EZ1)*U2(SX1:EX1,SY1-1:EY1-1,SZ1:EZ1))
    
    WT(SX1:EX1,SY1:EY1,SZ1:EZ1) =0.5D0*	&
		(DJAC(SX1:EX1,SY1:EY1,SZ1:EZ1)*U3(SX1:EX1,SY1:EY1,SZ1:EZ1)+DJAC(SX1:EX1,SY1:EY1,SZ1+1:EZ1+1)*U3(SX1:EX1,SY1:EY1,SZ1+1:EZ1+1))
    
    WB(SX1:EX1,SY1:EY1,SZ1:EZ1) =0.5D0*	&
		(DJAC(SX1:EX1,SY1:EY1,SZ1:EZ1)*U3(SX1:EX1,SY1:EY1,SZ1:EZ1)+DJAC(SX1:EX1,SY1:EY1,SZ1-1:EZ1-1)*U3(SX1:EX1,SY1:EY1,SZ1-1:EZ1-1))
!--------------------------------------------------------------------------------------------------------
!   CELL-FAVE VOLUME FLUX
    IF(FLAG) THEN
    AP(SX,SY1:EY1,SZ1:EZ1) = AP(SX1,SY1:EY1,SZ1:EZ1)
    AP(EX,SY1:EY1,SZ1:EZ1) = AP(EX1,SY1:EY1,SZ1:EZ1)
    AP(SX1:EX1,SY,SZ1:EZ1) = AP(SX1:EX1,EY1,SZ1:EZ1)  !AP(SX1:EX1,SY,SZ1:EZ1) = AP(SX1:EX1,SY1,SZ1:EZ1)
    AP(SX1:EX1,EY,SZ1:EZ1) = AP(SX1:EX1,SY1,SZ1:EZ1)  !AP(SX1:EX1,EY,SZ1:EZ1) = AP(SX1:EX1,EY1,SZ1:EZ1)
    AP(SX1:EX1,SY1:EY1,SZ) = AP(SX1:EX1,SY1:EY1,SZ1)
    AP(SX1:EX1,SY1:EY1,EZ) = AP(SX1:EX1,SY1:EY1,EZ1)
    
!   VOLUME FLUX NORMAL TO XI-CONSTANT SURFACE
    DO K=SZ1,EZ1
    DO J=SY1,EY1
    DO I=SX1,EX1-1
    UE(I,J,K) = UE(I,J,K) +0.25D0																						&
    *(																													&
    +2.D0*RE*XIX(I+1,J,K)/AP(I+1,J,K)																					&
    *(DJAC(I+2,J,K)*XIX(I+2,J,K)*P(I+2,J,K)-2.D0*DJAC(I+1,J,K)*XIX(I+1,J,K)*P(I+1,J,K)+DJAC(I,J,K)*XIX(I,J,K)*P(I,J,K)) &
    -2.D0*RE*XIX(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I+1,J,K)*XIX(I+1,J,K)*P(I+1,J,K)-2.D0*DJAC(I,J,K)*XIX(I,J,K)*P(I,J,K)+DJAC(I-1,J,K)*XIX(I-1,J,K)*P(I-1,J,K)) &
    +2.D0*RE*XIY(I+1,J,K)/AP(I+1,J,K)																					&
    *(DJAC(I+2,J,K)*XIY(I+2,J,K)*P(I+2,J,K)-2.D0*DJAC(I+1,J,K)*XIY(I+1,J,K)*P(I+1,J,K)+DJAC(I,J,K)*XIY(I,J,K)*P(I,J,K)) &
    -2.D0*RE*XIY(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I+1,J,K)*XIY(I+1,J,K)*P(I+1,J,K)-2.D0*DJAC(I,J,K)*XIY(I,J,K)*P(I,J,K)+DJAC(I-1,J,K)*XIY(I-1,J,K)*P(I-1,J,K)) &
	+2.D0*RE*XIZ(I+1,J,K)/AP(I+1,J,K)																					&
    *(DJAC(I+2,J,K)*XIZ(I+2,J,K)*P(I+2,J,K)-2.D0*DJAC(I+1,J,K)*XIZ(I+1,J,K)*P(I+1,J,K)+DJAC(I,J,K)*XIZ(I,J,K)*P(I,J,K)) &
    -2.D0*RE*XIZ(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I+1,J,K)*XIZ(I+1,J,K)*P(I+1,J,K)-2.D0*DJAC(I,J,K)*XIZ(I,J,K)*P(I,J,K)+DJAC(I-1,J,K)*XIZ(I-1,J,K)*P(I-1,J,K)) &
    )
    ENDDO
    ENDDO
    ENDDO

    DO K=SZ1,EZ1
    DO J=SY1,EY1
    DO I=SX1+1,EX1
    UW(I,J,K) = UW(I,J,K) +0.25D0																						&
    *(																													&
    +2.D0*RE*XIX(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I+1,J,K)*XIX(I+1,J,K)*P(I+1,J,K)-2.D0*DJAC(I,J,K)*XIX(I,J,K)*P(I,J,K)+DJAC(I-1,J,K)*XIX(I-1,J,K)*P(I-1,J,K)) &
    -2.D0*RE*XIX(I-1,J,K)/AP(I-1,J,K)																					&
    *(DJAC(I,J,K)*XIX(I,J,K)*P(I,J,K)-2.D0*DJAC(I-1,J,K)*XIX(I-1,J,K)*P(I-1,J,K)+DJAC(I-2,J,K)*XIX(I-2,J,K)*P(I-2,J,K)) &
    +2.D0*RE*XIY(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I+1,J,K)*XIY(I+1,J,K)*P(I+1,J,K)-2.D0*DJAC(I,J,K)*XIY(I,J,K)*P(I,J,K)+DJAC(I-1,J,K)*XIY(I-1,J,K)*P(I-1,J,K)) &
    -2.D0*RE*XIY(I-1,J,K)/AP(I-1,J,K)																					&
    *(DJAC(I,J,K)*XIY(I,J,K)*P(I,J,K)-2.D0*DJAC(I-1,J,K)*XIY(I-1,J,K)*P(I-1,J,K)+DJAC(I-2,J,K)*XIY(I-2,J,K)*P(I-2,J,K)) &
    +2.D0*RE*XIZ(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I+1,J,K)*XIZ(I+1,J,K)*P(I+1,J,K)-2.D0*DJAC(I,J,K)*XIZ(I,J,K)*P(I,J,K)+DJAC(I-1,J,K)*XIZ(I-1,J,K)*P(I-1,J,K)) &
    -2.D0*RE*XIZ(I-1,J,K)/AP(I-1,J,K)																					&
    *(DJAC(I,J,K)*XIZ(I,J,K)*P(I,J,K)-2.D0*DJAC(I-1,J,K)*XIZ(I-1,J,K)*P(I-1,J,K)+DJAC(I-2,J,K)*XIZ(I-2,J,K)*P(I-2,J,K)) &
    )
    ENDDO
    ENDDO
    ENDDO

!   VOLUME FLUX NORMAL TO ET-CONSTANT SURFACE
    DO K=SZ1,EZ1
    DO J=SY1,EY1-1
    DO I=SX1,EX1
    VN(I,J,K) = VN(I,J,K) +0.25D0																						&
    *(																													&
    +2.D0*RE*ETX(I,J+1,K)/AP(I,J+1,K)																					&
    *(DJAC(I,J+2,K)*ETX(I,J+2,K)*P(I,J+2,K)-2.D0*DJAC(I,J+1,K)*ETX(I,J+1,K)*P(I,J+1,K)+DJAC(I,J,K)*ETX(I,J,K)*P(I,J,K)) &
    -2.D0*RE*ETX(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I,J+1,K)*ETX(I,J+1,K)*P(I,J+1,K)-2.D0*DJAC(I,J,K)*ETX(I,J,K)*P(I,J,K)+DJAC(I,J-1,K)*ETX(I,J-1,K)*P(I,J-1,K)) &
    +2.D0*RE*ETY(I,J+1,K)/AP(I,J+1,K)																					&
    *(DJAC(I,J+2,K)*ETY(I,J+2,K)*P(I,J+2,K)-2.D0*DJAC(I,J+1,K)*ETY(I,J+1,K)*P(I,J+1,K)+DJAC(I,J,K)*ETY(I,J,K)*P(I,J,K)) &
    -2.D0*RE*ETY(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I,J+1,K)*ETY(I,J+1,K)*P(I,J+1,K)-2.D0*DJAC(I,J,K)*ETY(I,J,K)*P(I,J,K)+DJAC(I,J-1,K)*ETY(I,J-1,K)*P(I,J-1,K)) &
	+2.D0*RE*ETZ(I,J+1,K)/AP(I,J+1,K)																					&
    *(DJAC(I,J+2,K)*ETZ(I,J+2,K)*P(I,J+2,K)-2.D0*DJAC(I,J+1,K)*ETZ(I,J+1,K)*P(I,J+1,K)+DJAC(I,J,K)*ETZ(I,J,K)*P(I,J,K)) &
    -2.D0*RE*ETZ(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I,J+1,K)*ETZ(I,J+1,K)*P(I,J+1,K)-2.D0*DJAC(I,J,K)*ETZ(I,J,K)*P(I,J,K)+DJAC(I,J-1,K)*ETZ(I,J-1,K)*P(I,J-1,K)) &
    )
    ENDDO
    VN(I,EY1,K) = VN(I,EY1,K) +0.25D0																					&
    *(																													&
    +2.D0*RE*ETX(I,EY1+1,K)/AP(I,EY1+1,K)																				&
    *(DJAC(I,SY1+1,K)*ETX(I,SY1+1,K)*P(I,SY1+1,K)-2.D0*DJAC(I,EY1+1,K)*ETX(I,EY1+1,K)*P(I,EY1+1,K)+DJAC(I,EY1,K)*ETX(I,EY1,K)*P(I,EY1,K))	&
    -2.D0*RE*ETX(I,EY1,K)/AP(I,EY1,K)																										&
    *(DJAC(I,EY1+1,K)*ETX(I,EY1+1,K)*P(I,EY1+1,K)-2.D0*DJAC(I,EY1,K)*ETX(I,EY1,K)*P(I,EY1,K)+DJAC(I,EY1-1,K)*ETX(I,EY1-1,K)*P(I,EY1-1,K))	&
    +2.D0*RE*ETY(I,EY1+1,K)/AP(I,EY1+1,K)																									&
    *(DJAC(I,SY1+1,K)*ETY(I,SY1+1,K)*P(I,SY1+1,K)-2.D0*DJAC(I,EY1+1,K)*ETY(I,EY1+1,K)*P(I,EY1+1,K)+DJAC(I,EY1,K)*ETY(I,EY1,K)*P(I,EY1,K))	&
    -2.D0*RE*ETY(I,EY1,K)/AP(I,EY1,K)																										&
    *(DJAC(I,EY1+1,K)*ETY(I,EY1+1,K)*P(I,EY1+1,K)-2.D0*DJAC(I,EY1,K)*ETY(I,EY1,K)*P(I,EY1,K)+DJAC(I,EY1-1,K)*ETY(I,EY1-1,K)*P(I,EY1-1,K))	&
    +2.D0*RE*ETZ(I,EY1+1,K)/AP(I,EY1+1,K)																									&
    *(DJAC(I,SY1+1,K)*ETZ(I,SY1+1,K)*P(I,SY1+1,K)-2.D0*DJAC(I,EY1+1,K)*ETZ(I,EY1+1,K)*P(I,EY1+1,K)+DJAC(I,EY1,K)*ETZ(I,EY1,K)*P(I,EY1,K))	&
    -2.D0*RE*ETZ(I,EY1,K)/AP(I,EY1,K)																										&
    *(DJAC(I,EY1+1,K)*ETZ(I,EY1+1,K)*P(I,EY1+1,K)-2.D0*DJAC(I,EY1,K)*ETZ(I,EY1,K)*P(I,EY1,K)+DJAC(I,EY1-1,K)*ETZ(I,EY1-1,K)*P(I,EY1-1,K))	&
    )
    ENDDO
    ENDDO

    DO K=SZ1,EZ1
    DO J=SY1+1,EY1
    DO I=SX1,EX1
    VS(I,J,K) = VS(I,J,K) +0.25D0																						&
    *(																													&
    +2.D0*RE*ETX(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I,J+1,K)*ETX(I,J+1,K)*P(I,J+1,K)-2.D0*DJAC(I,J,K)*ETX(I,J,K)*P(I,J,K)+DJAC(I,J-1,K)*ETX(I,J-1,K)*P(I,J-1,K)) &
    -2.D0*RE*ETX(I,J-1,K)/AP(I,J-1,K)																					&
    *(DJAC(I,J,K)*ETX(I,J,K)*P(I,J,K)-2.D0*DJAC(I,J-1,K)*ETX(I,J-1,K)*P(I,J-1,K)+DJAC(I,J-2,K)*ETX(I,J-2,K)*P(I,J-2,K)) &
    +2.D0*RE*ETY(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I,J+1,K)*ETY(I,J+1,K)*P(I,J+1,K)-2.D0*DJAC(I,J,K)*ETY(I,J,K)*P(I,J,K)+DJAC(I,J-1,K)*ETY(I,J-1,K)*P(I,J-1,K)) &
    -2.D0*RE*ETY(I,J-1,K)/AP(I,J-1,K)																					&
    *(DJAC(I,J,K)*ETY(I,J,K)*P(I,J,K)-2.D0*DJAC(I,J-1,K)*ETY(I,J-1,K)*P(I,J-1,K)+DJAC(I,J-2,K)*ETY(I,J-2,K)*P(I,J-2,K)) &
    +2.D0*RE*ETZ(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I,J+1,K)*ETZ(I,J+1,K)*P(I,J+1,K)-2.D0*DJAC(I,J,K)*ETZ(I,J,K)*P(I,J,K)+DJAC(I,J-1,K)*ETZ(I,J-1,K)*P(I,J-1,K)) &
    -2.D0*RE*ETZ(I,J-1,K)/AP(I,J-1,K)																					&
    *(DJAC(I,J,K)*ETZ(I,J,K)*P(I,J,K)-2.D0*DJAC(I,J-1,K)*ETZ(I,J-1,K)*P(I,J-1,K)+DJAC(I,J-2,K)*ETZ(I,J-2,K)*P(I,J-2,K)) &
    )
    ENDDO
    VS(I,SY1,K) = VS(I,SY1,K) +0.25D0																										&
    *(																																		&
    +2.D0*RE*ETX(I,SY1,K)/AP(I,SY1,K)																										&
    *(DJAC(I,SY1+1,K)*ETX(I,SY1+1,K)*P(I,SY1+1,K)-2.D0*DJAC(I,SY1,K)*ETX(I,SY1,K)*P(I,SY1,K)+DJAC(I,SY1-1,K)*ETX(I,SY1-1,K)*P(I,SY1-1,K))	&
    -2.D0*RE*ETX(I,SY1-1,K)/AP(I,SY1-1,K)																									&
    *(DJAC(I,SY1,K)*ETX(I,SY1,K)*P(I,SY1,K)-2.D0*DJAC(I,SY1-1,K)*ETX(I,SY1-1,K)*P(I,SY1-1,K)+DJAC(I,EY1-1,K)*ETX(I,EY1-1,K)*P(I,EY1-1,K))	&
    +2.D0*RE*ETY(I,SY1,K)/AP(I,SY1,K)																										&
    *(DJAC(I,SY1+1,K)*ETY(I,SY1+1,K)*P(I,SY1+1,K)-2.D0*DJAC(I,SY1,K)*ETY(I,SY1,K)*P(I,SY1,K)+DJAC(I,SY1-1,K)*ETY(I,SY1-1,K)*P(I,SY1-1,K))	&
    -2.D0*RE*ETY(I,SY1-1,K)/AP(I,SY1-1,K)																									&
    *(DJAC(I,SY1,K)*ETY(I,SY1,K)*P(I,SY1,K)-2.D0*DJAC(I,SY1-1,K)*ETY(I,SY1-1,K)*P(I,SY1-1,K)+DJAC(I,EY1-1,K)*ETY(I,EY1-1,K)*P(I,EY1-1,K))	&
    +2.D0*RE*ETZ(I,SY1,K)/AP(I,SY1,K)																										&
    *(DJAC(I,SY1+1,K)*ETZ(I,SY1+1,K)*P(I,SY1+1,K)-2.D0*DJAC(I,SY1,K)*ETZ(I,SY1,K)*P(I,SY1,K)+DJAC(I,SY1-1,K)*ETZ(I,SY1-1,K)*P(I,SY1-1,K))	&
    -2.D0*RE*ETZ(I,SY1-1,K)/AP(I,SY1-1,K)																									&
    *(DJAC(I,SY1,K)*ETZ(I,SY1,K)*P(I,SY1,K)-2.D0*DJAC(I,SY1-1,K)*ETZ(I,SY1-1,K)*P(I,SY1-1,K)+DJAC(I,EY1-1,K)*ETZ(I,EY1-1,K)*P(I,EY1-1,K))	&
    )
    ENDDO
    ENDDO
    
	DO K=SZ1,EZ1-1
    DO J=SY1,EY1
    DO I=SX1,EX1
    WT(I,J,K) = WT(I,J,K) +0.25D0																						&
    *(																													&
    +2.D0*RE*ZTX(I,J,K+1)/AP(I,J,K+1)																					&
    *(DJAC(I,J,K+2)*ZTX(I,J,K+2)*P(I,J,K+2)-2.D0*DJAC(I,J,K+1)*ZTX(I,J,K+1)*P(I,J,K+1)+DJAC(I,J,K)*ZTX(I,J,K)*P(I,J,K)) &
    -2.D0*RE*ZTX(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I,J,K+1)*ZTX(I,J,K+1)*P(I,J,K+1)-2.D0*DJAC(I,J,K)*ZTX(I,J,K)*P(I,J,K)+DJAC(I,J,K-1)*ZTX(I,J,K-1)*P(I,J,K-1)) &
    +2.D0*RE*ZTY(I,J,K+1)/AP(I,J,K+1)																					&
    *(DJAC(I,J,K+2)*ZTY(I,J,K+2)*P(I,J,K+2)-2.D0*DJAC(I,J,K+1)*ZTY(I,J,K+1)*P(I,J,K+1)+DJAC(I,J,K)*ZTY(I,J,K)*P(I,J,K)) &
    -2.D0*RE*ZTY(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I,J,K+1)*ZTY(I,J,K+1)*P(I,J,K+1)-2.D0*DJAC(I,J,K)*ZTY(I,J,K)*P(I,J,K)+DJAC(I,J,K-1)*ZTY(I,J,K-1)*P(I,J,K-1)) &
    +2.D0*RE*ZTZ(I,J,K+1)/AP(I,J,K+1)																					&
    *(DJAC(I,J,K+2)*ZTZ(I,J,K+2)*P(I,J,K+2)-2.D0*DJAC(I,J,K+1)*ZTZ(I,J,K+1)*P(I,J,K+1)+DJAC(I,J,K)*ZTZ(I,J,K)*P(I,J,K)) &
    -2.D0*RE*ZTZ(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I,J,K+1)*ZTZ(I,J,K+1)*P(I,J,K+1)-2.D0*DJAC(I,J,K)*ZTZ(I,J,K)*P(I,J,K)+DJAC(I,J,K-1)*ZTZ(I,J,K-1)*P(I,J,K-1)) &
    )
    ENDDO
    ENDDO
    ENDDO

    DO K=SZ1+1,EZ1
    DO J=SY1,EY1
    DO I=SX1,EX1
    WB(I,J,K) = WB(I,J,K) +0.25D0																						&
    *(																													&
    +2.D0*RE*ZTX(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I,J,K+1)*ZTX(I,J,K+1)*P(I,J,K+1)-2.D0*DJAC(I,J,K)*ZTX(I,J,K)*P(I,J,K)+DJAC(I,J,K-1)*ZTX(I,J,K-1)*P(I,J,K-1)) &
    -2.D0*RE*ZTX(I,J,K-1)/AP(I,J,K-1)																					&
    *(DJAC(I,J,K)*ZTX(I,J,K)*P(I,J,K)-2.D0*DJAC(I,J,K-1)*ZTX(I,J,K-1)*P(I,J,K-1)+DJAC(I,J,K-2)*ZTX(I,J,K-2)*P(I,J,K-2)) &
    +2.D0*RE*ZTY(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I,J,K+1)*ZTY(I,J,K+1)*P(I,J,K+1)-2.D0*DJAC(I,J,K)*ZTY(I,J,K)*P(I,J,K)+DJAC(I,J,K-1)*ZTY(I,J,K-1)*P(I,J,K-1)) &
    -2.D0*RE*ZTY(I,J,K-1)/AP(I,J,K-1)																					&
    *(DJAC(I,J,K)*ZTY(I,J,K)*P(I,J,K)-2.D0*DJAC(I,J,K-1)*ZTY(I,J,K-1)*P(I,J,K-1)+DJAC(I,J,K-2)*ZTY(I,J,K-2)*P(I,J,K-2)) &
    +2.D0*RE*ZTZ(I,J,K)/AP(I,J,K)																						&
    *(DJAC(I,J,K+1)*ZTZ(I,J,K+1)*P(I,J,K+1)-2.D0*DJAC(I,J,K)*ZTZ(I,J,K)*P(I,J,K)+DJAC(I,J,K-1)*ZTZ(I,J,K-1)*P(I,J,K-1)) &
    -2.D0*RE*ZTZ(I,J,K-1)/AP(I,J,K-1)																					&
    *(DJAC(I,J,K)*ZTZ(I,J,K)*P(I,J,K)-2.D0*DJAC(I,J,K-1)*ZTZ(I,J,K-1)*P(I,J,K-1)+DJAC(I,J,K-2)*ZTZ(I,J,K-2)*P(I,J,K-2)) &
    )
    ENDDO
    ENDDO
    ENDDO
    
    ENDIF
!--------------------------------------------------------------------------------------------------------
!   BOUNDARY CONDITION FOR VOLUME FLUX
    UW(SX1,SY1:EY1,SZ1:EZ1) = 0.D0
!--------------------------------------------------------------------------------------------------------
	END SUBROUTINE CONTRA_VEL