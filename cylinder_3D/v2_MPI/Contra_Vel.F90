	SUBROUTINE CONTRA_VEL

	USE COMDAT_SHARED 
      
	IMPLICIT NONE
    SAVE
	
	INCLUDE 'mpif.h'
	
!--------------------------------------------------------------------------------------------------------
!   CONTRA-VARIENT VELOCITY
    U1(IS:IE,JS:JE,KS:KE) = XIX(IS:IE,JS:JE,KS:KE)*U(IS:IE,JS:JE,KS:KE) + XIY(IS:IE,JS:JE,KS:KE)*V(IS:IE,JS:JE,KS:KE) + XIZ(IS:IE,JS:JE,KS:KE)*W(IS:IE,JS:JE,KS:KE)
    U2(IS:IE,JS:JE,KS:KE) = ETX(IS:IE,JS:JE,KS:KE)*U(IS:IE,JS:JE,KS:KE) + ETY(IS:IE,JS:JE,KS:KE)*V(IS:IE,JS:JE,KS:KE) + ETZ(IS:IE,JS:JE,KS:KE)*W(IS:IE,JS:JE,KS:KE)
    U3(IS:IE,JS:JE,KS:KE) = ZTX(IS:IE,JS:JE,KS:KE)*U(IS:IE,JS:JE,KS:KE) + ZTY(IS:IE,JS:JE,KS:KE)*V(IS:IE,JS:JE,KS:KE) + ZTZ(IS:IE,JS:JE,KS:KE)*W(IS:IE,JS:JE,KS:KE)
!--------------------------------------------------------------------------------------------------------

!	MPI : COLUMN-DIRECTION BLOCK DIVISION
	CALL MPI_X(U1)
	CALL MPI_X(U2)
	CALL MPI_X(U3)

!--------------------------------------------------------------------------------------------------------
!   VOLUME FLUX : CONTRA-VARIENT VELOCITY MULTIPLIED BY J^-1
    UE(IS1:IE1,JS1:JE1,KS1:KE1) =0.5D0*	&
		(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*U1(IS1:IE1,JS1:JE1,KS1:KE1)+DJAC(IS1+1:IE1+1,JS1:JE1,KS1:KE1)*U1(IS1+1:IE1+1,JS1:JE1,KS1:KE1))
    
    UW(IS1:IE1,JS1:JE1,KS1:KE1) =0.5D0*	&
		(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*U1(IS1:IE1,JS1:JE1,KS1:KE1)+DJAC(IS1-1:IE1-1,JS1:JE1,KS1:KE1)*U1(IS1-1:IE1-1,JS1:JE1,KS1:KE1))
    
    VN(IS1:IE1,JS1:JE1,KS1:KE1) =0.5D0*	&
		(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*U2(IS1:IE1,JS1:JE1,KS1:KE1)+DJAC(IS1:IE1,JS1+1:JE1+1,KS1:KE1)*U2(IS1:IE1,JS1+1:JE1+1,KS1:KE1))
    
    VS(IS1:IE1,JS1:JE1,KS1:KE1) =0.5D0*	&
		(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*U2(IS1:IE1,JS1:JE1,KS1:KE1)+DJAC(IS1:IE1,JS1-1:JE1-1,KS1:KE1)*U2(IS1:IE1,JS1-1:JE1-1,KS1:KE1))
    
    WT(IS1:IE1,JS1:JE1,KS1:KE1) =0.5D0*	&
		(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*U3(IS1:IE1,JS1:JE1,KS1:KE1)+DJAC(IS1:IE1,JS1:JE1,KS1+1:KE1+1)*U3(IS1:IE1,JS1:JE1,KS1+1:KE1+1))
    
    WB(IS1:IE1,JS1:JE1,KS1:KE1) =0.5D0*	&
		(DJAC(IS1:IE1,JS1:JE1,KS1:KE1)*U3(IS1:IE1,JS1:JE1,KS1:KE1)+DJAC(IS1:IE1,JS1:JE1,KS1-1:KE1-1)*U3(IS1:IE1,JS1:JE1,KS1-1:KE1-1))
!--------------------------------------------------------------------------------------------------------
!   CELL-FAVE VOLUME FLUX
    IF(FLAG) THEN
    
!   VOLUME FLUX NORMAL TO XI-CONSTANT SURFACE
    DO K=KS1,KE1
    DO J=JS1,JE1
    DO I=IS1,IE1-1
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

    DO K=KS1,KE1
    DO J=JS1,JE1
    DO I=IS1+1,IE1
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
    DO K=KS1,KE1
    DO J=JS1,JE1-1
    DO I=IS1,IE1
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
    VN(I,JE1,K) = VN(I,JE1,K) +0.25D0																					&
    *(																													&
    +2.D0*RE*ETX(I,JE1+1,K)/AP(I,JE1+1,K)																				&
    *(DJAC(I,JS1+1,K)*ETX(I,JS1+1,K)*P(I,JS1+1,K)-2.D0*DJAC(I,JE1+1,K)*ETX(I,JE1+1,K)*P(I,JE1+1,K)+DJAC(I,JE1,K)*ETX(I,JE1,K)*P(I,JE1,K))	&
    -2.D0*RE*ETX(I,JE1,K)/AP(I,JE1,K)																										&
    *(DJAC(I,JE1+1,K)*ETX(I,JE1+1,K)*P(I,JE1+1,K)-2.D0*DJAC(I,JE1,K)*ETX(I,JE1,K)*P(I,JE1,K)+DJAC(I,JE1-1,K)*ETX(I,JE1-1,K)*P(I,JE1-1,K))	&
    +2.D0*RE*ETY(I,JE1+1,K)/AP(I,JE1+1,K)																									&
    *(DJAC(I,JS1+1,K)*ETY(I,JS1+1,K)*P(I,JS1+1,K)-2.D0*DJAC(I,JE1+1,K)*ETY(I,JE1+1,K)*P(I,JE1+1,K)+DJAC(I,JE1,K)*ETY(I,JE1,K)*P(I,JE1,K))	&
    -2.D0*RE*ETY(I,JE1,K)/AP(I,JE1,K)																										&
    *(DJAC(I,JE1+1,K)*ETY(I,JE1+1,K)*P(I,JE1+1,K)-2.D0*DJAC(I,JE1,K)*ETY(I,JE1,K)*P(I,JE1,K)+DJAC(I,JE1-1,K)*ETY(I,JE1-1,K)*P(I,JE1-1,K))	&
    +2.D0*RE*ETZ(I,JE1+1,K)/AP(I,JE1+1,K)																									&
    *(DJAC(I,JS1+1,K)*ETZ(I,JS1+1,K)*P(I,JS1+1,K)-2.D0*DJAC(I,JE1+1,K)*ETZ(I,JE1+1,K)*P(I,JE1+1,K)+DJAC(I,JE1,K)*ETZ(I,JE1,K)*P(I,JE1,K))	&
    -2.D0*RE*ETZ(I,JE1,K)/AP(I,JE1,K)																										&
    *(DJAC(I,JE1+1,K)*ETZ(I,JE1+1,K)*P(I,JE1+1,K)-2.D0*DJAC(I,JE1,K)*ETZ(I,JE1,K)*P(I,JE1,K)+DJAC(I,JE1-1,K)*ETZ(I,JE1-1,K)*P(I,JE1-1,K))	&
    )
    ENDDO
    ENDDO

    DO K=KS1,KE1
    DO J=JS1+1,JE1
    DO I=IS1,IE1
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
    VS(I,JS1,K) = VS(I,JS1,K) +0.25D0																										&
    *(																																		&
    +2.D0*RE*ETX(I,JS1,K)/AP(I,JS1,K)																										&
    *(DJAC(I,JS1+1,K)*ETX(I,JS1+1,K)*P(I,JS1+1,K)-2.D0*DJAC(I,JS1,K)*ETX(I,JS1,K)*P(I,JS1,K)+DJAC(I,JS1-1,K)*ETX(I,JS1-1,K)*P(I,JS1-1,K))	&
    -2.D0*RE*ETX(I,JS1-1,K)/AP(I,JS1-1,K)																									&
    *(DJAC(I,JS1,K)*ETX(I,JS1,K)*P(I,JS1,K)-2.D0*DJAC(I,JS1-1,K)*ETX(I,JS1-1,K)*P(I,JS1-1,K)+DJAC(I,JE1-1,K)*ETX(I,JE1-1,K)*P(I,JE1-1,K))	&
    +2.D0*RE*ETY(I,JS1,K)/AP(I,JS1,K)																										&
    *(DJAC(I,JS1+1,K)*ETY(I,JS1+1,K)*P(I,JS1+1,K)-2.D0*DJAC(I,JS1,K)*ETY(I,JS1,K)*P(I,JS1,K)+DJAC(I,JS1-1,K)*ETY(I,JS1-1,K)*P(I,JS1-1,K))	&
    -2.D0*RE*ETY(I,JS1-1,K)/AP(I,JS1-1,K)																									&
    *(DJAC(I,JS1,K)*ETY(I,JS1,K)*P(I,JS1,K)-2.D0*DJAC(I,JS1-1,K)*ETY(I,JS1-1,K)*P(I,JS1-1,K)+DJAC(I,JE1-1,K)*ETY(I,JE1-1,K)*P(I,JE1-1,K))	&
    +2.D0*RE*ETZ(I,JS1,K)/AP(I,JS1,K)																										&
    *(DJAC(I,JS1+1,K)*ETZ(I,JS1+1,K)*P(I,JS1+1,K)-2.D0*DJAC(I,JS1,K)*ETZ(I,JS1,K)*P(I,JS1,K)+DJAC(I,JS1-1,K)*ETZ(I,JS1-1,K)*P(I,JS1-1,K))	&
    -2.D0*RE*ETZ(I,JS1-1,K)/AP(I,JS1-1,K)																									&
    *(DJAC(I,JS1,K)*ETZ(I,JS1,K)*P(I,JS1,K)-2.D0*DJAC(I,JS1-1,K)*ETZ(I,JS1-1,K)*P(I,JS1-1,K)+DJAC(I,JE1-1,K)*ETZ(I,JE1-1,K)*P(I,JE1-1,K))	&
    )
    ENDDO
    ENDDO
    
	DO K=KS1,KE1-1
    DO J=JS1,JE1
    DO I=IS1,IE1
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

    DO K=KS1+1,KE1
    DO J=JS1,JE1
    DO I=IS1+1,IE1
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
	IF(RANK.EQ.0) THEN
	UW(1,JS1:JE1,KS1:KE1) = 0.D0
	ENDIF
!--------------------------------------------------------------------------------------------------------
	END SUBROUTINE CONTRA_VEL