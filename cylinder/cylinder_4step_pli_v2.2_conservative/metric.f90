	SUBROUTINE METRIC
	
	USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE

	ALLOCATE(XXI(0:NI,0:NJ,2:2),STAT=ALLOC_ERR)

	ALLOCATE(YXI(0:NI,0:NJ,2:2),STAT=ALLOC_ERR)

	ALLOCATE(XET(0:NI,0:NJ,2:2),STAT=ALLOC_ERR)

	ALLOCATE(YET(0:NI,0:NJ,2:2),STAT=ALLOC_ERR)
	
!----------------------------------------------------------
!   CELL-CENTER VALUE
	DO K=2,NKM1
	DO J=SY1,ey1
	DO I=sx1,ex1
	    XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
	    XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0

	    YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
	    YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
	ENDDO
	ENDDO
	ENDDO
!----------------------------------------------------------

!----- SOUTH FACE (I,0,K)----------
    DO K=2,NKM1
    DO I=SX1,EX1
        XXI(I,0,K)=(X(I+1,0,K)-X(I-1,0,K))/2.0D0
        XET(I,0,K)=(X(I,NJ,K)-X(I,NJ-2,K))/2.0D0 !(-3.0D0*X(I,0,K)+4.0D0*X(I,1,K)-X(I,2,K))/2.0D0
	
        YXI(I,0,K)=(Y(I+1,0,K)-Y(I-1,0,K))/2.0D0
        YET(I,0,K)=(Y(I,NJ,K)-Y(I,NJ-2,K))/2.0D0 !(-3.0D0*Y(I,0,K)+4.0D0*Y(I,1,K)-Y(I,2,K))/2.0D0
	ENDDO
	ENDDO  
!------ NORTH FACE (I,NJ,K) --------
	DO K=2,NKM1
	DO I=SX1,EX1
        XXI(I,NJ,K)=(X(I+1,NJ,K)-X(I-1,NJ,K))/2.0D0
        XET(I,NJ,K)=(X(I,2,K)-X(I,0,K))/2.0D0 !(3.0D0*X(I,NJ,K)-4.0D0*X(I,NJM1,K)+X(I,NJM2,K))/2.0D0
	
        YXI(I,NJ,K)=(Y(I+1,NJ,K)-Y(I-1,NJ,K))/2.0D0
        YET(I,NJ,K)=(Y(I,2,K)-Y(I,0,K))/2.0D0 !(3.0D0*Y(I,NJ,K)-4.0D0*Y(I,NJM1,K)+Y(I,NJM2,K))/2.0D0
	ENDDO
	ENDDO 
!----- WEST FACE (0,J,K) --------- 
	DO K=2,NKM1
	DO J=SY1,EY1
        XXI(0,J,K)=(-3.0D0*X(0,J,K)+4.0D0*X(1,J,K)-X(2,J,K))/2.0D0
        XET(0,J,K)=(X(0,J+1,K)-X(0,J-1,K))/2.0D0
	
        YXI(0,J,K)=(-3.0D0*Y(0,J,K)+4.0D0*Y(1,J,K)-Y(2,J,K))/2.0D0
        YET(0,J,K)=(Y(0,J+1,K)-Y(0,J-1,K))/2.0D0
    ENDDO
	ENDDO
!----- EAST FACE (NI,J,K)------------
    DO K=2,NKM1
	DO J=SY1,EY1
        XXI(NI,J,K)=(3.0D0*X(NI,J,K)-4.0D0*X(NIM1,J,K)+X(NIM2,J,K))/2.0D0
        XET(NI,J,K)=(X(NI,J+1,K)-X(NI,J-1,K))/2.0D0
	
        YXI(NI,J,K)=(3.0D0*Y(NI,J,K)-4.0D0*Y(NIM1,J,K)+Y(NIM2,J,K))/2.0D0
        YET(NI,J,K)=(Y(NI,J+1,K)-Y(NI,J-1,K))/2.0D0
	ENDDO
	ENDDO
!-----LINE(0,0,K)
	I=0
	J=0
	DO K=2,NKM1	
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0                                                                                                                                                  
        XET(I,J,K)=(X(I,NJ,K)-X(I,NJ-2,K))/2.0D0 !(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0

        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,NJ,K)-Y(I,NJ-2,K))/2.0D0 !(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
	ENDDO
!-----LINE(0,NJ,K)
	I=0
	J=NJ
	DO K=2,NKM1
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(X(I,2,K)-X(I,0,K))/2.0D0 !( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0

        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,2,K)-Y(I,0,K))/2.0D0 !( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
	ENDDO
!-----LINE(NI,0,K)
	I=NI
	J=0
	DO K=2,NKM1
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(X(I,NJ,K)-X(I,NJ-2,K))/2.0D0 !(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0

        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,NJ,K)-Y(I,NJ-2,K))/2.0D0 !(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
	ENDDO
!-----LINE(NI,NJ,K)
	I=NI
	J=NJ
	DO K=2,NKM1
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(X(I,2,K)-X(I,0,K))/2.0D0 !( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0

        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,2,K)-Y(I,0,K))/2.0D0 !( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
	ENDDO
!----------------------------------------------------------
!   DETERMINENT OF INVERSE JACOBIAN MATRIX
	DO K=2,NKM1		!DO K=1,NK
	DO J=SY,EY
	DO I=SX,EX
         DJAC(I,J,K)	&
         =  XXI(I,J,K)*YET(I,J,K)	&
           -YXI(I,J,K)*XET(I,J,K)	!&
        !  check the jacobian error
        IF (ABS(DJAC(I,J,K)).LT.1.0E-15) WRITE(6,*) I,J,K,DJAC(I,J,K) 
	ENDDO
	ENDDO
	ENDDO
!----------------------------------------------------------
 	DO K=2,NKM1		!DO K=1,NK
    DO J=SY,EY
	DO I=SX,EX
          XIX(I,J,K) =  YET(I,J,K)/DJAC(I,J,K)
          XIY(I,J,K) = -XET(I,J,K)/DJAC(I,J,K)

          ETX(I,J,K) = -YXI(I,J,K)/DJAC(I,J,K)
          ETY(I,J,K) =  XXI(I,J,K)/DJAC(I,J,K)
	ENDDO
	ENDDO
	ENDDO
	
	DEALLOCATE(XXI)

	DEALLOCATE(YXI)

	DEALLOCATE(XET)

	DEALLOCATE(YET)
	
!--- Q11,Q12,Q21,Q22 calculation ----------------------------------------------
 	DO K=2,NKM1		!	DO K=1,NK
	DO J=SY,EY
	DO I=SX,EX
          Q11(I,J,K)= XIX(I,J,K)**2.0D0     + XIY(I,J,K)**2.0D0    
          Q12(I,J,K)= XIX(I,J,K)*ETX(I,J,K) + XIY(I,J,K)*ETY(I,J,K) 

          Q21(I,J,K)= Q12(I,J,K)
          Q22(I,J,K)= ETX(I,J,K)**2.0D0     + ETY(I,J,K)**2.0D0
	ENDDO
	ENDDO
	ENDDO

!--- DE,DW,DN,DS calculation --------------------------------------------------
	DO K=2,NKM1
	DO J=SY1,EY1
 	DO I=SX1,EX1
 	    DE(I,J,K)=0.5D0*(DJAC(I,J,K)*Q11(I,J,K)+DJAC(I+1,J,K)*Q11(I+1,J,K))
	    DW(I,J,K)=0.5D0*(DJAC(I,J,K)*Q11(I,J,K)+DJAC(I-1,J,K)*Q11(I-1,J,K))
	    DN(I,J,K)=0.5D0*(DJAC(I,J,K)*Q22(I,J,K)+DJAC(I,J+1,K)*Q22(I,J+1,K))
	    DS(I,J,K)=0.5D0*(DJAC(I,J,K)*Q22(I,J,K)+DJAC(I,J-1,K)*Q22(I,J-1,K))
    ENDDO
	ENDDO
	ENDDO

	END SUBROUTINE METRIC