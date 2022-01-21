    SUBROUTINE ITER_PLI
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!---------------------------------------------------------------------------
!   TRI-DIAGONAL COMPONENTS OF LINEAR MATRIX EQUATION
!   AT THE CELL CENTER : PLI
!   THE NEIGHBOR CELLS : E,W,N,S
!   THE EXPLICIT TERMS : Su
    AP = DE+DW+DN+DS
    AE = DE
	AW = DW
	AN = DN
	AS = DS
	SU = 0
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
!   TDMA solver iteration
    NITER = 0 
804	NITER = NITER+1
    CALL DFI(PLI,DFIP)
    CALL DFE(PLI,DFEP) ! DIAGONAL DIF OF PLI : TREATED AS EXPLICIT

!--- SUMMATION OF EXPLICIT COMPONENTS --------------------------------------
    SU(SX1:EX1,SY1:EY1,2:NKM1) = -1.0D0/DTIME/DJAC(SX1:EX1,SY1:EY1,2:NKM1) *(		&
                        +UE(SX1:EX1,SY1:EY1,2:NKM1) -UW(SX1:EX1,SY1:EY1,2:NKM1)		&
                        +VN(SX1:EX1,SY1:EY1,2:NKM1) -VS(SX1:EX1,SY1:EY1,2:NKM1))	&
	                    +DFIP(SX1:EX1,SY1:EY1,2:NKM1) +DFEP(SX1:EX1,SY1:EY1,2:NKM1)

!--- LINE GAUSS-GEIDEL ITERATION METHOD ------------------------------------
    CALL BC_PLI
    
    DO K=2,NKM1
    DO J=SY1,EY1
    DO I=SX1,EX1
        A(I)=-AW(I,J,K)             ! A(1) SHOULD BE 0(ZERO)
        D(I)=+AP(I,J,K)-A(I)*B(I-1) ! A(1) SHOULD BE 0(ZERO)
        B(I)=-AE(I,J,K)/D(I)        ! B(NI2-1) SHOULD BE 0(ZERO)
        C(I)=AN(I,J,K)*PLI(I,J+1,K)+AS(I,J,K)*PLI(I,J-1,K)  &
!           +AT(I,J,K)*PLI(I,J,K+1)+AB(I,J,K)*PLI(I,J,K-1)	&
            +SU(I,J,K)
        PLI(I,J,K)=(C(I)-A(I)*PLI(I-1,J,K))/D(I) ! A(1) SHOULD BE 0(ZERO)
    ENDDO
    DO I=EX1,SX1,-1
        PLI(I,J,K)=PLI(I,J,K)-B(I)*PLI(I+1,J,K)
    ENDDO
    ENDDO
    ENDDO
    
    CALL BC_PLI
        
!---- PLI CONVERGE TEST : RESIDUAL OF PLI ----------------------------------
	RESIP = 0.0D0
    DUM = 0.D0
	DO K=2,NKM1
	DO J=SY1,EY1
	DO I=SX1,EX1
	    RESOR = AE(I,J,K)*PLI(I+1,J  ,K) + AW(I,J,K)*PLI(I-1,J  ,K) +   &
				AN(I,J,K)*PLI(I  ,J+1,K) + AS(I,J,K)*PLI(I  ,J-1,K) +	&
				SU(I,J,K) - AP(I,J,K)*PLI(I,J,K) 
		RESIP = RESIP + ABS(RESOR)
		DUM = DUM + ABS(AP(I,J,K)*PLI(I,J,K) )
	ENDDO
	ENDDO
	ENDDO
	RESIP = RESIP / DUM ! SCALED RESIDUAL
	IF(RESIP.GT.RESIMAX .AND. NITER.LT.NSWPP) GOTO 804
!---------------------------------------------------------------------------

    END SUBROUTINE ITER_PLI