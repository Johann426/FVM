    SUBROUTINE BC_PLI

    USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE

    DUM_3D=PLI
!-----------------------------------------
!   CYLINDER SURFACE AND INLET & OUTLET
!-----------------------------------------
	DO K=1,NKM1
	DO J=1,NJM1
	I=0
	! WALL BOUNDARY SURFACE
    !AP(I+1,J,K) = AP(I+1,J,K)-AW(I+1,J,K)
    AW(I+1,J,K) = 0.0D0
    PLI(I,J,K) = 1.D0/3.D0*((4.D0*PLI(I+1,J,k)-PLI(I+2,J,K))					&
                +Q12(I,J,K)/Q11(I,J,K) *0.5D0*(DUM_3D(I,J+1,K)-DUM_3D(I,J-1,K)))

    I=NI
	! OUTFLOW BOUNDARY SURFACE
    !AP(I-1,J,K) = AP(I-1,J,K)-AE(I-1,J,K)
    AE(I-1,J,K) = 0.0D0
    PLI(I,J,K) = 1.D0/3.D0*((4.D0*PLI(I-1,J,k)-PLI(I-2,J,K))					&
 				-Q12(I,J,K)/Q11(I,J,K) *0.5D0*(DUM_3D(I,J+1,K)-DUM_3D(I,J-1,K)))
	ENDDO
	ENDDO
!-------------------------------------------
!   BRANCH CUT REGION
!-------------------------------------------
	DO K=0,NK
	DO I=0,NI
    J=0
	    PLI(I,J,K) = PLI(I,NJ-1,K)
	J=NJ
		PLI(I,J,K) = PLI(I,1,K)
	ENDDO
	ENDDO
!-------------------------------------------
!   PERIODIC WALL OF BOUNDARY
!-------------------------------------------
	DO J=0,NJ
	DO I=0,NI
	K=0
!		AP(I,J,SZ1) = AP(I,J,SZ1)-AB(I,J,SZ1)
!		AB(I,J,SZ1) = 0.D0
!		PLI(I,J,K) = PLI(I,J,K+1)
		PLI(I,J,K) = PLI(I,J,NKM1)
	K=NK
!		AP(I,J,EZ1) = AP(I,J,EZ1)-AT(I,J,EZ1)
!		AT(I,J,EZ1) = 0.D0
!		PLI(I,J,K) = PLI(I,J,K-1)
		PLI(I,J,K) = PLI(I,J,1)
	ENDDO
	ENDDO
	
    END SUBROUTINE BC_PLI