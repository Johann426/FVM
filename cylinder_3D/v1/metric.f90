	SUBROUTINE METRIC
	
	USE COMDAT_SHARED

	IMPLICIT NONE
	SAVE
!----------------------------------------------------------
!   CENTERED DIFFERENCE VALUE
    DO K=1,NKM1
	DO J=1,NJM1
	DO I=1,NIM1
	    XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
	    XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
		XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0

	    YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
	    YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
		YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0

		ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
		ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
		ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
	ENDDO
	ENDDO
	ENDDO
!----------------------------------------------------------
!   ON THE BOUNDARY FACES
	DO K=1,NKM1
	DO J=1,NJM1
	I=0
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
	    XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
	    
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
		
		ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
		ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
	I=NI
        XXI(NI,J,K)=(3.0D0*X(NI,J,K)-4.0D0*X(NIM1,J,K)+X(NIM2,J,K))/2.0D0
        XET(NI,J,K)=(X(NI,J+1,K)-X(NI,J-1,K))/2.0D0
	    XZT(NI,J,K)=(X(NI,J,K+1)-X(NI,J,K-1))/2.0D0
	    
        YXI(NI,J,K)=(3.0D0*Y(NI,J,K)-4.0D0*Y(NIM1,J,K)+Y(NIM2,J,K))/2.0D0
        YET(NI,J,K)=(Y(NI,J+1,K)-Y(NI,J-1,K))/2.0D0
        YZT(NI,J,K)=(Y(NI,J,K+1)-Y(NI,J,K-1))/2.0D0
		
		ZXI(NI,J,K)=(3.0D0*Z(NI,J,K)-4.0D0*Z(NIM1,J,K)+Z(NIM2,J,K))/2.0D0
        ZET(NI,J,K)=(Z(NI,J+1,K)-Z(NI,J-1,K))/2.0D0
		ZZT(NI,J,K)=(Z(NI,J,K+1)-Z(NI,J,K-1))/2.0D0
	ENDDO
	ENDDO
	
    DO K=1,NKM1
    DO I=1,NIM1
    J=0
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(X(I,NJ,K)-X(I,NJ-2,K))/2.0D0 !(-3.0D0*X(I,0,K)+4.0D0*X(I,1,K)-X(I,2,K))/2.0D0
	    XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
	    
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(Y(I,NJ,K)-Y(I,NJ-2,K))/2.0D0 !(-3.0D0*Y(I,0,K)+4.0D0*Y(I,1,K)-Y(I,2,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
        
        ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,NJ,K)-Z(I,NJ-2,K))/2.0D0 !(-3.0D0*Z(I,0,K)+4.0D0*Z(I,1,K)-Z(I,2,K))/2.0D0
		ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
	J=NJ
        XXI(I,NJ,K)=(X(I+1,NJ,K)-X(I-1,NJ,K))/2.0D0
        XET(I,NJ,K)=(X(I,2,K)-X(I,0,K))/2.0D0 !(3.0D0*X(I,NJ,K)-4.0D0*X(I,NJM1,K)+X(I,NJM2,K))/2.0D0
	    XZT(I,NJ,K)=(X(I,NJ,K+1)-X(I,NJ,K-1))/2.0D0
	    
        YXI(I,NJ,K)=(Y(I+1,NJ,K)-Y(I-1,NJ,K))/2.0D0
        YET(I,NJ,K)=(Y(I,2,K)-Y(I,0,K))/2.0D0 !(3.0D0*Y(I,NJ,K)-4.0D0*Y(I,NJM1,K)+Y(I,NJM2,K))/2.0D0
        YZT(I,NJ,K)=(Y(I,NJ,K+1)-Y(I,NJ,K-1))/2.0D0
        
        ZXI(I,NJ,K)=(Z(I+1,NJ,K)-Z(I-1,NJ,K))/2.0D0
        ZET(I,NJ,K)=(Z(I,2,K)-Z(I,0,K))/2.0D0 !(3.0D0*Z(I,NJ,K)-4.0D0*Z(I,NJM1,K)+Z(I,NJM2,K))/2.0D0
		ZZT(I,NJ,K)=(Z(I,NJ,K+1)-Z(I,NJ,K-1))/2.0D0
	ENDDO
	ENDDO 

    DO J=1,NJM1
	DO I=1,NIM1
	K=0
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
		XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
	
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
		YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0

		ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
		ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    K=NK
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
		XZT(I,J,K)=(3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
	
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
		YZT(I,J,K)=(3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0

		ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
		ZZT(I,J,K)=(3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
	ENDDO
	ENDDO
	
!----------------------------------------------------------
!   ON THE BOUNDARY LINES
    DO I=1,NIM1
    J=0
    K=0
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(X(I,NJ,K)-X(I,NJ-2,K))/2.0D0 !(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
        
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(Y(I,NJ,K)-Y(I,NJ-2,K))/2.0D0 !(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0
        
        ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,NJ,K)-Z(I,NJ-2,K))/2.0D0 !(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    J=NJ
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(X(I,2,K)-X(I,0,K))/2.0D0 !(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
        
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(Y(I,2,K)-Y(I,0,K))/2.0D0 !(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0
        
        ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,2,K)-Z(I,0,K))/2.0D0 !(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    K=NK
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(X(I,2,K)-X(I,0,K))/2.0D0 !(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
        
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(Y(I,2,K)-Y(I,0,K))/2.0D0 !(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0
        
        ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,2,K)-Z(I,0,K))/2.0D0 !(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    J=0
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(X(I,NJ,K)-X(I,NJ-2,K))/2.0D0 !(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
        
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(Y(I,NJ,K)-Y(I,NJ-2,K))/2.0D0 !(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0
        
        ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,NJ,K)-Z(I,NJ-2,K))/2.0D0 !(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    ENDDO
    
    DO J=1,NJM1
    K=0
    I=0
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0

        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0

        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    K=NK
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0

        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0

        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    I=NI
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0

        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0

        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    K=0
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0

        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0

        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    ENDDO
    
	DO K=1,NKM1	
	I=0
	J=0
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0                                                                                                                                                  
        XET(I,J,K)=(X(I,NJ,K)-X(I,NJ-2,K))/2.0D0 !(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
        
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,NJ,K)-Y(I,NJ-2,K))/2.0D0 !(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
        
        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,NJ,K)-Z(I,NJ-2,K))/2.0D0 !(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
	J=NJ
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(X(I,2,K)-X(I,0,K))/2.0D0 !( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
        
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,2,K)-Y(I,0,K))/2.0D0 !( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
        
        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,2,K)-Z(I,0,K))/2.0D0 !( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
    I=NI
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(X(I,2,K)-X(I,0,K))/2.0D0 !( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
        
        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,2,K)-Y(I,0,K))/2.0D0 !( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
        
        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,2,K)-Z(I,0,K))/2.0D0 !( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
    J=0
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(X(I,NJ,K)-X(I,NJ-2,K))/2.0D0 !(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
        
        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,NJ,K)-Y(I,NJ-2,K))/2.0D0 !(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
        
        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,NJ,K)-Z(I,NJ-2,K))/2.0D0 !(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
	ENDDO

!----------------------------------------------------------
!   AT THE BOUNDARY POINTSS
    I=0
    J=0
    K=0
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(X(I,NJ,K)-X(I,NJ-2,K))/2.0D0 !(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
	
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,NJ,K)-Y(I,NJ-2,K))/2.0D0 !(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0
		  
        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,NJ,K)-Z(I,NJ-2,K))/2.0D0 !(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    J=NJ
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(X(I,2,K)-X(I,0,K))/2.0D0 !( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
	
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,2,K)-Y(I,0,K))/2.0D0 !( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0
		  
        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,2,K)-Z(I,0,K))/2.0D0 !( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    I=NI
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(X(I,2,K)-X(I,0,K))/2.0D0 !( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
	
        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,2,K)-Y(I,0,K))/2.0D0 !( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0
		  
        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,2,K)-Z(I,0,K))/2.0D0 !( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    J=0
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(X(I,NJ,K)-X(I,NJ-2,K))/2.0D0 !(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
	
        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,NJ,K)-Y(I,NJ-2,K))/2.0D0 !(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0
		  
        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,NJ,K)-Z(I,NJ-2,K))/2.0D0 !(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    K=NK
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(X(I,NJ,K)-X(I,NJ-2,K))/2.0D0 !(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
	
        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,NJ,K)-Y(I,NJ-2,K))/2.0D0 !(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0
		  
        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,NJ,K)-Z(I,NJ-2,K))/2.0D0 !(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    I=0
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(X(I,NJ,K)-X(I,NJ-2,K))/2.0D0 !(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
	
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,NJ,K)-Y(I,NJ-2,K))/2.0D0 !(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0
		  
        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,NJ,K)-Z(I,NJ-2,K))/2.0D0 !(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    J=NJ
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(X(I,2,K)-X(I,0,K))/2.0D0 !( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
	
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,2,K)-Y(I,0,K))/2.0D0 !( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0
		  
        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,2,K)-Z(I,0,K))/2.0D0 !( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    I=NI
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(X(I,2,K)-X(I,0,K))/2.0D0 !( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
	
        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,2,K)-Y(I,0,K))/2.0D0 !( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0
		  
        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,2,K)-Z(I,0,K))/2.0D0 !( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
!----------------------------------------------------------
        
!----------------------------------------------------------
!   DETERMINENT OF INVERSE JACOBIAN MATRIX : INVERSE OF JACOBIAN
	DO K=0,NK
	DO J=0,NJ
	DO I=0,NI
         DJAC(I,J,K)	&
         =  XXI(I,J,K)*(YET(I,J,K)*ZZT(I,J,K)-YZT(I,J,K)*ZET(I,J,K))	&
           -YXI(I,J,K)*(XET(I,J,K)*ZZT(I,J,K)-ZET(I,J,K)*XZT(I,J,K))	&
           +ZXI(I,J,K)*(XET(I,J,K)*YZT(I,J,K)-XZT(I,J,K)*YET(I,J,K))	
        !  CHECK THE JACOBIAN ERROR
        IF (DJAC(I,J,K).LT.1.0E-8) WRITE(6,*) I,J,K,'DJAC=',DJAC(I,J,K)
	ENDDO
	ENDDO
	ENDDO
!----------------------------------------------------------
	DO I=0,NI
    DO J=0,NJ
 	DO K=0,NK
          XIX(I,J,K) = (YET(I,J,K)*ZZT(I,J,K)-YZT(I,J,K)*ZET(I,J,K))/DJAC(I,J,K)
          XIY(I,J,K) = (ZET(I,J,K)*XZT(I,J,K)-ZZT(I,J,K)*XET(I,J,K))/DJAC(I,J,K)
          XIZ(I,J,K) = (XET(I,J,K)*YZT(I,J,K)-XZT(I,J,K)*YET(I,J,K))/DJAC(I,J,K)

          ETX(I,J,K) = (YZT(I,J,K)*ZXI(I,J,K)-YXI(I,J,K)*ZZT(I,J,K))/DJAC(I,J,K)
          ETY(I,J,K) = (ZZT(I,J,K)*XXI(I,J,K)-ZXI(I,J,K)*XZT(I,J,K))/DJAC(I,J,K)
          ETZ(I,J,K) = (XZT(I,J,K)*YXI(I,J,K)-XXI(I,J,K)*YZT(I,J,K))/DJAC(I,J,K)

          ZTX(I,J,K) = (YXI(I,J,K)*ZET(I,J,K)-YET(I,J,K)*ZXI(I,J,K))/DJAC(I,J,K)
          ZTY(I,J,K) = (ZXI(I,J,K)*XET(I,J,K)-ZET(I,J,K)*XXI(I,J,K))/DJAC(I,J,K)
          ZTZ(I,J,K) = (XXI(I,J,K)*YET(I,J,K)-XET(I,J,K)*YXI(I,J,K))/DJAC(I,J,K)
	ENDDO
	ENDDO
	ENDDO

!   Q11,Q12,Q21,Q22 CALCULATION
	DO I=0,NI
	DO J=0,NJ
 	DO K=0,NK
          Q11(I,J,K) = XIX(I,J,K)**2.0D0    +XIY(I,J,K)**2.0D0    +XIZ(I,J,K)**2.0D0
          Q12(I,J,K) = XIX(I,J,K)*ETX(I,J,K)+XIY(I,J,K)*ETY(I,J,K)+XIZ(I,J,K)*ETZ(I,J,K)
          Q13(I,J,K) = XIX(I,J,K)*ZTX(I,J,K)+XIY(I,J,K)*ZTY(I,J,K)+XIZ(I,J,K)*ZTZ(I,J,K)

          Q21(I,J,K) = Q12(I,J,K)
          Q22(I,J,K) = ETX(I,J,K)**2.0D0    +ETY(I,J,K)**2.0D0    +ETZ(I,J,K)**2.0D0
          Q23(I,J,K) = ETX(I,J,K)*ZTX(I,J,K)+ETY(I,J,K)*ZTY(I,J,K)+ETZ(I,J,K)*ZTZ(I,J,K)

          Q31(I,J,K) = Q13(I,J,K)
          Q32(I,J,K) = Q23(I,J,K)
          Q33(I,J,K) = ZTX(I,J,K)**2.0D0    +ZTY(I,J,K)**2.0D0    +ZTZ(I,J,K)**2.0D0
	ENDDO
	ENDDO
	ENDDO
		
!	DE,DW,DN,DS CALCULATION
 	DO I=1,NIM1
	DO J=1,NJM1
	DO K=1,NKM1
 	    DE(I,J,K)=0.5D0/DJAC(I,J,K)*(DJAC(I,J,K)*Q11(I,J,K)+DJAC(I+1,J,K)*Q11(I+1,J,K))
	    DW(I,J,K)=0.5D0/DJAC(I,J,K)*(DJAC(I,J,K)*Q11(I,J,K)+DJAC(I-1,J,K)*Q11(I-1,J,K))
	    DN(I,J,K)=0.5D0/DJAC(I,J,K)*(DJAC(I,J,K)*Q22(I,J,K)+DJAC(I,J+1,K)*Q22(I,J+1,K))
	    DS(I,J,K)=0.5D0/DJAC(I,J,K)*(DJAC(I,J,K)*Q22(I,J,K)+DJAC(I,J-1,K)*Q22(I,J-1,K))
        DT(I,J,K)=0.5D0/DJAC(I,J,K)*(DJAC(I,J,K)*Q33(I,J,K)+DJAC(I,J,K+1)*Q33(I,J,K+1))
		DB(I,J,K)=0.5D0/DJAC(I,J,K)*(DJAC(I,J,K)*Q33(I,J,K)+DJAC(I,J,K-1)*Q33(I,J,K-1))
    ENDDO
	ENDDO
	ENDDO
!------------------------------------------------------------------------------
	END SUBROUTINE METRIC