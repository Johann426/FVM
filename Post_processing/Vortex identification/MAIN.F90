	PROGRAM MAIN
	
	USE COMDAT_SHARED
	
	IMPLICIT NONE
	SAVE
	
!--------------------------------------------------------------------------------------------
!	FROM 2D DATA
	
	IS =0
	IS1=1
	IE=NI
	IE1=NI-1
	
	JS =0
	JS1=1
	JE=NJ
	JE1=NJ-1
	
	KS =0
	KS1=1
	KE=NK
	KE1=NK-1	
	
	OPEN(UNIT=10,FILE='081x081.dat',STATUS='OLD')
	K=1
	DO J=1,NJ
	DO I=1,NI
		READ(10,*) XG(I,J,K),YG(I,J,K)
	ENDDO
	ENDDO
	CLOSE(10)
	
	DO K=1,NK
	DO J=1,NJ
	DO I=1,NI
		XG(I,J,K) = XG(I,J,1)
		YG(I,J,K) = YG(I,J,1)
		ZG(I,J,K) = 3.141592*DBLE(K)/DBLE(NK)
	ENDDO
	ENDDO
	ENDDO
	
	OPEN(UNIT=10,FILE='test.dat',STATUS='OLD')
	K=0
	DO J=0,NJ
	DO I=0,NI
		READ(10,*) X(I,J,K),Y(I,J,K),U(I,J,K),V(I,J,K),DUM,DUM_A,DUM_B
	ENDDO
	ENDDO
	CLOSE(10)
	
	DO K=0,NK
	DO J=0,NJ
	DO I=0,NI
		X(I,J,K) = X(I,J,0)
		Y(I,J,K) = Y(I,J,0)
		Z(I,J,K) = 3.141592*DBLE(K)/DBLE(NK)
		U(I,J,K) = U(I,J,0)
		V(I,J,K) = V(I,J,0)
	ENDDO
	ENDDO
	ENDDO
	
	!   CENTERED DIFFERENCE VALUE
    DO K=KS1,KE1
	DO J=JS1,JE1
	DO I=IS1,IE1
	    XXI(I,J,K)=0.5D0*(X(I+1,J,K)-X(I-1,J,K))
	    XET(I,J,K)=0.5D0*(X(I,J+1,K)-X(I,J-1,K))
		XZT(I,J,K)=0.5D0*(X(I,J,K+1)-X(I,J,K-1))

	    YXI(I,J,K)=0.5D0*(Y(I+1,J,K)-Y(I-1,J,K))
	    YET(I,J,K)=0.5D0*(Y(I,J+1,K)-Y(I,J-1,K))
		YZT(I,J,K)=0.5D0*(Y(I,J,K+1)-Y(I,J,K-1))

		ZXI(I,J,K)=0.5D0*(Z(I+1,J,K)-Z(I-1,J,K))
		ZET(I,J,K)=0.5D0*(Z(I,J+1,K)-Z(I,J-1,K))
		ZZT(I,J,K)=0.5D0*(Z(I,J,K+1)-Z(I,J,K-1))
	ENDDO
	ENDDO
	ENDDO

	!   ON THE BOUNDARY FACES
	DO K=KS1,KE1
	DO J=JS1,JE1
	I=IS
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
	    XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
	    
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
		
		ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
		ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
	I=IE
        XXI(NI,J,K)=(3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(NI,J,K)=(X(NI,J+1,K)-X(NI,J-1,K))/2.0D0
	    XZT(NI,J,K)=(X(NI,J,K+1)-X(NI,J,K-1))/2.0D0
	    
        YXI(NI,J,K)=(3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(NI,J,K)=(Y(NI,J+1,K)-Y(NI,J-1,K))/2.0D0
        YZT(NI,J,K)=(Y(NI,J,K+1)-Y(NI,J,K-1))/2.0D0
		
		ZXI(NI,J,K)=(3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(NI,J,K)=(Z(NI,J+1,K)-Z(NI,J-1,K))/2.0D0
		ZZT(NI,J,K)=(Z(NI,J,K+1)-Z(NI,J,K-1))/2.0D0
	ENDDO
	ENDDO
	
    DO K=KS1,KE1
    DO I=IS1,IE1
    J=JS
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
	    XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
	    
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
        
        ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
		ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
	J=JE
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
	    XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
	    
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
        
        ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
		ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
	ENDDO
	ENDDO 

    DO J=JS1,JE1
	DO I=IS1,IE1
	K=KS
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
		XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
	
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
		YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0

		ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
		ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    K=KE
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
	

	!   ON THE BOUNDARY LINES
    DO I=IS1,IE1
    J=JS
    K=KS
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
        
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0
        
        ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    J=JE
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
        
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0
        
        ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    K=KE
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=(3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
        
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=(3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0
        
        ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=(3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    J=JS
        XXI(I,J,K)=(X(I+1,J,K)-X(I-1,J,K))/2.0D0
        XET(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
        
        YXI(I,J,K)=(Y(I+1,J,K)-Y(I-1,J,K))/2.0D0
        YET(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0
        
        ZXI(I,J,K)=(Z(I+1,J,K)-Z(I-1,J,K))/2.0D0
        ZET(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    ENDDO
    
    DO J=JS1,JE1
    K=KS
    I=IS
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0

        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0

        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    K=KE
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0

        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0

        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    I=IE
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(X(I,J+1,K)-X(I,J-1,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0

        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(Y(I,J+1,K)-Y(I,J-1,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0

        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=(Z(I,J+1,K)-Z(I,J-1,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    K=KS
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
    
	DO K=KS1,KE1	
	I=IS
	J=JS
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0                                                                                                                                                  
        XET(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
        
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
        
        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
	J=JE
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
        
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
        
        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
    I=IE
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
        
        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
        
        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
    J=JS
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=(X(I,J,K+1)-X(I,J,K-1))/2.0D0
        
        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=(Y(I,J,K+1)-Y(I,J,K-1))/2.0D0
        
        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=(Z(I,J,K+1)-Z(I,J,K-1))/2.0D0
	ENDDO


	!   AT THE BOUNDARY POINTSS
    I=IS
    J=JS
    K=KS
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
	
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0
		  
        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    J=JE
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
	
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0
		  
        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    I=IE
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
	
        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0
		  
        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    J=JS
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J,K+1)-X(I,J,K+2))/2.0D0
	
        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J,K+1)-Y(I,J,K+2))/2.0D0
		  
        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J,K+1)-Z(I,J,K+2))/2.0D0
    K=KE
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
	
        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0
		  
        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    I=IS
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I,J+1,K)-X(I,J+2,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
	
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I,J+1,K)-Y(I,J+2,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0
		  
        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I,J+1,K)-Z(I,J+2,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    J=JE
        XXI(I,J,K)=(-3.0D0*X(I,J,K)+4.0D0*X(I+1,J,K)-X(I+2,J,K))/2.0D0
        XET(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
	
        YXI(I,J,K)=(-3.0D0*Y(I,J,K)+4.0D0*Y(I+1,J,K)-Y(I+2,J,K))/2.0D0
        YET(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0
		  
        ZXI(I,J,K)=(-3.0D0*Z(I,J,K)+4.0D0*Z(I+1,J,K)-Z(I+2,J,K))/2.0D0
        ZET(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
    I=IE
        XXI(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I-1,J,K)+X(I-2,J,K))/2.0D0
        XET(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J-1,K)+X(I,J-2,K))/2.0D0
        XZT(I,J,K)=( 3.0D0*X(I,J,K)-4.0D0*X(I,J,K-1)+X(I,J,K-2))/2.0D0
	
        YXI(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I-1,J,K)+Y(I-2,J,K))/2.0D0
        YET(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J-1,K)+Y(I,J-2,K))/2.0D0
        YZT(I,J,K)=( 3.0D0*Y(I,J,K)-4.0D0*Y(I,J,K-1)+Y(I,J,K-2))/2.0D0
		  
        ZXI(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I-1,J,K)+Z(I-2,J,K))/2.0D0
        ZET(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J-1,K)+Z(I,J-2,K))/2.0D0
        ZZT(I,J,K)=( 3.0D0*Z(I,J,K)-4.0D0*Z(I,J,K-1)+Z(I,J,K-2))/2.0D0
	
	!   DETERMINENT OF INVERSE JACOBIAN MATRIX : INVERSE OF JACOBIAN
	DO K=KS,KE
	DO J=JS,JE
	DO I=IS,IE
         DJAC(I,J,K)	&
         =  XXI(I,J,K)*(YET(I,J,K)*ZZT(I,J,K)-YZT(I,J,K)*ZET(I,J,K))	&
           -YXI(I,J,K)*(XET(I,J,K)*ZZT(I,J,K)-ZET(I,J,K)*XZT(I,J,K))	&
           +ZXI(I,J,K)*(XET(I,J,K)*YZT(I,J,K)-XZT(I,J,K)*YET(I,J,K))	
        !  CHECK THE JACOBIAN ERROR
        IF (DJAC(I,J,K).LT.1.0E-8) WRITE(6,*) I,J,K,'DJAC=',DJAC(I,J,K)
	ENDDO
	ENDDO
	ENDDO

	DO K=KS,KE
	DO J=JS,JE
	DO I=IS,IE
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
!--------------------------------------------------------------------------------------------

	DO K=KS1,KE1
	DO J=JS1,JE1
	DO I=IS1,IE1
	
	DUDX=0.5D0*(XIX(I,J,K)*(U(I+1,J,K)-U(I-1,J,K))+ETX(I,J,K)*(U(I,J+1,K)-U(I,J-1,K))+ZTX(I,J,K)*(U(I,J,K+1)-U(I,J,K-1)))
	DVDX=0.5D0*(XIX(I,J,K)*(V(I+1,J,K)-V(I-1,J,K))+ETX(I,J,K)*(V(I,J+1,K)-V(I,J-1,K))+ZTX(I,J,K)*(V(I,J,K+1)-V(I,J,K-1)))
	DWDX=0.5D0*(XIX(I,J,K)*(W(I+1,J,K)-W(I-1,J,K))+ETX(I,J,K)*(W(I,J+1,K)-W(I,J-1,K))+ZTX(I,J,K)*(W(I,J,K+1)-W(I,J,K-1)))
	
	DUDY=0.5D0*(XIY(I,J,K)*(U(I+1,J,K)-U(I-1,J,K))+ETY(I,J,K)*(U(I,J+1,K)-U(I,J-1,K))+ZTY(I,J,K)*(U(I,J,K+1)-U(I,J,K-1)))
	DVDY=0.5D0*(XIY(I,J,K)*(V(I+1,J,K)-V(I-1,J,K))+ETY(I,J,K)*(V(I,J+1,K)-V(I,J-1,K))+ZTY(I,J,K)*(V(I,J,K+1)-V(I,J,K-1)))
	DWDY=0.5D0*(XIY(I,J,K)*(W(I+1,J,K)-W(I-1,J,K))+ETY(I,J,K)*(W(I,J+1,K)-W(I,J-1,K))+ZTY(I,J,K)*(W(I,J,K+1)-W(I,J,K-1)))
	
	DUDZ=0.5D0*(XIZ(I,J,K)*(U(I+1,J,K)-U(I-1,J,K))+ETZ(I,J,K)*(U(I,J+1,K)-U(I,J-1,K))+ZTZ(I,J,K)*(U(I,J,K+1)-U(I,J,K-1)))
	DVDZ=0.5D0*(XIZ(I,J,K)*(V(I+1,J,K)-V(I-1,J,K))+ETZ(I,J,K)*(V(I,J+1,K)-V(I,J-1,K))+ZTZ(I,J,K)*(V(I,J,K+1)-V(I,J,K-1)))
	DWDZ=0.5D0*(XIZ(I,J,K)*(W(I+1,J,K)-W(I-1,J,K))+ETZ(I,J,K)*(W(I,J+1,K)-W(I,J-1,K))+ZTZ(I,J,K)*(W(I,J,K+1)-W(I,J,K-1)))
	
!	SYMMETRIC PARTS OF THE VELOCITY GRADIENT TENSOR
	
	S1(1,1) = 0.5D0*(DUDX+DUDX)
	S1(2,1) = 0.5D0*(DVDX+DUDY)
	S1(3,1) = 0.5D0*(DWDX+DUDz)
	
	S1(1,2) = 0.5D0*(DUDY+DVDX)
	S1(2,2) = 0.5D0*(DVDY+DVDY)
	S1(3,2) = 0.5D0*(DWDY+DVDZ)
	
	S1(1,3) = 0.5D0*(DUDZ+DWDX)
	S1(2,3) = 0.5D0*(DVDZ+DWDY)
	S1(3,3) = 0.5D0*(DWDZ+DWDZ)
	
	O1(1,1) = 0.5D0*(DUDX-DUDX)
	O1(2,1) = 0.5D0*(DVDX-DUDY)
	O1(3,1) = 0.5D0*(DWDX-DUDz)
	
	O1(1,2) = 0.5D0*(DUDY-DVDX)
	O1(2,2) = 0.5D0*(DVDY-DVDY)
	O1(3,2) = 0.5D0*(DWDY-DVDZ)
	
	O1(1,3) = 0.5D0*(DUDZ-DWDX)
	O1(2,3) = 0.5D0*(DVDZ-DWDY)
	O1(3,3) = 0.5D0*(DWDZ-DWDZ)
	
	S2 = MATMUL(S1,S1)
	O2 = MATMUL(O1,O1)
		
!	S2(:,:) = 0.D0
!	O2(:,:) = 0.D0
!	
!	DO N=1,3
!	DO M=1,3
!	DO L=1,3
!		S2(M,N) = S2(M,N) + S1(L,N)*S1(M,L)
!		O2(M,N) = O2(M,N) + O1(L,N)*O1(M,L)
!	ENDDO
!	ENDDO
!	ENDDO
	
	S2O2(1:3,1:3) = S2(1:3,1:3) + O2(1:3,1:3)
	
	DUM_A = -1.D0
	DUM_B = +S2O2(1,1)+S2O2(2,2)+S2O2(3,3)
	DUM_C = -S2O2(1,1)*S2O2(2,2)-S2O2(2,2)*S2O2(3,3)-S2O2(3,3)*S2O2(1,1)	&
			+S2O2(1,2)*S2O2(2,1)+S2O2(2,3)*S2O2(3,2)+S2O2(3,1)*S2O2(1,3)
	DUM_D = +S2O2(1,1)*S2O2(2,2)*S2O2(3,3)+S2O2(1,2)*S2O2(2,3)*S2O2(3,1)+S2O2(1,3)*S2O2(2,1)*S2O2(3,2)	&
			-S2O2(1,1)*S2O2(2,3)*S2O2(3,2)-S2O2(2,2)*S2O2(1,3)*S2O2(3,1)-S2O2(3,3)*S2O2(1,2)*S2O2(2,1)
	
	CALL FIND_ROOT(DUM_A,DUM_B,DUM_C,DUM_D,LAMDA)
	
	LAMDA1 = MAX(LAMDA(1),LAMDA(2),LAMDA(3))
	LAMDA3 = MIN(LAMDA(1),LAMDA(2),LAMDA(3))
	
		DO N=1,3
			IF(LAMDA(N).LT.LAMDA1 .AND. LAMDA(N).GT.LAMDA3)	LAMDA2(I,J,K) = LAMDA(N)
		ENDDO
	
	ENDDO
	ENDDO
	ENDDO
	
	LAMDA2(IS,JS:JE,KS:KE) = 3.D0*LAMDA2(IS+1,JS:JE,KS:KE)-3.D0*LAMDA2(IS+2,JS:JE,KS:KE)+LAMDA2(IS+3,JS:JE,KS:KE)
	LAMDA2(IE,JS:JE,KS:KE) = 3.D0*LAMDA2(IE-1,JS:JE,KS:KE)-3.D0*LAMDA2(IE-2,JS:JE,KS:KE)+LAMDA2(IE-3,JS:JE,KS:KE)
	LAMDA2(IS:IE,JS,KS:KE) = 3.D0*LAMDA2(IS:IE,JS+1,KS:KE)-3.D0*LAMDA2(IS:IE,JS+2,KS:KE)+LAMDA2(IS:IE,JS+3,KS:KE)
	LAMDA2(IS:IE,JE,KS:KE) = 3.D0*LAMDA2(IS:IE,JE-1,KS:KE)-3.D0*LAMDA2(IS:IE,JE-2,KS:KE)+LAMDA2(IS:IE,JE-3,KS:KE)
	LAMDA2(IS:IE,JS:JE,KS) = 3.D0*LAMDA2(IS:IE,JS:JE,KS+1)-3.D0*LAMDA2(IS:IE,JS:JE,KS+2)+LAMDA2(IS:IE,JS:JE,KS+3)
	LAMDA2(IS:IE,JS:JE,KE) = 3.D0*LAMDA2(IS:IE,JS:JE,KE-1)-3.D0*LAMDA2(IS:IE,JS:JE,KE-2)+LAMDA2(IS:IE,JS:JE,KE-3)
	
	LAMDAG(1:NI,1:NJ,1:NK) =0.125D0*(																					&
	+LAMDA2(1:NI,1:NJ,1:NK)+LAMDA2(1-1:NI-1,1:NJ,1:NK)+LAMDA2(1:NI,1-1:NJ-1,1:NK)+LAMDA2(1-1:NI-1,1-1:NJ-1,1:NK)		&
	+LAMDA2(1:NI,1:NJ,0:NK-1)+LAMDA2(1-1:NI-1,1:NJ,0:NK-1)+LAMDA2(1:NI,1-1:NJ-1,0:NK-1)+LAMDA2(1-1:NI-1,1-1:NJ-1,0:NK-1))
		
	!   NODAL RESULTS
	OPEN(UNIT=100,FILE='VOR_STR.PLT',STATUS='UNKNOWN')
	WRITE(100,*) 'TITLE="VORTICAL STRUCTURE"'
	WRITE(100,*) 'variables=x,y,z,l'
    WRITE(100,24) NI,NJ,NK
24  FORMAT(TR1,'ZONE T="WHOLE", I=',I6,', j=',I6,', k=',I6)
	DO K=KS1,KE
    DO J=JS1,JE
    DO I=IS1,IE
		WRITE(100,*) XG(I,J,K),YG(I,J,K),ZG(I,J,K),LAMDAG(I,J,K)
	ENDDO
	ENDDO
	ENDDO
	CLOSE(100)
	
	END PROGRAM MAIN