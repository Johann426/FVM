    SUBROUTINE CALC_VOR
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE

	!	Cell center vorticity

    DO K=KS1,KE1
    DO J=JS1,JE1
    DO I=IS1,IE1
    
    VOR_Z(I,J,K)=XIX(I,J,K)*0.5D0*(V(I+1,J,K)-V(I-1,J,K)) &
                +ETX(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K)) &
                +ZTX(I,J,K)*0.5D0*(V(I,J,K+1)-V(I,J,K-1)) &
                -XIY(I,J,K)*0.5D0*(U(I+1,J,K)-U(I-1,J,K)) &
                -ETY(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K)) &
                -ZTY(I,J,K)*0.5D0*(U(I,J,K+1)-U(I,J,K-1)) 

    VOR_Y(I,J,K)=XIZ(I,J,K)*0.5D0*(U(I+1,J,K)-U(I-1,J,K)) &
                +ETZ(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K)) &
                +ZTZ(I,J,K)*0.5D0*(U(I,J,K+1)-U(I,J,K-1)) &
                -XIX(I,J,K)*0.5D0*(W(I+1,J,K)-W(I-1,J,K)) &
                -ETX(I,J,K)*0.5D0*(W(I,J+1,K)-W(I,J-1,K)) &
                -ZTX(I,J,K)*0.5D0*(W(I,J,K+1)-W(I,J,K-1)) 
                
    VOR_X(I,J,K)=XIY(I,J,K)*0.5D0*(W(I+1,J,K)-W(I-1,J,K)) &
                +ETY(I,J,K)*0.5D0*(W(I,J+1,K)-W(I,J-1,K)) &
                +ZTY(I,J,K)*0.5D0*(W(I,J,K+1)-W(I,J,K-1)) &
                -XIZ(I,J,K)*0.5D0*(V(I+1,J,K)-V(I-1,J,K)) &
                -ETZ(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K)) &
                -ZTZ(I,J,K)*0.5D0*(V(I,J,K+1)-V(I,J,K-1)) 
    ENDDO
    ENDDO
    ENDDO

    DO K=KS1,KE1
    DO J=JS1,JE1
    I=IS
        VOR_Z(I,J,K)=XIX(I,J,K)*0.5D0*(V(I+1,J,K)-(3.D0*V(I,J,K)-3.D0*V(I+1,J,K)+V(I+2,J,K)))	&
                    +ETX(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K))									&
					+ZTX(I,J,K)*0.5D0*(V(I,J,K+1)-V(I,J,K-1))									&
                    -XIY(I,J,K)*0.5D0*(U(I+1,J,K)-(3.D0*U(I,J,K)-3.D0*U(I+1,J,K)+U(I+2,J,K)))	&
                    -ETY(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))									&
					-ZTY(I,J,K)*0.5D0*(U(I,J,K+1)-U(I,J,K-1))
		
		VOR_Y(I,J,K)=XIZ(I,J,K)*0.5D0*(U(I+1,J,K)-(3.D0*U(I,J,K)-3.D0*U(I+1,J,K)+U(I+2,J,K)))	&
					+ETZ(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))									&
				    +ZTZ(I,J,K)*0.5D0*(U(I,J,K+1)-U(I,J,K-1))									&
					-XIX(I,J,K)*0.5D0*(W(I+1,J,K)-(3.D0*W(I,J,K)-3.D0*W(I+1,J,K)+W(I+2,J,K)))	&
	                -ETX(I,J,K)*0.5D0*(W(I,J+1,K)-W(I,J-1,K))									&
		            -ZTX(I,J,K)*0.5D0*(W(I,J,K+1)-W(I,J,K-1))
	
	    VOR_X(I,J,K)=XIY(I,J,K)*0.5D0*(W(I+1,J,K)-(3.D0*W(I,J,K)-3.D0*W(I+1,J,K)+W(I+2,J,K)))	&
		            +ETY(I,J,K)*0.5D0*(W(I,J+1,K)-W(I,J-1,K))									&
			        +ZTY(I,J,K)*0.5D0*(W(I,J,K+1)-W(I,J,K-1))									&
				    -XIZ(I,J,K)*0.5D0*(V(I+1,J,K)-(3.D0*V(I,J,K)-3.D0*V(I+1,J,K)+V(I+2,J,K)))	&
					-ETZ(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K))									&
	                -ZTZ(I,J,K)*0.5D0*(V(I,J,K+1)-V(I,J,K-1)) 
    I=IE
		VOR_Z(I,J,K)=XIX(I,J,K)*0.5D0*((3.D0*V(I,J,K)-3.D0*V(I-1,J,K)+V(I-2,J,K))-V(I-1,J,K))	&
					+ETX(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K))									&
					+ZTX(I,J,K)*0.5D0*(V(I,J,K+1)-V(I,J,K-1))									&
					-XIY(I,J,K)*0.5D0*((3.D0*U(I,J,K)-3.D0*U(I-1,J,K)+U(I-2,J,K))-U(I-1,J,K))	&
					-ETY(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))									&
                    -ZTY(I,J,K)*0.5D0*(U(I,J,K+1)-U(I,J,K-1))
		
		VOR_Y(I,J,K)=XIZ(I,J,K)*0.5D0*((3.D0*U(I,J,K)-3.D0*U(I-1,J,K)+U(I-2,J,K))-U(I-1,J,K))	&
					+ETZ(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))									&
				    +ZTZ(I,J,K)*0.5D0*(U(I,J,K+1)-U(I,J,K-1))									&
					-XIX(I,J,K)*0.5D0*((3.D0*W(I,J,K)-3.D0*W(I-1,J,K)+W(I-2,J,K))-W(I-1,J,K))	&
	                -ETX(I,J,K)*0.5D0*(W(I,J+1,K)-W(I,J-1,K))									&
		            -ZTX(I,J,K)*0.5D0*(W(I,J,K+1)-W(I,J,K-1))
	
	    VOR_X(I,J,K)=XIY(I,J,K)*0.5D0*((3.D0*W(I,J,K)-3.D0*W(I-1,J,K)+W(I-2,J,K))-W(I-1,J,K))	&
		            +ETY(I,J,K)*0.5D0*(W(I,J+1,K)-W(I,J-1,K))									&
			        +ZTY(I,J,K)*0.5D0*(W(I,J,K+1)-W(I,J,K-1))									&
				    -XIZ(I,J,K)*0.5D0*((3.D0*V(I,J,K)-3.D0*V(I-1,J,K)+V(I-2,J,K))-V(I-1,J,K))	&
					-ETZ(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K))									&
	                -ZTZ(I,J,K)*0.5D0*(V(I,J,K+1)-V(I,J,K-1)) 
    ENDDO
    ENDDO
    
    DO J=JS1,JE1
	DO I=IS1,IE1
    K=KS
	    VOR_Z(I,J,K)=XIX(I,J,K)*0.5D0*(V(I+1,J,K)-V(I-1,J,K))									&
		            +ETX(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K))									&
			        +ZTX(I,J,K)*0.5D0*(V(I,J,K+1)-(3.D0*V(I,J,K)-3.D0*V(I,J,K+1)+V(I,J,K+2)))	&
				    -XIY(I,J,K)*0.5D0*(U(I+1,J,K)-U(I-1,J,K))									&
					-ETY(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))									&
	                -ZTY(I,J,K)*0.5D0*(U(I,J,K+1)-(3.D0*U(I,J,K)-3.D0*U(I,J,K+1)+U(I,J,K+2)))

		VOR_Y(I,J,K)=XIZ(I,J,K)*0.5D0*(U(I+1,J,K)-U(I-1,J,K))									&
	                +ETZ(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))									&
		            +ZTZ(I,J,K)*0.5D0*(U(I,J,K+1)-(3.D0*U(I,J,K)-3.D0*U(I,J,K+1)+U(I,J,K+2)))	&
			        -XIX(I,J,K)*0.5D0*(W(I+1,J,K)-W(I-1,J,K))									&
				    -ETX(I,J,K)*0.5D0*(W(I,J+1,K)-W(I,J-1,K))									&
					-ZTX(I,J,K)*0.5D0*(W(I,J,K+1)-(3.D0*W(I,J,K)-3.D0*W(I,J,K+1)+W(I,J,K+2))) 
                
		VOR_X(I,J,K)=XIY(I,J,K)*0.5D0*(W(I+1,J,K)-W(I-1,J,K))									&
	                +ETY(I,J,K)*0.5D0*(W(I,J+1,K)-W(I,J-1,K))									&
		            +ZTY(I,J,K)*0.5D0*(W(I,J,K+1)-(3.D0*W(I,J,K)-3.D0*W(I,J,K+1)+W(I,J,K+2)))	&
			        -XIZ(I,J,K)*0.5D0*(V(I+1,J,K)-V(I-1,J,K))									&
				    -ETZ(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K))									&
					-ZTZ(I,J,K)*0.5D0*(V(I,J,K+1)-(3.D0*V(I,J,K)-3.D0*V(I,J,K+1)+V(I,J,K+2))) 
    K=KE
		VOR_Z(I,J,K)=XIX(I,J,K)*0.5D0*(V(I+1,J,K)-V(I-1,J,K))									&
	                +ETX(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K))									&
		            +ZTX(I,J,K)*0.5D0*((3.D0*V(I,J,K)-3.D0*V(I,J,K-1)+V(I,J,K-2))-V(I-1,J,K))	&
			        -XIY(I,J,K)*0.5D0*(U(I+1,J,K)-U(I-1,J,K))									&
				    -ETY(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))									&
					-ZTY(I,J,K)*0.5D0*((3.D0*U(I,J,K)-3.D0*U(I,J,K-1)+U(I,J,K-2))-U(I-1,J,K)) 

		VOR_Y(I,J,K)=XIZ(I,J,K)*0.5D0*(U(I+1,J,K)-U(I-1,J,K))									&
					+ETZ(I,J,K)*0.5D0*(U(I,J+1,K)-U(I,J-1,K))									&
				    +ZTZ(I,J,K)*0.5D0*((3.D0*U(I,J,K)-3.D0*U(I,J,K-1)+U(I,J,K-2))-U(I-1,J,K))	&
			        -XIX(I,J,K)*0.5D0*(W(I+1,J,K)-W(I-1,J,K))									&
		            -ETX(I,J,K)*0.5D0*(W(I,J+1,K)-W(I,J-1,K))									&
	                -ZTX(I,J,K)*0.5D0*((3.D0*W(I,J,K)-3.D0*W(I,J,K-1)+W(I,J,K-2))-W(I-1,J,K)) 
                
		VOR_X(I,J,K)=XIY(I,J,K)*0.5D0*(W(I+1,J,K)-W(I-1,J,K))									&
	                +ETY(I,J,K)*0.5D0*(W(I,J+1,K)-W(I,J-1,K))									&
		            +ZTY(I,J,K)*0.5D0*((3.D0*W(I,J,K)-3.D0*W(I,J,K-1)+W(I,J,K-2))-W(I-1,J,K))	&
			        -XIZ(I,J,K)*0.5D0*(V(I+1,J,K)-V(I-1,J,K))									&
				    -ETZ(I,J,K)*0.5D0*(V(I,J+1,K)-V(I,J-1,K))									&
					-ZTZ(I,J,K)*0.5D0*((3.D0*V(I,J,K)-3.D0*V(I,J,K-1)+V(I,J,K-2))-V(I-1,J,K)) 
    ENDDO
    ENDDO
    
    
    DO K=KS,KE
    DO I=IS,IE
    J=JS !branch cut region
        VOR_Z(I,J,K) = VOR_Z(I,NJ-1,K)
        VOR_Y(I,J,K) = VOR_Y(I,NJ-1,K)
        VOR_X(I,J,K) = VOR_X(I,NJ-1,K)
        
    J=JE !branch cut region
        VOR_Z(I,J,K) = VOR_Z(I,1,K)
        VOR_Y(I,J,K) = VOR_Z(I,1,K)
        VOR_X(I,J,K) = VOR_Z(I,1,K)
    ENDDO
    ENDDO

    END SUBROUTINE CALC_VOR