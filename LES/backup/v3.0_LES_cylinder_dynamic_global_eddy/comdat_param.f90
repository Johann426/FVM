	MODULE COMDAT_PARAM                                  

	IMPLICIT NONE
	SAVE

	INTEGER, PARAMETER :: NI	= 81	!201
	INTEGER, PARAMETER :: NJ	= 81	!201
	INTEGER, PARAMETER :: NK	= 11	!121
	
	INTEGER, PARAMETER :: NDIMJ = 6	!8
	INTEGER, PARAMETER :: NDIMK = 2	!4
	
	INTEGER, PARAMETER :: NDIMS = 2
	INTEGER, PARAMETER :: NSWPP = 20
	
	INTEGER, PARAMETER :: IRESTART = 0
	INTEGER, PARAMETER :: USE_LES = 1
	
	INTEGER :: IPOST	= 0
	
    INTEGER :: NITER_W   = 2E5
    INTEGER :: NITER_R   = 1E5
	
	DOUBLE PRECISION ::				RTIME_MAX = 500.D0
	
	DOUBLE PRECISION, PARAMETER ::	RTIME_END = 800.D0
    
	DOUBLE PRECISION, PARAMETER :: REF_A = 6.283185	!	REFERENCE AREA FOR CIRCULAR CYLINDER
!	DOUBLE PRECISION, PARAMETER :: REF_A = 6.236295	!	REFERENCE AREA FOR SKEWEDD CYLINDER WITH A=0.1
!	DOUBLE PRECISION, PARAMETER :: REF_A = 6.236295	!	REFERENCE AREA FOR SKEWEDD CYLINDER WITH A=0.05
    
    DOUBLE PRECISION, PARAMETER :: RESIMAX = 1.0D-4
    DOUBLE PRECISION, PARAMETER :: CFL_USE = 0.29D0
    
    DOUBLE PRECISION :: UIN
!	DOUBLE PRECISION, PARAMETER :: UIN		= 1.D0
	DOUBLE PRECISION, PARAMETER :: TIN		= 0.D0
    DOUBLE PRECISION, PARAMETER :: T_WALL	= 1.D0
    DOUBLE PRECISION, PARAMETER :: VEL_ROT	= 0.D0		!	ROTATING VELOCITY

    DOUBLE PRECISION :: DTIME =	5.0D-3	!2.0D-3
    
	DOUBLE PRECISION :: RE = 3900.D0		! REYNOLDS NUMBER
	DOUBLE PRECISION :: PR = 0.71D0			! PRANDTLE NUMBER
	DOUBLE PRECISION :: PE = 2769.D0		! PECLET NUMBER

	END MODULE COMDAT_PARAM