	MODULE ComDat_Param                                  

	IMPLICIT NONE
	SAVE

	INTEGER, PARAMETER :: NI	= 201
	INTEGER, PARAMETER :: NJ	= 201
	INTEGER, PARAMETER :: NK	= 121
	
	INTEGER, PARAMETER :: NDIMJ = 8
	INTEGER, PARAMETER :: NDIMK = 4
	
	INTEGER, PARAMETER :: NDIMS = 2
	INTEGER, PARAMETER :: NSWPU=1
	INTEGER, PARAMETER :: NSWPP=20
	INTEGER, PARAMETER :: NSWPT=1

    INTEGER :: NITER_W   = 1E6
    INTEGER :: NITER_R   = 2E4
    INTEGER :: RTIME_MAX = 300
    
   	DOUBLE PRECISION, PARAMETER :: REF_A = 6.283185307	!	REFERENCE AREA FOR CIRCULAR CYLINDER
!	DOUBLE PRECISION, PARAMETER :: REF_A = 6.283185307	!	REFERENCE AREA FOR WAVY CYLINDER
!	DOUBLE PRECISION, PARAMETER :: REF_A = 6.283185307	!	REFERENCE AREA FOR SKEWEDD CYLINDER
    
    DOUBLE PRECISION, PARAMETER :: RESIMAX = 1.0D-6
    DOUBLE PRECISION, PARAMETER :: CFL_USE = 0.29D0
    
    DOUBLE PRECISION :: UIN
!	DOUBLE PRECISION, PARAMETER :: UIN     = 1.D0
	DOUBLE PRECISION, PARAMETER :: TIN     = 0.D0
    DOUBLE PRECISION, PARAMETER :: T_WALL  = 1.D0

    DOUBLE PRECISION :: DTIME = 2.5D-3
!	DOUBLE PRECISION :: VEL_ROT
!--------------------------------------------------------    
	INTEGER :: IRESTART = 0
	INTEGER :: IPOST	= 0
	DOUBLE PRECISION :: RE = 300.D0     ! REYNOLDS NUMBER
	DOUBLE PRECISION :: PR = 0.7D0      ! PRANDTLE NUMBER
	DOUBLE PRECISION :: PE = 210.D0     ! PECLET NUMBER
!--------------------------------------------------------

	END MODULE ComDat_Param