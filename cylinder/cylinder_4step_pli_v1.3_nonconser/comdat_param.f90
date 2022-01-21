	MODULE ComDat_Param                                  

	IMPLICIT NONE                                         
	SAVE                                                  

	INTEGER, PARAMETER :: NI	= 80
	INTEGER, PARAMETER :: NJ	= 81
	INTEGER, PARAMETER :: NK	= 3

!	INTEGER, PARAMETER :: IMON	= 100
!	INTEGER, PARAMETER :: JMON	= 100
	INTEGER, PARAMETER :: KMON	=   2

	INTEGER, PARAMETER :: NSWPU=3
	INTEGER, PARAMETER :: NSWPV=3
	INTEGER, PARAMETER :: NSWPW=3
	INTEGER, PARAMETER :: NSWPP=5
	INTEGER, PARAMETER :: NSWPT=3

	DOUBLE PRECISION, PARAMETER :: RESIMAX = 1.0D-5
    DOUBLE PRECISION, PARAMETER :: CFL_USE = 0.25D0
    
	DOUBLE PRECISION, PARAMETER :: UIN     = 1.0D0
	DOUBLE PRECISION, PARAMETER :: VIN     = 0.0D0
	DOUBLE PRECISION, PARAMETER :: TIN     = 0.0D0
    DOUBLE PRECISION, PARAMETER :: T_WALL  = 1.0D0
    
!    DOUBLE PRECISION, PARAMETER :: RE      = 100.D0 ! REYNOLDS NUMBER
!    DOUBLE PRECISION, PARAMETER :: PR      = 1.D0   ! PRANDTL  NUMBER
!    DOUBLE PRECISION, PARAMETER :: PE      = 100.D0 ! PECLET   NUMBER
!    DOUBLE PRECISION, PARAMETER :: DTIME   = 1.0D-3
!    DOUBLE PRECISION, PARAMETER :: VEL_ROT = 0.0D0

    INTEGER :: NITER_W   = 1E6
    INTEGER :: RTIME_MAX = 5E2
    INTEGER :: IRESTART
    
    DOUBLE PRECISION :: RE
    DOUBLE PRECISION :: PR
    DOUBLE PRECISION :: PE
    DOUBLE PRECISION :: DTIME
    DOUBLE PRECISION :: VEL_ROT
        
	END MODULE ComDat_Param