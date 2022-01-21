    SUBROUTINE POST_PROCESS
    
    USE COMDAT_SHARED
    
    IMPLICIT NONE
    SAVE
!---------------------------------------------------------------------------------------------
!   POST-PROCESSING FOR NU,CD,CL
    IF(IPOST.EQ.0) THEN
		IF(RANK.EQ.0) THEN
		OPEN(UNIT=12,FILE='POST_NU_CD_CL.PLT',STATUS='UNKNOWN')
        WRITE(12,*) 'variables=TIME,NU,CD,CL,CL_AMP,CPB,ST,,RMS_CD,RMS_CL'
        CLOSE(12)
        ENDIF
		RTIME_MAX  = RTIME_MAX * 2.D0
		NITER_POS = 0
        IRESTART = 1
        TOGGLE = 1
        IPOST = 1
        SIGN = .TRUE.
    ELSE IF(IPOST.EQ.1.AND.IRESTART.EQ.1) THEN
        RTIME_MAX  = RTIME_MAX + (1.D0/ST)*1.D0 ! 1 CYCLE
        NITER_W = (1.D0/ST)/DTIME/24 ! 24 FRAME
        IRESTART = 0
        SIGN = .TRUE.
    ELSE IF(IPOST.EQ.1.AND.IRESTART.EQ.0) THEN
        SIGN = .FALSE.
    ENDIF
!---------------------------------------------------------------------------------------------
    END SUBROUTINE POST_PROCESS