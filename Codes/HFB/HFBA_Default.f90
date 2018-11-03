        subroutine Default()
        USE VAPHFB_Par

        Num_Rho1B = 0

        IntIMSRG  = 1
!    .................. PNP
        NFOM      = 1
!    .................. AMP
        NLEG_ALP  = 1
        NLEG_BET  = 8
        NLEG_GAM  = 1
!    .......... switches
       iME3B        = 0   ! with/without Rho3B
       Is_Read_Dens = 0   ! calc or read densities

        return
        end
