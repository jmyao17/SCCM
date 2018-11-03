        subroutine PNAMP_Default()
        USE VAPHFB_Par

        Num_Rho1B = 0
        Input%IntIMSRG  = 1
!    .................. PNP
        PNP%NFOM      = 7
!    .................. AMP
        AMP%NLEG_ALP  = 1
        AMP%NLEG_BET  = 8
        AMP%NLEG_GAM  = 1
!    .......... switches
       Input%iRho3B        = 0   ! with/without Rho3B
       Input%idens       = 0   ! calc or read densities

        return
        end
