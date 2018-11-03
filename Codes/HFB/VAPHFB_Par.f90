      MODULE VAPHFB_PAR
!     .............................................................
!     This module is mainly for the definition of global variables.
!     .............................................................
      implicit none

      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

      TYPE Optimizer
         character*10 method
         real*8 lr,mp,eps,bet1,bet2,gam1,gam2
      END TYPE Optimizer
      Type(Optimizer) Opt

      TYPE Shared_Var
           real*8 cf
      END Type Shared_Var
      TYPE(Shared_Var) shared


      character(500) :: INT_DIR
      integer, parameter :: lin=3 
!     .............. HO basis
      Type HOBasis
      real*8   b_osc,hw,hb0 
      integer  emax,tmax,nmax,lmax,twojmax,jmaxp5,ljmax,tnljmax,nljmax,NLEV,NLEV2      
      integer  NOrbit(0:1)
      character*21, dimension(:), allocatable :: tts 
      character*2 tp(2)
      character*1 tis(2),tl(0:30) 
      End Type HOBasis
      Type(HOBasis) HO

!     ............. Files
      TYPE FileName
        character(len=100) :: NN,p1p2,wf,out,f1b,f2bJ,f2bJTMT,f2bm,f3bJ,f3b22bJTMT
        character(len=200) :: IMSRG_Hme1b,IMSRG_Hme2b
      END Type FileName
      TYPE(FileName) File


      real*8, dimension(:,:,:,:,:,:), allocatable :: CG_Save
!    ..................... model space
      TYPE SPB_nlj
        INTEGER, dimension(:), allocatable :: e,n,l,twoj,lj  ! here n=0,1,2,3,...
        INTEGER, dimension(:,:), allocatable :: nlj
      END Type SPB_nlj
      TYPE(SPB_nlj) SPB

!   ......................in m-scheme 
      TYPE SPB_tnljm
        INTEGER, dimension(:), allocatable :: nlj,t,n,l,twoj,twom,lj  ! here n=1,2,3,..
        INTEGER, dimension(:,:,:,:), allocatable :: level ! (-jmax2:jmax2,0:ljmax,0:nmax,0:1) :: level
      END Type SPB_tnljm
        TYPE(SPB_tnljm) tnljm

!   .................... Hamiltonian
      TYPE Hamilton
          integer      iden,Aref
          integer      iabcd_max
          real*8       E0,ddd_kin,ddd_tbme
          INTEGER, dimension (:), allocatable :: ka,kb,kc,kd
          real*8, dimension(:,:), allocatable :: ME1BM
          real*8, dimension(:), allocatable   :: ME2BM
      END TYPE Hamilton
        TYPE(Hamilton) H

!    .................. HFB wave functions (U, V) and densities (rho, kappa)
      Type HFBWFS
          real*8, DIMENSION(:,:), allocatable :: U0,V0,U1,V1,GRAD_0,Thoul_0_work,ALcholes
          real*8, DIMENSION(:),   allocatable :: Thoul_0,Thoul_old
          real*8, DIMENSION(:,:), allocatable :: Ut,Vt
          real*8, DIMENSION(:,:), allocatable :: RO_0,Akapa10_0,Akapa01_0 
          real*8, DIMENSION(:,:), allocatable :: gamma_0,delta10,ham_0 ! HF field; delta10=Delta
          real*8, DIMENSION(:,:), allocatable :: it,ip
          integer iconv,ide(0:1)
          real*8  EFermi(0:1),alpha,Et 
      End Type HFBWFS
        Type(HFBWFS) HFB    

!      ............... 1B matrix elements (for constraint calculations) 
       Type ME1B
         real*8, dimension(:), allocatable :: Q2_2t,Q2_1t,Q20t,Q21t,Q22t,Q20m,Q22m,Q2_2m,Q21m,Q2_1m 
         real*8, dimension(:), allocatable :: Q40t,Q40m
         real*8, dimension(:), allocatable :: r2,Q30t
         real*8, dimension(:), allocatable :: AJX_ME,AJY_ME,AJZ_ME
         real*8, dimension(:), allocatable :: AJX2_ME,AJY2_ME,AJZ2_ME
         real*8, dimension(:), allocatable :: P_00_1m1_me,P_00_1p1_me,P_1m1_00_me,P_1p1_00_me, &
     &                                        P_00_10_me,P_10_00_me
       End Type ME1B
       Type(ME1B) cME1B
!      .............. for constraint calculations
      Type ConstPar
        real*8, dimension(:), allocatable :: hw,Etot,Q20t_mesh,Q22t_mesh,beta2t_mesh,gamma2t_mesh,P00_mesh
        Integer iQB,Ncons,iQBType,NQ
        real*8  Q2BA,Q2BP,Q2BN
      End Type ConstPar
      Type(ConstPar) Const

!     ............... Input Parameters
      Type ReaderFile
       character(len=:), allocatable :: cIntID,cValID,cFLow,ctpp
       character*6 chwHO
       Integer   KMix,NPMix,IsHFB,InME3B,IPMix,IntType,IntJT,IntIMSRG,Isospin,itermax,inwf,iscale
       Integer   ihwHO,iCOM,nprot,nneut
       real*8    alpha,eta,tolcons,tolgrad,bkick
       Integer   iRho3B,idens,iE2 
      End Type ReaderFile
      Type(ReaderFile) Input 

!     ............... Information about the nucleus of concerned
      Type NucleusInfo
        character*2 nucnam
        Integer nucleon(0:2),ncore(0:2)
      End Type NucleusInfo
      Type(NucleusInfo) Nucl

!     .................. PNP
      TYPE TypePNP
      real*8   phi_n, phi_p
      integer  NFOM
      END Type TypePNP
      TYPE(TypePNP) PNP
!     ................ scale factor for the matrix, the pfaffian of which is to be calc.ed
      Integer   iscale
      Real*8    Cscale


      Integer, parameter :: NCONSMAX=19 !15
      Integer iq2_const,iJ_const,ipair_const,ind_con
      Integer Icons_y_n(NCONSMAX)  ! 1 or 0
      real*8  qtc20,ACONS_MV_aux(NCONSMAX),ACONS_MV(NCONSMAX)
      real*8  beta2,gamma2,beta2_IV,gamma2_IV
      real*8, dimension(:), allocatable :: ANme,AZme,ACONS_ME !(NCONSMAX*NLEV*NLEV)
      real*8, dimension(:,:), allocatable :: t3me
      real*8, dimension(1:NCONSMAX) :: alag_0,alag !(NCONSMAX)

!   ............... constant parameter set
      Real*8,  PARAMETER :: CHOP = 1.d-9
      Real*8,  PARAMETER :: one = 1.d0
      Real*8,  PARAMETER :: two = 2.d0
      Real*8,  PARAMETER :: half= 0.5d0
      Real*8,  PARAMETER :: onem=-1.d0
      Real*8,  PARAMETER :: zero=0.d0
      Real*8,  PARAMETER :: third=1.d0/3.d0
      Real*8,  PARAMETER :: pi=4.d0*datan(1.d0)
      Integer, PARAMETER :: NOBS=10
      Integer, PARAMETER :: NOBS_PAV=15
      Integer, PARAMETER :: iCS = 1
      Real*8,  PARAMETER :: r0 = 1.2d0
      Integer, PARAMETER :: Niter_cons=40 !100

      Real*8, parameter ::  amu = 939.d0
      Real*8, parameter ::  hbc = 197.328284d0
      complex*16,PARAMETER :: zone=(1.d0,0.d0)
      complex*16,PARAMETER :: zzero=(0.d0,0.d0)
      complex*16,PARAMETER :: zimag=(0.d0,1.d0)
!     ................................. for MATH_LIB
      Integer, PARAMETER :: igfv  =1000
      INTEGER iv(-igfv:igfv),ibc(0:igfv,0:igfv)
      Real*8  sq(0:igfv),sqi(0:igfv),sqh(0:igfv),shi(0:igfv)
      Real*8  fak(-igfv:igfv),fad(0:igfv),fi(0:igfv),fdi(0:igfv),wf(0:igfv),wfi(0:igfv)
      Real*8  wfd(0:igfv),gm2(0:igfv),gmi(0:igfv),wg(0:igfv),wgi(0:igfv)
      INTEGER mit(-1:1),mip(-1:1)

      END MODULE VAPHFB_PAR
