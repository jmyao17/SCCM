      MODULE VAPHFB_PAR
      implicit none

      Integer, PARAMETER :: JJmax=6
      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
      Type HOBasis
      real*8   b_osc,hw,hb0 
      integer  emax,tmax,nmax,lmax,twojmax,jmaxp5,ljmax,nljmax,NLEV,NLEV2      
      integer  NOrbit(0:1),tnljmax
      character*14, dimension(:), allocatable :: tts 
      character*1  tp(2),tis(2),tl(0:30) 
      End Type HOBasis
      Type(HOBasis) HO

      Integer, parameter :: lou = 24
      INTEGER, PARAMETER :: Nuphimax=100
      TYPE FileName
        character(len=200) :: wf,out,FF
        CHARACTER*200 cwf_ini(Nuphimax),cwf_fin(Nuphimax)
        CHARACTER*200 cwf(Nuphimax)
        character*200      TD1B,wf1,wf2,elem,Rho1B,Rho2B,Rho3B
      END Type FileName
      TYPE(FileName) File

!    ..................... for model space
      TYPE SPB_nlj
        INTEGER, dimension(:),   allocatable :: e,n,l,twoj,lj  ! here n=0,1,2,3,...
        INTEGER, dimension(:,:), allocatable :: nlj         ! (0:nmax,0:ljmax)
        INTEGER, dimension(:,:,:), allocatable :: VType, CType ! (0:ljmax,0:nmax,0:1) :: VType
      END Type SPB_nlj
        TYPE(SPB_nlj) SPB

!   ......................in m-scheme 
      TYPE SPB_tnljm
        INTEGER, dimension(:), allocatable :: t,n,l,twoj,twom,lj  ! here n=1,2,3,..
        INTEGER, dimension(:,:,:,:), allocatable :: level ! (-jmax2:jmax2,0:ljmax,0:nmax,0:1) :: level
      END Type SPB_tnljm
        TYPE(SPB_tnljm) tnljm

      TYPE Hamilton
          integer      iden,Aref
          integer      iabcd_max
          real*8       E0,ddd_kin,ddd_tbme
          INTEGER, dimension (:), allocatable ::   ka,kb,kc,kd
          real*8, dimension(:,:), allocatable :: ME1BM
          real*8, dimension(:), allocatable   :: ME2BM
          complex*16, dimension(:,:), allocatable :: zME1BM
      END TYPE Hamilton
        TYPE(Hamilton) H

      Type HFBWFS
          real*8, DIMENSION(:,:), allocatable :: U0,V0,U1,V1,GRAD_0,Thoul_0_work,ALcholes
          real*8, DIMENSION(:),   allocatable :: Thoul_0
          real*8, DIMENSION(:,:), allocatable :: gamma_0,ham_0 ! HF field
          real*8, DIMENSION(:,:), allocatable :: RO_2,Akapa10_2,Akapa01_2 
          real*8, DIMENSION(:,:), allocatable :: RO_1,Akapa10_1,Akapa01_1 
          real*8, DIMENSION(:), allocatable :: it1,it2 
          integer, DIMENSION(:), allocatable :: ip1,ip2 
          real*8, DIMENSION(:),   allocatable :: vv1,vv2 
          real*8, DIMENSION(:),   allocatable :: kocc1,kocc2 
      End Type HFBWFS
        Type(HFBWFS) HFB    


!      ........ reduced  matrix elements
       TYPE RedME1B
         real*8, dimension(:,:,:), allocatable :: T !T0,T2  
       END TYPE RedME1B
       Type(RedME1B) rME1B
!      ............... cost lots of memory
       Type ME1B
         real*8, dimension(:), allocatable :: ANme,AZme,zANme,zAZme
         real*8, dimension(:), allocatable :: Q2_2t,Q2_1t,Q20t,Q21t,Q22t,Q20m,Q22m,Q2_2m,Q21m,Q2_1m 
         real*8, dimension(:), allocatable :: AJX_ME,AJY_ME,AJZ_ME
         real*8, dimension(:), allocatable :: AJX2_ME,AJY2_ME,AJZ2_ME
         complex*16, dimension(:), allocatable :: ZJXME,ZJYME,ZJZME
         complex*16, dimension(:), allocatable :: ZJX2ME,ZJY2ME,ZJZ2ME
         real*8, dimension(:), allocatable :: P_00_1m1_me,P_00_1p1_me,P_1m1_00_me,P_1p1_00_me, &
     &                                          P_00_10_me,P_10_00_me
       End Type ME1B
       Type(ME1B) cME1B

      integer Num_Rho1B  
!   ................... Input
      Type ReaderFile
       character(len=:), allocatable :: cIntID,cValID,cFlow
       character*12 vs4me3b
       Integer   NPMix,IntType,IntJT,IntIMSRG,Isospin,itermax,inwf,iscale
       Integer   ihwHO,IsHFB,iCOM,nprot,nneut
       real*8    en,ep,Cscale,eta,tolcons,tolgrad
       integer   NGCM,Nuphi_ini,Nuphi_fin,nnnini,nnnfin
       integer   nq0i,nq0f,nq1i,nq1f
       Integer   iTD_Calc,IE_Conv,NOS,icr
       Integer   iRho3B,idens,J0k 
       Integer   NOSJ(0:JJmax) 
       End Type ReaderFile 
       Type(ReaderFile) Input 

      Type NucleusInfo
        character*2 nucnam
        Integer nucleon(0:2),ncore(0:2)
      End Type NucleusInfo
      Type(NucleusInfo) Nucl

      Type ConstPar
        real*8, dimension(:), allocatable :: hw_mesh,Etot,Q20t_mesh,Q22t_mesh,beta2t_mesh,gamma2t_mesh,P00_mesh
        Integer chi,iQB,Ncons,iQBType,NQ
        real*8  Q2BA,Q2BP,Q2BN
      End Type ConstPar
      Type(ConstPar) Const

      Integer, parameter :: NCONSMAX=17 !15
      Integer Icons_y_n(NCONSMAX)
      Integer ipair_const,ind_con
      real*8  ACONS_MV_aux(NCONSMAX),ACONS_MV(NCONSMAX)
      real*8  beta2,gamma2,beta2_IV,gamma2_IV
      real*8, dimension(:), allocatable :: ACONS_ME !(NCONSMAX*NLEV*NLEV)
      real*8, dimension(:,:), allocatable :: t3me 
      integer, dimension(:,:), allocatable :: ipme 
      real*8, dimension(1:NCONSMAX) :: alag_0,alag !(NCONSMAX)

!    ......................... projection
!    .......... AMP
      TYPE TypeAMP
      integer   NNNP,NLEG_ALP, NLEG_BET, NLEG_GAM
      integer   IsAMP,i3DAMP
      INTEGER   JJMax   ! Maximum number of angular momentum
      real*8, dimension(:), allocatable :: alpha,beta,gamma,cosbet
      real*8, dimension(:), allocatable :: wwalp,wwbet,wwgam
      complex*16, dimension(:,:), allocatable :: ZROT_m1m2 
      END Type TypeAMP
      TYPE(TypeAMP) AMP

      real*8,dimension(:,:,:,:,:),allocatable ::  sdjmk

!     ........................
      TYPE TypePNP
      real*8       phi(1:2)
      complex*16   zemiphi(1:2)
      Integer      MPhi,NFOM          ! number of guage angles for PNP
      real*8       eps,eps1,eps2           ! cutoff 
      integer      npz           ! dim. of trunacted s.p. space
      END Type TypePNP
        TYPE(TypePNP) PNP
!     ................ scale factor for the matrix, the pfaffian of which is to be calc.ed
      Integer   iscale
      Real*8    Cscale

!     .................................... parameter for interation


!   ............... constant parameter set
      Real*8,  PARAMETER :: CHOP = 1.d-8
      Real*8,  PARAMETER :: one = 1.d0
      Real*8,  PARAMETER :: two = 2.d0
      Real*8,  PARAMETER :: half= 0.5d0
      Real*8,  PARAMETER :: onem=-1.d0
      Real*8,  PARAMETER :: zero=0.d0
      Real*8,  PARAMETER :: third=1.d0/3.d0
      Real*8,  PARAMETER :: pi=4.d0*datan(1.d0)
      Integer, PARAMETER :: NOBS=18
      Integer, PARAMETER :: NOBS_PAV=15
      Integer, PARAMETER :: iCS = 1
      Real*8,  PARAMETER :: r0 = 1.2d0
      Integer, PARAMETER :: Niter_cons=100

      Real*8, parameter ::  amu = 939.d0
      Real*8, parameter ::  hbc = 197.328284d0
      complex*16,PARAMETER :: zone=(1.d0,0.d0)
      complex*16,PARAMETER :: zzero=(0.d0,0.d0)
      complex*16,PARAMETER :: zimag=(0.d0,1.d0)
      complex*16,PARAMETER :: ztwo=(2.d0,0.d0)
!     ................................. for MATH_LIB
      Integer, PARAMETER :: igfv  =1000
      Integer, PARAMETER :: igfvbc=1000
      INTEGER iv(-igfv:igfv),ibc(0:igfvbc,0:igfvbc)
      Real*8  sq(0:igfv),sqi(0:igfv),sqh(0:igfv),shi(0:igfv)
      Real*8  fak(0:igfv),fad(0:igfv),fi(0:igfv),fdi(0:igfv),wf(0:igfv),wfi(0:igfv)
      Real*8  wfd(0:igfv),gm2(0:igfv),gmi(0:igfv),wg(0:igfv),wgi(0:igfv)
      INTEGER mip(-1:1)


!     .......... CG coefficient
!      real*8, dimension(:,:,:,:,:,:), allocatable :: CG_Save
!     ................. basis for two-particle
      INTEGER, PARAMETER :: aMaxMax=800 !1000
      INTEGER, PARAMETER :: bMax=100 ! 100 !150
      TYPE TPBC
          integer,dimension(:,:,:), allocatable :: a
          INTEGER,dimension(:,:,:), allocatable :: block !(0:JMax2,0:1,0:2)  ! JJ, PP, TT 
          INTEGER   J12(0:bMax,1:aMaxMax)
          INTEGER   P12(0:bMax,1:aMaxMax)
          INTEGER   n1(0:bMax,1:aMaxMax)   ! n=1,2,3,...
          INTEGER   t1(0:bMax,1:aMaxMax)
!          INTEGER   l1(0:bMax,1:aMaxMax)
          INTEGER   twol1(0:bMax,1:aMaxMax)
          INTEGER   j1(0:bMax,1:aMaxMax)
          INTEGER   n2(0:bMax,1:aMaxMax)
          INTEGER   t2(0:bMax,1:aMaxMax)
!          INTEGER   l2(0:bMax,1:aMaxMax)
          INTEGER   twol2(0:bMax,1:aMaxMax)
          INTEGER   j2(0:bMax,1:aMaxMax)
          INTEGER   iV(0:bMax,1:aMaxMax)
          INTEGER   aMax(0:bMax)
          INTEGER   cMax(0:bMax)
          INTEGER   btt(0:bMax)
          INTEGER   bPP(0:bMax)
          INTEGER   bJJ(0:bMax)
          INTEGER   bMax
          INTEGER   bJMax
      END TYPE TPBC
      TYPE(TPBC) :: TPB !(0:bMax,1:aMaxMax)      

!    ..........................
      TYPE Full_SPB_tnlj
!        INTEGER, dimension (0:nljmax*2)         :: t,e,n,l,twoj,lj  ! here n=1,2,3,...
!        INTEGER, dimension (0:ljmax,0:nmax,0:1) :: tnlj
        INTEGER, dimension (:), allocatable     :: t,e,n,l,twoj,lj  ! here n=1,2,3,...
        INTEGER, dimension (:,:,:), allocatable :: tnlj
        INTEGER                                 :: tnlj_fmax
      END Type Full_SPB_tnlj
        TYPE(Full_SPB_tnlj) FSPB

!   ............................................ 
!   structure for ME3B in valence space
!   ...........................................
      TYPE VS_nlj
        INTEGER, dimension (:), allocatable       :: n,l,twoj,lj ! here n=0,1,2,3,...
        INTEGER, dimension (:,:,:), allocatable   :: nlj
        INTEGER                                      nlj_vmax
      END Type VS_nlj

      TYPE VS_tnlj
        INTEGER, dimension (:), allocatable       :: t,n,l,twoj,lj,VType  ! here n=0,1,2,3,...
        INTEGER, dimension (:,:,:), allocatable   :: tnlj
        INTEGER                                      tnlj_vmax
      END Type VS_tnlj

      TYPE VS_tnljm
        INTEGER, dimension (:), allocatable                  :: t,n,l,twoj,twom,lj  ! here n=1,2,3,..
        INTEGER, dimension (:,:,:,:), allocatable :: tnljm
        INTEGER                                                 tnljm_vmax
      END Type VS_tnljm

      Type Struct_ME3J
           TYPE(VS_nlj)    :: Vnlj
           TYPE(VS_tnlj)   :: VSPB
           TYPE(VS_tnljm)  :: Vtnljm
           TYPE(TPBC)      :: VTPB
           INTEGER, dimension (:), allocatable :: J12,P12,t12,a12,t3,n3,lj3,J123
           INTEGER, dimension(:,:,:,:), allocatable :: idx
           INTEGER                             idx123_vmax
      END Type Struct_ME3J
      TYPE(Struct_ME3J) :: Rho3B

      
!    ......... for GCM
      Integer, PARAMETER :: NOQmax=100
      real*8, dimension(0:JJmax) :: weightJ 
      TYPE TypeKernel
         complex*16, dimension (1:NOQmax,1:NOQmax, 0:JJmax,-JJmax:JJmax,-JJmax:JJmax) :: njkk,hjkk,nn,zz,J2,&
     &                                             qp,qcp,q0p,q0pc 
         complex*16, dimension (1:NOQmax,1:NOQmax, 0:JJmax,-JJmax-2:JJmax+2,-JJmax-2:JJmax+2,-2:2) :: qpred,qcpred   
      END Type TypeKernel
      Type(TypeKernel) Kernel

      Type Struct_GCM
           integer kmax,NOQ, Jmax, Jmin, Jdf, nmaxdi(0:JJmax)
           integer qkmax(0:JJmax),ik(1:NOQmax*(JJmax*2+1),0:JJmax),iq(1:NOQmax*(JJmax*2+1),0:JJmax)
           integer NOS_J(1:NOQmax*JJmax,0:JJmax)
!           real*8, dimension (1:NOQmax*JJmax,1:NOQmax*JJmax,0:JJmax) :: njkkqq,hjkkqq
!          ... bug here
           real*8, dimension (0:JJmax,1:NOQmax*(JJmax*2+1),1:NOQmax*(JJmax*2+1)) :: njkkqq,hjkkqq
           real*8, dimension (1:NOQmax*(JJmax*2+1),0:JJmax,10) :: GJkq    
           real*8, dimension (1:NOQmax*(JJmax*2+1),0:JJmax,10) :: FqJk
           real*8, dimension (0:JJmax,10) :: EJk,mean,variance 
           real*8  eps
           integer, dimension(1:NOQmax*JJmax,1:NOQmax*(JJmax*2+1)) :: Iexst
!           integer, dimension (:,:), allocatable :: Iexst
      END Type Struct_GCM
      Type(Struct_GCM) GCM

      Type Struct_GCMDens
           real*8, dimension(:,:,:,:,:,:,:), allocatable :: TD1BSum
           real*8, dimension(:,:,:,:), allocatable :: Rho1BSum
           real*8, dimension(:,:,:), allocatable :: Rho2BSum
           real*8, DIMENSION(:,:), allocatable :: L3BSum,Rho3BSum  
      END Type Struct_GCMDens
      Type(Struct_GCMDens) GCMDens


      END MODULE VAPHFB_PAR
