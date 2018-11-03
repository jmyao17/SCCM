      MODULE VAPHFB_PAR
      implicit none


      integer, parameter :: lin=3
      character(500) :: INT_DIR
      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
      Type HOBasis
      real*8   b_osc,hw,hb0 
      integer  emax,tmax,nmax,lmax,twojmax,jmaxp5,ljmax,nljmax,NLEV,NLEV2      
      integer  NOrbit(0:1),tnljmax
      character*21, dimension(:), allocatable :: tts 
      character*1  tp(2),tis(2),tl(0:30) 
      End Type HOBasis
      Type(HOBasis) HO

      Integer, parameter :: lou = 24
      INTEGER, PARAMETER :: nosmax=100
      TYPE FileName
        character(len=100) :: NN,p1p2,f1b,f2bm,f2bJ,f2bJTMT,f3bJ,f3b22bJTMT
        character(len=100) :: IMSRG_Hme1b,IMSRG_Hme2b
        CHARACTER*100 cwf_ini(nosmax),cwf_fin(nosmax)
        CHARACTER*100 cwf(nosmax)
        CHARACTER*100 wf,out,FF 
        character*120 wf1,wf2,elem,TD1B,Rho1B,Rho2B,Rho3B
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
        INTEGER, dimension(:), allocatable :: nlj,t,n,l,twoj,twom,lj  ! here n=1,2,3,..
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
          integer ide(0:1),ide0(0:1),ide1(0:1)
      End Type HFBWFS
        Type(HFBWFS) HFB    

       Type ME1B
         real*8, dimension(:), allocatable :: ANme,AZme,zANme,zAZme
         real*8, dimension(:), allocatable :: r2,Q2_2t,Q2_1t,Q20t,Q21t,Q22t,Q20m,Q22m,Q2_2m,Q21m,Q2_1m 
         real*8, dimension(:), allocatable :: AJX_ME,AJY_ME,AJZ_ME
         real*8, dimension(:), allocatable :: AJX2_ME,AJY2_ME,AJZ2_ME
!         complex*16, dimension(:), allocatable :: zQ2_2t,zQ2_1t,zQ20t,zQ21t,zQ22t,zQ20m,zQ22m,zQ2_2m,zQ21m,zQ2_1m 
         complex*16, dimension(:), allocatable :: ZJXME,ZJYME,ZJZME
         complex*16, dimension(:), allocatable :: ZJX2ME,ZJY2ME,ZJZ2ME
         real*8, dimension(:), allocatable :: P_00_1m1_me,P_00_1p1_me,P_1m1_00_me,P_1p1_00_me, &
     &                                          P_00_10_me,P_10_00_me
       End Type ME1B
       Type(ME1B) cME1B

!      ........ reduced  matrix elements
       TYPE RedME1B
         real*8, dimension(:,:,:), allocatable :: T !T0,T2  
       END TYPE RedME1B 
       Type(RedME1B) rME1B


      integer Num_Rho1B  
!   ................... Input
      Type ReaderFile
       character(len=:), allocatable :: cIntID,cValID,cFLow,ctpp
       character*6 chwHO
       character*12 vs4me3b
       Integer   ICons,KMix,IpMix,IsHFB,InME3B,NPMix,IntType,IntJT,IntIMSRG,Isospin,itermax,inwf,iscale
       Integer   ihwHO,iCOM,nprot,nneut
       real*8    en,ep,Cscale,eta,tolcons,tolgrad
       integer   NLayer0,NLayer1,NGCM,Nuphi_ini,Nuphi_fin,nnnini,nnnfin
       integer   nq0i,nq0f,nq1i,nq1f
       Integer   iTD_Calc
       Integer   iRho3B,idens,iE2,icr 
      End Type ReaderFile
      Type(ReaderFile) Input 

      Type NucleusInfo
        character*2 nucnam
        Integer nucleon(0:2),ncore(0:2)
      End Type NucleusInfo
      Type(NucleusInfo) Nucl

      Type ConstPar
        real*8, dimension(:), allocatable :: hw_mesh,Etot,Q20t_mesh,Q22t_mesh,beta2t_mesh,gamma2t_mesh,P00_mesh
        Integer iQB,Ncons,iQBType,NQ
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
      real*8       NZmesh,Nmesh(400),Zmesh(400),phi(1:2)
      complex*16   zemiphi(1:2)
      Integer      NFOM,MPhi          ! number of guage angles for PNP
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
      Integer, PARAMETER :: NOBS=19
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
      INTEGER mit(-1:1),mip(-1:1)


!     .......... CG coefficient
      real*8, dimension(:,:,:,:,:,:), allocatable :: CG_Save
!     ................. basis for two-particle
!     eMax08: bMax=102 (temp.=106), aMaxMax=648 (see fort.100)

      INTEGER, PARAMETER :: aMaxMax= 650 !700 ! 1000 ! 500 ! 500 !1000
      INTEGER, PARAMETER :: bMax= 107 !140 !150 !100 !150
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
!        INTEGER, dimension (0:ljmax,0:nmax,0:1)   :: nlj
        INTEGER, dimension (:,:,:), allocatable   :: nlj
        INTEGER                                      nlj_vmax
      END Type VS_nlj

      TYPE VS_tnlj
        INTEGER, dimension (:), allocatable       :: t,n,l,twoj,lj,VType  ! here n=0,1,2,3,...
!        INTEGER, dimension (0:ljmax,0:nmax,0:1)   :: tnlj
        INTEGER, dimension (:,:,:), allocatable   :: tnlj
        INTEGER                                      tnlj_vmax
      END Type VS_tnlj

      TYPE VS_tnljm
        INTEGER, dimension (:), allocatable                  :: t,n,l,twoj,twom,lj  ! here n=1,2,3,..
!        INTEGER, dimension (-jmax2:jmax2,0:ljmax,0:nmax,0:1) :: tnljm
        INTEGER, dimension (:,:,:,:), allocatable :: tnljm
        INTEGER                                                 tnljm_vmax
      END Type VS_tnljm

!     .................. cost of lots of memory
      Type Struct_Rho2B
!          complex*16 ZRho2BJJ  (-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)
!          complex*16 ZTBDJJ    (-aMaxMax:aMaxMax,-aMaxMax:aMaxMax,0:bMax)
           complex*16, dimension (:,:,:), allocatable ::  ZTBDJJ ,ZRho2BJJ 
      END Type Struct_Rho2B
      Type(Struct_Rho2B) :: Dens

      Type Struct_ME3J
           TYPE(VS_nlj)    :: Vnlj
           TYPE(VS_tnlj)   :: VSPB
           TYPE(VS_tnljm)  :: Vtnljm
           TYPE(TPBC)      :: VTPB
           INTEGER, dimension (:), allocatable :: J12,P12,t12,a12,t3,n3,lj3,J123
!           INTEGER                             idx(0:bMax,1:aMaxMax,1:tnljmax,0:JMax2*3)
           INTEGER, dimension(:,:,:,:), allocatable :: idx
           INTEGER                             idx123_vmax
      END Type Struct_ME3J
      TYPE(Struct_ME3J) :: Rho3B

      
      Integer, PARAMETER :: JJmax=6
      Integer, PARAMETER :: JJmaxp=8
      real*8, dimension(0:JJmax) :: weightJ 
      TYPE TypeKernel
           INTEGER :: Nelem
           INTEGER, dimension(1:2000) ::  Nphi1,Nphi2
         complex*16, dimension (0:JJmax,-JJmax:JJmax,-JJmax:JJmax) :: njkk,hjkk,nn,zz,J2,&
     &                                             qp,qcp,q0p,q0pc 
         ! qpred(J, K1, K2) denoting <J1, K1, q1 || Q2 || J2, K2, q2>
         !complex*16, dimension (0:JJmax,0:JJmax,-JJmaxp:JJmaxp,-JJmax:JJmax) :: qpred,qcpred
         complex*16, dimension (0:JJmax,-JJmaxp-2:JJmaxp+2,-JJmax-2:JJmax+2,-2:2) :: qpred,qcpred
      END Type TypeKernel
      Type(TypeKernel) Kernel

      END MODULE VAPHFB_PAR
