#!/bin/bash
#PBS -l walltime=20:00:00
#PBS -l nodes=5:ppn=1,walltime=28:50:00,mem=20gb      ###nodes:ppn - how many nodes & cores per node (ppn)
#PBS -N Ti48-N8
#PBS -o Ti48-N8

# ........................ script for HFB calculation
N=4
itemax=300
ZZ=2 #32
NN=6 #44
Nucl=Ti48
iwf=2           # 2 or 1
IntJT=0      # 0 (in M-shceme); 1 (J2M); 2 (J2JT2M) 
Flow=s000 #s999
IntID=KB3G #${flow} #chi23bCa48_srg0625
Hs=H0
HFBdir=MPI_CHF2 #MPI_HFB7 
hw=11
eMax=eMax0${N}
lMax=lMax0${N}
hwHO=hwHO${hw}
Input=${Nucl}_magic_${eMax}.dat
ValID=pf #pf5g9

me1b=SM_KB3G_ME1B.dat 
me2b=SM_KB3G_ME2B.dat

PATH_WORK=./ 
PATH_INT=../${IntID}

cd ${PATH_WORK}



#cp ${GCM_ME_FILES}/${me1b} ${PATH_INT}/IMSRG_s000_SPE_${IntID}_${ValID}_${hwHO}.dat
#cp ${GCM_ME_FILES}/${me2b} ${PATH_INT}/IMSRG_${Flow}_${IntID}_${ValID}_J.dat
#cp ${GCM_ME_FILES}/${me2b} ${PATH_INT}/IMSRG_${Flow}_${IntID}_${ValID}_${hwHO}_J.dat


cat <<EOF > chi2b3b.int 
File4ME1B:    ${me1b}
File4ME2B:    ${me2b} 
EOF

cat <<EOF > hfb.dat 
NP-Mix    1                             ! np mixture: 0 (w/o); 1(w/)
IP-Mix    0                             ! parity mixture: 0 (w/o); 1(w/)
KMix      1                             ! Triaxiality/K-mixture: 0 (w/o); 1(w/)
IsHFB     1                             ! 0 (HF); 1(HFB)
IntType   0                             ! 0: shell-model; 1: Chiral
InME3B    0                             ! implicit ME3B: 0 (w/o); 1 (w/)
IntID     ${IntID}                      !
ValID     ${ValID}                        ! model space: emax0N 
hw        ${hw}                            ! i2,frequency of HO 
COM       2                             ! Center-of-Mass Correction: 1 (1B); 2 (1B+2B) 
FlowPars  ${Flow}                          ! s value
Int-JT    ${IntJT}                      ! 0(m-scheme); 1(JTTz-scheme); 2 (J2JT-scheme): 3(From Nathan,J2JT)
IntIMSRG  1                             ! IMSRG interation (1) or not (0)
IsoSpin   0                             ! (0) Int. does not have isospin symmetry; (1) Int. has isospin symm.
nprot     ${ZZ}                            ! neutron
nneut     ${NN}                            ! proton
Iter max  ${itemax}                           ! No. of iteractions 
Method    GDM                             ! GD: gradient descent method; GDM(with mom.); ADAM
eta1      7.e-3                         ! step size in gradient method
tolcons   1.e-4
tolgrad   1.e-4                         !  precision
INP_wf    ${iwf}                        ! (0)read wf from fort.10; (1) read from saved file; (2) from a random wf
PNP_phi   1                             !  1 (HFB) or 7 (VAPHFB)
........................................ constraint quantities
kick_step 0.0
B_G       1                            ! 0 (Q value); 1 (bet,gam)
NCONS     6
Q20t      1 0.d0                       ! 3: Q20 or beta
Q22t      1 0.d0                       ! 4: Q22 or gamma
Q21t      1 0.d0                       ! 5:
Q20m      0 0.d0                       ! 6: Q20m or beta
Q22m      1 0.d0                       ! 7: Q22m or gamma
Q21m      1 0.d0                       ! 8
JX        0 0.0d0                      ! 9  ! equivalent to choose diff. crank freq.
JY        0 0.d0                       ! 10
JZ        0 0.d0                       ! 11
P00_10    1 0.0d0                      ! 12  P_TMT_JMJ ! T=0, L=0, J=S=1, MJ=0
P1m1_00   0 0.d0                       ! 13
P1p1_00   0 0.d0                       ! 14
P10_00    0 0.d0                       ! 15
P00_1m1   0 0.0d0                      ! 16  P_TMT_JMJ ! T=0, L=0, J=S=1, MJ=-1
P00_1p1   0 0.0d0                      ! 17  P_TMT_JMJ ! T=0, L=0, J=S=1, MJ=+1
Q40t      0 0.d0                       ! 3: Q40t
Q40m      0 0.d0                       ! 3: Q40m
EOF

cat <<EOF > betgam.dat
 ............
 max
 bet     gam    P00(T=0,J=1)  crank.freq
 ----------------------------------------
 3
 0.0      0      0          0.0
 0.1      0      0          0.0
 0.2      0      0          0.0
EOF

#.... copy matrix elements from EOM-IMSRG calculation
cd ${PATH_WORK}

./HFB 

cat fort.9* >>E_beta_gam.dat
rm fort.*
