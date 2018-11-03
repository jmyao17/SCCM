#!/bin/bash
#PBS -l walltime=20:00:00
#PBS -l nodes=5:ppn=1,walltime=28:50:00,mem=20gb      ###nodes:ppn - how many nodes & cores per node (ppn)
#PBS -N Ti48-N8
#PBS -o Ti48-N8

# ........................ script for HFB calculation
N=4
itemax=400  
ZZ=1
NN=7
Nucl=Sc48
iwf=2           # 2 or 1
IntJT=0      # Interaction from IMSRG(by Heiko): choose 4(T2M) or 0(M); 2(J  to JT); 1(JT to M); 3(Nathan)
Hs=H0 #Ca48_H1 #H1
Flow=s000 #s999 #s999
IntID=KB3 #-48Ca 
HFBdir=MPI_CHF2 
eMax=eMax0${N}
lMax=lMax0${N}
Input=${Nucl}_KB3G_${eMax}.dat
ValID=pf

PATH_SM=/mnt/research/imsrg/jmyao/Int/KB3G/${Hs}
PATH_ME=${PATH_magic}/${IntID}_${eMax}

export GCM_ME_FILES=${PATH_SM}

me1b=KB3G_ME1B.dat #EOM_Ca48_chi2b3b_srg0625_eMax04_hwHO020_unnorm_ham_1b.dat
me2b=KB3G_ME2B_JTMT.dat #EOM_Ca48_chi2b3b_srg0625_eMax04_hwHO020_unnorm_ham_2b.dat

PATH_WORK=/mnt/home/yaojiang/GCM/ABGCM/${HFBdir}
PATH_INT=/mnt/home/yaojiang/GCM/ABGCM/Int

cd ${PATH_WORK}

cp ${PATH_SM}/${me1b} ${PATH_INT}/IMSRG_${Flow}_SPE_${IntID}_${ValID}.dat

#cp ${PATH_SM}/${me2b} ${PATH_INT}/IMSRG_${Flow}_${IntID}_${ValID}_JTMT.dat
cp ${PATH_SM}/${me2b} ${PATH_INT}/IMSRG_${Flow}_${IntID}_${ValID}_J.dat

cat <<EOF > chi2b3b.int 
File4ME1B:    ${me1b}
File4ME2B:    ${me2b} 
EOF

cat <<EOF > hfb.dat
NP-Mix    0                             ! np mixture: 0 (w/o); 1(w/)
IP-Mix    0                             ! parity mixture: 0 (w/o); 1(w/)
KMix      1                             ! Triaxiality/K-mixture: 0 (w/o); 1(w/)
IsHFB     0                             ! 0 (HF); 1(HFB)
IntType   0                             ! 0: shell-model; 1: Chiral
IntID     ${IntID}                      !
ValID     ${ValID}                        ! model space: emax0N 
FlowPars  ${Flow}                          ! s value
Int-JT    ${IntJT}                      ! 0(m-scheme); 1(JTTz-scheme); 2 (J2JT-scheme): 3(From Nathan,J2JT)
IntIMSRG  1                             ! IMSRG interation (1) or not (0)
IsoSpin   0                             ! (0) Int. does not have isospin symmetry; (1) Int. has isospin symm.
nprot     ${ZZ}                            ! neutron
nneut     ${NN}                            ! proton
Iter max  ${itemax}                           ! No. of iteractions 
eta1      7.e-3                         ! step size in gradient method
alpha     0.0                           ! momentum parameter: gamma [0,1] in gradient method
tolcons   1.e-4
tolgrad   1.e-4                         !  precision
INP_wf    ${iwf}                        ! (0)read wf from fort.10; (1) read from saved file; (2) from a random wf
PNP_phi   1                             !  1 (HFB) or 7 (VAPHFB)
........................................ constraint quantities
kick_step 0.0                          !
B_G       1                            ! 0 (Q value); 1 (bet,gam)
NCONS     4
Q20t      1 0.d0                       ! 3: Q20 or beta
Q22t      1 0.d0                       ! 4: Q22 or gamma
Q21t      1 0.d0                       ! 5:
Q20m      0 0.d0                       ! 6: Q20m or beta
Q22m      0 0.d0                       ! 7: Q22m or gamma
Q21m      0 0.d0                       ! 8
JX        0 0.d0                       ! 9
JY        0 0.d0                       ! 10
JZ        0 0.d0                       ! 11
P00_10    1 0.0d0                      ! 12  P_TMT_JMJ ! T=0, L=0, J=S=1, MJ=0
P1m1_00   0 0.d0                       ! 13
P1p1_00   0 0.d0                       ! 14
P10_00    0 0.d0                       ! 15
P00_1m1   0 0.0d0                      ! 16  P_TMT_JMJ ! T=0, L=0, J=S=1, MJ=-1
P00_1p1   0 0.0d0                      ! 17  P_TMT_JMJ ! T=0, L=0, J=S=1, MJ=+1

EOF

cat <<EOF > betgam.dat
 ............
 max
 bet     gam    P00(T=0,J=1)
 ---------------------
 5
 0.2     0       0    0
 0.2     0       0    0.5
 0.2     0       0    1.0
 0.2     0       0    1.5
 0.2     0       0    2.0
 15
-0.35    0       0
-0.3     0       0
-0.25    0       0
-0.2     0       0
-0.15    0       0
-0.1     0       0
 0.0     0       0
 0.1     0       0
 0.15    0       0
 0.2     0       0
 0.25    0       0
 0.3     0       0
 0.35    0       0
 0.4     0       0
 0.5     0       0
EOF

#.... copy matrix elements from EOM-IMSRG calculation
cd ${PATH_WORK}

./HFB < hfb.dat # ${Input}
