#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=08:50:00,mem=160gb      ###nodes:ppn - how many nodes & cores per node (ppn)
#PBS -N HFB4-Ge76-hw16 
#PBS -o Ge76-hw16-N8
# ........................ script for HFB calculation
N=6
itemax=800
hw=10
ZZ=4 #32
NN=16 #44
Nucl=Ge76
iwf=2           # 2 or 1
IntJT=0      # 0 (in M-shceme) or 1(J->M) or 2->1(J->JTMT->M) 
Flow=s000 #s999
IntID=GCN2850 #${flow} #chi23bCa48_srg0625
Hs=H0
#HFBdir=MPI_CHF2 #MPI_HFB7 
eMax=eMax0${N}
lMax=lMax0${N}
hwHO=hwHO${hw}
Input=${Nucl}_magic_${eMax}.dat
ValID=pf5g9

#PATH_SM=/mnt/research/imsrg/jmyao/Int #/GCN2850
#dir=${IntID}/${Hs}
#export GCM_ME_FILES=${PATH_SM}/$dir

me1b=SM_GCN2850_ME1B.dat #EOM_Ca48_chi2b3b_srg0625_eMax04_hwHO020_unnorm_ham_1b.dat
me2b=SM_GCN2850_ME2B.dat #EOM_Ca48_chi2b3b_srg0625_eMax04_hwHO020_unnorm_ham_2b.dat

#me1b=Ca48_chi2b3b400_srg0625J_eMax04_lMax04_hwHO020_evolved_unnorm_ham_1b.dat
#me2b=Ca48_chi2b3b400_srg0625J_eMax04_lMax04_hwHO020_evolved_unnorm_ham_2b.dat

#PATH_WORK=/global/homes/j/jmyao/GCM_Solver/CPCv5/${HFBdir}
#PATH_INT=/global/u2/j/jmyao/GCM_Solver/CPCv5/Int #/global/homes/j/jmyao/GCM_Solver/CPCv5/Int

PATH_WORK=./ #/mnt/home/yaojiang/GCM/ABGCM/${HFBdir}
PATH_INT=../${IntID}

cd ${PATH_WORK}

#cp ${GCM_ME_FILES}/${me1b} ${PATH_INT}/IMSRG_s000_SPE_${IntID}_${ValID}_${hwHO}.dat
#cp ${GCM_ME_FILES}/${me2b} ${PATH_INT}/IMSRG_${Flow}_${IntID}_${ValID}_J.dat
#cp ${GCM_ME_FILES}/${me2b} ${PATH_INT}/IMSRG_${Flow}_${IntID}_${ValID}_${hwHO}_J.dat

#echo ${PATH_INT}/IMSRG_s000_SPE_${IntID}_${ValID}.dat

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
Method    GD                             ! GD: gradient descent method; GDM(with mom.); ADAM
eta1      7.e-3                         ! step size in gradient method
tolcons   1.e-4
tolgrad   1.e-5                         !  precision
INP_wf    ${iwf}                        ! (0)read wf from fort.10; (1) read from saved file; (2) from a random wf
PNP_phi   1                             !  1 (HFB) or 7 (VAPHFB)
........................................ constraint quantities
kick_step 0.0
B_G       1                            ! 0 (Q value); 1 (bet,gam)
NCONS     7
Q20t      1 0.d0                       ! 3: Q20 or beta
Q22t      1 0.d0                       ! 4: Q22 or gamma
Q21t      1 0.d0                       ! 5:
Q20m      0 0.d0                       ! 6: Q20m or beta
Q22m      0 0.d0                       ! 7: Q22m or gamma
Q21m      0 0.d0                       ! 8
JX        1 0.d0                       ! 9
JY        1 0.d0                       ! 10
JZ        1 0.d0                       ! 11
P00_10    1 0.0d0                      ! 12  np isoscalar pairing ! T=0, L=0, J=S=1, MJ=0
P1m1_00   0 0.d0                       ! 13  pp isovector pairing
P1p1_00   0 0.d0                       ! 14  nn isovector pairing
P10_00    0 0.d0                       ! 15  np isovector pairing
P00_1m1   0 0.0d0                      ! 16  P_TMT_JMJ ! T=0, L=0, J=S=1, MJ=-1
P00_1p1   0 0.0d0                      ! 17  P_TMT_JMJ ! T=0, L=0, J=S=1, MJ=+1
Q40t      0 0.d0                       ! 18: Q40 
Q40m      0 0.d0                       ! 19: Q40 
EOF

nq=16 #17

cd ${PATH_WORK}

./HFB

#mpirun -n 1 ./HFB #>Ge76_MF_out
  
#mpirun -n ${nq} ./HFB  
