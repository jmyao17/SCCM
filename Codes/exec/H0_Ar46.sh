#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=18:50:00,mem=80gb      ###nodes:ppn - how many nodes & cores per node (ppn)
#PBS -N HFB2-Ar46-hw16-s000-1 
#PBS -o Ar46-hw16-s000-1

# ........................ script for HFB calculation
N=8
itemax=300
hw=16
ZZ=18
NN=28
Nucl=Ar46
iwf=1           # 2 or 1
IntJT=0      # 0 (in M-shceme) or choose 4(J->M) or 2->1(J->JTMT->M) 
IntID=magic-Ar46_srg0953 #chi23bCa48_srg0625
Hs=H0
HFBdir=MPI_CHF2 
Flow=s000
eMax=eMax0${N}
lMax=lMax0${N}
hwHO=hwHO${hw}
Input=${Nucl}_magic_${eMax}.dat

PATH_magic=/mnt/research/imsrg/jmyao/Int/magic
PATH_ME=${PATH_magic}/${IntID}_${eMax}

dir=${IntID}_${eMax}/${Hs}/${hwHO}
export GCM_ME_FILES=${PATH_magic}/$dir

#echo ${PATH_magic}/$dir

me1b=em1.8-2.0_ME1B.dat #EOM_Ca48_chi2b3b_srg0625_eMax04_hwHO020_unnorm_ham_1b.dat
me2b=em1.8-2.0_ME2B.dat #EOM_Ca48_chi2b3b_srg0625_eMax04_hwHO020_unnorm_ham_2b.dat

#me1b=Ca48_chi2b3b400_srg0625J_eMax04_lMax04_hwHO020_evolved_unnorm_ham_1b.dat
#me2b=Ca48_chi2b3b400_srg0625J_eMax04_lMax04_hwHO020_evolved_unnorm_ham_2b.dat

#PATH_WORK=/global/homes/j/jmyao/GCM_Solver/CPCv5/${HFBdir}
#PATH_INT=/global/u2/j/jmyao/GCM_Solver/CPCv5/Int #/global/homes/j/jmyao/GCM_Solver/CPCv5/Int

PATH_WORK=/mnt/home/yaojiang/GCM/ABGCM/${HFBdir}
PATH_INT=/mnt/home/yaojiang/GCM/ABGCM/Int

cd ${PATH_WORK}

#cp ${PATH_ME}/${Hs}/${hwHO}/${me1b} ${PATH_INT}/IMSRG_${Flow}_SPE_${IntID}_eMax0${N}_${hwHO}.dat

#cp ${PATH_ME}/${Hs}/${hwHO}/${me2b} ${PATH_INT}/IMSRG_${Flow}_${IntID}_tpp_eMax0${N}_${hwHO}_J.dat
#echo ${PATH_INT}/IMSRG_${Flow}_SPE_${IntID}_eMax0${N}_${hwHO}.dat

cat <<EOF > chi2b3b.int 
File4ME1B:    ${me1b}
File4ME2B:    ${me2b} 
EOF

cat <<EOF > hfb.dat 
NP-Mix    0                             ! np mixture: 0 (w/o); 1(w/)
IP-Mix    0                             ! parity mixture: 0 (w/o); 1(w/)
KMix      1                             ! Triaxiality/K-mixture: 0 (w/o); 1(w/)
IsHFB     1                             ! 0 (HF); 1(HFB)
IntType   1                             ! 0: shell-model; 1: Chiral
InME3B    0                             ! implicit ME3B: 0 (w/o); 1 (w/)
IntID     ${IntID}                      !
ValID     ${eMax}                        ! model space: emax0N 
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
tolgrad   1.e-4                         !  precision
INP_wf    ${iwf}                        ! (0)read wf from fort.10; (1) read from saved file; (2) from a random wf
PNP_phi   1                             !  1 (HFB) or 7 (VAPHFB)
........................................ constraint quantities
kick_step 0.0
B_G       1                            ! 0 (Q value); 1 (bet,gam)
NCONS     8
Q20t      1 0.d0                       ! 3: Q20 or beta
Q22t      1 0.d0                       ! 4: Q22 or gamma
Q21t      1 0.d0                       ! 5:
Q20m      0 0.d0                       ! 6: Q20m or beta
Q22m      1 0.d0                       ! 7: Q22m or gamma
Q21m      1 0.d0                       ! 8
JX        1 0.d0                       ! 9
JY        0 0.d0                       ! 10
JZ        1 0.d0                       ! 11
P00_10    1 0.0d0                      ! 12  P_TMT_JMJ ! T=0, L=0, J=S=1, MJ=0
P1m1_00   0 0.d0                      ! 13
P1p1_00   0 0.d0                      ! 14
P10_00    0 0.d0                       ! 15
P00_1m1   0 0.0d0                      ! 16  P_TMT_JMJ ! T=0, L=0, J=S=1, MJ=-1
P00_1p1   0 0.0d0                      ! 17  P_TMT_JMJ ! T=0, L=0, J=S=1, MJ=+1
Q40t      0 0.d0                       ! 18: Q40 
Q40m      0 0.d0                       ! 19: Q40 
EOF

nq=16 #17
cat <<EOF > betgam.dat
 ............
 max
 bet     gam    P00(T=0,J=1)    hbar*omega 
 ---------------------
 13 
-0.4     0       0   0
-0.3     0       0   0
-0.2     0       0   0
-0.1     0       0   0
 0.0     0       0   0
 0.1     0       0   0
 0.2     0       0   0 
 0.3     0       0   0
 0.4     0       0   0
 0.5     0       0   0
 0.6     0       0   0
 0.7     0       0   0
 0.8     0       0   0
EOF

#.... copy matrix elements from EOM-IMSRG calculation
cd ${PATH_WORK}

export OMP_NUM_THREADS=8
mpirun -n 4 ./HFB >Ar46_MF_out
#mpirun -n 6 ./HFB  
#mpirun -n ${nq} ./HFB  
