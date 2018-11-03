#!/bin/bash
#PBS -l nodes=1:ppn=4,walltime=43:50:00,mem=80gb      ###nodes:ppn - how many nodes & cores per node (ppn)
#PBS -o Ge76 
#PBS -N PNAMP-Ge76

IntID=GCN2850
Nucl=Ge76
ZZ=4
NN=16
Nsh=6
ValID=pf5g9
hw=10 #09 #20
eMax=eMax0${Nsh}
hwHO=hwHO${hw}

INT=0   # 0 or 4
HFB=1 #1   # 2(VAPHFB); 1(HFB); 0(HF)
Rho2B=0  # 1 (calc.); 0 (skip); 
Rho3B=0 #1  # 1 (calc.); 0 (skip)
E2=1     # 1 (calc.); 0 (skip)

Hs=H0 #H0
Flow=s000 #s999

Nphi=5 #5
Nalp=1  #12 #6   #6 #14 #6 #14 #6 #4 #14 #6 #12 #12 #6 #6 #8    # [0,pi]
Nbet=16  #12 #8   #8 #14 #8 #14 #14 #14 #8 #16 #8 #16 #8 #16 #8 #16    # [0,pi/2]
Ngam=1   #12 #14 #12 #14 #14 #6 #4 #14 #12 #12  #16   # [0,2pi] 

vs4me3b=pf5g9


PATH_WORK=./ #/mnt/home/yaojiang/GCM/ABGCM/${HFBdir}
PATH_INT=../${IntID}
me1b=SM_GCN2850_ME1B.dat #EOM_Ca48_chi2b3b_srg0625_eMax04_hwHO020_unnorm_ham_1b.dat
me2b=SM_GCN2850_ME2B.dat #EOM_Ca48_chi2b3b_srg0625_eMax04_hwHO020_unnorm_ham_2b.dat

cat <<EOF > chi2b3b.int 
File4ME1B:    ${me1b}
File4ME2B:    ${me2b} 
EOF

#file4me3b=me3b_psd.val #sdf.val #p1sdf7.val #sd.val #sd.val  #me3b_sd.val

#cp ../Int/${file4me3b} ../Int/me3b.val

cd ${PATH_WORK}
cat <<EOF> input.dat 
NP-Mix    1                             ! np mixture: 0 (w/o); 1(w/)
IP-Mix    1                             ! parity mixture: 0 (w/o); 1(w/)
KMix      1                             ! Triaxiality/K-mixture: 0 (w/o); 1(w/)
IsHFB     ${HFB}                             ! 0 (HF); 1(HFB)
IntType   0                             ! 0: shell-model; 1: Chiral
IntID     ${IntID}                      !magic18-20_srg0953            ! interaction
ValID     ${ValID}                      ! model space (sd), (pf), (pf5g9)
hw        ${hw}                         ! i2,frequency of HO 
COM       2                             ! Center-of-Mass Correction: 1 (1B); 2 (1B+2B)
IntIMSRG  1                             ! IMSRG interation (1) or not (0)
FlowPars  ${Flow}                         ! s value
Int-JT    ${INT}                        ! 0(m-scheme); 1(JTTz-scheme); 2 (Generate V from J-scheme to JT-scheme)
IsoSpin   0                             ! (0) Int. does not have isospin symmetry; (1) Int. has isospin symm.
nprot     ${ZZ}                             ! proton number of initial nucleus
nneut     ${NN}                            ! neutron number of intial nucleus
PNP       ${Nphi}                            ! PNP
Nalp      ${Nalp}                             ! alp for AMP
Nbet      ${Nbet}                            ! beta for AMP
Ngam      ${Ngam}                             ! gam for AMP
ISCALE    0                             ! scale factor for the matrix, the pfaffian of which needs to be calc.ed
Icons     1                             ! 0 (upper triangle); 1 (all)
nq0i      1 
nq0f      5 
nq1i      1                                < 1| O | 0>
nq1f      5  
Idens     ${Rho2B}                             ! 2 (skip 2B density); 1 (calc. density); 0 (read from file)
Rho3B     ${Rho3B}                             ! 1 (with Rho3B ); 0 (without Rho3B)
vs4me3b   ${vs4me3b}                                 ! model space for Rho3B (sd), (pf), (pf5g9)
iE2       ${E2}                             ! 1 (with E2 ); 0 (without E2)
crank     0                                 ! 2 (no any symmetry); 1 (cranking code); 0 (w/o cranking code) 
EOF

cd ${PATH_WORK}

mpirun -n 1  ./PNAMP #> Ge76_PNAMP.out 
