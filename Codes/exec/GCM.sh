#!/bin/bash
#PBS -l walltime=23:50:00
#PBS -l nodes=5:ppn=1      ###nodes:ppn - how many nodes & cores per node (ppn)
#PBS -o Ge76 
#PBS -N Ge76
### mem: amount of memory that the job will need
#PBS -l mem=5gb


IntID=GCN2850
# .... information about the initial nucleus
Nucl=Ge76
ZZ=4
NN=16
Nsh=6
hw=09
#ValID=eMax0${Nsh}
ValID=pf5g9
hwHO=hwHO${hw}

INT=0   # 0 or 4
HFB=1   # 2(VAPHFB); 1(HFB); 0(HF)
Rho3B=0  # 1 (calc.); 0 (skip)
E2=1     # 1 (calc.); 0 (skip)
idens=0

Nphi=5 #7
Nalp=1 #14 #6 #14 #6 #6 #1 #12 #6 #6
Nbet=16  #14 #8 #14 #8 #16 #8 #16 #8 #16 #8  #16
Ngam=1 #14 #12 #14 #12 #12 #1 #12  #12
vs4me3b=pf5g9 #sdf7
filename=${Nucl}_${eMax}_${hwHO}.dat

#file4me3b=me3b_p.val #me3b_sd.val

#cp ../Int/${file4me3b} ../Int/me3b.val


filename=${Nucl}_${eMax}_${hwHO}.dat

cat <<EOF > input.dat 
IntType   0                             ! 0: shell-model; 1: Chiral
IsHFB     ${HFB}                        ! 0 (HF); 1(HFB)
IntID     ${IntID}                     !magic18-20_srg0953            ! interaction
ValID     ${ValID}                        ! model space (sd), (pf), (pf5g9)
hwHO      ${hw}                            ! frequency 
FlowPars  s000                          ! s value
nprot     ${ZZ}                            ! proton number of initial nucleus
nneut     ${NN}                            ! neutron number of intial nucleus
PNP       ${Nphi}                       ! PNP
Nalp      ${Nalp}                       ! alp for AMP
Nbet      ${Nbet}                       ! beta for AMP
Ngam      ${Ngam}                       ! gam for AMP
Jmax      6                             ! max J from AMP
kmax      2                             ! kmax for a given J
Idens     ${idens}                      ! 1(calc. GCM density); 0 (read from file)
Rho3B     ${Rho3B}                      ! 1 (with Rho3B ); 0 (without Rho3B)
vs4me3b   ${vs4me3b}                                ! model space for Rho3B (sd), (pf), (pf5g9)
cutoff    1.e-5                         ! cut off in norm
J0k       1                               ! print out the w.f. of the k-th 0+ state 
NOSJ0     0                               ! 0 (based on variance of E); otherwise choose the value 
NOSJ2     0                               ! 
NOSJ4     0                               ! 
crank     0                               ! 1 (w/ cranking); 0 (w/o cranking)
EOF

./GCM 
