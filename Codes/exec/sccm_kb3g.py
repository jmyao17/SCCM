#"' This python file is made for compiling codes
#   and run them
#'"
# import modules
import sys, os, shutil
from shutil import copyfile

# path
os.system("source set_path.sh")
dir_dst="./"

# Options
OPT1 = [True,True,True]
OPT2 = [True,True,False]
OPT3 = [True,False,False]
OPT4 = [False,False,False]

#compile_hfb   = True #False
compile_hfb,compile_pnamp,compile_gcm =OPT3 

# parameters
run_hfb,run_pnamp,run_gcm = OPT3

# copy source *.f90 files for compiling
if compile_hfb:
   dir_src="../HFB/"
   for filename in os.listdir(dir_src):
       if filename.endswith('.f90'):
           shutil.copy(dir_src + filename, dir_dst)
   copyfile("makefile4HFB","makefile")
   os.system("make")

if compile_pnamp:
    dir_src="../PNAMP/"
    for filename in os.listdir(dir_src):
        if filename.endswith('.f90'):
            shutil.copy(dir_src + filename, dir_dst)
    copyfile("makefile4PNAMP","makefile")
    os.system("make")


if compile_gcm:
    dir_src="../GCM/"
    for filename in os.listdir(dir_src):
        if filename.endswith('.f90'):
            shutil.copy(dir_src + filename, dir_dst)
    copyfile("makefile4GCM","makefile")
    os.system("make")

# clean all the source files
for filename in os.listdir(dir_dst):
    if filename.endswith('.o') or filename.endswith('.f90'):
        os.remove(dir_dst + filename)


# run the code
if run_hfb:
    os.system("./Ti48_KB3G.sh")
if run_pnamp:
    os.system("./PNAMP.sh")
if run_gcm:
    os.system("./GCM.sh")

#os.remove("fort.*")

