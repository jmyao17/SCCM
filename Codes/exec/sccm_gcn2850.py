#"' This python file is made for compiling codes
#   and run them
#'"
# import modules
import sys, os, shutil
from shutil import copyfile

# path
os.system("./set_path")

# parameters

Compiler=[False,False,False]

#Compiler=[True,False,False]

#compile_hfb   = True #False
compile_hfb,compile_pnamp,compile_gcm = Compiler

# set path
dir_dst="./"

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
os.system("./HFB.sh")
#os.system("./PNAMP.sh")
#os.system("./GCM.sh")

os.remove("fort.*")

