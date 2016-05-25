#!/usr/bin/env python

import os, sys, time
import subprocess
import urllib
import zipfile
import tarfile
from sys import platform as _platform
import platform

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

def climbUpDirectoryTree(filepath,recursion_depth):
    if filepath:
        for i in range(recursion_depth):
            filepath = os.path.join(filepath,'..')
        resultant_path = os.path.abspath(filepath)
    else:
        resultant_path = []

    return resultant_path

def checkForExistingDirectoriesAndPrompt(dir_name,file_trigger,root_dir,recursion_depth):
    filepath = find(file_trigger,root_dir)

    preinstalled_path = climbUpDirectoryTree(filepath,recursion_depth)

    skip = None
    if preinstalled_path:
        print "\n  "+dir_name+" seems to be installed at this location: " + str(preinstalled_path)
        res = raw_input("  Shall I re-install "+dir_name+", possibly overwriting this location? [y/N] ")
        if res.lower() != 'y':
            skip = True
            print "  Skipping "+dir_name+" install.  I will use the pre-installed "+dir_name+" path."
        else:
            skip = False

    return skip, preinstalled_path

def checkIfArchiveExistsAndPrompt(dir_name,filenames,root_dir):
    # TODO: OS-specific checks
    flag_for_download = True
    if any([os.path.exists(os.path.join(root_dir, a)) for a in filenames]):
        res = raw_input("\n  "+dir_name+" archive found.  Overwrite? [y/N] ")
        flag_for_download = True if res.lower() == "y" else False

    return flag_for_download

if __name__ == "__main__":
    # Make sure we are not running on a network drive on Windows
    if os.path.abspath(__file__).startswith(r"\\"):
        print "ERROR: Sorry, this script cannot be run from a network drive."
        #             ... (because windows UNC paths seem to break things)
        print "Please move the reasyns folder a local disk and run this script again."
        print "(You can then move it back later.)" 
        print
        print "Press any key to quit..."
        raw_input()
        sys.exit(2)

    # Find the root directory
    root_dir = os.path.normpath(os.path.abspath(os.path.dirname(__file__)))
    lib_dir = os.path.join(root_dir, "lib")

    os.chdir(lib_dir)

    # -----------------
    #      SeDuMi 
    # -----------------
    print "\n-> Setting up SeDuMi..."

    # There could be several sedumi installations - single out the one we need by name
    skip_sedumi, preinstalled_sedumi_path = checkForExistingDirectoriesAndPrompt('SeDuMi','install_sedumi.m',lib_dir,1)

    if skip_sedumi:
        sedumi_dir = preinstalled_sedumi_path
    else:
        if checkIfArchiveExistsAndPrompt('SeDuMi',['sedumi.zip'],lib_dir):
            print ">>> Downloading SeDuMi..."
            urllib.urlretrieve("https://github.com/sqlp/sedumi/archive/master.zip", \
                os.path.join(lib_dir, "sedumi.zip"))

        print ">>> Unzipping SeDuMi..."
        with zipfile.ZipFile(os.path.join(lib_dir, "sedumi.zip")) as sedumi_zip:
            sedumi_zip.extractall(lib_dir)

        sedumi_dir = climbUpDirectoryTree(os.path.join(find('install_sedumi.m',lib_dir)),1)

    # -----------------
    #    Ellipsoids 
    # -----------------
    print "\n-> Setting up Ellipsoids..."

    skip_ellipsoids, preinstalled_ellipsoids_path = checkForExistingDirectoriesAndPrompt('Ellipsoids','ellipsoids_init.m',root_dir,1)

    if skip_ellipsoids:
        ellipsoids_dir = preinstalled_ellipsoids_path
    else:
        if checkIfArchiveExistsAndPrompt('Ellipsoids',['ellipsoids.zip'],lib_dir):
            print ">>> Downloading Ellipsoids..."
            urllib.urlretrieve("https://github.com/SystemAnalysisDpt-CMC-MSU/ellipsoids/releases/download/v1.1.3/elltbx_1_1_3_full.zip", \
                os.path.join(lib_dir, "ellipsoids.zip"))

        print ">>> Unzipping Ellipsoids..."
        with zipfile.ZipFile(os.path.join(lib_dir, "ellipsoids.zip")) as ellipsoids_zip:
            ellipsoids_zip.extractall(lib_dir)

        ellipsoids_dir = climbUpDirectoryTree(os.path.join(find('ellipsoids_init.m',root_dir)),1)

    # -----------------
    #       MPT 
    # -----------------
    print "\n-> Setting up MPT..."

    skip_mpt, preinstalled_mpt_path = checkForExistingDirectoriesAndPrompt('MPT','mpt_init.m',root_dir,1)

    if skip_mpt:
        mpt_dir = preinstalled_mpt_path
    else:
        if checkIfArchiveExistsAndPrompt('MPT',['mpt.zip'],lib_dir):        
            print ">>> Downloading MPT..."
            urllib.urlretrieve("http://people.ee.ethz.ch/~mpt/2/downloads/mpt_ver263.zip", \
                os.path.join(lib_dir, "mpt.zip"))

        print ">>> Unzipping MPT..."
        with zipfile.ZipFile(os.path.join(lib_dir, "mpt.zip")) as mpt_zip:
            mpt_zip.extractall(lib_dir)

        mpt_dir = climbUpDirectoryTree(os.path.join(find('mpt_init.m',root_dir)),1) 

    # -----------------
    #      Drake 
    # -----------------
    print "\n-> Setting up Drake..."

    skip_drake, preinstalled_drake_path = checkForExistingDirectoriesAndPrompt('Drake','balanceQuadForm.m',root_dir,3)

    if skip_drake:
        drake_dir = preinstalled_drake_path
    else:
        if checkIfArchiveExistsAndPrompt('Drake',['drake.zip','drake.tar.gz'],lib_dir):
            print ">>> Downloading Drake..."

            if _platform in ['win32','cygwin']:
                if platform.machine().endswith('64'):
                    urllib.urlretrieve("https://github.com/RobotLocomotion/drake/releases/download/v0.9.11/drake-v0.9.11-win64.zip", \
					os.path.join(lib_dir, "drake.zip"))
                else:
                    urllib.urlretrieve("https://github.com/RobotLocomotion/drake/releases/download/v0.9.11/drake-v0.9.11-129-g20a1c65-win32.zip", \
                        os.path.join(lib_dir, "drake.zip"))
            elif _platform in ['linux','linux2']:
                urllib.urlretrieve("https://github.com/RobotLocomotion/drake/releases/download/v0.9.11/drake-v0.9.11-linux.tar.gz", \
                    os.path.join(lib_dir, "drake.tar.gz"))
            elif _platform == "darwin":
                urllib.urlretrieve("https://github.com/RobotLocomotion/drake/releases/download/v0.9.11/drake-v0.9.11-mac.tar.gz", \
                    os.path.join(lib_dir, "drake.tar.gz"))
            else:
                raise Exception("Unrecognized OS!")

        print ">>> Unzipping Drake..."

        if os.path.exists(os.path.join(lib_dir,"drake.tar.gz")):
            with tarfile.open(os.path.join(lib_dir, "drake.tar.gz"), "r:gz") as drake_tar:
                drake_tar.extractall(lib_dir)
        else:        
            with zipfile.ZipFile(os.path.join(lib_dir, "drake.zip")) as drake_zip:
                drake_zip.extractall(lib_dir)

        drake_dir = climbUpDirectoryTree(os.path.join(find('balanceQuadForm.m',root_dir)),3)

    # -----------------
    #      Mosek 
    # -----------------
    print "\n-> Setting up Mosek (optional)...\n"
    skip_mosek = None
    while skip_mosek is None:
        skip_mosek_response = raw_input("  Skip this step [y/N]? ")
        if skip_mosek_response.lower() == "y":
            print "\n>>> Skipping."
            skip_mosek = True
        elif skip_mosek_response.lower() == "n" or skip_mosek_response == "":
            skip_mosek = False

    if not skip_mosek:

        skip_mosek, preinstalled_mosek_path = checkForExistingDirectoriesAndPrompt('Mosek','mosekopt.py',root_dir,4)

        if skip_mosek:
            mosek_dir = preinstalled_mosek_path
        else:
            if checkIfArchiveExistsAndPrompt('Mosek',['mosek.msi','mosek.tar.bz2'],lib_dir):
                print ">>> Downloading Mosek..."

                if _platform in ['win32','cygwin']:
                    if platform.machine().endswith('64'):
                        urllib.urlretrieve("http://download.mosek.com/stable/7.1.0.53/moseksetupwin64x86.msi", os.path.join(lib_dir, "mosek.msi"))
                    else:
                        urllib.urlretrieve("http://download.mosek.com/stable/7.1.0.53/moseksetupwin32x86.msi", os.path.join(lib_dir, "mosek.msi"))
                elif _platform in ['linux','linux2']:
                    if platform.machine().endswith('64'):
                        urllib.urlretrieve("http://download.mosek.com/stable/7.1.0.53/mosektoolslinux64x86.tar.bz2", os.path.join(lib_dir, "mosek.tar.bz2"))
                    else:
                        urllib.urlretrieve("http://download.mosek.com/stable/7.1.0.53/mosektoolslinux32x86.tar.bz2", os.path.join(lib_dir, "mosek.tar.bz2"))                
                elif _platform == "darwin":
                    urllib.urlretrieve("http://download.mosek.com/stable/7.1.0.53/mosektoolsosx64x86.tar.bz2", os.path.join(lib_dir, "mosek.tar.bz2"))
                else:
                    raise Exception("Unrecognized OS!")

            print ">>> Unpacking Mosek..."

            if os.path.exists(os.path.join(lib_dir,"mosek.tar.bz2")):
                with tarfile.open(os.path.join(lib_dir, "mosek.tar.bz2"), "r:bz2") as mosek_tar:
                    mosek_tar.extractall(lib_dir)
            else:
				os.system('msiexec /i %s /qn' % os.path.join(lib_dir, "mosek.msi"))
				#subprocess.call('msiexec /i %s /qn' % os.path.join(lib_dir, "mosek.msi"), shell=True)

            mosek_dir = climbUpDirectoryTree(os.path.join(find('mosekopt.m',root_dir)),4)

    # -----------------
    #      LTLMoP 
    # -----------------
    print "\n-> Setting up reasyns-compatible version of LTLMoP (optional)...\n"
    skip_ltlmop = None
    while skip_ltlmop is None:
        skip_ltlmop_response = raw_input("  Skip this step [y/N]? ")
        if skip_ltlmop_response.lower() == "y":
            print "\n>>> Skipping."
            skip_ltlmop = True
        elif skip_ltlmop_response.lower() == "n" or skip_ltlmop_response == "":
            skip_ltlmop = False

    if not skip_ltlmop:

        skip_ltlmop, preinstalled_ltlmop_path = checkForExistingDirectoriesAndPrompt('LTLMoP','specEditor.py',root_dir,1)
        
        if skip_ltlmop:
            ltlmop_dir = preinstalled_ltlmop_path
        else:
            if checkIfArchiveExistsAndPrompt('LTLMoP',['ltlmop.zip'],lib_dir):
                print ">>> Downloading LTLMoP..."
                urllib.urlretrieve("https://github.com/jadecastro/LTLMoP/archive/reasyns_fast.zip", os.path.join(lib_dir, "ltlmop.zip"))

            print ">>> Unzipping LTLMoP..."
            with zipfile.ZipFile(os.path.join(lib_dir, "ltlmop.zip")) as ltlmop_zip:
                ltlmop_zip.extractall(lib_dir)

            ltlmop_dir = climbUpDirectoryTree(os.path.join(find('specEditor.py',root_dir)),1)
    
    # --------------------
    #    Write to M-file 
    # --------------------  
    print "\n-> Writing Matlab path initialization file...\n"

    with open(os.path.join(root_dir,'reasyns_init.m'),'w') as matlab_path_file:
        matlab_path_file.write("%% This is an automatically-generated file.  See setup.py in the top-level directory. \n\n")
        matlab_path_file.write("restoredefaultpath\n\n")
        matlab_path_file.write("warning off\n\n")
        matlab_path_file.write("addpath(genpath('"+str(os.path.abspath(root_dir))+"'));\n")
        matlab_path_file.write("rmpath(genpath('"+str(os.path.abspath(os.path.join(drake_dir,'drake','examples')))+"'));\n")
        matlab_path_file.write("rmpath(genpath('"+str(os.path.abspath(os.path.join(ellipsoids_dir,'solvers','SeDuMi_1_1')))+"'));\n")
        matlab_path_file.write("rmpath(genpath('"+str(os.path.abspath(os.path.join(mpt_dir,'solvers','SeDuMi_1_3')))+"'));\n\n")
        matlab_path_file.write("run('"+str(os.path.abspath(os.path.join(drake_dir,'drake','addpath_drake')))+"');\n")
        matlab_path_file.write("run('"+str(os.path.abspath(os.path.join(drake_dir,'externals','spotless','spotless','spot_install')))+"');\n\n")
        matlab_path_file.write("cd('"+str(os.path.abspath(sedumi_dir))+"');\n")
        matlab_path_file.write("install_sedumi;\n")
        matlab_path_file.write("cd('"+str(os.path.abspath(root_dir))+"');\n\n")
        matlab_path_file.write("examplesPath = '"+str(os.path.abspath(os.path.join(root_dir,'examples')))+"';\n\n")
        matlab_path_file.write("warning on\n")

    print "\n>>> Done!\n"

    if os.name == "nt":
        raw_input("Press any key to quit...")
