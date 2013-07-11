#/bin/bash

#This is the path where you've put the PBSweet

export SWEETPATH=/stornext/snfs5/next-gen/scratch/english/Jelly/DevJelly/trunk

#for libsvm 
export LIBSVM=$SWEETPATH/lib/libsvm-3.17/
export PATH=$PATH:$LIBSVM/bin/
export PYTHONPATH=$PYTHONPATH:$LIBSVM/python/

#for python modules 
export PYTHONPATH=$PYTHONPATH:$SWEETPATH
#for executables 
export PATH=$PATH:$SWEETPATH/bin/

if [ -z "${SEYMOUR_HOME}" ]
then
    #This is where you've installed SMRTAnalysis
    SEYMOUR_HOME="/hgsc_software/pacbioSMRTAnalysis/smrtanalysis"
fi

#Adds PacBio SMRTAnalysis suite into the environment
source $SEYMOUR_HOME/etc/setup.sh
