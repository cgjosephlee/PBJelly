#/bin/bash

#This is the path where you've put the PBSweet
export SWEETPATH=/users/p-pacbio/english/Jelly/DevJelly/restructure/

#for libraries
export PYTHONPATH=$PYTHONPATH:$SWEETPATH

export PATH=$PATH:$SWEETPATH/bin/

if [ -z "${SEYMOUR_HOME}" ]
then
    #This is where you've installed SMRTAnalysis
    SEYMOUR_HOME="/hgsc_software/pacbioSMRTAnalysis/smrtanalysis"
fi

#Adds PacBio SMRTAnalysis suite into the environment
source $SEYMOUR_HOME/etc/setup.sh
