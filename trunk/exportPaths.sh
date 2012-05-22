#/bin/bash

#This is the path where you've put Jelly
export JELLYPATH=/users/p-pacbio/english/Jelly/Jelly/

export PATH=$PATH:$JELLYPATH
if [ -z "${SEYMOUR_HOME}" ]
then
    #This is where you've installed SMRTAnalysis
    SEYMOUR_HOME="/hgsc_software/pacbioSMRTAnalysis/smrtanalysis"
fi

#Adds PacBio SMRTAnalysis suite into the environment
source $SEYMOUR_HOME/etc/setup.sh
