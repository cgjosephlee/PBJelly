#/bin/bash

#This is the path where you've put Jelly
export JELLYPATH=/users/english/english/Jelly/Jelly_12.5.10/
export PATH=$PATH:$JELLYPATH
if [ -z "${SEYMOUR_HOME}" ]
then
    #This is where you've installed SMRTAnalysis
    SEYMOUR_HOME="/opt/smrtanalysis/"
fi

#Adds PacBio SMRTAnalysis suite into the environment
source $SEYMOUR_HOME/etc/setup.sh
