#/bin/bash

#source /hgsc_software/PBSuite/pbsuiteVirtualEnv/bin/activate

#This is the path where you've install the suite.
export SWEETPATH=/users/english/english/PBSuite/trunk
#for python modules 
export PYTHONPATH=$PYTHONPATH:$SWEETPATH
#for executables 
export PATH=$PATH:$SWEETPATH/bin/
