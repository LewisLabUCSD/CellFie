#!/bin/sh

# launch MATLAB
if [ "$ARCH" == "Linux" ]; then

    export CURRENTDIR=`pwd`

    # update the COBRA Toolbox
    cd $ARTENOLIS_DATA_PATH/scratch/cobratoolbox
    git pull origin master

    cd CURRENTDIR

    # launch the test suite
    $ARTENOLIS_SOFT_PATH/MATLAB/$MATLAB_VER/bin/./matlab -nodesktop -nosplash < test/testAll.m

fi

CODE=$?
exit $CODE
