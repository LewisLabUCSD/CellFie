#!/bin/sh

# launch MATLAB
if [ "$ARCH" == "Linux" ]; then

    export CURRENTDIR=`pwd`

    # update the COBRA Toolbox
    cd $ARTENOLIS_DATA_PATH/repos/cobratoolbox
    git pull origin master

    # change to the current directory
    echo $CURRENTDIR
    cd $CURRENTDIR

    # launch the test suite
    $ARTENOLIS_SOFT_PATH/MATLAB/$MATLAB_VER/bin/./matlab -nodesktop -nosplash < test/testAll.m

fi

CODE=$?
exit $CODE
