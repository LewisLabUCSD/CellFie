#!/bin/sh

# launch MATLAB
if [ "$ARCH" == "Linux" ]; then
    $ARTENOLIS_SOFT_PATH/$MATLAB_VER/bin/./matlab -nodesktop -nosplash -r "fprintf('Hello ARTENOLIS.\n'); quit();"
fi

CODE=$?
exit $CODE
