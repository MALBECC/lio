#!/bin/bash

# Gets this directory and sets it as LIOHOME.
INSTALL_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
export LIOHOME=$INSTALL_DIR

# Adds LIO libraries.
LIOLIBS=$LIOHOME/g2g:$LIOHOME/lioamber
export LIBRARY_PATH=$LIBRARY_PATH:$LIOLIBS
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIOLIBS

# Adds LIO tools and executables to path.
LIOEXES=$LIOHOME/liosolo:$LIOHOME/tools
export PATH=$PATH:$LIOEXES

