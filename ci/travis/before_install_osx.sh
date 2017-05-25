#!/bin/bash

brew update
brew install qt5
brew link qt5 --force

ENVFILE=qt-5.8.0.env
echo Create $ENVFILE
cat << EOF > $ENVFILE
export PATH="/usr/local/opt/qt/bin:$PATH"
EOF
