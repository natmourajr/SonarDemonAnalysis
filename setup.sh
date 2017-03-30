#!/bin/bash 

# SonarDemonAnalysis Setup Script
#
# Author: natmourajr@gmail.com
#

# Env Variables

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    # Ubuntu
    export INPUTDATAPATH=/home/natmourajr/Public/Marinha/Data
elif [[ "$OSTYPE" == "darwin"* ]]; then
    # Mac OSX
    export INPUTDATAPATH=/Users/natmourajr/Workspace/Doutorado/Data/Sonar/Demon
    # For matplotlib
    export LC_ALL=en_US.UTF-8
    export LANG=en_US.UTF-8
fi
