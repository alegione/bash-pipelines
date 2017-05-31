#!/bin/bash
# A few lines of code that can be handy to include at the start of any bash script.
# Mounts shared folders, sets terminal colour variables and maximises window size


#get CPU metrics
CORES=$(nproc)

# The below can be adjusted if you have a folder that requires mounting (which is not auto-mounted). Requires a text file
# 'mount_test.txt' in the target directory
if [ ! -e Desktop/Shared_folder/mount_test.txt ]; then
	sudo mount -t vboxsf Sequences Desktop/Shared_folder/
fi

#Set colour variables
RED='\033[1;31m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
PURPLE='\033[1;35m'
NOCOLOUR='\033[0m'

#get screen size and maximise
LINE=`xrandr -q | grep Screen`
WIDTH=`echo ${LINE} | awk '{ print $8 }'`
HEIGHT=`echo ${LINE} | awk '{ print $10 }' | awk -F"," '{ print $1 }'`
echo -e "\e[4;$HEIGHT;${WIDTH}t"

