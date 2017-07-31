#!/bin/bash
# A few lines of code that can be handy to include at the start of any bash script.
# Sets CPU and RAM variables, terminal colour variables and maximises window size


#!/bin/bash

#get CPU and RAM metrics
CORES=$(nproc)
RAM=$(free -g | tr -s "[:space:]" "\t" | cut -f11)

#Only useful if you have a shared drive you always want to check is mounted,
#requires a text file 'mount_test.txt' in the shared folder

#if [ ! -e Desktop/Shared_folder/mount_test.txt ]; then
#	sudo mount -t vboxsf Sequences Desktop/Shared_folder/
#fi

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
