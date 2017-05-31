#!/bin/bash

# renaming script

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]; then
	echo "usage: please enter target folder, a rename file (in that order) and extension"
	exit
fi

Folder=$(printf $1 | sed 's/\/$//')
renamefile=$2
extension=$3
find $Folder -type f | awk -F. '!a[$NF]++{print $NF}' > tmp.rename_extensions.txt


while read old new; do
	rename $old ${new} $Folder/${old}*${extension}*
done < $renamefile


rm tmp.rename_extensions.txt
