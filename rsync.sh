#!/bin/bash

usage (){
    echo "Usage:"
    echo "$1 [pull | push]";
    exit;
}

# no "/" at the end
dest="hyhy123@euler.phys.cmu.edu:/home/hyhy123/backup/src/3bdfxdat_public"

#echo $#
if [ $# == 0 ]
then
    usage $0
fi

key=$1

if [  $key == "pull" ]
then
    rsync -a -e "ssh" $dest/ .
fi


if [  $key == "push" ]
then
    rsync -a -e "ssh" ./ $dest
fi
