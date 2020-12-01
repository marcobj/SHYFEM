#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# runs shyfem model in mpi mode
#
#-----------------------------------------

fem3d=~/shyfem/fem3d
shyfem=$fem3d/shyfem

if [ $# -lt 2 ]; then
  echo "Usage: mpi_run.sh nproc [shyfem-options] str-file"
  exit 1
fi
nproc=$1
shift
str=$*

#-----------------------------------------

if [ $nproc = 0 ]; then
  $shyfem $str
  status=$?
elif [ $nproc = 1 ]; then
  $shyfem -mpi $str
  status=$?
else
  /usr/bin/mpirun -np $nproc $shyfem -mpi $str
  status=$?
fi

#-----------------------------------------

if [ $status -eq 99 ]; then
  exit 0
else
  echo "exit status is $status ... error"
  exit 3
fi

#-----------------------------------------

echo "command run with nproc=$nproc and status=$status"

#-----------------------------------------

