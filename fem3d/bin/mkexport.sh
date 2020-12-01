#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

################################################## make directory

echo "....................................................... making directory"

if [ ! -d export ]; then
  mkdir export
else
  rm -f export/*
fi

################################################## .f

echo "........................................................ copying f files"

cp *.f export

################################################## .h

echo "........................................................ copying h files"

cp *.h export

################################################## .F

echo "........................................................ copying F files"

cp *.F export

################################################## comment

echo "........................................................... uncommenting"

cd export
../bin/uncomment.pl *.[fFh]
rm -f *.bak
cd ..

################################################## end

echo "................................................................... done"

