#!/usr/bin/perl

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

while(<>) {

  #if( /nlv\s*=/ ) {
  if( /hlv\s*\(\s*\S+\s*\)\s*=/ ) {
    print "$ARGV: $_";
  }
}

