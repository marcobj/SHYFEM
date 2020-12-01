
/************************************************************************\
 *
 *    Copyright (C) 1995  Georg Umgiesser
 *
 *    This file is part of SHYFEM.
 *
 *    SHYFEM is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    SHYFEM is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with SHYFEM. Please see the file COPYING in the main directory.
 *    If not, see <http://www.gnu.org/licenses/>.
 *
 *    Contributions to this file can be found below in the revision log.
 *
\************************************************************************/


/************************************************************************\
 *
 * meshck.h - check routines for mesh
 *
 * revision log :
 *
 * 27.07.1995	ggu	routines written from scratch
 *
\************************************************************************/


#ifndef __GUH_MESHCK_
#define __GUH_MESHCK_


#include "nlist.h"
#include "heap.h"

int CheckInput( void );
void CheckBoundary( void );
void CheckConvex( NodeList hull, NodeList intern );
void CheckHeapPropertyUP( HeapTable HL, int up );
void CheckHeap( HeapTable HL );
void CheckCircumCircleProperty( void );
void CheckArea( void );
void CheckStatus( int id );
void CheckNeibor( int id );
void CheckCircumCircle( void );


#endif


