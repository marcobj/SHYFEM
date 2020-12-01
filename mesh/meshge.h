
/************************************************************************\
 *
 *    Copyright (C) 1995,1997,2011  Georg Umgiesser
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
 * meshge.h - geometric routines for mesh
 *
 * revision log :
 *
 * 01.08.1995	ggu	routines written from scratch
 * 15.10.1997	ggu	new routine InClosedLine()
 * 16.10.1997	ggu	new routine FindElement()
 * 17.10.1997	ggu	new routine TurnClosedLine()
 * 05.12.2011	ggu	new general routines
 *
\************************************************************************/


#ifndef __GUH_MESHGE_
#define __GUH_MESHGE_


#include "mesh.h"
#include "hash.h"
#include "list.h"
#include "nlist.h"

/*
#include "nlist.h"
#include "heap.h"
*/


void MakeCircumCircle( Elem_type *pe );
void ControlCircumCircle( Elem_type *pe );
float MakeInCircleRadius( Elem_type *pe );
void MakeCM( Elem_type *pe, float *xm, float *ym );

int InCircumCircle( Elem_type *pe , float x , float y );
int InElemCircumCircle( Elem_type *circle , Elem_type *pe , float fact );
int InConvex( int n , float *xe , float *ye , float x , float y );
int InElement( Hashtable_type H , Elem_type *pe , float x , float y );
void PolyMinMax( Hashtable_type H , Elem_type *pe , Rect *r );
void PolyMinMaxIndex( Hashtable_type H , int nvert , int *index , Rect *r );

Elem_type *FindElement( Hashtable_type H, float x, float y );
Elem_type *FindXYElement( ListTable L, float x, float y );

float AreaElement( Hashtable_type H , Elem_type *pe );
float AreaLine( Hashtable_type H , Line_type *pl );
float AreaConvex( NodeList list );
float LengthConvex( NodeList list );

double angle( float x1, float y1, float x2, float y2, float x3, float y3 );

int InClosedLine( int ndim , float *xl , float *yl , float x , float y );
int TurnClosedLine( int ndim , float *xl , float *yl );

int IsPointInLine( Line_type *pl , float x , float y );
int IsLineInLine( Line_type *plext , Line_type *plint );
int IsLineClosed( Line_type *pl );


#endif


