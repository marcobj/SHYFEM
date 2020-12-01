
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
 * nlist.h - NodeList utility routines
 *
 * revision log :
 *
 * 27.07.1995	ggu	routines written from scratch
 * 11.08.1995	ggu	Nodelist renamed to NodeList
 * ...		ggu	Nodelist_type * substituted by NodeList
 *
\************************************************************************/


#ifndef __GUH_NLIST_
#define __GUH_NLIST_


typedef struct nodelist_tag {
        int count;
        int *index;
} NodeList_type;

typedef NodeList_type *NodeList;

NodeList MakeNodeList( int count );
void InvertNodeList( NodeList L );
void FreeNodeList( NodeList L );


#endif

