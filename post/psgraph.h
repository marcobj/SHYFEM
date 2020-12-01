
/************************************************************************\
 *
 *    Copyright (C) 1992,1994,1996-1997,2004,2009  Georg Umgiesser
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
 * psgraph.h - graphic routines for postscript output
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 06.04.1994	ggu	copyright notice added to file
 * 23.02.1996	ggu	new routines added (tailored for PS)
 * 11.09.1997	ggu	minor modifications (PsGraphOpen,...)
 * 12.09.1997	ggu	PsAdjustScale becomes static, PsWindow -> PsSetWorld
 * 28.04.2004	ggu	new routines dash, rotate text and arc
 * 21.01.2009	ggu	new routine for text centering inserted
 * 12.06.2009	ggu	new routines for scale factor and no clipping
 *
\************************************************************************/


#ifndef __GUH_PSGRAPH_
#define __GUH_PSGRAPH_

/* #include "general.h" */

void PsGraphInit( char *file );
void PsGraphOpen( void );
void PsGraphOpenFile( char *file );
void PsGraphClose( void );
void PsStartPage( void );
void PsEndPage( void );
void PsNewPage( void );

void PsSetViewportNoClip( float left , float bottom , float right , float top );
void PsSetViewport( float left , float bottom , float right , float top );
void PsGetViewport( float *left , float *bottom , float *right , float *top );
void PsClearViewport( void );
void PsSetWorld( float xmin , float ymin , float xmax , float ymax );
void PsRectifyScale( void );
void PsFactorScale( float xfact , float yfact  );

void PsLine( float x1 , float y1 , float x2 , float y2 );
void PsMove( float x , float y );
void PsPlot( float x , float y );
void PsPoint( float x , float y );

void PsAreaFill( int ndim , float *x , float *y );
void PsRectFill( float x1 , float y1 , float x2 , float y2 );

void PsArc( float x0, float y0, float r, float ang1, float ang2 );
void PsArcFill( float x0, float y0, float r, float ang1, float ang2 );

void PsTextPointSize( int size );
void PsTextRealSize( float size );
void PsTextFont( char *font );
void PsTextRotate( float angle );
void PsText( float x , float y , char *s );
void PsTextDimensions( char *s , float *width , float *height );
void PsTextSetCenter( float hc , float vc );

void PsSetGray( float gray );
void PsSetRGB( float red , float green , float blue );
void PsSetHSB( float hue , float sat , float bri );
void PsSetHue( float hue );
void PsSetDefaultColorTable( int color_table );
void PsSetGenericColor( float color );
void PsPaintWhite( int ipaint );

void PsInitColorTable( int size );
void PsSetColorTable( int i , int type , float c1 , float c2 , float c3 );
void PsSetColorPen( int i );
void PsSetColorRange( float col );

void PsSetDashPattern( float fact, float offset, int n, float *array );
void PsResetDashPattern( void );
void PsSetLineWidth( float width );
void PsSetPointSize( float size );

void PsComment( char *s );
void PsRealXY( float vx , float vy , float *x , float *y );
void PsPageXY( float x , float y , float *vx , float *vy );
void PsCmLength( float *xcm, float *ycm );

void PsFlush( void );
void PsSync( void );
void PsSyncron( int sync );
int PsActScreen( void );
void PsSetTest( void );
void PsClearTest( void );
void PsGetMinMaxPix( int *xmin , int *ymin , int *xmax , int *ymax );
void PsSetMinMaxPix( int xmin , int ymin , int xmax , int ymax );


#endif /* __GUH_PSGRAPH_ */
