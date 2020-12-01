
/************************************************************************\
 *
 *    Copyright (C) 1992,1994,1997,2001,2003-2004,2009-2010  Georg Umgiesser
 *    Copyright (C) 2015-2016  Georg Umgiesser
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
 * psgraphf.c - POST interface for plotting Postscript from Fortran
 *
 * revision log :
 *
 * 01.01.1992	ggu	routines written from scratch
 * 11.02.1994	ggu	copyright notice added to all files
 * 21.03.1994	ggu	include xgraph.h
 * 02.05.1997	ggu	changed call to qrfill (no color anymore)
 * 02.05.1997	ggu	deleted qnewp,qagray,qahue
 * 01.12.1997	ggu	new routines qtsize (NEW)
 * 26.04.2001	ggu	new routine qpsize to set size of single point
 * 18.08.2003	ggu	new routines dealing with color table
 * 28.04.2004	ggu	new routines dash, rotate text and arc
 * 20.01.2009	ggu	routine for centering text
 * 21.01.2009	ggu	better dealing with fortran strings (all in c routine)
 * 27.01.2009	ggu	bug in convert_f2c()
 * 12.06.2009	ggu	new routines for scale factor and no clipping
 * 14.09.2009	ggu	new routine qcm to get length in cm
 * 23.02.2010	ggu	new routines qcolor and qtdef
 * 11.10.2015	ggu	new routines qopenfile to open with given file name
 * 14.09.2016	ggu	code to ignore underscore
 *
\************************************************************************/


/* 
 *	contents :
 *
 *	qopen()			opens output file / window
 *	qopenfile(file)		opens output file file
 *	qclose()		closes output file / window
 *	qstart()		opens a new page (clears screen)
 *	qend()			closes a page (flushes the buffer)
 *
 *	qgetvp(xmin,ymin,xmax,ymax)	gets viewport dimensions (in cm)
 *	qsetvp(xmin,ymin,xmax,ymax)	defines window for clipping (in cm)
 *	qsetvpnc(xmin,ymin,xmax,ymax)	defines window for clipping (no clip)
 *	qclrvp()			clears viewport
 *	qworld(xmin,ymin,xmax,ymax)	sets world coordinates
 *	qrcfy()				rectifies (adjusts) scale factor
 *	qfact(xfact,yfact)		multiplies scale with factor
 *	qcm(xcm,ycm)			computes 1 cm in world coordinates
 *	
 *	qline(x1,y1,x2,y2)	draws line from (1) to (2)
 *	qmove(x,y)		moves to (x,y)
 *	qplot(x,y)		draws line from actual position to (x,y)
 *	qpoint(x,y)		draws point at (x,y)
 *
 *	qafill(n,x,y)		fills x,y (n points) with actual color
 *	qrfill(x1,y1,x2,y2)	fills rectangle with actual color
 *
 *	qarc(x0,y0,r,ang1,ang2)		draws arc with radius r around x0,y0
 *	qarcf(x0,y0,r,ang1,ang2)	fills arc with radius r around x0,y0
 *
 *	qtxts(ip)		sets text size in points
 *	qtxtr(angle)		rotate text strings with angle
 *	qtxtcc(hc,vc)		sets horizontal/vertical text centering
 *	qtxtcr(hc,vc)		sets horizontal/vertical text centering (real)
 *	qfont(font)		sets font of text
 *	qtext(x,y,s)		writes text s at (x,y)
 *	qtsize(s,w,h)		computes size of text string
 *	qignu(flag)		convert underscore to blank in text (0/1)
 *
 *	qgray(gray)		sets shade of gray (c=[0,1], 0=black, 1=white)
 *	qrgb(red,green,blue)	sets red,green,blue [0,1]
 *	qhsb(hue,sat,bri)	sets hue,sat,bri [0,1]
 *	qhue(hue)		sets hue [0,1] with sat = bri = 1
 *	qcolor(color)		sets generic color
 *	qwhite(bpaint)		paint with white
 *	qtdef(ictab)		sets default color table
 *
 *	qinitct(isize)			make color table of size isize
 *	qsetct(i,itype,c1,c2,c3)	define color i in color table
 *	qsetcpen(i)			use color i [1-isize]
 *	qsetcrange(col)			use color col [0-1]
 *
 * 	qdash(fact,offset,n,array)	sets dash pattern
 * 	qdash0()			resets dash pattern
 *	qlwidth(width)			sets line width
 *	qpsize(size)			sets point size
 *
 *	qcomm(s)		writes comment to PS file
 *
 */

/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "general.h"
#include "psgraph.h"

#define BYTES_IN_INTEGER 4	/* size of fortran integer */

typedef int fint;		/* this must correspond to above */

static int ignore_underscore = 0;

/*****************************************************************/

static void testint( void )

{
	if( BYTES_IN_INTEGER != sizeof(fint) ) {
		Error("Erroneous size of integer");
	}
}

static char* convert_f2c( char *s , fint slen )

{
	static char string[257];
	static char blank = ' ';
	static char tab = '\t';
	static char underscore = '_';
	int i;

	if( slen > 256 ) slen = 256;
	strncpy(string,s,slen);
	string[slen] = '\0';

	for( i=slen-1; i>=0; i-- ) {
	  if( string[i] != blank && string[i] != tab ) break;
	}
	string[i+1] = '\0';

	if( ignore_underscore > 0 ) {
	  for( i=slen-1; i>=0; i-- ) {
	    if( string[i] == underscore ) string[i] = blank;
	  }
	}

	/* fprintf(stderr,"converting string: |%s|\n",string); */

	return string;
}

/*****************************************************************/

void qopen_( void )

{
	testint();
	PsGraphInit("plot.ps");
}

void qopenfile_( char *s , fint slen )

{
	testint();
	PsGraphInit( convert_f2c(s,slen) );
}

void qclose_( void )

{
	PsGraphClose();
}

void qstart_( void )

{
	PsStartPage();
}

void qend_( void )

{
	PsEndPage();
}

/*****************************************************************/

void qgetvp_( float *xmin , float *ymin , float *xmax , float *ymax )

{
	PsGetViewport(xmin,ymin,xmax,ymax);
}

void qsetvp_( float *xmin , float *ymin , float *xmax , float *ymax )

{
	PsSetViewport(*xmin,*ymin,*xmax,*ymax);
}

void qsetvpnc_( float *xmin , float *ymin , float *xmax , float *ymax )

{
	PsSetViewportNoClip(*xmin,*ymin,*xmax,*ymax);
}

void qclrvp_( void )

{
	PsClearViewport();
}

void qworld_( float *xmin , float *ymin , float *xmax , float *ymax )

{
	PsSetWorld( *xmin , *ymin , *xmax , *ymax );
}

void qrcfy_( void )

{
	PsRectifyScale();
}

void qfact_( float *xfact , float *yfact )

{
	PsFactorScale( *xfact , *yfact );
}

void qcm_( float *xcm , float *ycm )

{
	PsCmLength( xcm , ycm );
}

/*****************************************************************/

void qline_( float *x1 , float *y1 , float *x2 , float *y2 )

{
	PsLine( *x1 , *y1 , *x2 , *y2 );
}

void qmove_( float *x , float *y )

{
	PsMove( *x , *y );
}

void qplot_( float *x , float *y )

{
	PsPlot( *x , *y );
}

void qpoint_( float *x , float *y )

{
	PsPoint( *x , *y );
}

/*****************************************************************/

void qafill_( fint *n , float *x , float *y )

{
	PsAreaFill( (int) *n , x , y );
}

void qrfill_( float *x1 , float *y1 , float *x2 , float *y2 )

{
	PsRectFill( *x1 , *y1 , *x2 , *y2 );
}

/*****************************************************************/

void qarc_( float *x0 , float *y0 , float *r , float *ang1, float *ang2 )

{
        PsArc( *x0 , *y0 , *r , *ang1, *ang2 );
}

void qarcf_( float *x0 , float *y0 , float *r , float *ang1, float *ang2 )

{
        PsArcFill( *x0 , *y0 , *r , *ang1, *ang2 );
}

/*****************************************************************/

void qtxts_( fint *ip )

{
	PsTextPointSize( (int) *ip );
}

void qtxtr_( float *angle )

{
        PsTextRotate( *angle );
}

void qtxtcc_( fint *hc , fint *vc )

{
	float hcr,vcr;

	hcr = *hc;
	vcr = *vc;

	PsTextSetCenter( hcr , vcr );
}

void qtxtcr_( float *hc , float *vc )

{
	PsTextSetCenter( *hc , *vc );
}

void qfont_( char *font , fint slen )

{
	PsTextFont( convert_f2c(font,slen) );
}

void qtext_( float *x , float *y , char *s , fint slen )

{
	PsText( *x , *y , convert_f2c(s,slen) );
}

void qtsize_( char *s , float *w , float *h , fint slen )

{
	PsTextDimensions( convert_f2c(s,slen) , w , h );
}

void qignu_( fint *flag )

{
	ignore_underscore = *flag;
}

/*****************************************************************/

void qgray_( float *gray )

{
	PsSetGray( *gray );
}

void qrgb_( float *red , float *green , float *blue )

{
	PsSetRGB( *red , *green , *blue );
}

void qhsb_( float *hue , float *sat , float *bri )

{
	PsSetHSB( *hue , *sat , *bri );
}

void qhue_( float *hue )

{
	PsSetHue( *hue );
}

void qcolor_( float *color )

{
	PsSetGenericColor( *color );
}

void qwhite_( fint *bpaint )

{
	PsPaintWhite( *bpaint );
}

void qtdef_( fint *ictab )

{
	PsSetDefaultColorTable( *ictab );
}

/*****************************************************************/

void qinitct_( fint *size )

{
	PsInitColorTable( *size );
}

void qsetct_( fint *i , fint *type , float *c1 , float *c2 , float *c3 )

{
	PsSetColorTable( *i-1 , *type , *c1 , *c2 , *c3 );
}

void qsetcpen_( fint *i )

{
	PsSetColorPen( *i-1 );
}

void qsetcrange_( float *col )

{
	PsSetColorRange( *col );
}

/*****************************************************************/

void qdash_( float *fact, float *offset, int *n, float *array )

{
        PsSetDashPattern( *fact, *offset, *n, array );
}

void qdash0_( void )

{
        PsResetDashPattern();
}

void qlwidth_( float *width )

{
	PsSetLineWidth( *width );
}

void qpsize_( float *size )

{
	PsSetPointSize( *size );
}

/*****************************************************************/

void qcomm_( char *s , fint slen )

{
	PsComment( convert_f2c(s,slen) );
}

/*****************************************************************/

