
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! revision log :
!
! 10.10.2015	ggu	changed VERS_7_3_2
! 12.10.2015	ggu	changed VERS_7_3_3
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 23.09.2019	ggu	for resolution make histogram

c***************************************************************
c***************************************************************
c***************************************************************

	subroutine bas_resolution

c handles horizontal resolution

	use basin
	use basutil

	implicit none

	real cvres(nkn)

	if( .not. breadbas ) then
	  write(6,*) 'for -resol we need a bas file'
	  stop 'error stop bas_resolution: need a bas file'
	end if

	call compute_resolution(cvres)
	call write_resolution(cvres)

	end

c***************************************************************

	subroutine write_resolution(cvres)

c computes horizontal resolution

	use basin
	use mod_depth

	implicit none

	real cvres(nkn)

	integer nb,ierr
	integer ilhkv(1)
	real rmin,rmax
	real hlv(1)
	character*80 file,title

	integer ifileo

	ilhkv(1) = 1
	hlv(1) = 10000.

        call mkname(' ','basres','.nos',file)
        write(6,*) 'writing file ',file(1:50)
        nb = ifileo(0,file,'unform','new')
        if( nb .le. 0 ) goto 98

	title = 'basin resolution'
        call wfnos(nb,3,nkn,nel,1,1,title,ierr)
        if( ierr .ne. 0 ) goto 97
        call wsnos(nb,ilhkv,hlv,hev,ierr)
        if( ierr .ne. 0 ) goto 97

	call mima(cvres,nkn,rmin,rmax)
	write(6,*) 'min/max resolution: ',rmin,rmax

        call wrnos(nb,0,334,1,ilhkv,cvres,ierr)    !aver
        if( ierr .ne. 0 ) goto 99

	write(6,*)
	write(6,*) 'data written to file basres.nos'
	write(6,*)

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   97   continue
        write(6,*) 'error writing header'
        stop 'error stop nosaver'
   98   continue
        write(6,*) 'error opening file'
        stop 'error stop nosaver'
   99   continue
        write(6,*) 'error writing data'
        stop 'error stop nosaver'
	end

c***************************************************************

	subroutine compute_resolution(cvres)

c computes horizontal resolution

	use basin
	use evgeom

	implicit none

	real cvres(nkn)

	integer k,ie,ii,iii,k1,k2
	real dd
	real cvaux(nkn)

	real dist

	cvres = 0.
	cvaux = 0.

	!write(6,*) 'isphe: ',isphe_ev
	
	do ie=1,nel
	  do ii=1,3
	    iii = 1 + mod(ii,3)
	    k1 = nen3v(ii,ie)
	    k2 = nen3v(iii,ie)
	    dd = dist(k1,k2)
	    cvres(k1) = cvres(k1) + dd
	    cvaux(k1) = cvaux(k1) + 1.
	    cvres(k2) = cvres(k2) + dd
	    cvaux(k2) = cvaux(k2) + 1.
	  end do
	end do

	do k=1,nkn
	  if( cvaux(k) .gt. 0. ) cvres(k) = cvres(k) / cvaux(k)
	end do

	call make_histo(nkn,cvres)

	end

c***************************************************************


	function dist(k1,k2)

c computes distance between nodes

	use basin

	implicit none

	real dist
	integer k1,k2

	real x1,x2,y1,y2,dx,dy

	x1 = xgv(k1)
	y1 = ygv(k1)
	x2 = xgv(k2)
	y2 = ygv(k2)

	call compute_distance(x1,y1,x2,y2,dx,dy)

	dist = sqrt( dx*dx + dy*dy )

	end

c***************************************************************

	subroutine make_histo(n,a)

	implicit none

	integer n
	real a(n)

	integer i,nr,ir
	integer, allocatable :: ic(:)
	real rmin,rmax
	real rrmin,rrmax
	real rtot,dr,r,perc

	real rnext

	rmin = minval(a)
	rmax = maxval(a)

	write(6,*) 'resolution min/max: ',rmin,rmax

	rrmin = rnext(rmin,+2)
	rrmax = rnext(rmax,-2)

	write(6,*) 'resolution min/max: ',rrmin,rrmax

	rtot = rmax - rmin
	dr = rnext(rtot/20.,2)
	nr = rmax / dr

	write(6,*) 'resolution: ',dr,nr,nr*dr

	allocate(ic(0:nr+1))
	ic = 0

	do i=1,n
	  ir = a(i)/dr
	  ic(ir) = ic(ir) + 1
	end do

	rtot = sum(ic)

	open(1,file='histo.txt',status='unknown',form='formatted')
	do i=0,nr
	  r = (dr/2.) * (i+i+1)
	  perc = 100. * ic(i) / rtot
	  write(6,*) i,r,ic(i),perc
	  write(1,*) i,r,ic(i),perc
	end do
	close(1)
	write(6,*) 'histogram data written to histo.txt'

	end

c***************************************************************
c***************************************************************
c***************************************************************

