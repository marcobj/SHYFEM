
!--------------------------------------------------------------------------
!
!    Copyright (C) 1996,2001,2010,2014-2015,2019  Georg Umgiesser
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

c explicitly sets the velocities (transports) to prescribed value
c
c revision log :
c
c 06.06.1996	ggu	written (from sp159f)
c 09.11.2001	ggu	compiler directives removed
c 23.03.2010	ggu	changed v6.1.1
c 19.12.2014	ggu	changed VERS_7_0_10
c 23.12.2014	ggu	changed VERS_7_0_11
c 19.01.2015	ggu	changed VERS_7_1_3
c 17.07.2015	ggu	changed VERS_7_1_80
c 20.07.2015	ggu	changed VERS_7_1_81
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	changed VERS_7_5_61
c 21.05.2019	ggu	changed VERS_7_5_62
c
c******************************************************************

	subroutine chao

	use mod_hydro_baro
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin

	implicit none

	real umax,dz,fact
	integer ie,l,last,iex
	integer ichange

	ichange = 0
	umax = 0.1
	dz = 3.

	do ie=1,nel

	  iex = ipev(ie)
	  last = mod(iex,100)

	  if( last .le. 2 ) then
	    ichange = ichange + 1
	    do l=1,nlv
		fact = - 0.5
		if( l .eq. 1 ) fact = 1.
		if( l .eq. 2 ) fact = 0.5
		utlnv(l,ie) = umax * dz * fact
		vtlnv(l,ie) = 0.
	    end do
	    unv(ie) = 0.
	    vnv(ie) = 0.
	  end if

	end do

	write(6,*) 'chao: ',umax,nlv,ichange

	return
	end

c******************************************************************

