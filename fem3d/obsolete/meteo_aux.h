
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
! 19.01.2015	ggu	changed VERS_7_1_3
! 10.07.2015	ggu	changed VERS_7_1_50
! 16.02.2019	ggu	changed VERS_7_5_60

        real wxv(nkndim),wyv(nkndim)
        common /wxv/wxv,/wyv/wyv
	real ppv(nkndim)
	common /ppv/ppv

        real metrad(nkndim),methum(nkndim)
        real mettair(nkndim),metcc(nkndim)
        common /metrad/metrad, /methum/methum
        common /mettair/mettair, /metcc/metcc

        real metrain(nkndim)
        common /metrain/metrain

        real tauxnv(nkndim),tauynv(nkndim)
        common /tauxnv/tauxnv,/tauynv/tauynv

        real metwbt(nkndim),metws(nkndim)
        common /metwbt/metwbt, /metws/metws

	save /metwbt/,/metws/,/metrain/
	save /ppv/,/wxv/,/wyv/,/tauxnv/,/tauynv/
	save /metrad/,/methum/,/mettair/,/metcc/

