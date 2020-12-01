
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
! 10.07.2015	ggu	changed VERS_7_1_50
! 20.07.2015	ggu	changed VERS_7_1_81
! 29.09.2015	ggu	changed VERS_7_2_5
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

!==================================================================
	module mod_bndo
!==================================================================

	implicit none

	integer, private, save  :: nrb_bndo = 0
	integer, private, save  :: ngr_bndo = 0

	integer, save :: nbndo = 0	!total number of OB nodes
	integer, save :: ndebug = 0	!unit number for debug messages

	integer, save :: kbcdim = 0	!to be deleted later
	integer, save :: kopdim = 0	!to be deleted later

	integer, allocatable, save :: nopnod(:)	!number of nodes close to OB
	integer, allocatable, save :: ibcnod(:)	!number of boundary
	integer, allocatable, save :: kbcnod(:)	!number of boundary node
	integer, allocatable, save :: itynod(:)	!type of boundary
	integer, allocatable, save :: nopnodes(:,:)	!nodes close to OB
	real, allocatable, save :: xynorm(:,:)		!normal direction
	real, allocatable, save :: wopnodes(:,:)	!weights

!==================================================================
	contains
!==================================================================

	subroutine mod_bndo_init(ngr,nrb)

	integer ngr,nrb

	integer nlk,naux

        if( ngr == ngr_bndo .and. nrb == nrb_bndo ) return

        if( ngr_bndo > 0 ) then
          deallocate(nopnod)
          deallocate(ibcnod)
          deallocate(kbcnod)
          deallocate(itynod)
          deallocate(nopnodes)
          deallocate(xynorm)
          deallocate(wopnodes)
        end if

        ngr_bndo = ngr
        nrb_bndo = nrb

	kbcdim = nrb		!to be deleted later
	kopdim = ngr		!to be deleted later

	if( ngr == 0 ) return

	naux = max(1,nrb)

        allocate(nopnod(naux))
        allocate(ibcnod(naux))
        allocate(kbcnod(naux))
        allocate(itynod(naux))
        allocate(nopnodes(ngr,naux))
        allocate(xynorm(2,naux))
        allocate(wopnodes(ngr,naux))

	end subroutine mod_bndo_init

	subroutine mod_bndo_info

	write(6,*) 'mod_bndo_info ================='
	write(6,*) 'ngr_bndo: ',ngr_bndo
	write(6,*) 'nrb_bndo: ',nrb_bndo
	write(6,*) 'nbndo: ',nbndo
	write(6,*) 'mod_bndo_info end ================='

	end subroutine mod_bndo_info

!==================================================================
	end module mod_bndo
!==================================================================

