
!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
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
! 25.05.2015    ggu     module introduced

!==================================================================
        module levels
!==================================================================

        implicit none

        integer, save :: nlv = 0
        integer, save :: nlvdi = 0

        integer, save, private :: nkn_levels = 0
        integer, save, private :: nel_levels = 0
        integer, save, private :: nlv_levels = 0

        integer, save, allocatable :: ilhv(:)
        integer, save, allocatable :: ilhkv(:)
        integer, save, allocatable :: ilmv(:)
        integer, save, allocatable :: ilmkv(:)

        real, save, allocatable :: hlv(:)
        real, save, allocatable :: hldv(:)

!==================================================================
        contains
!==================================================================

	subroutine levels_init(nkn,nel,nl)

	integer nkn,nel,nl

        if( nkn == nkn_levels .and. nel == nel_levels .and.
     +		nl == nlv_levels ) return

        if( nkn > 0 .or. nel > 0 .or. nl > 0 ) then
          if( nkn == 0 .or. nel == 0 .or. nl == 0 ) then
            write(6,*) 'nkn,nel,nlv: ',nkn,nel,nl
            stop 'error stop levels_init: incompatible params'
          end if
        end if

	if( nlv_levels > 0 ) then
	  deallocate(ilhv)
	  deallocate(ilmv)
	  deallocate(ilhkv)
	  deallocate(ilmkv)
	  deallocate(hlv)
	  deallocate(hldv)
	end if

	nlvdi = nl
	nlv_levels = nl
	nkn_levels = nkn
	nel_levels = nel

	if( nl == 0 ) return

	allocate(ilhv(nel))
	allocate(ilmv(nel))
	allocate(ilhkv(nkn))
	allocate(ilmkv(nkn))
	allocate(hlv(nl))
	allocate(hldv(nl))
	
	hlv = 0.
	hldv = 0.

	!write(6,*) 'levels allocated: ',nkn,nel,nl

	end subroutine levels_init

!******************************************************************

	subroutine levels_init_2d(nkn,nel)

	integer nkn,nel,nl

	call levels_init(nkn,nel,1)

	ilhv = 1
	ilmv = 1
	ilhkv = 1
	ilmkv = 1
	hlv = 10000.
	hldv = 10000.

	end subroutine levels_init_2d

!******************************************************************

	subroutine levels_hlv_init(nl)

! allocates only hlv

	integer nl

	if( nlv_levels == nl ) return

        if( nlv_levels > 0 ) then
          deallocate(hlv)
        end if

	nlvdi = nl
        nlv_levels = nl

        if( nl == 0 ) return

        allocate(hlv(nl))

	end subroutine levels_hlv_init

!******************************************************************

	subroutine levels_reinit(nl)

! re-allocates arrays depending on nl

	integer nl

	integer n
	real, allocatable :: hlv_aux(:)
	real, allocatable :: hldv_aux(:)

	if( nlv_levels == nl ) return

	n = min(nl,nlv_levels)
	write(6,*) 'levels_reinit: ',nl,nlv_levels,n

	if( nlv_levels > 0 ) then
	  if( nl > 0 ) then
	    allocate(hlv_aux(nl))
	    allocate(hldv_aux(nl))
	    hlv_aux = 0.
	    hldv_aux = 0.
	    hlv_aux(1:n) = hlv(1:n)
	    hldv_aux(1:n) = hldv(1:n)
	  end if
	  deallocate(hlv)
	  deallocate(hldv)
	end if

	write(6,*) 'levels_reinit: ',nl,nlv_levels,n
	nlvdi = nl
	nlv_levels = nl

	if( nl > 0 ) then
	  allocate(hlv(nl))
	  allocate(hldv(nl))
	  hlv=0.
	  hldv=0.
	  if( n > 0 ) then
	    hlv(1:n) = hlv_aux(1:n)
	    hldv(1:n) = hldv_aux(1:n)
	    deallocate(hlv_aux)
	    deallocate(hldv_aux)
	  end if
	end if
	
	end subroutine levels_reinit

!******************************************************************

	subroutine levels_get_dimension(nl)

	integer nl

	nl = nlvdi

	end subroutine levels_get_dimension

!==================================================================
        end module levels
!==================================================================

	function get_max_node_level(k)

	use levels

	implicit none

	integer get_max_node_level
	integer k

	get_max_node_level = ilhkv(k)

	end function get_max_node_level

!******************************************************************

