
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
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
! 26.04.2018	ggu	changed VERS_7_5_46
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 11.11.2020	clr	changed time to shympitime
! 11.03.2021	clr&ggu	integrated into petsc branch

!***************************************************************
!
! handles timing information
!
!***************************************************************

!===============================================================
	module shympi_time
!===============================================================

	implicit none

	integer, parameter :: ndim = 10

	double precision, save :: shympitime(ndim) = 0.

        integer, parameter :: shympi_t_solve           =1
        integer, parameter :: shympi_t_init_solver     =2
!===============================================================
	end module shympi_time
!===============================================================

	subroutine shympi_time_reset(id)

	use shympi_time

	implicit none

	integer id

	if( id < 1 .or. id > ndim ) call shympi_time_error(id)

	shympitime(id) = 0.

	end

!***************************************************************

	subroutine shympi_time_accum(id,dt)

	use shympi_time

	implicit none

	integer id
	double precision dt

	if( id < 1 .or. id > ndim ) call shympi_time_error(id)

	shympitime(id) = shympitime(id) + dt

	end

!***************************************************************

	subroutine shympi_time_get(id,at)

	use shympi_time

	implicit none

	integer id
	double precision at

	if( id < 1 .or. id > ndim ) call shympi_time_error(id)

	at = shympitime(id)

	end

!***************************************************************

	subroutine shympi_time_error(id)

	use shympi_time

	implicit none

	integer id

	write(6,*) 'id = ',id,'  ndim = ',ndim
	stop 'error stop shympi_time_error: id out of range'

	end

!***************************************************************

