
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016-2017,2019-2020  Georg Umgiesser
!    Copyright (C) 2018  Christian Ferrarin
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

! netcdf routines for shyelab: elab_nc
!
! revision log :
!
! 07.06.2016	ggu	changed VERS_7_5_12
! 11.07.2017	ggu	changed VERS_7_5_30
! 05.12.2017	ggu	changed VERS_7_5_39
! 25.10.2018	ccf	write field already on regular grid
! 16.02.2019	ggu	changed VERS_7_5_60
! 14.05.2019	ggu	changes in nc_output_record()
! 16.05.2019	ggu	new aux variable to write nc_records (rncaux)
! 16.05.2019	ggu	also hydro output fixed
! 29.01.2020	ggu	always write variables as 3d
!
!************************************************************

	subroutine nc_output_init(ncid,title,nvar,ivars)

	use basin
	use levels
	use mod_depth
	use shyelab_out
        use elabutil

	implicit none

	integer ncid			!id of file (return)
	character*(*) title		!name of simulation
	integer nvar			!total number of variables to be written
	integer ivars(nvar)		!variable id of SHYFEM

	logical bhydro
	integer date0,time0
	integer lmax,i,iztype,idim,ivar,ncnlv
	real, save :: ncflag = -999.

	call dts_get_date(date0,time0)
	call compute_iztype(iztype)
	call nc_set_quiet(bquiet)

        ncnlv = nlv
	if ( b2d ) ncnlv = 1

	if( breg ) then
	  allocate(value2d(nxreg,nyreg))
	  allocate(value3d(ncnlv,nxreg,nyreg))
	  allocate(vnc3d(nxreg,nyreg,ncnlv))
	  call get_lmax_reg(nxreg,nyreg,fmreg,ilhv,lmax)
	  lmax = max(1,lmax)	!at least one layer
	  if ( b2d ) lmax = 1
	  lmaxreg = lmax
	  call nc_open_reg(ncid,nxreg,nyreg,lmax
     +				,ncflag,date0,time0,iztype)
	else
	  allocate(var3d(ncnlv*nkn))
	  call nc_open_fem(ncid,nkn,nel,ncnlv,date0,time0,iztype)
	end if

	call nc_global(ncid,title)

	bhydro = ivars(1) == 1		!this is a hydro file
	allocate(var_ids(nvar))
	var_ids = 0

        do i=1,nvar
	  ivar = ivars(i)
	  idim = 3
	  if( bhydro ) then
	    if( i == 1 ) idim = 2
	    if( i == 2 ) cycle	!do not write second level
	    if( i == 3 ) ivar = 2
	  end if
          call nc_init_variable(ncid,breg,idim,ivar,ncflag,var_ids(i))
        end do

        call nc_end_define(ncid)

        if( breg ) then
          call nc_write_coords_reg(ncid,nxreg,nyreg,ncnlv
     +					,xlon,ylat,hcoord,hlv)
        else
          call nc_write_coords_fem(ncid,nkn,nel,ncnlv,xgv,ygv
     +					,hkv,nen3v,hlv)
        end if

	end

!********************************************************************

	subroutine nc_output_record(ncid,var_id,np,sv)

	use basin
	use levels
	use shyelab_out

	implicit none

	integer ncid
	integer var_id
	integer np
	!real cv3(nlvdi,nkn)
	real sv(nlvdi,np)

	integer lmax,nx,ny,iwrite
	real, allocatable :: rncaux(:,:,:)

	iwrite = iwrite_nc

        if( breg ) then
	  lmax = lmaxreg
	  nx = nxreg
	  ny = nyreg
	  allocate(rncaux(nx,ny,lmax))
	  if( nx*ny /= np ) stop 'error stop nc_output_record: breg'
          call nc_rewrite_3d_reg(nlvdi,lmax,nx,ny,sv,rncaux)
          call nc_write_data_3d_reg(ncid,var_id,iwrite
     +					,lmax,nx,ny,rncaux)
	else
	  if( nkn /= np ) stop 'error stop nc_output_record: nkn'
          call nc_compact_3d(nlv,nlv,nkn,sv,var3d)
          call nc_write_data_3d(ncid,var_id,iwrite,nlv,nkn,var3d)
        end if

	end

!********************************************************************

	subroutine nc_output_hydro(ncid,znv,uprv,vprv)

	use basin
	use levels
	use shyelab_out

	implicit none

	integer ncid
	real znv(nkn)
	real uprv(nlvdi,nkn)
	real vprv(nlvdi,nkn)

	integer lmax,nx,ny,var_id,iwrite
	real, allocatable :: rncaux(:,:,:)

	iwrite = iwrite_nc

        if( breg ) then
	  lmax = lmaxreg
	  nx = nxreg
	  ny = nyreg
	  allocate(rncaux(nx,ny,lmax))
	  !write(6,*) 'nc_output_hydro: ',nx,ny,lmax,nlv,nlvdi

	  var_id = var_ids(1)
          call fm2am2d(znv,nx,ny,fmreg,value2d)
          call nc_write_data_2d_reg(ncid,var_id,iwrite,nx,ny,value2d)

	  var_id = var_ids(3)
          call fm2am3d(nlv,ilhv,uprv,nlvdi,nx,ny,fmreg,value3d)
          call nc_rewrite_3d_reg(nlvdi,lmax,nx,ny,value3d,rncaux)
          call nc_write_data_3d_reg(ncid,var_id,iwrite
     +				,lmax,nx,ny,rncaux)

	  var_id = var_ids(4)
          call fm2am3d(nlv,ilhv,vprv,nlvdi,nx,ny,fmreg,value3d)
          call nc_rewrite_3d_reg(nlvdi,lmax,nx,ny,value3d,rncaux)
          call nc_write_data_3d_reg(ncid,var_id,iwrite
     +				,lmax,nx,ny,rncaux)
        else
	  var_id = var_ids(1)
          call nc_write_data_2d(ncid,var_id,iwrite,nkn,znv)

	  var_id = var_ids(3)
          call nc_compact_3d(nlv,nlv,nkn,uprv,var3d)
          call nc_write_data_3d(ncid,var_id,iwrite,nlv,nkn,var3d)

	  var_id = var_ids(4)
          call nc_compact_3d(nlv,nlv,nkn,vprv,var3d)
          call nc_write_data_3d(ncid,var_id,iwrite,nlv,nkn,var3d)
        end if

	end

!********************************************************************

	subroutine nc_output_time(ncid,dtime)

	use shyelab_out

	implicit none

	integer ncid
	double precision dtime

	integer, save :: icall = 0
	double precision, save :: dtime_old = 0

	if( icall == 0 ) dtime_old = dtime - 1.
	icall = icall + 1
	if( dtime_old == dtime ) return			!already written

	iwrite_nc = iwrite_nc + 1

	call nc_write_dtime(ncid,iwrite_nc,dtime)

	end

!********************************************************************

	subroutine nc_output_final(ncid)

	implicit none

	integer ncid

	call nc_close(ncid)

	end

!********************************************************************

