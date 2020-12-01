
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013-2015,2019  Georg Umgiesser
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
! 13.06.2013	ggu	changed VERS_6_1_65
! 05.11.2014	ggu	changed VERS_7_0_5
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_52
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62

!**************************************************************************

	program check_debug

c checks two files written with check_debug from ht

	implicit none

	integer ndim
	parameter (ndim=2000000)

	character*60 name_one,name_two
	character*80 text1,text2,text
	logical bcheck,bstop,bverbose
	integer it1,it2,it
	integer nt1,nt2,nt
	integer nf1,nf2,nf
	integer nrec
	integer i,idiff,idiff_tot
	integer nc

	real val1(ndim)
	real val2(ndim)

	bstop = .false.			!stop on error
	bstop = .true.			!stop on error
	bverbose = .true.		!only write error
	bverbose = .false.		!only write error
	bcheck = .true.			!check for differences

	nc = command_argument_count()
	if( nc .ne. 2 ) then
	  write(6,*) 'Usage: check_debug file1 file2'
	  stop 'error stop check_debug: no files given'
	end if

	call get_command_argument(1,name_one)
	call get_command_argument(2,name_two)

	open(1,file=name_one,status='old',form='unformatted')
	open(2,file=name_two,status='old',form='unformatted')

	write(6,*) 'file 1: ',trim(name_one)
	write(6,*) 'file 2: ',trim(name_two)

	idiff_tot = 0

	do while(.true.)

	  read(1,end=9) it1
	  read(2,end=9) it2
	  if( it1 .ne. it2 ) goto 99
	  it = it1
	  write(6,*) 'time = ',it

	  nrec = 0
	  do while(.true.)
	    read(1) nt1,nf1
	    read(2) nt2,nf2
	    nrec = nrec + 1
	    if( nt1 .ne. nt2 ) goto 98
	    if( nf1 .ne. nf2 ) goto 98
	    nt = nt1
	    nf = nf1

	    if( nt .eq. 0 ) exit
	    if( nt .gt. ndim ) goto 97

	    read(1,end=9) text1
	    read(2,end=9) text2
	    if( text1 .ne. text2 ) goto 96
	    text = text1

	    read(1) (val1(i),i=1,nt)
	    read(2) (val2(i),i=1,nt)

	    idiff = 0
	    !bcheck = it .ge. 87300
	    !bcheck = bcheck .and. nrec .ne. 5
	    if( bcheck ) then
	      call check_val(it,nrec,nt,nf,val1,val2,idiff)
	      if( idiff > 0 .or. bverbose ) then
	        write(6,*) trim(text),nrec,nt,nf,idiff
	      end if
	      idiff_tot = idiff_tot + idiff
	    end if

	  end do

	  if( bstop .and. idiff_tot > 0 ) exit
	  if( bverbose ) write(6,*) 'nrecs checked: ',nrec
	end do

    9	continue

	close(1)
	close(2)

	write(6,*) 'total differences found: ',idiff_tot

	stop
   99	continue
	write(6,*) it1,it2
	stop 'error stop check_debug: time mismatch'
   98	continue
	write(6,*) nt1,nt2,nf1,nf2
	stop 'error stop check_debug: size mismatch'
   97	continue
	write(6,*) it,nrec,nt,ndim
	stop 'error stop check_debug: dimension'
   96	continue
	write(6,*) trim(text1),' - ',trim(text2)
	stop 'error stop check_debug: text mismatch'
	end

c*******************************************************************

	subroutine check_val(it,nrec,nt,nf,val1,val2,idiff)

	implicit none

	integer it,nrec
	integer nt,nf,idiff
	real val1(nt)
	real val2(nt)

	integer i,k,l

	idiff = 0

	do i=1,nt
	  if( val1(i) .ne. val2(i) ) then
	    k = 1 + (i-1)/nf
	    l = 1 + mod(i-1,nf)
	    write(77,*) it,nrec,k,l,val1(i),val2(i)
	    idiff = idiff + 1
	  end if
	end do

	end

c*******************************************************************

