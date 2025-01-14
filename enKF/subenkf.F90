!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
!********************************************************

  subroutine rst_read(rstname,atimea)

  use mod_restart
  use levels, only: nlvdi,nlv,hlv,ilhv,ilhkv
  use shympi
  implicit none

  character(len=*), intent(in) :: rstname
  double precision, intent(in) :: atimea

  integer ierr,iflag
  double precision atimef

  integer, save :: icall = 0
  real*4 ibarcl4,iconz4,imerc4,iturb4,iwvert_rst4,ieco_rst4,zero4

  zero4 = 0.

  open(24,file=trim(rstname),status='old',form='unformatted',iostat=ierr)
  if (ierr /= 0) error stop 'rst_read: Error opening file'

  89  call rst_read_record(24,atimef,iflag,ierr)
     if (ierr /= 0) goto 90
     if (atimef /= atimea) goto 89
  close(24)

  if ( icall == 0 ) then

     hlv = hlvrst
     ilhv = ilhrst
     ilhkv = ilhkrst
     ibarcl4 = ibarcl_rst
     iwvert_rst4 = iwvert_rst
     ieco_rst4 = ieco_rst
     iconz4 = iconz_rst
     imerc4 = imerc_rst
     iturb4 = iturb_rst

     call addpar('ibarcl',ibarcl4)
     call addpar('iconz',iconz4)
     call addpar('iwvert',iwvert_rst4) !maybe not
     call addpar('ieco',ieco_rst4) !maybe not
     call addpar('ibio',zero4)
     call addpar('ibfm',zero4)
     call addpar('imerc',imerc4)
     call addpar('iturb',iturb4)

     call daddpar('date',0.)
     call daddpar('time',0.)

     write(*,*) 'SHYFEM flags from restart:'
     write(*,*) 'ibarcl	= ',ibarcl_rst
     write(*,*) 'iconz	= ',iconz_rst
     write(*,*) 'iwvert	= ',iwvert_rst
     write(*,*) 'ieco	= ',ieco_rst
     write(*,*) 'imerc	= ',imerc_rst
     write(*,*) 'iturb	= ',iturb_rst
  end if

  icall = icall + 1

  return

 90 write(*,*) 'Error in the restart file. Are you sure that the analysis step'
    write(*,*) 'has a time present in the restart records?'
    error stop


  end subroutine rst_read

!********************************************************

  subroutine rst_write(rstname,atimea)

  implicit none

  character(len=*), intent(in) :: rstname
  double precision, intent(in) :: atimea
  real*4 :: svar

  open(34,file=rstname,form='unformatted')
  call rst_write_record(atimea,34)
  close(34)

  end subroutine rst_write

!********************************************************

  subroutine find_node(x,y,ie,ik)

  use shyfile
  use basin

  implicit none

  real, intent(in) :: x,y
  integer, intent(in) :: ie
  integer, intent(out) :: ik
  
  integer i,ikmin
  real d,dmin

  ikmin = -999
  dmin=10e10
  do i = 1,3
     ik = nen3v(i,ie)
     d = ( x - xgv(ik) )**2 + ( y - ygv(ik) )**2
     if (d <= dmin) then
       dmin = d
       ikmin = ik
     end if
  end do
  ik = ikmin

  end subroutine find_node

!********************************************************

   subroutine layer_thick(nelem)
   ! Set zenv from znv
   ! A value bigger than the depth hm3v is set
   use mod_hydro
   use basin

   implicit none
   integer,intent(in) :: nelem
   integer ie,ii,k
   real*4 z,h
   real*4 hmin

   hmin = 0.03
   
   do ie = 1,nelem 
     do ii = 1,3
        k = nen3v(ii,ie)
        zenv(ii,ie) = znv(k)

        z = zenv(ii,ie)
        h = hm3v(ii,ie)
        if (z + h < hmin) then
           z = -h + hmin
           zenv(ii,ie) = z
        end if

     end do
   end do
   
   end subroutine layer_thick

!********************************************************

  subroutine num2str(num,str)

  implicit none

  integer, intent(in) :: num
  character(len=5), intent(out) :: str

  if ((num >= 0).and.(num < 10)) then
    write(str,'(a4,i1)') '0000',num
  elseif ((num >= 10).and.(num < 100)) then
    write(str,'(a3,i2)') '000',num
  elseif ((num >= 100).and.(num < 1000)) then
    write(str,'(a2,i3)') '00',num
  elseif ((num >= 1000).and.(num < 10000)) then
    write(str,'(a1,i4)') '0',num
  elseif ((num >= 10000).and.(num < 100000)) then
    write(str,'(i5)') num
  else
    error stop 'num2str: num out of range'
  end if

  end subroutine num2str

!********************************************************

  subroutine random_vec(v,vdim)

  use m_random

  implicit none
  integer vdim
  real v(vdim),vaux(vdim-1)
  real aaux,ave
  integer n

  call random(vaux,vdim-1)
  ! remove outlayers
  do n = 1,vdim-1
     aaux = vaux(n)
     if (abs(aaux) >= 3.) then
       aaux = aaux/abs(aaux) * (abs(aaux)-floor(abs(aaux)) + 1.) 
     end if
     vaux(n) = aaux
  end do

  ! set mean eq to zero
  ave = sum(vaux)/float(vdim-1)
  vaux = vaux - ave

  v(1) = 0.
  v(2:vdim) = vaux

  end subroutine random_vec

!********************************************************

  subroutine make_0Dpert(vflag,n,na,id,vec,t,tau)

  implicit none

  character(len=5), intent(in) :: vflag
  integer, intent(in) :: n,na,id
  real, intent(out) :: vec(n)
  double precision, intent(in) :: t,tau
  character(len=5) :: nal,idl
  character(len=21) :: pfile
  logical bfile
  integer nf
  real, allocatable :: vec_old(:)
  double precision t_old
  double precision alpha,dt
  character(len=1) :: vfl

  ! make a new perturbation
  !
  call random_vec(vec,n)

  ! white noise if tau is lower than 0
  !
  if (tau <= 0.) return

  if (vflag == '0DLEV') then
     vfl = 'z'
  else if (vflag == '0DTEM') then
     vfl = 't'
  else if (vflag == '0DSAL') then
     vfl = 's'
  end if

  ! if exist load an old perturbation and merge
  !
  call num2str(na-1,nal)
  call num2str(id,idl)
  pfile = vfl // 'pert_' // nal // '_' // idl // '.bin'
  inquire(file=pfile,exist=bfile)
  if (bfile) then
     open(22,file=pfile,status='old',form='unformatted')
     read(22) nf
     if (nf /= n) error stop 'make_0Dpert: dimension mismatch'
     read(22) t_old
     allocate(vec_old(nf))
     read(22) vec_old
     close(22)

     dt = t - t_old
     alpha = 1. - (dt/tau) 

     ! when the old time is far compared to tau
     if (alpha < 0.) alpha = 0.

     vec = alpha * vec_old + sqrt(1 - alpha**2) * vec
  else
     !write(*,*) 'No old file with perturbations: ',pfile
     continue
  end if

  ! save the last perturbation
  !
  call num2str(na,nal)
  pfile = vfl // 'pert_' // nal // '_' // idl // '.bin'
  open(32,file=pfile,form='unformatted')
  write(32) n
  write(32) t
  write(32) vec
  close(32)

  end subroutine make_0Dpert

!********************************************************

  subroutine make_2Dpert(vec,n,nens)

  use basin
  use m_sample2D
  use mod_para

  implicit none

  integer, intent(in) :: n,nens
  real, intent(out) :: vec(n,nens)

  integer :: nx,ny	!number of grid points in x and y direction for the regular grid
  real x1,x2,y1,y2,xlength,ylength
  real dx,dy,rx,ry
  real*4 sdx,sdy,sx0,sy0
  real*4 :: sflag = -999.

  real, allocatable :: mat(:,:,:)
  real*4, allocatable :: mat4(:,:),vec4fem(:)

  integer ne

  !----------------------------------------------------
  ! Computes some geometric parameters from the grid
  !----------------------------------------------------
  x1 = minval(xgv)
  x2 = maxval(xgv)
  y1 = minval(ygv)
  y2 = maxval(ygv)

  x1 = floor(x1)
  y1 = floor(y1)
  x2 = ceiling(x2)
  y2 = ceiling(y2)

  xlength = x2 - x1
  ylength = y2 - y1

  if ((xlength > 180.).or.(ylength > 90.)) error stop 'Coordinates are not geographical'

  if (xlength < 4) then
     dx = 0.05
     dy = 0.05
     rx = 2
     ry = 2
     nx = xlength/dx + 1
     ny = ylength/dy + 1
  else
     dx = 0.1
     dy = 0.1
     rx = 4
     ry = 4
     nx = xlength/dx + 1
     ny = ylength/dy + 1
  end if

  if (verbose) then
    write(*,'(a20,2f8.4,i5,f8.4)') 'x1,xlength,nx,dx: ',x1,xlength,nx,dx
    write(*,'(a20,2f8.4,i5,f8.4)') 'y1,ylength,ny,dy: ',y1,ylength,ny,dy
    write(*,'(a14,2f10.4,1x,i3)') 'rx,ry,fmult: ',rx,ry,fmult_init
  end if

  !----------------------------------------------------
  ! creates the sample
  !----------------------------------------------------
  allocate(mat(nx,ny,nens))
  call sample2D(mat,nx,ny,nens,fmult_init,dx,dy,rx,ry,theta_init,sample_fix_init,verbose)

  !----------------------------------------------------
  ! Interpolates over the FEM grid
  !----------------------------------------------------
  write(*,*) 'Interpolating 2D field over the FEM grid...'
  sdx = dx
  sdy = dy
  sx0 = x1
  sy0 = y1
  call setgeo(sx0,sy0,sdx,sdy,sflag)

  allocate(mat4(nx,ny),vec4fem(n))
  do ne = 1,nens
    mat4 = mat(:,:,ne)
    call am2av(mat4,vec4fem,nx,ny)
    vec(:,ne) = vec4fem
  end do
  deallocate(mat,mat4,vec4fem)

  end subroutine make_2Dpert

!********************************************************

  subroutine find_el_node(x,y,iie,ik)

  use basin

  implicit none

  real, intent(in) :: x,y
  integer, intent(out) :: iie,ik
  real*4 x4,y4
  real dst, dstmax
  integer iik,ii
  integer k

  !-----------
  ! Finds the grid element and the node 
  ! nearest to the observation (the subroutine is in real4)
  !-----------
  x4 = x
  y4 = y
  call find_element(x4,y4,iie)

  dstmax = 1e15
  ! Inside the grid
  if (iie /= 0) then
     dst = 0.
     do ii = 1,3
        iik = nen3v(ii,iie)
        dst = sqrt( (xgv(iik)-x4)**2 + (ygv(iik)-y4)**2 )
        if (dst < dstmax) then
          dstmax = dst
          ik = iik
        end if
     end do
  ! outside the grid
  else
     !write(*,*) 'Warning! Observation outside the grid'
     dst = 0.
     do k = 1,nkn
        dst = sqrt( (xgv(k)-x4)**2 + (ygv(k)-y4)**2 )
	if (dst < dstmax) then
          dstmax = dst
          ik = k
        end if
     end do
     call find_element(xgv(ik),ygv(ik),iie)
     !write(*,*) 'Distance, node, element: ',dstmax, ik, iie
  end if

  end subroutine find_el_node

!*************************************************************

  subroutine find_weight_GC(rho,dst,w)
  implicit none
  real, intent(in) :: rho,dst
  real, intent(out) :: w
  real s

  s = dst/rho

  ! Taper function of Gaspari-Cohn
  if ((s >=0) .and. (s<1)) then
     w = 1 - (5./3. * s**2) + (5./8. * s**3)  + (1./2. * s**4) - (1./4. * s**5)
  else if ((s>=1) .and. (s<2)) then
     w = - (2./3. * s**-1) + 4 - (5. * s) + (5./3. * s**2) + (5./8. * s**3)& 
         - (1./2. * s**4) + (1./12. * s**5)
  else
     w = 0.
  end if

  end subroutine find_weight_GC


!********************************************************
  
  subroutine superobs_horiz_el(no,x,y,ostatus,val1,val2)

  use basin

  implicit none

  integer,intent(in) :: no
  real,intent(in) :: x(no),y(no)
  integer,intent(inout) :: ostatus(no)
  real,intent(inout) :: val1(no),val2(no)
  integer nobs(nel)
  integer, allocatable :: ieobs(:)
  integer, allocatable :: oindex(:,:)
  integer omax
  real*4 x4,y4
  integer n,ie,nn
  real avval

  ! find the index of the elements containing obs
  !
  allocate(ieobs(no))
  ieobs = -999
  nobs = 0
  do n = 1,no
     if (ostatus(n) > 0) cycle

     x4 = x(no)
     y4 = y(no)
     call find_element(x4,y4,ie)
     ieobs(n) = ie
     nobs(ie) = nobs(ie) + 1
  end do

  ! find the maximum number of obs inside an element
  !
  omax = 0
  do ie = 1,nel
     omax = max(omax,nobs(ie))
  end do

  ! find the final obs indexes
  !
  allocate(oindex(0:omax,nel))
  oindex = 0
  do n = 1,no
     if (ostatus(n) > 0) cycle

     ie = ieobs(n)
     nn = oindex(0,ie)
     nn = nn + 1
     oindex(nn,ie) = n
     oindex(0,ie) = nn
  end do

  ! average the values and assign the status
  !
  do ie = 1,nel
     nn = oindex(0,ie)
     if (nn == 0) cycle
     write(*,*) 'making super-observation...'
     avval = sum(val1(oindex(1:nn,ie)))/nn
     val1(oindex(1,ie)) = avval
     avval = sum(val2(oindex(1:nn,ie)))/nn
     val2(oindex(1,ie)) = avval
     ostatus(oindex(1,ie)) = 1
     ostatus(oindex(2:nn,ie)) = 2
  end do
  
  end subroutine superobs_horiz_el

!********************************************************

  subroutine save_X5(alabel,tt)

  implicit none
  character, intent(in) :: alabel*6
  double precision, intent(in) :: tt

  integer :: nrens
  integer :: nrobs
  real, allocatable :: X5(:,:)
  real, allocatable :: X3(:,:)
  real, allocatable :: S(:,:)
  character(len=2) :: tag2

  ! read X5.uf depending on the mult method (X5 or X3,S)
  open(unit=25,file='X5.uf',form='unformatted',status='old')
   read(25) tag2
   rewind 25
   if (tag2 == 'X3') then
      read(25) tag2,nrens,nrobs
      allocate(X3(nrens,nrens),S(nrobs,nrens))
      rewind 25
      read(25) tag2,nrens,nrobs,X3,S
   else
      read(25) tag2,nrens
      allocate(X5(nrens,nrens))
      rewind 25
      read(25) tag2,nrens,X5
   end if
  close(25)
  
  ! write the total X5_tot.uf adding time and type of analysis (global or local)
  write(*,*) 'saving X5 matrix for enKS'
  open(unit=35,file='X5_tot.uf',form='unformatted',status='unknown',access='append')
   write(35) tt,alabel,tag2
   if (tag2 == 'X3') then
      write(35) nrens,nrobs
      write(35) X3,S
   else
      write(35) nrens
      write(35) X5
   end if
  close(35)

  end subroutine save_X5


