!
! Paralution solver routines
!
! revision log :
!
! 27.11.2015    mbj     First version
!
!*************************************************************************

      subroutine para_init_system

! Initialize vector and matrix      

	use mod_system
	use basin

      implicit none
      include 'param.h'

      integer n

      integer icall_coo
      data icall_coo /0/
      save icall_coo

      do n=1,nkn
         rvec(n) = 0.
         raux(n) = 0.
      end do

      if (icall_coo.eq.0) then		! only first time
	 write(6,*) 'SOLVER: Paralution'
         call coo_init(nel,nkn,mbw,nen3v,csrdim,nnzero,ijp,icoo,jcoo)
         print*, 'coo-matrix initialisation...'
         print*, 'Number of non-zeros: ',nnzero
         icall_coo=1
      end if

      do n=1,nnzero
         coo(n) = 0.
      end do

      end

!*************************************************************************

 subroutine para_solve_system(nin)

! Solver routine with Paralution iterative methods.

 use mod_system
 use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR, C_NULL_CHAR, C_LOC

 implicit none

 integer, intent(in) :: nin

 interface
   subroutine paralution_fortran_solve_coo( n, m, nnz, solver, mformat, preconditioner, pformat,    &
    &                                        rows, cols, rval, rhs, atol, rtol, div, maxiter, basis, &
    &                                        p, q, x, iter, resnorm, ierr ) BIND(C)

      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR

      integer(kind=C_INT), value, intent(in)  :: n, m, nnz, maxiter, basis, p, q
      real(kind=C_DOUBLE), value, intent(in)  :: atol, rtol, div
      integer(kind=C_INT),        intent(out) :: iter, ierr
      real(kind=C_DOUBLE),        intent(out) :: resnorm
      type(C_PTR),         value, intent(in)  :: rows, cols, rval, rhs
      type(C_PTR),         value              :: x
      character(kind=C_CHAR)                  :: solver, mformat, preconditioner, pformat

    end subroutine paralution_fortran_solve_coo
 end interface

 integer(kind=C_INT)   :: n, m, nnz, iter, ierr
 real(kind=C_DOUBLE)   :: resnorm

 integer(kind=C_INT), allocatable, target :: rows(:), cols(:)
 real(kind=C_DOUBLE), allocatable, target :: rval(:)
 real(kind=C_DOUBLE), allocatable, target :: rhs(:), x(:)

 n = nin
 m = n	!n row and col are the same
 nnz = nnzero

 allocate( rhs(n), x(n) )
 allocate( rows(nnz), cols(nnz) )
 allocate( rval(nnz) )
 rhs = rvec
 x = raux
 rows = icoo
 cols = jcoo
 rval = coo

  ! Run paralution C function for COO matrices
  ! Doing a GMRES with MultiColored ILU(1,2) preconditioner
  ! Check paralution documentation for a detailed argument explanation
  call paralution_fortran_solve_coo( n, m, nnz,                                          &
  &                                  'CG' // C_NULL_CHAR,                                &
  &                                  'CSR' // C_NULL_CHAR,                               &
  &                                  'MultiColoredILU' // C_NULL_CHAR,                   &
  &                                  'CSR' // C_NULL_CHAR,                               &
  &                                  C_LOC(rows), C_LOC(cols), C_LOC(rval), C_LOC(rhs),  &
  &                                  1e-15_C_DOUBLE, 1e-8_C_DOUBLE, 1e+8_C_DOUBLE, 5000, &
  &                                  30, 0, 1, C_LOC(x), iter, resnorm, ierr )

 rvec = x

 deallocate( rows, cols, rval, rhs, x )


 end subroutine para_solve_system

!*************************************************************************
