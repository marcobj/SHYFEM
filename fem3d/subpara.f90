!
! Paralution solver routines
!
! revision log :
!
! 28.11.2015    mbj     First version
!
!*************************************************************************

      subroutine para_init_system

! Initialize vector and matrix      

      use mod_system
      use basin
      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR, C_NULL_CHAR, C_LOC

      implicit none
      include 'param.h'

      integer icall_coo
      data icall_coo /0/
      save icall_coo

      interface
        subroutine paralution_init() BIND(C)
           use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR
        end subroutine
      end interface
      
      rvec = 0.
      raux = 0.
      coo = 0.

      if (icall_coo.eq.0) then		! only first time
	 write(6,*) 'SOLVER: Paralution'
         call coo_init(nel,nkn,mbw,nen3v,csrdim,nnzero,ijp,icoo,jcoo)
         print*, 'coo-matrix initialisation...'
         print*, 'Number of non-zeros: ',nnzero
         call paralution_init()
         icall_coo=1
      end if


      end

!*************************************************************************

 subroutine para_solve_system_coo(nin,xguess)

! Solver routine with Paralution iterative methods.

 use mod_system
 use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR, C_NULL_CHAR, C_LOC

 implicit none

 integer, intent(in) :: nin
 real, intent(in) :: xguess(nin)

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

 !Paralution parameters
 character(len=80) :: Solver,Op_mat_form,Prec,Prec_mat_form

 n = nin
 m = n	!n row and col are the same
 nnz = nnzero

 if( .not. allocated(x) ) then
   allocate( rhs(n), x(n) )
   allocate( rows(nnz), cols(nnz) )
   allocate( rval(nnz) )
 end if

 rhs = rvec
 !x = 0._C_DOUBLE
 x = xguess
 rows = icoo
 cols = jcoo
 rval = coo

! Solver (CG,BiCGStab,GMRES,Fixed-Point)
 Solver = 'BiCGStab'	!fastest
! Operation matrix format(DENSE,CSR,MCSR,COO,DIA,ELL,HYB)
 Op_mat_form = 'MCSR'	!fastest
! Preconditioner (None,Jacobi,MultiColoredGS,MultiColoredSGS,ILU,MultiColoredILU)
 Prec = 'Jacobi'
 !Prec = 'ILU'	!much slower with big matrices
! Preconditioner matrix format (DENSE,CSR,MCSR,COO,DIA,ELL,HYB)
 Prec_mat_form = 'MCSR'
 
 
 ! Run paralution C function for COO matrices
 call paralution_fortran_solve_coo( n, m, nnz,                                          &
 &                                  trim(Solver) // C_NULL_CHAR,                        &
 &                                  trim(Op_mat_form) // C_NULL_CHAR,                   &
 &                                  trim(Prec) // C_NULL_CHAR,                          &
 &                                  trim(Prec_mat_form) // C_NULL_CHAR,                 &
 &                                  C_LOC(rows), C_LOC(cols), C_LOC(rval), C_LOC(rhs),  &
 &                                  1e-8_C_DOUBLE, 1e-8_C_DOUBLE, 1e+8_C_DOUBLE, 3000,  &
 &                                  10, 0, 1, C_LOC(x), iter, resnorm, ierr )

 rvec = x

 !deallocate( rows, cols, rval, rhs, x )

 end subroutine para_solve_system_coo


!*************************************************************************

 subroutine para_solve_system_csr(nin,xguess)

! Solver routine with Paralution iterative methods.

 use mod_system
 use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR, C_NULL_CHAR, C_LOC

 implicit none

 integer, intent(in) :: nin
 real, intent(in) :: xguess(nin)

 interface
   subroutine paralution_fortran_solve_csr( n, m, nnz, solver, mformat, preconditioner, pformat,    &
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

    end subroutine paralution_fortran_solve_csr
 end interface

 integer(kind=C_INT)   :: n, m, nnz, iter, ierr
 real(kind=C_DOUBLE)   :: resnorm

 integer(kind=C_INT), allocatable, target :: rows(:), cols(:)
 real(kind=C_DOUBLE), allocatable, target :: rval(:)
 real(kind=C_DOUBLE), allocatable, target :: rhs(:), x(:)

 !Paralution parameters
 character(len=80) :: Solver,Op_mat_form,Prec,Prec_mat_form

 integer, allocatable :: iwork(:)

 n = nin
 m = n	!n row and col are the same
 nnz = nnzero

 if( .not. allocated(x) ) then
   allocate( rhs(n), x(n) )
   allocate( rows(n+1), cols(nnz) )
   allocate( rval(nnz) )
   allocate(iwork(max(n+1,2*nnz)))
 end if

! convert from coo to csr
 call coocsr(n,nnz,coo,icoo,jcoo,rval,cols,rows)
 call csort (n,rval,cols,rows,iwork,.true.)

 rhs = rvec
 !x = 0._C_DOUBLE
 x = xguess

! Solver (CG,BiCGStab,GMRES,Fixed-Point)
 Solver = 'BiCGStab'	!fastest
! Operation matrix format(DENSE,CSR,MCSR,COO,DIA,ELL,HYB)
 Op_mat_form = 'MCSR'	!fastest
! Preconditioner (None,Jacobi,MultiColoredGS,MultiColoredSGS,ILU,MultiColoredILU)
 Prec = 'Jacobi'
 !Prec = 'MultiColoredILU'	!much slower with big matrices
! Preconditioner matrix format (DENSE,CSR,MCSR,COO,DIA,ELL,HYB)
 Prec_mat_form = 'MCSR'
 
 ! Run paralution C function for CSR matrices
 call paralution_fortran_solve_csr( n, m, nnz,                                          &
 &                                  trim(Solver) // C_NULL_CHAR,                        &
 &                                  trim(Op_mat_form) // C_NULL_CHAR,                   &
 &                                  trim(Prec) // C_NULL_CHAR,                          &
 &                                  trim(Prec_mat_form) // C_NULL_CHAR,                 &
 &                                  C_LOC(rows), C_LOC(cols), C_LOC(rval), C_LOC(rhs),  &
 &                                  1e-8_C_DOUBLE, 1e-8_C_DOUBLE, 1e+8_C_DOUBLE, 3000,  &
 &                                  20, 0, 1, C_LOC(x), iter, resnorm, ierr )

 rvec = x

 !deallocate( rows, cols, rval, rhs, x , iwork)

 end subroutine para_solve_system_csr
