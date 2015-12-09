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
      
      integer, save :: icall_coo = 0

      interface
        subroutine paralution_init() BIND(C)
           use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR
        end subroutine
      end interface
      
      if (icall_coo.eq.0) then		! only first time

         call coo_init_new

         if( n2zero > n2max ) then
           stop 'error stop para_init_system: non zero 2d max'
         end if
         if( n3zero > n3max ) then
           stop 'error stop para_init_system: non zero 3d max'
         end if

         call paralution_init()

	 write(6,*) 'SOLVER: Paralution'
         print*, 'coo-matrix initialisation...'
         print*, 'Number of non-zeros: ',n2zero,n2max

         icall_coo=1
      end if

      rvec = 0.
      raux = 0.
      c2coo = 0.
      c3coo = 0.

      end

!*************************************************************************

 subroutine para_solve_system_coo(nnz,nin,xguess)

! Solver routine with Paralution iterative methods.

 use mod_system
 use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR, C_NULL_CHAR, C_LOC

 implicit none

 integer, intent(in) :: nnz,nin
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

 integer(kind=C_INT)   :: n, m, iter, ierr
 real(kind=C_DOUBLE)   :: resnorm

 integer(kind=C_INT), allocatable, target :: rows(:), cols(:)
 real(kind=C_DOUBLE), allocatable, target :: rval(:)
 real(kind=C_DOUBLE), allocatable, target :: rhs(:), x(:)

 !Paralution parameters
 character(len=80) :: Solver,Op_mat_form,Prec,Prec_mat_form

 n = nin
 m = n	!n row and col are the same

 if( .not. allocated(x) ) then
   allocate( rhs(n), x(n) )
   allocate( rows(nnz), cols(nnz) )
   allocate( rval(nnz) )
 end if

 rhs = rvec
 !x = 0._C_DOUBLE
 x = xguess
 rows = i2coo
 cols = j2coo
 rval = c2coo

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

 subroutine para_solve_system_csr(nnz,nin,xguess)

! Solver routine with Paralution iterative methods.

 use mod_system
 use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR, C_NULL_CHAR, C_LOC

 implicit none

 integer, intent(in) :: nnz, nin
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

 integer(kind=C_INT)   :: n, m, iter, ierr
 real(kind=C_DOUBLE)   :: resnorm

 integer(kind=C_INT), allocatable, target :: rows(:), cols(:)
 real(kind=C_DOUBLE), allocatable, target :: rval(:)
 real(kind=C_DOUBLE), allocatable, target :: rhs(:), x(:)

 !Paralution parameters
 character(len=80) :: Solver,Op_mat_form,Prec,Prec_mat_form

 integer, allocatable :: iwork(:)

 n = nin
 m = n	!n row and col are the same

 if( .not. allocated(x) ) then
   allocate( rhs(n), x(n) )
   allocate( rows(n+1), cols(nnz) )
   allocate( rval(nnz) )
   allocate(iwork(max(n+1,2*nnz)))
 end if

! convert from coo to csr
 call coocsr(n,nnz,c2coo,i2coo,j2coo,rval,cols,rows)
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
