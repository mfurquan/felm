module shapefunc
!==============================================================================
! Shape function types and operations
!------------------------------------------------------------------------------
! M. Furquan | Aug 2016
!==============================================================================
! Primary shape functions implemented
!------------------------------------------------------------------------------
! qtype     :  Description
!
! Lagr1DMN  :  1D Lagrangian of order (10*M+N)
!==============================================================================
   use global
   use polynomial
   implicit none

   type :: shapefn
      type(polynomial),allocatable,private :: F(:,:)
      integer,private :: ndim, npar, nord
      contains
         procedure :: set_shapefn

         procedure :: eval
         procedure :: jac
   end type shapefn

   contains

      subroutine set_shapefn(this,qtype)
         class(shapefn),intent(inout) :: this
         character(len=shap_len) :: qtype

         read(qtype(5)  ,'(I)' ) this%ndim
         read(qtype(7:8),'(I2)') this%nord

         select case(qtype(1:4))
         case('Lagr')
            this%npar = this%ndim
            allocate(F(this%nord,this%nord*this%ndim))
         end select
      end subroutine set_shapefn
end module shapefunc
