module polynomial
!==========================================================================
! Polynomial type and operations
!--------------------------------------------------------------------------
! M. Furquan | May 2016
!==========================================================================
   use global
   implicit none

   type :: polynomial
      real(kind=wp),allocatable,private :: coeff(:)
      real(kind=wp),private :: const
      integer,allocatable,private :: monomial(:)
      integer,private :: nvar, nterms
      character(len=var_len),allocatable,private :: var(:)
      contains
         generic :: set_polynomial => set_poly,  set_const
         procedure :: eval
         procedure :: deriv
   end type polynomial

   contains

      subroutine set_poly(this,tcoeff,tmono,tvar,val_const)
         class(polynomial),intent(inout) :: this
         real(kind=wp),intent(in) :: tcoeff(:)
         character(len=*),intent(in):: tvar(:)
         integer,intent(in) :: tmono(:)

         if(size(tmono)/=size(coeff)*size(tvar)) &
            error stop "in set_polynomial: incompatible input data"

         this%nvar = size(tvar)
         this%nterms = size(tcoeff)
         this%const = val_const

         ! set coeff
         if(allocated(this%coeff)) then
            if(size(coeff) /= this%nterms) &
               deallocate(this%coeff)
         end if
         this%coeff = tcoeff

         ! set monomial
         if(allocated(this%monomial)) then
            if(size(monomial) /= this%nvar * this%nterms) &
               deallocate(this%monomial
         end if
         this%monomial = reshape(tmono,[this%nvar * this%nterms]))

         !set var
         if(allocated(this%var)) then
            if(size(var) /= this%nvar) &
               deallocate(this%var)
         end if
         this%var = tvar
      end subroutine set_poly

      subroutine set_const(this,val_const)
         class(polynomial),intent(inout) :: this
         real(kind=wp),intent(in) :: val_const

         this%nvar = 0
         this%nterms = 0
         this%const = val_const

         if(allocated(coeff)) deallocate(coeff)
         if(allocated(monomial)) deallocate(monomial)
         if(allocated(var)) deallocate(var)
      end subroutine set_const

      pure function eval(this,x)
         class(polynomial),intent(in) :: this
         real(kind=wp),intent(in) :: x(:)
         real(kind=wp) :: eval, tmp
         integer :: i, j

         ! simple concurrent evaluation
         ! should perhaps try Horner's method
         eval=0.
         do concurrent (i = 1:this%nterms)
            tmp = 1
            do concurrent (j = 1:this%nvar)
               tmp = x(j)**monomial((i-1)*this%nvar +j)
            end do
            eval = eval + tmp
         end do
      end function eval

      pure function deriv(this,tvar)
         class(polynomial),intent(in) :: this
         class(polynomial),intent(out) :: deriv
         character(len=*),intent(in) :: tvar

         if(ANY(tvar /= this%var(:))) then
            allocate()
      end function deriv
end module polynomial
