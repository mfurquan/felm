module polynomials
   use param
   implicit none

   type :: polynom(deg)
      integer,len :: deg
      real(kind=rp) :: roots(deg), factor
   contains
      procedure :: eval
      procedure :: deriv
   end type polynom

contains
   elemental function eval(this,x)
      type(polynom),intent(in) :: this
      real(kind=rp),intent(in) :: x
      real(kind=rp) :: eval
      
      eval = PRODUCT(x-this%roots)*factor
   end function eval

   elemental function deriv(this,x)
      type(polynom),intent(in) :: this
end module polynomials
