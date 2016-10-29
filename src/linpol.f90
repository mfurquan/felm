#include "basic_types.h"

module linpol
   use global
   use polynomial_type
   implicit none
   contains
      elemental function linpol(var,const)
         character(len=*),intent(in) :: var
         longint,intent(in) :: const
         type(polynomial) :: linpol
         double :: coeff
         longint :: mono

         coeff = 1.0
         mono  = 1
         call linpol%(coeff,mono,var,const)
      end function linpol
end module linpol
