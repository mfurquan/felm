#include "../include/basic_types.h"

module spec_poly
   use global
   use polynomial_type
   implicit none
   contains
      elemental function linpol(var,const)
         character(len=*),intent(in) :: var
         double,intent(in) :: const
         type(polynomial) :: linpol
         character(len=var_len) :: varay(1)
         double :: coeff(1)
         longint :: mono(1,1)

         varay = var
         coeff = 1.0
         mono  = 1
         call linpol%set_polynomial(coeff,mono,varay,const)
      end function linpol
end module spec_poly
