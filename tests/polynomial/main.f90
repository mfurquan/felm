#include "basic_types.h"

program test_polynomial
   use global
   use polynomial_type
   use linpol
   implicit none

   type(polynomial) :: p, dp, q
   double :: a(3), c
   longint :: e(1,3), e1(2,3)
   character(len=var_len) :: x(1), x1(2)

   a = [5.0, 1.0, 7.2]
   c = 5.0
   e(1,:) = [3, 2, 1]
   x = 'x'

   call p%set_polynomial(a,e,x,c)
   call p%prnt
   dp=p%deriv('x')
   call dp%prnt

   e1(1,:) = e(1,:)
   e1(2,:) = [2, 1, 0]
   x1      = ['x', 'y']

   call p%set_polynomial(a,e1,x1,c)
   call p%prnt
   dp=p%deriv('x')
   call dp%prnt
   dp=p%deriv('y')
   call dp%prnt
   dp=dp%deriv('x')
   call dp%prnt

   write(*,*) 'dp(3,2)=',dp%eval([3.0_wp,2.0_wp])

   q = linpol('x',-2.0)
   call q%prnt
   q = q*p
   call q%prnt
end program test_polynomial
