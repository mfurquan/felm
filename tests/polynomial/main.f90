#include "../../include/basic_types.h"

program test_polynomial
   use global
   use polynomial_type
   use spec_poly
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
   write(*,*) '*** Testing single variable polynomial ***'
   write(*,'(A7)',advance='no') 'p(x) = '
   call p%prnt
   dp=p%deriv('x')
   write(*,'(A8)',advance='no') 'dp/dx = '
   call dp%prnt

   e1(1,:) = e(1,:)
   e1(2,:) = [2, 1, 0]
   x1      = ['x', 'y']

   call p%set_polynomial(a,e1,x1,c)
   write(*,*) '*** Testing two variable polynomial *** '
   write(*,'(A9)',advance='no') 'p(x,y) = '
   call p%prnt
   write(*,'(A8)',advance='no') 'dp/dx = '
   dp=p%deriv('x')
   call dp%prnt
   write(*,'(A8)',advance='no') 'dp/dy = '
   dp=p%deriv('y')
   call dp%prnt
   write(*,'(A11)',advance='no') 'd2p/dydx = '
   dp=dp%deriv('x')
   call dp%prnt

   write(*,*) 'd2p/dydx|(3,2) =',dp%eval([3.0_rp,2.0_rp])

   write(*,'(A7)',advance='no') 'q(x) = '
   q = linpol('x',-2.0_rp)
   call q%prnt
   write(*,'(A13)',advance='no') 'q(x)p(x,y) = '
   q = q*p
   call q%prnt
end program test_polynomial
