program test_polynomial
   use global
   use polynomial_type
   implicit none

   type(polynomial) :: p
   real(kind=wp) :: a(3), c
   integer :: e(1,3)
   character(len=var_len) :: x(1)

   a = [5, 1, 7.2]
   c = 5
   e = [3, 2, 1]
   x = 'x'

   call p%set_poly(a,e,x,c)
   call p%prnt
end program polynomial
