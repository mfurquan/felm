module shapefunc
   use global
   implicit none

   type :: lagrange
      logical,allocatable,private :: node(:)
      integer,private :: ndim
      contains
         procedure :: set_lagrange

         procedure :: is_node
         procedure :: eval
         procedure :: jac
   end type lagrange
end module shapefunc
