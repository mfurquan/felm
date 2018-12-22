module element_type
   use global
   use integration
   implicit none

   type,abstract :: element
      type(polynomial),allocatable,private :: shapefn(:), dshapefn(:,:)
      type(integration),private :: elemcub
      contains
         procedure(set_elem),deferred :: set_element
         procedure :: eval, deriv
   end type element

   abstract interface
      pure subroutine set_elem(this,elshape,order)
         class(element),intent(inout) :: this
         character(len=ele_len),intent(in) :: elshape
         shortint,intent(in) :: order
   end interface

   contains

      elemental function eval(this,)
end module element_type
