module element_type
   use global
   use polynomial
   use integration
   implicit none

   type :: element
      type(polynomial),allocatable,private :: shapefn(:), dshapefn(:,:)
      type(integration),private :: elemcub
      contains
         procedure :: set_element
         procedure :: eval, deriv
   end type element

   contains

      subroutine set_element(this,code)
         class(this),intent(inout) :: this
         character(len=*),intent(in) :: code
         character(len=2) :: family
         character :: dummy
         integer :: order, dimen

         read(code,'(2A2IAI)') family,order,dummy,dimen
         select case(family)
            case('LR')
               if(dimen == '1') then
                  allocate(shapefn(order + 1))
               end if
         end select
      end subroutine set_element
end module element_type
