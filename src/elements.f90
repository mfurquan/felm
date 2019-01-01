module element_type
   use param
   use integration
   implicit none

   type :: element
      integer :: nqd, npr, nsh
      real(kind=rp) :: sh(:,:), dsh(:,:,:), wq(:)
      type(polynomial),allocatable,private :: shapefn(:), dshapefn(:,:)
      contains
   end type element
end module element_type
