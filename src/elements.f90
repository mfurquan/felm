module element_type
   use param
   use integration
   implicit none

   type :: element
      integer :: nqd, npr, nsh
      real(kind=rp) :: sh(:,:), dsh(:,:,:), wq(:)
   end type element
end module element_type
