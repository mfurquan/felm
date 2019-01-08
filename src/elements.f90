module element_type
   use param
   use utils
   implicit none

   type :: element
      integer :: nqd, npr, nsh
      real(kind=rp) :: sh(:,:), dsh(:,:,:), wq(:)
   end type element

contains
   pure function tensorprod_element(A,B) result(C)
      type(element),intent(in) :: A, B
      type(element) :: C

      C%nqd = A%nqd*B%nqd
      C%npr = A%npr+B%npr
      C%nsh = A%nsh*B%nsh

      C% sh = RESHAPE(A%sh.otimes.B%sh,[C%nqd,C%nsh])
      do concurrent (i = 1:A%npr)
         C%dsh(:,:,i) = A%dsh(:,:,i).otimes.B%sh
      end do
      do concurrent (i = 1:B%npr)
         C%dsh(:,:,A%npr+i) = A%sh.otimes.B%dsh(:,:,i)
      end do

      C%wq = A%wq.otimes.B%wq
   end function tensorprod_element
end module element_type
