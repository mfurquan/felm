module utils
   use param
   implicit none

   interface operator (.otimes.)
      module procedure tesnorprod_array
   end interface operator (.otimes.)
contains
   pure function tensorprod_array(A,B) result (C)
      real(kind=rp),intent(in) :: A(:), B(:)
      real(kind=rp) :: C(:,:)

      C = SPREAD(A,2,SIZE(B))*SPREAD(B,1,SIZE(A))
   end function tensorprod_array
end module utils
