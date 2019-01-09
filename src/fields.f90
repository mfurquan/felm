module field_MOD
   use param
   implicit none

   type :: field
      real(kind=rp),allocatable :: d(:,:)
   end type field
end field_MOD
