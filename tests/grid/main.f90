#include "../../include/basic_types.h"

program main
   use global
   use grid_type
   implicit none

   type(grid) :: simple_grid
   longint :: nn, nsd, ne, cnen_max
   double,allocatable :: x(:,:)
   longint,allocatable :: ien(:,:)
   longint :: i, j, nen

   open(10,file='test.nod')
   read(10,*) nn, nsd
   allocate(x(nsd,nn))
   do j = 1, nn
      read(10,*) (x(i,j), i = 1, nsd)
   end do
   close(10)

   open(10,file='test.ele')
   read(10,*) ne, cnen_max
   allocate(ien(cnen_max,ne))
   do j = 1, ne
      read(10,*) nen, (ien(i,j), i = 1, nen)
   end do
   close(10)

   call simple_grid%set_grid(nsd,nn,ien,x)
   do i = 1, ne
      write(*,*) simple_grid%get_vals(i)
   end do
   call simple_grid%set_val(2,3,[2.0_rp, 2.0_rp, 1.0_rp])
   write(*,*) 'after changing a value'
   write(*,*) simple_grid%get_aval(2,3)
end program main
