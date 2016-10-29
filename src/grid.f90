#include "basic_types.h"

module grid
!===========================================================================
! defined type for 'grid of nodal values'
!---------------------------------------------------------------------------
! M. Furquan | Oct 2016
!===========================================================================
   use global
   implicit none

   type :: grid
      longint,                   protected :: ne, nn, ndim
      longint,allocatable,         private :: ien(:), nen(:)
      double,   allocatable,target,protected :: val(:)
      contains
         procedure :: set_grid, set_val
         procedure :: get_nen,  get_val
   end type grid

   contains

!***************************** grid-bound procedures ***********************
   subroutine set_grid(this,numdof,numnp,cells,nodes)
      class(grid),intent(inout) :: this
      double,       intent(in)    :: nodes(:,:)
      longint,    intent(in)    :: cells(:,:), numnp, numdof
      longint                   :: ncn, nen_max, m, i, j

      if(present(nodes)) then
         if(numdof /= size(nodes,1)) &
            error stop "in set_grid: incompatible input data 'nodes' and 'numdof'"
         if(numnp /= size(nodes,2)) &
            error stop "in set_grid: incompatible input data 'nodes' and 'numnp'"
      end if
      
      this%ndim = numdof
      this%nn   = numnp
      this%ne   = size(cells,2)
      ncn       = count(cells/=0)
      nen_max   = size(cells,1)

      allocate(this%val(this%ndim*this%nn))
      allocate(this%ien(ncn))
      allocate(this%nen(this%ne + 1))

      if(present(nodes)) &
         val = reshape(nodes,[1, this%ndim*this%nn])

      m = 0
      do i = 1, this%ne
         do j = 1, this%nen_max
            if(cells(j,i) /= 0) then
               m = m + 1
               this%ien(m) = (cells(j,i) - 1)*this%ndim + 1
               if(i == 1) this%nen(i) = m
            end if
         end do
      end do
      this%nen(this%ne + 1) = ncn + 1
   end subroutine set_grid

   elemental function get_nen(this,iele)
      class(this),intent(in) :: this
      shortint               :: get_nen
      longint,intent(in)     :: iele

      get_nen = this%nen(iele + 1) - this%nen(iele)
   end function get_nen

   pure function get_val(this,iele,inod)
      class(this),intent(in) :: this
      longint,    intent(in) :: iele, inod
      double,pointer           :: get_val
      longint                :: n

      n = this%nen(iele + inod - 1)
      get_val => this%val(n : n + this%ndim -1)
   end function get_val

   function set_val(this,iele,inod,nodvalue)
      class(this),intent(in) :: this
      longint,    intent(in) :: iele, inod
      double                   :: nodvalue(:)
      longint                :: n

      n = this%nen(iele + inod - 1)
      this%val(n : n + this%ndim -1) = nodvalue
   end function set_value
