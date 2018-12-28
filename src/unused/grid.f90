#include "../include/basic_types.h"

module grid_type
!===========================================================================
! defined type for 'grid of nodal values'
!---------------------------------------------------------------------------
! M. Furquan | Oct 2016
!===========================================================================
   use global
   implicit none

   type :: elenods
      double :: val(nsd_max,nen_max)
   end type elenods

   type :: grid
      longint,            private :: ne, nn, ndim
      longint,allocatable,private :: ien(:), nen(:)
      double, allocatable,private :: val(:)
      contains
         procedure :: set_grid, set_val
         procedure :: get_nen,  get_aval, get_vals
   end type grid

   contains

!***************************** grid-bound procedures ***********************
   pure subroutine set_grid(this,numdof,numnp,cells,nodes)
      class(grid),    intent(inout) :: this
      double,optional,intent(in)    :: nodes(:,:)
      longint,        intent(in)    :: cells(:,:), numnp, numdof
      longint                       :: ncn, lnen_max, m, i, j

      if(present(nodes)) then
         if(numdof /= size(nodes,1)) &
            error stop &
            "in set_grid: incompatible input data 'nodes' and 'numdof'"
         if(numnp /= size(nodes,2)) &
            error stop &
            "in set_grid: incompatible input data 'nodes' and 'numnp'"
      end if
      
      this%ndim = numdof
      this%nn   = numnp
      this%ne   = size(cells,2)
      ncn       = count(cells/=0)
      lnen_max  = size(cells,1)

      allocate(this%val(this%ndim*this%nn))
      allocate(this%ien(ncn))
      allocate(this%nen(this%ne + 1))

      if(present(nodes)) &
         this%val = reshape(nodes,[this%ndim*this%nn])

      m = 0
      do i = 1, this%ne
         do j = 1, lnen_max
            if(cells(j,i) /= 0) then
               m = m + 1
               this%ien(m) = (cells(j,i) - 1)*this%ndim + 1
               if(j == 1) this%nen(i) = m
            end if
         end do
      end do
      this%nen(this%ne + 1) = ncn + 1
   end subroutine set_grid

   elemental function get_nen(this,iele)
      class(grid),intent(in)  :: this
      shortint                :: get_nen
      longint,    intent(in)  :: iele

      get_nen = this%nen(iele + 1) - this%nen(iele)
   end function get_nen

   pure function get_aval(this,iele,inod)
      class(grid),intent(in)  :: this
      longint,    intent(in)  :: iele, inod
      double                  :: get_aval(this%ndim)
      longint                 :: n

      n = this%ien(this%nen(iele) + inod - 1)
      get_aval = this%val(n : n + this%ndim - 1)
   end function get_aval

   pure function get_vals(this,iele)
      class(grid),  intent(in)  :: this
      longint,      intent(in)  :: iele
      type(elenods)             :: get_vals
      longint                   :: inod, isd, m, n, lnen

      lnen = this%get_nen(iele)
      m = this%nen(iele)
      do concurrent (inod = 1:lnen)
         n = this%ien(m + inod - 1)
         do concurrent (isd = 1:this%ndim)
            get_vals%val(isd,inod) = this%val(n + isd - 1)
         end do
      end do
      get_vals%val(:,lnen + 1:) = 0.0_rp
   end function get_vals

   pure subroutine set_val(this,iele,inod,nodvalue)
      class(grid),intent(inout) :: this
      longint,    intent(in)    :: iele, inod
      double,     intent(in)    :: nodvalue(:)
      longint                   :: n

      n = this%nen(iele) + inod - 1
      n = this%ien(n)
      this%val(n : n + this%ndim - 1) = nodvalue
   end subroutine set_val

end module grid_type
