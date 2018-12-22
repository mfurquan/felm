#include "../include/basic_types.h"

module array_utils
!========================================================================
! Some utilities for handling dynamic arrays
!------------------------------------------------------------------------
! M. Furquan | Jun 2016
!========================================================================
   use global
   implicit none

   interface del_elem
      module procedure del_elemf, del_elemi, del_elema
   end interface del_elem

   interface del_rows
      module procedure del_rowf, del_rowi
   end interface del_rows

   interface del_cols
      module procedure del_colf, del_coli
   end interface del_cols

   contains

      pure subroutine del_elemf(array,elem)
         double,allocatable,intent(inout) :: array(:)
         double,allocatable :: tmp(:)
         longint,intent(in) :: elem(:)
         longint :: i, j, m, n

         m = size(elem)
         if(m > 0) then
            n = size(array)
            allocate(tmp(n - m))
            j = 0
            do i = 1, n
               if(ANY(elem /= i)) then
                  j = j + 1
                  tmp(j) = array(i)
               end if
            end do
            call move_alloc(tmp,array)
         end if
      end subroutine del_elemf

      pure subroutine del_elemi(array,elem)
         longint,allocatable,intent(inout) :: array(:)
         longint,allocatable :: tmp(:)
         longint,intent(in) :: elem(:)
         longint :: i, j, m, n

         m = size(elem)
         if(m > 0) then
            n = size(array)
            allocate(tmp(n - m))
            j = 0
            do i = 1, n
               if(ANY(elem /= i)) then
                  j = j + 1
                  tmp(j) = array(i)
               end if
            end do
            call move_alloc(tmp,array)
         end if
      end subroutine del_elemi

      pure subroutine del_elema(array,elem)
         character(len=*),allocatable,intent(inout) :: array(:)
         character(len=LEN(array(1))),allocatable :: tmp(:)
         longint,intent(in) :: elem(:)
         longint :: i, j, m, n

         m = size(elem)
         if(m > 0) then
            n = size(array)
            allocate(tmp(n - m))
            j = 0
            do i = 1, n
               if(ANY(elem /= i)) then
                  j = j + 1
                  tmp(j) = array(i)
               end if
            end do
            call move_alloc(tmp,array)
         end if
     end subroutine del_elema

     pure subroutine del_rowf(array,elem)
         double,allocatable,intent(inout) :: array(:,:)
         double,allocatable :: tmp(:,:)
         longint,intent(in) :: elem(:)
         longint :: i, j, m, n

         m = size(elem)
         if(m > 0) then
            n = size(array,1)
            allocate(tmp(n - m,size(array,2)))
            j = 0
            do i = 1, n
               if(ANY(elem /= i)) then
                  j = j + 1
                  tmp(j,:) = array(i,:)
               end if
            end do
            call move_alloc(tmp,array)
         end if
      end subroutine del_rowf

      pure subroutine del_colf(array,elem)
         double,allocatable,intent(inout) :: array(:,:)
         double,allocatable :: tmp(:,:)
         longint,intent(in) :: elem(:)
         longint :: i, j, m, n

         m = size(elem)
         if (m > 0) then
            n = size(array,2)
            allocate(tmp(size(array,1), n - m))
            j = 0
            do i = 1, n
               if(ANY(elem /= i)) then
                  j = j + 1
                  tmp(:,j) = array(:,i)
               end if
            end do
            call move_alloc(tmp,array)
         end if
      end subroutine del_colf

      pure subroutine del_rowi(array,elem)
         longint,allocatable,intent(inout) :: array(:,:)
         longint,allocatable :: tmp(:,:)
         longint,intent(in) :: elem(:)
         longint :: i, j, m, n

         m = size(elem)
         if(m > 0) then
            n = size(array,1)
            allocate(tmp(n - m,size(array,2)))
            j = 0
            do i = 1, n
               if(ANY(elem /= i)) then
                  j = j + 1
                  tmp(j,:) = array(i,:)
               end if
            end do
            call move_alloc(tmp,array)
         end if
      end subroutine del_rowi

      pure subroutine del_coli(array,elem)
         longint,allocatable,intent(inout) :: array(:,:)
         longint,allocatable :: tmp(:,:)
         longint,intent(in) :: elem(:)
         longint :: i, j, m, n

         m = size(elem)
         if(m > 0) then
            n = size(array,2)
            allocate(tmp(size(array,1), n - m))
            j = 0
            do i = 1, n
               if(ANY(elem /= i)) then
                  j = j + 1
                  tmp(:,j) = array(:,i)
               end if
            end do
            call move_alloc(tmp,array)
         end if
      end subroutine del_coli
end module array_utils
