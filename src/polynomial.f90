#include "../include/basic_types.h"

module polynomial_type
!==========================================================================
! Polynomial type and operations
!--------------------------------------------------------------------------
! M. Furquan | May 2016
!==========================================================================
   use global
   use array_utils
   implicit none

   type :: polynomial
      double,allocatable,private :: coeff(:)
      double,private :: const
      longint,allocatable,private :: monomial(:,:)
      longint,private :: nvar, nterms
      character(len=var_len),allocatable,private :: var(:)
      contains
         procedure,private :: set_poly
         procedure,private :: set_const
         procedure,private :: assign_poly
         procedure,private :: multiply_poly
         generic :: set_polynomial => set_poly,  set_const
         generic :: assignment(=) => assign_poly
         generic :: operator(*) => multiply_poly
         procedure :: get_coeff, get_monomial, get_const
         procedure :: get_nvar, get_nterms, get_var
         procedure :: clean
         procedure :: eval
         procedure :: deriv

         procedure :: prnt
   end type polynomial

   contains

      pure subroutine set_poly(this,tcoeff,tmono,tvar,val_const)
         class(polynomial),intent(inout) :: this
         double,intent(in) :: tcoeff(:), val_const
         character(len=*),intent(in):: tvar(:)
         longint,intent(in) :: tmono(:,:)

         if(size(tmono)/=size(tcoeff)*size(tvar)) &
            error stop "in set_polynomial: incompatible input data"

         this%nvar = size(tvar)
         this%nterms = size(tcoeff)
         this%const = val_const

         ! set coeff
         if(allocated(this%coeff)) then
            if(size(this%coeff) /= this%nterms) &
               deallocate(this%coeff)
         end if
         this%coeff = tcoeff

         ! set monomials
         if(allocated(this%monomial)) then
            if(ALL(shape(this%monomial) /= shape(tmono))) &
               deallocate(this%monomial)
         end if
         this%monomial = tmono

         !set var
         if(allocated(this%var)) then
            if(size(this%var) /= this%nvar) &
               deallocate(this%var)
         end if
         this%var = tvar
      end subroutine set_poly

      pure subroutine set_const(this,val_const)
         class(polynomial),intent(inout) :: this
         double,intent(in) :: val_const

         this%nvar = 0
         this%nterms = 0
         this%const = val_const

         if(allocated(this%coeff))    deallocate(this%coeff)
         if(allocated(this%monomial)) deallocate(this%monomial)
         if(allocated(this%var))      deallocate(this%var)
      end subroutine set_const

      pure function get_coeff(this)
         class(polynomial),intent(in) :: this
         double :: get_coeff(this%nterms)

         get_coeff = this%coeff
      end function get_coeff

      pure function get_monomial(this)
         class(polynomial),intent(in) :: this
         longint :: get_monomial(this%nvar,this%nterms)

         get_monomial = this%monomial
      end function get_monomial

      elemental function get_var(this,ivar)
         ! procedure made elemental to avoid gfortran internal error
         ! should return whole var array, otherwise
         class(polynomial),intent(in) :: this
         shortint,intent(in) :: ivar
         character(len=var_len) :: get_var

         get_var = this%var(ivar)
      end function get_var

      pure function get_nvar(this)
         class(polynomial),intent(in) :: this
         longint :: get_nvar

         get_nvar = this%nvar
      end function get_nvar

      pure function get_nterms(this)
         class(polynomial),intent(in) :: this
         longint :: get_nterms

         get_nterms = this%nterms
      end function get_nterms

      pure function get_const(this)
         class(polynomial),intent(in) :: this
         double :: get_const

         get_const = this%const
      end function get_const

      subroutine assign_poly(lvalue, rvalue)
         class(polynomial),intent(in) :: rvalue
         class(polynomial),intent(out) :: lvalue

         if(rvalue%nterms>0) then
            call lvalue%set_polynomial(rvalue%coeff,rvalue%monomial, &
               rvalue%var,rvalue%const)
         else
            call lvalue%set_polynomial(rvalue%const)
         end if
      end subroutine assign_poly

      pure function eval(this,x)
         class(polynomial),intent(in) :: this
         double,intent(in) :: x(:)
         double :: eval, tmp
         longint :: i, j

         ! simple concurrent evaluation
         ! should perhaps try Horner's method
         eval=0.
         do concurrent (i = 1:this%nterms)
            tmp = 1
            do concurrent (j = 1:this%nvar)
               tmp = tmp * x(j)**this%monomial(j,i)
            end do
            eval = eval + tmp * this%coeff(i)
         end do
      end function eval

      function deriv(this,tvar)
         class(polynomial),intent(in) :: this
         type(polynomial) :: deriv
         character(len=*),intent(in) :: tvar
         longint :: ivar

         if(this%nvar == 0) then
            call deriv%set_polynomial(0.0_rp)
         else if(ALL(tvar /= this%var(:))) then
            call deriv%set_polynomial(0.0_rp)
         else
            do ivar = 1,this%nvar
               if(tvar == this%var(ivar)) exit
            end do
            deriv = this
            deriv%const = 0.0_rp
            deriv%coeff = this%coeff * this%monomial(ivar,:)
            deriv%monomial(ivar,:) = this%monomial(ivar,:) - 1
         end if
         call deriv%clean
      end function deriv

      ! removes extra allocation, 
      ! stores polynomial in its simplest form
      subroutine clean(this)
         class(polynomial),intent(inout) :: this
         double,allocatable :: ccoeff(:)
         double :: cconst
         longint,allocatable :: cmono(:,:), drow(:), dcol(:)
         character(len=var_len),allocatable :: cvar(:)
         logical :: extra_var(this%nvar), extra_trm(this%nterms)
         longint :: i, j, ivar, imono

         ccoeff = this%coeff
         cmono = this%monomial
         cvar = this%var
         cconst = this%const
         ! checking for superfluous variables
         extra_var = .FALSE.
         do concurrent (ivar = 1 : this%nvar)
            if(ALL(this%monomial(ivar,:) == 0)) &
               extra_var(ivar) = .TRUE.
         end do

         ! checking for superfluous monomials
         extra_trm = .FALSE.
         where(this%coeff == 0.)
            extra_trm = .TRUE.
         end where

         do concurrent (imono = 1 : this%nterms)
            if(ALL(this%monomial(:,imono)==0)) then
               cconst = cconst + this%coeff(imono)
               extra_trm(imono) = .TRUE.
            end if
         end do

         ! deleting extra rows/columns
         allocate(drow(count(extra_var)))
         allocate(dcol(count(extra_trm)))
         j = 0
         do i = 1,size(extra_trm)
            if(extra_trm(i)) then
               j = j + 1
               dcol(j) = i
            end if
         end do
         j = 0
         do i = 1,size(extra_var)
            if(extra_var(i)) then
               j = j + 1
               drow(j) = i
            end if
         end do

         call del_elem(ccoeff,dcol)
         call del_elem(cvar,  drow)
         call del_cols(cmono, dcol)
         call del_rows(cmono, drow)

         deallocate(drow); deallocate(dcol)

         call this%set_polynomial(ccoeff,cmono,cvar,cconst)
      end subroutine clean

      subroutine prnt(this)
         class(polynomial) :: this
         longint :: itrm, ivar
         character(len=4) :: afmt

         do itrm = 1, this%nterms
            if(this%coeff(itrm) > 0. .AND. itrm > 1) &
               write(*,'(a2)',advance='no') ' +'
            write(*,'(f0.5a)',advance='no') this%coeff(itrm),' '
            do ivar = 1, this%nvar
               write(afmt,'(a2i0a)') '(a', len_trim(this%var(ivar)), ')'
               if(this%monomial(ivar,itrm)/=0) then
                  write(*,afmt,advance='no') this%var(ivar)
                  if(this%monomial(ivar,itrm)/=1) &
                     write(*,'(ai0)',advance='no') '^',this%monomial(ivar,itrm)
               end if
            end do
         end do
         if(this%const/=0) then
            if(this%const>0) write(*,'(af0.5)',advance='no') '+',this%const
            if(this%const<0) write(*,'(f0.5)',advance='no') this%const
         end if
         write(*,*)
      end subroutine prnt

      pure function multiply_poly(poly1,poly2) result (poly3)
         class(polynomial),intent(in) :: poly1, poly2
         type(polynomial) :: poly3
         double,allocatable :: coeff2(:), coeff3(:)
         double :: const2, const3
         longint,allocatable :: mono2(:,:), mono3(:,:)
         character(len=var_len),allocatable :: var2(:), var3(:)
         shortint :: nterms2, nvar2, nterms3, nvar3
         shortint :: comvar(poly1%nvar)
         shortint :: i, j, k, l, m, n

         coeff2  = poly2%get_coeff()
         const2  = poly2%get_const()
         mono2   = poly2%get_monomial()
         var2    = poly2%get_var([(i, i = 1, poly2%get_nvar())])
         nterms2 = poly2%get_nterms()
         nvar2   = poly2%get_nvar()

         comvar  = common_var(poly1%var, var2)
         nterms3 = poly1%nterms * nterms2
         nvar3   = poly1%nvar + nvar2 - count(comvar /= 0)
         const3  = poly1%const * const2

         allocate(coeff3(nterms3))
         allocate(mono3(nvar3,nterms3))
         allocate(var3(nvar3))

         n = 0
         do i = 1, poly1%nvar
            if(comvar(i) == 0) then
               n = n + 1
               var3(n) = poly1%var(i)
            end if
         end do
         var3(n + 1:nvar2) = var2
         
         mono3 = 0
         do concurrent (i = 1:poly1%nterms, j = 1:nterms2)
            k = (i - 1)*nterms2 + j
            m = 0
            do l = 1, poly1%nvar
               if(comvar(l) == 0) then
                  m = m + 1
                  mono3(m,k) = poly1%monomial(l,i)
               else
                  mono3(n + comvar(l),k) = poly1%monomial(l,i)
               end if
            end do
            mono3(n + 1:nvar2,k) = mono3(n + 1:nvar2,k) + mono2(:,j)
            coeff3(k) = poly1%coeff(i)*coeff2(j)
         end do

         call poly3%set_polynomial(coeff3,mono3,var3,const3)

         contains
            pure function common_var(var_str1,var_str2)
               character(len=*),intent(in) :: var_str1(:), var_str2(:)
               shortint :: common_var(poly1%nvar)
               shortint :: i, j, m, n

               m = poly1%nvar; n = nvar2
               common_var = 0
               do concurrent (i = 1:m, j = 1:n)
                  if(var_str1(i) == var_str2(j)) common_var(i) = j
               end do
            end function common_var
      end function multiply_poly
end module polynomial_type
