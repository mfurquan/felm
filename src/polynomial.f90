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
      real(kind=wp),allocatable,private :: coeff(:)
      real(kind=wp),private :: const
      integer,allocatable,private :: monomial(:,:)
      integer,private :: nvar, nterms
      character(len=var_len),allocatable,private :: var(:)
      contains
         procedure,private :: set_poly
         procedure,private :: set_const
         procedure,private :: assign_poly
         generic :: set_polynomial => set_poly,  set_const
         generic :: assignment(=) => assign_poly
         procedure :: clean
         procedure :: eval
         procedure :: deriv

         procedure :: prnt
   end type polynomial

   contains

      subroutine set_poly(this,tcoeff,tmono,tvar,val_const)
         class(polynomial),intent(inout) :: this
         real(kind=wp),intent(in) :: tcoeff(:), val_const
         character(len=*),intent(in):: tvar(:)
         integer,intent(in) :: tmono(:,:)

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

      subroutine set_const(this,val_const)
         class(polynomial),intent(inout) :: this
         real(kind=wp),intent(in) :: val_const

         this%nvar = 0
         this%nterms = 0
         this%const = val_const

         if(allocated(this%coeff)) deallocate(this%coeff)
         if(allocated(this%monomial)) deallocate(this%monomial)
         if(allocated(this%var)) deallocate(this%var)
      end subroutine set_const

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
         real(kind=wp),intent(in) :: x(:)
         real(kind=wp) :: eval, tmp
         integer :: i, j

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
         integer :: ivar

         if(this%nvar == 0) then
            call deriv%set_polynomial(0.0_wp)
         else if(ALL(tvar /= this%var(:))) then
            call deriv%set_polynomial(0.0_wp)
         else
            do ivar = 1,this%nvar
               if(tvar == this%var(ivar)) exit
            end do
            deriv = this
            deriv%const = 0.0_wp
            deriv%coeff = this%coeff * this%monomial(ivar,:)
            deriv%monomial(ivar,:) = this%monomial(ivar,:) - 1
         end if
         call deriv%clean
      end function deriv

      ! removes extra allocation, 
      ! stores polynomial in its simplest form
      subroutine clean(this)
         class(polynomial),intent(inout) :: this
         real(kind=wp),allocatable :: ccoeff(:)
         real(kind=wp) :: cconst
         integer,allocatable :: cmono(:,:), drow(:), dcol(:)
         character(len=var_len),allocatable :: cvar(:)
         logical :: extra_var(this%nvar), extra_trm(this%nterms)
         integer :: i, j, ivar, imono

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
         integer :: itrm, ivar
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
end module polynomial_type
