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
         generic :: set_polynomial => set_poly,  set_const
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
               tmp = x(j)**this%monomial(j,i)
            end do
            eval = eval + tmp
         end do
      end function eval

      pure function deriv(this,tvar)
         class(polynomial),intent(in) :: this
         type(polynomial) :: deriv
         character(len=*),intent(in) :: tvar
         integer :: ivar

         if(this%nvar == 0) then
            this%set_polynomial(0.0)
         else if(ANY(tvar /= this%var(:))) then
            this%set_polynomial(0.0)
         else
            do ivar = 1,this%nvar
               if(tvar == this%var(ivar)) exit
            end do
            this%coeff = this%coeff * this%monomial(ivar,:)
            this%monomial(ivar,:) = this%monomial(ivar,:) - 1
         end if
      end function deriv

      ! removes extra allocation, 
      ! stores polynomial in its simplest form
      subroutine clean(this)
         class(polynomial),intent(in) :: this
         real(kind=wp),allocatable :: ccoeff(:)
         real(kind=wp) :: cconst
         integer,allocatable :: cmono(:,:), drow(:), dcol(:)
         character(len=var_len),allocatable :: cvar(:)
         logical :: extra_var(this%nvar), extra_trm(this%nterms)
         integer :: i, j, ivar

         ccoeff = this%coeff
         cmono = this%monomial
         cvar = this%var
         cconst = this%const
         ! checking for superfluous variables
         extra_var = .FALSE.
         do concurrent (ivar = 1 : this%nvar)
            if(ALL(this%monomial(ivar,:) == 0)) then
               cconst = cconst + this%coeff(ivar)
               extra_var(ivar) = .TRUE.
            end if
         end do

         ! checking for superfluous monomials
         extra_trm = .FALSE.
         where(this%coeff == 0.)
            extra_trm = .TRUE.
         end where

         ! deleting extra rows/columns
         allocate(drow(count(extra_var)))
         allocate(dcol(count(extra_trm)))
         j = 0
         do i = 1,size(extra_trm)
            j = j + 1
            if(extra_trm(i)) dcol(j) = i
         end do
         j = 0
         do i = 1,size(extra_var)
            j = j + 1
            if(extra_var(i)) drow(j) = i
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
         integer :: itrm, i, j

         do itrm = 1, this%nterms
            if(this%coeff(i) > 0. .AND. itrm > 1) &
               write(*,'(a)',advance='no') '+'
            write(*,'(f8.5)',advance='no') this%coeff(i)
            do j = 1, this%nvar
               write(*,'(5a)',advance='no') this%var(j)
               if(this%monomial(j,i)/=1) then
                  write(*,'(ai2)',advance='no') '^',this%monomial(j,i)
               end if
            end do
         end do
         write(*,*)
      end subroutine prnt
end module polynomial_type
