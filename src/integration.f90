module integration
!==========================================================================
! Cubature rules types and integration function
!--------------------------------------------------------------------------
! M. Furquan | Apr 2016
!==========================================================================
   use global
   use gauss_cubature
   implicit none

   type :: cubrule
      real(kind=wp),allocatable,private :: cpts(:,:), cwts(:)
      contains
         procedure :: set_crule

         procedure :: get_cpts
         procedure :: get_cwts
         procedure :: get_ncub
         procedure :: get_ndim
   end type cubrule

   interface cubrule
      module procedure exc_sim, gaus_lin
   end interface cubrule

   interface operator (.otimes.)
      module procedure tensorprod_cubrule
   end interface operator (.otimes.)

   contains

      subroutine set_crule(this,coor,wts)
         implicit none
         class(cubrule),intent(inout) :: this
         real(kind=wp),intent(in) :: coor(:,:),wts(:)

         if(size(coor,2)/=size(wts)) &
            error stop "ERROR in set_crule: incompatible input data"

         if(.NOT.allocated(this%cpts)) &
            allocate(this%cpts,mold=coor)

         if(size(this%cpts,2)/=size(coor,2)) then
            deallocate(this%cpts)
            allocate(this%cpts,mold=coor)
         end if

         this%cpts=coor

         if(.NOT.allocated(this%cwts)) &
            allocate(this%cwts,mold=wts)

         if(size(this%cwts)/=size(wts)) then
            deallocate(this%cwts)
            allocate(this%cwts,mold=wts)
         end if

         this%cwts=wts
      end subroutine set_crule

      elemental function get_cpts(this,ind)
         implicit none
         class(cubrule),intent(in) :: this
         integer,intent(in) :: ind
         real(kind=wp),intent(out) :: get_cpts(:)

         get_cpts=this%cpts(:,ind)
      end function get_cpts

      elemental function get_cwts(this,ind)
         implicit none
         class(cubrule),intent(in) :: this
         integer,intent(in) :: ind
         real(kind=wp),intent(out) :: get_cwts

         get_cwts=this%cwts(ind)
      end function get_cwts

      elemental function get_ncub(this)
         implicit none
         class(cubrule),intent(in) :: this
         integer,intent(out) :: get_ncub

         get_ncub=size(this%cwts)
      end function get_ncub

      elemental function get_ndim(this)
         implicit none
         class(cubrule),intent(in) :: this
         integer,intent(out) :: get_ndim

         get_ndim=size(1,cpts)
      end function get_ndim

      pure function tensorprod_cubrule(ruleA,ruleB) result (ruleC)
      !==========================================================
      ! Implements tensor product of two cubature rules
      !==========================================================
         implicit none
         class(cubrule),intent(in) :: ruleA,ruleB
         class(cubrule),intent(out) :: ruleC
         real(kind=wp),allocatable :: tcor(:,:),twts(:)
         integer :: i,j,k
         integer :: ncubA,ncubB,ncubC
         integer :: ndimA,ndimB,ndimC
         integer,allocatable :: ind(:)

         ncubA=ruleA%get_ncub()
         ncubB=ruleB%get_ncub()

         ndimA=ruleA%get_ndim()
         ndimB=ruleB%get_ndim()

         ncubC=ncubA*ncubB
         ndimC=ndimA+ndimB

         allocate(ind(max(ncubC,ndimC)))
         ind=[(i,i=1,size(ind))]

         allocate(tcor(ndimC,ncubC))
         allocate(twts(ncubC))

         k=0
         do i=1,ncubA
            do j=1,ncubB
               k=k+1
               tcor(1:ndimA,k)      =ruleA%get_cpts(i)
               tcor(ndimA+1:ndimB,k)=ruleB%get_cpts(j)
               twts(k)=ruleA%get_cwts(i)*ruleB%get_cwts(j)
            end do
         end do

         call ruleC%set_crule(tcor,twts)
      end function tensorprod_cubrule

      pure function exc_sim(nord,sdim)
      !================================================================
      ! Exact cubature rule for integrating cubic polynomial
      ! over a sdim-th dimensional simplex.
      ! Ref: Hammer & Stroud "Numerical integration over simplexes"
      !================================================================
         implicit none
         integer,intent(in) :: nord,sdim
         class(cubrule),intent(out) :: exc_sim
         real(kind=wp),allocatable :: excor(:,:),exwt(:)
         real(kind=wp),allocatable :: vertex(:,:)
         real(kind=wp) :: vol

         allocate(vertex(sdim,sdim+1))
         allocate(excor(sdim,sdim+2))
         allocate(exwt(sdim+2))

         vertex=0.
         vol=1.
         do i=2,sdim+1
            vertex(i-1,i)=1.
            vol=vol/(i-1)
         end do

         excor(:,sdim+2)=sum(vertex)/(sdim+1)
         exwt(sdim+2)=-(sdim+1)**2*vol/(4*(sdim+2))

         excor(:,1:sdim+1)=(2./(sdim+3))*vertex &
                          +((sdim+1)/(sdim+3))*excor(:,sdim+2)

         do concurrent (i=1:sdim+1)
            exwt(i)=(i+2)**2*vol/(4*i*(i+1))
         end do

         call exc_sim%set_crule(excor,exwt)
      end function exc_sim

      pure function gaus_lin(npts)
        implicit none
        integer,intent(in) :: npts
        integer :: loc
        class(cubrule),intent(out) :: gaus_lin

        loc=npts*(npts-1)/2
        call gaus_lin%set_crule(gaus_pts(loc:loc+npts),gaus_wts(loc:loc+npts))
      end function gaus_lin

      pure function integrate(integrand,crule)
         implicit none
         real(kind=wp),intent(out) :: integrate
         class(cubrule),intent(in) :: crule
         integer,allocatable :: ind(:)
         integer :: i
         interface
            elemental function integrand(xi)
               real(kind=wp),intent(in) :: xi(:)
               real(kind=wp),intent(out) :: integrand
            end function integrand
         end interface

         allocate(ind(crule%get_ncub()))
         ind=[(i,i=1,size(ind))]

         integrate=sum(crule%get_cwts(ind)*integrand(crule%get_cpts(ind)))
      end function integrate
end module integration
