module integration
!==========================================================================
! Cubature rules type and integration function
!--------------------------------------------------------------------------
! M. Furquan | Apr 2016
!==========================================================================
   use global
   use gauss_cubature
   implicit none

   type :: cubrule
      real(kind=wp),allocatable,private :: cpts(:), cwts(:)
      integer :: ndim
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

      subroutine set_crule(this, coor, wts, dimn)
         implicit none
         class(cubrule),intent(inout) :: this
         real(kind=wp),intent(in) :: coor(:), wts(:)
         integer,intent(in) :: dimn
         integer :: npts

         if(size(coor) /= size(wts)*dimn) &
            error stop "in set_crule: incompatible input data"

         this%ndim = dimn
         npts = size(wts)

         if(.NOT.allocated(this%cpts)) &
            allocate(this%cpts(dimn*npts))

         if(size(this%cpts) /= npts*dimn) then
            deallocate(this%cpts)
            allocate(this%cpts(dimn*npts))
         end if

         this%cpts = coor

         if(.NOT.allocated(this%cwts)) &
            allocate(this%cwts(npts))

         if(size(this%cwts) /= npts) then
            deallocate(this%cwts)
            allocate(this%cwts(npts))
         end if

         this%cwts = wts
      end subroutine set_crule

      pure function get_cpts(this,ind)
         implicit none
         class(cubrule),intent(in) :: this
         integer,intent(in) :: ind
         real(kind=wp) :: get_cpts(this%ndim)

         get_cpts=this%cpts((ind-1)*this%ndim+1: &
                                ind*this%ndim)
      end function get_cpts

      pure function get_cwts(this,ind)
         implicit none
         class(cubrule),intent(in) :: this
         integer,intent(in) :: ind
         real(kind=wp) :: get_cwts

         get_cwts=this%cwts(ind)
      end function get_cwts

      pure function get_ncub(this)
         implicit none
         class(cubrule),intent(in) :: this
         integer :: get_ncub

         get_ncub = size(this%cwts)
      end function get_ncub

      pure function get_ndim(this)
         implicit none
         class(cubrule),intent(in) :: this
         integer :: get_ndim

         get_ndim = this%ndim
      end function get_ndim

      function tensorprod_cubrule(ruleA,ruleB) result (ruleC)
      !==========================================================
      ! Implements tensor product of two cubature rules
      !==========================================================
         implicit none
         type(cubrule),intent(in) :: ruleA, ruleB
         type(cubrule) :: ruleC
         real(kind=wp),allocatable :: tcor(:,:), twts(:)
         integer :: i, j, k
         integer :: ncubA, ncubB, ncubC
         integer :: ndimA, ndimB, ndimC

         ncubA = ruleA%get_ncub()
         ncubB = ruleB%get_ncub()

         ndimA = ruleA%get_ndim()
         ndimB = ruleB%get_ndim()

         ncubC = ncubA * ncubB
         ndimC = ndimA + ndimB

         allocate(tcor(ndimC,ncubC))
         allocate(twts(ncubC))
         k=0
         do i=1,ncubA
            do j=1,ncubB
               k = k + 1
               tcor(      1:ndimA,k) = ruleA%get_cpts(i)
               tcor(ndimA+1:ndimC,k) = ruleB%get_cpts(j)
               twts(k) = ruleA%get_cwts(i) * ruleB%get_cwts(j)
            end do
         end do

         call ruleC%set_crule(reshape(tcor,[size(tcor)]),twts,ndimC)
      end function tensorprod_cubrule

      function exc_sim(nord,sdim)
      !================================================================
      ! Exact cubature rule for integrating cubic polynomial
      ! over a sdim-th dimensional simplex.
      ! Ref: Hammer & Stroud "Numerical integration over simplexes"
      !================================================================
         implicit none
         integer,intent(in) :: nord, sdim
         type(cubrule) :: exc_sim
         real(kind=wp),allocatable :: excor(:,:), exwt(:)
         real(kind=wp),allocatable :: vertex(:,:)
         real(kind=wp) :: vol
         integer :: i, j

         allocate(vertex(sdim,sdim + 1))
         allocate(excor(sdim,sdim + 2))
         allocate(exwt(sdim + 2))

         vertex = 0.
         vol = 1.
         do concurrent (i = 2:sdim+1)
            vertex(i-1,i)=1.
            vol=vol/(i-1)
         end do

         excor(:,sdim+2) = sum(vertex,2)/(sdim+1)
         exwt(sdim+2) = -(sdim+1)**2*vol/(4*(sdim+2))

         do concurrent (i = 1:sdim, j = 1:sdim+1)
            excor(i,j) = ((sdim+1.)/(sdim+3.)) * excor(i,sdim+2)
         end do
         excor(:,1:sdim+1) = excor(:,1:sdim+1) + (2./(sdim+3)) * vertex
         exwt(1:sdim+1)    = (sdim+3)**2*vol/(4*(sdim+1)*(sdim+2))

         call exc_sim%set_crule(reshape(excor,[size(excor)]),exwt,sdim)
      end function exc_sim

      function gaus_lin(npts)
        implicit none
        integer,intent(in) :: npts
        integer :: loc
        type(cubrule) :: gaus_lin

        loc = npts * (npts-1)/2
        call gaus_lin%set_crule(gaus_pts(loc+1:loc+npts), &
                                gaus_wts(loc+1:loc+npts),1)
      end function gaus_lin

      pure function integrate(integrand,crule)
         implicit none
         real(kind=wp) :: integrate
         type(cubrule),intent(in) :: crule
         integer :: i, ncub
         interface
            pure function integrand(xi)
               use global
               real(kind=wp),intent(in) :: xi(:)
               real(kind=wp) :: integrand
            end function integrand
         end interface

         ncub = crule%get_ncub()
         integrate = 0.
         do concurrent (i = 1:ncub)
            integrate = integrate &
                      + crule%get_cwts(i)*integrand(crule%get_cpts(i))
         end do
      end function integrate
end module integration
