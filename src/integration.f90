module integration
!==========================================================================
! Quadrature rules type and integration function
!--------------------------------------------------------------------------
! M. Furquan | Dec 2018
!==========================================================================
   use global
   use quadrature_tables
   implicit none

   type :: qrule
      integer :: nq
      real(kind=rp),allocatable :: pts(:), wts(:)
      contains
         procedure :: init => init_qrule
   end type qrule

   interface operator (.otimes.)
      module procedure tensorprod_qrule
   end interface operator (.otimes.)

   contains

      subroutine init_qrule(this, coor, weights)
         implicit none
         class(qrule),intent(inout) :: this
         real(kind=rp),intent(in) :: coor(:,:), weights(:)

         if(allocated(this%pts)) deallocate(this%pts)
         if(allocated(this%wts)) deallocate(this%wts)
         this%pts = coor
         this%wts = weights
      end subroutine set_qrule

      function tensorprod_qrule(ruleA,ruleB) result (ruleC)
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

         call ruleC%init_qrule(tcor,twts)
      end function tensorprod_qrule

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

         call exc_sim%init_qrule(excor,exwt)
      end function exc_sim

      function gauss_1D(npts)
        implicit none
        integer,intent(in) :: npts
        integer :: loc
        type(cubrule) :: gaus_lin

        loc = npts * (npts-1)/2
        call gaus_lin%init_qrule(RESHAPE(gaus_pts(loc+1:loc+npts),[npts+1], &
                                         gaus_wts(loc+1:loc+npts))
      end function gauss_1D
end module integration
