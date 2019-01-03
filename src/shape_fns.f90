module shape_fns
   use param
   use elements
   implicit none
contains
   ! tensor product of shape_fns

   ! 1D Lagrange polynomial
   pure function lagrng_poly(nnod,qr) result (L)
      integer,    intent(in) :: nnod
      type(qrule),intent(in) :: qr
      type(element) :: L
      integer       :: ia(nnod), i
      real(kind=rp) :: xia(nnod)
      
      xia  = [(-1._rp+2._rp*i/(nnod-1), i = 0,nnod-1)]
      ia   = [(i, i = 1,nnod)]

      L%wq = qr%wts
      L%sh = lagr_poly(qr%pts,nnod)
      L%dsh = lagr_poly_deriv(qr%pts,nnod)

      allocate(L% sh(nq,nnod))
      allocate(L%dsh(nq,1,nnod))
      do concurrent (i = 1:nnod)
         L% sh(:,i)   = lagr   (i)
         L%dsh(:,1,i) = lagr_xi(i)
      end do
   contains
      pure function lagr(xq,nnode)
         integer,intent(in)        :: ind
         type(qscalar) :: lagr
         logical                   :: am(nnod), bm(nnod,nnod-1)
         real(kind=rp)             :: xdiff(nnod), divsor
         integer                   :: i


         call lagr%init(SIZE(xq),1,1)
         am = ia /= ind
         bm = RESHAPE([(am.AND.CSHIFT(am,i), i = 1,nnod-1)],[nnod,nnod-1])

         do concurrent (i = 1:lagr%nqd)
            xdiff          = xq(i) - xia
            divsor         = PRODUCT(xia(ind) - xia, am)
            lagr%v(i,1)    = PRODUCT(xdiff,am)/divsor
            lagr%dv(i,1,1) = SUM(PRODUCT(SPREAD(xdiff,2,nnod-1),1,bm))   &
                           / divsor
         end do
      end function lagr
   end function lagrng_poly

   ! Simplex shape functions
   pure function simplex(this,norder,nq,nsd)
      integer,intent(in) :: norder, ninteg, nsd
      type(qvar) :: simplex
      type(tuple(nsd+1)),allocatable :: ia(:)
      type(qvar) ,allocatable :: T(nq,1,1)
      type(quad_rule) :: qr
      integer :: i, j, isd

      ! generate indices for nodes of the simplex:
      ia = AP(1,1,norder+1).raised.(nsd+1)
      ia = PACK(ia,tsum(ia)==nsd+norder+1)

      ! import quadrature rule:
      call qr%simplex_quadrature(nsd,nq)
      this%wt = qr%wt

      ! calculate T factors and derivatives:
      call calc_T   (T,   qr%xq,norder+1,nq)
      call calc_T_xi(T_xi,qr%xq,norder+1,nq)

      ! initialize shape_fn:
      this%nsd = nsd
      this%nqd = nq
      this%nen = SIZE(indices,2)
      this%f   = combine_T(ia)
      do concurrent (isd = 1:nsd)
         this%Df(:,isd) = combine_T_xi(ia)
      end do
   contains
      elemental function tsum(a)
         type(tuple(nsd+1)),intent(in) :: a
         integer :: tsum

      end function tsum

      pure subroutine calc_T(T,xi,n,nquad)
         real(kind=rp),intent(inout) :: T(:,:)
         real(kind=rp),intent(in)    :: xi(:)
         integer      ,intent(in)    :: n, nquad
         real(kind=rp) :: xia(n)
         integer       :: i, j

         allocate(T(nquad,n))
         calc_T(:,1) = 1._rp
         xia = [((1._rp/(n-1))*j, j = 0,n-1)]
         do concurrent (i = 1:nquad, j = 2:n)
            T(i,j) = PRODUCT(xia(1:j-1)-xi(i))                           &
                   / PRODUCT(xia(1:j-1)-xia(j))
         end do
      end subroutine calc_T

      pure subroutine calc_T_xi(T_xi,xi,n,nquad)
         real(kind=rp),intent(inout) :: T_xi(:,:)
         real(kind=rp),intent(in)    :: xi(:)
         integer      ,intent(in)    :: n, nquad
         real(kind=rp) :: xia(n)
         integer       :: i, j

         T(:,1) = 0._rp
         xia = [((1._rp/(n-1))*j, j = 0,n-1)]
         do concurrent (i = 1:nquad, j = 2:n)
            do concurrent (k = 2:n)
               T_xi(i,j) = T_xi(i,j) + PRODUCT(xia(1:k-2)-xi(i))         &
                                     * PRODUCT(xia(k:j-1))
            end do
            T_xi(i,j) = T_xi(i,j)/ PRODUCT(xia(1:j-1)-xia(j))
         end do
      end subroutine calc_T_xi

      elemental function combine_T(ind)
         type(tuple(nsd+1)),intent(in) :: ind
         type(rtuple(nq)) :: combine_T
         integer :: i

         do concurrent (i = 1:nq)
            combine_T(i) = PRODUCT(T(i,ind%val))
         end do
      end function combine_T

      elemental function combine_T_xi(ind)
         type(tuple(nsd+1)),intent(in) :: ind
         type(rtuple(nq)) :: combine_T_xi
         integer :: i, p(nsd)

         do concurrent (i = 1:nq)
            p = [ind%val(1:isd-1),ind%val(isd+1:nsd+1)]
            q =  ind%val(1:nsd)
            combine_T_xi(i) = PRODUCT(T(i,p))*T_xi(i,isd)                &
                            - PRODUCT(T(i,q))*T_xi(i,nsd+1)
         end do
      end function combine_T_xi
   end subroutine simplex

end module shape_fns
