module element_library
   use param
   use elements
   implicit none
contains
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
      do concurrent (i = 1:nnod)
         xib = PACK(xia,ia/=i)
         L% sh(:,i) = eval_poly      (qr%pts,xib)
         L%dsh(:,i) = eval_poly_deriv(qr%pts,xib)
      end do
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

!====================== UTILITY PROCEDURES ===============================
   pure function eval_poly(x,roots)
      real(kind=rp),intent(in) :: x(:), roots(:)
      real(kind=rp) :: eval_poly(SIZE(x))

      eval_poly = PRODUCT(SPREAD(x    ,1,SIZE(roots))                    &
                -         SPREAD(roots,2,SIZE(x))    ,1)
   end function eval_poly

   pure function eval_poly_deriv(x,roots)
      real(kind=rp),intent(in) :: x(:), roots(:)
      real(kind=rp) :: eval_poly_deriv(SIZE(x))
      integer :: i, n

      n = SIZE(roots)
      do concurrent (i = 1,n)
         eval_poly_deriv = eval_poly_deriv                               &
                         + eval_poly(x,[roots(1:i-1),roots(i+1:n)])
      end do
   end function eval_poly_deriv
end module element_library
