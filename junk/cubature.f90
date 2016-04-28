module cubature
   use global
   implicit none

!integrates a function over a reference element
   function integ(func,eltyp,ncub)
      integer,parameter :: maxint
      character(len=elty_len),intent(in) :: eltyp
      integer,intent(in) :: ncub, i, j, n, zint=maxint**(1/3)
      real(kind=wp) :: wt(maxint), xi(nsd,maxint), z(zint)
      interface
         function func(zi)
            real(kind=wp) :: zi(nsd),func(ndf)
         end function
      end interface

      xi=0.; wt=0.
      select case(eltyp(1:2))
      case('LN')
         select case(ncub)
         case(1)
            wt(1)=2.; xi(:,1)=[0., 0., 0.]
         case(2)
            wt(1)=1.; xi(:,1)=[-1./sqrt(3.), 0., 0.]
            wt(2)=1.; xi(:,2)=[ 1./sqrt(3.), 0., 0.]
         case(3)
            wt(1)=5./9.; xi(:,1)=[-sqrt(3./5.), 0., 0.]
            wt(2)=8./9.; xi(:,2)=[0.,0.,0.]
            wt(3)=5./9.; xi(:,3)=[ sqrt(3./5.), 0., 0.]
         case default
            write(*,*) 'No rule for',ncub, &
                       'cubature points in',eltyp,'element'
         end select

      case('TR')
         select case(ncub)
         case(1)
            wt(1)=1.; xi(:,1)=[1./3., 1./3., 0.]
         case(3)
            wt(1)=1./3.; xi(:,1)=[1./6., 1./6., 0.]
            wt(2)=1./3.; xi(:,2)=[2./3., 1./6., 0.]
            wt(3)=1./3.; xi(:,3)=[1./6., 2./3., 0.]
         case(4)
            wt(1)=-27./48.; xi(:,1)=[1./3., 1./3., 0.]
            wt(2)= 25./48.; xi(:,2)=[1./5., 1./5., 0.]
            wt(3)= 25./48.; xi(:,3)=[3./5., 1./5., 0.]
            wt(4)= 25./48.; xi(:,4)=[1./5., 3./5., 0.]
         case default
            write(*,*) 'No rule for',ncub, &
               'cubature points in',eltyp,'element'
         end select

      case('QD')
         select case(ncub)
         case(1)
            wt(1)=4.; xi(:,1)=[0., 0., 0.]
         case(4)
            wt(1)=1.; xi(:,1)=[-1./sqrt(3.), -1./sqrt(3.), 0.]
            wt(2)=1.; xi(:,2)=[ 1./sqrt(3.), -1./sqrt(3.), 0.]
            wt(3)=1.; xi(:,3)=[ 1./sqrt(3.),  1./sqrt(3.), 0.]
            wt(4)=1.; xi(:,4)=[-1./sqrt(3.),  1./sqrt(3.), 0.]
         case(9)
            n=0
            zw(1:3)=[5./9., 8./9., 5./9.]
            z=
            do i=1,3
               do j=1,3
                  n=n+1
                  wt(n)=zw(i)*zw(j)
                  xi(:,n)=[z(i), z(j), 0.]
               end do
            end do
   end function integ
end module cubature
