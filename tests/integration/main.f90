program test_integration
   use global
   use integration
   implicit none

   type(cubrule) :: lrule,trule,sqrule,tetrule,wdgrule,hexrule
   real(kind=wp) :: integral

   ! testing 2 point linear rule
   lrule   = cubrule(3)
   integral = integrate(func1,lrule)
   write(*,*) 'integral over line =',integral

   ! testing cubic triangular rule
   trule   = cubrule(3,2)
   integral = integrate(func2,trule)
   write(*,*) 'integral over triangle =',integral

   ! testing 9 point square rule
   sqrule  = lrule.otimes.lrule
   integral = integrate(func2,sqrule)
   write(*,*) 'integral over square =',integral

   ! testing cubic tetrahedral rule
   tetrule = cubrule(3,3)
   integral = integrate(func3,tetrule)
   write(*,*) 'integral over tetrahedron =',integral

   ! testing cubic wedge rule
   wdgrule = trule.otimes.lrule
   integral = integrate(func3,wdgrule)
   write(*,*) 'integral over wedge =',integral

   ! testing 8 point hexahedral rule
   hexrule = sqrule.otimes.lrule
   integral = integrate(func3,hexrule)
   write(*,*) 'integral over hexahedron =',integral

   contains

      pure function func1(x)
         implicit none
         real(kind=wp),intent(in) :: x(:)
         real(kind=wp) :: func1

         func1=x(1)**3
      end function func1

      pure function func2(x)
         implicit none
         real(kind=wp),intent(in) :: x(:)
         real(kind=wp) :: func2

         func2=x(1)**2+x(2)**2
      end function func2

      pure function func3(x)
         implicit none
         real(kind=wp),intent(in) :: x(:)
         real(kind=wp) :: func3

         func3=product(x)
      end function func3
end program test_integration
