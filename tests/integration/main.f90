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

   ! testing 4 point square rule
   sqrule  = lrule.otimes.lrule
   integral = integrate(func1,sqrule)
   write(*,*) 'integral over square =',integral

   ! testing cubic tetrahedral rule
   tetrule = cubrule(3,3)
   integral = integrate(func1,tetrule)
   write(*,*) 'integral over tetrahedron =',integral

   ! testing cubic wedge rule
   wdgrule = trule.otimes.lrule
   integral = integrate(func1,wdgrule)
   write(*,*) 'integral over wedge =',integral

   ! testing 8 point hexahedral rule
   hexrule = sqrule.otimes.lrule
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

         func2=x(1)**2*x(2)
      end function func2
end program test_integration
