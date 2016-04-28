program test_integration
   use global
   use integration
   implicit none

   type(cubrule) :: lrule,trule,sqrule,tetrule,wdgrule,hexrule

   ! testing 2 point linear rule
   lrule   = cubrule(3)
   integrate(func1,lrule)
   ! testing cubic triangular rule
   trule   = cubrule(3,2)
   ! testing 4 point square rule
   sqrule  = lrule.otimes.lrule
   ! testing cubic tetrahedral rule
   tetrule = cubrule(3,3)
   ! testing cubic wedge rule
   wdgrule = trule.otimes.lrule
   ! testing 8 point hexahedral rule
   hexrule = sqrule.otimes.lrule

   contains

      function func1(x)
         implicit none
         real(kind=wp) :: x(:),func1

         func1=product(x)**3
      end function func1
end program test_integration
