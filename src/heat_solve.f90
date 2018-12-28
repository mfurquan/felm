program heat_solve
   use mesh_MOD, only : mesh
   use fe_space_MOD, only : fe_space
   implicit none

   type(mesh),target     :: plate
   type(fe_space),target :: Vh
   type(test_field)  :: w
   type(field)       :: u, f, h
   integer :: i

   call plate%read_file('holed_plate.msh',1,[(i,i=2,6)],2)
   call Vh%set(plate,'P1',1,1)
   call u%assoc(Vh)
   call f%assoc(Vh)

   call u%solve(int_dV(grad(w).ddot.grad(u)) &
                    .equals. int_dV(w.dot.f)+int_dS(w.dot.h))
   call u%write2file
end program heat_solve
