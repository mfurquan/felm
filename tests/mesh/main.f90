program main
   use global, only : fnam_max
   use mesh, only : zone
   implicit none

   type(zone) :: fluid
   character(len=fnam_max) :: nodfile,elemfile,bndryfile,outfile

   nodfile="test.nod"
   elemfile="test.ele"
   bndryfile="test.bnd"
   outfile="test.out"
   call fluid%read_zone(nodfile,elemfile,bndryfile)
   call fluid%print_zone(outfile)
end program main
