program main
   use global, only : fnam_max
   use model, only : zone
   implicit none

   type(zone) :: fluid
   character(len=fnam_max) :: nodfile,elemfile,outfile

   nodfile="test.nod"
   elemfile="test.ele"
   outfile="test.out"
   call fluid%read_zone(nodfile,elemfile,5,2)
   call fluid%print_zone(outfile)
end program main
