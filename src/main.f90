program poiss
   use mesh
   use global
   implicit none

   type(zone) :: domain
   character(len=fnam_max) :: nodfile,elemfile

   call read_inp(nodfile,elemfile)
   call domain%read_zone(nodfile,elemfile)
end program poiss
