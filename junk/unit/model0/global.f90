module global
   use iso_fortran_env
   implicit none

   integer,parameter :: wp=REAL64, &
                        nsd=3, ndf=1, &
                        fnam_max=20, elty_len=4
   logical,parameter :: debugOn=.FALSE.

end module global
