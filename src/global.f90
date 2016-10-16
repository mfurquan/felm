module global
   use iso_fortran_env
   implicit none

   integer,parameter :: rp = REAL64, &
                        si = selected_int_kind(2), &
                        li = selected_int_kind(8), &
                        var_len = 4
   logical,parameter :: debugOn=.FALSE.

end module global
