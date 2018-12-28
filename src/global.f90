module global
   use iso_fortran_env, only : REAL64
   implicit none

   integer,parameter :: rp = REAL64, &
!                        si = selected_int_kind(2), &
!                        li = selected_int_kind(8), &
                        nen_max =  4, &
                        nsd_max =  3, &
                        var_len =  4, &
                        ele_len =  4

   logical,parameter :: debugOn=.TRUE.

end module global
