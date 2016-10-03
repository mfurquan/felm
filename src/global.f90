module global
   use iso_fortran_env
   implicit none

   integer,parameter :: wp = REAL64, &
                        nsd = 3, &
                        fnam_max = 20, celltype_len = 4, nen_max = 20, &
                        var_len = 4, shap_len = 8
   logical,parameter :: debugOn=.FALSE.

end module global
