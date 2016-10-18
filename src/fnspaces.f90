#include "basic_types.h"

module fnspaces
   use global
   use grid
   use ploynomial_type
   implicit none

   type :: fnspace
      type(grid),   pointer,private :: nod
      type(element),pointer,private :: master(:)
   end type fnspace
end module fnspaces
