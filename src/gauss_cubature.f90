module gauss_cubature
   use global
   implicit none

   integer,private,parameter :: maxint = 6
   real(kind=wp),public,parameter :: gaus_pts(maxint) = [  &
      ! 1 point rule
       0.0,                                                &
      ! 2 point rule
      -1.0/sqrt(3.0),                                      &
       1.0/sqrt(3.0),                                      &
      !3 pointy rule 
      -sqrt(3.0/5.0),                                      &
       0.0,                                                &
       sqrt(3.0/5.0) ],                                    &
!
!     weights
      gaus_wts(maxint) = [                                 &
      !
      2.0,                                                 &
      !
      1.0,                                                 &
      1.0,                                                 &
      !
      5.0/9.0,                                             &
      8.0/9.0,                                             &
      5.0/9.0 ]
end module gauss_cubature
