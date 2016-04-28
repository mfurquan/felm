program read_val
   integer v
   character(len=4) string
   string='TR03'
   read (string(3:4),'(I2)') v
   print *, v
end program read_val
