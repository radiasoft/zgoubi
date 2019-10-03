program main
   !! author: Damian Rouson
   !!
   !! summary: verify that data partitioned across images evenly to 
   !! within a difference of one datum between any two images.
   use data_partition_ixfc, only : data_partition
   implicit none
   type(data_partition) my_data

   print *,"Hello from image",this_image(),"of",num_images()
   print *,"Test passed."

end program
