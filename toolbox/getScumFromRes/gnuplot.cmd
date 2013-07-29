
set tit "KICKH"
plot "getScumFromRes.out_KICKH" u :2 w p tit "zgoubi" ,\
     "/home/meot/zgoubi/struct/bnl/rhic/K.75/110509_includingLuoHVOrbits/filesFromLuo/list_blue_hkicker.dat" u :3 w p tit "Luo" 

pause 2

set tit "KICKV"
plot "getScumFromRes.out_KICKV" u :2 w p tit "zgoubi" ,\
     "/home/meot/zgoubi/struct/bnl/rhic/K.75/110509_includingLuoHVOrbits/filesFromLuo/list_blue_vkicker.dat" u :3 w p tit "Luo" 

pause 80

exit
