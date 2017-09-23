plot 'animation.fai' u ($10):($38==i?($11):1/0) w p title sprintf("Turn=%i",i)
print 'i=',i
pause 1
i=i+10
if (i < n) reread