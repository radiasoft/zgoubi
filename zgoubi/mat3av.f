      subroutine MAT3AV(var,nr,v_aver,s_dev)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer*4 i,nr
      DIMENSION  var(*)
c
      v_aver=0
      s_dev=0.0
c
      Do i=1,nr
         v_aver=v_aver+var(i)
         enddo
c
         v_aver=v_aver/Float(nr)
c
      Do i=1,nr
         s_dev=s_dev+(var(i)-v_aver)*(var(i)-v_aver)
         enddo
c
	ffnum=Float(nr)
	ffnum=s_dev/ffnum
         s_dev=DSQRT(ffnum) 
c
      return
      end
