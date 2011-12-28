      subroutine bbrsp(Ptsl,thispart,a,b,c,sa,sb,npt)
      implicit double precision (a-h,o-z)
      parameter (mxpt=10000)
      dimension  Ptsl(9,mxpt) 
      integer  thispart
      double precision  a,b,c,sa,sb,phi
      dimension  rotmat(3,3) 
      double precision  tmp1,tmp2,tmp3

      rotmat(1,1) = 1-sa*(b**2+c**2)
      rotmat(1,2) = a*b*sa+c*sb
      rotmat(1,3) = a*c*sa-b*sb
      rotmat(2,1) = a*b*sa-c*sb
      rotmat(2,2) = 1-sa*(a**2+c**2)
      rotmat(2,3) = b*c*sa+a*sb
      rotmat(3,1) = a*c*sa+b*sb
      rotmat(3,2) = b*c*sa-a*sb
      rotmat(3,3) = 1-sa*(a**2+b**2)

      tmp1 = Ptsl(7,thispart)
      tmp2 = Ptsl(8,thispart)
      tmp3 = Ptsl(9,thispart)

      Ptsl(7,thispart) = 
     > rotmat(1,1)*tmp1+rotmat(1,2)*tmp2+rotmat(1,3)*tmp3
      Ptsl(8,thispart) = 
     > rotmat(2,1)*tmp1+rotmat(2,2)*tmp2+rotmat(2,3)*tmp3
      Ptsl(9,thispart) = 
     > rotmat(3,1)*tmp1+rotmat(3,2)*tmp2+rotmat(3,3)*tmp3
      
      return
      end 
