C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  François Méot
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 51 Franklin Street, Fifth Floor,
C  Boston, MA  02110-1301  USA
C
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE SYMPL(R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(6,*) , T(6,6,*)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG

      DIMENSION S(14)

C----- First order. Cours INSTN G. Leleux, p5-14.
C      conditions equivalent to T(R).S.R = S
      S(1) = R(1,1)*R(2,2)-R(1,2)*R(2,1)+R(3,1)*R(4,2)-R(3,2)*R(4,1)-1D0
      S(2) = R(1,3)*R(2,4)-R(1,4)*R(2,3)+R(3,3)*R(4,4)-R(3,4)*R(4,3)-1D0
      S(3)= R(2,4)*R(1,1)-R(1,4)*R(2,1) + R(4,4)*R(3,1)-R(3,4)*R(4,1)
      S(4)= R(2,4)*R(1,2)-R(1,4)*R(2,2) + R(4,4)*R(3,2)-R(3,4)*R(4,2)
      S(5)=-R(2,3)*R(1,1)+R(1,3)*R(2,1) - R(4,3)*R(3,1)+R(3,3)*R(4,1)
      S(6)=-R(2,3)*R(1,2)+R(1,3)*R(2,2) - R(4,3)*R(3,2)+R(3,3)*R(4,2)
      
      WRITE(NRES,100) (S(J),J=1,6)
 100  FORMAT(/,5X,
     >' First order symplectic conditions (expected values = 0) :',
     >/,5X,1P,6G14.4)
      RETURN

      ENTRY SYMPL2(R,T)

C----- Second order. These J.L. Laclare, p18.
      S(1) =    R(1,1)*(T(2,1,2)+T(2,2,1)) + 2.D0*R(2,2)*T(1,1,1) 
     > - 2.D0*R(1,2)*T(2,1,1) -    R(2,1)*(T(1,1,2)+T(1,2,1))

      S(2) = 2.D0*R(1,1)*T(2,2,2) +    R(2,2)*(T(1,1,2)+T(1,2,1))
     > -    R(1,2)*(T(2,1,2)+T(2,2,1)) - 2.D0*R(2,1)*T(1,2,2)

      S(3) =    R(1,1)*(T(2,2,6)+T(2,6,2)) +  R(2,2)*(T(1,1,6)+T(1,6,1))
     > -    R(1,2)*(T(2,1,6)+T(2,6,1)) -    R(2,1)*(T(1,2,6)+T(1,6,2))

      S(4) = 2.D0*R(1,1)*T(2,3,3) +    R(4,3)*(T(3,1,3)+T(3,3,1)) 
     > -    R(3,3)*(T(4,1,3)+T(4,3,1)) - 2.D0*R(2,1)*T(1,3,3)

      S(5) =    R(1,1)*(T(2,3,4)+T(2,4,3)) +  R(4,3)*(T(3,1,4)+T(3,4,1))
     > -    R(3,3)*(T(4,1,4)+T(4,4,1)) -    R(2,1)*(T(1,3,4)+T(1,4,3))

      S(6) =    R(1,1)*(T(2,3,4)+T(2,4,3)) +  R(4,4)*(T(3,1,3)+T(3,3,1))
     > -    R(3,4)*(T(4,1,3)+T(4,3,1)) -    R(2,1)*(T(1,3,4)+T(1,4,3))

      S(7) = 2.D0*R(1,1)*T(2,4,4) +    R(4,4)*(T(3,1,4)+T(3,4,1))
     > -    R(3,4)*(T(4,1,4)+T(4,4,1)) - 2.D0*R(2,1)*T(1,4,4)

      S(8) =    R(1,2)*(T(2,3,4)+T(2,4,3)) +  R(4,4)*(T(3,2,3)+T(3,3,2))
     > -    R(3,4)*(T(4,2,3)+T(4,3,2)) -    R(2,2)*(T(1,3,4)+T(1,4,3))

      S(9) = 2.D0*R(1,2)*T(2,4,4) +    R(4,4)*(T(3,2,4)+T(3,4,2)) 
     > -    R(3,4)*(T(4,2,4)+T(4,4,2)) - 2.D0*R(2,2)*T(1,4,4)

C FM 07/97      S(10)= 2.D0*R(1,2)*T(2,3,3) +    R(2,1)*(T(3,2,3)+T(3,3,2)) 
      S(10)= 2.D0*R(1,2)*T(2,3,3) +    R(4,3)*(T(3,2,3)+T(3,3,2)) 
     > -    R(3,3)*(T(4,2,3)+T(4,3,2)) - 2.D0*R(2,2)*T(1,3,3)

C FM 07/97       S(11)= 2.D0*R(1,2)*(T(2,3,4)+T(2,4,3))+ R(2,1)*(T(3,2,4)+T(3,4,2)) 
      S(11)=    R(1,2)*(T(2,3,4)+T(2,4,3))+ R(4,3)*(T(3,2,4)+T(3,4,2)) 
     > -    R(3,3)*(T(4,2,4)+T(4,4,2)) -    R(2,2)*(T(1,3,4)+T(1,4,3))

      S(12)=    R(3,3)*(T(4,1,4)+T(4,4,1)) +  R(4,4)*(T(3,1,3)+T(3,3,1)) 
     > -    R(3,4)*(T(4,1,3)+T(4,3,1)) -    R(4,3)*(T(3,1,4)+T(3,4,1))

      S(13)=    R(3,3)*(T(4,2,4)+T(4,4,2)) +  R(4,4)*(T(3,2,3)+T(3,3,2)) 
     > -    R(3,4)*(T(4,2,3)+T(4,3,2)) -    R(4,3)*(T(3,2,4)+T(3,4,2))

      S(14)=    R(3,3)*(T(4,4,6)+T(4,6,4)) +  R(4,4)*(T(3,3,6)+T(3,6,3)) 
     > -    R(3,4)*(T(4,3,6)+T(4,6,3)) -    R(4,3)*(T(3,4,6)+T(3,6,4))

      WRITE(NRES,101) (S(J),J=1,14)
 101  FORMAT(/,5X,
     >' Second order symplectic conditions (expected values = 0) :',
     >/,1P,2(/,5X,7G14.4) )

      RETURN
      END
