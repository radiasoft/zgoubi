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
C  François Méot <fmeot@bnl.gov>
C  Brookhaven National Laboratory  
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE mulerr(noel,irr,iseed,BM, 
     >kPOL,TYPERR,TYPAR,TYPDIS,ERRCEN,ERRSIG,ERRCUT,
     >                                     DB,dpos,tilt)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension bm(*)
      PARAMETER (MXTA=45)
      PARAMETER (MXERR=MXTA)
      PARAMETER (MPOL=10)
      DIMENSION kpol(mxerr,mpol)
      character(2) TYPERR(mxerr,mpol)
      character(1) TYPAR(mxerr,mpol),TYPDIS(mxerr,mpol)
      DIMENSION ERRCEN(mxerr,mpol),ERRSIG(mxerr,mpol),ERRCUT(mxerr,mpol)

      INCLUDE 'MXLD.H'
      dimension db(MXL,mpol),dpos(MXL,mpol,3),tilt(MXL,mpol,3)

      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QEL,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN

          call rndm2(iseed)      !!!!!!!!!! seed

      do i = 1, mpol
        if(kpol(irr,i) .eq.1) then
          if    (typdis(irr,i).eq.'G') then
            SM = EXP(-(errcut(irr,i)*errsig(irr,i))**2/2.D0)
            derr = (2.d0*RNDM() -1.d0)*SM
          elseif(typdis(irr,i).eq.'U') then
            derr = errsig(irr,i)* 2.d0*(rndm()-0.5d0) 
          endif
          if    (typerr(irr,i)(1:1).eq.'B') then
            if    (typar(irr,i).eq.'A') then
C              Absolute error
              db(noel,i) = errcen(irr,i) + derr
            elseif(typar(irr,i).eq.'R') then
C              Relative error
              db(noel,i) = errcen(irr,i) + derr*BM(I)
            endif
          elseif(typerr(irr,i)(2:2).eq.'S') then
              dpos(noel,i,3) = 0.d0
          elseif(typerr(irr,i)(2:2).eq.'R') then
              tilt(noel,i,3) = 0.d0
          else
            call endjob('SBR mulerr. No such option for typerr',-99)
          endif
        endif
      enddo      

      RETURN      
      END
