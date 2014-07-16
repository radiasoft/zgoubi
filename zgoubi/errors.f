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
      SUBROUTINE ERRORS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ********************************
C     READS DATA FOR PROCEDURE 'ERRORS'
C     ********************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      COMMON/DONT/ TA(MXL,MXTA)

      CHARACTER(132) TXT132
      LOGICAL STRCON, CMMNT, FITFNL, ok
      CHARACTER(40) STRA(10)
 
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLERR
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LBL1, LBL2
C                  XR, ZS... 
      CHARACTER(2) TYPERR
C                  A or R   G or U
      CHARACTER(1) TYPAR,   TYPDIS
      LOGICAL EMPTY
      INTEGER DEBSTR, FINSTR 

      DATA LBL1, LBL2 / 2*' ' /

C on/off switch  (1/0), number of lines to follow (each line sets a particular error)
      iop = nint(A(NOEL,1) )
      nbr = nint (A(NOEL,2) )
      iseed = nint (A(NOEL,3) )

      if (iop .eq. 0) then 
C        Switch off all possible earlier error settings
        call MULTP4
      endif

      if(nres.gt.0) then
        write(nres,fmt='(/25x,''--- SETTING ERRORS ---'',/)') 
        write(nres,fmt='(/15x,''On/off, number, random seed :'',3I9)') 
     >  iop, nbr, iseed
        write(nres,fmt='(/15x,''Errors to be introduced : '')')
        DO IRR = 1, NBR
          write(nres,fmt='(20x,a)')
     >    TA(NOEL,IRR)(debstr(TA(NOEL,IRR)):finstr(TA(NOEL,IRR)))
        enddo
      endif

      IF(NBR.GT.MXTA) CALL ENDJOB('SBR rerror. Number of instructions '
     >//' cannot exceed ',MXTA)

C Example of an error assignment line : 
C          MULTIPOL{lbl1,lbl2} 1, XR, R, G, center, sigma, cut
C {lbl1,lbl2} is optional, can be {,lbl2}, {lbl1}  
C 1 is the  pole # (dipole). can be 1-10 for dipole-20_pole
C XR is roll. Other possibilities : YR, ZR, XS, YS, ZS, BP (B_pole)
C R is for relative, A for absolute
C G for gaussian, U for uniform
C Case U : "sigma" stands for half-width
C cut is in units of sigma
      DO IRR = 1, NBR
        Txt132 = TA(NOEL,IRR)(debstr(TA(NOEL,IRR)):finstr(TA(NOEL,IRR)))
C         Get possible label1 and/or label2
        ok = strcon(txt132,'{',
     >                         is)
        if(ok) then 
          READ(txt132(1:is-1),*) klerr
        else
          READ(txt132,*) klerr
        endif
        if(klerr.eq.'MULTIPOL') then
          txt132 = txt132(9:finstr(txt132))
          if(ok) then 
            ok = strcon(txt132,'{',
     >                             is)
            ok = strcon(txt132,'}',
     >                             is2)
            ok = strcon(txt132(is:is2),',',
     >                                     is3)
            if(is+1.lt.is2-1) then
              if(.not. empty(txt132(is+1:is2-1))) 
     >             read(txt132(is+1:is2-1),*) lbl1
            endif
            if(ok) then
              if(.not. empty(txt132(is3+1:is2-1))) 
     >          read(txt132(2:is2+1),*) lbl2
            endif
            txt132 = txt132(is2+1:finstr(txt132))
          endif
C          Get the rest of the arguments
          txt132 = txt132(debstr(txt132):finstr(txt132))
          CALL STRGET(TXT132,99,
     >                          NSTR,STRA)
C          write(*,fmt='(20a)') ' sbr errors ',(stra(ii),ii=1,nstr)

          read(stra(1),*) ipol     ! Pole to which the error applies (1-10)
          read(stra(2),*) typerr   ! Error type : BP (B_pole), XR,YR,ZR,XS,YS.ZS
          read(stra(3),*) typar    ! Relative or absolute : R, A
          read(stra(4),*) typdis   ! Type of density : G or U,  or 0 to switch off ERRORS
          if(typdis.eq.'0') then
            call MULTP4
          else
            read(stra(5),*) errcen
            read(stra(6),*) errsig   ! sigma for G, half-width for U
            read(stra(7),*) errcut   ! in units of errsig for G, unused for U
            call MULTP2(irr,iseed,ipol,typerr,typar,typdis,
     >      errcen,errsig,errcut,lbl1,lbl2)          
          endif
        endif
      ENDDO

      RETURN
      END
