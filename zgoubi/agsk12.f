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
      SUBROUTINE AGSK12(NOEL,X10,BK1,BK2,MOD,
     >                                       XL,BM,ANGMM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BK1(*),BK2(*),BM(*)
      INCLUDE 'MXFS.H'
      INCLUDE 'MXLD.H'
      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LABEL
      COMMON /LABEL/ LABEL(MXL,2)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF)
      SAVE YSHFT
      INTEGER FINSTR
      CHARACTER(LBLSIZ) LBL1
      PARAMETER (CM2M=1.D-2)

      lbl1 = label(noel,1)
        last = finstr(lbl1)
C        last = finstr(label(noel,2))

        if(MOD .eq. 1) then
C centered dipole model
          sag = -0.d0  !21906d0 * 2.54d0   /16.d0*18.d0
          sagS = -0.d0  !15472d0 * 2.54d0  /16.d0*18.d0
          sgn = -1.d0
          fb1 = 1.d0   !/1.0028d0
          fb2 = 1.d0/fb1
          if    (lbl1(last-1:last) .eq. 'BF') then
            sag = sagS
            xlmag = 79.d0 * 0.0254d0
            b3 = bk2(1) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(1) 
            b2 = bk1(1) - bk2(1)*(sag/2.d0)/100.d0
          elseif(lbl1(last-1:last) .eq. 'CD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(2) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(2) 
            b2 = bk1(2) - bk2(2)*(sag/2.d0)/100.d0 
          elseif(lbl1(last-1:last) .eq. 'AF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(3) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(3) 
            b2 = bk1(3) - bk2(3)*(sag/2.d0)/100.d0 
          elseif(lbl1(last-1:last) .eq. 'BD') then  
            sag = sagS
            xlmag = 79.d0 * 0.0254d0
            b3 = bk2(4) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(4) 
            b2 = bk1(4) - bk2(4)*(sag/2.d0)/100.d0 
          elseif(lbl1(last-1:last) .eq. 'CF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(5) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(5) 
            b2 = bk1(5) - bk2(5)*(sag/2.d0)/100.d0 
          elseif(lbl1(last-1:last) .eq. 'AD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(6) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(6) 
            b2 = bk1(6) - bk2(6)*(sag/2.d0)/100.d0 
          endif
          yce = 0.d0
        elseif(MOD .eq. 2) then
C long-shifted dipole model
          fyce = 1.d0
          sgn = 1.d0
          b1 = 0.d0
          fb1 = 1.d0
          fb2 = 1.d0
          if    (lbl1(last-1:last) .eq. 'BF') then
            xlmag = 79.d0 * 0.0254d0
            b3 = -bk2(1) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (23.17d0 )
            b2 = bk1(1) - bk2(1)*yce/100.d0
          elseif(lbl1(last-1:last) .eq. 'CD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = -bk2(2) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-24.07d0 )
            b2 = bk1(2) - bk2(2)*yce/100.d0 
          elseif(lbl1(last-1:last) .eq. 'AF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = -bk2(3) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (23.d0 )
            b2 = bk1(3) - bk2(3)*yce/100.d0 
          elseif(lbl1(last-1:last) .eq. 'BD') then  
            xlmag = 79.d0 * 0.0254d0
            b3 = -bk2(4) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-23.93d0 )
            b2 = bk1(4) - bk2(4)*yce/100.d0 
          elseif(lbl1(last-1:last) .eq. 'CF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = -bk2(5) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (23.d0 )
            b2 = bk1(5) - bk2(5)*yce/100.d0 
          elseif(lbl1(last-1:last) .eq. 'AD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = -bk2(6) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-24.07d0 )
            b2 = bk1(6) - bk2(6)*yce/100.d0 
          endif
c          fyce = 1.d0
c          sgn = 1.d0
c          fb1 = 1.d0
c          fb2 = 1.d0
c          if    (lbl1(last-1:last) .eq. 'BF') then
c            xlmag = 79.d0 * 0.0254d0
c            b3 = -bk2(1) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (23.17d0 )
c            b2 = bk1(1) 
c          elseif(lbl1(last-1:last) .eq. 'CD') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(2) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (-24.07d0 )
c            b2 = bk1(2) 
c          elseif(lbl1(last-1:last) .eq. 'AF') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(3) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (23.d0 )
c            b2 = bk1(3) 
c          elseif(lbl1(last-1:last) .eq. 'BD') then  
c            xlmag = 79.d0 * 0.0254d0
c            b3 = -bk2(4) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (-23.93d0 )
c            b2 = bk1(4) 
c          elseif(lbl1(last-1:last) .eq. 'CF') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(5) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (23.d0 )
c            b2 = bk1(5) 
c          elseif(lbl1(last-1:last) .eq. 'AD') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(6) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (-24.07d0 )
c            b2 = bk1(6) 
c          endif
c          b1 = b1 + b2*yce/100.d0 - b3 * yce*yce/1.d4
c          yce = -yce
        elseif(MOD .eq. 3) then
C short-shifted dipole model
          fyce = 1.d0
          sgn = -1.d0
          fb1 = 1.d0
          fb2 = 1.d0
c              write(*,*) ' agsk12 b3 bk2 avant ',b3,bk2(1)
          if    (lbl1(last-1:last) .eq. 'BF') then
            xlmag = 79.d0 * 0.0254d0
            b3 = bk2(1) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.3942767240d0)
            b2 = bk1(1) 
c              write(*,*) ' agsk12 b3 apres ',b3
          elseif(lbl1(last-1:last) .eq. 'CD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(2) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.5539466648d0)
            b2 = bk1(2) 
          elseif(lbl1(last-1:last) .eq. 'AF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(3) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.5589609130d0)
            b2 = bk1(3) 
          elseif(lbl1(last-1:last) .eq. 'BD') then  
            xlmag = 79.d0 * 0.0254d0
            b3 = bk2(4) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.3917462578d0)
            b2 = bk1(4) 
          elseif(lbl1(last-1:last) .eq. 'CF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(5) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.5589406517d0)
            b2 = bk1(5) 
          elseif(lbl1(last-1:last) .eq. 'AD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(6) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.5539271156d0)
            b2 = bk1(6) 
          endif
        else
          write(*,*) ' No such option  MOD = ',MOD
          stop
        endif

        XL = XLMAG / CM2M
        if    (lbl1(last-1:last) .eq. 'BF') then
          ANGMM = 0.02350230d0
        elseif(lbl1(last-1:last) .eq. 'BD') then  
          ANGMM = 0.02350230d0
        else
          ANGMM = 0.02796503d0
        endif
        BM(1) = -1.D3 / (XL / (2.D0*SIN(ANGMM/2.D0)))
        BM(2) = FB2*B2
c         write(*,*) ' agsk12 bm b2,b3 ',(bm(i),i=1,3),b2,b3
        BM(3) = B3*10.D0*SGN
c         write(*,*) ' agsk12 bm b2,b3 ',(bm(i),i=1,3),b2,b3
        YSHFT = YCE
               
      RETURN

      ENTRY AGSK11(
     >             YSHFTO)
      YSHFTO = YSHFT
      RETURN

      ENTRY AGSK13(MODI,NOEL,
     >                       YSHFTO)
      lbl1 = label(noel,1)
        last = finstr(lbl1)

        if(MODI .eq. 1) then
C centered dipole model
          yce = 0.d0
        elseif(MODI .eq. 2) then
C long-shifted dipole model
          fyce = 1.d0
          if    (lbl1(last-1:last) .eq. 'BF') then
            yce = fyce * (23.17d0 )
          elseif(lbl1(last-1:last) .eq. 'CD') then  
            yce = fyce * (-24.07d0 )
          elseif(lbl1(last-1:last) .eq. 'AF') then  
            yce = fyce * (23.d0 )
          elseif(lbl1(last-1:last) .eq. 'BD') then  
            yce = fyce * (-23.93d0 )
          elseif(lbl1(last-1:last) .eq. 'CF') then  
            yce = fyce * (23.d0 )
          elseif(lbl1(last-1:last) .eq. 'AD') then  
            yce = fyce * (-24.07d0 )
          endif
c          fyce = 1.d0
c          sgn = 1.d0
c          fb1 = 1.d0
c          fb2 = 1.d0
c          if    (lbl1(last-1:last) .eq. 'BF') then
c            xlmag = 79.d0 * 0.0254d0
c            b3 = -bk2(1) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (23.17d0 )
c            b2 = bk1(1) 
c          elseif(lbl1(last-1:last) .eq. 'CD') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(2) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (-24.07d0 )
c            b2 = bk1(2) 
c          elseif(lbl1(last-1:last) .eq. 'AF') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(3) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (23.d0 )
c            b2 = bk1(3) 
c          elseif(lbl1(last-1:last) .eq. 'BD') then  
c            xlmag = 79.d0 * 0.0254d0
c            b3 = -bk2(4) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (-23.93d0 )
c            b2 = bk1(4) 
c          elseif(lbl1(last-1:last) .eq. 'CF') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(5) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (23.d0 )
c            b2 = bk1(5) 
c          elseif(lbl1(last-1:last) .eq. 'AD') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(6) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (-24.07d0 )
c            b2 = bk1(6) 
c          endif
c          b1 = b1 + b2*yce/100.d0 - b3 * yce*yce/1.d4
c          yce = -yce
        else
C short-shifted dipole model
          fyce = 1.d0
          if    (lbl1(last-1:last) .eq. 'BF') then
            yce = fyce * (-0.3942767240d0)
          elseif(lbl1(last-1:last) .eq. 'CD') then  
            yce = fyce * (-0.5539466648d0)
          elseif(lbl1(last-1:last) .eq. 'AF') then  
            yce = fyce * (-0.5589609130d0)
          elseif(lbl1(last-1:last) .eq. 'BD') then  
            yce = fyce * (-0.3917462578d0)
          elseif(lbl1(last-1:last) .eq. 'CF') then  
            yce = fyce * (-0.5589406517d0)
          elseif(lbl1(last-1:last) .eq. 'AD') then  
            yce = fyce * (-0.5539271156d0)
          endif
        endif
        YSHFTO = YCE
      RETURN

      END
