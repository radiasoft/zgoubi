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
C  Upton, NY, 11973, USA
C  -------
      BLOCK DATA BLOCK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
C----- Pick-ups
      PARAMETER (MXPUD=9,MXPU=5000)
       INCLUDE "C.CO.H"     ! COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
C----- CONSTANTES
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
C--------
      INCLUDE "C.DEPL.H"     ! COMMON/DEPL/ XF(3),DXF(3),DQBRO,DTAR
C      PARAMETER (MDR=9)
      INCLUDE "C.DROITE_2.H"     ! COMMON/DROITE/ AM(MDR),BM(MDR),CM(MDR),IDRT
      INCLUDE "C.EFBS.H"     ! COMMON/EFBS/ AFB(2), BFB(2), CFB(2), IFB
      INCLUDE "MAXTRA.H"
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.GASC.H"     ! COMMON/GASC/ AI, DEN, KGA
      CHARACTER(LEN=1) KAR(41)
      INCLUDE "C.KAR.H"     ! COMMON/KAR/ KAR
      INCLUDE 'MXLD.H'
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) LABEL
      INCLUDE "C.LABEL.H"     ! COMMON/LABEL/ LABEL(MXL,2)
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      INCLUDE 'C.PDATA.H'       ! COMMON /PDATA/ AMLEC,GLEC,AMPRO,GPRO,AMMU,GMU,TAUMU,AM3HE,G3HE,
                                ! AMDEU,GDEU
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AMASS,Q,G,TO
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      INCLUDE "C.SCAL.H"     ! COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
      LOGICAL TSPCH
      INCLUDE "C.SPACECHA.H" ! COMMON/SPACECHA/ TLAMBDA,RBEAM(2),XAVE(2),EMITT(2),TAVE,BUNCH_LEN,
C                             >                EMITTZ, BTAG, SCKX, SCKY, TSPCH
      INCLUDE "C.STEP.H"     ! COMMON/STEP/ TPAS(3), KPAS
      INCLUDE "C.SYNRA.H"     ! COMMON/SYNRA/ KSYN
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
C----- CONVERSION DES COORD. (CM,MRD) -> (M,RD)
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)
      PARAMETER (MXV=60)
      INCLUDE "C.VARY.H"  ! COMMON/VARY/ NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
                          !     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,27)

      PARAMETER (MDR3= 3*MDR)

      DATA LF,LST / 2 * 0 /
      DATA KCO / 0 /
      DATA XF,DXF,DQBRO / 3*0.D0, 3*0.D0, 0.D0 /
      DATA ZERO, UN / 0.D0, 1.D0 /
      DATA IDRT, AM, BM, CM / 0, MDR3*0.D0 /
      DATA IFB / 0 /
      DATA IFDES / 0 /
      DATA KGA / 0 /
      DATA AMQLU,PABSLU/6*.FALSE./
      DATA (KAR(I),I=1,41) /
     > 'O','A','B','C','D','E','F','G','H','I','J','K','L','M','N'
     >,'P','Q','R','U','V','W','X','Y','Z','2','3','4','5','6'
     >,'7','8','9','0','(',')','+','-','/','=','"','*'/
      DATA KZOB, KOBJ / 0, 0 /
      DATA KFLD,MG,LC,ML,ZSYM/ 1,1,2,3,.TRUE./
      DATA IDS, KORD / 4, 2 /
      DATA AMASS, Q / 0.D0, 1D0 /
      DATA NRBLT,IPASS/ 0, 1/
Changed dpref to dp/dp_ref - integer part. FM Mar 2018
C     DATA DPREF / 1.D0 /
      DATA DPREF, HDPRF / 0.D0, 1.D0 /
      DATA KPAS / 0 /
      DATA CPAR / 1620*0.D0 /
      DATA IDMAX / 1 /
      DATA KSYN / 0 /

C ----- Fundamental Physical Constants -----
C Updated 08/2018 - FM; 09/2018 - DTA
C Unless noted otherwise, these data are taken from the database
C of 2014 CODATA Recommended Values, available online at
C   NIST, https://physics.nist.gov/cuu/Constants/index.html
C Parentheses delimit the standard uncertainty in the last digits.
C
C The magnetic moment anomaly we require is the quantity (g - 2) / 2
C traditionally denoted by either a (leptons) or G (baryons) in the
C well-known Thomas-BMT equation of spin dynamics. Here g denotes a
C so-called g-factor, which is a dimensionless constant that measures
C the extent to which a particle's gyromagnetic ratio (the ratio of
C magnetic moment to spin angular momentum) differs from the classical
C value of q / (2 m). More specifically, one writes the gyromagnetic
C ratio, denoted γ, in the form
C   γ = g (q / (2 m)),
C where q denotes the (signed) particle charge and m the particle mass.
C One may therefore (I'm leaving out a few steps here!) compute the
C required g-factor using the formula
C   g = (u / u_N) * (m / m_p) / (Z * S / h-bar),
C where
C         u   = particle magnetic moment,
C         u_N = nuclear magneton = (e * h-bar) / (2 * m_p),
C         m   = particle mass,
C         m_p = proton mass,
C         Z   = particle charge in units of elementary charge e,
C   S / h-bar = particle spin (max. proj. of S_z) in units of h-bar.
C
C Significant confusion may arise from the fact that the database of
C CODATA Recommended Values lists not the g-factor we require, but a
C different g-factor. A useful reference on the latter is the article
C   PJ Mohr, DB Newell, and BN Taylor, “CODATA recommended values of
C   the fundamental physical constants: 2014”, Rev. Modern Phys.,
C   88(3), July–Sept. 2016; DOI: 10.1103/revmodphys.88.035009.
C One may download this article from
C   https://physics.nist.gov/cuu/Constants/article2014.html
C See, in particular, the introductory portions of sections V and VI.
C
C For a lepton---electron, muon, tau---the differing definitions
C affect only the sign of g, and one may compute the desired magnetic
C moment anomaly as (|g| - 2) / 2. Because of the importance of these
C values to an understanding of QED, both electron amd muon magnetic
C anomalies are included in the CODATA Recommended Values.

C For nucleons, the g-factor (here denoted g_n) used by the nuclear
C physics community is defined by the rule
C   u = g_n * (e / 2 * m_p) * S = g_n * u_N * (S / h-bar),
C It follows that the two different g-factors are related according to
C   g_s = g_n * (m_n / m_p) * (1 / Z).
C Here g_s denotes the g-factor we need for computing spin dynamics
C with the Thomas-BMT equation.

C
C speed of light in vacuum / m.s^-1 (exact)
      DATA CL / 2.99792458D+08 /
C
C elementary charge / C (98)
      DATA QEL / 1.602176487D-19 /   ! CODATA 2006
C     DATA QEL / 1.6021766208D-19 /  ! CODATA 2014
C
C atomic mass unit energy equivalent u.c^2 / MeV (57)
      PARAMETER (AMU = 931.4940954D0)
C
C electron, spin +1/2
C electron mass energy equivalent = 0.5109989461(31) MeV
C electron mass = 548.579909070(16) x 10^-6 u
      PARAMETER (XMLEC = 548.579909070D-06 * AMU) ! = 0.5109989461|537738
      DATA AMLEC / XMLEC /
C electron g-factor = 2.00231930436182(52)
C electron magnetic moment anomaly / 1 (26) :: a = (g-2)/2
      DATA GLEC / 1.159652181D-3 /     ! c. CODATA 2006
C     DATA GLEC / 1.15965218091D-03 /  ! CODATA 2014
C
C muon, spin +1/2
C muon mass energy equivalent = 105.6583745(24) MeV
C muon mass = 0.1134289257(25) u
      PARAMETER (XMMU = 0.1134289257D0 * AMU) ! = 105.6583745|371153
      DATA AMMU / XMMU /
C muon g-factor = 2.0023318418(13)
C muon magnetic moment anomaly / 1 (63) :: a = (g-2)/2
      DATA GMU / 1.16592089D-03 /
C muon lifetime / s (22)
C Ref: http://pdg.lbl.gov/2018/tables/rpp2018-sum-leptons.pdf
C According to this reference, the ratio of lifetimes for positive and
C negative muons is very nearly unity: tau_mu+ / tau_mu- = 1.00002(08).
C In other words, experiment cannot yet say that the anti-particle
C lifetime differs at all from the particle lifetime.
      DATA TAUMU / 2.1969811D-06 /
C
C proton, spin +1/2
C proton mass energy equivalent = 938.2720813(58) MeV
C proton mass = 1.007276466879(91) u
      PARAMETER (XMPRO = 1.00727646688D0 * AMU)  ! c. CODATA 2014
C     PARAMETER (XMPRO = 1.007276466879D0 * AMU) ! = 938.2720813|33162
      DATA AMPRO / XMPRO /
      DATA AMPROT / XMPRO /
C proton g-factor = 5.585694702(17)
C proton magnetic moment anomaly / 1 (85) :: G = (g-2)/2
      DATA GPRO / 1.79284735D0 /    ! c. CODATA 2006
C     DATA GPRO / 1.7928473508D0 /  ! CODATA 2014
C
C deuteron, spin +1 (boson)
C deuteron mass energy equivalent = 1875.612928(12) MeV
C deuteron mass = 2.013553212745(40) u
      PARAMETER (XMDEU = 2.013553212745D0 * AMU) ! = 1875.612928|445668
      DATA AMDEU / XMDEU /
C deuteron g-factor = 1.7140254555(98)
C deuteron magnetic moment anomaly / 1 (49) :: G = (g-2)/2
      DATA GDEU / -0.14301D0 /
C     DATA GDEU / -0.1429872722D0 /  ! DTA
C
C helion, spin +1/2
C helion mass energy equivalent = 2808.391586(17) MeV
C helion mass = 3.01493224673(12) u
      PARAMETER (XM3HE = 3.01493224673D0 * AMU) ! = 2808.391585|86005
      DATA AM3HE / XM3HE /
C helion g-factor = -6.368307372(74)
C helion magnetic moment anomaly / 1 (37) :: G = (g-2)/2
      DATA G3HE / -4.1841538D0 /
C     DATA G3HE / -4.184153686D0 /  ! DTA

      DATA TSPCH / .FALSE. /

C----- To yield MKSA units :
C                                1      2     3     4     5    6     7
C                                Y      T     Z     P     S    D    time
C                                m     rad    m    rad    m    1     s
      DATA (UNIT(I),I=1,MXJ) / 1.D-2,1.D-3,1.D-2,1.D-3,1.D-2,1.D0,1.D-6/
      DATA (NTIM(I),I=1,MXF) / MXF * 0 /
      PARAMETER (MXL2=MXL*2)
      DATA LABEL / MXL2*' ' /
      END
