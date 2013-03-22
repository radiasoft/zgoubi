      FUNCTION AGSMMA(DEV)
C     ------------------------------------------
C     Convert rigidity to amperes in main coils.
C     R. Thern, E. Blesser, AGS/AD/Tech.Note 424
C     ------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (AM4=-1.070694D3, AM3=+5.836957D2, AM2=-1.119693D2, 
     >AM1=+1.282979D0, A0=8.224488D0, A1=-3.408321D-4, A2=+7.313631D-6, 
     >A3=-5.642907D-8,A4=+1.966953D-10,A5=-3.143464D-13,A6=1.897916D-16)

      BL = BR * DEV
      BL1 = 1.D0 / BL
      AGSMMA = AM1 +(AM2 + (AM3 + AM4*BL1)*BL1)*BL1 + 
     >(A0 + (A1 + (A2 + (A3 + (A4 + (A5 + A6*BL)*BL)*BL)*BL)*BL)*BL)*BL
      RETURN
      END
