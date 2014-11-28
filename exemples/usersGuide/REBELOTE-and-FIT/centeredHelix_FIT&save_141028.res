Centering 5 helical orbits in the AGS warm helical snake 3-D OPERA map.                                       
 'OBJET'                    ! This data list may be copy-pasted and run, as is.                               1
7.2051782956D3                                                                                                
2                                                                                                             
1 1                                                                                                           
-2.2  0. 0. 0. 0.  1. 'o'                                                                                     
1                                                                                                             
 'PARTICUL'                                                                                                   2
938.27203 1.602176487E-19 1.7928474 0 0                                                                       
                                                                                                              
 'SPNTRK'                                                                                                     3
4.1                                                                                                           
0. 0. 1.                                                                                                      
                                                                                                              
 'FAISCEAU'                                                                                                   4
 'SPNPRT'                                                                                                     5
                                                                                                              
 'TOSCA'                                                                                                      6
0  20                                                                                                         
1.e1   100. 100. 100.                                                                                         
HEADER_0 wsnake                                                                                               
801 29 29 12.1                                                                                                
warmSnake.map                                                                                                 
0 0 0 0                                                                                                       
2                                                                                                             
.1                                                                                                            
2  0.  .0  0.  0.                                                                                             
                                                                                                              
 'FAISCEAU'                                                                                                   7
                                                                                                              
 'FIT'                                                                                                        8
2  save                     ! Save updated FIT variables in zgoubi.FITVALS.out.                               
1 30 0  [-3,3]                                                                                                
1 31 0  [-3,3]                                                                                                
3    1E-2                                                                                                     
3.1 1 2 6 0. 1. 0                                                                                             
3.1 1 3 6 0. 1. 0                                                                                             
7.3 1 2 6 0. 1. 0                                                                                             
                                                                                                              
 'SPNPRT'  PRINT            ! zgoubi.SPNPRT.Out reamins open in the presence of                               9
! REBELOTE, output data will be stacked.                                                                      
                                                                                                              
 'SYSTEM'                   ! Save zgoubi.FITVALS.out data                                                   10
1                                                                                                             
cat zgoubi.FITVALS.out >> zgoubi.FITVALS.out_cat                                                              
                                                                                                              
 'REBELOTE'                                                                                                  11
4  0.1  0 1                                                                                                   
1                                                                                                             
OBJET  35  1.3872739973  2.1368296674  4.8261190694  11.015241321                                             
                                                                                                              
 'FAISCEAU'                                                                                                  12
                                                                                                              
 'SYSTEM'                   ! Save zgoubi.SPNPRT.Out                                                         13
1                                                                                                             
cp zgoubi.SPNPRT.Out zgoubi.SPNPRT.Out_save                                                                   
                                                                                                              
 'END'                                                                                                       14

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


     Particle  properties :

                     Mass          =    938.272        MeV/c2
                     Charge        =   1.602176E-19    C     
                     G  factor     =    1.79285              

              Reference  data :
                    mag. rigidity (kG.cm)   :   7205.1783      =p/q, such that dev.=B*L/rigidity
                    mass (MeV/c2)           :   938.27203    
                    momentum (MeV/c)        :   2160.0581    
                    energy, total (MeV)     :   2355.0383    
                    energy, kinetic (MeV)   :   1416.7663    
                    beta = v/c              :  0.9172072069    
                    gamma                   :   2.509973905    
                    beta*gamma              :   2.302166155    
                    G*gamma                 :   4.500000190    
                    electric rigidity (MeV) :   1981.220867      =T[eV]*(gamma+1)/gamma, such that dev.=E*L/rigidity
  
 I, AMQ(1,I), AMQ(2,I)/QE, P/Pref, v/c, time :
  
     1   9.38272030E+02  1.00000000E+00  1.00000000E+00  9.17207207E-01  0.00000000E+00

************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              All spins entered particle by particle
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.0000    -2.200     0.000     0.000     0.000       0.0000    1.0000   -2.200    0.000    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -2.2000E-02    0.000           1       0    0.000    B-Dim     1      1
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      1
    0.000        0.000        0.000        0.000        1417.           1       0    0.000    B-Dim     3      1


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                            



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000              0.000000    0.000000    1.000000    1.000000    0.000000    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000           0.000000  0.000000  1.000000  1.000000      2.5100    0.0000    0.0000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.00000E+00  2.50997E+00      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


                OPEN FILE zgoubi.plt                                                                      
                FOR PRINTING TRAJECTORIES


     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.I is 12. 1


           New field map(s) now used, cartesian mesh (MOD.le.19) ; 
           name(s) of map data file(s) : 

          warmSnake.map                                                                   

warmSnake.map map,  FORMAT type :  regular            

 HEADER  (0 lines) : 
  SBR FMAPW/FMAPR3 : completed job of reading map.


     Min/max fields in map       :  -2.934230E-01                  /   3.119730E-01
       @  X(CM),  Y(CM), Z(CM) :  -1.04     -6.000E-02  6.500E-02  /   1.05      5.500E-02  7.000E-02
      given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     Min/max normalised fields (kG)  : -2.934230E+00                      3.119730E+00

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm


     Element  is  mis-aligned  wrt.  the  optical  axis
          Center  of  entrance  EFB  is  at    X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD

  A    1  1.0000    -2.200     0.000     0.000     0.000          599.500    -2.182     0.000    -0.055     0.000            1


                CONDITIONS  DE  MAXWELL  (     4004.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -3.3328E-12      -1.0454E-13       -1.5080E-14
                                      1.5768E-12       -9.5122E-13
                                     -2.7052E-13       -6.0426E-12
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.0000    -2.200     0.000     0.000     0.000       0.0000    1.0000   -2.182    0.090   -0.055    0.009   4.003418E+02     1
               Time of flight (mus) :  1.45593751E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -2.1817E-02   9.0396E-05       1       0    0.000    B-Dim     1      1
    0.000        0.000        0.000      -5.4741E-04   8.6225E-06       1       0    0.000    B-Dim     2      1
    0.000        0.000        0.000       1.4559E-02    1417.           1       0    0.000    B-Dim     3      1


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      8  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           constraint #            1       IR =            6 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            6 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            6 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will be saved in zgoubi.FITVALS.out                                                              


 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    1    1     30   -3.00       -1.84       -1.840000000       3.00      2.000E-02 OBJET      *          *         
    1    2     31   -3.00       6.000E-02   6.0000000000E-02   3.00      2.000E-02 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     6      0.0000000E+00     1.0000E+00    8.1418231E-03   7.2565E-03   TOSCA      *          *          *   0 : 
   3   1   3     6      0.0000000E+00     1.0000E+00    8.1629180E-02   7.2942E-01   TOSCA      *          *          *   0 : 
   7   1   2     6      0.0000000E+00     1.0000E+00   -4.9046288E-02   2.6333E-01   TOSCA      *          *          *   0 : 
 Fit reached penalty value   9.1352E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              All spins entered particle by particle
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.0000    -1.840     0.060     0.000     0.000       0.0000    1.0000   -1.840    0.060    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -1.8400E-02   6.0000E-05       1       0    0.000    B-Dim     1      1
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      1
    0.000        0.000        0.000        0.000        1417.           1       0    0.000    B-Dim     3      1


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                            



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000              0.000000    0.000000    1.000000    1.000000    0.000000    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000           0.000000  0.000000  1.000000  1.000000      2.5100    0.0000    0.0000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.00000E+00  2.50997E+00      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.I is 12. 1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           warmSnake.map                                                                   
 SBR TOSCAC,  restored mesh coordinates for field map #   1,  name : warmSnake.map


     Min/max fields in map       :  -2.934230E-01                  /   3.119730E-01
       @  X(CM),  Y(CM), Z(CM) :  -1.04     -6.000E-02  6.500E-02  /   1.05      5.500E-02  7.000E-02
      given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     Min/max normalised fields (kG)  : -2.934230E+00                      3.119730E+00

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm


     Element  is  mis-aligned  wrt.  the  optical  axis
          Center  of  entrance  EFB  is  at    X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD

  A    1  1.0000    -1.840     0.060     0.000     0.000          599.500    -1.832     0.000    -0.045     0.000            1


                CONDITIONS  DE  MAXWELL  (     4004.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -2.7838E-12      -1.0454E-13        6.5642E-15
                                      1.5768E-12       -9.5291E-13
                                     -2.7052E-13       -6.0517E-12
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    8.00000000     m ;  Time  (for ref. rigidity & particle) =   2.909389E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.0000    -1.840     0.060     0.000     0.000       0.0000    1.0000   -1.832   -0.022   -0.045    0.007   4.003421E+02     1
               Time of flight (mus) :  1.45593861E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -1.8319E-02  -2.1629E-05       1       0    0.000    B-Dim     1      1
    0.000        0.000        0.000      -4.5196E-04   6.8008E-06       1       0    0.000    B-Dim     2      1
    0.000        0.000        0.000       1.4559E-02    1417.           1       0    0.000    B-Dim     3      1


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************

      8   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
      9  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000             -0.000966    0.208322    0.978060    1.000000   12.024168    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000966  0.208322  0.978060  1.000000      2.5100   12.0242   12.0240    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  1.00000E+00  2.50997E+00      1   1

************************************************************************************************************************************
     10  Keyword, label(s) :  SYSTEM                            

     Number  of  commands :   1,  as  follows : 

 cat zgoubi.FITVALS.out >> zgoubi.FITVALS.out_cat

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        1 through the optical structure 

                     Total of          1 particles have been launched

     Multiple pass, 
          from element #     1 : OBJET     /label1=          /label2=           to REBELOTE /label1=          /label2=          
          ending at pass #       5 at element #    11 : REBELOTE  /label1=          /label2=          


     Parameter #   35 in element #    1 will be modified at each pass. 
     list of requested parameter values :
                   1   1.38727400E+00
                   2   2.13682967E+00
                   3   4.82611907E+00
                   4   1.10152413E+01

 Pgm rebel. At pass #    1/   5.  In element #    1,  parameter # 35  changed to    1.38727400E+00   (was    1.00000000E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    1    1     30   -3.00       -1.30       -1.300000000       3.00      2.000E-02 OBJET      *          *         
    1    2     31   -3.00       0.120       0.1200000000       3.00      2.000E-02 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     6      0.0000000E+00     1.0000E+00    3.6561969E-02   2.2023E-01   TOSCA      *          *          *   0 : 
   3   1   3     6      0.0000000E+00     1.0000E+00    5.9165583E-02   5.7671E-01   TOSCA      *          *          *   0 : 
   7   1   2     6      0.0000000E+00     1.0000E+00    3.5107390E-02   2.0306E-01   TOSCA      *          *          *   0 : 
 Fit reached penalty value   6.0699E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              All spins entered particle by particle
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.3873    -1.300     0.120     0.000     0.000       0.0000    1.3873   -1.300    0.120    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -1.3000E-02   1.2000E-04       1       0    0.000    B-Dim     1      2
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      2
    0.000        0.000        0.000        0.000        2202.           1       0    0.000    B-Dim     3      2


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                            



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000              0.000000    0.000000    1.000000    1.000000    0.000000    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000           0.000000  0.000000  1.000000  1.000000      3.3466    0.0000    0.0000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  1.38727E+00  3.34663E+00      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.I is 12. 1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           warmSnake.map                                                                   
 SBR TOSCAC,  restored mesh coordinates for field map #   1,  name : warmSnake.map


     Min/max fields in map       :  -2.934230E-01                  /   3.119730E-01
       @  X(CM),  Y(CM), Z(CM) :  -1.04     -6.000E-02  6.500E-02  /   1.05      5.500E-02  7.000E-02
      given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     Min/max normalised fields (kG)  : -2.934230E+00                      3.119730E+00

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm


     Element  is  mis-aligned  wrt.  the  optical  axis
          Center  of  entrance  EFB  is  at    X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD

  A    1  1.3873    -1.300     0.120     0.000     0.000          599.500    -1.263     0.000    -0.038     0.000            1


                CONDITIONS  DE  MAXWELL  (     4002.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -1.4120E-12      -5.6219E-14       -3.7321E-15
                                      1.0172E-12       -6.7881E-13
                                     -1.4808E-13       -4.3143E-12
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    8.00000000     m ;  Time  (for ref. rigidity & particle) =   2.909389E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   1.3873    -1.300     0.120     0.000     0.000       0.0000    1.3873   -1.263    0.061   -0.038    0.002   4.001769E+02     1
               Time of flight (mus) :  1.39875067E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -1.2634E-02   6.0834E-05       1       0    0.000    B-Dim     1      2
    0.000        0.000        0.000      -3.7661E-04   2.1120E-06       1       0    0.000    B-Dim     2      2
    0.000        0.000        0.000       1.3988E-02    2202.           1       0    0.000    B-Dim     3      2


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************

      8   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
      9  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000             -0.000671    0.196147    0.980574    1.000000   11.311791    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000671  0.196147  0.980574  1.000000      3.3466   11.3118   11.3117    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  1.38727E+00  3.34663E+00      1   1

************************************************************************************************************************************
     10  Keyword, label(s) :  SYSTEM                            

     Number  of  commands :   1,  as  follows : 

 cat zgoubi.FITVALS.out >> zgoubi.FITVALS.out_cat

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        2 through the optical structure 

                     Total of          2 particles have been launched

 Pgm rebel. At pass #    2/   5.  In element #    1,  parameter # 35  changed to    2.13682967E+00   (was    1.38727400E+00)

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              All spins entered particle by particle
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   2.1368    -1.300     0.120     0.000     0.000       0.0000    2.1368   -1.300    0.120    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -1.3000E-02   1.2000E-04       1       0    0.000    B-Dim     1      3
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      3
    0.000        0.000        0.000        0.000        3772.           1       0    0.000    B-Dim     3      3


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                            



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000              0.000000    0.000000    1.000000    1.000000    0.000000    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000           0.000000  0.000000  1.000000  1.000000      5.0199    0.0000    0.0000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  2.13683E+00  5.01995E+00      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.I is 12. 1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           warmSnake.map                                                                   
 SBR TOSCAC,  restored mesh coordinates for field map #   1,  name : warmSnake.map


     Min/max fields in map       :  -2.934230E-01                  /   3.119730E-01
       @  X(CM),  Y(CM), Z(CM) :  -1.04     -6.000E-02  6.500E-02  /   1.05      5.500E-02  7.000E-02
      given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     Min/max normalised fields (kG)  : -2.934230E+00                      3.119730E+00

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm


     Element  is  mis-aligned  wrt.  the  optical  axis
          Center  of  entrance  EFB  is  at    X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD

  A    1  2.1368    -1.300     0.120     0.000     0.000          599.500    -1.249     0.000    -0.030     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -9.0275E-13      -2.3726E-14        1.4132E-15
                                      7.2217E-13       -4.4137E-13
                                     -1.0949E-13       -2.8051E-12
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   4.364084E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   2.1368    -1.300     0.120     0.000     0.000       0.0000    2.1368   -1.249    0.135   -0.030    0.000   4.000743E+02     1
               Time of flight (mus) :  1.36179758E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -1.2488E-02   1.3504E-04       1       0    0.000    B-Dim     1      3
    0.000        0.000        0.000      -2.9790E-04  -3.2808E-08       1       0    0.000    B-Dim     2      3
    0.000        0.000        0.000       1.3618E-02    3772.           1       0    0.000    B-Dim     3      3


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      8  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1


************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              All spins entered particle by particle
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   2.1368    -1.300     0.120     0.000     0.000       0.0000    2.1368   -1.300    0.120    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -1.3000E-02   1.2000E-04       1       0    0.000    B-Dim     1      3
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      3
    0.000        0.000        0.000        0.000        3772.           1       0    0.000    B-Dim     3      3


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                            



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000              0.000000    0.000000    1.000000    1.000000    0.000000    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000           0.000000  0.000000  1.000000  1.000000      5.0199    0.0000    0.0000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  2.13683E+00  5.01995E+00      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.I is 12. 1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           warmSnake.map                                                                   
 SBR TOSCAC,  restored mesh coordinates for field map #   1,  name : warmSnake.map


     Min/max fields in map       :  -2.934230E-01                  /   3.119730E-01
       @  X(CM),  Y(CM), Z(CM) :  -1.04     -6.000E-02  6.500E-02  /   1.05      5.500E-02  7.000E-02
      given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     Min/max normalised fields (kG)  : -2.934230E+00                      3.119730E+00

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm


     Element  is  mis-aligned  wrt.  the  optical  axis
          Center  of  entrance  EFB  is  at    X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD

  A    1  2.1368    -1.300     0.120     0.000     0.000          599.500    -1.249     0.000    -0.030     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -9.0275E-13      -2.3726E-14        1.4132E-15
                                      7.2217E-13       -4.4137E-13
                                     -1.0949E-13       -2.8051E-12
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   2.1368    -1.300     0.120     0.000     0.000       0.0000    2.1368   -1.249    0.135   -0.030    0.000   4.000743E+02     1
               Time of flight (mus) :  1.36179758E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -1.2488E-02   1.3504E-04       1       0    0.000    B-Dim     1      3
    0.000        0.000        0.000      -2.9790E-04  -3.2808E-08       1       0    0.000    B-Dim     2      3
    0.000        0.000        0.000       1.3618E-02    3772.           1       0    0.000    B-Dim     3      3


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      8  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           constraint #            1       IR =            6 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            6 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            6 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will be saved in zgoubi.FITVALS.out                                                              


 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    1    1     30   -3.00      -0.880      -0.8800000000       3.00      2.000E-02 OBJET      *          *         
    1    2     31   -3.00       0.120       0.1200000000       3.00      2.000E-02 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     6      0.0000000E+00     1.0000E+00    4.2566376E-02   4.4138E-01   TOSCA      *          *          *   0 : 
   3   1   3     6      0.0000000E+00     1.0000E+00    2.8341262E-02   1.9567E-01   TOSCA      *          *          *   0 : 
   7   1   2     6      0.0000000E+00     1.0000E+00   -3.8600183E-02   3.6296E-01   TOSCA      *          *          *   0 : 
 Fit reached penalty value   4.1051E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              All spins entered particle by particle
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   2.1368    -0.880     0.120     0.000     0.000       0.0000    2.1368   -0.880    0.120    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -8.8000E-03   1.2000E-04       1       0    0.000    B-Dim     1      3
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      3
    0.000        0.000        0.000        0.000        3772.           1       0    0.000    B-Dim     3      3


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                            



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000              0.000000    0.000000    1.000000    1.000000    0.000000    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000           0.000000  0.000000  1.000000  1.000000      5.0199    0.0000    0.0000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  2.13683E+00  5.01995E+00      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.I is 12. 1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           warmSnake.map                                                                   
 SBR TOSCAC,  restored mesh coordinates for field map #   1,  name : warmSnake.map


     Min/max fields in map       :  -2.934230E-01                  /   3.119730E-01
       @  X(CM),  Y(CM), Z(CM) :  -1.04     -6.000E-02  6.500E-02  /   1.05      5.500E-02  7.000E-02
      given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     Min/max normalised fields (kG)  : -2.934230E+00                      3.119730E+00

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm


     Element  is  mis-aligned  wrt.  the  optical  axis
          Center  of  entrance  EFB  is  at    X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD

  A    1  2.1368    -0.880     0.120     0.000     0.000          599.500    -0.837     0.000    -0.027     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -6.0593E-13      -2.3726E-14        1.0384E-16
                                      7.2217E-13       -4.4136E-13
                                     -1.0949E-13       -2.8078E-12
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    8.00000000     m ;  Time  (for ref. rigidity & particle) =   2.909389E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   2.1368    -0.880     0.120     0.000     0.000       0.0000    2.1368   -0.837    0.092   -0.027    0.000   4.000743E+02     1
               Time of flight (mus) :  1.36179768E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -8.3743E-03   9.1659E-05       1       0    0.000    B-Dim     1      3
    0.000        0.000        0.000      -2.7487E-04   1.0380E-07       1       0    0.000    B-Dim     2      3
    0.000        0.000        0.000       1.3618E-02    3772.           1       0    0.000    B-Dim     3      3


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************

      8   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
      9  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000             -0.000530    0.188579    0.982058    1.000000   10.869907    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000530  0.188579  0.982058  1.000000      5.0199   10.8699   10.8699    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  2.13683E+00  5.01995E+00      1   1

************************************************************************************************************************************
     10  Keyword, label(s) :  SYSTEM                            

     Number  of  commands :   1,  as  follows : 

 cat zgoubi.FITVALS.out >> zgoubi.FITVALS.out_cat

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        3 through the optical structure 

                     Total of          3 particles have been launched

 Pgm rebel. At pass #    3/   5.  In element #    1,  parameter # 35  changed to    4.82611907E+00   (was    2.13682967E+00)

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              All spins entered particle by particle
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   4.8261    -0.880     0.120     0.000     0.000       0.0000    4.8261   -0.880    0.120    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -8.8000E-03   1.2000E-04       1       0    0.000    B-Dim     1      4
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      4
    0.000        0.000        0.000        0.000        9529.           1       0    0.000    B-Dim     3      4


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                            



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000              0.000000    0.000000    1.000000    1.000000    0.000000    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000           0.000000  0.000000  1.000000  1.000000     11.1554    0.0000    0.0000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  4.82612E+00  1.11554E+01      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.I is 12. 1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           warmSnake.map                                                                   
 SBR TOSCAC,  restored mesh coordinates for field map #   1,  name : warmSnake.map


     Min/max fields in map       :  -2.934230E-01                  /   3.119730E-01
       @  X(CM),  Y(CM), Z(CM) :  -1.04     -6.000E-02  6.500E-02  /   1.05      5.500E-02  7.000E-02
      given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     Min/max normalised fields (kG)  : -2.934230E+00                      3.119730E+00

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm


     Element  is  mis-aligned  wrt.  the  optical  axis
          Center  of  entrance  EFB  is  at    X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD

  A    1  4.8261    -0.880     0.120     0.000     0.000          599.500    -0.832     0.000    -0.013     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -2.6796E-13      -1.0505E-14        4.5234E-17
                                      3.1975E-13       -1.9949E-13
                                     -4.8479E-14       -1.2702E-12
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   4.364084E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   4.8261    -0.880     0.120     0.000     0.000       0.0000    4.8261   -0.832    0.118   -0.013    0.000   4.000145E+02     1
               Time of flight (mus) :  1.33969849E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -8.3234E-03   1.1788E-04       1       0    0.000    B-Dim     1      4
    0.000        0.000        0.000      -1.3446E-04  -4.6479E-07       1       0    0.000    B-Dim     2      4
    0.000        0.000        0.000       1.3397E-02    9529.           1       0    0.000    B-Dim     3      4


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      8  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1


************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              All spins entered particle by particle
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   4.8261    -0.880     0.120     0.000     0.000       0.0000    4.8261   -0.880    0.120    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -8.8000E-03   1.2000E-04       1       0    0.000    B-Dim     1      4
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      4
    0.000        0.000        0.000        0.000        9529.           1       0    0.000    B-Dim     3      4


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                            



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000              0.000000    0.000000    1.000000    1.000000    0.000000    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000           0.000000  0.000000  1.000000  1.000000     11.1554    0.0000    0.0000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  4.82612E+00  1.11554E+01      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.I is 12. 1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           warmSnake.map                                                                   
 SBR TOSCAC,  restored mesh coordinates for field map #   1,  name : warmSnake.map


     Min/max fields in map       :  -2.934230E-01                  /   3.119730E-01
       @  X(CM),  Y(CM), Z(CM) :  -1.04     -6.000E-02  6.500E-02  /   1.05      5.500E-02  7.000E-02
      given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     Min/max normalised fields (kG)  : -2.934230E+00                      3.119730E+00

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm


     Element  is  mis-aligned  wrt.  the  optical  axis
          Center  of  entrance  EFB  is  at    X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD

  A    1  4.8261    -0.880     0.120     0.000     0.000          599.500    -0.832     0.000    -0.013     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -2.6796E-13      -1.0505E-14        4.5234E-17
                                      3.1975E-13       -1.9949E-13
                                     -4.8479E-14       -1.2702E-12
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   4.8261    -0.880     0.120     0.000     0.000       0.0000    4.8261   -0.832    0.118   -0.013    0.000   4.000145E+02     1
               Time of flight (mus) :  1.33969849E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -8.3234E-03   1.1788E-04       1       0    0.000    B-Dim     1      4
    0.000        0.000        0.000      -1.3446E-04  -4.6479E-07       1       0    0.000    B-Dim     2      4
    0.000        0.000        0.000       1.3397E-02    9529.           1       0    0.000    B-Dim     3      4


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      8  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           constraint #            1       IR =            6 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            6 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            6 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will be saved in zgoubi.FITVALS.out                                                              


 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    1    1     30   -3.00      -0.400      -0.4000000000       3.00      2.000E-02 OBJET      *          *         
    1    2     31   -3.00       0.120       0.1200000000       3.00      2.000E-02 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     6      0.0000000E+00     1.0000E+00    4.5659894E-02   7.5534E-01   TOSCA      *          *          *   0 : 
   3   1   3     6      0.0000000E+00     1.0000E+00    1.2231785E-02   5.4207E-02   TOSCA      *          *          *   0 : 
   7   1   2     6      0.0000000E+00     1.0000E+00   -2.2927628E-02   1.9045E-01   TOSCA      *          *          *   0 : 
 Fit reached penalty value   2.7601E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              All spins entered particle by particle
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   4.8261    -0.400     0.120     0.000     0.000       0.0000    4.8261   -0.400    0.120    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -4.0000E-03   1.2000E-04       1       0    0.000    B-Dim     1      4
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      4
    0.000        0.000        0.000        0.000        9529.           1       0    0.000    B-Dim     3      4


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                            



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000              0.000000    0.000000    1.000000    1.000000    0.000000    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000           0.000000  0.000000  1.000000  1.000000     11.1554    0.0000    0.0000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  4.82612E+00  1.11554E+01      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.I is 12. 1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           warmSnake.map                                                                   
 SBR TOSCAC,  restored mesh coordinates for field map #   1,  name : warmSnake.map


     Min/max fields in map       :  -2.934230E-01                  /   3.119730E-01
       @  X(CM),  Y(CM), Z(CM) :  -1.04     -6.000E-02  6.500E-02  /   1.05      5.500E-02  7.000E-02
      given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     Min/max normalised fields (kG)  : -2.934230E+00                      3.119730E+00

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm


     Element  is  mis-aligned  wrt.  the  optical  axis
          Center  of  entrance  EFB  is  at    X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD

  A    1  4.8261    -0.400     0.120     0.000     0.000          599.500    -0.354     0.000    -0.013     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -1.0848E-13      -5.5780E-15       -5.7161E-16
                                      3.2984E-13       -1.9945E-13
                                     -5.3047E-14       -1.2710E-12
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    8.00000000     m ;  Time  (for ref. rigidity & particle) =   2.909389E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1   4.8261    -0.400     0.120     0.000     0.000       0.0000    4.8261   -0.354    0.108   -0.013    0.000   4.000145E+02     1
               Time of flight (mus) :  1.33969849E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -3.5434E-03   1.0777E-04       1       0    0.000    B-Dim     1      4
    0.000        0.000        0.000      -1.2806E-04  -3.7453E-07       1       0    0.000    B-Dim     2      4
    0.000        0.000        0.000       1.3397E-02    9529.           1       0    0.000    B-Dim     3      4


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************

      8   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
      9  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000             -0.000434    0.184148    0.982898    1.000000   10.611513    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000434  0.184148  0.982898  1.000000     11.1554   10.6115   10.6115    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  4.82612E+00  1.11554E+01      1   1

************************************************************************************************************************************
     10  Keyword, label(s) :  SYSTEM                            

     Number  of  commands :   1,  as  follows : 

 cat zgoubi.FITVALS.out >> zgoubi.FITVALS.out_cat

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        4 through the optical structure 

                     Total of          4 particles have been launched


      Next  pass  is  #     5 and  last  pass  through  the  optical  structure


 Pgm rebel. At pass #    4/   5.  In element #    1,  parameter # 35  changed to    1.10152413E+01   (was    4.82611907E+00)

     WRITE statements to zgoubi.res are re-established from now on.

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              All spins entered particle by particle
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1  11.0152    -0.400     0.120     0.000     0.000       0.0000   11.0152   -0.400    0.120    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -4.0000E-03   1.2000E-04       1       0    0.000    B-Dim     1      5
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      5
    0.000        0.000        0.000        0.000       2.2874E+04       1       0    0.000    B-Dim     3      5


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                            



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000              0.000000    0.000000    1.000000    1.000000    0.000000    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000           0.000000  0.000000  1.000000  1.000000     25.3786    0.0000    0.0000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  1.10152E+01  2.53786E+01      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.I is 12. 1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           warmSnake.map                                                                   
 SBR TOSCAC,  restored mesh coordinates for field map #   1,  name : warmSnake.map


     Min/max fields in map       :  -2.934230E-01                  /   3.119730E-01
       @  X(CM),  Y(CM), Z(CM) :  -1.04     -6.000E-02  6.500E-02  /   1.05      5.500E-02  7.000E-02
      given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     Min/max normalised fields (kG)  : -2.934230E+00                      3.119730E+00

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm


     Element  is  mis-aligned  wrt.  the  optical  axis
          Center  of  entrance  EFB  is  at    X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD

  A    1 11.0152    -0.400     0.120     0.000     0.000          599.500    -0.353     0.000    -0.006     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -4.7510E-14      -2.4439E-15       -2.3047E-16
                                      1.4452E-13       -8.7736E-14
                                     -2.3242E-14       -5.5918E-13
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   4.364084E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1  11.0152    -0.400     0.120     0.000     0.000       0.0000   11.0152   -0.353    0.115   -0.006    0.000   4.000028E+02     1
               Time of flight (mus) :  1.33530270E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -3.5287E-03   1.1542E-04       1       0    0.000    B-Dim     1      5
    0.000        0.000        0.000      -5.7228E-05  -2.2002E-07       1       0    0.000    B-Dim     2      5
    0.000        0.000        0.000       1.3353E-02   2.2874E+04       1       0    0.000    B-Dim     3      5


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      8  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1


************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              All spins entered particle by particle
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1  11.0152    -0.400     0.120     0.000     0.000       0.0000   11.0152   -0.400    0.120    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -4.0000E-03   1.2000E-04       1       0    0.000    B-Dim     1      5
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      5
    0.000        0.000        0.000        0.000       2.2874E+04       1       0    0.000    B-Dim     3      5


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                            



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000              0.000000    0.000000    1.000000    1.000000    0.000000    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000           0.000000  0.000000  1.000000  1.000000     25.3786    0.0000    0.0000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  1.10152E+01  2.53786E+01      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.I is 12. 1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           warmSnake.map                                                                   
 SBR TOSCAC,  restored mesh coordinates for field map #   1,  name : warmSnake.map


     Min/max fields in map       :  -2.934230E-01                  /   3.119730E-01
       @  X(CM),  Y(CM), Z(CM) :  -1.04     -6.000E-02  6.500E-02  /   1.05      5.500E-02  7.000E-02
      given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     Min/max normalised fields (kG)  : -2.934230E+00                      3.119730E+00

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm


     Element  is  mis-aligned  wrt.  the  optical  axis
          Center  of  entrance  EFB  is  at    X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD

  A    1 11.0152    -0.400     0.120     0.000     0.000          599.500    -0.353     0.000    -0.006     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -4.7510E-14      -2.4439E-15       -2.3047E-16
                                      1.4452E-13       -8.7736E-14
                                     -2.3242E-14       -5.5918E-13
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1  11.0152    -0.400     0.120     0.000     0.000       0.0000   11.0152   -0.353    0.115   -0.006    0.000   4.000028E+02     1
               Time of flight (mus) :  1.33530270E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -3.5287E-03   1.1542E-04       1       0    0.000    B-Dim     1      5
    0.000        0.000        0.000      -5.7228E-05  -2.2002E-07       1       0    0.000    B-Dim     2      5
    0.000        0.000        0.000       1.3353E-02   2.2874E+04       1       0    0.000    B-Dim     3      5


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      8  Keyword, label(s) :  FIT                               

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           constraint #            1       IR =            6 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            6 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            6 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will be saved in zgoubi.FITVALS.out                                                              


 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    1    1     30   -3.00      -0.220      -0.2200000000       3.00      2.000E-02 OBJET      *          *         
    1    2     31   -3.00       0.120       0.1200000000       3.00      2.000E-02 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     6      0.0000000E+00     1.0000E+00    4.7047888E-02   2.2979E-01   TOSCA      *          *          *   0 : 
   3   1   3     6      0.0000000E+00     1.0000E+00    4.9788712E-03   2.5734E-03   TOSCA      *          *          *   0 : 
   7   1   2     6      0.0000000E+00     1.0000E+00   -8.5991763E-02   7.6764E-01   TOSCA      *          *          *   0 : 
 Fit reached penalty value   9.6329E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing a last run using variable values from FIT. 

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =       7205.178 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


************************************************************************************************************************************
      3  Keyword, label(s) :  SPNTRK                            


                SPIN  TRACKING  REQUESTED  

                          PARTICLE  MASS          =    938.2720     MeV/c2
                          GYROMAGNETIC  FACTOR  G =    1.792847    

                          INITIAL SPIN CONDITIONS TYPE  4 :
                              All spins entered particle by particle
                              Particles # 1 to      10 may be subjected to spin matching using FIT procedure


                          PARAMETRES  DYNAMIQUES  DE  REFERENCE :
                               BORO   =      7205.178 KG*CM
                               BETA   =    0.917207
                               GAMMA*G =   4.500000


                          POLARISATION  INITIALE  MOYENNE  DU  FAISCEAU  DE        1  PARTICULES :
                               <SX> =     0.0000
                               <SY> =     0.0000
                               <SZ> =     1.0000
                               <S>  =     1.0000

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1  11.0152    -0.220     0.120     0.000     0.000       0.0000   11.0152   -0.220    0.120    0.000    0.000   0.000000E+00     1
               Time of flight (mus) :   0.0000000     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -2.2000E-03   1.2000E-04       1       0    0.000    B-Dim     1      5
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      5
    0.000        0.000        0.000        0.000       2.2874E+04       1       0    0.000    B-Dim     3      5


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      5  Keyword, label(s) :  SPNPRT                            



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000              0.000000    0.000000    1.000000    1.000000    0.000000    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000           0.000000  0.000000  1.000000  1.000000     25.3786    0.0000    0.0000    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  1.10152E+01  2.53786E+01      1   1

************************************************************************************************************************************
      6  Keyword, label(s) :  TOSCA                             


     zgoubi.plt                                                                      
      already open...

     NDIM =   3 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.I is 12. 1

          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           warmSnake.map                                                                   
 SBR TOSCAC,  restored mesh coordinates for field map #   1,  name : warmSnake.map


     Min/max fields in map       :  -2.934230E-01                  /   3.119730E-01
       @  X(CM),  Y(CM), Z(CM) :  -1.04     -6.000E-02  6.500E-02  /   1.05      5.500E-02  7.000E-02
      given normalisation coeffs on field, x, y, z :   1.000000E+01   1.000000E+02   1.000000E+02   1.000000E+02
     Min/max normalised fields (kG)  : -2.934230E+00                      3.119730E+00

     Length of element,  XL =  4.000000E+02 cm 
                                               from  XI =   1.995000E+02 cm 
                                               to    XF =   5.995000E+02 cm 

     Nbr of nodes in X = 801;  nbr of nodes in Y =   29
     X-size of mesh =  5.000000E-01 cm ; Y-size =  5.000000E-01 cm
     nber of nodes in Z =   57 ; Step in Z =  0.500000E+00 cm

                     OPTION  DE  CALCUL  :  GRILLE  3-D  A  3*3*3  POINTS , INTERPOLATION  A  L'ORDRE  2

                    Integration step :  0.1000     cm


     Element  is  mis-aligned  wrt.  the  optical  axis
          Center  of  entrance  EFB  is  at    X =   0.000     CM   Y =   0.000     cm,  tilt  angle =   0.00000     RAD

  A    1 11.0152    -0.220     0.120     0.000     0.000          599.500    -0.173     0.000    -0.006     0.000            1


                CONDITIONS  DE  MAXWELL  (     4001.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                    -2.1972E-14       4.3083E-17        1.7902E-16
                                      1.2338E-13       -8.7761E-14
                                     -1.9032E-14       -5.5921E-13
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    8.00000000     m ;  Time  (for ref. rigidity & particle) =   2.909389E-08 s 

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISCEAU                          

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

o  1  11.0152    -0.220     0.120     0.000     0.000       0.0000   11.0152   -0.173    0.115   -0.006    0.000   4.000028E+02     1
               Time of flight (mus) :  1.33530271E-02 mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000      -1.7295E-03   1.1502E-04       1       0    0.000    B-Dim     1      5
    0.000        0.000        0.000      -5.6525E-05  -2.1974E-07       1       0    0.000    B-Dim     2      5
    0.000        0.000        0.000       1.3353E-02   2.2874E+04       1       0    0.000    B-Dim     3      5


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN           NaN
           AlpZ, BetZ =           NaN           NaN

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************

      8   Keyword FIT[2] is skipped since this is the (end of) last run following the fitting procedure.

          Now carrying on beyond FIT keyword.


************************************************************************************************************************************
      9  Keyword, label(s) :  SPNPRT      PRINT                 



                          Average  over  particles at this pass ;   beam with       1  particles :

                   INITIAL                                           FINAL

        <SX>        <SY>        <SZ>        <S>                   <SX>       <SY>         <SZ>         <S>  <(SI,SF)>  sigma_(SI,SF)
                                                                                                               (deg)       (deg)
    0.000000    0.000000    1.000000    1.000000             -0.000419    0.183281    0.983060    1.000000   10.560969    0.000000


                Spin  components  of  each  of  the      1  particles,  and  rotation  angle :

                   INITIAL                                           FINAL

           SX        SY        SZ        |S|               SX        SY        SZ        |S|        GAMMA    (Si,Sf)   (Si,Sf_x)
                                                                                                              (deg.)     (deg.)
                                                                                           (Sf_x : projection of Sf on plane x=0)

 o  1  0.000000  0.000000  1.000000  1.000000          -0.000419  0.183281  0.983060  1.000000     25.3786   10.5610   10.5609    1



                Min/Max  components  of  each  of  the      1  particles :

  SX_mi       SX_ma       SY_mi       SY_ma       SZ_mi       SZ_ma       |S|_mi      |S|_ma      p/p_0        GAMMA          I  IEX

 -9.6557E-04  0.0000E+00  0.0000E+00  2.0832E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.0000E+00  1.10152E+01  2.53786E+01      1   1

************************************************************************************************************************************
     10  Keyword, label(s) :  SYSTEM                            

     Number  of  commands :   1,  as  follows : 

 cat zgoubi.FITVALS.out >> zgoubi.FITVALS.out_cat

************************************************************************************************************************************
     11  Keyword, label(s) :  REBELOTE                          


                         ****  End  of  'REBELOTE'  procedure  ****

      There  has  been          5  passes  through  the  optical  structure 

                     Total of          5 particles have been launched

 That's where this REBELOTE option ends,  for the time being.

End of job !

  
