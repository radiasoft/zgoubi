Centering 5 helical orbits in the AGS warm helical snake 3-D OPERA map.                                       
 'OBJET'                                                                                                      1
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
                                                                                                              
 'TOSCA'                                                                                                      4
0  2                                                                                                          
1.e1   100. 100. 100.                                                                                         
HEADER_0 wsnake                                                                                               
801 29 29 12.1                                                                                                
warmSnake.map                                                                                                 
0 0 0 0                                                                                                       
2                                                                                                             
.1                                                                                                            
2  0.  .0  0.  0.                                                                                             
                                                                                                              
 'FIT'                                                                                                        5
2                                                                                                             
1 30 0  [-3,3]                                                                                                
1 31 0  [-3,3]                                                                                                
3    1E-2                                                                                                     
3.1 1 2 4 0. 1. 0                                                                                             
3.1 1 3 4 0. 1. 0                                                                                             
7.3 1 2 4 0. 1. 0                                                                                             
                                                                                                              
 'SPNPRT'  PRINT                                                                                              6
                                                                                                              
 'REBELOTE'                                                                                                   7
4  0.1  0 1                                                                                                   
1                                                                                                             
OBJET  35  1.3872739973  2.1368296674  4.8261190694  11.015241321                                             
                                                                                                              
 'END'                                                                                                        8

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
      4  Keyword, label(s) :  TOSCA                             


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
      5  Keyword, label(s) :  FIT                               

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            4 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save' command

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    1    1     30   -3.00       -1.84       -1.840000000       3.00      2.000E-02 OBJET      *          *         
    1    2     31   -3.00       6.000E-02   6.0000000000E-02   3.00      2.000E-02 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     4      0.0000000E+00     1.0000E+00    8.1418231E-03   7.2565E-03   TOSCA      *          *          *   0 : 
   3   1   3     4      0.0000000E+00     1.0000E+00    8.1629180E-02   7.2942E-01   TOSCA      *          *          *   0 : 
   7   1   2     4      0.0000000E+00     1.0000E+00   -4.9046288E-02   2.6333E-01   TOSCA      *          *          *   0 : 
 Fit reached penalty value   9.1352E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing final run using FIT variable values. 

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
      4  Keyword, label(s) :  TOSCA                             


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

      5   Keyword FIT is skipped since this is the (end of) final run following the fitting procedure.

          Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  SPNPRT      PRINT                 



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

 -9.6557E-04 -9.6557E-04  2.0832E-01  2.0832E-01  9.7806E-01  9.7806E-01  1.0000E+00  1.0000E+00  1.00000E+00  2.50997E+00      1   1

************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                          


     Multiple pass, 
          from element #     1 : OBJET     /label1=          /label2=           to REBELOTE /label1=          /label2=          
          ending at pass #       5 at element #     7 : REBELOTE  /label1=          /label2=          


     Parameter #   35 in element #    1 will be modified at each pass, 
     list of requested parameter values :
                   1   1.38727400E+00
                   2   2.13682967E+00
                   3   4.82611907E+00
                   4   1.10152413E+01

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.840     0.060     0.000     0.000          599.500    -1.801     0.000    -0.046     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.780     0.060     0.000     0.000          599.500    -1.744     0.000    -0.045     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.780     0.120     0.000     0.000          599.500    -1.720     0.000    -0.044     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.720     0.180     0.000     0.000          599.500    -1.640     0.000    -0.043     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.660     0.240     0.000     0.000          599.500    -1.559     0.000    -0.042     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.600     0.300     0.000     0.000          599.500    -1.479     0.000    -0.041     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.540     0.360     0.000     0.000          599.500    -1.398     0.000    -0.040     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.480     0.420     0.000     0.000          599.500    -1.317     0.000    -0.039     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.420     0.480     0.000     0.000          599.500    -1.237     0.000    -0.038     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.360     0.540     0.000     0.000          599.500    -1.156     0.000    -0.037     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.300     0.600     0.000     0.000          599.500    -1.076     0.001    -0.036     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.360     0.540     0.000     0.000          599.500    -1.156     0.000    -0.037     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.300     0.540     0.000     0.000          599.500    -1.099     0.000    -0.037     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.420     0.540     0.000     0.000          599.500    -1.213     0.000    -0.038     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.360     0.600     0.000     0.000          599.500    -1.133     0.001    -0.037     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.360     0.480     0.000     0.000          599.500    -1.180     0.000    -0.038     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.360     0.420     0.000     0.000          599.500    -1.203     0.000    -0.038     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.360     0.360     0.000     0.000          599.500    -1.227     0.000    -0.038     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.360     0.300     0.000     0.000          599.500    -1.250     0.000    -0.038     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.360     0.240     0.000     0.000          599.500    -1.274     0.000    -0.038     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.360     0.180     0.000     0.000          599.500    -1.297     0.000    -0.038     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.360     0.120     0.000     0.000          599.500    -1.321     0.000    -0.038     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.360     0.180     0.000     0.000          599.500    -1.297     0.000    -0.038     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.300     0.180     0.000     0.000          599.500    -1.240     0.000    -0.038     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.300     0.240     0.000     0.000          599.500    -1.217     0.000    -0.037     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

  A    1  1.3873    -1.240     0.060     0.000     0.000          599.500    -1.230     0.000    -0.037     0.000            1

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

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
      4  Keyword, label(s) :  TOSCA                             


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

 Cumulative length of optical axis =    4.00000000     m ;  Time  (for ref. rigidity & particle) =   1.454695E-08 s 

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    1    1     30   -3.00       -1.30       -1.300000000       3.00      2.000E-02 OBJET      *          *         
    1    2     31   -3.00       0.120       0.1200000000       3.00      2.000E-02 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     4      0.0000000E+00     1.0000E+00    3.6561969E-02   2.2023E-01   TOSCA      *          *          *   0 : 
   3   1   3     4      0.0000000E+00     1.0000E+00    5.9165583E-02   5.7671E-01   TOSCA      *          *          *   0 : 
   7   1   2     4      0.0000000E+00     1.0000E+00    3.5107390E-02   2.0306E-01   TOSCA      *          *          *   0 : 
 Fit reached penalty value   6.0699E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing final run using FIT variable values. 

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
      4  Keyword, label(s) :  TOSCA                             


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

      5   Keyword FIT is skipped since this is the (end of) final run following the fitting procedure.

          Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  SPNPRT      PRINT                 



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

 -9.6557E-04 -6.7114E-04  1.9615E-01  2.0832E-01  9.7806E-01  9.8057E-01  1.0000E+00  1.0000E+00  1.38727E+00  3.34663E+00      1   1

************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                          


 SBR rebel. At pass #    2.  In element #    1,  changed value of parameter # 35  to :   2.13682967E+00   (was :   1.38727400E+00)

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
      4  Keyword, label(s) :  TOSCA                             


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
      5  Keyword, label(s) :  FIT                               


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
      4  Keyword, label(s) :  TOSCA                             


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
      5  Keyword, label(s) :  FIT                               

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            4 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save' command

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    1    1     30   -3.00      -0.880      -0.8800000000       3.00      2.000E-02 OBJET      *          *         
    1    2     31   -3.00       0.120       0.1200000000       3.00      2.000E-02 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     4      0.0000000E+00     1.0000E+00    4.2566376E-02   4.4138E-01   TOSCA      *          *          *   0 : 
   3   1   3     4      0.0000000E+00     1.0000E+00    2.8341262E-02   1.9567E-01   TOSCA      *          *          *   0 : 
   7   1   2     4      0.0000000E+00     1.0000E+00   -3.8600183E-02   3.6296E-01   TOSCA      *          *          *   0 : 
 Fit reached penalty value   4.1051E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing final run using FIT variable values. 

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
      4  Keyword, label(s) :  TOSCA                             


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

      5   Keyword FIT is skipped since this is the (end of) final run following the fitting procedure.

          Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  SPNPRT      PRINT                 



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

 -9.6557E-04 -5.2989E-04  1.8858E-01  2.0832E-01  9.7806E-01  9.8206E-01  1.0000E+00  1.0000E+00  2.13683E+00  5.01995E+00      1   1

************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                          


 SBR rebel. At pass #    3.  In element #    1,  changed value of parameter # 35  to :   4.82611907E+00   (was :   2.13682967E+00)

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
      4  Keyword, label(s) :  TOSCA                             


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
      5  Keyword, label(s) :  FIT                               


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
      4  Keyword, label(s) :  TOSCA                             


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
      5  Keyword, label(s) :  FIT                               

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            4 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save' command

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    1    1     30   -3.00      -0.400      -0.4000000000       3.00      2.000E-02 OBJET      *          *         
    1    2     31   -3.00       0.120       0.1200000000       3.00      2.000E-02 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     4      0.0000000E+00     1.0000E+00    4.5659894E-02   7.5534E-01   TOSCA      *          *          *   0 : 
   3   1   3     4      0.0000000E+00     1.0000E+00    1.2231785E-02   5.4207E-02   TOSCA      *          *          *   0 : 
   7   1   2     4      0.0000000E+00     1.0000E+00   -2.2927628E-02   1.9045E-01   TOSCA      *          *          *   0 : 
 Fit reached penalty value   2.7601E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing final run using FIT variable values. 

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
      4  Keyword, label(s) :  TOSCA                             


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

      5   Keyword FIT is skipped since this is the (end of) final run following the fitting procedure.

          Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  SPNPRT      PRINT                 



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

 -9.6557E-04 -4.3378E-04  1.8415E-01  2.0832E-01  9.7806E-01  9.8290E-01  1.0000E+00  1.0000E+00  4.82612E+00  1.11554E+01      1   1

************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                          


                                -----  REBELOTE  -----

     End of pass #        4 through the optical structure 

                     Total of          4 particles have been launched

********************************************************************************************************************************

********************************************************************************************************************************

********************************************************************************************************************************


      Next  pass  is  #     5 and  last  pass  through  the  optical  structure


 SBR rebel. At pass #    4.  In element #    1,  changed value of parameter # 35  to :   1.10152413E+01

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
      4  Keyword, label(s) :  TOSCA                             


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
      5  Keyword, label(s) :  FIT                               


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
      4  Keyword, label(s) :  TOSCA                             


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
      5  Keyword, label(s) :  FIT                               

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            1 ,   ok.
           variable #            2       IP =           31 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.
           constraint #            3       IR =            4 ,   ok.

                    FIT  variables  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save' command

 STATUS OF VARIABLES  (Iteration #     0 /     90 max.)
 LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL         MAXIMUM      STEP     NAME       LBL1     LBL2
    1    1     30   -3.00      -0.220      -0.2200000000       3.00      2.000E-02 OBJET      *          *         
    1    2     31   -3.00       0.120       0.1200000000       3.00      2.000E-02 OBJET      *          *         
 STATUS OF CONSTRAINTS (Target penalty =   1.0000E-02)
 TYPE  I   J  LMNT#       DESIRED           WEIGHT         REACHED         KI2       NAME       LBL1     LBL2         *  Parameter(s) 
   3   1   2     4      0.0000000E+00     1.0000E+00    4.7047888E-02   2.2979E-01   TOSCA      *          *          *   0 : 
   3   1   3     4      0.0000000E+00     1.0000E+00    4.9788712E-03   2.5734E-03   TOSCA      *          *          *   0 : 
   7   1   2     4      0.0000000E+00     1.0000E+00   -8.5991763E-02   7.6764E-01   TOSCA      *          *          *   0 : 
 Fit reached penalty value   9.6329E-03

************************************************************************************************************************************

           MAIN PROGRAM :  FIT completed.  Now doing final run using FIT variable values. 

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
      4  Keyword, label(s) :  TOSCA                             


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

      5   Keyword FIT is skipped since this is the (end of) final run following the fitting procedure.

          Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  SPNPRT      PRINT                 



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

 -9.6557E-04 -4.1869E-04  1.8328E-01  2.0832E-01  9.7806E-01  9.8306E-01  1.0000E+00  1.0000E+00  1.10152E+01  2.53786E+01      1   1

************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                          


                         ****  End  of  'REBELOTE'  procedure  ****

      There  has  been          5  passes  through  the  optical  structure 

