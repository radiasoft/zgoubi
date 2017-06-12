Uniform field sector
 'OBJET'                                                                                                      1
64.62444403717985                                   ! 200keV
2
4 1
12.9248888074 0.  0.   0.   0.  1. 'm'              ! 200keV. R=Brho/B=*/.5
28.9070891209 0.  0.   0.   0.  2.23654451125 'm'   ! 1 MeV. R=Brho/B=*/.5
50.           0.  0.   0.   0.  3.86850523397 'o'   ! at RM  (B*rho=0.5*0.5=0.25T.m, 2.9885 MeV)
64.7070336799 0.  0.   0.   0.  5.0063899693  'M'   ! 5 MeV. R=Brho/B=*/.5
1 1 1 1
 
 'PARTICUL'   ! This is required only because we want to get the time-of-flight,                              2
PROTON       ! otherwise zgoubi only requires rigidity.
 
 'TOSCA'                                                                                                      3
0  0
1.  1. 1. 1.
HEADER_8
106 121 1 22.1 1.       IZ=1 -> 2D ; MOD=22 -> polar map ; .MOD2=.1 -> one map file
geneSectorMap.out
0 0 0 0
2
.1
2
0. 0. 0. 0.
 
 'FAISCEAU'                                                                                                   4
 
 'FIT'                                                                                                        5
2   nofinal
1 30 0 [12.,65.]     ! Variable : Y_0
3 50 0 .5            ! Variable : step size in DIPOLE
2  2e-6  199         ! Penalty; max numb of calls to function
3.1 1 2 #End 0. 1. 0    ! Constraint :  Y_final=Y_0
3.1 1 3 #End 0. 1. 0    ! Constraint :  T_final=T_0
 
 'FAISTORE'                                                                                                   6
zgoubi.fai
1
 
 'REBELOTE'                                                                                                   7
20 0.2  0 1        ! 20 different rigidities ;  ; coordinates as found in OBJET ; change parameter(s)
1
OBJET 35     1:5.0063899693    ! 0.2 MeV to 5 MeV
 
 'END'                                                                                                        8

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       4 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


     Particle  properties :
     PROTON
                     Mass          =    938.272        MeV/c2
                     Charge        =   1.602176E-19    C     
                     G  factor     =    1.79285              
                     COM life-time =   1.000000E+99    s     

              Reference  data :
                    mag. rigidity (kG.cm)   :   64.624444      =p/q, such that dev.=B*L/rigidity
                    mass (MeV/c2)           :   938.27203    
                    momentum (MeV/c)        :   19.373921    
                    energy, total (MeV)     :   938.47203    
                    energy, kinetic (MeV)   :  0.20000000    
                    beta = v/c              :  2.0644111178E-02
                    gamma                   :   1.000213158    
                    beta*gamma              :  2.0648511631E-02
                    G*gamma                 :   1.793229509    
                    electric rigidity (MeV) :  0.3999573775      =T[eV]*(gamma+1)/gamma, such that dev.=E*L/rigidity
  
 I, AMQ(1,I), AMQ(2,I)/QE, P/Pref, v/c, time, s :
  
     1   9.38272030E+02  1.00000000E+00  1.00000000E+00  2.06441112E-02  0.00000000E+00  0.00000000E+00
     2   9.38272030E+02  1.00000000E+00  2.23654451E+00  4.61321482E-02  0.00000000E+00  0.00000000E+00
     3   9.38272030E+02  1.00000000E+00  3.86850523E+00  7.96252494E-02  0.00000000E+00  0.00000000E+00
     4   9.38272030E+02  1.00000000E+00  5.00638997E+00  1.02826545E-01  0.00000000E+00  0.00000000E+00

************************************************************************************************************************************
      3  Keyword, label(s) :  TOSCA                                                 


     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00

          New field map(s) now used, polar mesh (MOD .ge. 20) ;  name(s) of map data file(s) are : 
          geneSectorMap.out

   ----
   Map file number    1 ( of 1) :

     geneSectorMap.out map,  FORMAT type : regular.  Field multiplication factor :  1.00000000E+00

 HEADER  (8 lines) : 
        10.000000000000000       0.50000000000000000       0.57142857142857140      -4.97067035767729219E-042      !  Rmi/cm,
      # Field map generated using geneSectorMap.f                                                                            
     # AT/rd,  AT/deg, Rmi/cm, Rma/cm, RM/cm, NR, dR/cm, NX, dX/cm, dA/rd :                                                  
     #   1.04719755E+00   6.00000000E+01   1.00000000E+01   7.00000000E+01   5.00000000E+01 121   5.00000000E-01 106   4.9866
      #                                                                                                                      
      # R, Z, X, BY, BZ, BX                                                                                                  
      # cm cm rd kG  kG  kG                                                                                                  
      #                                                                                                                      
 
R_min (cm), DR (cm), DTTA (deg), DZ (cm) :   1.000000E+01  5.000000E-01  5.714286E-01 -4.970670E-42



     Field map limits, angle :  min, max, max-min (rad) :  -0.52858543E+00   0.51861212E+00   0.10471976E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+02   0.70000000E+02   0.60000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      10.0     -0.00      /   0.00      10.0     -0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  106/ 121/   1
     Node  distance  in   X/Y/Z :   9.973310E-03/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :  0.1000     cm   (i.e.,   2.5000E-03 rad  at mean radius RM =    40.00    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            0.519    12.925    -0.000     0.000     0.000            1
  A    1  2.2365    28.907     0.000     0.000     0.000            0.519    28.907    -0.000     0.000     0.000            2
  A    1  3.8685    50.000     0.000     0.000     0.000            0.519    50.000     0.000     0.000     0.000            3
  A    1  5.0064    64.707     0.000     0.000     0.000            0.519    64.707    -0.000     0.000     0.000            4

 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  4 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    0.0000   12.925   -0.000    0.000    0.000   1.353491E+01     1
               Time of flight (mus) :  2.18694843E-02 mass (MeV/c2) :   938.272    
m  1   2.2365    28.907     0.000     0.000     0.000       0.0000    1.2365   28.907   -0.000    0.000    0.000   3.027143E+01     2
               Time of flight (mus) :  2.18881269E-02 mass (MeV/c2) :   938.272    
o  1   3.8685    50.000     0.000     0.000     0.000       0.0000    2.8685   50.000    0.000    0.000    0.000   5.235988E+01     3
               Time of flight (mus) :  2.19344684E-02 mass (MeV/c2) :   938.272    
M  1   5.0064    64.707     0.000     0.000     0.000       0.0000    4.0064   64.707   -0.000    0.000    0.000   6.776105E+01     4
               Time of flight (mus) :  2.19813400E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   2.0205E-11  -1.3291E+00   6.0779E+09   3.913475E-01  -3.417563E-11        4        4    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        4        4    1.000      (Z,P)         1
   9.3007E-13  -2.7283E+08   6.3577E+03   2.191835E-02   2.297156E+00        4        4    1.000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   1.977092E-01
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   5.410451E-11

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   4.338446E-05
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   1.861731E+00


  Beam  sigma  matrix : 

   3.908893E-02   8.547660E-12   0.000000E+00   0.000000E+00
   8.547660E-12   2.927298E-21   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     6.431366E-12    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      5  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            3 ,   ok.
           variable #            2       IP =           50 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.

                    FIT  variables  and  constraints  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        12.9       12.924889       65.0      0.177      OBJET      -                    -                   
   3   2    50   5.000E-02   0.100      0.10000000      0.150      3.333E-04  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    9.425420E-10    5.45E-05 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    1.277014E-07    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.6309E-14

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


                OPEN FILE zgoubi.fai                                                                      
                FOR PRINTING COORDINATES 

               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        1 through the optical structure 

                     Total of          4 particles have been launched

     Multiple pass, 
          from element #     1 : OBJET     /label1=                    /label2=                    
                             to  REBELOTE  /label1=                    /label2=                    
     ending at pass #      21 at element #     7 : REBELOTE  /label1=                    /label2=                    


     Parameter #   35 in element #    1 will be modified at each pass. 
     list of requested parameter values :
                   1   1.00000000E+00
                   2   1.21086263E+00
                   3   1.42172526E+00
                   4   1.63258789E+00
                   5   1.84345052E+00
                   6   2.05431315E+00
                   7   2.26517578E+00
                   8   2.47603841E+00
                   9   2.68690104E+00
                  10   2.89776367E+00
                  11   3.10862630E+00
                  12   3.31948893E+00
                  13   3.53035156E+00
                  14   3.74121419E+00
                  15   3.95207682E+00
                  16   4.16293945E+00
                  17   4.37380208E+00
                  18   4.58466471E+00
                  19   4.79552734E+00
                  20   5.00638997E+00

 Pgm rebel. At pass #    1/  21.  In element #    1,  parameter # 35  changed to    1.00000000E+00   (was    1.00000000E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        12.9       12.924889       65.0      0.177      OBJET      -                    -                   
   3   2    50   5.000E-02   0.100      0.10000000      0.150      3.333E-04  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    9.425420E-10    5.45E-05 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    1.277014E-07    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.6309E-14

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        2 through the optical structure 

                     Total of          8 particles have been launched

 Pgm rebel. At pass #    2/  21.  In element #    1,  parameter # 35  changed to    1.21086263E+00   (was    1.00000000E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        15.7       15.650257       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14982213      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    3.978185E-06    8.16E-05 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    4.402741E-04    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.9386E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        3 through the optical structure 

                     Total of         12 particles have been launched

 Pgm rebel. At pass #    3/  21.  In element #    1,  parameter # 35  changed to    1.42172526E+00   (was    1.21086263E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        18.4       18.375625       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14980018      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    7.960114E-06    1.13E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    7.503044E-04    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   5.6302E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        4 through the optical structure 

                     Total of         16 particles have been launched

 Pgm rebel. At pass #    4/  21.  In element #    1,  parameter # 35  changed to    1.63258789E+00   (was    1.42172526E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        21.1       21.100993       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14978372      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.194157E-05    1.48E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    9.802098E-04    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   9.6095E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        5 through the optical structure 

                     Total of         20 particles have been launched

 Pgm rebel. At pass #    5/  21.  In element #    1,  parameter # 35  changed to    1.84345052E+00   (was    1.63258789E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        23.8       23.826361       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14977138      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.592280E-05    1.89E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    1.157503E-03    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.3401E-06

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        6 through the optical structure 

                     Total of         24 particles have been launched

 Pgm rebel. At pass #    6/  21.  In element #    1,  parameter # 35  changed to    2.05431315E+00   (was    1.84345052E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        26.6       26.551729       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14976132      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.990390E-05    2.35E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    1.298392E-03    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.6862E-06

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        7 through the optical structure 

                     Total of         28 particles have been launched

 Pgm rebel. At pass #    7/  21.  In element #    1,  parameter # 35  changed to    2.26517578E+00   (was    2.05431315E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        29.3       29.277097       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14975354      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    2.388492E-05    2.86E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    1.413046E-03    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.9973E-06

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        8 through the optical structure 

                     Total of         32 particles have been launched

 Pgm rebel. At pass #    8/  21.  In element #    1,  parameter # 35  changed to    2.47603841E+00   (was    2.26517578E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        32.0       32.002546       65.0      2.693E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14975141      0.150      5.081E-08  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.252425E-05    3.41E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    6.778411E-04    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   4.5963E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        9 through the optical structure 

                     Total of         36 particles have been launched

 Pgm rebel. At pass #    9/  21.  In element #    1,  parameter # 35  changed to    2.68690104E+00   (was    2.47603841E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        34.7       34.727914       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14975857      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    8.543284E-06    4.02E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    4.260954E-04    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.8163E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       10 through the optical structure 

                     Total of         40 particles have been launched

 Pgm rebel. At pass #   10/  21.  In element #    1,  parameter # 35  changed to    2.89776367E+00   (was    2.68690104E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        37.5       37.453282       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14976132      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    4.562339E-06    4.67E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    2.109883E-04    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   4.4537E-08

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       11 through the optical structure 

                     Total of         44 particles have been launched

 Pgm rebel. At pass #   11/  21.  In element #    1,  parameter # 35  changed to    3.10862630E+00   (was    2.89776367E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        40.2       40.178650       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14978784      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    5.814088E-07    5.38E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    2.506381E-05    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   6.2853E-10

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       12 through the optical structure 

                     Total of         48 particles have been launched

 Pgm rebel. At pass #   12/  21.  In element #    1,  parameter # 35  changed to    3.31948893E+00   (was    3.10862630E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        42.9       42.904019       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14979012      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    3.399512E-06    6.13E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    1.372395E-04    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.8846E-08

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       13 through the optical structure 

                     Total of         52 particles have been launched

 Pgm rebel. At pass #   13/  21.  In element #    1,  parameter # 35  changed to    3.53035156E+00   (was    3.31948893E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        45.6       45.629387       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14981939      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    7.380424E-06    6.94E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    2.801543E-04    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   7.8541E-08

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       14 through the optical structure 

                     Total of         56 particles have been launched

 Pgm rebel. At pass #   14/  21.  In element #    1,  parameter # 35  changed to    3.74121419E+00   (was    3.53035156E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        48.4       48.354755       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14984911      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.136133E-05    7.79E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    4.069590E-04    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.6574E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       15 through the optical structure 

                     Total of         60 particles have been launched

 Pgm rebel. At pass #   15/  21.  In element #    1,  parameter # 35  changed to    3.95207682E+00   (was    3.74121419E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        51.1       51.080123       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14988157      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.534223E-05    8.69E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    5.202322E-04    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.7088E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       16 through the optical structure 

                     Total of         64 particles have been launched

 Pgm rebel. At pass #   16/  21.  In element #    1,  parameter # 35  changed to    4.16293945E+00   (was    3.95207682E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        53.8       53.805491       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14988797      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.932313E-05    9.64E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    6.220303E-04    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   3.8730E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       17 through the optical structure 

                     Total of         68 particles have been launched

 Pgm rebel. At pass #   17/  21.  In element #    1,  parameter # 35  changed to    4.37380208E+00   (was    4.16293945E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        56.5       56.530859       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14989255      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    2.330403E-05    1.06E-03 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    7.140128E-04    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   5.1036E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       18 through the optical structure 

                     Total of         72 particles have been launched

 Pgm rebel. At pass #   18/  21.  In element #    1,  parameter # 35  changed to    4.58466471E+00   (was    4.37380208E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        59.3       59.256227       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14995610      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    2.728493E-05    1.17E-03 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    7.975342E-04    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   6.3681E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       19 through the optical structure 

                     Total of         76 particles have been launched

 Pgm rebel. At pass #   19/  21.  In element #    1,  parameter # 35  changed to    4.79552734E+00   (was    4.58466471E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        62.0       61.981595       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14993919      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    3.126582E-05    1.28E-03 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    8.737106E-04    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   7.6435E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       20 through the optical structure 

                     Total of         80 particles have been launched


      Next  pass  is  #    21 and  last  pass  through  the  optical  structure


 Pgm rebel. At pass #   20/  21.  In element #    1,  parameter # 35  changed to    5.00638997E+00   (was    4.79552734E+00)

     WRITE statements to zgoubi.res are re-established from now on.

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       4 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  TOSCA                                                 


     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.52858543E+00   0.51861212E+00   0.10471976E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+02   0.70000000E+02   0.60000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      10.0     -0.00      /   0.00      10.0     -0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  106/ 121/   1
     Node  distance  in   X/Y/Z :   9.973310E-03/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :  0.1499     cm   (i.e.,   3.7485E-03 rad  at mean radius RM =    40.00    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  5.0064    61.982     0.000     0.000     0.000            0.519    63.301     0.036     0.000     0.000            1
  A    1  2.2365    28.907     0.000     0.000     0.000            0.519    28.907    -0.000     0.000     0.000            2
  A    1  3.8685    50.000     0.000     0.000     0.000            0.519    50.000    -0.000     0.000     0.000            3
  A    1  5.0064    64.707     0.000     0.000     0.000            0.519    64.707    -0.000     0.000     0.000            4

 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  4 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   5.0064    61.982     0.000     0.000     0.000       0.0000    4.0064   63.301   36.485    0.000    0.000   6.540022E+01     1
               Time of flight (mus) :  2.12155011E-02 mass (MeV/c2) :   938.272    
m  1   2.2365    28.907     0.000     0.000     0.000       0.0000    1.2365   28.907   -0.000    0.000    0.000   3.027143E+01     2
               Time of flight (mus) :  2.18881269E-02 mass (MeV/c2) :   938.272    
o  1   3.8685    50.000     0.000     0.000     0.000       0.0000    2.8685   50.000   -0.000    0.000    0.000   5.235988E+01     3
               Time of flight (mus) :  2.19344684E-02 mass (MeV/c2) :   938.272    
M  1   5.0064    64.707     0.000     0.000     0.000       0.0000    4.0064   64.707   -0.000    0.000    0.000   6.776105E+01     4
               Time of flight (mus) :  2.19813400E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   6.3153E-03  -5.2509E-01   1.0275E+01   5.172884E-01   9.121198E-03        4        4    1.000      (Y,T)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        4        4    1.000      (Z,P)        21
   1.4733E-03   4.7674E-01   2.0908E-04   2.175486E-02   3.497156E+00        4        4    1.000      (t,K)        21

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   1.437163E-01
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   1.579838E-02

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   3.131375E-04
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   1.659177E+00


  Beam  sigma  matrix : 

   2.065439E-02   1.055542E-03   0.000000E+00   0.000000E+00
   1.055542E-03   2.495888E-04   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     2.010207E-03    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      5  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1


                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       4 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                                              


************************************************************************************************************************************
      3  Keyword, label(s) :  TOSCA                                                 


     NDIM =   2 ;  Number of data file sets used is   1 ;  Stored in field array # IMAP =    1 ;  
     Value of MOD.MOD2 is 22.1

           3-D map. MOD=22 or 23.  Single field map, with field coefficient value : 
                 1.000000E+00
          No  new  map  file  to  be  opened. Already  stored.
          Skip  reading  field  map  file :           geneSectorMap.out
 Pgm toscap,  restored mesh coordinates for field map #   1,  name : geneSectorMap.out


     Field map limits, angle :  min, max, max-min (rad) :  -0.52858543E+00   0.51861212E+00   0.10471976E+01
     Field map limits, radius :  min, max, max-min (cm) :   0.10000000E+02   0.70000000E+02   0.60000000E+02
      Min / max  fields  drawn  from  map  data :    0.0000                     /    5.0000    
       @  X-node, Y-node, z-node :                0.00      10.0     -0.00      /   0.00      10.0     -0.00    
     Normalisation  coeff.  BNORM   :   1.00000    
     Field  min/max  normalised  :                0.0000        /    5.0000    
     Nbre  of  nodes  in  X/Y/Z :  106/ 121/   1
     Node  distance  in   X/Y/Z :   9.973310E-03/  0.500000    /   0.00000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :  0.1499     cm   (i.e.,   3.7485E-03 rad  at mean radius RM =    40.00    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  5.0064    61.982     0.000     0.000     0.000            0.519    63.301     0.036     0.000     0.000            1
  A    1  2.2365    28.907     0.000     0.000     0.000            0.519    28.907    -0.000     0.000     0.000            2
  A    1  3.8685    50.000     0.000     0.000     0.000            0.519    50.000    -0.000     0.000     0.000            3
  A    1  5.0064    64.707     0.000     0.000     0.000            0.519    64.707    -0.000     0.000     0.000            4

 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  4 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   5.0064    61.982     0.000     0.000     0.000       0.0000    4.0064   63.301   36.485    0.000    0.000   6.540022E+01     1
               Time of flight (mus) :  2.12155011E-02 mass (MeV/c2) :   938.272    
m  1   2.2365    28.907     0.000     0.000     0.000       0.0000    1.2365   28.907   -0.000    0.000    0.000   3.027143E+01     2
               Time of flight (mus) :  2.18881269E-02 mass (MeV/c2) :   938.272    
o  1   3.8685    50.000     0.000     0.000     0.000       0.0000    2.8685   50.000   -0.000    0.000    0.000   5.235988E+01     3
               Time of flight (mus) :  2.19344684E-02 mass (MeV/c2) :   938.272    
M  1   5.0064    64.707     0.000     0.000     0.000       0.0000    4.0064   64.707   -0.000    0.000    0.000   6.776105E+01     4
               Time of flight (mus) :  2.19813400E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   6.3153E-03  -5.2509E-01   1.0275E+01   5.172884E-01   9.121198E-03        4        4    1.000      (Y,T)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        4        4    1.000      (Z,P)        21
   1.4733E-03   4.7674E-01   2.0908E-04   2.175486E-02   3.497156E+00        4        4    1.000      (t,K)        21

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   1.437163E-01
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   1.579838E-02

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   3.131375E-04
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   1.659177E+00


  Beam  sigma  matrix : 

   2.065439E-02   1.055542E-03   0.000000E+00   0.000000E+00
   1.055542E-03   2.495888E-04   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     2.010207E-03    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      5  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           variable #            2       IR =            3 ,   ok.
           variable #            2       IP =           50 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.

                    FIT  variables  and  constraints  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        64.7       64.706963       65.0      8.078E-05  OBJET      -                    -                   
   3   2    50   5.000E-02   0.150      0.14992593      0.150      1.524E-07  TOSCA      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-06)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    3.524671E-05    1.39E-03 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    9.434700E-04    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   8.9138E-07

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai                                                                      
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      7  Keyword, label(s) :  REBELOTE                                              


                         ****  End  of  'REBELOTE'  procedure  ****

      There  has  been         21  passes  through  the  optical  structure 

                     Total of         84 particles have been launched

************************************************************************************************************************************
      8  Keyword, label(s) :  END                                                   


************************************************************************************************************************************
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  12-06-2017,  at  07:10:18 
  JOB  ENDED  ON    12-06-2017,  AT  07:10:20 

   CPU time, total :     2.0561279999999997     
