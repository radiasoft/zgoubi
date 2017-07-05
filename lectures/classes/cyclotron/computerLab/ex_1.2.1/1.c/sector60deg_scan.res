Uniform field sector
 'OBJET'                                                                                                      1
64.62444403717985                                   ! 200keV
2
1 1
12.9248888074 0.  0.   0.   0.  1. 'm'              ! 200keV. R=Brho/B=*/.5
1
 
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
.01         !  5.  0.1 .03  .01    Need step < rho*AT ~ 13+cm
2
0. 0. 0. 0.
 
 'FAISCEAU'                                                                                                   4
 
 'FIT'                                                                                                        5
1   nofinal
1 30 0 [12.,65.]        ! Variable : Y_0 (cm)
2  2e-10  199           ! Penalty; max numb of calls to function
3.1 1 2 #End 0. 1. 0    ! Constraint :  Y_final=Y_0
3   1 3 #End 0. 1. 0    ! Constraint :  T_final=0 (by symmetry)
 
 'FAISCEAU'                                                                                                   6
 
! 'FAISTORE'                                                                                                   6
! zgoubi.fai_stepSize13cm
! 1
 
! 'FAISTORE'                                                                                                   6
! zgoubi.fai_stepSize5cm
! 1
 
! 'FAISTORE'                                                                                                   6
! zgoubi.fai_stepSize0.1cm
! 1
 
! 'FAISTORE'                                                                                                   7
! zgoubi.fai_stepSize0.03cm
! 1
 
 'FAISTORE'                                                                                                   7
zgoubi.fai_stepSize0.01cm
1
 
 'REBELOTE'                                                                                                   8
20 0.2  0 1        ! 20 different rigidities ;  ; coordinates as found in OBJET ; change parameter(s)
1
OBJET 35     1:5.0063899693    ! 0.2 MeV to 5 MeV
 
 'END'                                                                                                        9

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



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

                    Integration step :  1.0000E-02 cm   (i.e.,   2.5000E-04 rad  at mean radius RM =    40.00    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            0.519    12.925     0.000     0.000     0.000            1

 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    0.0000   12.925    0.000    0.000    0.000   1.353491E+01     1
               Time of flight (mus) :  2.18694843E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01   2.388097E-12        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   2.186948E-02   2.000000E-01        1        1    1.000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      5  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.

                    FIT  variables  and  constraints  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        12.9       12.924889       65.0      0.177      OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.818456E-11    5.80E-05 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    2.388097E-09    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   5.7033E-18

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    0.0000   12.925    0.000    0.000    0.000   1.353491E+01     1
               Time of flight (mus) :  2.18694843E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01   2.388097E-12        1        1    1.000      (Y,T)         1
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         1
   0.0000E+00   0.0000E+00   1.0000E+00   2.186948E-02   2.000000E-01        1        1    1.000      (t,K)         1

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


                OPEN FILE zgoubi.fai_stepSize0.01cm                                                       
                FOR PRINTING COORDINATES 

               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        1 through the optical structure 

                     Total of          1 particles have been launched

     Multiple pass, 
          from element #     1 : OBJET     /label1=                    /label2=                    
                             to  REBELOTE  /label1=                    /label2=                    
     ending at pass #      21 at element #     8 : REBELOTE  /label1=                    /label2=                    


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
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.818456E-11    5.80E-05 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00    2.388097E-09    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   5.7033E-18

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    0.0000   12.925    0.000    0.000    0.000   1.353491E+01     1
               Time of flight (mus) :  2.18694843E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.292489E-01   2.388097E-12        1        1    1.000      (Y,T)         2
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         2
   0.0000E+00   0.0000E+00   1.0000E+00   2.186948E-02   2.000000E-01        1        1    1.000      (t,K)         2

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        2 through the optical structure 

                     Total of          2 particles have been launched

 Pgm rebel. At pass #    2/  21.  In element #    1,  parameter # 35  changed to    1.21086263E+00   (was    1.00000000E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        15.7       15.650265       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    8.244212E-09    8.16E-05 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -9.124008E-07    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   8.3254E-13

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.2109    15.650     0.000     0.000     0.000       0.0000    0.2109   15.650   -0.000    0.000    0.000   1.638892E+01     1
               Time of flight (mus) :  2.18716567E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.565026E-01  -9.124008E-10        1        1    1.000      (Y,T)         3
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         3
   0.0000E+00   0.0000E+00   1.0000E+00   2.187166E-02   2.932231E-01        1        1    1.000      (t,K)         3

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        3 through the optical structure 

                     Total of          3 particles have been launched

 Pgm rebel. At pass #    3/  21.  In element #    1,  parameter # 35  changed to    1.42172526E+00   (was    1.21086263E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        18.4       18.375641       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.650539E-08    1.13E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -1.555863E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   2.4210E-12

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.4217    18.376     0.000     0.000     0.000       0.0000    0.4217   18.376   -0.000    0.000    0.000   1.924293E+01     1
               Time of flight (mus) :  2.18742432E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   1.837564E-01  -1.555863E-09        1        1    1.000      (Y,T)         4
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         4
   0.0000E+00   0.0000E+00   1.0000E+00   2.187424E-02   4.042166E-01        1        1    1.000      (t,K)         4

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        4 through the optical structure 

                     Total of          4 particles have been launched

 Pgm rebel. At pass #    4/  21.  In element #    1,  parameter # 35  changed to    1.63258789E+00   (was    1.42172526E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        21.1       21.101017       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    2.476805E-08    1.48E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -2.033091E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   4.1341E-12

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.6326    21.101     0.000     0.000     0.000       0.0000    0.6326   21.101   -0.000    0.000    0.000   2.209693E+01     1
               Time of flight (mus) :  2.18772437E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   2.110102E-01  -2.033091E-09        1        1    1.000      (Y,T)         5
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         5
   0.0000E+00   0.0000E+00   1.0000E+00   2.187724E-02   5.329741E-01        1        1    1.000      (t,K)         5

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        5 through the optical structure 

                     Total of          5 particles have been launched

 Pgm rebel. At pass #    5/  21.  In element #    1,  parameter # 35  changed to    1.84345052E+00   (was    1.63258789E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        23.8       23.826393       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    3.302813E-08    1.89E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -2.401203E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   5.7669E-12

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   1.8435    23.826     0.000     0.000     0.000       0.0000    0.8435   23.826   -0.000    0.000    0.000   2.495094E+01     1
               Time of flight (mus) :  2.18806580E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   2.382639E-01  -2.401203E-09        1        1    1.000      (Y,T)         6
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         6
   0.0000E+00   0.0000E+00   1.0000E+00   2.188066E-02   6.794884E-01        1        1    1.000      (t,K)         6

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        6 through the optical structure 

                     Total of          6 particles have been launched

 Pgm rebel. At pass #    6/  21.  In element #    1,  parameter # 35  changed to    2.05431315E+00   (was    1.84345052E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        26.6       26.551769       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    4.129071E-08    2.35E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -2.693688E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   7.2577E-12

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   2.0543    26.552     0.000     0.000     0.000       0.0000    1.0543   26.552   -0.000    0.000    0.000   2.780495E+01     1
               Time of flight (mus) :  2.18844858E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   2.655177E-01  -2.693688E-09        1        1    1.000      (Y,T)         7
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         7
   0.0000E+00   0.0000E+00   1.0000E+00   2.188449E-02   8.437511E-01        1        1    1.000      (t,K)         7

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        7 through the optical structure 

                     Total of          7 particles have been launched

 Pgm rebel. At pass #    7/  21.  In element #    1,  parameter # 35  changed to    2.26517578E+00   (was    2.05431315E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        29.3       29.277145       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    4.955440E-08    2.86E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -2.931691E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   8.5973E-12

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   2.2652    29.277     0.000     0.000     0.000       0.0000    1.2652   29.277   -0.000    0.000    0.000   3.065895E+01     1
               Time of flight (mus) :  2.18887271E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   2.927715E-01  -2.931691E-09        1        1    1.000      (Y,T)         8
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         8
   0.0000E+00   0.0000E+00   1.0000E+00   2.188873E-02   1.025753E+00        1        1    1.000      (t,K)         8

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        8 through the optical structure 

                     Total of          8 particles have been launched

 Pgm rebel. At pass #    8/  21.  In element #    1,  parameter # 35  changed to    2.47603841E+00   (was    2.26517578E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        32.0       32.002521       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    5.781399E-08    3.41E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -3.129201E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   9.7952E-12

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   2.4760    32.003     0.000     0.000     0.000       0.0000    1.4760   32.003   -0.000    0.000    0.000   3.351296E+01     1
               Time of flight (mus) :  2.18933815E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   3.200252E-01  -3.129201E-09        1        1    1.000      (Y,T)         9
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)         9
   0.0000E+00   0.0000E+00   1.0000E+00   2.189338E-02   1.225484E+00        1        1    1.000      (t,K)         9

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #        9 through the optical structure 

                     Total of          9 particles have been launched

 Pgm rebel. At pass #    9/  21.  In element #    1,  parameter # 35  changed to    2.68690104E+00   (was    2.47603841E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        34.7       34.727897       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    6.607563E-08    4.02E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -3.295759E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.0866E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   2.6869    34.728     0.000     0.000     0.000       0.0000    1.6869   34.728   -0.000    0.000    0.000   3.636697E+01     1
               Time of flight (mus) :  2.18984487E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   3.472790E-01  -3.295759E-09        1        1    1.000      (Y,T)        10
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        10
   0.0000E+00   0.0000E+00   1.0000E+00   2.189845E-02   1.442932E+00        1        1    1.000      (t,K)        10

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       10 through the optical structure 

                     Total of         10 particles have been launched

 Pgm rebel. At pass #   10/  21.  In element #    1,  parameter # 35  changed to    2.89776367E+00   (was    2.68690104E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        37.5       37.453273       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    7.434068E-08    4.67E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -3.437921E-06    1.00E+00 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.1825E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   2.8978    37.453     0.000     0.000     0.000       0.0000    1.8978   37.453   -0.000    0.000    0.000   3.922098E+01     1
               Time of flight (mus) :  2.19039286E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   3.745327E-01  -3.437921E-09        1        1    1.000      (Y,T)        11
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        11
   0.0000E+00   0.0000E+00   1.0000E+00   2.190393E-02   1.678085E+00        1        1    1.000      (t,K)        11

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       11 through the optical structure 

                     Total of         11 particles have been launched

 Pgm rebel. At pass #   11/  21.  In element #    1,  parameter # 35  changed to    3.10862630E+00   (was    2.89776367E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        40.2       40.178649       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    8.259841E-08    5.38E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -3.560985E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.2687E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   3.1086    40.179     0.000     0.000     0.000       0.0000    2.1086   40.179   -0.000    0.000    0.000   4.207498E+01     1
               Time of flight (mus) :  2.19098207E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   4.017865E-01  -3.560985E-09        1        1    1.000      (Y,T)        12
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        12
   0.0000E+00   0.0000E+00   1.0000E+00   2.190982E-02   1.930931E+00        1        1    1.000      (t,K)        12

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       12 through the optical structure 

                     Total of         12 particles have been launched

 Pgm rebel. At pass #   12/  21.  In element #    1,  parameter # 35  changed to    3.31948893E+00   (was    3.10862630E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        42.9       42.904025       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    9.086330E-08    6.13E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -3.668315E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.3465E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   3.3195    42.904     0.000     0.000     0.000       0.0000    2.3195   42.904   -0.000    0.000    0.000   4.492899E+01     1
               Time of flight (mus) :  2.19161248E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   4.290403E-01  -3.668315E-09        1        1    1.000      (Y,T)        13
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        13
   0.0000E+00   0.0000E+00   1.0000E+00   2.191612E-02   2.201454E+00        1        1    1.000      (t,K)        13

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       13 through the optical structure 

                     Total of         13 particles have been launched

 Pgm rebel. At pass #   13/  21.  In element #    1,  parameter # 35  changed to    3.53035156E+00   (was    3.31948893E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        45.6       45.629402       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    9.912580E-08    6.94E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -3.762796E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.4168E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   3.5304    45.629     0.000     0.000     0.000       0.0000    2.5304   45.629   -0.000    0.000    0.000   4.778300E+01     1
               Time of flight (mus) :  2.19228405E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   4.562940E-01  -3.762796E-09        1        1    1.000      (Y,T)        14
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        14
   0.0000E+00   0.0000E+00   1.0000E+00   2.192284E-02   2.489639E+00        1        1    1.000      (t,K)        14

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       14 through the optical structure 

                     Total of         14 particles have been launched

 Pgm rebel. At pass #   14/  21.  In element #    1,  parameter # 35  changed to    3.74121419E+00   (was    3.53035156E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        48.4       48.354778       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.074013E-07    7.79E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -3.846504E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.4807E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   3.7412    48.355     0.000     0.000     0.000       0.0000    2.7412   48.355   -0.000    0.000    0.000   5.063700E+01     1
               Time of flight (mus) :  2.19299673E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   4.835478E-01  -3.846504E-09        1        1    1.000      (Y,T)        15
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        15
   0.0000E+00   0.0000E+00   1.0000E+00   2.192997E-02   2.795471E+00        1        1    1.000      (t,K)        15

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       15 through the optical structure 

                     Total of         15 particles have been launched

 Pgm rebel. At pass #   15/  21.  In element #    1,  parameter # 35  changed to    3.95207682E+00   (was    3.74121419E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        51.1       51.080154       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.156497E-07    8.69E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -3.921606E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.5392E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   3.9521    51.080     0.000     0.000     0.000       0.0000    2.9521   51.080   -0.000    0.000    0.000   5.349101E+01     1
               Time of flight (mus) :  2.19375050E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   5.108015E-01  -3.921606E-09        1        1    1.000      (Y,T)        16
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        16
   0.0000E+00   0.0000E+00   1.0000E+00   2.193751E-02   3.118931E+00        1        1    1.000      (t,K)        16

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       16 through the optical structure 

                     Total of         16 particles have been launched

 Pgm rebel. At pass #   16/  21.  In element #    1,  parameter # 35  changed to    4.16293945E+00   (was    3.95207682E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        53.8       53.805530       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.239115E-07    9.64E-04 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -3.988912E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.5927E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   4.1629    53.806     0.000     0.000     0.000       0.0000    3.1629   53.806   -0.000    0.000    0.000   5.634502E+01     1
               Time of flight (mus) :  2.19454531E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   5.380553E-01  -3.988912E-09        1        1    1.000      (Y,T)        17
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        17
   0.0000E+00   0.0000E+00   1.0000E+00   2.194545E-02   3.460003E+00        1        1    1.000      (t,K)        17

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       17 through the optical structure 

                     Total of         17 particles have been launched

 Pgm rebel. At pass #   17/  21.  In element #    1,  parameter # 35  changed to    4.37380208E+00   (was    4.16293945E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        56.5       56.530906       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.321859E-07    1.06E-03 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -4.049643E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.6417E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   4.3738    56.531     0.000     0.000     0.000       0.0000    3.3738   56.531   -0.000    0.000    0.000   5.919903E+01     1
               Time of flight (mus) :  2.19538112E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   5.653091E-01  -4.049643E-09        1        1    1.000      (Y,T)        18
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        18
   0.0000E+00   0.0000E+00   1.0000E+00   2.195381E-02   3.818666E+00        1        1    1.000      (t,K)        18

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       18 through the optical structure 

                     Total of         18 particles have been launched

 Pgm rebel. At pass #   18/  21.  In element #    1,  parameter # 35  changed to    4.58466471E+00   (was    4.37380208E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        59.3       59.256282       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.404454E-07    1.17E-03 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -4.104924E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.6870E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   4.5847    59.256     0.000     0.000     0.000       0.0000    3.5847   59.256   -0.000    0.000    0.000   6.205303E+01     1
               Time of flight (mus) :  2.19625787E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   5.925628E-01  -4.104924E-09        1        1    1.000      (Y,T)        19
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        19
   0.0000E+00   0.0000E+00   1.0000E+00   2.196258E-02   4.194901E+00        1        1    1.000      (t,K)        19

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       19 through the optical structure 

                     Total of         19 particles have been launched

 Pgm rebel. At pass #   19/  21.  In element #    1,  parameter # 35  changed to    4.79552734E+00   (was    4.58466471E+00)

     WRITE statements to zgoubi.res are inhibeted from now on.

************************************************************************************************************************************

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        62.0       61.981658       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.486871E-07    1.28E-03 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -4.155469E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.7290E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   4.7955    61.982     0.000     0.000     0.000       0.0000    3.7955   61.982   -0.000    0.000    0.000   6.490704E+01     1
               Time of flight (mus) :  2.19717552E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   6.198166E-01  -4.155469E-09        1        1    1.000      (Y,T)        20
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        20
   0.0000E+00   0.0000E+00   1.0000E+00   2.197176E-02   4.588686E+00        1        1    1.000      (t,K)        20

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                                -----  REBELOTE  -----

     End of pass #       20 through the optical structure 

                     Total of         20 particles have been launched


      Next  pass  is  #    21 and  last  pass  through  the  optical  structure


 Pgm rebel. At pass #   20/  21.  In element #    1,  parameter # 35  changed to    5.00638997E+00   (was    4.79552734E+00)

     WRITE statements to zgoubi.res are re-established from now on.

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



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

                    Integration step :  1.0000E-02 cm   (i.e.,   2.5000E-04 rad  at mean radius RM =    40.00    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  5.0064    61.982     0.000     0.000     0.000            0.519    63.301     0.036     0.000     0.000            1

 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   5.0064    61.982     0.000     0.000     0.000       0.0000    4.0064   63.301   36.484    0.000    0.000   6.540028E+01     1
               Time of flight (mus) :  2.12155188E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   6.330129E-01   3.648395E-02        1        1    1.000      (Y,T)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        21
   0.0000E+00   0.0000E+00   1.0000E+00   2.121552E-02   5.000000E+00        1        1    1.000      (t,K)        21

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      5  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1


                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.



************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                                                 

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



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

                    Integration step :  1.0000E-02 cm   (i.e.,   2.5000E-04 rad  at mean radius RM =    40.00    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  5.0064    61.982     0.000     0.000     0.000            0.519    63.301     0.036     0.000     0.000            1

 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      4  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      3)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   5.0064    61.982     0.000     0.000     0.000       0.0000    4.0064   63.301   36.484    0.000    0.000   6.540028E+01     1
               Time of flight (mus) :  2.12155188E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   6.330129E-01   3.648395E-02        1        1    1.000      (Y,T)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        21
   0.0000E+00   0.0000E+00   1.0000E+00   2.121552E-02   5.000000E+00        1        1    1.000      (t,K)        21

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      5  Keyword, label(s) :  FIT                                                   

     FIT procedure launched. Method is 1

           variable #            1       IR =            1 ,   ok.
           variable #            1       IP =           30 ,   ok.
           constraint #            1       IR =            4 ,   ok.
           constraint #            1       I  =            1 ,   ok.
           constraint #            2       IR =            4 ,   ok.
           constraint #            2       I  =            1 ,   ok.

                    FIT  variables  and  constraints  in  good  order,  FIT  will proceed. 

                    Final FIT status will NOT be saved. For so, use the 'save [FileName]' command

 STATUS OF VARIABLES  (Iteration #     0 /    199 max.)
LMNT VAR PARAM  MINIMUM    INITIAL         FINAL         MAXIMUM     STEP        NAME   LBL1                 LBL2
   1   1    30    12.0        64.7       64.707034       65.0      3.324E-07  OBJET      -                    -                   
 STATUS OF CONSTRAINTS (Target penalty =   2.0000E-10)
TYPE  I   J LMNT#     DESIRED          WEIGHT         REACHED         KI2     NAME   LBL1                 LBL2      Nb param. [value]
  3   1   2     4    0.000000E+00    1.000E+00    1.569583E-07    1.39E-03 FAISCEAU   -                    -                    0
  3   1   3     4    0.000000E+00    1.000E+00   -4.201552E-06    9.99E-01 FAISCEAU   -                    -                    0
 Fit reached penalty value   1.7678E-11

   Last run following FIT[2] is skipped,as requested.  Now carrying on beyond FIT keyword.

************************************************************************************************************************************
      6  Keyword, label(s) :  FAISCEAU                                              

0                                             TRACE DU FAISCEAU
                                           (follows element #      5)
                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(cm)     T(mr)     Z(cm)     P(mr)       S(cm)       D-1     Y(cm)    T(mr)    Z(cm)    P(mr)      S(cm)

m  1   5.0064    64.707     0.000     0.000     0.000       0.0000    4.0064   64.707   -0.000    0.000    0.000   6.776105E+01     1
               Time of flight (mus) :  2.19813401E-02 mass (MeV/c2) :   938.272    


------
  Characteristics of concentration ellipse (Surface, ALP, BET, <X>, <XP>, #prtcls, #prtcls inside ellips, ratio, space, pass#) : 

   0.0000E+00   0.0000E+00   1.0000E+00   6.470703E-01  -4.201552E-09        1        1    1.000      (Y,T)        21
   0.0000E+00   0.0000E+00   1.0000E+00   0.000000E+00   0.000000E+00        1        1    1.000      (Z,P)        21
   0.0000E+00   0.0000E+00   1.0000E+00   2.198134E-02   5.000000E+00        1        1    1.000      (t,K)        21

(Y,T)  space (units : (cm,rd)   ) :  
      sigma_Y = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_T = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(Z,P)  space (units : (cm,rd)   ) :  
      sigma_Z = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_P = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00

(t,K)  space (units : (mu_s,MeV)) :  
      sigma_t = sqrt(Surface/pi * BET) =   0.000000E+00
      sigma_K = sqrt(Surface/pi * (1+ALP^2)/BET) =   0.000000E+00


  Beam  sigma  matrix : 

   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

      sqrt(det_Y), sqrt(det_Z) :     0.000000E+00    0.000000E+00    (Note :  sqrt(determinant) = ellipse surface / pi)

************************************************************************************************************************************
      7  Keyword, label(s) :  FAISTORE                                              


     zgoubi.fai_stepSize0.01cm                                                       
      already open...
               Print will occur at element[s] labeled : 


************************************************************************************************************************************
      8  Keyword, label(s) :  REBELOTE                                              


                         ****  End  of  'REBELOTE'  procedure  ****

      There  has  been         21  passes  through  the  optical  structure 

                     Total of         21 particles have been launched

************************************************************************************************************************************
      9  Keyword, label(s) :  END                                                   


************************************************************************************************************************************
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi, author's dvlpmnt version.
  Job  started  on  05-07-2017,  at  14:22:01 
  JOB  ENDED  ON    05-07-2017,  AT  14:22:02 

   CPU time, total :     1.4680000000000000     

                    An updated version of the input data file, with variables in FIT'ed state, has been saved in zgoubi.FIT.out.dat.


