Uniform field sector
 'OBJET'                                                                                                      1
64.62444403717985                       ! 200keV proton
2
1 1
12.9248888074 0.  0.   0.   0.  1. 'm'
1
 
 'PARTICUL'              ! This is required only because we want to get the time-of-flight                    2
938.27203D0 1.602176487D-19  1.79284735D0 0. 0.  ! otherwise zgoubi only requires rigidity.
! 'PARTICUL'
! PROTON
 
 'FAISTORE'                                                                                                   3
zgoubi.fai   #End
1
 'TOSCA'                                                                                                      4
0  2
1.  1. 1. 1.
HEADER_8
629 151 1 22.1 1.       IZ=1 -> 2D ; MOD=22 -> polar map ; .MOD2=.1 -> one map file
geneSectorMap.out
0 0 0 0
2
1.
2
0. 0. 0. 0.
 
 'FAISCEAU'  #End                                                                                             5
 'END'                                                                                                        6

************************************************************************************************************************************
      1  Keyword, label(s) :  OBJET                             

                          MAGNETIC  RIGIDITY =         64.624 kG*cm

                                         TRAJECTOIRY SETTING UP

                              OBJET  (2)  BUILT  UP  FROM       1 POINTS 



************************************************************************************************************************************
      2  Keyword, label(s) :  PARTICUL                          


     Particle  properties :

                     Mass          =    938.272        MeV/c2
                     Charge        =   1.602176E-19    C     
                     G  factor     =    1.79285              

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
  
 I, AMQ(1,I), AMQ(2,I)/QE, P/Pref, v/c, time :
  
     1   9.38272030E+02  1.00000000E+00  1.00000000E+00  2.06441112E-02  0.00000000E+00

************************************************************************************************************************************
      3  Keyword, label(s) :  FAISTORE                          


                OPEN FILE zgoubi.fai                                                                      
                FOR PRINTING COORDINATES 

               Print will occur at element[s] labeled : 
                    #End      

************************************************************************************************************************************
      4  Keyword, label(s) :  TOSCA                             


                OPEN FILE zgoubi.plt                                                                      
                FOR PRINTING TRAJECTORIES


     NDIM = 2 ;   Value of MOD is  22 ;  Number of data file sets used is   1


          New field map(s) now used, polar mesh (MOD .ge. 20) ; 
           name(s) of map data file(s) are : 
          geneSectorMap.out                                                               

geneSectorMap.out map,  FORMAT type :  regular            
 IMAP :            1
 HEADER  (8 lines) : 
        1.0000000000000000       0.50000000000000000       0.57324840764331209       2.12199579145933796E-314      !  Rmi/cm,
      # Field map generated using geneSectorMap.f                                                                            
     # AT/rd,  AT/deg, Rmi/cm, Rma/cm, RM/cm, NR, dR/cm, NX, dX/cm, dA/rd :                                                  
     #   6.28318531E+00   3.60000000E+02   1.00000000E+00   7.60000000E+01   5.00000000E+01 151   5.00000000E-01 629   5.0025
      # For TOSCA:          629         151  1 22.1 1.  !IZ=1 -> 2D ; MOD=22 -> polar map ; .MOD2=.1 -> one map file         
      # R*cosA (A:0->360), Z==0, R*sinA, BY, BZ, BX                                                                          
      # cm                 cm    cm      kG  kG  kG                                                                          
      #                                                                                                                      
 
R0, DR, DTTA, DZ :   1.000000E+00  5.000000E-01  5.732484E-01  2.121996-314



     Field map limits, angle, min, max, max-min (rad) :  -0.314159E+01   0.314159E+01   0.628319E+01
     Field map limits, radius, min, max, max-min (cm) :   0.100000E+01   0.760000E+02   0.750000E+02
      Min / max  fields  drawn  from  map  data :    0.00                       /    5.00    
       @  x-node, y-node, z-node :                0.00      1.00      0.00      /   0.00      1.00      0.00    
     Normalisation  coeff.  BNORM   :   1.000    
     Field  min/max  normalised  :                0.00          /    5.00    
     Nbre  of  nodes  in  x/y/z :  629/ 151/   1
     Node  distance  in   x/y/z :   1.0005E-02/  0.5000    /   0.000    

                     Option  for  interpolation : 2
                     Smoothing  using  9  points 

                    Integration step :   1.000     cm   (i.e.,   2.5974E-02 rad  at mean radius RM =    38.50    )

                               KPOS = 2.  Position of reference orbit on mechanical  faces
                                         at entrance    RE =   0.00000     cm  TE =   0.00000     rad
                                         at exit        RS =   0.00000     cm  TS =   0.00000     rad

  A    1  1.0000    12.925     0.000     0.000     0.000            3.142    12.925     0.000     0.000     0.000            1


                CONDITIONS  DE  MAXWELL  (       82.  PAS )  :
                       DIV(B)        LAPLACIEN(B)     ROTATIONNEL(B)
                      0.000            0.000             0.000    
                                       0.000             0.000    
                                       0.000             0.000    
                       LAPLACIEN SCALAIRE =   0.000    



 Cumulative length of optical axis =    0.00000000     m ;  Time  (for ref. rigidity & particle) =    0.00000     s 

************************************************************************************************************************************
      5  Keyword, label(s) :  FAISCEAU    #End                  

0                                             TRACE DU FAISCEAU

                                                  1 TRAJECTOIRES

                                   OBJET                                                  FAISCEAU

          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)          D       Y(CM)    T(MR)    Z(CM)    P(MR)      S(CM)

m  1   1.0000    12.925     0.000     0.000     0.000       0.0000    1.0000   12.925    0.000    0.000    0.000   8.120937E+01     1
               Time of flight (mus) :  0.13121675     mass (MeV/c2) :   938.272    


  Beam  characteristics   (EMIT,ALP,BET,XM,XPM,NLIV,NINL,RATIN) : 

    0.000        0.000        0.000       0.1292       2.4120E-08       1       0    0.000    B-Dim     1      1
    0.000        0.000        0.000        0.000        0.000           1       0    0.000    B-Dim     2      1
    0.000        0.000        0.000       0.1312       0.2000           1       0    0.000    B-Dim     3      1


  Beam  characteristics   SIGMA(4,4) : 

           Ex, Ez =  0.000000E+00  0.000000E+00
           AlpX, BetX =           NaN      Infinity
           AlpZ, BetZ =           NaN           NaN

  4.057292E-37  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00

************************************************************************************************************************************
      6  Keyword, label(s) :  END                               


                             1 particles have been launched
                     Made  it  to  the  end :      1

************************************************************************************************************************************

           MAIN PROGRAM : Execution ended upon key  END       

************************************************************************************************************************************
   
            ZGOUBI RUN COMPLETED. 

  Zgoubi version 6.0 16-Jun-2015.
  Job  started  on  06-02-2018,  at  07:09:01 
  JOB  ENDED  ON    06-02-2018,  AT  07:09:01 

   CPU time, total :    0.43602599999999997     
