      parameter (G=1.7928474, am = 938.27203e6, c=2.99792458e8)

      qz =  8.76346

      T1 = (5./G - 1. ) * am
      e1 = T1 + am
      g1 = e1/am
      p1 =  sqrt(T1 * (T1 + 2.*am))
      bta1 = p1 / e1
      bg1 = bta1 * g1
      epsi0 = 20./bg1

        write(*,*) ' Imperfection resonances '
      do i = 5, 51, 1
        gG = i
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        if(T.ge.T1 .and. T.le.26.e9) 
     >    write(*,*) 'gG, T(GeV), p(GeV/c) : ',
     >    gG, T/1e9, p/1e9,turn,sqrt(bg1/bg),bg1/bg*epsi0,brho
     >    ,bg1*epsi0
      enddo

      write(*,*) ' '
      write(*,*) ' Intrinsic resonances, random, k - qz  (qz=',qz,')'
        M=12

        k=0
        gG = k*M + qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        if(T.ge.T1 .and. T.le.26.e9) 
     >            write(*,*) ' kM+qz gG, T(GeV), p(GeV/c) : ',
     >  k*M,'+qz  ',gG, T/1e9, p/1e9,turn,sqrt(bg1/bg),bg1/bg*epsi0,brho
        k=1
        gG = k*M + qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        if(T.ge.T1 .and. T.le.26.e9) 
     >            write(*,*) '  kM+qz  gG, T(GeV), p(GeV/c) : ',
     >  k*M,'+qz  ',gG, T/1e9, p/1e9,turn,sqrt(bg1/bg),bg1/bg*epsi0,brho
        k=1
        gG = k*M - qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        if(T.ge.T1 .and. T.le.26.e9) 
     >            write(*,*) '  kM-qz gG, T(GeV), p(GeV/c) : ',
     >  k*M,'-qz  ',gG, T/1e9, p/1e9,turn,sqrt(bg1/bg),bg1/bg*epsi0,brho
        k=2
        gG = k*M + qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        if(T.ge.T1 .and. T.le.26.e9) 
     >            write(*,*) '  kM+qz gG, T(GeV), p(GeV/c) : ',
     >     k*M,'+qz  ',gG, T/1e9, p/1e9,turn,sqrt(bg1/bg),bg1/bg*epsi0
        k=2
        gG = k*M - qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        if(T.ge.T1 .and. T.le.26.e9) 
     >            write(*,*) '  kM-qz gG, T(GeV), p(GeV/c) : ',
     >  k*M,'-qz  ',gG, T/1e9, p/1e9,turn,sqrt(bg1/bg),bg1/bg*epsi0,brho
        k=3
        gG = k*M + qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        if(T.ge.T1 .and. T.le.26.e9) 
     >            write(*,*) '  kM+qz gG, T(GeV), p(GeV/c) : ',
     >     k*M,'+qz  ',gG, T/1e9, p/1e9,turn,sqrt(bg1/bg),bg1/bg*epsi0
        k=3
        gG = k*M - qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        if(T.ge.T1 .and. T.le.26.e9) 
     >            write(*,*) '  kM-qz gG, T(GeV), p(GeV/c) : ',
     >     k*M,'-qz  ',gG, T/1e9, p/1e9,turn,sqrt(bg1/bg),bg1/bg*epsi0
        k=4
        gG = k*M + qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        k=4
        gG = k*M - qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        if(T.ge.T1 .and. T.le.26.e9) 
     >            write(*,*) '  kM-qz gG, T(GeV), p(GeV/c) : ',
     >     k*M,'-qz  ',gG, T/1e9, p/1e9,turn,sqrt(bg1/bg),bg1/bg*epsi0
        k=5
        gG = k*M + qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        k=5
        gG = k*M - qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        if(T.ge.T1 .and. T.le.26.e9) 
     >            write(*,*) '  kM-qz gG, T(GeV), p(GeV/c) : ',
     >     k*M,'-qz  ',gG, T/1e9, p/1e9,turn,sqrt(bg1/bg),bg1/bg*epsi0


      write(*,*) ' '
      write(*,*) ' Intrinsic resonances, random, k - qz  (qz=',qz,')'
      do i = -1, 42, 1
        gG = i + qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        if(T.ge.T1 .and. T.le.26.e9) 
     >            write(*,*) ' k+qz, gG, T(GeV), p(GeV/c) : ',
     >    i,'+qz  ', gG, T/1e9, p/1e9,turn,sqrt(bg1/bg),bg1/bg*epsi0
      enddo

      write(*,*) ' '
      write(*,*) ' Intrinsic resonances, random, k - qz  (qz=',qz,')'
      nTurn = 4000
      do i = 4, 59, 1
        gG = i - qz
        gamma = gG / G
        e = gamma * am
        T = e - am
        p = sqrt(T * (T + 2.*am))
        brho = p/c
        bta = p / e
        bg = bta * gamma
        gag = (T+am)/am * G / gG
        turn=(T)/157836.
        if(T.ge.T1 .and. T.le.26.e9) then
          turn1 = turn - nTurn
          T1 = turn1 * 2809.91999999
          p1 = sqrt(T1 * (T1 + 2.*am))          
          brho1 = p1/c
                 write(*,*) ' k-qz, gG, T(GeV), p(GeV/c), turn# : ',
     >    i,'-qz  ', gG, T/1e9, p/1e9,turn,turn1,T1/1e9,BRho1
        endif
      enddo

      stop
      end

