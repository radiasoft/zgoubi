
      implicit double precision (a-h,o-z)

      character(200) txt200
      parameter (mxa = 600, mxr = 400, mxz = 100)
      dimension aa(mxa), rr(mxr), zz(mxz)
      dimension ba(mxa,mxr,mxz),br(mxa,mxr,mxz),bz(mxa,mxr,mxz)
      parameter (rad = (4.d0 * atan(1.d0)) /180.d0)

      open(unit=1,file='zgoubi.impdev.out')
      open(unit=2,file='impdev2FieldMap.out')

C read 2-line header
          read(1,fmt='(a)')
          read(1,fmt='(a)')

      AT = 30.d0  * rad    ! rad
      ia = 567    ! number of mesh nodes, longitudinal
      dtta = AT/dble(ia-1)   ! rad

C Filed map extends ~ 440 - 520 cm
      dz = .1d0 ;  dR = 0.4 ! cm
      dz = .5d0 ;  dR = 1. ! cm
      dz = 1.d0 ;  dR = 2. ! cm

      jr = 201 ; kz = 41  ! has to be odd
      jr = 81 ; kz = 9  ! has to be odd
      jr = 61 ; kz = 5  ! has to be odd

      write(*,*) 'Number of mesh nodes in X/Y/Z : ',ia,'/',jr,'/',kz
      write(*,*) 'type enter to continue'
      read(*,*)
      write(*,*) 'Ok, now busy, making map ... wait'
      write(*,*) 'Map will be in impdev2FieldMap.out'

      if(ia .gt. mxa) stop ' Too many angles'
      if(jr .gt. mxr) stop ' Too many radii'
      if(kz .gt. mxz) stop ' Too many Zs'

      do i = 1, ia            
        aa(i) = dble(i-1) * dtta
      enddo

      bami = 1d10
      bama = -1d10
      brmi = 1d10
      brma = -1d10
      bzmi = 1d10
      bzma = -1d10

      do j = 1, jr
        jj = (jr+1)/2 + j/2*(-1)**j 
        do k = 1, kz
          kk = (kz+1)/2 + k/2*(-1)**k
          do i = 1, ia
             read(1,*,err=10,end=10) D, Y, T, Z, P, X, S, BX, BY, BZZ
             if(i.eq.ia) read(1,*,err=10,end=10) ! because integr writes twice at the last integration step
            rr(jj) = Y
            zz(kk) = Z

            ba(i,jj,kk) = BX
            br(i,jj,kk) = BY
            bz(i,jj,kk) = BZZ

            if(bami .gt. bx) bami = bx
            if(bama .lt. bx) bama = bx
            if(brmi .gt. by) brmi = by
            if(brma .lt. by) brma = by
            if(bzmi .gt. bzz) bzmi = bzz
            if(bzma .lt. bzz) bzma = bzz

c             write(*,*) i,jj,kk,aa(i),rr(j),zz(k),
c     >         ba(i,jj,kk), br(i,jj,kk), bz(i,jj,kk)

          enddo
        enddo
      enddo

 10   continue

      write(*,*) ' Done reading ! '

      RM = rr(1 + jr/2)    ! cm
      xpas = dtta * RM   ! cm
      write(*,*) ' Mesh step at RM = ',RM,' is ',xpas,' cm'

      write(*,*) ' Min - max field components :  '
      write(*,*) '      ba :  ',bami,'/',bama
      write(*,*) '      br :  ',brmi,'/',brma
      write(*,*) '      bz :  ',bzmi,'/',bzma

      Rmin = rr(1)
      write(2,fmt='(f8.2,f8.2,f20.14,f8.2,a)') Rmin , DR, DTTA/rad, DZ,
     >'       ! R_min (cm), DR (cm), DTTA (deg), DZ (cm) '
      write(2,fmt='(a)') '# Field map generated using impdev2FieldMap'
      write(2,fmt='(a)') 
     >'# AT/deg   RM/cm   xpas/cm  dR/cm  dZ/cm  ia jr  kz : '
      write(2,fmt='(a,1p,5(e14.6,1x),3(i5,1x))') 
     >'#', AT/((4.d0 * atan(1.d0)) /180.d0), RM,xpas,dR,dZ,ia,jr,kz
      write(2,fmt='(a)') 
     >'#  theta/rad    R/cm    Z/cm      B_theta   B_R    B_Z'

      do j = 1, jr
        do k = 1, kz
          do i = 1, ia            
            write(2,fmt='(1p,e18.10,1x,
     >      0p,2(f11.5,1x), 
     >      1p,3(e18.10,1x),3(i5,1x))') aa(i),rr(j),zz(k),
     >      ba(i,j,k), br(i,j,k), bz(i,j,k),i,j,k
          enddo
        enddo
      enddo

      stop
      end
















