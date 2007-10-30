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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
C-------------------------------------------------------------------------------
C
C                             BIBLIOTHEQUE RPSLD
C 
C                            LECTURE FORMAT LIBRE
C
C                                 J.L.Hamel
C 
C-------------------------------------------------------------------------------
C                
C
C*****EXTRACT: Extraction d'une chaine de caracteres dans TAMP*****
C
      SUBROUTINE EXTRAC(RECORD,IER)
      CHARACTER TAMP*81,RECORD *(*),KR*1
      LOGICAL IDEC,SEPAR,CTRL,UPPER 
      COMMON/EDIT/TAMP
      COMMON/EDIS/ITAMP,IDEC,ITYPE,CTRL,IERLOC,UPPER

      SEPAR(KR)=KR.EQ.' '.OR.((ICHAR(KR).LT.32).AND.(ICHAR(KR).GE.0))

      L=LEN(RECORD)
      RECORD=' '
      IER=-32767
      IF(ITAMP.NE.0) THEN
1        CONTINUE 
         IF(IER.EQ.0) THEN    
            IF(ITAMP.EQ.81) THEN
               I2=81
               ITAMP=0
            ELSE
               IF(SEPAR(TAMP(ITAMP:ITAMP))) THEN 
                  I2=ITAMP-1 
               ELSE             
                  K=ITAMP-I1
                  ITAMP=ITAMP+1
                  IF(K.LT.L) THEN
                     GOTO 1 
                  ELSE
                     I2=ITAMP-1
                  ENDIF
               ENDIF
            ENDIF
         ELSE    
            IF(ITAMP.EQ.81) THEN 
               ITAMP=0
            ELSE
               IF(.NOT.SEPAR(TAMP(ITAMP:ITAMP))) THEN
                  IER=0
                  I1=ITAMP
               ENDIF  
               ITAMP=ITAMP+1
               GOTO 1
            ENDIF
         ENDIF
         IF(IER.EQ.0) RECORD=TAMP(I1:I2)  
      ENDIF
      RETURN
      END
C
C*****RPSLDP:impression de texte et introduction dans la pile*****
C
      SUBROUTINE RPSLDP(NCAR)
      CHARACTER*4 TZERO
C      PARAMETER (TZERO=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0))
      integer*4 itzero
      equivalence (tzero, itzero)
      data itzero/0/
      LOGICAL IDEC,CTRL,UPPER,FLG
      LOGICAL PSLDP
      CHARACTER*81 TAMP,BUF
      CHARACTER*80 ITAB
      COMMON/EDIT/TAMP 
      COMMON/EDIS/ITAMP,IDEC,ITYPE,CTRL,IERLOC,UPPER
      SAVE BUF

      PSLDP=.TRUE.
      CALL PRINT_ONLY(TAMP,NCAR)
      GOTO 100
C
C+++++RPSLDS:impression question et edition de la reponse+++++
C
      ENTRY RPSLDS(NCAR)
      PSLDP=.FALSE.
100   CONTINUE
      CTRL=.FALSE.
C      IERLOC=0 (===> sur demande de Bedfer, pas de RAZ IERLOC)
      IF (.NOT.IDEC .AND. .NOT.PSLDP) THEN
         CALL LINE_EDIT(TAMP,NCAR)
         ITAMP=1
      ELSE IF(IDEC) THEN
         TAMP=BUF
         IF (NCAR.EQ.0) THEN
            I=INDEX(TAMP(1:80),':')
            IF(I.NE.0) THEN
               ITAMP=I+1
            ELSE
               ITAMP=1
            END IF
         ELSE IF (NCAR.GT.0.AND.NCAR.LT.80) THEN
            ITAMP=NCAR
         ELSE
            WRITE(6,'(''*ERREUR DANS RPSLDS (INDEX ABERRANT)*'')')
            STOP
         ENDIF
      ENDIF
      RETURN
C
C+++++RPSUPP:blocage du clavier en majuscules+++++
C
      ENTRY RPSUPP(FLG)
      UPPER=FLG
      RETURN
C
C+++++RPSDEC:entree depuis un tampon en memoire+++++
C
      ENTRY RPSDEC(ITAB)
      IDEC=(ITAB(1:4).NE.TZERO)
      IF (IDEC) THEN
         BUF=ITAB 
         ITAMP=1
      END IF  
      CTRL=.FALSE.
      RETURN

      END
C
C*****RPSLDE:lecture d'un flottant simple precision*****
C
      SUBROUTINE RPSLDE(VAL)
      CHARACTER*20 RECORD

      CALL EXTRAC(RECORD,IER)
      IF(IER.EQ.0) READ(RECORD,*,IOSTAT=IER)VAL
      CALL SAVIER(IER)    
      RETURN
      END      
C
C*****RPSLDD:lecture d'un flottant double precision*****
C
      SUBROUTINE RPSLDD(VAL)
      DOUBLE PRECISION VAL
      CHARACTER*20 RECORD
 
      CALL EXTRAC(RECORD,IER)
      IF(IER.EQ.0) READ(RECORD,*,IOSTAT=IER)VAL
      CALL SAVIER(IER) 
      RETURN
      END   
C
C*****RPSLDI:lecture d'un entier*****
C
      SUBROUTINE RPSLDI(IVAL)
      CHARACTER*20 RECORD

      CALL EXTRAC(RECORD,IER)
      IF(IER.EQ.0) READ(RECORD,*,IOSTAT=IER)IVAL
      CALL SAVIER(IER) 
      RETURN     
      END  
C
C*****RPSLDM:lecture d'un mot de 4 car.*****
C
      SUBROUTINE RPSLDM(VAL)
      CHARACTER RECORD*20
      CHARACTER*4 CVAL
      EQUIVALENCE (VVAL,CVAL)

      CALL EXTRAC(RECORD,IER)
      CVAL=' '
      VAL=VVAL
      IF(IER.EQ.0) CVAL=RECORD(1:4)
      VAL=VVAL
      CALL SAVIER(IER)
      RETURN      
      END
C
C*****RPSLDL:lecture d'un mot de 8car.*****
C
      SUBROUTINE RPSLDL(VAL)
      DOUBLE PRECISION VVAL,VAL
      CHARACTER*8 CVAL
      CHARACTER RECORD*20
      EQUIVALENCE (VVAL,CVAL)

      CALL EXTRAC(RECORD,IER)
      CVAL=' '
      VAL=VVAL
      IF (IER.EQ.0) CVAL=RECORD(1:8)
      VAL=VVAL
      CALL SAVIER(IER)
      RETURN
      END              
C
C*****RPSLDC:lecture d'un complexe simple precision*****
C
      SUBROUTINE RPSLDC(COMPL)
      CHARACTER RECORD*20
      COMPLEX COMPL
      REAL RE,IM

      CALL EXTRAC(RECORD,IER)
      IF (IER.EQ.0) READ(RECORD,*,IOSTAT=IER) RE
      IF (IER.EQ.0) THEN
         CALL EXTRAC(RECORD,IER)
         IF (IER.EQ.0) READ(RECORD,*,IOSTAT=IER) IM
      END IF
      CALL SAVIER(IER)
      COMPL=CMPLX(RE,IM)
      RETURN
      END    
C
C*****RPSLDT:lecture d'un mot de taille quelconque (type CHARACTER)*****
C
      SUBROUTINE RPSLDT(REC)
      CHARACTER*81 RECORD
      CHARACTER REC*(*)

      CALL EXTRAC(RECORD,IER)
      REC=' '
      IF (IER.EQ.0) REC=RECORD
      CALL SAVIER(IER)
      RETURN
      END
C
C*****RPSLDX:lecture d'un reel simple precision ou d'un mot de 4 car.*****
C
      SUBROUTINE RPSLDX(VAL)
      CHARACTER CINI*1,RECORD*20
      CHARACTER*4 CVAL
      LOGICAL IDEC,CTRL,UPPER
      COMMON/EDIS/ITAMP,IDEC,ITYPE,CTRL,IERLOC,UPPER
      EQUIVALENCE (IBLANC,BLANC),(VVAL,CVAL)

      CALL EXTRAC(RECORD,IER)
      CINI=RECORD(1:1)
      IF (((CINI.GE.'0').AND.(CINI.LE.'9'))
     X     .OR.(CINI.EQ.'+').OR.(CINI.EQ.'-').OR.(CINI.EQ.'.')) THEN
         RECORD=' '//RECORD
         IF(IER.EQ.0) READ(RECORD,*,IOSTAT=IER)VAL 
         ITYPE=1  
      ELSE
         CVAL=' '
         VAL=VVAL
         ITYPE=0
         IF(IER.EQ.0) THEN
            CVAL=RECORD(1:4)
            VAL=VVAL
            ITYPE=2
         ENDIF
      ENDIF
      CALL SAVIER(IER)
      RETURN
      END 
C
C*****JRPSLD:lecture du type de donnee lu par RPSLDX*****
C
      SUBROUTINE JRPSLD(JTYPE)
      LOGICAL IDEC,CTRL,UPPER
      COMMON/EDIS/ITAMP,IDEC,ITYPE,CTRL,IERLOC,UPPER

      JTYPE=ITYPE
      RETURN
      END
C
C*****RPSLDZ:lecture d'un entier en hexadecimal*****
C
      SUBROUTINE RPSLDZ(IVAL)
      CHARACTER RECORD*20
      LOGICAL ARRET

      CALL EXTRAC(RECORD,IER)
      IF (IER.EQ.0) THEN
         ITEMP=0
         ARRET=.FALSE.
         I=1
         DO WHILE(.NOT.ARRET .AND. I.LE.8)
            IC=ICHAR(RECORD(I:I))
            I=I+1
            IF (IC.GE.ICHAR('0').AND.IC.LE.ICHAR('9')) THEN
               IC=IC-ICHAR('0')
            ELSE IF (IC.GE.ICHAR('A').AND.IC.LE.ICHAR('F')) THEN
               IC=IC-ICHAR('A')+10
            ELSE IF (IC.GE.ICHAR('a').AND.IC.LE.ICHAR('f')) THEN
               IC=IC-ICHAR('a')+10
            ELSE
               ARRET=.TRUE.
            ENDIF
            IF (.NOT.ARRET) ITEMP=ITEMP*16+IC
         ENDDO
         IF (ARRET.AND.IC.NE.ICHAR(' ')) THEN
            IER=1
         ELSE
            IVAL=ITEMP
         ENDIF
      ENDIF
      CALL SAVIER(IER)
      RETURN
      END    
C
C*****RPSLDH:lecture d'une ligne de 80 car.*****
C
      SUBROUTINE RPSLDH(ITABEX)
      DIMENSION ITABEX(20)
      CHARACTER*81 TAMP
      COMMON/EDIT/TAMP
      LOGICAL IDEC,CTRL,UPPER
      COMMON/EDIS/ITAMP,IDEC,ITYPE,CTRL,IERLOC,UPPER

      READ(TAMP,'(20A4)') ITABEX
      RETURN
      END 
C
C*****SAVIER:sauvegarde de l'erreur*****
C
      SUBROUTINE SAVIER(IER)
      LOGICAL IDEC,CTRL,UPPER
      COMMON/EDIS/ITAMP,IDEC,ITYPE,CTRL,IERLOC,UPPER

      IF(IERLOC.LE.0) THEN
         IF(IER.EQ.-32767) THEN
            IERLOC=IERLOC-1
         ELSE
            IER=IABS(IER)
            IERLOC=IER
         ENDIF
      ENDIF
      IF (.NOT.CTRL .AND. IER.GT.0) THEN
         WRITE(6,*)' *** ERREUR PSLD ',IER,' ***'
C         STOP
      ENDIF
      RETURN
C
C+++++IRPSLD:armement et lecture de l'erreur*****
C
      ENTRY IRPSLD(IER)
      CTRL=.TRUE.
      IER=IERLOC
      IERLOC=0
      RETURN

      END   
C
C*****LINE_EDIT:edition de ligne*****
C
      SUBROUTINE LINE_EDIT(TAMP,CURMIN)
      IMPLICIT INTEGER(A-Z)
      LOGICAL IDEC,CTRL,UPPER,FTSTTY,SBUF
      COMMON/EDIS/ITAMP,IDEC,ITYPE,CTRL,IERLOC,UPPER
      CHARACTER*81 TAMP
      CHARACTER *(*) PRKEY,LINSTA
      CHARACTER*80 BUF,SAVBUF,PFKEY(4)
      CHARACTER*1 CAR,ESC,BKS,CFORM,RC
      CHARACTER*5 WRON,WROFF
      LOGICAL ERREUR,INS,PRONLY,PFKYON
      PARAMETER (ESC='')
      PARAMETER (BKS='')
      PARAMETER (RC='')
C      PARAMETER (LF=CHAR(8))
      PARAMETER (LRING=20)
      PARAMETER (LTRING=LRING/4)
      PARAMETER (MRING=LRING-1)
      PARAMETER (WRON ='[?7h')
      PARAMETER (WROFF='[?7l')
C      CHARACTER*1 CHZERO
C      PARAMETER (CHZERO=CHAR(0))
      CHARACTER*80 RING(0:MRING)
      INTEGER*4 TRING(LTRING)
      EQUIVALENCE (RING, TRING)
      SAVE BUF,CODE,CURC,CURX,INS,RING,NRING,PFKEY,PFKYON
      DATA NRING/0/
      DATA PFKYON/.FALSE./
      DATA TRING/LTRING*0/
      DATA PFKEY/'PF1','PF2','PF3','PF4'/

C+++++ENTREE LINE_EDIT
      PRONLY=.FALSE.
      GOTO 100

C+++++ENTREE PRINT_ONLY
      ENTRY PRINT_ONLY(TAMP,CURMIN)
      PRONLY=.TRUE.

C+++++RECUPERER LA QUESTION A ECRIRE
100   CONTINUE
      CALL FBGTXT
      ERREUR=(CURMIN.LT.0.OR.CURMIN.GT.80)
      IF (.NOT.ERREUR) THEN
         CFORM=TAMP(1:1)
         BUF=TAMP(2:81)
      END IF
C+++++TROUVER LA POSITION INITIALE DU CURSEUR
      IF (.NOT.ERREUR) THEN
         IF (CURMIN.EQ.0) THEN
            I=INDEX(BUF,':')
            IF(I.NE.0) THEN
               CURM=I+1
            ELSE
               CURM=1
            END IF
         ELSE
            CURM=CURMIN
         END IF
C-----ECRITURE DU FORMAT S'IL Y A LIEU
C         IF(CURM.NE.1) THEN
         IF(FTSTTY(0) .AND. .NOT.PRONLY) THEN
            IF(CFORM.EQ.'1') CALL SCRCLR(1,1)
            IF(CFORM.EQ.'0') CALL PUT_CHAR(RC//CHAR(8),0)
         ENDIF
C+++++TRAITEMENT ENTREE PAR FICHIER
         IF(.NOT.FTSTTY(0) .AND. .NOT.PRONLY) THEN
            READ(5,'(A)') TAMP
            CALL PUT_CHAR(RC,0)
            IF(TAMP(1:1).EQ.'#') THEN
               DO WHILE (TAMP(1:1).EQ.'#')
C...pour tests...                  CALL PUT_CHAR(TAMP//RC//CHAR(8),0)
                  READ(5,'(A)') TAMP
               END DO
            END IF
            IF(CURM.EQ.1) THEN
               CALL PUT_CHAR(BUF,0)
               CALL PUT_CHAR(RC//CHAR(8)//TAMP//RC,0)
            ELSE
               CALL PUT_CHAR(BUF(1:CURM-1),0)
               CALL PUT_CHAR(TAMP//RC,0)                
            END IF
            IF(TAMP.EQ.' ') TAMP=BUF(CURM:80)
            RETURN
         ENDIF
C+++++SUITE DU TRAITEMENT INTERACTIF: AFFICHAGE TEXTE 
         WRITE(6,'($)') 
         CALL PUT_CHAR(CHAR(13),1)
         IF (CURM.GT.1) CALL PUT_CHAR(BUF(1:CURM-1),0)
         DO I=CURM,79
            IF(BUF(I:I).EQ.CHAR(9)) BUF(I:I)=' '
         END DO
         CALL PUT_CHAR(BUF(CURM:80),1)
      END IF
C-----PREPARER LE BUFFER DE TRAVAIL
      IF (.NOT.ERREUR) THEN
         BUF=BUF(CURM:80)
         CODE=0
         CURC=1
         CURX=81-CURM
         INS =.FALSE.
         SBUF=.FALSE.
         KR=0
      ENDIF
C+++++SORTIE SI PRINT_ONLY
      IF(PRONLY .AND. .NOT.ERREUR) GOTO 200
C+++++EDITION
      IF (.NOT.ERREUR) THEN
         CALL open_key
         CALL PUT_CHAR(WROFF,0)
         CALL CUR_EXT(BUF(1:CURX),CURC,0)
         DO WHILE (CODE.NE.13.AND.CODE.NE.3)
C-----LECTURE CARACTERE COURANT
            CALL get_c(CAR)
            CODE=ICHAR(CAR)
            IF(CODE.LT.0) CODE=CODE+256
            IF(CODE.EQ.10) CODE=13
            IF (CODE.EQ.27) THEN
               CALL get_c(CAR)
               IF (CAR.EQ.'[') CODE=-1
               IF (CAR.EQ.'O') CODE=-2
               IF (CODE.LT.0)  CALL get_c(CAR)
            ELSE IF (CODE.EQ.3) THEN
               CALL PUT_CHAR(WRON,0)
               CALL close_key
               CALL send_int
            END IF
C-----CARACTERES ORDINAIRES (code ISO 8859 inclus)
            IF ( (CODE.GE.32 .AND. CODE.LT.127) .OR.
     *           (CODE.GE.161) ) THEN
               IF (UPPER .AND. ( (CODE.GE.97.AND.CODE.LE.122) .OR.
     *                           (CODE.GE.224.AND.CODE.LE.253) ) 
     *            ) CODE=CODE-32
               CAR=CHAR(CODE)
               IF (CURC.LE.CURX) THEN
                  CALL PUT_CHAR(CAR,0)
               END IF
               IF (INS.AND.CURC.LE.CURX) THEN
                  I=CURX-CURC+4
                  CALL PUT_CHAR(BUF(CURC:CURX-1),1)
                  BUF=BUF(1:CURC-1)//CAR//BUF(CURC:CURX)
               ELSE 
                  BUF(CURC:CURC)=CAR
               END IF
               IF (CURC.LT.CURX) THEN
                  CURC=CURC+1
               END IF
C-----CURSEUR GAUCHE
            ELSE IF (CAR.EQ.'D'.AND.CURC.GT.1) THEN
               CALL PUT_CHAR(BKS,0)
               CURC=CURC-1
C-----CURSEUR DROIT
            ELSE IF (CAR.EQ.'C'.AND.CURC.LT.CURX) THEN
               CALL PUT_CHAR(ESC//'[C',0)
               CURC=CURC+1
C-----CURSEUR HAUT OU BAS
            ELSE IF (CAR.EQ.'A' .OR. CAR.EQ.'B') THEN
               IF(.NOT.SBUF) THEN
                  SAVBUF=BUF(1:CURX)
                  NRSAV=NRING
                  SBUF=.TRUE.
               ENDIF
               IF(CAR.EQ.'A') THEN
                  IR=-1
                  IF(KR.EQ.0) THEN
                     NR=NRING
                  ELSE
                     NR=MOD(NRING+LRING-1,LRING)
                  ENDIF
               ELSE
                  IR=1
                  NR=MOD(NRING+1,LRING)
               ENDIF
               IF(RING(NR)(1:1).NE.CHAR(0) .OR. KR.EQ.-IR) THEN
                  KR=KR+IR
                  IF(KR.EQ.0) THEN
                     BUF(1:CURX)=SAVBUF
                     NRING=NRSAV
                  ELSE
                     NRING=NR
                     BUF(1:CURX)=RING(NR)
                  ENDIF
                  CALL AFF_LINE(BUF(1:CURX),CURC)
               ENDIF
C-----INSERT/OVER
            ELSE IF (CODE.EQ.1) THEN
               INS=.NOT.INS
C-----DELETE
            ELSE IF (CODE.EQ.127.AND.CURC.GT.1) THEN
               I=CURX-CURC+7
               CALL PUT_CHAR(BKS,0)  
               CALL PUT_CHAR(BUF(CURC:CURX)//' ',1)
               BUF=BUF(1:CURC-2)//BUF(CURC:CURX)//' '
               CURC=CURC-1
C-----CTRL H (DEBUT DE LIGNE)
            ELSE IF (CODE.EQ.8) THEN
               CALL CUR_EXT(BUF(1:CURX),CURC,0)
C-----CTRL E (FIN DE LIGNE)
            ELSE IF (CODE.EQ.5) THEN
               CALL CUR_EXT(BUF(1:CURX),CURC,1)
C-----TAB
            ELSE IF (CODE.EQ.9 .AND. INDEX(BUF(CURC:),' ').NE.0) THEN
               CALL CUR_TAB(BUF,CURC)
C-----CLEFS PF1 A PF4
            ELSE IF (PFKYON .AND. CAR.GE.'P' .AND. CAR.LE.'S' 
     *                          .AND. CODE.EQ.-2) THEN
                I=ICHAR(CAR)-ICHAR('P')+1
                BUF(1:CURX)=PFKEY(I)
                CALL AFF_LINE(BUF(1:CURX),CURC)
                CODE=13
C-----PF1 (DELETE SOUS LE CURSEUR)
C            ELSE IF (CODE.EQ.256.AND.CURC.LT.CURX) THEN
C               I=CURX-CURC+5
C               CALL PUT_CHAR(
C     1            BUF(CURC+1:CURX)//' ',1)
C               BUF=BUF(1:CURC-1)//BUF(CURC+1:CURX)//' '
            END IF
         END DO
         IF (CODE.NE.3) THEN
            CALL PUT_CHAR(WRON,0) 
            CALL close_key
         END IF
      END IF
C+++++TRAITEMENT ERREUR
      IF (ERREUR) THEN
         WRITE(6,'('' *ERREUR DANS LINE_EDIT*'')')
         STOP
      END IF
C+++++FIN DE LINE_EDIT
200   CONTINUE
      CALL put_c(CHAR(13))
      WRITE(6,'()')
      TAMP=BUF(1:CURX)
C.....Si le tampon envoye est dans la pile, il y reste...
C.....Sinon le nouveau tampon est mis en pile
      IF(KR.EQ.0) NRING=MOD(NRING+1,LRING)
      RING(NRING)=TAMP
C.....Si les cles PF sont armees, les desarmer, et RAZ ligne d'etat
      IF(PFKYON) THEN
         CALL WRSTAT(' ')
         PFKYON=.FALSE.
      ENDIF

      RETURN

C+++++ENTREE INPFKY:programmation des cles PF
      ENTRY INPFKY(NKEY,PRKEY)
      PFKEY(NKEY)=PRKEY
      RETURN

C+++++ENTREE ENPFKY:armement des cles PF et ecriture ligne d'etat
      ENTRY ENPFKY(LINSTA)
      PFKYON=.TRUE.
      CALL WRSTAT(LINSTA)
      RETURN

      END
C
C*****AFF_LINE:affichage d'une ligne de la pile*****
C
      SUBROUTINE AFF_LINE(TEXT,CURS)
      INTEGER CURS
      CHARACTER *80 TEXT
      CHARACTER*5 TEMP

      IF(CURS.GT.1) THEN
         WRITE(TEMP,'(A2,I2,A1)') CHAR(27)//'[',CURS-1,'D'
         IF(TEMP(3:3).EQ.' ') TEMP(3:3)='0'
         CALL PUT_CHAR(TEMP//TEXT,1)
      ELSE
         CALL PUT_CHAR(TEXT,1)
      ENDIF
      RETURN
      END
C
C*****CUR_EXT:positionnement curseur sur 1er ou dernier item de la ligne*****
C
      SUBROUTINE CUR_EXT(TEXT,CURS,SENS)
      INTEGER CURS,SENS
      CHARACTER *(*) TEXT
      CHARACTER *1 C
      CHARACTER *5 TEMP

C+++++TRAITEMENT DEBUT ET FIN DE LIGNE
      L=LEN(TEXT)
      IF(SENS.EQ.0) THEN
         I=1
         DO WHILE(TEXT(I:I).EQ.' ' .AND. I.NE.L)
            I=I+1
         END DO
         IF(I.EQ.L) I=1
      ELSE
         I=L
         DO WHILE(TEXT(I:I).EQ.' ' .AND. I.NE.1)
            I=I-1
         END DO
         IF(I.EQ.1) THEN
            I=L
         ELSE
            DO WHILE(TEXT(I:I).NE.' ' .AND. I.NE.1)
               I=I-1
            END DO
            IF(TEXT(I:I).EQ.' ') I=I+1
         ENDIF
      ENDIF
      GOTO 100

C+++++TRAITEMENT TAB
      ENTRY CUR_TAB(TEXT,CURS)
      L=LEN(TEXT)
      I=CURS
      IF(TEXT(I:I).NE.' ') THEN
         J=INDEX(TEXT(I:L),' ')
         IF(J.EQ.0) RETURN
         I=I+J-1
      ENDIF
      DO WHILE(TEXT(I:I).EQ.' ' .AND. I.LT.L)
         I=I+1
      END DO
      IF(TEXT(I:I).EQ.' ') RETURN

C-----POSITIONNEMENT DU CURSEUR
100   CONTINUE
      IF(I.EQ.CURS) RETURN
      N=CURS-I
      IF(N.GT.0) THEN
         C='D'
      ELSE
         C='C'
         N=-N
      ENDIF
      WRITE(TEMP,'(A2,I2,A1)') CHAR(27)//'[',N,C
      IF(TEMP(3:3).EQ.' ') TEMP(3:3)='0'
      CALL PUT_CHAR(TEMP,0)
      CURS=I
      RETURN
      END
C
C*****PUT_CHAR:sortie de car. avec restauration du curseur*****
C
      SUBROUTINE PUT_CHAR(CHAINE,RESCUR)
      IMPLICIT INTEGER (A-Z)   
      CHARACTER *(*) CHAINE
      CHARACTER*1 BKS,ESC
      CHARACTER*2 SCUR,RCUR
      PARAMETER (BKS='')
      PARAMETER (ESC='')
      PARAMETER (SCUR='7')
      PARAMETER (RCUR='8')

      IF (RESCUR.NE.0) THEN
          CALL put_c(ESC)
          CALL put_c('7')
      ENDIF

      L=LEN(CHAINE)
      DO I=1,L
         CALL put_c(CHAINE(I:I))
      END DO

      IF (RESCUR.NE.0) THEN
          CALL put_c(ESC)
          CALL put_c('8')
      ENDIF

C         DO I=1,L
C            CALL put_c(BKS)
C         END DO

      RETURN
      END

C--------------------------------------------------------C
C                  FUNCTION IDLG                         C
C                                                        C
C   ARGUMENTS:                                           C
C                                                        C
C     -FORME:TEXTE A IMPRIMER ,PASSE DANS UNE VARIABLE   C
C            OU EN CONSTANTE                             C
C            Format: (' texte:'[)]                       C
C                   de 80 caracteres au plus             C
C                                                        C
C     -ITAB :CHAINE DE CARACTERES                        C
C            CONTENANT LES REPONSES (DE 4 CAR.)          C
C                                                        C
C     -NREP :DIMENSION DU TABLEAU ITAB                   C
C                                                        C
C   RETOUR   :                                           C
C                                                        C
C     -IDLG :RANG DANS LE TABLEAU ITAB DE LA REPONSE;    C
C            PAR DEFAUT,LA VALEUR 1                      C
C--------------------------------------------------------C
      FUNCTION IDLG(FORME,ITAB,NREP)
      CHARACTER * (*) ITAB,FORME
      CHARACTER * 4 REP
      CHARACTER*81 TAMP,BUF
      LOGICAL IDLGDEF
      COMMON/EDIT/TAMP

      BUF=TAMP
      IF (FORME(1:1).NE.CHAR(0)) WRITE(BUF,FMT=FORME)
1     CONTINUE
      TAMP=BUF
      CALL RPSLDS(0)
      BUF(1:1)=' '
      CALL RPSLDM(XREP)
      WRITE(REP,'(A4)') XREP
      IF(NREP.GE.1 .AND. .NOT.(IDLGDEF() .AND. REP.EQ.'    ')) THEN
         DO 2 I=1,NREP
            J=(I-1)*4+1
            IF(REP.EQ.ITAB(J:J+3))GOTO 3
2        CONTINUE
         CALL PUT_CHAR(CHAR(27)//'[A',0)
         GOTO 1
      ELSE
         I=1
      ENDIF 
3     CONTINUE
      IDLG=I
      RETURN
      END
 
C--------------------------------------------------------C
C                  FUNCTION IDLGA                        C
C                                                        C
C   ARGUMENTS:                                           C
C                                                        C
C     -FORME:TEXTE A IMPRIMER ,PASSE DANS UNE VARIABLE   C
C            OU EN CONSTANTE                             C
C            Format: (' texte:'[)]                       C
C                   de 80 caracteres au plus             C
C                                                        C
C     -ITAB :TABLEAU D'INTEGER OU DE REAL                C
C            CONTENANT LES REPONSES (DE 4 CAR.)          C
C                                                        C
C     -NREP :DIMENSION DU TABLEAU ITAB                   C
C                                                        C
C   RETOUR   :                                           C
C                                                        C
C     -IDLG :RANG DANS LE TABLEAU ITAB DE LA REPONSE;    C
C            PAR DEFAUT,LA VALEUR 1                      C
C--------------------------------------------------------C
      FUNCTION IDLGA(FORME,ITAB,NREP)
      CHARACTER * (*) FORME
      DIMENSION ITAB(1)
      EQUIVALENCE (XREP,IREP)
      CHARACTER*4 REP
      LOGICAL IDLGDEF
      CHARACTER*81 TAMP,BUF
      COMMON/EDIT/TAMP
   
      BUF=TAMP
      IF (FORME(1:1).NE.CHAR(0)) WRITE(BUF,FMT=FORME)
1     CONTINUE
      TAMP=BUF
      CALL RPSLDS(0)
      BUF(1:1)=' '
      CALL RPSLDM(XREP)
      WRITE(REP,'(A4)') XREP
      IF(NREP.GE.1 .AND. .NOT.(IDLGDEF() .AND. REP.EQ.'    ')) THEN
         DO 2 I=1,NREP
            IF(IREP.EQ.ITAB(I))GOTO 3
2        CONTINUE
         CALL PUT_CHAR(CHAR(27)//'[A',0)
         GOTO 1
      ELSE
         I=1
      ENDIF 
3     CONTINUE
      IDLGA=I
      RETURN
      END
C--------------------------------------------------------C
C   LOGICAL FUNCTION IDLGDEF                             C
C   TESTE SI 1ERE REPONSE PAR DEFAUT                     C
C                                                        C
C   SUBROUTINE DLGDEF(LOGICAL_VARIABLE)                  C
C   ARME/DESARME 1ERE REPONSE PAR DEFAUT                 C
C--------------------------------------------------------C

      LOGICAL FUNCTION IDLGDEF()
      LOGICAL DEFCOU, DEF
      SAVE DEFCOU
      DATA DEFCOU/.FALSE./

      IDLGDEF = DEFCOU
      RETURN

      ENTRY DLGDEF(DEF)
      DEFCOU = DEF
      RETURN

      END
