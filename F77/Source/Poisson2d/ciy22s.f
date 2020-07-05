C
C   Copyright (c) 1997 Silvano Bonazzola
C
C    This file is part of LORENE.
C
C    LORENE is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    LORENE is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with LORENE; if not, write to the Free Software
C    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
       SUBROUTINE CIY22S(NDEG,NDIMY,NN64,ITCH,IDR,IND,C64,CC,CS,DEN)

		implicit double precision (a-h,o-z)

C
C           ROUTINE POUR LE CALCUL DES TRANSFORMEES DE FOURIER
C           E DE TCHEBYTCHEV D'UN TABLEAU A 2 DIMENSIONS. CETTE
C           ROUTINE A ETE' SPECIALEMENT REALISEE POUR ETRE APPELLE
C           PAR FCIR2S ET FGIR2S. IL S'AGIT PARCONSEQUENT D'UNE
C           ROUTINE HAUTEMENT SPECIALISEE.
C
C
C           ROUTINE COMPLETEMENT CRAYTINIZEE.
C
C           ARGUMENTS DE LA ROUTINE:
C
C               NDEG = TABLEAU, NDEG(2) CONTENANT LES DE-
C                    GREES DE LIBERTE DES TRANSFORMEES
C                    A EFFECTUER, NDEG(1) CONCERNE LE PRE-
C                    MIER INDICE DE LA MATRICE, NDEG(2)
C                    LE 2ME INDICE DE LA MATRICE.
C                    MATRICE.     
C                    NDEG DOIT IMPERATIVEMENT ETRE DE LA
C                    FORME 2**m*3**p*5**q+1 (m,p,q NOMBRES 
C                    INTIERS).
C
C           NDIMY   =DIMESION DU TABLEAU YY(LR,LY,LZ).
C           POUR DES RAISONS DE CRAYTINISATION NDIMY ET NDIMZ NE DOIT PA
C  S
C           ETRE UN MULTIPLE DE 8.
C
C           NN64 = PARAMETRE DE LA VECTORIZATION, PAR EXEMPLE
C                NN64=64 SIGNIFIE QUE 64 FONCTIONS A TRANSFORMER
C                SONT VECTORIZEE.
C
C           ITCH =PARAMETRE, SII ITCH.EQ.1 LA TRANSFORMEE
C                DE FOURIER EST EFFECTUEE, SI ITCH=2 LA SUBROU-
C                TINE EFFECTUE LA TRANSFORMEE DE TCHEBYTCHEV.
C
C
C           IDR =PRAMETRE: VOIR SUBROUTINE FCIR2S
C
C           IN1 = PARAMETRE, SI IN1=1 LA TRANSFORMEE EST
C                EFFECTUEE SUR LE PREMIER INDICE, SI IN1=2 SUR LE
C                DEUXIEME.
C                PAR EXEMPLE ITCH=2, IN1=1 SIGNIFIE
C                QU'ON EFFECTUE LA TRANSFORMATION DE TCHEBY-
C                TCHEV SR LE PREMIER INDICE DU TABLEAU YY.
C
C           IND = PARAMETRE: EN OUTPUT ON A LES COEFFICIENTS
C                 DE FOURIER (OU DE TCHEBYTCHEV SELON ITCH) DES DERI-
C                 VEES (PREMIERES OU DEUXIEMES SELON IDR) SI IND=1.
C                 SI IND=2 ON A LES DERIVEES DES FONCTIONS.
C
C           C64,CC,CS= TABLEAUX DE TRAVAIL: DIMENSION MINIME=
C                  (NN64+1)*((MAX(NDEG(1),NDEG(2))+3)
C
C           DEN =TABLEAU A 2 DIMENSIONS CONTENANT LA FONCTION
C                A TRANSFORMER EN IMPUT, ET LA TRANSFORMEE EN
C                OUTPUT. 
C
C           ROUTINE ayant testee avec le protocol ordinaire le 10/1/1987
C
C
C $Id: ciy22s.f,v 1.2 2012/03/30 12:12:43 j_novak Exp $
C $Log: ciy22s.f,v $
C Revision 1.2  2012/03/30 12:12:43  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.1  1997/10/23  08:39:49  eric
c Initial revision
c
C
C $Header: /cvsroot/Lorene/F77/Source/Poisson2d/ciy22s.f,v 1.2 2012/03/30 12:12:43 j_novak Exp $
C
C
	character*120 header
	data header/'$Header: /cvsroot/Lorene/F77/Source/Poisson2d/ciy22s.f,v 1.2 2012/03/30 12:12:43 j_novak Exp $'/


       DIMENSION DEN(NDIMY,*),C64(*),CC(*),CS(*),NDEG(3)
       DATA NDY,NDZ,NEQ,INITI/0,0,0,0/

	save	NDY,NDZ,NEQ,NZ,NY,NY64,MULTY,NY65,NYR64
	save	NYR65

C
C           INITIALISATION.
C
       NY1=NDEG(1)
       NZ1=NDEG(2)
C
       IF(NDY.EQ.NY1.AND.NDZ.EQ.NZ1.AND.NEQ.EQ.NN64) GO TO 300
C
C           DETERMINATION DES PARAMETRES NECESSAIRES POUR LA TRANSFORMAT
C  ION
C           DU 2ME INDICE.
C
       NDY=NY1
       NDZ=NZ1
       NEQ=NN64
       NZ=NZ1-1
       NY=NY1-1
C
       NY64=NN64
       MULTY=NY1/NY64
       IF(MULTY.EQ.0) THEN 
       NY64=NY1
       MULTY=1
       ENDIF
C
       NY65=NY64
       IF((NY64/8)*8.EQ.NY64) NY65=NY64+1
       NYR64=NY1-MULTY*NY64
       NYR65=NYR64
       IF((NYR64/8)*8.EQ.NYR64)NYR65=NYR64+1
  300  CONTINUE
C
C
CC2222222222222222222222222222222222222222222222222222222222222222222222
C  2222
CC2222222222222222222222222222222222222222222222222222222222222222222222
C  222.
C22222222222222222222222222222222222222222222222222222222222222222222222
C  2.
C
       IF(NZ1.LT.2) RETURN
       MM=ITCH-1
C
C      CALCUL DES COEFFICIENTS DE LA DERIVEE 2ME
C
       IF(IDR.EQ.2) THEN
C
       DO 507 LY=1,NY1
       DEN(LY,1)=0
       DEN(LY,2)=-DEN(LY,2)
  507  CONTINUE              
C
       DO 509 LZ=3,NZ1
       DE2=-(LZ-1)**2
       DO 511 LY=1,NY1
       DEN(LY,LZ)=DEN(LY,LZ)*DE2
  511  CONTINUE              
  509  CONTINUE              
       IF(IND.EQ.1)RETURN
       ENDIF
C
C           CALCUL DES COEFFICIENTS DE LA DERIVEE PREMIERE
C
       IF(IDR.EQ.7) THEN
       IF(ITCH.EQ.1) THEN
       DO 95 LZ=1,NZ1
       FACT=-(LZ-1)
       DO 96 LY=1,NY1
       DEN(LY,LZ)=DEN(LY,LZ)*FACT
  96   CONTINUE
  95   CONTINUE
C
       IF(IND.EQ.1) RETURN
       ENDIF
C
       IF(ITCH.EQ.2) THEN
       DO 97 LZ=1,NZ1
       FACT=LZ-1
       DO 98 LY=1,NY1
       DEN(LY,LZ)=DEN(LY,LZ)*FACT
  98   CONTINUE
  97   CONTINUE
       ENDIF
       IF(IND.EQ.1) RETURN
       ENDIF
C
C           CAS NN64.GE.NY1
C
             IF(NY64.EQ.NY1) THEN
C
C
       IF(IDR.LT.5) THEN
       IDR1=IDR-2
       IDR2=IDR
       IF(IDR.GT.2) IDR2=1
       ENDIF
C
       LN65=0
       DO 66 LZ=1,NZ1
       DO 67 LY=1,NY64
       CC(LY+LN65)=DEN(LY,LZ)
  67   CONTINUE                    
       LN65=LN65+NY65
  66   CONTINUE                    
C
       IF(IDR.EQ.0.AND.ITCH.EQ.1) CALL CHIYMS(NZ,NY64,CC,CS,C64)
       IF(IDR.EQ.0.AND.ITCH.EQ.2) CALL CIY2MS(NZ,NY64,CC,CS,C64)
C
       IF(IDR.EQ.1) THEN
       IF(ITCH.EQ.1) CALL DECFMS(NZ,NY64,CC,CS,C64)
       IF(ITCH.EQ.2) CALL DESFMS(NZ,NY64,CC,CS,C64)
       ENDIF
C
       IF(IDR.EQ.2) THEN
       IF(ITCH.EQ.1) CALL CHIYMS(NZ,NY64,CC,CS,C64)
       IF(ITCH.EQ.2) CALL CIY2MS(NZ,NY64,CC,CS,C64)
       ENDIF
C
       IF(IDR.EQ.3.OR.IDR.EQ.4) THEN
       CALL DISEMS(NZ,NY64,IDR1,ITCH,CC,CS,C64)
       ENDIF
C
       IF(IDR.EQ.5.OR.IDR.EQ.6) CALL CHELES(NZ,ITCH,NY64,MM,CC,CS,C64)
       IF(IDR.EQ.8) CALL DIS2MS(NZ,NY64,CC,CS,C64)
       IF(IND.EQ.2.AND.IDR2.EQ.1)THEN
       IF(ITCH.EQ.1) CALL CHIYMS(NZ,NY64,C64,CS,CC)
       IF(ITCH.EQ.2) CALL CIY2MS(NZ,NY64,C64,CS,CC)
C
       LN65=0
       DO  68 LZ=1,NZ1
       DO  69 LY=1,NY64
       DEN(LY,LZ)=CC(LY+LN65)
  69   CONTINUE                    
       LN65=LN65+NY65
  68   CONTINUE                    
       RETURN
       ENDIF
C
       LN65=0
       DO 70 LZ=1,NZ1
       DO 71 LY=1,NY64
       DEN(LY,LZ)=C64(LY+LN65)
  71   CONTINUE                          
       LN65=LN65+NY65
  70   CONTINUE                          
       RETURN
       ENDIF
C
C           TRANSFORMATION DANS LE CAS NN64< NY1
       L164=1
       L264=NY64
       LL64=0
       DO 83 LMUL=1,MULTY
       LN65=LL64
C
       LN65=LL64
       DO 77 LZ=1,NZ1
       DO 78 LY=L164,L264
       CC(LY+LN65)=DEN(LY,LZ)
  78   CONTINUE                    
       LN65=LN65+NY65
  77   CONTINUE                    
C
       IF(IDR.EQ.0.AND.ITCH.EQ.1) CALL CHIYMS(NZ,NY64,CC,CS,C64)
       IF(IDR.EQ.0.AND.ITCH.EQ.2) CALL CIY2MS(NZ,NY64,CC,CS,C64)
C
       IF(IDR.EQ.1) THEN
       IF(ITCH.EQ.1) CALL DECFMS(NZ,NY64,CC,CS,C64)
       IF(ITCH.EQ.2) CALL DESFMS(NZ,NY64,CC,CS,C64)
       ENDIF
C
       IF(IDR.EQ.2) THEN
       IF(ITCH.EQ.1) CALL CHIYMS(NZ,NY64,CC,CS,C64)
       IF(ITCH.EQ.2) CALL CIY2MS(NZ,NY64,CC,CS,C64)
       ENDIF
C
       IF(IDR.EQ.3.OR.IDR.EQ.4) THEN
       CALL DISEMS(NZ,NY64,IDR1,ITCH,CC,CS,C64)
       ENDIF
C
       IF(IDR.EQ.5.OR.IDR.EQ.6) CALL CHELES(NZ,ITCH,NY64,MM,CC,CS,C64)
       IF(IDR.EQ.8) CALL DIS2MS(NZ,NY64,CC,CS,C64)
       IF(IND.EQ.2.AND.IDR2.EQ.1)THEN
       IF(ITCH.EQ.1) CALL CHIYMS(NZ,NY64,C64,CS,CC)
       IF(ITCH.EQ.2) CALL CIY2MS(NZ,NY64,C64,CS,CC)
C
       LN65=LL64
       DO 79 LZ=1,NZ1
       DO 80 LY=L164,L264
       DEN(LY,LZ)=CC(LY+LN65)
  80   CONTINUE                    
       LN65=LN65+NY65
  79   CONTINUE                    
       GO TO 801
       END IF
C
       LN65=LL64
C
       DO 81 LZ=1,NZ1
       DO 82 LY=L164,L264
       DEN(LY,LZ)=C64(LY+LN65)
  82   CONTINUE                    
       LN65=LN65+NY65
  81   CONTINUE                    
C
  801  CONTINUE
       L164=L264+1
       L264=L264+NY64
       LL64=-L164+1
  83   CONTINUE                    
C
C           LE CALCUL EST TERMINE' SI NY1 EST UN MLTIPLE DE NN64
C
             IF(NYR64.EQ.0) RETURN
C
C           FIN DU CALCUL SI NY1 N'EST PAS UN MULTIPLE DE NN64
C
       LN65=LL64
        DO 89 LZ=1,NZ1
       DO 90 LY=L164,NY1
        CC(LY+LN65)=DEN(LY,LZ)
  90   CONTINUE                    
       LN65=LN65+NYR65
  89   CONTINUE                    
C
       IF(IDR.EQ.0.AND.ITCH.EQ.1) CALL CHIYMS(NZ,NYR64,CC,CS,C64)
       IF(IDR.EQ.0.AND.ITCH.EQ.2) CALL CIY2MS(NZ,NYR64,CC,CS,C64)
C
       IF(IDR.EQ.1) THEN
       IF(ITCH.EQ.1) CALL DECFMS(NZ,NYR64,CC,CS,C64)
       IF(ITCH.EQ.2) CALL DESFMS(NZ,NYR64,CC,CS,C64)
       ENDIF
C
       IF(IDR.EQ.2) THEN
       IF(ITCH.EQ.1) CALL CHIYMS(NZ,NYR64,CC,CS,C64)
       IF(ITCH.EQ.2) CALL CIY2MS(NZ,NYR64,CC,CS,C64)
       ENDIF
C
       IF(IDR.EQ.3.OR.IDR.EQ.4) THEN
       CALL DISEMS(NZ,NYR64,IDR1,ITCH,CC,CS,C64)
       ENDIF
C
       IF(IDR.EQ.5.OR.IDR.EQ.6) CALL CHELES(NZ,ITCH,NYR64,MM,CC,CS,C64)
       IF(IDR.EQ.8) CALL DIS2MS(NZ,NYR64,CC,CS,C64)
       IF(IND.EQ.2.AND.IDR2.EQ.1)THEN
       IF(ITCH.EQ.1) CALL CHIYMS(NZ,NYR64,C64,CS,CC)
       IF(ITCH.EQ.2) CALL CIY2MS(NZ,NYR64,C64,CS,CC)
C
       LN65=LL64
       DO 91 LZ=1,NZ1
       DO 92 LY=L164,NY1
       DEN(LY,LZ)=CC(LY+LN65)
  92   CONTINUE                    
       LN65=LN65+NYR65
  91   CONTINUE                    
       RETURN
       END IF
C
       LN65=LL64
       DO 93 LZ=1,NZ1
       DO 94 LY=L164,NY1
       DEN(LY,LZ)=C64(LY+LN65)
  94   CONTINUE                    
       LN65=LN65+NYR65
  93   CONTINUE                    
C
C      FIN DE LA TRANSFORMATION DES COLONNES DE LA MATRICE DEN
C      ......................................................
C
  100  FORMAT(1X,'FOU',9D12.4)
  101  FORMAT(1X,'FOU')
  200  FORMAT(1X,'FOU',20I5)
       RETURN
       END
C
