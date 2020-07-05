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
       SUBROUTINE CERX3S(NDEG,NDIMR,NDIMY,NNN64,ITCH,C64,CC,CS,DEN)

		implicit double precision (a-h,o-z)

C
C           ROUTINE POUR LE CALCUL DES TRANSFORMEES DE FOURIER
C           E DE TCHEBYTCHEV D'UN TABLEAU A 3 DIMENSIONS POUR LA PARTIE
C           RADIALE DES FONCTIONS EN COORDONNES SPHERIQUES.
C           CETTE ROUTINE EST HAUTEMENT SPECIALISEE' ET DOIT ETRE
C           APPELLEE PAR FCE3S OU FGE3RS.
C
C
C           ROUTINE COMPLETEMENT CRAYTINIZEE.
C
C           ARGUMENTS DE LA ROUTINE:
C
C               NDEG = TABLEAU, NDEG(3) CONTENANT LES DE-
C                    GREES DE LIBERTE DES TRANSFORMEES
C                    A EFFECTUER, NDEG(1) CONCERNE LE PRE-
C                    MIER INDICE DE LA MATRICE, NDEG(2)
C                    LE 2ME INDICE, NDEG(3) LE 3me INDICE DE LA
C                    MATRICE.     
C                    NDEG DOIT IMPERATIVEMENT ETRE DE LA
C                    FORME 2**m*3**p*5**q POUR LES TRANS-
C                    FORMEEES DE FOURIER (m,p,q NOMBRES 
C                    INTIERS)
C                    ET 2**p*3**p*5*q+1 POUR LES TRANSFOR-
C                    MEES DE TCHEBYTCHEV.
C
C           NDIMY,NDIMZ  =DIMESION DU TABLEAU YY(LR,LY,LZ).
C           POUR DES RAISONS DE CRAYTINISATION NDIMY ET NDIMZ NE DOIT PA
C  S
C           ETRE UN MULTIPLE DE 8.
C
C           NNN64 = PARAMETRE DE LA VECTORIZATION, PAR EXEMPLE
C                NNN64=64 SIGNIFIE QUE 64 FONCTIONS A TRANSFORMER
C                SONT VECTORIZEE.
C
C           ITCH =PARAMETRE:IL DOIT ETRE=0 SI LA FONCTION A TRANSFORMER
C                EST SYMMETRIQUE PAR RAPPORT r=0, = 1 DANS LE CAS
C                CONTRAIRE.
C
C           C64,CC,CS= TABLEAUX DE TRAVAIL: DIMENSION MINIME=
C                  (NNN64+1)*((MAX(NDEG(1),NDEG(2))+3)
C
C           DEN =TABLEAU A 3 DIMENSIONS CONTENANT LA FONCTION
C                A TRANSFORMER EN IMPUT, ET LA TRANSFORMEE EN
C                OUTPUT. LES COEFFICIENTS DE FOURIER (DANS LE CAS
C                D'UNE TRANSFORMATION DE FOURIER) SONT STOCKES
C                DANS LA FACON SUIVANTE: (EXEMPLE DANS LE CAS
C                DE LA TRANSFORMATION DU 1ER INDICE) DANS DEN(1,L,M)
C                IL-Y-A LE COEFFICIENT DE FOURIER COSINUS DE LA 
C                FREQUENCE ZERO, DANS LES COEFFICIENTS PAIRES 
C                (DEN(2,L,M),DEN(4,L,M)...DEN(2*N,L,M) LES COEFI-
C                CIENTS EN COSINUS, DANS LES COEFFICIENTS IMPAIRES
C                DEN(3,L,M),DEN(5,L,M),...DEN(2*N-1) LES TERMES EN 
C                SINUS. PARCONSEQUENT LE TERMES DEN(2,L,M) ET DEN(3,L,M)
C                SONT LES TERMES DU DEVELOPPEMENT DE FOURIER CORRESPON-
C                DANTS A LA MEME FREQUENCE. EN TOTALE IL-Y-A 2*N
C                DGRES DE LIBERTE.
C
C           Routine ayant testee avec le protocol ordinaire le 10/12/198
C  6
C
C
C
C $Id: cerx3s.f,v 1.2 2012/03/30 12:12:42 j_novak Exp $
C $Log: cerx3s.f,v $
C Revision 1.2  2012/03/30 12:12:42  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.1  1997/10/23  08:29:48  eric
c Initial revision
c
C
C $Header: /cvsroot/Lorene/F77/Source/Poisson2d/cerx3s.f,v 1.2 2012/03/30 12:12:42 j_novak Exp $
C
C
	character*120 header
	data header/'$Header: /cvsroot/Lorene/F77/Source/Poisson2d/cerx3s.f,v 1.2 2012/03/30 12:12:42 j_novak Exp $'/


       DIMENSION DEN(NDIMR,NDIMY,*),C64(*),CC(*),CS(*),NDEG(3)
       DATA NDY,NDR,NDZ/0,0,0/
       DATA NEQ/0/

	save	NDR,NDY,NDZ,NEQ,NR,NZ,NY,NNN65,NN64R,NN65R,NY1Z1,MULTR
	save	IRESTR,N6R,MULT1,IREEE2,II,I2,LRR1,LRR2,LDR,N64R,N65R,N6565R
	save	N365R,IRR2R

       NR1=NDEG(1)
       NY1=NDEG(2)
       NZ1=NDEG(3)
C
C           INITIALISATION.
C
       IF(NDY.EQ.NY1.AND.NDZ.EQ.NZ1.AND.NDR.EQ.NR1.AND.NEQ.EQ.NNN64) GO 
     &TO 800
C           PRINT*,'INITIALISATION'
       NDR=NR1
       NDY=NY1
       NDZ=NZ1
       NEQ=NNN64
       NR=NR1-1
       NZ=NZ1-1
       NY=NY1-1
       NNN65=NNN64
       IF((NNN64/8)*8.EQ.NNN64)NNN65=NNN64+1
C
C               PREPARATION DES QUANTITES NCESSAIRES POUR LA TRANSFOR-
C           MATION DU 1er INDEX DU TABLEAU.
C
       NN64R=NNN64
       NN65R=NNN65
       NY1Z1=NY1*NZ1
       MULTR=(NY1Z1)/NN64R
       IRESTR=NY1Z1-MULTR*NN64R
C
C           OPTIMISATION DE NN64R: ON CHERCHE UNE VALEUR DE NN64R QUI
C           SOIT< NNN64 MAIS UN MULTIPLE DE NY1. CELA EN VUE DE REDUIRE
C           LE NOMBRE D'OPERATION DANS LE TRANSFERT DES VALEURS DE DEN
C           DANS C64.
C
       IF(MULTR.GT.0.AND.NNN64.GT.NY1) THEN
       N6R=(NNN64/NY1)*NY1
       MULT1=NY1Z1/N6R
       IREEE2=NY1Z1-MULT1*N6R
       II=0
       I2=0
       IF(IRESTR.GT.0)II=1
       IF(IREEE2.GT.0)I2=1
       IF(MULT1+I2.LE.MULTR+II) THEN
       MULTR=MULT1
       NN64R=N6R
       NN65R=NN64R
       IF((NN64R/8)*8.EQ.NN64R)NN65R=NN64R+1
       IRESTR=IREEE2
       ENDIF
       ENDIF
C
       IF(MULTR.EQ.0) THEN
       NN64R=NY1Z1
       NN65R=NN64R
       IF((NN64R/8)*8.EQ.NN64R)NN65R=NN64R+1
       MULTR=1
       IRESTR=0
       ENDIF
       LRR1=1
       LRR2=NN64R/NY1
       LDR=LRR2-LRR1
       IF(IRESTR.GT.0) THEN
       N64R=IRESTR
       N65R=N64R
       IF((N64R/8)*8.EQ.N64R)N65R=N64R+1 
       ENDIF
C
       N6565R=NN65R+NN65R
       N365R=N6565R+NN65R
       IRR2R=NN64R-(NN64R/NY1)*NY1
C
 800   CONTINUE
C
C           TRANSFORMATION DU PREMIER INDICE
C
C 1111111111111111111111111111111111111111111111111111111111111111111111
C  111
C 1111111111111111111111111111111111111111111111111111111111111111111111
C  111
C
       IF(IRR2R.EQ.0)THEN
C
C           ON EFFECTUE LA TRANSFORMATION DANS LE CAS NN64R MULTIPLE DE 
C  NY1.
C      
C
       LR1=LRR1
       LR2=LRR2
C
       DO 40 LMU=1,MULTR
       JR=0
       DO 1 LR=1,NR1
       JZ=JR 
       DO 2 LZ=LR1,LR2
       DO 3 LY=1,NY1
       C64(LY+JZ)=DEN(LR,LY,LZ)
   3   CONTINUE
       JZ=JZ+NY1
   2   CONTINUE
       JR=JR+NN65R
   1   CONTINUE
C
        CALL CERAMS(NR,NN64R,ITCH,C64,CS,CC)
C
C           ON STOCKE LES COEFFICIENTS DE LA TRANSFORMATION DANS 
C           DEN(LR,LY,LZ).
C
C
       JR=0
       DO 13 LR=1,NR1
       JZ=JR 
       DO 14 LZ=LR1,LR2
       DO 15 LY=1,NY1
       DEN(LR,LY,LZ)=CC(LY+JZ)
  15   CONTINUE
       JZ=JZ+NY1
  14   CONTINUE
       JR=JR+NN65R
  13   CONTINUE
C
       LR1=LR2+1
       LR2=LR1+LDR
  40   CONTINUE
C
C           LE CALCUL EST CONTINUE' SI NY1*NR1 N'EST PAS UN MULTIPLE
C           DE N64Y.
C
       IF(IRESTR.GT.0) THEN
C
       JR=0
       DO 16 LR=1,NR1
       JZ=JR
       DO 17 LZ=LR1,NZ1
       DO 18 LY=1,NY1
       C64(LY+JZ)=DEN(LR,LY,LZ)
  18   CONTINUE
       JZ=JZ+NY1
  17   CONTINUE
       JR=JR+N65R
  16   CONTINUE
C
        CALL CERAMS(NR,N64R,ITCH,C64,CS,CC)
C
C           ON REINTRODUIT LES COEFF. DE FOURIER DANS DEN(LR,LY,LZ)
C
       JR=0
       DO 28 LR=1,NR1
       JZ=JR
       DO 29 LZ=LR1,NZ1
       DO 30 LY=1,NY1
       DEN(LR,LY,LZ)=CC(LY+JZ)
  30   CONTINUE
       JZ=JZ+NY1
  29   CONTINUE
       JR=JR+N65R
  28   CONTINUE
       ENDIF
C
       RETURN
       ENDIF
C
C           CALCUL DE LA TF DANS LE CAS OU N64Y N'EST PAS UN MULTIPLE DE
C   
C           NY1.
C
       IF(IRR2R.GT.0) THEN
C           
       NR64R=NN64R
       N63Y=NY1-NR64R+1
       NR65R=NN65R
       MU2=NR64R
       LZ1=1
       LY1=1
       LRR=1
C
       DO 99 LM=1,MULTR
       IF(LRR.GT.NY1) LRR=1
       LRR1=0
       LZ1=(MU2-NR64R)/NY1+1
       LZ2=.99999+FLOAT(MU2)/NY1
       LY1=+MU2+N63Y-LZ1*NY1
       LY2=MU2-(LZ2-1)*NY1
       JR=0
C
       DO 46 LR=1,NR1
       JJZ=JR
       IF(LZ2.GT.LZ1) THEN
       JZ=JR-LRR+1
C
       DO 41 LY=LRR,NY1
       C64(LY+JZ)=DEN(LR,LY,LZ1)
  41   CONTINUE
       JJZ=NY1+JZ
       LRR1=1
       ENDIF
C
       IF(LZ2.GT.LZ1+1) THEN
       JZ=JJZ
       DO 42 LZ=LZ1+1,LZ2-1
       DO 43 LY=1,NY1
       C64(LY+JZ)=DEN(LR,LY,LZ)
  43   CONTINUE
       JZ=JZ+NY1
  42   CONTINUE
       JJZ=JZ
       ENDIF
C
       IF(LY2.GE.LY1.OR.LRR1.EQ.1) THEN
C
       IF(LRR1.EQ.0) THEN
       JZ=JJZ+1-LY1
       DO 44 LY=LY1,LY2
       C64(LY+JZ)=DEN(LR,LY,LZ2)
  44   CONTINUE
       ENDIF
       IF(LRR1.EQ.1) THEN
       JZ=JJZ
       DO 45 LY=1,LY2
       C64(LY+JZ)=DEN(LR,LY,LZ2)
  45   CONTINUE
       ENDIF
C
       LSY=LY2+1
       ENDIF
       JR=JR+NR65R
  46   CONTINUE   
C
C
        CALL CERAMS(NR,NR64R,ITCH,C64,CS,CC)
C
       JR=0
       DO 71 LR=1,NR1
       JJZ=JR
       IF(LZ2.GT.LZ1) THEN
C
       JZ=JR-LRR+1
       DO 66 LY=LRR,NY1
       DEN(LR,LY,LZ1)=CC(LY+JZ)
  66   CONTINUE
       JJZ=NY1+JZ
       LRR1=1
       ENDIF
C
       IF(LZ2.GT.LZ1+1) THEN
       JZ=JJZ
       DO 67 LZ=LZ1+1,LZ2-1
       DO 68 LY=1,NY1
       DEN(LR,LY,LZ)=CC(LY+JZ)
  68   CONTINUE
       JZ=JZ+NY1
  67   CONTINUE
       JJZ=JZ
       ENDIF
C
       IF(LY2.GE.LY1.OR.LRR1.EQ.1) THEN
       IF(LRR1.EQ.0) THEN
       JZ=JJZ+1-LY1
       DO 69 LY=LY1,LY2
       DEN(LR,LY,LZ2)=CC(LY+JZ)
  69   CONTINUE
       ENDIF
C
       IF(LRR1.EQ.1) THEN
       JZ=JJZ
       DO 70 LY=1,LY2
       DEN(LR,LY,LZ2)=CC(LY+JZ)
  70   CONTINUE
       ENDIF
C
       LSY=LY2+1
       ENDIF
       JR=JR+NR65R
  71   CONTINUE   
C
       LRR=LSY 
       MU2=MU2+NR64R
  99   CONTINUE
       ENDIF
C
       IF(IRESTR.EQ.0) RETURN
C
       JR=0
       LY1=LRR
       LZ1=LZ2
       LZ3=LZ2
       IF(LRR.GT.NY1) THEN
       LY1=1
       LZ3=LZ2+1
       ENDIF
       NR64R=IRESTR
       NR65R=NR64R
       IF((NR64R/8)*8.EQ.NR64R)NR65R=NR64R+1
       N265=NR65R+NR65R
       JR=-LY1+1
C
       DO 76 LR=1,NR1
C
       JZ=JR
       IF(LRR.LE.NY1) THEN
       DO 73 LY=LY1,NY1
       C64(LY+JZ)=DEN(LR,LY,LZ2)
  73   CONTINUE
       LZ3=1+LZ2
       JZ=JZ+NY1
       ENDIF
C
       IF(LZ3.LE.NZ1)THEN
       DO 74 LZ=LZ3,NZ1
       DO 75 LY=1,NY1
       C64(LY+JZ)=DEN(LR,LY,LZ)
  75   CONTINUE
       JZ=JZ+NY1
  74   CONTINUE
       ENDIF  
       JR=JR+NR65R
  76   CONTINUE
C      
        CALL CERAMS(NR,NR64R,ITCH,C64,CS,CC)
C           
C           ON STOCKES LES COEFFICIENTS DE FOURIER DANS DEN(LR,LY,LZ).
C
       JR=-LY1+1
C
       DO 93 LR=1,NR1
C
       JZ=JR
       IF(LRR.LE.NY1) THEN
       DO 90 LY=LY1,NY1
       DEN(LR,LY,LZ2)=CC(LY+JZ)
  90   CONTINUE
       LZ3=1+LZ2
       JZ=JZ+NY1
       ENDIF
C
       IF(LZ3.LE.NZ1)THEN
       DO 91 LZ=LZ3,NZ1
       DO 92 LY=1,NY1
       DEN(LR,LY,LZ)=CC(LY+JZ)
  92   CONTINUE
       JZ=JZ+NY1
  91   CONTINUE
       ENDIF  
       JR=JR+NR65R
  93   CONTINUE
C
       RETURN
C      ENDIF
 1000  FORMAT(1X,10E12.4)
 1010  FORMAT(1X,' ')
       END
