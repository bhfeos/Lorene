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
       SUBROUTINE CERAMS(N,N64,INDEX,F,CS,CC)

		implicit double precision (a-h,o-z)

C           
C           ROUTINE POUR LE CALCUL DES COEFFICIENTS DE
C           TCHEBYTCHEV DE N64 FONCTIONS  AVEC ECHANTILLONAGE
C           RAREFIE' A L'ORIGINE. L'ECHANTILLONAGE EST EFFECTUE'
C           DANS LES POINTS XL=-DCOS((L-1)*PI/(N*N)+PI/2) DE L'IN-
C           TERVALLE 0<X<1.
C      
C               ARGUMENTS DE LA ROUTINE:
C
C           N   = NOMBRE DES DEGRES DE LIBERTE'
C           N64 = NOMBRE DES FONCTIONS A TRANSFORMER.
C           INDEX  = PARAMETRE: SI LA FONCTION EST PAIRE INDEX
C                 DOIT ETERE=0, SI LA FONCTION EST IMPAIRE
C                INDEX=1
C           F   =TABLEAU CONTENANT LES N64 FONCTIONS QU'ON VEUT
C                TRANSFORMER. LES TABLEAU DOIT CONTENIR LES ECHANTILLONS
C                DES FONCTION EN 'PARALLEL'. (Voir la routine tfms)
C           CS  = TABLEAU DE TRAVAIL
C           CC  = OUTPUT.(COEEF. DE TCHEBITCHEV DE N64 FONCTIONS F
C                     EN 'PARALLEL'.
C               DIMENSIONS MINIMALES DES TABLEAUX = (N64+1)*(N+4)
C......................................................................
C
C      Routine testee le 16/10/1986.
C
C
C $Id: cerams.f,v 1.2 2012/03/30 12:12:42 j_novak Exp $
C $Log: cerams.f,v $
C Revision 1.2  2012/03/30 12:12:42  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.1  1997/10/23  08:38:23  eric
c Initial revision
c
C
C $Header: /cvsroot/Lorene/F77/Source/Poisson2d/cerams.f,v 1.2 2012/03/30 12:12:42 j_novak Exp $
C
C
	character*120 header
	data header/'$Header: /cvsroot/Lorene/F77/Source/Poisson2d/cerams.f,v 1.2 2012/03/30 12:12:42 j_novak Exp $'/


       DIMENSION CC(*),CS(*),F(*),SEN(520),COSE(520),FN11(257)
C
       DATA NCONTR,NFON /0,0/

	save	N520,NC257,X0,X1,PI,PIT,N1,N2,AN2,N21
	save	SEN,COSE,N65,NFON,NCONTR,NM65
	save	NM655,NM265,NM259,NM261,N63,N66,N654
	save	N6565,N356,N357,N365,N6555

       IF(N.EQ.NCONTR) GO TO 999
       N520=520
C
C           PREPARATION DES FONCTIONS AUXILIAIRES.
C
       IF(N.LT.1025) GO TO 8
       PRINT 200,N
 200   FORMAT(10X,'DIMENSION DU TABLEAU SEN INSUFISANTES DANS LA'
     , ,' ROUTINE CHEIM: N=',I4)
       CALL EXIT
  8    CONTINUE
       NC257=257
       IF(N64.GT.257) THEN
       PRINT*,'DIMENSIONS INSUFFISANTES DANS LA ROUT,CERAMS,N64=',N64
       CALL EXIT
       ENDIF
C
       X0=0
       X1=1
       PI=2*ACOS(X0)/N
       PIT=PI/2
       N1=N+1
       N2=N/2
       AN2=X1/N2
       N21=N2+1
       DO 7 L=1,N1
       SEN(L)=.5+SIN((N2-L)*PI)
       COSE(L)=SIN((L-1)*PIT)
  7    CONTINUE
 999   CONTINUE
C
       IF(N64.EQ.NFON.AND.N.EQ.NCONTR) GO TO 998
       N65=N64
       NFON=N64
       NCONTR=N
       IF((N64/8)*8.EQ.N64) N65=N64+1
       NM65=N*N65
       NM655=NM65+N65
       NM265=N2*N65
       NM259=NM265-N65
       NM261=NM265+N65
       N63=N64-1
       N66=N65+1
       N654=N65+N64
       N6565=N65+N65
       N356=N66+N6565
       N357=N356+N63
       N365=N6565+N65
       N6555=5*N65
  998  CONTINUE
C
       IF(INDEX.EQ.0) GO TO 10
C
C           LA FONCTION IMPAIRE EST RENDUE PAIRE PAR MULTIPLICA-
C           TION AVEC SIN(TETA).
C
       JM1=1
       JM2=JM1+N63
       DO 1 L=1,N1
       COSEN=COSE(L)
       DO 2 M=JM1,JM2
       F(M)=F(M)*COSEN
  2    CONTINUE
C
       JM1=JM1+N65
       JM2=JM1+N63
  1    CONTINUE
C
  10   CONTINUE
C
C
C           TRASFORMATION DE TCHEBYTCHEV DE LA FONCTION PAIRE.
C
       DO 3 M=1,N64
       FN11(M)=(F(M)-F(M+NM65))*.5
  3    CONTINUE   
C
C ****   PONDERIZATION DE LA PARIE ASSYMETRIQUE  PAR LA FONCTION SINUS.
C
       DO 30 L=1,NM655
       CC(L)=0
  30   CONTINUE
C
       JM1=NM261
       JM2=NM259
       JM3=JM2+1
       JM4=JM2+N64
       JM5=JM1+1
       JM6=JM5+N64
       DO 4 L =1,N2
       SENN=SEN(L)
       DO 5 M=1,N64
       CC(M)=(F(M+JM2)-F(M+JM1))*SENN
  5    CONTINUE   
C
       DO 6 M=JM3,JM4
       F(M)=F(M)-CC(M-JM2)
  6    CONTINUE   
       DO 11 M=JM5,JM6
       F(M)=F(M)+CC(M-JM1)
  11   CONTINUE   
       JM1=JM1+N65
       JM2=JM2-N65
       JM3=JM2+1
       JM4=JM2+N64
       JM5=JM1+1
       JM6=JM1+N64
  4    CONTINUE   
C
C
C ***   CALCUL DE COEFF. DE FOURIER DE LA FONTION PONDEREE.
C
C
C ***  REARENGEMENT DE COEFF DE FOURIER POUR OBTENIR LES COEFF. DE TCHEB
C
       CALL TFMXS(N,N64,F,CC)
C
       DO 12 M=N356,N357
       F(M)=CC(M)
  12   CONTINUE
C
       JL2=N6555+1
       JL3=JL2+N63
       DO 13 L=6,N1,2
       DO 14 M=JL2,JL3
       CS(M)=CC(M)+F(M-N6565)
  14   CONTINUE
       DO 15 M=JL2,JL3
       F(M)=CS(M)
  15   CONTINUE
       JL2=JL2+N6565
       JL3=JL2+N63
  13   CONTINUE
C
       DO 16 L=1,N64
       CS(L)=F(L+N365)
  16   CONTINUE
C
       JL=N6555
       DO 17 L=6,N1,2
       DO 18 M=1,N64
       CS(M)=CS(M)+F(M+JL)
  18   CONTINUE
       JL=JL+N6565
  17   CONTINUE
C
       DO 19 M=1,N64
       F(M+N65)=(FN11(M)-CS(M))*AN2
  19   CONTINUE
C
       DO 20 M=N66,N654
       CS(M)=F(M)
  20   CONTINUE
C
       JL1=N65*3+1
       JL2=JL1+N63
       JL3=N66-JL1
       DO 21 L=4,N,2
       DO 22 M=JL1,JL2
       CC(M)=-(F(M)+CS(M+JL3))
  22   CONTINUE
       JL1=JL1+N6565
       JL2=JL1+N63
       JL3=N66-JL1
  21   CONTINUE
C
       DO 23 M=N66,N64+N65
       CC(M)=-F(M)
  23   CONTINUE
C 
       IF(INDEX.EQ.0) RETURN
C
       DO 24 L=2,NM655
       CS(L)=CC(L)
  24   CONTINUE
C
C
C           A PARTIR DE LA TRANSFORMME DE LA FONCTION PAIRE ON
C           OBTIENT LES COEFF. DE LA FONCTION IMPAIRE.
C
       DO 25 M=1,N64
       CC(M)=-CC(M)
  25   CONTINUE
C
       JM1=N66
       JM2=JM1+N63
       DO 26 L=2,N
       DO 27 M=JM1,JM2
       F(M)=-CC(M-N65)-2*CS(M)
  27   CONTINUE
C
       DO 28 M=JM1,JM2
       CC(M)=F(M)
  28   CONTINUE
       JM1=JM1+N65
       JM2=JM1+N63
  26   CONTINUE
C
       DO 29 M=NM65+1,NM655
       CC(M)=0
  29   CONTINUE
C 
100    FORMAT(10X,10E12.4)
101    FORMAT(1X,'CERAMS')
       RETURN
       END
