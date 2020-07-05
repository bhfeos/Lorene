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

       SUBROUTINE TFMZS(N,N64,Y,CC)
C
	IMPLICIT NONE

C           ROUTINE POUR LE CALCUL DE LES TF DE N64 FONCTIONS
C           A LA FOIS.
C
C           ARGUMENTS DE LA ROUTINE:
C
C           N   = NOMBRE DES DEGRES DE LIBERTEE. N DOIT
C                 ETRE= A 2**p*3**q*5**r AVEC p,q,r NOMBRS
C                 ENTIERS (p.NE.0).
C
C           N64 = NOMBRES DES FONCTIONS DONT ON VEUT CALCULER
C                 LES TF.
C           Y   = TABLEAU CONTENANT EN IMPUT LES ECHANTILLONS DES FON-
C                 CTIONS. Y DOIT AVOIR (N+4)*(N64+1)  DIMENSIONS.
C                 LES N ECHANTILLONS DES N64 FONCTIONS SONT
C                 STOCKES EN PARALLEL OU EN SERIE. EN PARALEL
C                 ON A:
C                 Y(1),Y(2),...Y(N64) CONTIENT LES VALEURES
C                 DES N64 FONCTIONS DANS LE POINT X=0,
C                 Y(1+N64),Y(2+N64),...Y(N64+N64)
C               LES VALEURES DES FONCTIONS POUR X=2*PI/N
C               Y(1+2*N64),Y(2+2*N64),...Y(N64+2*N64)
C               LES VALEURES POUR X=(2*PI/N)*2 ET AINSI DE 
C                JUSQU' A X=(2*PI/N)*(N-1).
C
C           CC  = OUTPUT: LES COEEFICIENTS DE CC IMPAIRES CONTIENENT
C                LES COSINUS ET LES COEFFICIENTS PAIRES LES SINUS.
C                (CC(2)=0 PAR DEFINITION. L'ARGUMENT DE CC EST DONC
C                COMPRIS ENTRE 1 ET N+1. 
C
C           Routine testee le 14/9/1985, modifiee le 26 juin 1990 
C
C
C $Id: tfmzs.f,v 1.2 2012/03/30 12:12:44 j_novak Exp $
C $Log: tfmzs.f,v $
C Revision 1.2  2012/03/30 12:12:44  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
C LORENE
C
c Revision 1.1  1997/10/23  08:29:01  eric
c Initial revision
c
C
C $Header: /cvsroot/Lorene/F77/Source/Poisson2d/tfmzs.f,v 1.2 2012/03/30 12:12:44 j_novak Exp $
C
C
	character*120 header
	data header/'$Header: /cvsroot/Lorene/F77/Source/Poisson2d/tfmzs.f,v 1.2 2012/03/30 12:12:44 j_novak Exp $'/

	REAL*8
     1  CC,Y,TRIGS
C
	INTEGER N1040,N,N64,NDIM,IFAX,NFON,N65,N63,N6565,NM65
     1 ,N66,L,JM1,JM2,M 
C
       DIMENSION CC(*),Y(*),IFAX(64),TRIGS(1600)
       DATA NDIM/0/
       DATA NFON/0/

	save	NFON,NDIM,N65,N63,N6565,NM65,N66,ifax,trigs

       N1040=1040
       IF(N.GT.N1040) THEN
       PRINT 800,N
  800  FORMAT(10X,'DIMENSION INSUFFISANTES DANS LA ROUTINE'
     , ,' TFM, N=',I5,' > A LA DIMENSION DE TRIGS*2/3=',I5)
             ENDIF
C
       IF(N.EQ.NDIM) GO TO 1
       CALL FAX(IFAX,N,3)
       CALL FFTRIG(TRIGS,N,3)
  1    CONTINUE
C
       IF(NFON.EQ.N64.AND.NDIM.EQ.N) GO TO 4
       NFON=N64
       NDIM=N
       N65=N64
       IF((N64/8)*8.EQ.N64)N65=N64+1
       N63=N64-1
       N6565=N65+N65
       NM65=N*N65
       N66=N65+1
  4    CONTINUE
C
C
C
       DO 10 L=1,NM65
       CC(L)=Y(L)*2
10     CONTINUE
C
C           CALCUL DE LA TF EN PARALLEL
C
       CALL FFT991(CC,Y,TRIGS,IFAX,N65,1,N,N64,-1)
       JM1=N66
       JM2=JM1+N63
       DO 3 L=2,N,2
       DO 2 M=JM1,JM2
       CC(M)=-CC(M)
  2    CONTINUE
       JM1=JM1+N6565
       JM2=JM1+N63
  3    CONTINUE
C
  101  FORMAT(1X,'TFM')
  100  FORMAT(1X,'TFM',10D12.4)
C
       RETURN
       END
C
