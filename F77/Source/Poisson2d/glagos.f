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

	SUBROUTINE GLAGOS(N,NDIM,IND,COEF,BB)

C
	IMPLICIT NONE
C
C		ROUTINE POUR LA SOLUTION DE L'EQUATION DE LAPLACE  
C	
C		COEF(1)*D2/Dr2 + COEF(2)/(COEF(4)+COEF(5)*x)*D/Dr +
C		COEF(3)/(COEF(4)+COEF(5)*x)**2
C
C		DANS UNE COQUILLE. (r=COEF(4)+COEF(5)*x
C		LA SOLUTION EST AUSSI TROUVE SI LE DEVELOPPEMENT EN 1/r
C		EST DEMANDE. (VOIR PARAMETRE IND.)
C
C		
C		LA SOLUTION EST TROUVEE EN INVERSANT L'OPERATEUR
C
C		(COEF(4)+COEF(5)*x)*(COEF(1)*D2/dr**2+COEF(1)*D/Dr)=
C		(COEF(4)+COEF(5)*x)*SOURCE 
C
C		SI COEF(3)=0 ET L'OPERATEUR
C	
C		(COEF(4)+COEF(5)*x)**2*COEF(1)*D2/d/r**2+
C		(COEF(4)+COEF(5)*x)*COEF(1)*D/Dr+COEF(3)=
C		(COEF(4)+COEF(5)*x)**2)*SOURCE 
C
C		SI COEF(3).NE.0.
C
C		CETTE ROUTINE TRANSFORME LES COEFFICIENTS DE L'EQUATION
C		POUR QU'ILS PUISSENT SERVIR COMME IMPUT DE LA ROUTINE MAGL3S.
C
C		ARGUMENTS DE LA ROUTINE:
C
C		N	= NOMBRE DE DEGRES DE LIBERTE-1
C		NDIM    = DIMENSIONS DES TABLEAUX COMME DEFINIT DANS LE PRO-
C			  GRAMMA APPELLANT.
C
C		IND.	= DRAPEAU: SI IND=1 LE CALCUL DU POTENTIEL EST
C		          EFFECTUE DANS LA COQUILLE EN DEVELOPPEMENT EN r
C		          SI IND=2, LE CACUL EST EFFECTUE DANS UNE COQUILLE
C		          COMACTIFIEE, AVEC DEVELOPPEMNT EN 1/r
C			  SI IND=3 LE DEVELOPPEMENT EST EFFECTUE DANS UNE CO-
C		          QUILLE EN DEVELOPPENT EN 1/r ET LES CONDITIONS
C		          AU CONTOUR SONT CELLE DU VIDE.
C
C		COEF	= COEFFICIENTS DE L'OPERATEUR
C
C		BB	=MATRICE DE L'OPERATEUR DIFFERENTIEL: DIMENSION
C			 .GE. 7
C
C		ATTENTION ! LES DIMENSION MINIMES DES TABLEAUX SONT N1+2
C               ---------
C			    DE MEME DANS LES ROUTINE ILGCOS ET MAGL3S QUI
C			    SONT NECESSAIREMENT ASSOCIEES A CETTE ROUTINE.
C
C		Subroutine modifiee le 20/11/94 en ayant generalise
C		le cas pour IND+@ (COEF(2).NE.0))
C
C
C $Id: glagos.f,v 1.2 2012/03/30 12:12:43 j_novak Exp $
C $Log: glagos.f,v $
C Revision 1.2  2012/03/30 12:12:43  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.2  1997/05/23  11:42:28  hyc
c *** empty log message ***
c
C Revision 1.1  1997/05/08 07:34:16  hyc
C Initial revision
C
C
C $Header: /cvsroot/Lorene/F77/Source/Poisson2d/glagos.f,v 1.2 2012/03/30 12:12:43 j_novak Exp $
C
C
	character*120 header
	data header/'$Header: /cvsroot/Lorene/F77/Source/Poisson2d/glagos.f,v 1.2 2012/03/30 12:12:43 j_novak Exp $'/

	INTEGER N,N1,NN,NN1,NDIM,IND

	double PRECISION COEF,BB,CC,R1,R2,X0

	DIMENSION COEF(*),BB(NDIM,*),CC(10)
	N1=N+1
C
	R1=COEF(4)
	R2=COEF(5)
C
	NN=N
	NN1=NN+1
C
	IF(IND.GT.1) GO TO 667
C
C		CAS R2 < 10*R1 (COQUILLE PAS TROP EPAISSE)
C
	IF(COEF(3).EQ.0.) THEN
C
	CC(1)=0
	CC(2)=1/R2
	CC(3)=R1/R2**2
	CC(4)=0
	CC(5)=COEF(2)/R2
	CC(6)=0
C
	BB(NN1,4)=COEF(4)
	BB(NN1,5)=COEF(5)
	BB(NN1,6)=316.1
C
	CALL MAGL3S(NN,NDIM,CC,BB)
	RETURN
	ENDIF
C
	X0=0
C
	IF(COEF(3).NE.X0) THEN
	CC(1)=COEF(1)
	CC(2)=2*R1*COEF(1)/R2
	CC(3)=(R1/R2)**2*COEF(1)
	CC(4)=COEF(2)
	CC(5)=COEF(2)*R1/R2
	CC(6)=COEF(3)
C
	BB(NN1,4)=COEF(4)
	BB(NN1,5)=COEF(5)
	BB(NN1,6)=317.1
	CALL MAGL3S(NN,NDIM,CC,BB)
	RETURN
	ENDIF
C
  667	CONTINUE
C
C		MATRICE DE ZONES COMPACTIFIEES'
C
C		CALCUL DANS LA ZONE COMPACTIFIEE : DEVELOPPEMENT EN u=1/r
C		ET SANS CONDITION AU CONTOUR: LA SOLUTION EST 
C		DEVELOPPEE EN SERIE DE u**2*Tn(u)
C
	IF(IND.EQ.2) THEN
	NN=N-2
	NN1=NN+1
	CC(1)=COEF(1)
	CC(2)=2*R1*COEF(1)/R2
	CC(3)=(R1/R2)**2*COEF(1)
C
	CC(4)=4*COEF(1)+COEF(2)                 !!!! MODIF LE 20/11/94
	CC(5)=4*COEF(1)*R1/R2+COEF(2)*R1/R2     !!!! MODIF LE 20/11/94
	CC(6)=COEF(3)+2*COEF(1)+2*COEF(2)       !!!! MODIF LE 20/11/94
C
	BB(NN1,4)=COEF(4)
	BB(NN1,5)=COEF(5)
	BB(NN1,6)=315.1
	CALL MAGL3S(NN,NDIM,CC,BB)
	RETURN
	ENDIF
C
C		CALCUL DU POTENTIELL EN EXPANSION EN u=1/r.
C		LA ROUTINE PREPARE LA MATRICE POUR RESOUDRE L'EQUATION
C		(u**2*D2/du2 - l*(l+1)) y = Source
C
	IF(IND.EQ.3) THEN
C
	CC(1)=COEF(1)
	CC(2)=2*R1*COEF(1)/R2
	CC(3)=(R1/R2)**2*COEF(1)
C
	CC(4)=0
	CC(5)=0
	CC(6)=COEF(3)
C
	NN=N
	NN1=N1
C
	BB(NN1,4)=COEF(4)
	BB(NN1,5)=COEF(5)
	BB(NN1,6)=318.1
	CALL MAGL3S(NN,NDIM,CC,BB)
	RETURN
	ENDIF
C
  100	FORMAT(1X,10E11.4)
  101	FORMAT(1X,' ')
	END 
C
