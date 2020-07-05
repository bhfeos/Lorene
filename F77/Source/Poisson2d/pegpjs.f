C
C   Copyright (c) 1998 Silvano Bonazzola
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

	SUBROUTINE PEGPJS(NDEG,NDR,NDT,IPAR,JDIM,IND,R1,R2,CC,DEN2
     1 	,BB,DEN1,SOLH,CDEN)

C
C## Routine modifiee le 09/04/96 : teste' pour ipar=4,5
C
C		ROUTINE POUR LA SOLUTION DE POISSON D' UN SCALAIRE 
C		DANS UN ESPACE A 2 OU A 3 DIMENSIONS EN COORDONNES
C		SPHERIQUES 3 DIMENSIONS AVEC PLUSIEURS COQUILLES.
C		DANS LE CAS SYMETRIQUE PAR RAPPORT AU PLAN z=0,
C		OU SUPERSYMETRIQUE CETTE ROUTINE CALCULE UNE SOLUTION
C		PARTICULIERE DE L' EQUATION DE POISSON DANS UN ESPACE
C		A 2 OU 3 DIMENSIONS (VOIR LE DRAPEAU JDIM).
C		LA ROUTINE PEUT AUSSI CALCULER LE POTENTIEL DANS
C		UNE COQUILLE COMPACTIFIEE,(AVEC DEVELOPPEMENT DE la
C		SOLUTION EN PUISSANCE DE u=1/r, OU DANS UNE COQUILLE
C		NORMALE AVEC UNE SOLUTION EN PUISSANCES DE u=1/r.
C
C		LA SOLUTION CALCULEE EST LA SOLUTION DE L'EQUATION
C
C		(D2/dr2+1/r*D/dr-l**2/r**2)PSI = SOURCE     SI JDIM=2
C 
C		OU
C		
C		(D2/dr2+2/r*D/dr-l*(l+1)/r**2)PSI = SOURCE  SI JDIM=3
C
C		L'INPUT DOIT ETRE PAR CONSEQUENT TCHEBYTCTCHEV EN r FOURIER
C		EN theta DANS LE PREMIER CAS, ET TCHEBYTCHEV EN r LEGENDRE
C		EN theta DANS LE 2me CAS.
C
C     N.B.	SI JDIM=1 LA SOLUTION POUR l=1 DANS LA COQUILLE COMPACTIFIEE
C     ----	N'EST PAS CALCULEE (OPERATEUR SINGULIER)
C
C	ARGUMENTS DE LA ROUTINE:
C
C		NDEG	=TABLEAU CONTENANT LE NOMBRE DE DEGREES DE LIBERTE'
C			 POUR LES COORDONNES r, POUR TETA, ET FI
C		NDR	=DIMENSION DU PREMIER INDICE DES TABLEAUX.
C		NDT	=DIMENSION DU 2ME INDICE DES TABLEAUX.
C
C		IPAR	=DRAPEAU: SI
C			
C			IPAR=0 CAS SENS AUCUNE SYMETRIE
C			IPAR=2 CAS SYMETRIQUE PAR RAPPORT AU PLAN z=0 SENS
C			       SUPERSYMETRIE
C		        IPAR=3 CAS ANTISYMERIQUE PAR RAPPORT AU PLAN z=0	
C			       SENSSUPERSIMMETRIES
C			IPAR=4 CAS SYMETRIQUE PAR RAPPORT AU PLAN z=0 AVEC
C			       SUPERSYMETRIES
C			IPAR=5 CAS ANTISYMETRIQUE PAR RAPPORT AU PLAN z=0
C			       AVEC SUPERSYMETREIES
C				LE TERME "SUPERSYMETRIES SIGNIFIE INVARIANCE
C			       DE LA SOURCE PAR RAPPORT A LA TRANSFORMATION
C			       x,y, -> -x,-y
C
C		JDIM	=DRAPEAU: SI JDIM=2 LA ROUTINE CALCULE LA SOLUTION
C			 DE L'EQUATION DE POISSON DANS UN ESPACE A 2 DI-
C		         MENSIONS, SI JDIM=3 L'EQUATION DE POISSON DANS UN
C                       ESPCE A 3 DIMENSIONS.
C		IND	=DRAPEAU: SI IND=1 LA SOLUTION EST CALCULEE EN 
C			 PUISSANCE DE r, SI IND=2 LA SOLUTION EST CALCULEE
C		         EN PUISSANCE DE 1/r DANS UNE GRILLE COMPACTIFIEE,
C			 SI IND=3 LA SOLUTION EST CALCULEE DANS UNE GRILLE 
C			 ORDINAIRE MAIS LA SOLUTION EST EN PUISSANCES DE 1/r
C
C
C		R1	= RAYONS INTERNE  DE LA COQUILLE.
C		R2	=COEFFICIENT DEFINISSANT LE RAYON EXTERNE DE LA
C		         COQUILLE:  RAY=R1+2*R2
C			
C		BB	=TABLEAU DE TRAVAIL DE DIMENSION .GE.(NR1+1)*6
C		
C
C		CC,DEN1,DEN2=TABLEAUX DE TRAVAIL. DIMENSIONS ((NR1+1)*NF)
C
C		SOLH	=TABLEAU (SOLH(NDR,NDT,2) CONTENENT LES 2 SOLUTIONS
C			 HOMOGENES r**l ET 1/r**(l+1)
C		CDEN	=TABLEAU CONTENANT EN IMPUT LES COEFFICIENTS DE LA 
C			 DENSITE' (FOURIER EN FI (TROISIEME INNDICE,LEGENDRE
C		         EN TETA(2me INDICE) ET TCHEBITCHEF EN r (1ere INDICE)
C			 DEN(LR,LY,LF) ET EN OUTPUT LA SOLUTION. ETANT DONNE'
C			 QUE LE POTENTIEL GRAVITATIONNEL A UN COEFFICIENT DE 
C			 PLUS QUE LE TERME SOURCE, LES DIMENSIONS DE DEN DIVENT
C			 ETRE AU MOINS DEN(NR1+4,NT1,NF)
C			 LE STOCAGE DES COEEF. EST LE MEME QUE DANS
C			 FCIR3S POUR IND=7. C'EST A DIRE: POUR LE 3me
C			 INDICE( PARTIE EN PHI) DANS LF=1 IL-Y-A
C			 LE COEFFICIENT CORRESPONDENT A LA FREQUENCE ZERO,
C			 DANS LF =2,3 LES COEFF. CORR. A LA FREQUENCE 1
C			 (COS ET SIN) ET AINSI DE SUITE.
C
C			 DANS LE CAS JDIM=3 LE STOCAGE DE FONCTIONS 
C			ASSOCIEES DE LEGENDRE  EST LE SUIVANT:
C
C			 SI m EST PAIRE LES 2 COEFFICIENTS (COS ET SIN) DE
C							   m	
C        		 LA FONCTION ASSOCIEE DE LEGENDRE P (theta,fi)
C					                   j
C	                 SE TROUVENT DANS DEN(j+1-m,2*m+1),DEN(j+1-m,2*m+1)
C   			 SI m EST IMPAIRE DANS DEN(j-m,2*m),DEN(j-m,2*m+1)
C
C			 LE NOMBRE QUANTIQUE l EST OBTENU
C			 DANS LA FACON SUIVANTE:
C
C				DO LF=1,NF
C				LF2=LF/2
C				LEF=LF2-(MOD(LF2,2)
C				DO LT=1,NT
C				ELLE=LT+LEF-1
C                               .............
C				ENDDO
C				ENDDO
C
C		ATTENTION !   LES DIMENSION MINIMES DES TABLEAUX SONT
C               ----------    NR1+3
C
C		N.B.     POUR IND=2 L'IMPUT ES SUPPOSE TAVOIR ESTE DIVISE
C               ----     PAR u**4, POUR IND=3 PAR u**2. AUCU ALIASING N'EST
C                        EXECUTE A L'INTERIEUR DE LA ROUTINE
C
		IMPLICIT NONE
C
C $Id: pegpjs.f,v 1.2 2012/03/30 12:12:44 j_novak Exp $
C $Log: pegpjs.f,v $
C Revision 1.2  2012/03/30 12:12:44  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.9  1998/09/02  15:26:11  eric
c Correction erreur NR=5 dans l'appel DIRCMS dans la derniere zone pour
c ipar = 3.
c
c Revision 1.8  1998/08/25  10:13:42  eric
c Retour a la version 1.6 fournie par Silvano (traitement cas ipar = 3).
c
c Revision 1.7  1998/07/28  12:54:21  eric
c Retour a la version 1.5
c
c Revision 1.6  1998/07/27  15:32:31  eric
c Nouvelle version fournie par Silvano:
c   traite le cas ipar = 3
c
c N'annule plus solh avant utilisation
c
c Revision 1.5  1997/10/10  08:46:13  eric
c Initialisation a zero de SOLH suivant la valeur de IPAR.
c
C Revision 1.4  1997/10/09 12:59:19  eric
C Correction erreur initialisation SOLH a zero
C
C Revision 1.3  1997/08/07 18:16:02  eric
C Corrige erreur.
C
C Revision 1.2  1997/05/23 11:38:42  hyc
C *** empty log message ***
C
C Revision 1.1  1997/05/07 16:41:08  hyc
C Initial revision
C
C
C $Header: /cvsroot/Lorene/F77/Source/Poisson2d/pegpjs.f,v 1.2 2012/03/30 12:12:44 j_novak Exp $
C
C
	character*120 header
	data header/'$Header: /cvsroot/Lorene/F77/Source/Poisson2d/pegpjs.f,v 1.2 2012/03/30 12:12:44 j_novak Exp $'/

	INTEGER NRL,N257,NDR,NDT,NDEG,NT1,NT,NF2,NR,LF
     1	,JDIM,NR1,LT,LR,JD3,LEL,LELLE,LF1,LF2,ILT,NF,IPAR
     1	,JF,IND,MM,JD1,LY,LF21,IPA, JJ

	double PRECISION COEF,R1,R2,RAY,CDEN,CC1,CS,BB,DEN1
     1	,CC,DEN2,SOLH,RAP,S1,S2,OUT

	PARAMETER (NRL=70,N257=262)
C
	DIMENSION CC(*),CDEN(NDR,NDT,*),COEF(6),SOLH(NDR,NDT,*)
	DIMENSION DEN1(NDR,*),BB(NDR,*),DEN2(NDR,*),ILT(NRL)
	DIMENSION CS(N257),CC1(N257)
	DIMENSION NDEG(3)
C
	NR1=NDEG(1)
	NT1=NDEG(2)
	NF= NDEG(3)
C
	IF(NF.GT.NRL) THEN
	PRINT*,'DIMENSIONS INSUFF. DANS LA ROUTINE PEGPJS'
	PRINT*,'NR1=,NT1=,NF=',NR1,NT1,NF
	STOP
	ENDIF
C
C		LES DIMENSIONS MIMIMES DOIVENT ETRE  = NR1+5 PARCE QUE
C		ON VEUT ELIMINER LE TERME SOURCE DU r**l*LOG(r).
C		SANS CELA ELLES AURAIT PU ETRE NR1+4
C
	IF(NDR.LT.NR1+3) THEN
	PRINT*,'1ERE DIMENSION INSUFFISANTE DANS LA ROUTINE PEGPJS.'
	STOP
	ENDIF
C
	IF(IND.GT.1.AND.NR1.LT.9) THEN
	PRINT*,'ROUTINE PEGPJS: LE NOMBRE DES DEGRES DE LBERTE DOIT ETRE'
	PRINT*,' > 8:,IND,NR1=',IND,NR1
	STOP
	ENDIF
C
	IF(IND.EQ.1.AND.NR1.LT.7) THEN
	PRINT*,'ROUTINE PEGPJS: LE NOMBRE DES DEGRES DE LBERTE DOIT ETRE'
	PRINT*,' > 6 :,IND,NR1=',IND,NR1
	STOP
	ENDIF
C
	IF(NR1.GT.N257) THEN
	PRINT*,'ROUTINE PEGPJS: DIMENSIONS DE CS1 ET CS2 INSUFFISANTES'
	PRINT*,'N257,NR1=',N257,NR1
	STOP
	ENDIF
C
	IF(NDT.LT.NT1) THEN
	PRINT*,'2ME DIMENSION INSUFFISANTE DANS LA ROUTINE PEGPJS'
	STOP
	ENDIF
C
	NR=NR1-1
	NT=NT1-1
	NF2=NF/2
C
	IF((NR1+1)*NF.GT.NDR*NDT) THEN
	PRINT*,'ROUTINE PEGPJS: DIMENSIONS INTERNES INSUFFISANTES'
	STOP
	ENDIF
C
C.............................................................................
C
C			LELLE REPRESENTE LE NOMBRE QUANTIQUE l.
C
	IF(IND.EQ.1) THEN
	RAY=R1+2*R2
	ELSE
	RAY=R1
	ENDIF
C
	S1=R1/RAY
	S2=R2/RAY
C
	COEF(4)=R1
	COEF(5)=R2
C
C		INVERSION DE L'OPERATEUR
C
	COEF(1)=1
	COEF(2)=JDIM-1
	JD1=JDIM-2
	IF(IND.GT.1.AND.JDIM.EQ.3) COEF(2)=0
C
	DO LR=1,NR1
	CC1(LR)=0
	ENDDO
	CC1(1)=2
C
	JJ=2
	IF((IPAR.EQ.2).OR.(IPAR.EQ.3)) JJ=4
	DO LF=1,JJ
	DO LT=1,NT1
	DO LR=1,NR1
	SOLH(LR,LT,LF)=0
	ENDDO
	ENDDO
	ENDDO
C
C		RELATION ENTRE LES VALEURS DE LF ET LE NOMBRE QUANTIQUE l
C
	IF(IPAR.EQ.0) THEN
C
	JF=1
	IF(JDIM.EQ.3.AND.IND.GT.1) THEN
C
C		CALCUL DU POTENTIEL DANS LE CAS 3 DIM. POUR l=0
C		ET DEVELOPPEMENT EN 1/r
C
	RAP=R2**2
	DO LR=1,NR1
	CC(LR)=CDEN(LR,1,1)*RAP
	ENDDO
C
	CALL PRIMS(NR,CC,CS,1,DEN1)
	CALL PRIMS(NR,CS,CC,1,DEN1)
C
	DO LR=1,NR1
	CDEN(LR,1,1)=CC(LR)
	ENDDO
	DO LR=1,5
	CC(LR)=0
	ENDDO
	CC(1)=2*(R1+R2)
	CC(2)=-R2
	CALL DIRCMS(4,N257,1,0,S1,S2,CC,CC1)
	JF=2
	SOLH(1,1,1)=CC(1)
	SOLH(2,1,1)=CC(2)
	SOLH(1,1,2)=2
	ENDIF
C
	DO LEL=JF,NT1
	LELLE=LEL-1
	IF(JD1.EQ.0.AND.IND.EQ.2.AND.LELLE.EQ.1) THEN
	RAP=R2**2
	DO LF=1,MIN0(NF,3)
	DO LR=1,NR1
	CC(LR)=CDEN(LR,2,LF)*RAP
	ENDDO
	CC(NR1)=CC(NR1)*.5
	CC(NR1+1)=0
	CC(NR1+2)=0
	CC(NR1+3)=0
C
	CALL PRIMS(NR+2,CC,CS,1,DEN1)
	CALL EXTM1S(NR+2,1,0,CS,OUT)
	CS(1)=CS(1)-2*OUT
	CALL DIRCMS(NR+2,N257,1,0,R1,R2,CS,CC)
	CALL PRIMS(NR+2,CC,CS,1,DEN1)
	CALL EXTM1S(NR+2,1,0,CS,OUT)
	CS(1)=CS(1)-2*OUT
	CALL DIRCMS(NR+2,N257,1,1,R1,R2,CS,CC)
	CC(NR1)=CC(NR1)*2
C
	DO LR=1,NR1
	CDEN(LR,2,LF)=CC(LR)
	ENDDO
	ENDDO
	GO TO 34
	ENDIF
C
	LF1=LEL*2-1
	IF(LF1.GT.NF) LF1=NF
	DO  LF=1,LF1
	LF2=LF/2
	LT=LEL-LF2
	IF((LF2/2)*2.NE.LF2) LT=LT+1
C
	ILT(LF)=LT
 	ENDDO
C
	COEF(3)=-LELLE*(LELLE+JD1)
C
	CALL GLAGOS(NR,NDR,IND,COEF,BB)
C
	DO LF=1,LF1
	DO LR=1,NR1
	DEN1(LR,LF)=CDEN(LR,ILT(LF),LF)
	ENDDO
	ENDDO
C
	CALL ILGGOS(NDR,LF1,CC,BB,DEN1,DEN2)
C
	DO LF=1,LF1
	DO LR=1,NR1
	CDEN(LR,ILT(LF),LF)=DEN2(LR,LF)
	ENDDO
	ENDDO
 34	CONTINUE
C
	DO LR=1,NR1
	SOLH(LR,LEL,1)=CC1(LR)
	ENDDO
C
	CALL DIRCMS(NR,N257,1,0,S1,S2,CC1,CC)
C
	IF(IND.NE.2) THEN
	LF2=LF1+1
	DO LR=1,NR1
	SOLH(LR,LEL,2)=DEN2(LR,LF2)
	ENDDO
	ENDIF
C
	DO LR=1,NR1
	CC1(LR)=CC(LR)
	ENDDO
C
	ENDDO
C
		RETURN
C
	ENDIF
C
C		CAS AVEC SYMETRIES SEULEMENT PAR RAPPORT AU PLAN z=0
C----------------------------------------------------------------------------
C
	IF(IPAR.EQ.2.OR.IPAR.EQ.3) THEN
C
C		CAS AVEC SYMETRIE PAR RAPPORT AU PLAN z=0
C
C@@@@@
	IPA=IPAR-2
C
C		SOUSCAS l = IPA
C
	LELLE=IPA
	LEL=LELLE*(LELLE+JD1)
!
	IF(JDIM.EQ.3.AND.LEL.EQ.0.AND.IND.GT.1) THEN
C
	RAP=R2**2
	DO LR=1,NR1
	CC(LR)=CDEN(LR,1,1)*RAP
	ENDDO
C
	CALL PRIMS(NR,CC,CS,1,DEN2)
	CALL PRIMS(NR,CS,CC,1,DEN2)
C
	DO LR=1,NR1
	CDEN(LR,1,1)=CC(LR)
	ENDDO
C
	SOLH(1,1,2)=2
	ELSE
C
	COEF(3)=-LEL
C
	CALL GLAGOS(NR,NDR,IND,COEF,BB)
C
	DO LR=1,NR1
	DEN1(LR,1)=CDEN(LR,1,1)
	ENDDO
C
	LF1=1
	
	CALL ILGGOS(NDR,LF1,CC,BB,DEN1,DEN2)
C
	DO LR=1,NR1
	CDEN(LR,1,1)=DEN2(LR,1)
	ENDDO
C
	IF(IND.NE.2) THEN
	DO LR=1,NR1
	SOLH(LR,1,2)=DEN2(LR,2)
	ENDDO
	ELSE
	DO LR=1,NR1
	SOLH(LR,1,2)=0
	ENDDO
C
	ENDIF
	ENDIF
C
C		CALCUL POUR LES CAS l PAIRE SI IPAR=2, ET POUR LE CAS
C		IMPAIRE SI IPAR=3
C
	DO 777 LY=2,NT1
	LEL=LY+LY-2
	LELLE=LEL+IPA
C
	DO LR=1,NR1
	DEN1(LR,1)=CDEN(LR,LY,1)
	ENDDO
C
	COEF(3)=-LELLE*(LELLE+JD1)
	CALL GLAGOS(NR,NDR,IND,COEF,BB)
	CALL ILGGOS(NDR,1,CC,BB,DEN1,DEN2)
C
	DO LR=1,NR1
	CDEN(LR,LY,1)=DEN2(LR,1)
	ENDDO
C
	IF(IND.NE.2) THEN
	DO LR=1,NR1
	SOLH(LR,LY,2)=DEN2(LR,2)
	ENDDO
	ELSE
	DO LR=1,NR1
	SOLH(LR,LY,2)=0
	ENDDO
	ENDIF
C
	IF(NF.LT.4) GO TO 776
	ILT(1)=1
	LF2=1
C
	LF1=MIN0(LEL+LEL+1,NF)
!
	IF(LF1.GE.4) THEN
	DO LF=4,LF1,4
	MM=LF/2
	LT=LY-LF/4
	DO JF=LF,MIN0(LF+1,LF1)
	LF2=LF2+1
	ILT(LF2)=LT
	DO LR=1,NR1
	DEN1(LR,LF2)=CDEN(LR,LT,JF)
	ENDDO
	ENDDO
	ENDDO
C
	COEF(3)=-LELLE*(LELLE+JD1)
	CALL GLAGOS(NR,NDR,IND,COEF,BB)
	CALL ILGGOS(NDR,LF2,CC,BB,DEN1,DEN2)
C
	LF2=1
	DO LF=4,LF1,4
	DO JF=LF,MIN0(LF+1,LF1)
	LF2=LF2+1
	DO LR=1,NR1
	CDEN(LR,ILT(LF2),JF)=DEN2(LR,LF2)
	ENDDO
	ENDDO
	ENDDO
	ENDIF
	LF21=LF2+1
C
776	CONTINUE
C
 777	CONTINUE
C
C		CALCUL D'UNE SOLUTION HOMOGENE
C
	DO LR=1,NR1
	CC1(LR)=0
	ENDDO
C
	DO LR=1,NR1
	CC(LR)=0
	ENDDO
!
	CC1(1)=2
	CC(1)=2
	IF(IPAR.EQ.3) THEN
	CC1(1)=2*(R1+R2)
	CC1(2)=-R2
C
	CC(1)=CC1(1)
	CC(2)=CC1(2)
	ENDIF
C
	IF(IND.GT.1) THEN
C
	CALL DIRCMS(4,N257,1,0,S1,S2,CC,CC1)
C
	DO LR=1,5
	CC(LR)=CC1(LR)
	ENDDO
!
	ENDIF
C
	DO LR=1,5
	SOLH(LR,1,1)=CC1(LR)
	ENDDO
C
	DO LT=2,NT1
	CALL DIRCMS(NR,N257,1,0,S1,S2,CC1,CS)
	CALL DIRCMS(NR,N257,1,0,S1,S2,CS,CC1)
C
	DO LR=1,NR1
	SOLH(LR,LT,1)=CC1(LR)
	ENDDO
	ENDDO
C
	IF(NF.EQ.1) RETURN
C
	IF(IPAR.EQ.2) THEN
	DO LR=1,NR1
	CC1(LR)=0
	ENDDO
C
	CALL DIRCMS(NR,5,1,0,S1,S2,CC,CC1)
!
	DO LR=1,4
	SOLH(LR,1,3)=CC1(LR)
	ENDDO
C
	DO LT=2,NT1
	CALL DIRCMS(NR,N257,1,0,S1,S2,CC1,CS)
	CALL DIRCMS(NR,N257,1,0,S1,S2,CS,CC1)
C
	DO LR=1,NR1
	SOLH(LR,LT,3)=CC1(LR)
	ENDDO
	ENDDO
C
	ELSE
	DO LR=1,NR1
	CC(LR)=0
	CC1(LR)=0
	ENDDO
	CC(1)=2
	CC1(1)=2
	IF(IND.EQ.2) THEN
	CALL DIRCMS(5,N257,1,0,S1,S2,CC1,CC)
	ENDIF
!
	DO LR=1,NR1
	SOLH(LR,1,3)=CC(LR)
	ENDDO
	DO LT=2,NT1
	CALL DIRCMS(NR,N257,1,0,S1,S2,CC,CS)
	CALL DIRCMS(NR,N257,1,0,S1,S2,CS,CC)
	DO LR=1,NR1
	SOLH(LR,LT,3)=CC(LR)
	ENDDO
	ENDDO	
	ENDIF
C
C		CAS l IMPAIRE SI IPAR=2 ET CAS PAIRE SI IPAR=3
C
	JJ=1
	IF(IND.EQ.2) JJ=1+IPA
	DO 778 LT=JJ,NT1
	LEL=LT+LT-1
	LELLE=LEL-IPA
C
	LF1=MIN0(LEL+LEL+1,NF-1)
C
	LF2=0
	DO LF=2,LF1,4
	MM=LF/2
	LY=LT-LF/4
	DO JF=LF,MIN0(LF+1,LF1)
	LF2=LF2+1
	ILT(LF2)=LY
	DO LR=1,NR1
	DEN1(LR,LF2)=CDEN(LR,LY,JF)
	ENDDO
	ENDDO
	ENDDO
C
	COEF(3)=-LELLE*(LELLE+JD1)
	CALL GLAGOS(NR,NDR,IND,COEF,BB)
	CALL ILGGOS(NDR,LF2,CC,BB,DEN1,DEN2)
C
	LF2=0
	DO LF=2,LF1,4
	DO JF=LF,MIN0(LF+1,LF1)
	LF2=LF2+1
	DO LR=1,NR1
	CDEN(LR,ILT(LF2),JF)=DEN2(LR,LF2)
	ENDDO
	ENDDO
	ENDDO
C
	LF21=LF2+1
	IF(IND.NE.2) THEN
	DO LR=1,NR1
	SOLH(LR,LT,4)=DEN2(LR,LF21)
	ENDDO
	ELSE
	DO LR=1,NR1
	SOLH(LR,LT,4)=0
	ENDDO
	ENDIF
 778	CONTINUE
C@@@
C
	RETURN
	ENDIF
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
	IF(IPAR.EQ.4) THEN
	JD3=1
C
	IF(JDIM.EQ.3.AND.IND.GT.1) THEN
C
C		CALCUL DU POTENTIEL DANS LE CAS 3 DIM. POUR l=0
C		ET DEVELOPPEMENT EN 1/r
C
	JD3=2
C
	RAP=R2**2
	DO LR=1,NR1
	CC(LR)=CDEN(LR,1,1)*RAP
	ENDDO
C
	CALL PRIMS(NR,CC,CS,1,DEN1)
	CALL PRIMS(NR,CS,CC,1,DEN1)
C
	DO LR=1,NR1
	CDEN(LR,1,1)=CC(LR)
	ENDDO
C
	CALL DIRCMS(4,N257,1,0,S1,S2,CC1,CC)
	DO LR=1,4
	CC1(LR)=CC(LR)
	ENDDO
	ENDIF
C
	DO 3 LY=JD3,NT1
	LEL=LY+LY-1
	LELLE=LEL-1
C
	LF1=LEL
	IF(LF1.GT.NF) LF1=NF
C
C		RELATION ENTRE LES VALEURS DE LF ET LE NOMBRE QUANTIQUE l
C
	DO 2 LF=1,LF1
	MM=(LF/2)*2
	LT=(LEL-MM)/2+1
C
	ILT(LF)=LT
  2	CONTINUE
C
	COEF(3)=-LELLE*(LELLE+JD1)
C
	CALL GLAGOS(NR,NDR,IND,COEF,BB)
C
	DO 4 LF=1,LF1
	DO 5 LR=1,NR1
	DEN1(LR,LF)=CDEN(LR,ILT(LF),LF)
   5	CONTINUE
   4	CONTINUE
C
	CALL ILGGOS(NDR,LF1,CC,BB,DEN1,DEN2)
C
	DO LF=1,LF1
	DO LR=1,NR1
	CDEN(LR,ILT(LF),LF)=DEN2(LR,LF)
	ENDDO
	ENDDO
C
	LF2=LF1+1
	DO LR=1,NR1
	SOLH(LR,LY,1)=CC1(LR)
	ENDDO
C
	IF(IND.NE.2) THEN
	DO LR=1,NR1
	SOLH(LR,LY,2)=DEN2(LR,LF2)
	ENDDO
	ENDIF
C
	CALL DIRCMS(NR,N257,1,0,S1,S2,CC1,CC)
	CALL DIRCMS(NR,N257,1,0,S1,S2,CC,CC1)
C
  3	CONTINUE
C
	IF(IND.EQ.1.OR.JD1.EQ.0) RETURN
C
C		LES SOLUTIONS HOMOGENES SONT MULTIPLIEE PAR r**2
C
	CALL DIRCMS(NR,NDR,NT1,0,S1,S2,SOLH,DEN1)
	CALL DIRCMS(NR,NDR,NT1,0,S1,S2,DEN1,SOLH)
	SOLH(1,1,1)=2*(R1+R2)
	SOLH(2,1,1)=-R2
	SOLH(1,1,2)=2
	RETURN
C
	ENDIF
C
Cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
C
C		CAS IPAR=5, AVEC SUPERSYMETRIE MAIS ANTISYMETRIQUE PAR
C                           RAPPORT AU PLAN z=0
C
	IF(IPAR.EQ.5) THEN
C
	IF(IND.EQ.1) THEN
	CC1(1)=2*(R1+R2)
	CC1(2)=-R2
	ENDIF
C
	DO 33 LY=1,NT1
	LEL=LY+LY
	LELLE=LEL-1
C
	IF(JD1.EQ.0.AND.IND.EQ.2.AND.LELLE.EQ.1) THEN
C
	RAP=R2**2
	DO LF=1,MIN0(NF,3)
	DO LR=1,NR1
	CC(LR)=CDEN(LR,1,LF)*RAP
	ENDDO
	CC(NR1)=CC(NR1)*.5
	CC(NR1+1)=0
	CC(NR1+2)=0
	CC(NR1+3)=0
C
	CALL PRIMS(NR+2,CC,CS,1,DEN1)
	CALL EXTM1S(NR+2,1,0,CS,OUT)
	CS(1)=CS(1)-2*OUT
	CALL DIRCMS(NR+2,N257,1,0,R1,R2,CS,CC)
	CALL PRIMS(NR+2,CC,CS,1,DEN1)
	CALL EXTM1S(NR+2,1,0,CS,OUT)
	CS(1)=CS(1)-2*OUT
	CALL DIRCMS(NR+2,N257,1,1,R1,R2,CS,CC)
	CC(NR1)=CC(NR1)*2
C
	DO LR=1,NR1
	CDEN(LR,1,LF)=CC(LR)
	ENDDO
	ENDDO
	GO TO 35
	ENDIF
C
	LF1=LEL-1
	IF(LF1.GT.NF) LF1=NF
C
C		RELATION ENTRE LES VALEURS DE LF ET LE NOMBRE QUANTIQUE l
C

	DO 22 LF=1,LF1
	MM=(LF/2)*2
	LT=(LEL-MM)/2
C
	ILT(LF)=LT
  22	CONTINUE
C
	COEF(3)=-LELLE*(LELLE+JD1)
	CALL GLAGOS(NR,NDR,IND,COEF,BB)
C
	DO 44 LF=1,LF1
	DO 55 LR=1,NR1
	DEN1(LR,LF)=CDEN(LR,ILT(LF),LF)
   55	CONTINUE
   44	CONTINUE
C
	CALL ILGGOS(NDR,LF1,CC,BB,DEN1,DEN2)
C
	DO LF=1,LF1
	DO LR=1,NR1
	CDEN(LR,ILT(LF),LF)=DEN2(LR,LF)
	ENDDO
	ENDDO
C
	LF2=LF1+1
C
	IF(IND.NE.2) THEN
	DO LR=1,NR1
	SOLH(LR,LY,2)=DEN2(LR,LF2)
	ENDDO
	ENDIF
  35	CONTINUE
C
	DO LR=1,NR1
	SOLH(LR,LY,1)=CC1(LR)
	ENDDO
C
	CALL DIRCMS(NR,N257,1,0,S1,S2,CC1,CC)
	CALL DIRCMS(NR,N257,1,0,S1,S2,CC,CC1)
C
  33	CONTINUE
C
	IF(IND.EQ.1) RETURN
C

C		LES SOLUTIONS HOMOGENES SONT MULTIPLIEE PAR r**2
C
	CALL DIRCMS(NR,NDR,NT1,0,S1,S2,SOLH,DEN1)
C
	IF(JD1.EQ.0) THEN
	DO LY=1,NT1
	DO LR=1,NR1
	SOLH(LR,LY,1)=DEN1(LR,LY)
	ENDDO
	ENDDO
	RETURN
	ENDIF
C
	CALL DIRCMS(NR,NDR,NT1,0,S1,S2,DEN1,SOLH)
	RETURN
C
	ENDIF
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
C
C  100	FORMAT(1X,10E10.3)
C  101	FORMAT(1X,' ')
C  200	FORMAT(1X,10I4)
C  201	FORMAT(10X,'LT=',I3,' N64=',I3)
C  202	FORMAT(10X,'LT=',I3,' N64=',I3)
C  203	FORMAT(10X,'LT=',I3,' N64=',I3,' N65=',I3,' LF1=',I3,' LF0=',I3)
	RETURN
	END
!
