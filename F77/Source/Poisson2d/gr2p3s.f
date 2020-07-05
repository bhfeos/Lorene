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
C
C $Id: gr2p3s.f,v 1.3 2014/03/26 10:44:19 j_novak Exp $
C $Log: gr2p3s.f,v $
C Revision 1.3  2014/03/26 10:44:19  j_novak
C Minor modifications to be compatible with intel 12 compiler.
C
C Revision 1.2  2012/03/30 12:12:43  j_novak
C Cleaning of fortran files
C
C Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
C LORENE
C
c Revision 1.1  1998/06/22  10:31:25  eric
c Initial revision
c
C
C $Header: /cvsroot/Lorene/F77/Source/Poisson2d/gr2p3s.f,v 1.3 2014/03/26 10:44:19 j_novak Exp $
C
C

	SUBROUTINE GR2P3S(NDL1,NDR,NDT,NDF,INDD,IPARR,CC,C64,BB,DEN1
     1	,DEN2,UGRAV,SOLH,ERRE,SOLHH,POTEN)

C
C## Routine modifiee le 09.04.1996 : Teste pour les cas ipar=4,ipar=5
C
C
C		ATTENTION !!
C			     LES DIMENSIONS DU TABLEAU SOLHH DOIVENT
C			     ETRE GE. A
C			     SOLHH(NDR,NDT,8,NDZ)
C			     LES TABLEAUX DEN1,DEN2 DOIVENT ETRE DIMENIONES
C		             DEN1(NDR,NDT,MAX(8,NF)),DEN2(NDR,NDT,MAX(8,NF))
C			     LES TABLEAUX UGRAV,SOLH (NDR,NDT,2)	
C
C		ROUTINE POUR LA SOLUTION DE POISSON D' UN SCALAIRE 
C		DANS UN ESPACE A 3 DIMENSIONS EN COORDONNES SPHERIQUES
C		3 DIMENSIONS AVEC PLUSIEURS COQUILLES.
C		L'ESPACE PEUT ETRE COMPACTIFIE. (VOIR PARAMETRE IND).
C		DANS LE CAS NON COMACTIFIE  LES CONDITIONS AU CONTOUR
C		DE LA DERNIERE COQUILLE SONT CELLES DU VIDE. LA ROUTINE
C		EST ADADPTEE A TRAVAILLER SELON LES DIFFERENTS SYMMETRIES
C		DU PROBLEME.(VOIR DRAPEAU IPAR). LA SOLUTION DANS CHAQUE
C		COQUILLE EST APPROCHEE
C		PAR UN DEVELOPPEMENT EN POLINOMES DE TCHEBYTCHEV EN r OU
C		OU EN 1/r (VOIR DRAPEAU IND)
C
C		SI INDD=1,3 LE POTENTIELLE EST CALCULE EN RESOLVANT L'EQ.
C
C		(r*D2/dr +2*D/dr) POT  =Source *r (SI l=0)  ET
C
C		(r**2*D2/dr +2*r*D/dr -l*(l+1) ) =SOURCE * r**2  SI l > 0
C	
C		(r = R1+R2*(1-x), x=1-COS(THETA))
C
C		POUR INDD=0,2 ON RESOUT L'EQ.
C
C		(D2/dr +2/r*D/dr -l*(l+1)/r**2 ) =SOURCE
C
C			ARGUMENTS DE LA ROUTINE:
C
C		NDL1	=TABLEAU CONTENENT LES PARAMETRES DES COQUILLES:
C			 DANS NDL1(1) IL-Y-A LE NOMBRE DES COQUILLES, DANS
C			 NDL1(2),NDL1(3),...NDL1(NZON+1) LES DEGRES DE LIBERTE
C			 EN r DES DIFFERENTS COQUILLES,DANS NDL1(NZON+2) LE
C			 NOMBRE DES DEGRES DE LIBERTE EN THETA, ET DANS 
C			 NDL1(NZON+3) LE DEGRE DE LIBERTE EN FI
C
C		NDR	=PREMIERE DIMENSION DES TABLEAUX DE TRAVAIL
C			 EXACTEMENT COMME DECLARE DANS LE PROGRAMME APPELLANT.
C		NDT,NDF	=2EME ET 3EME DIMENSION DES TOUS LES TABLEAUX
C			 EXACTEMENT COMME DECLARE DANS LE PROGRAMME APPELLANT.
C
C		INDD	=TABLEAU CONTENENT LES DRAPEAUX POUR CHAQUE COQUILLE:
C                        A EXEPTION DE LA PREMIERE ET DERNIERE COQUILLE,
C		         SI :
C
C		 	           POUR LA PREMIERE COQUILLE, SI:
C	                    INDD=0 LA PREMIERE ZONE EST SUPPOSEE ETRE   
C			           UN NOYAU SPHERIQUE ET LA SOLUTION EST
C				   CALCULLE AVEC UN ECHANTILLONAGE RAREFIE' A
C				   L'ORIGIONE
C			    INDD=1 LA REMIERE ZONE EST UNE COQUILLE ET LA
C			           SOLUTION EST EN FONCTION DE r, ET UNE
C			           CONDITION AU CONTOUR DOIT ETRE IMPOSOEE
C			           AUR LA PARTIE GAUCHE DE LA COQUILLE. SI
C			    INDD=2 LA PREMIERE ZONE PEUT ETRE UN NOYAU SPHE-
C			           RIQUE, LA SOLUTION EST EN r MAIS AVEC E-
C				   CHANTILLONAGE DENSE A L'ORIGINE, LA CONDI-
C			           TION AU CONTOUR A GAUCHE EST AUTOMATIQUEMENT 
C			           REMPLACEE PAR UNE CONDITION DE REGULARITE.
C				(CE CAS N'EST PAS ENCORE IMPLEMENTE)
C					SI:
C                           INDD=3 LA PREMIERE ZONE EST NECESSERAIMENT UNE
C                                  COQUILLE ET LA SOLUTION EST EN u=1/r.
C		
C			    POUR LES COQUILLES INTERMEDIAIRES IND DOIT ETRE
C		            NECESSERAIMENT =1 OU =3. (DEVELLOPEMNT EN r OU
C			    EN 1/r.
C
C				     POUR LA DERNIERE ZONE, SI:
C                           INDD=1 LA SOLUTION EST EN FONCTION DE r, LA
C			           CONDITION AU CONTOUR (SURFACE A DROITE DE
C                                  COQUILLE EST CELLE U VIDE. SI
C                           INDD=2  LA COQUILLE PEUT ETRE COMPACTIFIEE, LA
C			           SOLUTION EST EN u=1/r, SI
C			    INDD=3  LA COQUILLE N'EST POAS COMPACTIFIEE, LA
C			            SOLUTION EST EN u=1/r ET LES CONDITIONS AU
C				   CONTOUR SONT CELLE DU VIDE (A DROITE)
C			
C		IPARR	= DRAPEAU DEFINNISSANT LES SYMMETRIES DU PROBLEME:
C
C			  IPARR=0 AUCUNE SYMMETRIE AST SUPPOSEE EXISTER
C			  IPARR=2 LA SOLUTION EST SYMMETRIQUE PAR RAPPORT
C			         LE PLAN EQUATORIALE z=0
C		          IPARR=3 LA SOLUTION EST ANTISYMMETRIQUE PAR RAPPORT
C                                LE PLAN EQUATORIALE z=0
C
C                         IPARR=4,5 COMME POUR LES CAS IPARR=2,3, MAIS AVEC EN
C                                  PLUS UNE SYMMETRIE PAR RAPPORT LA TRANSFOR-
C                                  MATION x,y -> -x,-y
C
C		CC,C64,BB,DEN2,DEN1,SOLH,UGRAV=TABLEAUX DE TRAVAIL.
C				LE DIMENSIONS MINIMALES DE CES TABLEAUX DE
C			        TRAVAIL SONT POUR CC(NDR,ND),NDR .GE MAX(NR1+2)
C			        ,ND GT. MAX(NF+1,NZON,NDT)
C				POUR C64 GE.(NRMAX+2)*NDT
C				POUR BB(NDR,ND), ND.GE.6,DEN1,
C			        POUR DEN2(NDR,ND,NDZ), ND.GE(NT1,NZON+NZON),
C			        POUR SOLH (NDR,NDT,ND) ND.GE.MAX(NZON,8)
C
C                           OU NRMAX EST LA VALEUR MAXIMUN DE TOUS LES NR1.
C
C		ERRE	= TABLEAU CONTENENT LA VARIABLE r ECHANTILLONNEE AUX
C			  BONS POINTS DANS CHAQUE COQUILLE.
C		SOLHH	=TABLEAU OUTPUT CONTENENT LES 2 SOLUTIONS HOMOGENES
C			 ONRACCORDEES A TRAVERS LES ZONES, PERMETTANT DE SATI-
C			 LES CONDITIONS AU CONTOUR. POUR IPARR=0,4,5
C			 DANS SOLHH(LR,LT,1,LZON) IL-Y-A LA SOLUTION HOMOGENE
C			 RELATIVE A LA COQUILLE LZONieme 
C			 r**l OU u**(l+1), u=1/r (SELON LA VALEUR DE IND(LZON))
C			 DANS SOLLH(LR,LT,2,LZON) IL-Y-A LA SOLUTION HOMOGENE
C			 EN 1/r**(l+1), SI IND(LZON)=1 OU LA SOL. u**l SI
C		         IND(LZON)=3 l ETANT LE NOMBRE QUANTIQUE BIEN CONNU
C			 l=LT-1 SI IPARR=0, l=LT+LT-2 SI IPARR=4, l=LT+LT-1
C			 SI IPARR=5, LES DIMENSIOND DE SOLHH, DOIVENT TER .GE.
C		         AUX VALEURS SUIVANTE: NDR .GE. A LA VALEUR MAX.
C			 DES DEGRES DE LIBERTE EN r DE CHAQUE COQUILLE,
C			 NDT .GE. AU NOMBRE DES DEGRES EN theta, 
C			 NDF .GE.MAX(NF,8), NDZ, (SOLHH(NDR,NDT,NDF,NDZ).GE.
C			 AU NOMBRE DES COQUILLES.
C
C			SI IPARR=2 ON A DANS SOURHH(LR,LT,1,LZON) ET 
C			SOLHH(LR,LT,2,LZON), REASPECTIVEMENT LES SOLUTIONS
C			HOMOGENES EN r**l et 1/r**(l+1) (SI IND(LZON)=1,3,VOIR
C			PLUS HAUT) POUR LES VALEURSE DE l PAIRES (l=LT+LT-2)
C			 ,DANS SOULH(LR,LT,3,LZON) ET SOURHH(LR,LT,4,LZON)
C			LES SOLUTIONS HOMOGENES EN r**l, et 1/r**(l+1), POUR
C			l IMPAIRE, l=LT+LT-1
C
C			SI IPARR=3 ON A DANS SOURHH(LR,LT,1,LZON) ET 
C			SOURHH(LR,LT,2,LZON), REASPECTIVEMENT LES SOLUTIONS
C			HOMOGENES EN r**l et 1/r**(l+1) (SI IND(LZON)=1,3,VOIR
C			PLUS HAUT) POUR LES VALEURSE DE l IMPPAIRES (l=LT+LT-1)
C			 ,DANS SOURHH(LR,LT,3,LZON) ET SOURHH(LR,LT,4,LZON)
C			LES SOLUTIONS HOMOGENES EN r**l, et 1/r**(l+1), POUR
C			l PAIRE, l=LT+LT
C			
C		POTEN	=TABLEAU IMPUT CONTENENT LES COEFFICIENTS DE FOURIERR
C 			 EN PHI,DE TCHEBYTCHEV EN THETA ET EN r DE LA SOURCE DU 
C			 POTENTIEL DANS LE DIFFERENTES COQUILLES. 
C			  EN OUTPUT ON A LE POTENTIEL DANS L'ESPACE DE TCHE
C			 BITCHEV, LEGENDRE ET FOURIER.
C
	IMPLICIT NONE

	character*120 header
	data header/'$Header: /cvsroot/Lorene/F77/Source/Poisson2d/gr2p3s.f,v 1.3 2014/03/26 10:44:19 j_novak Exp $'/

	INTEGER N257,NDR,NDT,NDF,NZON,NT1,NF,NDZ,IPARR
     1	,NDL1,NDEG,NR1,NR,LR,LZON,INT,IND
     1	,LF,LT,LELLE,LELL1,LEL,LF1,LF2,ILT,J,L,J2,NZ2
     1	,LZ1,IPAR,IPA,IDEAL,INDD,MM,JJ,IPA2
     1  ,IQAR,JLZ,INT1,IPP

	double PRECISION ERRE,BB,RR1,RR2
     1	,DEN1,CC,UGRAV,C64,DEN2,POTEN,SOLHH,R1,RA1
     1	,SOLH,R2,VA1,DE1,RAP

		PARAMETER (INT=11,N257=260)
C
	DIMENSION SOLHH(NDR,NDT,8,*),POTEN(NDR,NDT,NDF,*)
	DIMENSION CC(NDR,*),C64(*),BB(NDR,*),ERRE(NDR,*),RR1(INT)
	DIMENSION SOLH(NDR,NDT,*),UGRAV(NDR,NDT,*),RR2(INT)
	DIMENSION DEN1(NDR,NDT,*),DEN2(NDR,NDT,*)
	DIMENSION NDL1(*),NDEG(3),INDD(*),ILT(N257)
	DIMENSION IQAR(N257),INT1(N257)
	DIMENSION VA1(1), DE1(1)
C
	NZON=NDL1(1)
C
		IF (IPARR.EQ.3) THEN
		PRINT*,'ROUTINE GR2P3S: LE CAS IPARR=3 NE MARCHE PAS'
                PRINT*,'ENCORE. Sorry !'
		STOP
		ENDIF
C
	IF (NZON.GT.INT) THEN
	PRINT*,'ROUTINE GR2P3S: DIMENSIONS DES TABLEAUX ILT,NDL INSUFFISANTES'
	STOP
	ENDIF
C
	NT1=NDL1(NZON+2)
	NF=NDL1(NZON+3)
	NDZ=NDL1(NZON+4)
C
	IF (NT1.GT.NDT) THEN
	PRINT*,'ROUTINE GR2P3S: TABLEAUX INSUFFISAMMENT DIMENSIONNES'
	PRINT*,'NDT,NT1=',NDT,NT1
	STOP
	ENDIF
C
	IF (NF.GT.NDF) THEN
	PRINT*,'ROUTINE GR2P3S: TABLEAUX INSUFFISAMMENT DIMENSIONNES'
	PRINT*,'NDF,NF=',NDF,NF
	STOP
	ENDIF
C
	IF (NF.GT.N257.OR.NT1.GT.N257) THEN
	PRINT*,'ROUTINE GR2P3S: TABLEAUX INTERNES INSUFFISAMMENT'
	PRINT*,' DIMENSIONNES,NF,NT1,N257=',NF,NT1,N257
	STOP
	ENDIF
C
	IPA2=IPARR/2
	IPA=MOD(IPARR,2)
C
	IF (IPARR.EQ.0) THEN
C
	INT1(1)=NT1
	IF (NF.GE.2) THEN
	DO LF=2,NF,4
	J2=LF/2-1
	DO J=LF,MIN(LF+1,NF)
	INT1(J)=NT1-J2
	ENDDO
	ENDDO
	ENDIF
C
	IF (NF.GE.4) THEN
	DO LF=4,NF,4			
	DO J=LF,MIN(LF+1,NF)
	INT1(J)=NT1-LF/2
	ENDDO
	ENDDO
	ENDIF
	ENDIF
C	
	IF (IPA2.EQ.1) THEN
	INT1(1)=NT1-IPA
	IQAR(1)=IPA
	IF (NF.GE.4) THEN
	DO LF=4,NF,4			
	DO J=LF,MIN(LF+1,NF)
	INT1(J)=NT1-(LF/2+IPA)/2
	IQAR(J)=IPA
	ENDDO
	ENDDO
	ENDIF
C
	IF (NF.GE.2) THEN
	LEL=MOD(IPARR-1,2)
C
	DO LF=2,NF,4
	DO J=LF,MIN(LF+1,NF)
	INT1(J)=NT1-(LF/2+LEL)/2+IPA
	IQAR(J)=LEL
	ENDDO
	ENDDO
	ENDIF
C
	ENDIF
C
	IF (IPA2.EQ.2) THEN
	DO LF=1,NF
	INT1(LF)=MIN(NT1-LF/2+IPA,NT1)
	IQAR(LF)=IPA
	ENDDO
	ENDIF
C
	NDEG(2)=NT1
	NDEG(3)=NF
C
	LZ1=1
	RR2(1)=ERRE(NDL1(2),1)
	IF (INDD(1).EQ.0) LZ1=2
	DO LZON=LZ1,NZON
	NR1=NDL1(LZON+1)
	NDEG(1)=NR1
C
	IF (NR1.GT.NDR) THEN
	PRINT*,'ROUTINE GR2P3S: TABLEAUX INSUFFISAMMENT DIMENSIONNES'
	PRINT*,'NDR,NR1=',NDR,NR1
	STOP
	ENDIF
C
	IF (NR1.GT.N257) THEN
	PRINT*,'ROUTINE GR2P3S: TABLEAUX INTERNES INSUFFISAMMENT'
	PRINT*,' DIMENSIONNES,NR1,N257=',NR1,N257
	STOP
	ENDIF

	RR1(LZON)=ERRE(1,LZON)
	RR2(LZON)=(ERRE(NR1,LZON)-RR1(LZON))*.5
	ENDDO
C	
C------------------------------------------------------------------------------
C
C		                      SI IDEAL=1 ON DESIALISE
	IDEAL=1
C
C------------------------------------------------------------------------------
C
C
	IPAR=IPARR
C
	IF (NF.EQ.1.AND.IPAR/2.EQ.1) IPAR=MIN(2*IPARR,5)
	IPA2=IPAR/2
C
	IF (NZON.GT.NDZ) THEN
	PRINT 101
	PRINT*,'LA DIMENSION NDZ EST INSUFFISANTE DANS LA ROUTINE GR2P3S:'
	PRINT*,'NDZ,NZON=',NDZ,NZON
	STOP
	ENDIF
C
	IF (NT1.GT.NDT) THEN
	PRINT*,'ROUTINE GR2P3S:'
	PRINT*,'LA 2ME DIMENSION DES TABLEAUX DEN2,DEN2,UGRAV,SOLH,SOLHH'
	PRINT*,'POTEN EST INSUFFISANTE: NT1,NDT=',NT1,NDT
	STOP
	ENDIF
C
	IF (INT.LT.NF.OR.INT.LT.NZON) THEN
	PRINT*,'ROUTINE GR2P3S: NF OU NZON > INT, NF,NZON,INT=',INT,NF,NZON
	STOP
	ENDIF
C
	IF (IPARR.EQ.2.OR.IPARR.EQ.3) THEN
	IF (NT1.LT.NF/2) THEN
	PRINT*,'ROUTINE GR2P3S'
	PRINT*,'POUR UN CALCUL CORRECTE DU POTENTIEL'
	PRINT*,'NF+NF DOIT ETRE < A NT1,IPARR=',NF+NF,NT1,IPARR
	STOP
	ENDIF	
	ELSE
	IF (NT1.LT.NF) THEN
	PRINT*,'ROUTINE GR2P3S'
	PRINT*,'POUR UN CALCUL CORRECTE DU POTENTIEL'
	PRINT*,'NF DOIT ETRE < A NT1'
	STOP
	ENDIF	
	ENDIF
C	
	DO LZON=1,NZON
	NR1=NDL1(LZON+1)
	DO J=1,8
	DO LT=1,NT1
	DO LR=1,NR1
	SOLHH(LR,LT,J,LZON)=0
	ENDDO
	ENDDO
	ENDDO
	ENDDO
C
	NR1=NDL1(2)
C	
	IPA=0
	IF ((IPAR/2)*2.NE.IPAR) IPA=1
C
	IF (INDD(1).EQ.2) THEN
	PRINT*,'LE NOYAU INTERNE N ADMET PAS DE DEVELOPEMENT EN 1/r'
	STOP
	ENDIF
C	
	IF (INDD(1).EQ.0) THEN
C
	NR1=NDL1(2)
	NR=NR1-1
	R2=RR2(1)
	DO LF=1,NF
	DO LT=1,NT1
	DO LR=1,NR1
	DEN1(LR,LT,LF)=POTEN(LR,LT,LF,1)
	ENDDO
	ENDDO
	ENDDO
C
	NDEG(1)=NR1
C
C		DESIALISATION DANS LE NOYAU AVEC ECHANTILLONAGE RAREFIE
C
	CALL GRGP2S(NDEG,NDR,NDT,IPAR,CC,DEN2,BB,SOLH,DEN1)
C
	DO LF=1,NF
	DO LT=1,NT1
	DO LR=1,NR1
	POTEN(LR,LT,LF,1)=DEN1(LR,LT,LF)*R2**2
	ENDDO
	ENDDO
	ENDDO
C
	DO LT=1,NT1
	DO LR=1,NR1
	SOLHH(LR,LT,1,1)=SOLH(LR,LT,1)
	ENDDO
	ENDDO
	IF (IPA2.EQ.1) THEN
	DO LT=1,NT1
	DO LR=1,NR1
	SOLHH(LR,LT,3,1)=SOLH(LR,LT,2)	
	ENDDO
	ENDDO
	ENDIF
C
	IF (NZON.EQ.1) RETURN
	ENDIF
C
	LZ1=1
	IF (INDD(1).EQ.0) LZ1=2
C
	DO LZON=LZ1,NZON

	NR1=NDL1(LZON+1)
	NDEG(1)=NR1
	NR=NR1-1
	DO LF=1,NF
	DO LT=1,NT1
	DO LR=1,NR1
	DEN1(LR,LT,LF)=POTEN(LR,LT,LF,LZON)
	ENDDO
	ENDDO
	ENDDO
C
C		SONT SUPRIMEES EN ANNULANT LES 2
C		DERNIER COEFFICIENT DU DEVELLOPEMENT DE TCHEBYTCHEV
C		EN r (OU EN 1/r) ET MODIFIANT LES 2 AVANTDERNIERS COEFF.
C		DE TEL SORTE QUE LES VALEURS SUR LES BORDS DES COQUILLES RESTE
C		INCHANGEES. 
C
	DO LF=1,NF
	DO LT=1,NT1
	DEN1(NR1-2,LT,LF)=DEN1(NR1-2,LT,LF)+.5*DEN1(NR1,LT,LF)
 	DEN1(NR1-3,LT,LF)=DEN1(NR1-3,LT,LF)+   DEN1(NR ,LT,LF)
	DEN1(NR,LT,LF)=0
	DEN1(NR1,LT,LF)=0
	ENDDO
	ENDDO
C
	IND=INDD(LZON)	
	R1=RR1(LZON)
	R2=RR2(LZON)
C
	IND=INDD(LZON)
C
	NDEG(1)=NR1
	R1=RR1(LZON)
	R2=RR2(LZON)
C
	CALL PEGPJS(NDEG,NDR,NDT,IPAR,2,IND,R1,R2,C64,CC
     1	,BB,DEN2,SOLH,DEN1)
C
C		LES SOLUTIONS PARTICULIERES SONT STOCKEES DANS POTEN
C
	IF (IND.EQ.2) THEN
	DO LT=1,NT1
	DO LR=1,NR1
	SOLH(LR,LT,2)=0
	SOLH(LR,LT,4)=0
	ENDDO
	ENDDO
	ENDIF
C
	DO LF=1,NF
	DO LT=1,NT1
	DO LR=1,NR1
	POTEN(LR,LT,LF,LZON)=DEN1(LR,LT,LF)
	ENDDO
	ENDDO
	ENDDO
C 
C		LES SOLUTIONS HOMOGENES SONT STOCKEES DANS SOLHH
C
	J2=2
	IF (IPAR.EQ.2.OR.IPAR.EQ.3) J2=4
C
	DO J=1,J2
	DO LT=1,NT1
	DO LR=1,NR1
	SOLHH(LR,LT,J,LZON)=SOLH(LR,LT,J)
	ENDDO
	ENDDO
	ENDDO
	ENDDO
C
C		CALCUL DES QUANTITES NECESSAIRES AU RACCORDEMENT
C		DANS DEN1(LT,1,LZON) IL-Y-A LES VALEURS DES SOL.
C		HOMOGENES A GAUCHE DE L'INTERVALLE, DANS DEN1(LT,2,LZON)
C		LES VAL. A DROTE DE L'INTERVALLE POUR LA PREMIERE SOLUTION 
C		HOMOGENE, DANS DEN1(LT,3,LZON) ET DANS DEN1(LT,LZON,4) LA MEMME
C		CHOSE POUR LA 2EME SOLUTION HOMOGENE. DANS LE CAS IPAR=2
C		IPARR=3, LA MEME CHOSE POUR LES SOLUTIONS CORRESPONDANTES A 
C		UN l DE LA MEME PARITE, ET DANS DEN1(LT,5,6,7,8,LZON) L'ERQUIVA
C		LENT POU LES SOLUTIONS CORRESPONDANT A l DE LA PARITE OPPOSE
C
C		DANS DEN2 IL-Y-A LES VALEURS DES DERIVES SUR LES BORDES, 
C		STOCKEES COMME POUR DEN1
C
	DO J=1,8
	DO LT=1,NT1
	DEN1(LT,J,1)=0
	DEN2(LT,J,1)=0
	ENDDO
	ENDDO
C
	J2=2
	IF (IPA2.EQ.1) J2=4
C
	JLZ=1
	IF (INDD(1).GT.0) GO TO 334
	NR1=NDL1(2)
	NR=NR1-1
	JLZ=2
C
	RAP=1/ERRE(NR1,1)
C
	IPA=MOD(IPAR,2)
C
	IF (IPA2.EQ.0) THEN

	DO LT=1,NT1
	DEN1(LT,1,1)=0
	DEN2(LT,1,1)=0
	ENDDO
C
	DO LT=1,NT1
	DEN1(LT,2,1)=1
	ENDDO
C
	DO LT=1,NT1
	DEN2(LT,2,1)=(LT-1+IPA)*RAP
	ENDDO
	ENDIF
C	
	IF (IPA2.EQ.1) THEN
	DO LT=1,NT1
	DEN1(LT,1,1)=0
	DEN2(LT,1,1)=0
	ENDDO
	DO LT=1,NT1
	DEN1(LT,2,1)=1
	DEN2(LT,2,1)=(LT+LT-2+IPA)*RAP
	ENDDO
C
	DO J=3,5
	DO LT=1,NT1
	DEN1(LT,J,1)=0
	DEN2(LT,J,1)=0
	ENDDO
	ENDDO
C
	DO LT=1,NT1
	DEN1(LT,6,1)=1
	DEN2(LT,6,1)=(LT+LT-1+IPA)*RAP
	ENDDO	
C
	DO J=7,8
	DO LT=1,NT1
	DEN1(LT,J,1)=0
	DEN2(LT,J,1)=0
	ENDDO
	ENDDO
	ENDIF
C
	IF (IPA2.EQ.2) THEN
C
	DO LT=1,NT1
	DEN1(LT,1,1)=0
	DEN2(LT,1,1)=0
	ENDDO
C
	DO LT=1,NT1
	DEN1(LT,2,1)=1
	ENDDO
C
	DO LT=1,NT1
	DEN2(LT,2,1)=(LT+LT-2+IPA)*RAP
	ENDDO
	ENDIF
C
 334	CONTINUE
C
	DO LZON=JLZ,NZON
	NR1=NDL1(LZON+1)
	NDEG(1)=NR1
	NDEG(2)=NT1
	NR=NR1-1
	NDEG(2)=NT1
	RAP=1/RR2(LZON)
C
	RA1=RAP
	IF (INDD(LZON).GT.1) RA1=-ERRE(1,LZON)**2*RAP
C
	DO J=1,J2
	JJ=J+J
	DO LT=1,NT1
	DO LR=1,NR1
	CC(LR,LT)=SOLHH(LR,LT,J,LZON)
	ENDDO
	ENDDO
C	
	CALL EXTM2S(NDEG,NDR,1,0,CC,C64)
	CALL EXTM2S(NDEG,NDR,1,1,CC,UGRAV)
C
	DO LT=1,NT1
	DEN1(LT,JJ,LZON)=C64(LT)
	ENDDO
C
	DO LT=1,NT1
	DEN2(LT,JJ,LZON)=UGRAV(LT,1,1)*RA1
	ENDDO
C
	CALL EXTM2S(NDEG,NDR,0,0,CC,C64)
	CALL EXTM2S(NDEG,NDR,0,1,CC,UGRAV)
C
	DO LT=1,NT1
	DEN1(LT,JJ-1,LZON)=C64(LT)
	ENDDO
	DO LT=1,NT1
	DEN2(LT,JJ-1,LZON)=UGRAV(LT,1,1)*RA1
	ENDDO
	ENDDO
	ENDDO
C
C		CALCUL DES VALEURS SUR LES BORDES DES COQUILLES DES SOLUTIONS
C		PARTICULIERES, LES VALEURS SONT STOCKEES DANS SOLHH DANS LA-
C		CON SUIVANTE: CAS IPARR=0,4,5: DANS SOLHH(LT,LF,3,LZON) ET
C		SOLHH(LT,LF,4,LZON) IL-Y-A LA VALEUR DE LA SOLUTION ET DE SA
C		DERIVEE SUR LE BORD INTERIEUR DE LA COQUILLE, DANS
C
C		SOLHH(LR,LT,5,LZN) ET SOLHH(LR,LT,6,LZON) LES VALEURS DE
C		LA SOLUTION ET DE SA DERIVE SUR LE BORD EXTERIEUR DE LA 
C		COQUILLE.
C		SI IPAR=2,3 LA MEMEME ME ILFAUT REMPLCER 3,4,5,6,PAR
C		5,6,7,8
C
	J2=0
	IF (IPA2.EQ.1) J2=2
C
	JLZ=1
	NR1=NDL1(2)
	NR=NR1-1
C
	RAP=1/RR2(1)
	IF (INDD(1).GT.0) GO TO 333
	JLZ=2
C
	IF (IPARR.GT.0) THEN
	DO LF=1,NF
	LF1=INT1(LF)
	IPP=IQAR(LF)
	DO LT=1,LF1
	DO LR=1,NR1
	CC(LR,LT)=POTEN(LR,LT,LF,1)
	ENDDO
	ENDDO
C
	CALL EXRM1S(NR,NDR,LF1,1,0,IPP,CC,C64)
	CALL EXRM1S(NR,NDR,LF1,1,1,IPP,CC,UGRAV)
C
	DO LT=1,LF1
	SOLHH(LT,LF,J2+5,1)=C64(LT)
	ENDDO
C
	DO LT=1,LF1
	SOLHH(LT,LF,J2+6,1)=UGRAV(LT,1,1)*RAP
	ENDDO
	ENDDO
C
	ELSE
C
C		CAS IPARR=0
C
	DO LF=1,NF
	LF1=INT1(LF)
	DO LT=1,LF1
	IPP=MOD(LT-1,2)
	DO LR=1,NR1
	CC(LR,1)=POTEN(LR,LT,LF,1)
	ENDDO
C
	CALL EXRM1S(NR,NDR,1,1,0,IPP,CC,VA1)
	CALL EXRM1S(NR,NDR,1,1,1,IPP,CC,DE1)
C
	SOLHH(LT,LF,J2+5,1)=VA1(1)
	SOLHH(LT,LF,J2+6,1)=DE1(1)*RAP
	ENDDO
	ENDDO
	ENDIF
C
  333	CONTINUE
C
C		CALCUL DES VALEUERES EXTREMES DANS LES COQUILLES EXTERIEURES
C		OU DANS LE NOYAU AVEC ECHANTILLONNAGE FIN
C
	IF (NZON.GE.JLZ) THEN
	DO LZON=JLZ,NZON
C
	NR1=NDL1(LZON+1)
	NDEG(1)=NR1
	NR=NR1-1
	NDEG(1)=NR1
	RAP=1/RR2(LZON)
C
	DO LF=1,NF
	DO LT=1,NT1
	DO LR=1,NR1
	CC(LR,LT)=POTEN(LR,LT,LF,LZON)
	ENDDO
	ENDDO
C
	CALL EXTM2S(NDEG,NDR,1,0,CC,C64)
	CALL EXTM2S(NDEG,NDR,1,1,CC,UGRAV)
C
	DO LT=1,NT1
	SOLHH(LT,LF,J2+5,LZON)=C64(LT)
	ENDDO
C
	RA1=RAP
	IF (INDD(LZON).GT.1) RA1=-RAP*ERRE(NR1,LZON)**2
	DO LT=1,NT1
	SOLHH(LT,LF,J2+6,LZON)=UGRAV(LT,1,1)*RA1
	ENDDO
C
	CALL EXTM2S(NDEG,NDR,0,0,CC,C64)
	CALL EXTM2S(NDEG,NDR,0,1,CC,UGRAV)
C
	DO LT=1,NT1
	SOLHH(LT,LF,J2+3,LZON)=C64(LT)
	ENDDO
	RA1=RAP
	IF (INDD(LZON).GT.1) RA1=-RAP*ERRE(1,LZON)**2
	DO LT=1,NT1
	SOLHH(LT,LF,J2+4,LZON)=UGRAV(LT,1,1)*RA1
	ENDDO
	ENDDO
	ENDDO
	ENDIF
C
	IF (NZON.EQ.1) THEN
C
C		CAS AVEC UNE SEULE COQUILLE: ON IMPOSE A DROITE LES CONDITIONS
C		AU CONTOUR DU VIDE (SI INDD(1).NE.0,OU NE. 2)
C
	IF (IPAR.EQ.0) THEN
	DO LF=1,NF
	MM=LF/2
	DO LT=1,NT1-MM
	LEL=LT+MM
	LELLE=LEL-1
	RAP=-(SOLHH(LT,LF,5,1)*LEL+SOLHH(LT,LF,6,1)*ERRE(NR1,1))/
     1	(DEN1(LEL,2,1)*LEL+DEN2(LEL,2,1)*ERRE(NR1,1))
	DO LR=1,NR1
	POTEN(LR,LT,LF,1)=POTEN(LR,LT,LF,1)+RAP*SOLHH(LR,LEL,1,1)
	ENDDO
	ENDDO
	ENDDO
	IF (INDD(1).EQ.0.OR.INDD(1).EQ.2) THEN
	RETURN
	ELSE
C
C		RACCORDEMENT AVEC LE VIDE
C
	DO LT=1,NT1
	RAP=(DEN1(LT,2,1)*LT+DEN2(LT,2,1)*ERRE(NR1,1))/(DEN1(LT,4,1)*
     1	LT+ERRE(NR1,1)*DEN2(LT,4,1))
	DO LR=1,NR1
	SOLHH(LR,LT,1,1)=SOLHH(LR,LT,1,1)-RAP*SOLHH(LR,LT,2,1)
	ENDDO
	ENDDO
C
	ENDIF
	ENDIF
C
	IF (IPAR.EQ.2.OR.IPAR.EQ.3) THEN
C
	J2=IPA+1
	DO LT=1,NT1
	LELLE=LT+LT-2+IPA
	RAP=-(SOLHH(LT,1,7,1)*(LELLE+1)+SOLHH(LT,1,8,1)*ERRE(NR1,1))/
     1	(DEN1(LT,2,1)*(LELLE+1)+DEN2(LT,2,1)*ERRE(NR1,1))
	DO LR=1,NR1
	POTEN(LR,LT,1,1)=POTEN(LR,LT,1,1)+RAP*SOLHH(LR,LT,1,1)	
	ENDDO
	ENDDO
C
	IF (NF.GT.1) THEN
C
	DO LF=2,NF,4
	MM=(LF/2)
	DO J=LF,MIN(LF+1,NF)
	DO LT=1,NT1-MM
	LEL=LT+IPA+MM-1
	LELLE=LT+LT-2+IPA+MM
	RAP=-(SOLHH(LT,J,7,1)*(LELLE+1)+SOLHH(LT,J,8,1)*ERRE(NR1,1))/
     1	(DEN1(LEL,6,1)*(LELLE+1)+DEN2(LEL,6,1)*ERRE(NR1,1))
	DO LR=1,NR1
	POTEN(LR,LT,J,1)=POTEN(LR,LT,J,1)+RAP*SOLHH(LR,LEL,3,1)
	ENDDO
	ENDDO
	ENDDO
	ENDDO
C
	DO LF=4,NF,4
	MM=(LF/2)
	DO J=LF,MIN(LF+1,NF)
	DO LT=1,NT1-MM
	LEL=LT+IPA+MM-1
	LELLE=LT+LT-2+IPA+MM
	RAP=-(SOLHH(LT,J,7,1)*(LELLE+1)+SOLHH(LT,J,8,1)*ERRE(NR1,1))/
     1	(DEN1(LEL,2,1)*(LELLE+1)+DEN2(LEL,2,1)*ERRE(NR1,1))
	DO LR=1,NR1
	POTEN(LR,LT,J,1)=POTEN(LR,LT,J,1)+RAP*SOLHH(LR,LEL,1,1)
	ENDDO
C
	ENDDO
	ENDDO
	ENDDO
	ENDIF
C
	IF (INDD(1).EQ.0.OR.INDD(1).EQ.2) THEN
	RETURN
	ELSE
C
	DO LT=1,NT1
	LEL=LT+LT-1+IPA
	RAP=(DEN1(LT,2,1)*LEL+DEN2(LT,2,1)*ERRE(NR1,1))/(DEN1(LT,4,1)*
     1	LEL+ERRE(NR1,1)*DEN2(LT,4,1))
C
	DO LR=1,NR1
	SOLHH(LR,LT,2,1)=SOLHH(LR,LT,1,1)-RAP*SOLHH(LR,LT,2,1)
	ENDDO
	ENDDO
C
	IF (NF.EQ.1) RETURN
C
	DO LT=1,NT1-1
	LEL=LT+LT+IPA
	RAP=(DEN1(LT,6,1)*LEL+DEN2(LT,6,1)*ERRE(NR1,1))/(DEN1(LT,8,1)*
     1	LEL+ERRE(NR1,1)*DEN2(LT,8,1))
	DO LR=1,NR1
	SOLHH(LR,LT,4,1)=SOLHH(LR,LT,3,1)-RAP*SOLHH(LR,LT,4,1)
	ENDDO
	ENDDO
	ENDIF
	ENDIF
C
	IF (IPAR.EQ.4.OR.IPAR.EQ.5) THEN
C
	DO LF=1,NF
	MM=(LF/2)
	DO LT=1,NT1-MM
	LEL=LT+MM
	LELLE=LT+LT-2+IPA
	RAP=-(SOLHH(LT,LF,5,1)*(LELLE+1)+SOLHH(LT,LF,6,1)*ERRE(NR1,1))/
     1	(DEN1(LEL,2,1)*(LELLE+1)+DEN2(LEL,2,1)*ERRE(NR1,1))
	DO LR=1,NR1
	POTEN(LR,LT,LF,1)=POTEN(LR,LT,LF,1)+RAP*SOLHH(LR,LEL,1,1)
	ENDDO
	ENDDO
	ENDDO
C
	IF (INDD(1).EQ.0.OR.INDD(1).EQ.2) THEN
	RETURN
	ELSE
C
	DO LT=1,NT1
	LEL=LT+LT-1+IPA
	RAP=(DEN1(LT,2,1)*LEL+DEN2(LT,2,1)*ERRE(NR1,1))/
     1	(DEN1(LT,4,1)*LEL+DEN2(LT,4,1)*ERRE(NR1,1))
	DO LR=1,NR1
	SOLHH(LR,LT,2,1)=SOLHH(LR,LT,1,1)-RAP*SOLHH(LR,LT,2,1)
	ENDDO
C
	ENDDO
	ENDIF
	ENDIF
	RETURN
	ENDIF
C
C*****************************************************************************
C
C		CAS NZON.GT.1 RACCORDEMENT ENTRE LES COQUILLE ET EVENTUELLEMNT
C		AVEC LE VIDE (SI INDD(NZON).NE.2)
C	
C		RELATION ENTRE LES VALEURS DE LF ET LE NOMBRE QUANTIQUE l
C
C		FORMATION DE LA MATRICE DANS LES CAS IPARR=0,4,5 NECESSAIRE
C		AU CALCUL DES COEFFICIENTS
C		NECESSAIRES AU COMBINAISON LINEAIRES POUR EFFECTUER LE
C		RACORDEMENT DE SOLUTIONS
C
	IF (IPA2.EQ.1) GO TO 999
C
	NZ2=NZON+NZON-2
	NR1=NDL1(2)
	LZ1=NZ2
	IF (INDD(NZON).NE.2) LZ1=NZ2+1
C
C		PREPARATION DE LAS MATRICE QUI CALCULE LES COEFFICIENTS
C		NECESSAIRES AU COMBINAISON LINEAIRES POUR EFFECTUER LE
C		RACORDEMENT
C
	JLZ=1
	IF (INDD(NZON).EQ.2.AND.IPARR.NE.5) JLZ=2
	DO LT=JLZ,NT1
	LELL1=LT
	DO LR=1,5
	DO J=1,NZ2+1
	BB(J,LR)=0
	ENDDO
	ENDDO
C
	BB(1,3)= DEN1(LT,2,1)
	BB(1,4)=-DEN1(LT,1,2)
	BB(1,5)=-DEN1(LT,3,2)
C
	BB(2,2)= DEN2(LT,2,1)
	BB(2,3)=-DEN2(LT,1,2)
	BB(2,4)=-DEN2(LT,3,2)
C
	DO LZON=2,NZON-1
	J=LZON+LZON-1
	BB(J,2)=DEN1(LT,2,LZON)
	BB(J,3)  =DEN1(LT,4,LZON)
	BB(J,4)=-DEN1(LT,1,LZON+1)
	BB(J,5)=-DEN1(LT,3,LZON+1)
C		
	BB(J+1,1)=DEN2(LT,2,LZON)
	BB(J+1,2)  =DEN2(LT,4,LZON)
	BB(J+1,3)=-DEN2(LT,1,LZON+1)
	BB(J+1,4)=-DEN2(LT,3,LZON+1)
	ENDDO
C
C
	IF (INDD(NZON).NE.2) THEN
C
C		CAS AVEC RACCORDEMENT AVEC LE VIDE
C
	IF (IPAR.GE.4) LELL1=LT+LT-1+IPA
C
C
	BB(NZ2+1,2)=DEN1(LT,2,NZON)*LELL1+DEN2(LT,2,NZON)*ERRE(NR1,NZON)
	BB(NZ2+1,3)=DEN1(LT,4,NZON)*LELL1+DEN2(LT,4,NZON)*ERRE(NR1,NZON)
	BB(NZ2+1,1)=0
	ENDIF
C
C		CAS IPARR=0,4,5
C
	IF (IPAR.EQ.0) THEN
	LF1=MIN(LT*2-1,NF)
C
	L=LT
	IPA=0
	DO LF=1,LF1
	LF2=LF/2
	LEL=LT-LF2+MOD(LF2,2)
C
	ILT(LF)=LEL
C
C		PREPARATION 2EME MEMBRE DU SYSTEME POUR LE CAS IPAR=0
C
	DO LZON=1,NZON-1
	J=LZON+LZON-1
	CC(J,LF)=-SOLHH(LEL,LF,5,LZON)+SOLHH(LEL,LF,3,LZON+1)
	CC(J+1,LF)=-SOLHH(LEL,LF,6,LZON)+SOLHH(LEL,LF,4,LZON+1)
	ENDDO
	ENDDO		
	ENDIF
C
C		CAS SUPERSYMMETRIQUE (IPAR=4,5)
C
	IF (IPA2.EQ.2) THEN
C
	L=LT+LT-1
	LF1=MIN(L,NF)
C
	DO LF=1,LF1
	LF2=(LF/2)*2
	LEL=(L-LF2)/2+1
	ILT(LF)=LEL
C
C		PREPARATION 2EME MEMBRE DU SYSTEME POUR LES CAS IPAR=4,5
C
	DO LZON=1,NZON-1
	J=LZON+LZON-1
	CC(J,LF)=-SOLHH(LEL,LF,5,LZON)+SOLHH(LEL,LF,3,LZON+1)
	CC(J+1,LF)=-SOLHH(LEL,LF,6,LZON)+SOLHH(LEL,LF,4,LZON+1)
	ENDDO
	ENDDO		
	ENDIF
C
	IF (INDD(NZON).NE.2) THEN
	DO LF=1,LF1
	CC(LZ1,LF)=-(SOLHH(ILT(LF),LF,5,NZON)*(L+IPA)+
     1	SOLHH(ILT(LF),LF,6,NZON)*ERRE(NR1,NZON))
	ENDDO
	ELSE
	DO LF=1,LF1+1
	CC(NZ2+1,LF)=0
	ENDDO
	ENDIF
C
	LF2=LF1+1
	DO LR=1,LZ1
	CC(LR,LF2)=0
	ENDDO
C
	CC(1,LF2)=-DEN1(LT,4,1)
	CC(2,LF2)=-DEN2(LT,4,1)
	DO LF=1,LF2
	ENDDO
C
	CALL LEQT1S(BB,LZ1,2,2,NDR,CC,LF2,NDR,0,UGRAV,J)
C
C		ON EFFECTUE LES COMBINAISON LINEAIRES POUR QUE LA SOLUTION
C		DE POISSON SOIT CONTINUE AVEC SES DERIVEES A TRAVERS LES
C		BORDS DES COQUILLES
C
	NR1=NDL1(2)
	NDEG(1)=NR1
	DO LF=1,LF1
	DO LR=1,NR1
	POTEN(LR,ILT(LF),LF,1)=POTEN(LR,ILT(LF),LF,1)+
     1	CC(1,LF)*SOLHH(LR,LT,1,1)
	ENDDO
	ENDDO
C
	DO LZON=2,NZON
	J=LZON+LZON-2
	NR1=NDL1(LZON+1)
	NDEG(1)=NR1
	DO LF=1,LF1
	DO LR=1,NR1
	POTEN(LR,ILT(LF),LF,LZON)=POTEN(LR,ILT(LF),LF,LZON)
     1	+CC(J,LF)*SOLHH(LR,LT,1,LZON)+CC(J+1,LF)*SOLHH(LR,LT,2,LZON)
	ENDDO
	ENDDO
	ENDDO
C
	NR1=NDL1(2)
C
	DO LR=1,NR1
	SOLHH(LR,LT,2,1)=SOLHH(LR,LT,2,1)+CC(1,LF2)*SOLHH(LR,LT,1,1)
	ENDDO
C
	DO LZON=2,NZON
	NR1=NDL1(LZON+1)
	NDEG(1)=NR1
	J=LZON+LZON-2
	DO LR=1,NR1
	SOLHH(LR,LT,2,LZON)=SOLHH(LR,LT,1,LZON)*CC(J,LF2)+
     1	CC(J+1,LF2)*SOLHH(LR,LT,2,LZON)
	ENDDO
	ENDDO
C
	ENDDO
C	
	IF (MOD(IPAR,2).EQ.0.AND.INDD(NZON).EQ.2) THEN
	NR1=NDL1(NZON+1)
	NR=NR1-1
	DO LR=1,NR1
	C64(LR)=POTEN(LR,1,1,NZON)
	ENDDO
	CALL EXTM1S(NR,1,0,C64,VA1(1))
	DO LZON=1,NZON
	POTEN(1,1,1,LZON)=POTEN(1,1,1,LZON)-2*VA1(1)
	ENDDO
	ENDIF
C
	RETURN
  999	CONTINUE
C
C
C		CAS SEULEMENT SYMMETRIQUE (IPAR=2,3)
C
C@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
	NZ2=NZON+NZON-2
	NR1=NDL1(2)
	LZ1=NZ2
	IF (INDD(NZON).NE.2) LZ1=NZ2+1
C
	DO LT=1,NT1
C
	LELL1=LT
	DO LR=1,5
	DO J=1,NZ2+1
	BB(J,LR)=0
	ENDDO
	ENDDO
C
	BB(1,3)= DEN1(LT,2,1)
	BB(1,4)=-DEN1(LT,1,2)
	BB(1,5)=-DEN1(LT,3,2)
C
	BB(2,2)= DEN2(LT,2,1)
	BB(2,3)=-DEN2(LT,1,2)
	BB(2,4)=-DEN2(LT,3,2)
C
	DO LZON=2,NZON-1
	J=LZON+LZON-1
	BB(J,2)=DEN1(LT,2,LZON)
	BB(J,3)  =DEN1(LT,4,LZON)
	BB(J,4)=-DEN1(LT,1,LZON+1)
	BB(J,5)=-DEN1(LT,3,LZON+1)
C		
	BB(J+1,1)=DEN2(LT,2,LZON)
	BB(J+1,2)  =DEN2(LT,4,LZON)
	BB(J+1,3)=-DEN2(LT,1,LZON+1)
	BB(J+1,4)=-DEN2(LT,3,LZON+1)
	ENDDO
C
	IF (INDD(NZON).NE.2) THEN
	LELL1=LT+LT-1+IPA
	BB(NZ2+1,2)=DEN1(LT,2,NZON)*LELL1+DEN2(LT,2,NZON)*ERRE(NR1,NZON)
	BB(NZ2+1,3)=DEN1(LT,4,NZON)*LELL1+DEN2(LT,4,NZON)*ERRE(NR1,NZON)
	BB(NZ2+1,1)=0
	ENDIF
C
	LF=1
	DO LZON=1,NZON-1
	J=LZON+LZON-1
	CC(J,1)=-SOLHH(LT,1,7,LZON)+SOLHH(LT,1,5,LZON+1)
	CC(J+1,1)=-SOLHH(LT,1,8,LZON)+SOLHH(LT,1,6,LZON+1)
	ENDDO		
	ILT(1)=LT
	LF1=1
C
	LF2=1
	L=LT+LT-1
	IF (LT.GT.1) THEN
	MM=L-1	
	LF1=MIN(2*MM+1,NF)
	DO LF=4,LF1,4
	LEL=LT-LF/4
	DO JJ=LF,MIN(LF+1,LF1)
	LF2=LF2+1
	ILT(LF2)=LEL
	DO LZON=1,NZON-1
	J=LZON+LZON-1
	CC(J,LF2)=-SOLHH(LEL,JJ,7,LZON)+SOLHH(LEL,JJ,5,LZON+1)
	CC(J+1,LF2)=-SOLHH(LEL,JJ,8,LZON)+SOLHH(LEL,JJ,6,LZON+1)
	ENDDO
	ENDDO		
	ENDDO
	ENDIF
C
	IF (INDD(NZON).NE.2) THEN
	CC(LZ1,1)=-(SOLHH(LT,1,7,NZON)*(L+IPA)+
     1	SOLHH(LT,1,8,NZON)*ERRE(NR1,NZON))
C
	LF2=1
	IF (LF1.GE.4) THEN
	DO LF=4,LF1,4
	DO JJ=LF,MIN(LF+1,LF1)
	LF2=LF2+1
	CC(LZ1,LF2)=-(SOLHH(ILT(LF2),JJ,7,NZON)*(L+IPA)+
     1	SOLHH(ILT(LF2),JJ,8,NZON)*ERRE(NR1,NZON))
	ENDDO
	ENDDO
	ENDIF
	ELSE
	DO LF=1,NF+1
	CC(LZ1+1,LF)=0
	ENDDO
	ENDIF
C
	J2=LF2+1
	DO J=1,LZ1
	CC(J,J2)=0
	ENDDO
	CC(1,J2)=-DEN1(LT,4,1)
	CC(2,J2)=-DEN2(LT,4,1)
C
C	
	CALL LEQT1S(BB,LZ1,2,2,NDR,CC,J2,NDR,0,UGRAV,J)
C
	NR1=NDL1(2)
	DO LR=1,NR1
	POTEN(LR,LT,1,1)=POTEN(LR,LT,1,1)+
     1	CC(1,1)*SOLHH(LR,LT,1,1)
	ENDDO
C
	LF2=1
	DO LZON=2,NZON
	J=LZON+LZON-2
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	POTEN(LR,LT,1,LZON)=POTEN(LR,LT,1,LZON)
     1	+CC(J,1)*SOLHH(LR,LT,1,LZON)+CC(J+1,1)*SOLHH(LR,LT,2,LZON)
	ENDDO
	ENDDO
C
	IF (LF1.GE.4) THEN
	NR1=NDL1(2)
	DO LF=4,LF1,4
	DO JJ=LF,MIN(LF+1,LF1)
	LF2=LF2+1
C
	DO LR=1,NR1
	POTEN(LR,ILT(LF2),JJ,1)=POTEN(LR,ILT(LF2),JJ,1)+
     1	CC(1,LF2)*SOLHH(LR,LT,1,1)
	ENDDO
C
	DO LZON=2,NZON
	J=LZON+LZON-2
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	POTEN(LR,ILT(LF2),JJ,LZON)=POTEN(LR,ILT(LF2),JJ,LZON)
     1	+CC(J,LF2)*SOLHH(LR,LT,1,LZON)+CC(J+1,LF2)*SOLHH(LR,LT,2,LZON)
	ENDDO
	ENDDO
	ENDDO
	ENDDO
	ENDIF
C
C		ON RACCORDE LES SOLUTIONS HOMOGENES
C
	NR1=NDL1(2)
	NR=NR1-1
	DO LR=1,NR1
	SOLHH(LR,LT,2,1)=SOLHH(LR,LT,2,1)+CC(1,J2)*SOLHH(LR,LT,1,1)
	ENDDO
C	
	DO LZON=2,NZON
	J=LZON+LZON-2
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	SOLHH(LR,LT,2,LZON)=SOLHH(LR,LT,1,LZON)*CC(J,J2)+
     1	CC(J+1,J2)*SOLHH(LR,LT,2,LZON)
	ENDDO
	ENDDO
	ENDDO
C
	DO LT=1,NT1-1
C
C		ON TRAITE LE CA m IMPAIRE
C
	DO LR=1,5
	DO J=1,NZ2+1
	BB(J,LR)=0
	ENDDO
	ENDDO
C
	BB(1,3)= DEN1(LT,6,1)
	BB(1,4)=-DEN1(LT,5,2)
	BB(1,5)=-DEN1(LT,7,2)
C
	BB(2,2)= DEN2(LT,6,1)
	BB(2,3)=-DEN2(LT,5,2)
	BB(2,4)=-DEN2(LT,7,2)
C
	DO LZON=2,NZON-1
	J=LZON+LZON-1
C
	BB(J,2)=DEN1(LT,6,LZON)
	BB(J,3)  =DEN1(LT,8,LZON)
	BB(J,4)=-DEN1(LT,5,LZON+1)
	BB(J,5)=-DEN1(LT,7,LZON+1)
C
	BB(J+1,1)=DEN2(LT,6,LZON)
	BB(J+1,2)  =DEN2(LT,8,LZON)
	BB(J+1,3)=-DEN2(LT,5,LZON+1)
	BB(J+1,4)=-DEN2(LT,7,LZON+1)
	ENDDO
C
	IF (INDD(NZON).NE.2) THEN
	LELL1=LT+LT+IPA
	BB(NZ2+1,2)=DEN1(LT,6,NZON)*LELL1+DEN2(LT,6,NZON)*ERRE(NR1,NZON)
	BB(NZ2+1,3)=DEN1(LT,8,NZON)*LELL1+DEN2(LT,8,NZON)*ERRE(NR1,NZON)
	BB(NZ2+1,1)=0
	ENDIF
C
	L=LT+LT-1
	LF1=MIN(2*L+1,NF)
C
	LF2=0
	DO LF=2,LF1,4
	MM=LF/2
	LEL=LT-LF/4
	DO JJ=LF,MIN(LF+1,LF1)
	LF2=LF2+1
	ILT(LF2)=LEL
	DO LZON=1,NZON-1
	J=LZON+LZON-1
	CC(J,LF2)=-SOLHH(LEL,JJ,7,LZON)+SOLHH(LEL,JJ,5,LZON+1)
	CC(J+1,LF2)=-SOLHH(LEL,JJ,8,LZON)+SOLHH(LEL,JJ,6,LZON+1)
	ENDDO
	ENDDO		
	ENDDO
C
	IF (INDD(NZON).NE.2) THEN
	LELL1=L+1+IPA
	LF2=0
	DO LF=2,LF1,4
	DO JJ=LF,MIN(LF+1,LF1)
	LF2=LF2+1
	CC(LZ1,LF2)=-(SOLHH(ILT(LF2),JJ,7,NZON)*LELL1+
     1	SOLHH(ILT(LF2),JJ,8,NZON)*ERRE(NR1,NZON))
	ENDDO
	ENDDO
	ELSE
	DO LF=1,NF+1
	CC(LZ1+1,LF)=0
	ENDDO
	ENDIF
C
	J2=LF2+1
	DO J=1,LZ1
	CC(J,J2)=0
	ENDDO
C
	CC(1,J2)=-DEN1(LT,8,1)
	CC(2,J2)=-DEN2(LT,8,1)
C
	CALL LEQT1S(BB,LZ1,2,2,NDR,CC,J2,NDR,0,UGRAV,J)
C
	NR1=NDL1(2)
	LF2=0
	DO LF=2,LF1,4
	DO JJ=LF,MIN(LF+1,LF1)
	LF2=LF2+1
C
	DO LR=1,NR1
	POTEN(LR,ILT(LF2),JJ,1)=POTEN(LR,ILT(LF2),JJ,1)+
     1	CC(1,LF2)*SOLHH(LR,LT,3,1)
	ENDDO
C
	DO LZON=2,NZON
	J=LZON+LZON-2
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	POTEN(LR,ILT(LF2),JJ,LZON)=POTEN(LR,ILT(LF2),JJ,LZON)
     1	+CC(J,LF2)*SOLHH(LR,LT,3,LZON)+CC(J+1,LF2)*SOLHH(LR,LT,4,LZON)
	ENDDO
	ENDDO
	ENDDO
	ENDDO
C
C		RACCORDEMENT FONCTIONS HOMOGENES
C
	NR1=NDL1(2)
	DO LR=1,NR1
	SOLHH(LR,LT,4,1)=SOLHH(LR,LT,4,1)+CC(1,J2)*SOLHH(LR,LT,3,1)
	ENDDO
	DO LZON=2,NZON
	J=LZON+LZON-2
	NR1=NDL1(LZON+1)
	DO LR=1,NR1
	SOLHH(LR,LT,4,LZON)=CC(J,J2)*SOLHH(LR,LT,3,LZON)+
     1	CC(J+1,J2)*SOLHH(LR,LT,4,LZON)
	ENDDO
	ENDDO
C
	ENDDO
C
	IF (IPARR.EQ.2.AND.INDD(NZON).EQ.2) THEN
	NR1=NDL1(NZON+1)
	NR=NR1-1
	DO LR=1,NR1
	C64(LR)=POTEN(LR,1,1,NZON)
	ENDDO
	CALL EXTM1S(NR,1,0,C64,VA1(1))
	DO LZON=1,NZON
	POTEN(1,1,1,LZON)=POTEN(1,1,1,LZON)-2*VA1(1)
	ENDDO
	ENDIF
C
C  100	FORMAT(1X,10E10.3)
  101	FORMAT(1X,' ')
C  200	FORMAT(1X,16I4)
	RETURN
	END
C

C