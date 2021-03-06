/*
 * Methods of the class Eos_incomp.
 *
 * (see file eos.h for documentation).
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */


 

/*
 * $Id: eos_incomp.C,v 1.7 2016/12/05 16:17:51 j_novak Exp $
 * $Log: eos_incomp.C,v $
 * Revision 1.7  2016/12/05 16:17:51  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:52:53  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:13:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2002/10/16 14:36:35  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.3  2002/04/09 14:32:15  e_gourgoulhon
 * 1/ Added extra parameters in EOS computational functions (argument par)
 * 2/ New class MEos for multi-domain EOS
 *
 * Revision 1.2  2001/12/04 21:27:53  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2001/02/07  09:49:08  eric
 * Suppression de la fonction derent_ent_p.
 * Ajout des fonctions donnant les derivees de l'EOS:
 *      der_nbar_ent_p
 *      der_ener_ent_p
 *      der_press_ent_p
 *
 * Revision 2.3  2000/02/14  14:49:48  eric
 * Modif affichage.
 *
 * Revision 2.2  2000/02/14  14:33:51  eric
 * Ajout du constructeur par lecture de fichier formate.
 *
 * Revision 2.1  2000/01/21  15:18:15  eric
 * Ajout des operateurs de comparaison == et !=
 *
 * Revision 2.0  2000/01/18  16:11:25  eric
 * *** empty log message ***
 *
 *
 * $Header: /cvsroot/Lorene/C++/Source/Eos/eos_incomp.C,v 1.7 2016/12/05 16:17:51 j_novak Exp $
 *
 */


// Headers C
#include <cstdlib>
#include <cstring>
#include <cmath>

// Headers Lorene
#include "eos.h"
#include "cmp.h"
#include "utilitaires.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor with ent0 = 1
// ---------------------------------
namespace Lorene {
Eos_incomp::Eos_incomp(double rho_c) : 
	Eos("EOS for relativistic incompressible matter"), 
	rho0(rho_c), ent0( double(-1.e-6) ) {}

// Standard constructor with ent0 specified
// ---------------------------------------
Eos_incomp::Eos_incomp(double rho_c, double ent_c) : 
	Eos("EOS for relativistic incompressible matter"), 
	rho0(rho_c), ent0( ent_c ) {

    assert( ent_c <= double(0) ) ; 
	
}
  
// Copy constructor
// ----------------
Eos_incomp::Eos_incomp(const Eos_incomp& eosi) : 
	Eos(eosi), 
	rho0(eosi.rho0), ent0(eosi.ent0) {}
  

// Constructor from a binary file
// ------------------------------
Eos_incomp::Eos_incomp(FILE* fich) : 
	Eos(fich) {
        
    fread_be(&rho0, sizeof(double), 1, fich) ;		
    fread_be(&ent0, sizeof(double), 1, fich) ;		

}

// Constructor from a formatted file
// ---------------------------------
Eos_incomp::Eos_incomp(ifstream& fich) : 
	Eos(fich) {
        
    char blabla[80] ;
        
    fich >> rho0 ; fich.getline(blabla, 80) ;
    fich >> ent0 ; fich.getline(blabla, 80) ;

}
			//--------------//
			//  Destructor  //
			//--------------//

Eos_incomp::~Eos_incomp(){
    
    // does nothing
        
}
			//--------------//
			//  Assignment  //
			//--------------//

void Eos_incomp::operator=(const Eos_incomp& eosi) {
    
    set_name(eosi.name) ; 
    
    rho0 = eosi.rho0 ; 
    ent0 = eosi.ent0 ; 

}

			//------------------------//
			//  Comparison operators  //
			//------------------------//

bool Eos_incomp::operator==(const Eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_incomp !" << endl ; 
	resu = false ; 
    }
    else{
	
	const Eos_incomp& eos = dynamic_cast<const Eos_incomp&>( eos_i ) ; 

	if (eos.rho0 != rho0) {
	    cout 
	    << "The two Eos_incomp have different rho0 : " << rho0 << " <-> " 
		<< eos.rho0 << endl ; 
	    resu = false ; 
	}

	if (eos.ent0 != ent0) {
	    cout 
	    << "The two Eos_incomp have different ent0 : " << ent0 << " <-> " 
		<< eos.ent0 << endl ; 
	    resu = false ; 
	}

    }
    
    return resu ; 
    
}

bool Eos_incomp::operator!=(const Eos& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}



			//------------//
			//  Outputs   //
			//------------//

void Eos_incomp::sauve(FILE* fich) const {

    Eos::sauve(fich) ; 
    
    fwrite_be(&rho0, sizeof(double), 1, fich) ;	
    fwrite_be(&ent0, sizeof(double), 1, fich) ;	
   
}

ostream& Eos_incomp::operator>>(ostream & ost) const {
    
    ost << "EOS of class Eos_incomp (relativistic incompressible matter) : " 
	<< endl ; 
    ost << "   Constant density : " << rho0 << " rho_nuc" << endl ; 
    ost << "   Log-enthalpy threshold for non-zero density : " << ent0 
	<< " c^2" <<  endl ; 
    
    return ost ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy 
//------------------------------

double Eos_incomp::nbar_ent_p(double ent, const Param* ) const {
    
    if ( ent >= ent0 ) {

	return rho0 ;
    }
    else{
	return 0 ;
    }
}

// Energy density from enthalpy 
//------------------------------

double Eos_incomp::ener_ent_p(double ent, const Param* ) const {
    
    if ( ent >= ent0 ) {

	return rho0 ;
    }
    else{
	return 0 ;
    }
}

// Pressure from enthalpy 
//------------------------

double Eos_incomp::press_ent_p(double ent, const Param* ) const {
    
    if ( ent >= ent0 ) {

	return rho0 * (exp(ent) - double(1)) ;
    }
    else{
	return 0 ;
    }
}


// dln(n)/ln(H) from enthalpy 
//---------------------------

double Eos_incomp::der_nbar_ent_p(double ent, const Param* ) const {
    
    if ( ent >= ent0 ) {
    
	return 0 ;
    }
    else{
	return 0 ;
    }
}

// dln(e)/ln(H) from enthalpy 
//---------------------------

double Eos_incomp::der_ener_ent_p(double ent, const Param* ) const {
    
    if ( ent >= ent0 ) {

	return 0 ;
    }
    else{
	return 0 ;
    }
}

// dln(p)/ln(H) from enthalpy 
//---------------------------

double Eos_incomp::der_press_ent_p(double ent, const Param* ) const {
    
    if ( ent >= ent0 ) {

	return ent / (double(1) - exp(-ent)) ;

    }
    else{
	return 0 ;
    }
}



}
