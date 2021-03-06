/*
 *  Methods of class Eos_CompOSE
 *
 *  (see file eos_compose.h for documentation).
 *
 */

/*
 *   Copyright (c) 2010, 2014 Jerome Novak
 *             (c) 2019 Lorenzo Sala
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
 * $Id: eos_compose.C,v 1.10 2019/03/28 13:41:02 j_novak Exp $
 * $Log: eos_compose.C,v $
 * Revision 1.10  2019/03/28 13:41:02  j_novak
 * Improved managed of saved EoS (functions sauve and constructor form FILE*)
 *
 * Revision 1.9  2019/03/25 14:21:24  j_novak
 * Improved constructor from a FILE*, following patch by Lorenzo Sala.
 *
 * Revision 1.8  2016/12/05 16:17:51  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2015/12/04 16:27:05  j_novak
 * Correction of constructor calling.
 *
 * Revision 1.6  2015/08/04 14:41:29  j_novak
 * Back to previous version for Eos_CompOSE. Enthalpy-consistent EoS can be accessed using Eos_consistent class (derived from Eos_CompOSE).
 *
 * Revision 1.5  2015/01/27 14:22:38  j_novak
 * New methods in Eos_tabul to correct for EoS themro consistency (optional).
 *
 * Revision 1.4  2014/10/13 08:52:52  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/07/01 09:26:21  j_novak
 * Improvement of comments
 *
 * Revision 1.2  2014/06/30 16:13:18  j_novak
 * New methods for reading directly from CompOSE files.
 *
 * Revision 1.1  2014/03/06 15:53:34  j_novak
 * Eos_compstar is now Eos_compOSE. Eos_tabul uses strings and contains informations about authors.
 *
 * Revision 1.1  2010/02/03 14:56:45  j_novak
 * *** empty log message ***
 *
 *
 * $Header: /cvsroot/Lorene/C++/Source/Eos/eos_compose.C,v 1.10 2019/03/28 13:41:02 j_novak Exp $
 *
 */

#include <string>

// Headers Lorene
#include "headcpp.h"
#include "eos.h"
#include "tbl.h"
#include "utilitaires.h"
#include "unites.h"

			//----------------------------//
			//   	Constructors	      //
			//----------------------------//

// Standard constructor
// --------------------			
namespace Lorene {

// Standard constructor with only the name of the file
//------------------------------------------------------
  Eos_CompOSE::Eos_CompOSE(const char* file_name)
    : Eos_tabul("CompOSE EoS"), format(0) {
    tablename = file_name ;
  }
  

// Constructor from binary file
// ----------------------------
  Eos_CompOSE::Eos_CompOSE(FILE* fich) : Eos_tabul("CompOSE Eos"), format(0) {
    
    fread(name, sizeof(char), 100, fich) ;
    char tmp_string[160] ;
    fread(tmp_string, sizeof(char), 160, fich) ;
    tablename = tmp_string;

    int format_read = int(fread_be(&format, sizeof(int), 1, fich)) ;

    if (format_read != 1) format = 0 ;

    if (format == 0) read_table() ;
    else read_compose_data() ;
  }
  


// Constructor from a formatted file
// ---------------------------------
  Eos_CompOSE::Eos_CompOSE(ifstream& fich) : Eos_tabul(fich), format(0)
  {}


// Constructor from CompOSE data files
// ------------------------------------
  Eos_CompOSE::Eos_CompOSE(const string& files) :
    Eos_tabul("CompOSE Eos"), format(1) {
    
    tablename = files ;
    read_compose_data() ;
  }


// Reading the data files in CompOSE format
//-----------------------------------------
  void Eos_CompOSE::read_compose_data() {

    using namespace Unites ;
    
    // Files containing data and a test
    //---------------------------------
    string file_nb = tablename + ".nb" ;
    string file_thermo = tablename + ".thermo" ;

    ifstream in_nb(file_nb.data()) ;
    if (!in_nb) {
      cerr << "Eos_CompOSE::read_compose_data(): " << endl ;
      cerr << "Problem in opening the EOS file!" << endl ;
      cerr << "While trying to open " << file_nb << endl ;
      cerr << "Aborting..." << endl ;
      abort() ;
    }

    // obtaining the size of the tables for memory allocation
    //-------------------------------------------------------
    int index1, index2 ;
    in_nb >> index1 >> index2 ;
    int nbp = index2 - index1 + 1 ;
    assert(nbp > 0) ;

    press = new double[nbp] ;
    nb    = new double[nbp] ;
    ro    = new double[nbp] ;
    
    logh = new Tbl(nbp) ;
    logp = new Tbl(nbp) ;
    dlpsdlh = new Tbl(nbp) ;
    lognb = new Tbl(nbp) ;
    dlpsdlnb = new Tbl(nbp) ;
    
    logh->set_etat_qcq() ;
    logp->set_etat_qcq() ;
    dlpsdlh->set_etat_qcq() ;
    lognb->set_etat_qcq() ;
    dlpsdlnb->set_etat_qcq() ;
    
    // Variables and conversion
    //-------------------------
    double nb_fm3, rho_cgs, p_cgs, p_over_nb_comp, eps_comp ;
    double dummy_x ;
    int dummy_n ;
    
    double rhonuc_cgs = rhonuc_si * 1e-3 ;
    double c2_cgs = c_si * c_si * 1e4 ;
    double m_neutron_MeV, m_proton_MeV ;
    
    ifstream in_p_rho (file_thermo.data()) ;
    if (!in_p_rho) {
      cerr << "Eos_CompOSE::read_compose_data(): " << endl ;
      cerr << "Problem in opening the EOS file!" << endl ;
      cerr << "While trying to open " << file_thermo << endl ;
      cerr << "Aborting..." << endl ;
      abort() ;
    }
    in_p_rho >> m_neutron_MeV >> m_proton_MeV ; //Neutron and proton masses
    in_p_rho.ignore(1000, '\n') ;
    
    double p_convert = mev_si * 1.e45 * 10. ; // Conversion from MeV/fm^3 to cgs
    double eps_convert = mev_si * 1.e42 / (c_si*c_si) ; //From meV/fm^3 to g/cm^3
    
    // Main loop reading the table
    //----------------------------
    for (int i=0; i<nbp; i++) {
      in_nb >> nb_fm3 ;
      in_p_rho >> dummy_n >> dummy_n >> dummy_n >> p_over_nb_comp ;
      in_p_rho >> dummy_x >> dummy_x >> dummy_x >> dummy_x >> dummy_x >> eps_comp ;
      in_p_rho.ignore(1000, '\n') ;
      p_cgs = p_over_nb_comp * nb_fm3 * p_convert ;
      rho_cgs = ( eps_comp + 1. ) * m_neutron_MeV * nb_fm3 * eps_convert ;
      
      if ( (nb_fm3<0) || (rho_cgs<0) || (p_cgs < 0) ){
	cout << "Eos_CompOSE::read_compose_data(): " << endl ;
	cout << "Negative value in table!" << endl ;
	cout << "nb = " << nb_fm3 << ", rho = " << rho_cgs <<
	  ", p = " << p_cgs << endl ;
	cout << "Aborting..." << endl ;
	abort() ;
      }
      
      press[i] = p_cgs / c2_cgs ;
      nb[i]    = nb_fm3 ;
      ro[i]    = rho_cgs ;
    }
    
    double ww = 0. ;
    for (int i=0; i<nbp; i++) {
      double h = log( (ro[i] + press[i]) / (10 * nb[i] * rhonuc_cgs) ) ;
      
      if (i==0) { ww = h ; }
      h = h - ww + 1.e-14 ;
      
      logh->set(i) = log10( h ) ;
      logp->set(i) = log10( press[i] / rhonuc_cgs ) ;
      dlpsdlh->set(i) = h * (ro[i] + press[i]) / press[i] ;
      lognb->set(i) = log10(nb[i]) ;
    }
    
    // Computation of dpdnb
    //---------------------
    double p0, p1, p2, n0, n1, n2, dpdnb;
    
    // special case: i=0
    p0 = log(press[0]);
    p1 = log(press[1]);
    p2 = log(press[2]);
    
    n0 = log(nb[0]);
    n1 = log(nb[1]);
    n2 = log(nb[2]);
    
    dpdnb = p0*(2*n0-n1-n2)/(n0-n1)/(n0-n2) +
      p1*(n0-n2)/(n1-n0)/(n1-n2) +
      p2*(n0-n1)/(n2-n0)/(n2-n1) ;
    
    dlpsdlnb->set(0) = dpdnb ;
    
    for(int i=1;i<nbp-1;i++) {
      
      p0 = log(press[i-1]);
      p1 = log(press[i]);
      p2 = log(press[i+1]);
      
      n0 = log(nb[i-1]);
      n1 = log(nb[i]);
      n2 = log(nb[i+1]);
      
      dpdnb = p0*(n1-n2)/(n0-n1)/(n0-n2) +
	p1*(2*n1-n0-n2)/(n1-n0)/(n1-n2) +
	p2*(n1-n0)/(n2-n0)/(n2-n1) ;
      
      dlpsdlnb->set(i) = dpdnb ;
      
    }

    // special case: i=nbp-1
    p0 = log(press[nbp-3]);
    p1 = log(press[nbp-2]);
    p2 = log(press[nbp-1]);
    
    n0 = log(nb[nbp-3]);
    n1 = log(nb[nbp-2]);
    n2 = log(nb[nbp-1]);
    
    dpdnb = p0*(n2-n1)/(n0-n1)/(n0-n2) +
      p1*(n2-n0)/(n1-n0)/(n1-n2) +
      p2*(2*n2-n0-n1)/(n2-n0)/(n2-n1) ;
    
    dlpsdlnb->set(nbp-1) = dpdnb ;
    
    hmin = pow( double(10), (*logh)(0) ) ;
    hmax = pow( double(10), (*logh)(nbp-1) ) ;
    
    // Cleaning
    //---------
    delete [] press ;
    delete [] nb ;
    delete [] ro ;

  }
  
			//--------------//
			//  Destructor  //
			//--------------//

  Eos_CompOSE::~Eos_CompOSE(){

    // does nothing

  }


			//------------------------//
			//  Comparison operators  //
			//------------------------//


  bool Eos_CompOSE::operator==(const Eos& eos_i) const {

    bool resu = true ;
    
    if ( eos_i.identify() != identify() ) {
      cout << "The second EOS is not of type Eos_CompOSE !" << endl ;
      resu = false ;
    }
    
    return resu ;
    
  }

  bool Eos_CompOSE::operator!=(const Eos& eos_i) const {
    
    return !(operator==(eos_i)) ;
    
  }

			//------------//
			//  Outputs   //
			//------------//
  void Eos_CompOSE::sauve(FILE* fich) const {

    Eos_tabul::sauve(fich) ;

    fwrite_be(&format, sizeof(int), 1, fich) ;
  }
  

ostream& Eos_CompOSE::operator>>(ostream & ost) const {

    ost << "EOS of class Eos_CompOSE." << endl ;
    ost << "Built from file " << tablename << endl ;
    ost << ((format == 0) ? "Old LORENE format" : "CompOSE format") << endl ; 
    ost << "Authors : " << authors << endl ;
    ost << "Number of points in file : " << logh->get_taille() << endl ;
    return ost ;

}

			
}
