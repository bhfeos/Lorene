/*
 * Main code for computing stationary axisymmetric rotating stars. 
 * 
 */

/*
 *   Copyright (c) 2010 Eric Gourgoulhon
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
 * $Id: nrotstar.C,v 1.17 2017/11/06 11:03:51 m_bejger Exp $
 * $Log: nrotstar.C,v $
 * Revision 1.17  2017/11/06 11:03:51  m_bejger
 * Update the gyoto output file: 4-acceleration on the surface
 *
 * Revision 1.16  2017/10/20 13:56:20  j_novak
 * Adapted to be run with nt=1.
 *
 * Revision 1.15  2017/04/11 12:58:54  m_bejger
 * Updated output file for gyoto
 *
 * Revision 1.14  2016/12/05 16:18:25  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.13  2014/10/13 08:53:58  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.12  2013/12/15 17:48:15  e_gourgoulhon
 * Added output file for GYOTO
 *
 * Revision 1.11  2012/01/02 13:57:58  j_novak
 * Reverting to version 1.9
 *
 * Revision 1.7  2010/03/29 14:36:13  e_gourgoulhon
 * Change of name of the output file calcul.d --> result.txt
 *
 * Revision 1.6  2010/03/29 14:07:27  e_gourgoulhon
 * Added file outputs prof_*.d for various radial profile plots.
 *
 * Revision 1.5  2010/01/31 18:27:38  e_gourgoulhon
 * Better reading of parameter files.
 *
 * Revision 1.4  2010/01/26 16:52:37  e_gourgoulhon
 * Added the graphical outputs at the end.
 *
 * Revision 1.3  2010/01/26 10:54:07  e_gourgoulhon
 * First complete version.
 *
 * Revision 1.2  2010/01/25 22:34:15  e_gourgoulhon
 * On avance...
 *
 * Revision 1.1  2010/01/25 18:16:39  e_gourgoulhon
 * First version.
 *
 *
 * $Header: /cvsroot/Lorene/Codes/Nrotstar/nrotstar.C,v 1.17 2017/11/06 11:03:51 m_bejger Exp $
 *
 */

// C headers
#include <cstdlib>
#include <cmath>
#include <cstring>

// Lorene headers
#include "star_rot.h"
#include "eos.h"
#include "utilitaires.h"
#include "graphique.h"
#include "nbr_spx.h"
#include "unites.h"	    

namespace Lorene{
// Local prototype (for drawings only)
Scalar raccord_c1(const Scalar& uu, int l1) ; 
}
//******************************************************************************

using namespace Lorene ;

int main(){

    using namespace Unites ; 

    //------------------------------------------------------------------
    //	    Parameters of the computation 
    //------------------------------------------------------------------

    int relat_i, mer_max, mer_rot, mer_change_omega, mer_fix_omega, 
	delta_mer_kep, mer_mass, mermax_poisson, graph, nz, nzet, nzadapt,
	nt, np, mer_triax ; 
    double ent_c, freq_si, fact_omega, mbar_wanted, precis, freq_ini_si, 
	   thres_adapt, aexp_mass, relax, relax_poisson, ampli_triax, 
	   precis_adapt ;  
    
    ifstream fpar("par_rot.d") ;
    if ( !fpar.good() ) {
        cerr << "Problem in opening the file par_rot.d ! " << endl ;
        abort() ;
    }
    

    fpar.ignore(1000,'\n') ;    // skip title
    fpar >> relat_i ; fpar.ignore(1000,'\n') ;
    bool relat = (relat_i == 1) ; 
    fpar >> ent_c ; fpar.ignore(1000,'\n') ;
    fpar >> freq_si ; fpar.ignore(1000,'\n') ;
    fpar >> fact_omega ; fpar.ignore(1000,'\n') ;
    fpar >> mbar_wanted ; fpar.ignore(1000,'\n') ;
    mbar_wanted *= msol ; 
    fpar.ignore(1000,'\n') ;	// skip title
    fpar >> mer_max ; fpar.ignore(1000,'\n') ;
    fpar >> precis ; fpar.ignore(1000,'\n') ;
    fpar >> mer_rot ; fpar.ignore(1000,'\n') ;
    fpar >> freq_ini_si ; fpar.ignore(1000,'\n') ;
    fpar >> mer_change_omega ; fpar.ignore(1000,'\n') ;
    fpar >> mer_fix_omega ; fpar.ignore(1000,'\n') ;
    fpar >> delta_mer_kep ; fpar.ignore(1000,'\n') ;
    fpar >> thres_adapt ; fpar.ignore(1000,'\n') ;
    fpar >> mer_triax ; fpar.ignore(1000,'\n') ;
    fpar >> ampli_triax ; fpar.ignore(1000,'\n') ;
    fpar >> mer_mass ; fpar.ignore(1000,'\n') ;
    fpar >> aexp_mass ; fpar.ignore(1000,'\n') ;
    fpar >> relax ; fpar.ignore(1000,'\n') ;
    fpar >> mermax_poisson ; fpar.ignore(1000,'\n') ;
    fpar >> relax_poisson ; fpar.ignore(1000,'\n') ;
    fpar >> precis_adapt ; fpar.ignore(1000,'\n') ;
    fpar >> graph ; fpar.ignore(1000,'\n') ;
    fpar.ignore(1000,'\n') ; // skip title
    fpar >> nz ; fpar.ignore(1000,'\n') ;
    fpar >> nzet; fpar.ignore(1000,'\n') ;
    fpar >> nzadapt; fpar.ignore(1000,'\n') ;
    fpar >> nt; fpar.ignore(1000,'\n') ;
    fpar >> np; fpar.ignore(1000,'\n') ;

    int* nr = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
     
    fpar.ignore(1000,'\n'); 	// skip title
    for (int l=0; l<nz; l++) {
	fpar >> nr[l]; 
	fpar >> bornes[l]; fpar.ignore(1000,'\n') ;
	np_tab[l] = np ; 
	nt_tab[l] = nt ; 
    }
    bornes[nz] = __infinity ;

    Tbl ent_limit(nzet) ;
    ent_limit.set_etat_qcq() ;
    ent_limit.set(nzet-1) = 0 ; 	// enthalpy at the stellar surface
    for (int l=0; l<nzet-1; l++) {
    	fpar >> ent_limit.set(l) ; fpar.ignore(1000,'\n') ;
    }

    fpar.close();

    // Particular cases
    // ----------------

    // Initial frequency = final frequency
    if ( freq_ini_si < 0 ) {
	freq_ini_si = freq_si ; 
	mer_change_omega = mer_rot ; 
	mer_fix_omega = mer_rot + 1 ;  
    }

    
    //-----------------------------------------------------------------------
    //		Equation of state
    //-----------------------------------------------------------------------

    fpar.open("par_eos.d") ;
    if ( !fpar.good() ) {
        cerr << "Problem in opening the file par_eos.d ! " << endl ;
        abort() ;
    }

    Eos* peos = Eos::eos_from_file(fpar) ;
    Eos& eos = *peos ;

    fpar.close() ;


    // Special treatment of crust - liquid core boundary in the case
    //  of Eos_strange
    if (eos.identify() == 6) {
    	assert( nzet == 2 ) ;    	
    	const Eos_strange_cr* peos_cr = dynamic_cast<const Eos_strange_cr*>(peos) ;
    	if (peos_cr == 0x0) {
    	       cout << "nrotstar: problem : peos is not of type Eos_strange_cr !" << endl ;
    	       abort() ;
    	}
    	
    	ent_limit.set(0) = peos_cr->get_ent_nd() ;  // enthalpy at core/crust transition

    }

    //-----------------------------------------------------------------------
    //		Construction of the multi-grid and the mapping
    //-----------------------------------------------------------------------

    // Rescale of bornes in the case where there more than 1 domain inside
    //   the star

    for (int l=0; l<nzet-1; l++) {

    	bornes[l+1] = bornes[nzet] * sqrt(1 - ent_limit(l) / ent_c) ;

    }

    // Type of r sampling :
    int* type_r = new int[nz];
    type_r[0] = RARE ; 
    for (int l=1; l<nz-1; l++) {
	type_r[l] = FIN ; 
    }
    type_r[nz-1] = UNSURR ; 
    
    // Type of sampling in theta and phi :
    int type_t = SYM ; 
    int type_p = SYM ; 
    
    Mg3d mg(nz, nr, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_et mp(mg, bornes) ;
   
    // Cleaning
    // --------

    delete [] nr ; 
    delete [] nt_tab ; 
    delete [] np_tab ; 
    delete [] type_r ; 
    delete [] bornes ; 
       


    cout << endl 
	 << "==========================================================" << endl
	 << "                    Physical parameters                   " << endl
	 << "=========================================================="
	 << endl ; 
    cout << endl ;

    cout << endl << "Equation of state : " 
	 << endl << "=================   " << endl ;
    cout << eos << endl ; 

    cout << "Central enthalpy : " << ent_c << " c^2" << endl ; 
    cout << "Rotation frequency : " << freq_si << " Hz" << endl ; 
    if ( abs(mer_mass) < mer_max ) {
	cout << "Required Baryon mass [M_sol] : " 
	     << mbar_wanted / msol << endl ; 
    }
    
    cout << endl 
	 << "==========================================================" << endl
	 << "               Computational parameters                   " << endl
	 << "=========================================================="
	 << endl << endl ; 

    cout << "Maximum number of steps in the main iteration : " 
	 << mer_max << endl ; 
    cout << "Relaxation factor in the main iteration  : " 
	 << relax << endl ; 
    cout << "Threshold on the enthalpy relative change for ending the computation : " 
	 << precis << endl ; 
    cout << "Maximum number of steps in Map_et::poisson : " 
	 << mermax_poisson << endl ; 
    cout << "Relaxation factor in Map_et::poisson : " 
	 << relax_poisson << endl ; 
    cout << "Step from which the baryon mass is forced to converge : " 
	 << mer_mass << endl ; 
    cout << "Exponent for the increase factor of the central enthalpy : " 
	 << aexp_mass << endl ; 
    cout << 
    "Threshold on |dH/dr|_eq / |dH/dr|_pole for the adaptation of the mapping"
    << endl << thres_adapt << endl ; 


    cout << endl << "Multi-grid : " 
	 << endl << "==========" << endl << mg << endl ; 
    cout << "Mapping : " 
	 << endl << "=======" << endl << mp << endl ; 


    //-----------------------------------------------------------------------
    //		Construction of the star
    //-----------------------------------------------------------------------
    
    Star_rot star(mp, nzet, relat, eos) ; 
    
    if ( star.is_relativistic() ) {
	cout << "========================" << endl ;
	cout << "Relativistic computation" << endl ;
	cout << "========================" << endl ;
    }
    else {
	cout << "=====================" << endl ;
	cout << "Newtonian computation" << endl ;
	cout << "=====================" << endl ;
    }
 
    //-----------------------------------------------------------------------
    //		Initialization of the enthalpy field
    //-----------------------------------------------------------------------


    const Coord& r = mp.r ;
    double ray0 = mp.val_r(nzet-1, 1., 0., 0.) ;  
    Scalar ent0(mp) ; 
    ent0 = ent_c * ( 1 - r*r / (ray0*ray0) ) ; 
    ent0.annule_domain(nz-1) ; 
    ent0.std_spectral_base() ; 
    star.set_enthalpy(ent0) ;  
    
    // Initialization of (n,e,p) from H
    star.equation_of_state() ; 

    // Initialization of (E,S,U,etc...) (quantities relative to the Eulerian obs)
    star.hydro_euler() ; 

    cout << endl << "Initial star : " 
	 << endl << "============   " << endl ;

    cout << star << endl ; 
     
    //-----------------------------------------------------------------------
    //		Computation of the rotating equilibrium
    //-----------------------------------------------------------------------

    double omega = 2 * M_PI * freq_si / f_unit ; 
    double omega_ini = 2 * M_PI * freq_ini_si / f_unit ; 

    Itbl icontrol(8) ;
    icontrol.set_etat_qcq() ; 
    icontrol.set(0) = mer_max ; 
    icontrol.set(1) = mer_rot ; 
    icontrol.set(2) = mer_change_omega ; 
    icontrol.set(3) = mer_fix_omega ; 
    icontrol.set(4) = mer_mass ; 
    icontrol.set(5) = mermax_poisson ; 
    icontrol.set(6) = mer_triax ; 
    icontrol.set(7) = delta_mer_kep ; 
    
    Tbl control(7) ; 
    control.set_etat_qcq() ; 
    control.set(0) = precis ; 
    control.set(1) = omega_ini ; 
    control.set(2) = relax ; 
    control.set(3) = relax_poisson ; 
    control.set(4) = thres_adapt ; 
    control.set(5) = ampli_triax ; 
    control.set(6) = precis_adapt ; 

    Tbl diff(8) ;     

    if (nt > 1) 
      star.equilibrium(ent_c, omega, fact_omega, nzadapt, ent_limit,
		       icontrol, control, mbar_wanted, aexp_mass, diff) ;
    else {
      cout << "************** Warning ****************" << endl ;
      cout <<" with nt = 1 some features are not available" << endl ;
      cout << "(e.g. convergence to a given mass)." << endl ;
      star.equilibrium_spher(ent_c, precis, &ent_limit) ;
    }

    cout << endl << "Final star : " 
	 << endl << "==========   " << endl ;

    cout.precision(10) ; 
    cout << star << endl ;

    double rho_c = star.get_ener().val_grid_point(0,0,0,0) ;

    cout << "r_p/r_eq :" << star.aplat() << endl ;
    cout << "Omega rho0^{-1/2} : " << star.get_omega_c() /
    	sqrt( ggrav * rho_c ) << endl ;

    cout << "M rho0^{1/2} : " << star.mass_g() * pow(ggrav,1.5) *
    	sqrt( rho_c ) << endl ;

    cout << "M_B rho0^{1/2} : " << star.mass_b() * pow(ggrav,1.5) *
    	sqrt( rho_c ) << endl ;

    cout << "R_circ rho0^{1/2} : " << star.r_circ() *
    	sqrt( ggrav * rho_c ) << endl ;
    	
    cout << "J rho0 : " << star.angu_mom() * ggrav * ggrav * rho_c  << endl ;
    	
    cout << "Z_p :      " << star.z_pole() << endl ;
    cout << "Z_eq^f :   " << star.z_eqf() << endl ;
    cout << "Z_eq^b :   " << star.z_eqb() << endl ;

    cout << "GRV2: " << star.grv2() << endl ;
    cout << "GRV3: " << star.grv3() << endl ;

    double vit_triax = 0. ;
    if (nt > 1) vit_triax = diff(7) ;

    //-----------------------------------------------
    //  General features of the final configuration
    //  saved in a file
    //-----------------------------------------------

    ofstream fichfinal("result.txt") ;
    fichfinal.precision(10) ; 
    
    if ( star.is_relativistic() ) {
	fichfinal << "Relativistic computation" << endl ;
    }
    else {
	fichfinal << "Newtonian computation" << endl ;
    }
    
    fichfinal << star.get_eos() << endl ;
    
    fichfinal << endl << "Total CPU time  : " << endl ;
    fichfinal << "Memory size : " << endl << endl ; 

    fichfinal << endl << endl ; 
    fichfinal << "Grid : " << endl ; 
    fichfinal << "------ " << endl ; 
    fichfinal << *(star.get_mp().get_mg()) << endl ; 
    fichfinal << endl << "Physical characteristics : " << endl ; 
    fichfinal	  << "-------------------------" << endl ; 
    fichfinal << star << endl ;
    fichfinal << "Growing rate of triaxial perturbation: " << vit_triax 
	      << endl ; 

    fichfinal << endl <<
    "===================================================================" 
    << endl ;
    if (nt > 1) 
      fichfinal << "Diff_ent : " << diff(0) << endl ; 
    fichfinal << "Relative error on the virial theorem GRV2 : "
	      << star.grv2() << endl ;   
    fichfinal << "Relative error on the virial theorem GRV3 : "
	      << star.grv3() << endl ;   
    
    fichfinal << endl <<
    "================================================================" << endl ;
    fichfinal <<
    "   PARAMETERS USED FOR THE COMPUTATION (file parrot.d) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    system("cat par_rot.d >> result.txt") ; 

    fichfinal.open("result.txt", ios::app) ;
    fichfinal << endl <<
    "================================================================" << endl ;
    fichfinal <<
    "	           EOS PARAMETERS (file par_eos.d) : " << endl ;
    fichfinal <<
    "================================================================" << endl ;
    fichfinal.close() ;
    system("cat par_eos.d >> result.txt") ;

    // Saveguard of the whole configuration
    // ------------------------------------

	FILE* fresu = fopen("resu.d", "w") ;

	star.get_mp().get_mg()->sauve(fresu) ;		// writing of the grid
	star.get_mp().sauve(fresu) ;                // writing of the mapping
	star.get_eos().sauve(fresu) ;  				// writing of the EOS
	star.sauve(fresu) ;                         // writing of the star
	
	fclose(fresu) ;
	
    // Drawings
    // --------

    ofstream fichdes("prof_n.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       N(theta=0)	N(theta=pi/2)" << endl ; 
    int npd = 300 ;
    double r_max = 3 ;
    double h = r_max/double(npd-1) ;
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << star.get_nn().val_point(rr,0.,0.) << "   " 
	  << star.get_nn().val_point(rr, M_PI/2.,0.) << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_omega.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       omega(theta=0)	  omega(theta=pi/2)" << endl ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << star.get_nphi().val_point(rr,0.,0.) << "   " 
	  << star.get_nphi().val_point(rr, M_PI/2.,0.) << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_omega.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       omega(theta=0)/Omega	  omega(theta=pi/2)/Omega" << endl ; 
    double omega_rot = star.get_omega_c() ; 
    if (omega_rot == double(0)) omega_rot = 1 ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << star.get_nphi().val_point(rr,0.,0.) / omega_rot << "   " 
	  << star.get_nphi().val_point(rr, M_PI/2.,0.) / omega_rot << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_a.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       A(theta=0)	  A(theta=pi/2)" << endl ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << sqrt( star.get_a_car().val_point(rr,0.,0.) ) << "   " 
	  << sqrt( star.get_a_car().val_point(rr, M_PI/2.,0.) ) << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_bma.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       B-A(theta=0)	  B-A(theta=pi/2)" << endl ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << sqrt( star.get_b_car().val_point(rr,0.,0.) ) 
	   - sqrt( star.get_a_car().val_point(rr,0.,0.) ) << "   " 
	  << sqrt( star.get_b_car().val_point(rr, M_PI/2.,0.) ) 
	   - sqrt( star.get_a_car().val_point(rr, M_PI/2.,0.) ) << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_ener.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       ener(theta=0)	  ener(theta=pi/2)" << endl ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << star.get_ener().val_point(rr,0.,0.) << "   " 
	  << star.get_ener().val_point(rr, M_PI/2.,0.) << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_h.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       H(theta=0)	  H(theta=pi/2)" << endl ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << star.get_ent().val_point(rr,0.,0.) << "   " 
	  << star.get_ent().val_point(rr, M_PI/2.,0.) << endl ; 
    }
    fichdes.close() ; 

    fichdes.open("prof_u.d") ;
    fichdes.precision(10) ; 
    fichdes << "#   r [km]       U(theta=0)/c	  U(theta=pi/2)/c" << endl ; 
    for (int i=0; i<npd; i++) {
	double rr = h * i ;
	fichdes << rr*10 << "  " 
	  << star.get_uuu().val_point(rr,0.,0.) << "   " 
	  << star.get_uuu().val_point(rr, M_PI/2.,0.) << endl ; 
    }
    fichdes.close() ; 

    
   if (graph == 1) {

	char title[80] ;
	char bslash[2] = {92, '\0'} ;  // 92 is the ASCII code for backslash 

	// Scalar defining the surface of the star (via the enthalpy field)
	Scalar surf( star.get_ent() ); 
	Scalar surf_ext(mp) ; 
	surf_ext = - 0.2 * surf.val_grid_point(0, 0, 0, 0) ; 
	surf_ext.annule(0, star.get_nzet()-1) ; 
	surf.annule(star.get_nzet(), mg.get_nzone()-1) ; 
	surf = surf + surf_ext ;
	surf = raccord_c1(surf, star.get_nzet()) ; 

	int nzdes = star.get_nzet() ; 

	des_coupe_y(star.get_ent(), 0., nzdes, "Log-enthalpy", &surf) ; 

	cout << endl << "Plot of the coefficients of the cos(j theta) expansion of the function G(theta) defining the stellar surface:" << endl ; 
	des_map_et(mp, star.get_nzet()-1) ; 

	des_coupe_y(star.get_ener(), 0., nzdes, "Fluid Proper energy density", &surf) ; 

	if (mer_triax < mer_max) { 
	    des_coupe_z(star.get_ent(), 0., nzdes, "Log-enthalpy (equatorial plane)", 
			&surf) ; 
	}
	    
	char partit[] = {92, 'g', 'n', '\0'} ; 
	strcpy(title, "Gravitational potential ") ; 
	strcat(title, partit) ; 

	des_coupe_y(star.get_logn(), 0., nzdes, title, &surf) ; 
	
	
	strcpy(title, "Azimuthal shift N") ; 
	strcat(title, bslash) ; 
	strcat(title, "u") ; 
	strcat(title, bslash) ; 
	strcat(title, "gf") ; 
	des_coupe_y(star.get_nphi(), 0., nzdes, title, &surf) ; 
	
	strcpy(title, "Metric potential ") ; 
	strcat(title, bslash) ; 
	strcat(title, "gz") ; 
	des_coupe_y(star.get_dzeta(), 0., nzdes, title, &surf) ; 
	
	strcpy(title, "Metric potential (NB-1) r sin") ; 
	strcat(title, bslash) ; 
	strcat(title, "gh") ; 
	des_coupe_y(star.get_tggg(), 0., nzdes, title, &surf) ; 
	
	char debtit[] = {'A', 92, 'u', '2', 92, 'd', ' ', 'K', 92, 'u', '\0'} ; 
	strcpy(title, debtit) ; 
	strcat(title, "ij") ; 
	strcat(title, bslash) ; 
	strcat(title, "d K") ; 
	strcat(title, bslash) ; 
	strcat(title, "dij") ; 
	strcat(title, bslash) ; 
	strcat(title, "u") ; 

	des_coupe_y(star.get_ak_car(), 0., nzdes, title, &surf) ; 

    }

    ofstream seq("seq.d", ios::app) ; 
    seq << star.get_ent().val_grid_point(0,0,0,0) << "   " << star.mass_g() / msol 
    << "    " << qpig/(4.*M_PI) * star.mass_g() / star.r_circ() << endl ; 
   
    seq.close() ; 
	
    
    // Output file for GYOTO
    //----------------------
    
    // Shift vector on spherical triad --> beta
    Scalar beta_phi = - star.get_nphi() ; // beta^phi = - N^phi
    Vector beta(star.get_mp(), CON, star.get_mp().get_bvect_spher()) ;
    beta.set(1) = 0 ;
    beta.set(2) = 0 ;
    Scalar tmp = beta_phi ; 
    tmp.mult_rsint() ; 
    beta.set(3) = tmp ; 
    
    // Extrinsic curvature --> kk
    Scalar nn = star.get_nn() ;
    Scalar bb2 = star.get_b_car() ; 
    Sym_tensor kk(star.get_mp(), COV, star.get_mp().get_bvect_spher()) ;
    kk.set(1,1) = 0 ; 
    kk.set(1,2) = 0 ; 
    tmp = 0.5 * bb2 * beta_phi.dsdr() / nn ;
    tmp.mult_rsint() ;
    kk.set(1,3) = tmp ;
    kk.set(2,2) = 0 ; 
    tmp = 0.5 * bb2 * beta_phi.dsdt() / nn ;
    tmp.mult_sint() ;
    tmp.inc_dzpuis(2); 
    kk.set(2,3) = tmp ;
    kk.set(3,3) = 0 ;

    // File for GYOTO
    //-----------
    Valeur ns_surf(mg.get_angu()) ;
    ns_surf.annule_hard() ;
    for (int k=0; k<np; k++)
      for (int j=0; j<nt; j++)  {
      ns_surf.set(nzet-1, k, j, 0) = mp.val_r_jk(star.l_surf()(k, j), star.xi_surf()(k, j), j, k) ;
      }

    // Horizon 
    //--------
    Valeur ah(mg.get_angu()) ;
    ah.annule_hard() ;

    Metric gamma = star.get_gamma() ;

    // Lorentz facor and 3-velocity
    //-----------------------------
    Vector velo = star.get_u_euler();
    velo.change_triad(mp.get_bvect_spher());
    Vector u_euler = velo.up_down(gamma);

    // 4-acceleration vector 
    //----------------------

    Scalar loggam = log(star.get_gam_euler()) ; 
    loggam.std_spectral_base() ; 

    Vector a = (star.get_logn()).derive_cov( mp.flat_met_spher() ) 
             - loggam.derive_cov( mp.flat_met_spher() ); 
    a.std_spectral_base() ; 

    for (int i=nzet; i<nz; i++) 
      a.annule_domain(i); 

    // Scalar g = sqrt(contract( a, 0, a.up_down(gamma), 0 )) ;

    // File for GYOTO
    FILE* file_out = fopen("resu_gyoto.d", "w") ;
    double total_time = 0. ; // for compatibility
    fwrite_be(&total_time, sizeof(double), 1, file_out) ;    
    star.get_mp().get_mg()->sauve(file_out) ;
    star.get_mp().sauve(file_out) ;
    star.get_nn().sauve(file_out) ;
    beta.sauve(file_out) ;
    star.get_gamma().cov().sauve(file_out) ;
    star.get_gamma().con().sauve(file_out) ;
    kk.sauve(file_out) ;

    // Lorentz factor and 3-velocity save 
    star.get_gam_euler().sauve(file_out) ;
    u_euler.sauve(file_out) ;

    // NS surface save 
    ns_surf.get_mg()->sauve(file_out) ;
    ns_surf.sauve(file_out) ;

    // 4-acceleration vector save 
    a.sauve(file_out) ;   
 
    // Horizon save 
    ah.get_mg()->sauve(file_out) ;
    ah.sauve(file_out) ;

    //  aa = c J / (G M^2)
    double aa = star.angu_mom()/( qpig / (4* M_PI) * pow(star.mass_g(), 2.) ) ;
    fwrite_be(&aa, sizeof(double), 1, file_out) ;
    fclose(file_out) ;    

    // Cleaning
    // --------

    delete peos ;    

    exit(EXIT_SUCCESS) ; 
    
    return EXIT_SUCCESS ; 
   
}
