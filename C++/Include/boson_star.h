/*
 *  Definition of Lorene class Boson_star
 *
 */

/*
 *   Copyright (c) 2012 Claire Some, Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

#ifndef __BOSON_STAR_H_ 
#define __BOSON_STAR_H_ 

/*
 * $Id: boson_star.h,v 1.4 2014/10/13 08:52:32 j_novak Exp $
 * $Log: boson_star.h,v $
 * Revision 1.4  2014/10/13 08:52:32  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2012/12/03 15:26:14  c_some
 * Added data member m2
 *
 * Revision 1.2  2012/11/23 15:42:14  c_some
 * Small changes
 *
 * Revision 1.1  2012/11/22 16:03:16  c_some
 * New class Boson_star
 *
 *
 * $Header: /cvsroot/Lorene/C++/Include/boson_star.h,v 1.4 2014/10/13 08:52:32 j_novak Exp $
 *
 */


// Headers Lorene
#include "compobj.h"


			//--------------------------//
			//    class Boson_star      //
			//--------------------------//

namespace Lorene {
/**
 * Class for stationary axisymmetric boson stars (***under development***). 
 * \ingroup(compactobjects)
 *
 */
class Boson_star : public Star_QI {

    // Data : 
    // -----
    protected:

	/** Real part of the scalar field Phi
	 */
	Scalar rphi ;

	/** Imaginary part of the scalar field Phi
	 */
	Scalar iphi ;

	/**  Coefficient omega in the time dependence of Phi
	 */
	 double omega ; 
	 
	/**  Coefficient kkk in the azimuthal dependence of Phi
	 */
	 int kkk ; 

	/**  Boson mass 
	 */
	 double mmm ; 
	 	 
	/**  Boson mass squared
	 */
	 double m2 ; 
	 	 

    // Derived data : 
    // ------------
    protected:
	
    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star is contructed
	 * @param m Boson mass  //## which unit ?
	 * @param k Coefficient \e k in the azimuthal dependence of Phi
     *
	 */
	Boson_star(Map& mp_i, double m, int k) ;			
	
	
	Boson_star(const Boson_star& ) ;		///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE*) ). 
	 * 
	 * @param mp_i Mapping on which the star is constructed
	 * @param fich	input file (must have been created by the function
	 *	\c Boson_star::sauve )
	 */
	Boson_star(Map& mp_i, FILE* fich) ;    		

	virtual ~Boson_star() ;			///< Destructor


    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0  all the pointers on derived quantities
	virtual void set_der_0x0() const ; 

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Boson_star 
	void operator=(const Boson_star& ) ;	
	
    // Accessors
    // ---------
    public:

	/** Returns the real part of the scalar field
	 */
	const Scalar& get_rphi() const {return rphi;} ;

	/** Returns the imaginary part of the scalar field
	 */
	const Scalar& get_iphi() const {return iphi;} ;

	/** Sets a value to the real part of the scalar field 
	 */
	Scalar& set_rphi() {del_deriv(); return rphi;} ;

	/** Sets a value to the imaginary part of the scalar field 
	 */
	Scalar& set_iphi() {del_deriv(); return iphi;} ;


    // Outputs
    // -------
    public:
	virtual void sauve(FILE* ) const ;	    ///< Save in a file
    
    protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    

    // Global quantities
    // -----------------
    public:
		
    
    // Computational routines
    // ----------------------
    public: 

	/** Computes the 3+1 components of the energy-momentum tensor (E, P_i and S_{ij})
	 *  from the values of the scalar field and the metric 
	 */
	void update_ener_mom() ; 
	
	/** Solves the equation satisfied by the scalar field
	 */
	// void solve_phi() ; 
	
	/** Computes an equilibrium configuration.
	 *  
	 *  @param rphi_c  [input] Central value of the real part of the scalar field
	 *  @param iphi_c  [input] Central value of the imaginary part of the scalar field
	 *  @param nzadapt  [input] Number of (inner) domains where the mapping 
	 *			    adaptation to an iso-enthalpy surface
	 *			    should be performed
	 *  @param ent_limit [input] 1-D \c Tbl  of dimension \c nzet  which
	 *				defines the enthalpy at the outer boundary
	 *				of each domain
	 *  @param icontrol [input] Set of integer parameters (stored as a
	 *			    1-D \c Itbl  of size 8) to control the 
	 *			    iteration: 
	 *	\li \c icontrol(0) = mer_max  : maximum number of steps 
	 *	\li \c icontrol(1) = mer_rot  : step at which the rotation is 
	 *				      switched on 
	 *	\li \c icontrol(2) = mer_change_omega  : step at which the rotation
	 *			  velocity is changed to reach the final one  
	 *	\li \c icontrol(3) = mer_fix_omega  :  step at which the final
	 *			    rotation velocity must have been reached  
	 *	\li \c icontrol(4) = mer_mass  : the absolute value of 
	 *			    \c mer_mass  is the step from which the 
	 *			    baryon mass is forced to converge, 
	 *			    by varying the central enthalpy 
	 *			    (\c mer_mass>0 ) or the angular 
	 *			    velocity (\c mer_mass<0 ) 
	 *	\li \c icontrol(5) = mermax_poisson  : maximum number of steps in 
	 *				\c Map_et::poisson  
	 *	\li \c icontrol(6) = mer_triax  : step at which the 3-D 
	 *				perturbation is switched on 
	 *	\li \c icontrol(7) = delta_mer_kep  : number of steps
	 *			    after \c mer_fix_omega  when \c omega 
	 *			    starts to be increased by \c fact_omega 
	 *			    to search for the Keplerian velocity
	 * 	 
	 *  @param control [input] Set of parameters (stored as a 
	 *			    1-D \c Tbl  of size 7) to control the 
	 *			    iteration: 
	 *	\li \c control(0) = precis  : threshold on the enthalpy relative 
	 *				change for ending the computation 
	 *	\li \c control(1) = omega_ini  : initial angular velocity, 
	 *			    switched on only if \c mer_rot<0 , 
	 *			    otherwise 0 is used  
	 *	\li \c control(2) = relax  : relaxation factor in the main 
	 *				   iteration  
	 *	\li \c control(3) = relax_poisson  : relaxation factor in 
	 *				   \c Map_et::poisson  
	 *	\li \c control(4) = thres_adapt  :  threshold on dH/dr for 
	 *			    freezing the adaptation of the mapping 
	 *	\li \c control(5) = ampli_triax  :  relative amplitude of 
	 *			    the 3-D perturbation 
	 *	\li \c control(6) = precis_adapt  : precision for 
	 *			    \c Map_et::adapt 
	 *
	 *  @param diff [output]   1-D \c Tbl  of size 7 for the storage of 
	 *			    some error indicators : 
	 *	    \li \c diff(0)  : Relative change in the enthalpy field
	 *			      between two successive steps 
	 *	    \li \c diff(1)  : Relative error in the resolution of the
	 *			    Poisson equation for \c nuf    
	 *	    \li \c diff(2)  : Relative error in the resolution of the
	 *			    Poisson equation for \c nuq    
	 *	    \li \c diff(3)  : Relative error in the resolution of the
	 *			    Poisson equation for \c dzeta    
	 *	    \li \c diff(4)  : Relative error in the resolution of the
	 *			    Poisson equation for \c tggg    
	 *	    \li \c diff(5)  : Relative error in the resolution of the
	 *			    equation for \c shift  (x comp.)   
	 *	    \li \c diff(6)  : Relative error in the resolution of the
	 *			    equation for \c shift  (y comp.)   
	 */
	virtual void equilibrium(double rphi_c, double iphi_c,
			 int nzadapt, const Tbl& ent_limit,
			 const Itbl& icontrol, const Tbl& control,
			 Tbl& diff, Param* = 0x0) ;
			
};


}
#endif
