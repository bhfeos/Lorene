/*
 *  Definition of Lorene class Etoile_bin_nsbh
 *
 */

/*
 *   Copyright (c) 2003 Keisuke Taniguchi
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

#ifndef __ET_BIN_NSBH_H_ 
#define __ET_BIN_NSBH_H_ 

/*
 * $Id: et_bin_nsbh.h,v 1.9 2014/10/13 08:52:34 j_novak Exp $
 * $Log: et_bin_nsbh.h,v $
 * Revision 1.9  2014/10/13 08:52:34  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2006/09/25 10:01:45  p_grandclement
 * Addition of N-dimensional Tbl
 *
 * Revision 1.7  2006/06/01 12:47:50  p_grandclement
 * update of the Bin_ns_bh project
 *
 * Revision 1.6  2006/04/25 07:21:54  p_grandclement
 * Various changes for the NS_BH project
 *
 * Revision 1.5  2005/10/18 13:12:31  p_grandclement
 * update of the mixted binary codes
 *
 * Revision 1.4  2005/08/29 15:10:12  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.3  2004/06/09 06:32:51  k_taniguchi
 * Introduce set_n_auto() and set_confpsi_auto().
 *
 * Revision 1.2  2003/10/24 12:25:19  k_taniguchi
 * Introduce the method of update metric for NS-BH
 *
 * Revision 1.1  2003/10/21 11:46:13  k_taniguchi
 * Definition of class Et_bin_nsbh
 *
 *
 * $Header: /cvsroot/Lorene/C++/Include/et_bin_nsbh.h,v 1.9 2014/10/13 08:52:34 j_novak Exp $
 *
 */

// Lorene headers
#include "etoile.h"
#include "bhole.h"

namespace Lorene {
/**
 * Class for a star in a NS-BH binary system
 * 
 * This class is a derived class from {\tt Etoile\_bin}
 *
 * @version #$Id: et_bin_nsbh.h,v 1.9 2014/10/13 08:52:34 j_novak Exp $#
 */
class Et_bin_nsbh : public Etoile_bin {

    // Data : 
    // -----
    protected:
    /// Part of the lapse {\it N} generated principaly by the star.
    Tenseur n_auto ;

    /// Part of the lapse {\it N} generated principaly by the companion star. 
    Tenseur n_comp ;

    /** Gradient of {\tt n\_auto}
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    Tenseur d_n_auto ;

    /** Gradient of {\tt n\_comp}
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    Tenseur d_n_comp ;

    /// Total conformal factor $\Psi$
    Tenseur confpsi ;

    /// Part of the conformal factor $\Psi$ generated principaly by the star.
    Tenseur confpsi_auto ;

    /** Part of the conformal factor $\Psi$ generated principaly
     *  by the companion star.
     */
    Tenseur confpsi_comp ;

    /** Gradient of {\tt confpsi\_auto}
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    Tenseur d_confpsi_auto ;

    /** Gradient of {\tt confpsi\_comp}
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    Tenseur d_confpsi_comp ;

    /** Part of the extrinsic curvature tensor 
     *  $\tilde A^{ij} = 2 N K^{ij}$
     *  generated by {\tt shift\_auto}. 
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    Tenseur_sym taij_auto ;

    /** Part of the extrinsic curvature tensor 
     *  $\tilde A^{ij} = 2 N K^{ij}$
     *  generated by {\tt shift\_comp}. 
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    Tenseur_sym taij_comp ;

    /** Total extrinsic curvature tensor 
     *  $\tilde A^{ij} = 2 N K^{ij}$
     *  generated by {\tt shift\_auto} and {\tt shift\_comp}. 
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    Tenseur_sym taij_tot ;

    /** Part of the extrinsic curvature tensor $K^{ij}$
     *  generated by {\tt shift\_auto}. 
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    Tenseur_sym tkij_auto ;

    /** Total extrinsic curvature tensor $K^{ij}$
     *  generated by {\tt shift\_auto} and {\tt shift\_comp}. 
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    Tenseur_sym tkij_tot ;

    /** Effective source at the previous step for the resolution of 
     *  the Poisson equation for {\tt n\_auto} by means of
     *  {\tt Map\_et::poisson}.
     */
    Cmp ssjm1_lapse ;

    /** Effective source at the previous step for the resolution of 
     *  the Poisson equation for {\tt confpsi\_auto} by means of
     *  {\tt Map\_et::poisson}.
     */
    Cmp ssjm1_confpsi ;

    // Constructors - Destructor
    // -------------------------
    public:
    /** Standard constructor.
     *
     * @param mp_i Mapping on which the star will be defined
     * @param nzet_i Number of domains occupied by the star
     * @param relat should be {\tt true} for a relativistic star,
     *                        {\tt false} for a Newtonian one
     * @param eos_i Equation of state of the stellar matter
     * @param irrot should be {\tt true} for an irrotational star,
     *                        {\tt false} for a corotating one
     * @param ref_triad_i  Reference triad ("absolute frame"),
     *        with respect to which the components of all the member
     *        {\tt Tenseur}'s are defined, except for {\tt w\_shift}
     *        and {\tt ssjm1\_wshift} whose components are defined
     *        with respect to the mapping {\tt mp} Cartesian triad.
     */
    Et_bin_nsbh(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i,
		bool irrot, const Base_vect& ref_triad_i) ;

    Et_bin_nsbh(const Et_bin_nsbh& ) ;		/// Copy constructor

    /** Constructor from a file (see {\tt sauve(FILE* )})
     *
     * @param mp_i Mapping on which the star will be defined
     * @param eos_i Equation of state of the stellar matter
     * @param ref_triad_i  Reference triad ("absolute frame"),
     *        with respect to which the components of all the member
     *        {\tt Tenseur}'s are defined, except for {\tt w\_shift}
     *        and {\tt ssjm1\_wshift} whose components are defined
     *        with respect to the mapping {\tt mp} Cartesian triad.
     * @param fich  input file (must have been created by the function
     *        {\tt sauve})
     */
    Et_bin_nsbh(Map& mp_i, const Eos& eos_i, const Base_vect& ref_triad_i,
		FILE* fich, bool old = false) ;    		

    virtual ~Et_bin_nsbh() ;			/// Destructor
 

    // Mutators / assignment
    // ---------------------
    public:
    /// Assignment to another {\tt Et\_bin\_nsbh}
    void operator=(const Et_bin_nsbh&) ;	
	
    /** Read/write the lapse {\it N} generated principaly
     *  by the star.
     */
    Tenseur& set_n_auto() ;

    /** Read/write the lapse {\it N} generated principaly
     *  by the companion star.
     */
    Tenseur& set_n_comp() ;

    /** Read/write the conformal factor $\Psi$ generated principaly
     *  by the star.
     */
    Tenseur& set_confpsi_auto() ;

    /** Read/write the conformal factor $\Psi$ generated principaly
     *  by the companion star.
     */
    Tenseur& set_confpsi_comp() ;

    // Accessors
    // ---------
    public:
    /** Returns the part of the lapse {\it N} generated principaly
     *  by the star. 
     */
    const Tenseur& get_n_auto() const {return n_auto;} ;

    /** Returns the part of the lapse {\it N} generated principaly
     *  by the companion star. 
     */
    const Tenseur& get_n_comp() const {return n_comp;} ;

    /** Returns the gradient of {\tt n\_auto}
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    const Tenseur& get_d_n_auto() const {return d_n_auto;} ;

    /** Returns the gradient of {\tt n\_comp}
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    const Tenseur& get_d_n_comp() const {return d_n_comp;} ;

    /// Returns the part of the conformal factor $\Psi$
    const Tenseur& get_confpsi() const {return confpsi;} ;

    /** Returns the part of the conformal factor $\Psi$
     *  generated principaly by the star. 
     */
    const Tenseur& get_confpsi_auto() const {return confpsi_auto;} ;

    /** Returns the part of the conformal factor $\Psi$
     *  generated principaly by the companion star. 
     */
    const Tenseur& get_confpsi_comp() const {return confpsi_comp;} ;

    /** Returns the gradient of {\tt confpsi\_auto}
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    const Tenseur& get_d_confpsi_auto() const {return d_confpsi_auto;} ;

    /** Returns the gradient of {\tt confpsi\_comp} 
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    const Tenseur& get_d_confpsi_comp() const {return d_confpsi_comp;} ;

    /** Returns the part of the extrinsic curvature tensor 
     *  $\tilde A^{ij} = 2 N K^{ij}$ generated by {\tt shift\_auto}. 
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    const Tenseur_sym& get_taij_auto() const {return taij_auto;} ;

    /** Returns the part of the extrinsic curvature tensor 
     *  $\tilde A^{ij} = 2 N K^{ij}$ generated by {\tt shift\_comp}. 
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    const Tenseur_sym& get_taij_comp() const {return tkij_comp;} ;

    /** Returns the total extrinsic curvature tensor 
     *  $\tilde A^{ij} = 2 N K^{ij}$
     *  generated by {\tt shift\_auto} and {\tt shift\_comp}. 
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    const Tenseur_sym& get_taij_tot() const {return taij_tot;} ;

    /** Returns the part of the extrinsic curvature tensor $K^{ij}$
     *  generated by {\tt shift\_auto}. 
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    const Tenseur_sym& get_tkij_auto() const {return tkij_auto;} ;

    /** Returns the total extrinsic curvature tensor $K^{ij}$
     *  generated by {\tt shift\_auto} and {\tt shift\_comp}. 
     *  (Cartesian components with respect to {\tt ref\_triad})
     */
    const Tenseur_sym& get_tkij_tot() const {return tkij_tot;} ;


    // Outputs
    // -------
    public:
    virtual void sauve(FILE *) const ;	    /// Save in a file

    protected:
    /// Operator >> (virtual function called by the operator <<).
    virtual ostream& operator>>(ostream& ) const ;


    // Global quantities
    // -----------------

    // Computational routines
    // ----------------------
    public:
        
    void fait_taij_auto() ;
    /** Computes metric coefficients from known potentials,
     *  when the companion is a black hole.
     *
     *  @param comp companion black hole
     *
     */
    void update_metric(const Bhole& comp) ;

    /** Computes the derivative of metric functions related to the
     *  companion black hole.
     *
     *  @param comp companion BH.
     *
     */
    void update_metric_der_comp(const Bhole& comp) ;

    /** Computes an equilibrium configuration in a NS-BH binary system.
     * 
     *  The values of {\tt n\_comp}, {\tt confpsi\_comp}, {\tt pot\_centri}
     *  are held fixed during the iteration. 
     *  
     *  @param ent_c  [input] Central enthalpy
     *  @param mermax [input] Maximum number of steps 
     *  @param mermax_poisson [input]   Maximum number of steps in 
     *				    Map\_et::poisson
     *  @param relax_poisson [input]  Relaxation factor in Map\_et::poisson
     *  @param mermax_potvit [input]  Maximum number of steps in 
     *				  Map\_radial::poisson\_compact
     *  @param relax_potvit [input]   Relaxation factor in 
     *				  Map\_radial::poisson\_compact
     *  @param thres_adapt  [input]   Threshold on dH/dr for the adaptation 
     *				  of the mapping
     *  @param fact [input]    1-D {\tt Tbl} for the input
     *                          of some factors : \\
     *          {\tt fact(0)} : A resizing factor for the first shell
     *  @param diff [output]   1-D {\tt Tbl} for the storage of some
     *			    error indicators : \\
     *	    {\tt diff(0)} : Relative change in the enthalpy field
     *			      between two successive steps \\
     *	    {\tt diff(1)} : Relative error returned by the routine
     *				{\tt Etoile\_bin::velocity\_potential}  
     *	    {\tt diff(2)} : Relative error in the resolution of the
     *			    Poisson equation for {\tt n\_auto} \\  
     *	    {\tt diff(3)} : Relative error in the resolution of the
     *			    Poisson equation for {\tt confpsi\_auto} \\  
     *	    {\tt diff(4)} : Relative error in the resolution of the
     *			    equation for {\tt shift\_auto} (x comp.) \\  
     *	    {\tt diff(5)} : Relative error in the resolution of the
     *			    equation for {\tt shift\_auto} (y comp.) \\  
     *	    {\tt diff(6)} : Relative error in the resolution of the
     *			    equation for {\tt shift\_auto} (z comp.)   
     */
    virtual void equilibrium_nsbh(double ent_c, int mermax,
				  int mermax_poisson, double relax_poisson,
				  int mermax_potvit, double relax_potvit,
				  double thres_adapt,
				  const Tbl& fact, Tbl& diff) ;

    void equilibrium_nsbh (bool, double, int&, int, int, double, int, double, Tbl&) ;
    
    /** Computes the quantities \c bsn  and \c pot_centri .
	 * 
	 *  The calculation is performed starting from the quantities
	 *  \c nnn , \c shift ,  \c a_car ,  
	 *  which are supposed to be up to date.  
	 * 
	 *  @param omega  angular velocity with respect to an asymptotically 
	 *		  inertial observer
	 *  @param x_axe  absolute X coordinate of the rotation axis
	 */
	virtual void kinematics(double omega, double x_axe) ; 
	
	double compute_angul() const ;
	double compute_axe(double) const ;
	
    friend class Bin_ns_bh ;

};

}
#endif