/*
 *  Definition of Lorene class Et_rot_bifluid
 *
 */

/*
 *   Copyright (c) 2001 Jerome Novak
 *             (c) 2015 Aurelien Sourie
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


#ifndef __ET_ROT_BIFLUID_H_ 
#define __ET_ROT_BIFLUID_H_ 

/*
 * $Id: et_rot_bifluid.h,v 1.21 2017/10/06 12:36:33 a_sourie Exp $
 * $Log: et_rot_bifluid.h,v $
 * Revision 1.21  2017/10/06 12:36:33  a_sourie
 * Cleaning of tabulated 2-fluid EoS class + superfluid rotating star model.
 *
 * Revision 1.20  2015/06/26 14:10:08  j_novak
 * Modified comments.
 *
 * Revision 1.19  2015/06/11 13:50:18  j_novak
 * Minor corrections
 *
 * Revision 1.18  2015/06/10 14:39:17  a_sourie
 * New class Eos_bf_tabul for tabulated 2-fluid EoSs and associated functions for the computation of rotating stars with such EoSs.
 *
 * Revision 1.17  2014/10/13 08:52:34  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.16  2013/11/25 13:50:55  j_novak
 * The inheritance from Etoile_rot is no longer virtual.
 *
 * Revision 1.15  2011/10/06 14:55:36  j_novak
 * equation_of_state() is now virtual to be able to call to the magnetized
 * Eos_mag.
 *
 * Revision 1.14  2004/09/01 10:56:05  r_prix
 * added option of converging baryon-mass to equilibrium_bi()
 *
 * Revision 1.13  2004/03/22 13:12:41  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.12  2003/12/04 14:13:32  r_prix
 * added method get_typeos {return typeos}; and fixed some comments.
 *
 * Revision 1.11  2003/11/20 14:01:45  r_prix
 * changed member names to better conform to Lorene coding standards:
 * J_euler -> j_euler, EpS_euler -> enerps_euler, Delta_car -> delta_car
 *
 * Revision 1.10  2003/11/18 18:32:36  r_prix
 * added new class-member: EpS_euler := ener_euler + s_euler
 * has the advantage of a nice Newtonian limit -> rho
 * (ener_euler is no longer used in this class!)
 *
 * Revision 1.9  2003/11/13 12:02:03  r_prix
 * - adapted/extended some of the documentation
 * - changed xxx2 -> Delta_car
 * - added members J_euler, sphph_euler, representing 3+1 components of Tmunu
 *   (NOTE: these are not 2-fluid specific, and should ideally be moved into Class Etoile!)
 *
 * Revision 1.8  2003/09/17 08:27:50  j_novak
 * New methods: mass_b1() and mass_b2().
 *
 * Revision 1.7  2002/10/09 07:54:29  j_novak
 * Et_rot_bifluid and Et_rot_mag inheritate virtually from Etoile_rot
 *
 * Revision 1.6  2002/09/13 09:17:33  j_novak
 * Modif. commentaires
 *
 * Revision 1.5  2002/04/05 09:09:36  j_novak
 * The inversion of the EOS for 2-fluids polytrope has been modified.
 * Some errors in the determination of the surface were corrected.
 *
 * Revision 1.4  2002/01/16 15:03:28  j_novak
 * *** empty log message ***
 *
 * Revision 1.3  2002/01/08 14:43:53  j_novak
 * better determination of surfaces for 2-fluid stars
 *
 * Revision 1.2  2002/01/03 15:30:27  j_novak
 * Some comments modified.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2001/10/03  09:49:06  novak
 * *** empty log message ***
 *
 * Revision 1.2  2001/08/28 14:14:10  novak
 * overrided l_surf function
 *
 * Revision 1.1  2001/06/22 15:38:52  novak
 * Initial revision
 *
 *
 * $Header: /cvsroot/Lorene/C++/Include/et_rot_bifluid.h,v 1.21 2017/10/06 12:36:33 a_sourie Exp $
 *
 */

// Headers Lorene
#include "eos_bifluid.h"
#include "etoile.h"

namespace Lorene {

  // Local prototype (for determining the surface)
  Cmp prolonge_c1(const Cmp& uu, const int nzet) ;
  
  /**
   * Class for two-fluid rotating relativistic stars. \ingroup (star)
   * 
   * This is a child class of \c Etoile_rot , with the same metric
   * and overloaded member functions. 
   * 
   * There are two number-density fields \c nbar  and \c nbar2  
   * (and 2 log-enthalpies, see class \c Eos_bifluid ), as well as two 
   * velocity fields, with phi-components (with respect to the Eulerian observer)
   * \c uuu  and \c uuu2 .
   *
   * Fluid 1 can be considered to correspond to the (superfluid) neutrons, whereas 
   * fluid 2 would consist of the protons (and electrons)
   *.
   * The quantity \c u_euler  of the \c class Etoile  is 
   * \b not \b used  in this class!
   * Only the "3+1" components of \f${T^\mu}_\nu\f$ should be used outside
   * of \c hydro_euler() , namely \c s_euler, \c sphph_euler, \c j_euler and 
   * \c ener_euler.
   *
   */
  class Et_rot_bifluid : public Etoile_rot {
    
    // Data : 
    // -----
  protected:
    const Eos_bifluid& eos ; ///< Equation of state for two-fluids model
    
    double omega2 ; ///< Rotation angular velocity for fluid 2 (\c [f_unit] )
    
    // Fluid quantities with respect to the fluid frame
    // ------------------------------------------------
    
    /// Log-enthalpy for the second fluid
    Tenseur ent2 ;
    
    Tenseur nbar2 ; ///< Baryon density in the fluid frame, for fluid 2
    
    Tenseur K_nn ;  ///< Coefficient \f$K_{nn}\f$
    Tenseur K_np ;  ///< Coefficient \f$K_{np}\f$
    Tenseur K_pp ;  ///< Coefficient \f$K_{pp}\f$
    Tenseur alpha_eos ;  ///< Coefficient \f$\alpha\f$ relative to entrainment effects

    // Fluid quantities with respect to the Eulerian frame
    // ---------------------------------------------------
    
    // FIXME: the following three variables are not specific to 2-fluid stars
    //  and should ideally be moved to class Etoile!
    
    /// The component \f$S^\varphi_\varphi\f$ of the stress tensor \f${S^i}_j\f$.
    Tenseur sphph_euler;
    
    /** Total angular momentum (flat-space!) 3-vector \f$J_\mathrm{euler}\f$,
     *  which is related to \f$J^i\f$ of the "3+1" decomposition, but expressed 
     *  in a flat-space triad. In axisymmetric circular cases, only 
     * \f$J_\mathrm{euler}(\varphi)=r \sin\theta\, J^\varphi\f$ is nonzero.
     */
    Tenseur j_euler;
    
    Tenseur j_euler1; ///< To compute \f$J_n\f$
    Tenseur j_euler2; ///< To compute \f$J_p\f$
    
    /// the combination \f$E+S_i^i\f$: useful because in the Newtonian limit \f$\rightarrow \rho\f$.
    Tenseur enerps_euler;
    
    /// Norm of the (fluid no.2) 3-velocity with respect to the eulerian observer
    Tenseur uuu2 ;
    
    /// Lorentz factor between the fluid 2 and Eulerian observers   
    Tenseur gam_euler2 ;
    
    /**
     * The "relative velocity" (squared) \f$\Delta^2\f$ of the two fluids.
     * See Prix et al.(2003) and see also \c Eos_bifluid . 
     */
    Tenseur delta_car ; 
    
    // Derived data : 
    // ------------
  protected:
    /// Coordinate radius at \f$\phi=0\f$, \f$\theta=\pi/2\f$. 
    mutable double* p_ray_eq2 ; 
    
    /// Coordinate radius at \f$\phi=\pi/2\f$, \f$\theta=\pi/2\f$. 
    mutable double* p_ray_eq2_pis2 ;
    
    /// Coordinate radius at \f$\phi=\pi\f$, \f$\theta=\pi/2\f$. 
    mutable double* p_ray_eq2_pi ;
    
    /// Coordinate radius at \f$\theta=0\f$. 
    mutable double* p_ray_pole2 ;
    
    /** Description of the surface of fluid 2: 2-D \c Itbl  containing the 
     *	values of the domain index \e l  on the surface at the 
     *	collocation points in \f$(\theta', \phi')\f$
     */
    mutable Itbl* p_l_surf2 ; 
    
    /** Description of the surface of fluid 2: 2-D \c Tbl  containing the 
     *	values of the radial coordinate \f$\xi\f$ on the surface at the 
     *	collocation points in \f$(\theta', \phi')\f$
     */
    mutable Tbl* p_xi_surf2 ; 
    
    mutable double* p_r_circ2 ;	///< Circumferential radius of fluid no.2
    mutable double* p_area2 ;	///< Surface area of fluid no.2
    mutable double* p_aplat2 ;	///< Flatening r_pole/r_eq of fluid no.2
    
    mutable double* p_mass_b1 ;	///< Baryon mass of fluid 1
    mutable double* p_mass_b2 ;	///< Baryon mass of fluid 2
    
    mutable double* p_angu_mom_1;	///< Angular momentum of fluid 1
    mutable double* p_angu_mom_2;	///< Angular momentum of fluid 2
   
    /// Quantities used to describe the different couplings between the fluids
    mutable double* p_coupling_mominert_1;	///< \f$\tilde{I}_n\f$ 
    mutable double* p_coupling_mominert_2;	///< \f$\tilde{I}_p\f$    
    mutable double* p_coupling_entr;		///< \f$\tilde{\varepsilon}_n\tilde{I}_n = \tilde{\varepsilon}_p\tilde{I}_p\f$ 
    mutable double* p_coupling_LT_1;		///< \f$\tilde{\omega}_n\tilde{I}_n\f$ 
    mutable double* p_coupling_LT_2;		///< \f$\tilde{\omega}_p\tilde{I}_p\f$ 

    
    // Constructors - Destructor
    // -------------------------
  public:
    
    Et_rot_bifluid(Map& mp_i, int nzet_i, bool relat, 
		   const Eos_bifluid& eos_i) ;     ///< Standard constructor
    
    Et_rot_bifluid(const Et_rot_bifluid& ) ;       ///< Copy constructor
    
    /** Constructor from a file (see \c sauve(FILE*) ) Works only for 
     *  relativistic stars.
     *  This has to be improved....
     */ 
    Et_rot_bifluid(Map& mp_i, const Eos_bifluid& eos_i, FILE* fich) ;    
    
    virtual ~Et_rot_bifluid() ;			///< Destructor
    
    // Memory management
    // -----------------
  protected:
    
    /// Deletes all the derived quantities
    virtual void del_deriv() const ; 
    
    /// Sets to \c 0x0  all the pointers on derived quantities
    virtual void set_der_0x0() const ; 
    
    /** Sets to \c ETATNONDEF  (undefined state) the hydrodynamical 
     *  quantities relative to the Eulerian observer.
     */
    virtual void del_hydro_euler() ; 
    
    // Mutators / assignment
    // ---------------------
  public:
    
    /// Assignment to another Et_rot_bifluid
    void operator=(const Et_rot_bifluid&) ;	
    
    /// Sets both enthalpy profiles
    void set_enthalpies(const Cmp&, const Cmp&) ;
    
    /** Computes a spherical static configuration. 
     * 
     *  @param ent_c [input] central value of the enthalpy 1
     *  @param ent_c2 [input] central value of the enthalpy 2
     *  @param precis [input] threshold in the relative difference between 
     *	the enthalpy fields of two consecutive steps
     *	to stop the iterative procedure (default value: 1.e-14)
     */
    void equilibrium_spher_bi(double ent_c, double ent_c2, 
			      double precis = 1.e-14) ;
    
    /** Computes a spherical static configuration. 
     *  The sources for Poisson equations are regularized
     *  by extracting analytical diverging parts.
     * 
     *  @param ent_c [input] central value of the enthalpy 1
     *  @param ent_c2 [input] central value of the enthalpy 2
     *  @param precis [input] threshold in the relative difference between 
     *	the enthalpy fields of two consecutive steps
     *	to stop the iterative procedure (default value: 1.e-14)
     */
    void equil_spher_regular(double ent_c, double ent_c2, 
			     double precis = 1.e-14) ;
    
    // Accessors
    // ---------
  public:
    
    /// Returns the equation of state
    const Eos_bifluid& get_eos() const {return eos; } ;
    
    /// Returns the enthalpy field for fluid 2 
    const Tenseur& get_ent2() const {return ent2 ; } ;
    
    /// Returns the proper baryon density for fluid 2
    const Tenseur& get_nbar2() const {return nbar2 ; } ;
    
    /// Returns the coefficient \f$K_{nn}\f$
    const Tenseur& get_K_nn() const {return K_nn ; } ;
    
    /// Returns the coefficient \f$K_{np}\f$
    const Tenseur& get_K_np() const {return K_np ; } ;
    
    /// Returns the coefficient \f$K_{pp}\f$
    const Tenseur& get_K_pp() const {return K_pp ; } ;
    
    /// Returns the coefficient \f$\alpha\f$ (entrainment parameter)
    const Tenseur& get_alpha_eos() const {return alpha_eos ; } ;    

    /// Returns the "relative velocity" (squared) \f$\Delta^2\f$ of the two fluids
    const Tenseur& get_delta_car() const {return delta_car ; } ;
    
    /// Returns the Lorentz factor between the fluid 2 and Eulerian observers
    const Tenseur& get_gam_euler2() const {return gam_euler2 ; } ;
    
    /// Returns the rotation angular velocity of fluid 2(\c [f_unit] ) 
    double get_omega2() const {return omega2 ; } ;

    /// Returns the norm of the fluid 2 3-velocity with respect to the eulerian frame
    const Tenseur& get_uuu2() const {return uuu2 ; } ;
    
    // Outputs
    // -------
  public:
    virtual void sauve(FILE *) const ;	    ///< Save in a file
    
    /// Operator >> (virtual function called by the operator <<). 
    virtual ostream& operator>>(ostream& ) const ;    
    
    /// Printing of some informations, excluding all global quantities
    virtual void partial_display(ostream& ) const ;    
    
    // Global quantities
    // -----------------
  public:
    
    /** Description of the surface of fluid 1: returns a 2-D \c Itbl  
     *	containing the 
     *	values of the domain index \e l  on the surface at the 
     *	collocation points in \f$(\theta', \phi')\f$.
     *	This surface is defined as the location where
     *	the density 1 (member \c nbar ) vanishes.
     */
    virtual const Itbl& l_surf() const ; 
    
    /** Description of the surface of fluid 2: returns a 2-D \c Itbl  
     *	containing the 
     *	values of the domain index \e l  on the surface at the 
     *	collocation points in \f$(\theta', \phi')\f$.
     *	This surface is defined as the location where
     *	the density 2 (member \c nbar2 ) vanishes.
     */
    const Itbl& l_surf2() const ; 
    
    /** Description of the surface of fluid 2: returns a 2-D \c Tbl  
     *	containing the values of the radial coordinate \f$\xi\f$ 
     *	on the surface at the 
     *	collocation points in \f$(\theta', \phi')\f$. 
     *	This surface is defined as the location where
     *	the density 2 (member \c nbar2 ) vanishes.
     */
    const Tbl& xi_surf2() const ; 

    /// Coordinate radius for fluid 2 at \f$\phi=0\f$, \f$\theta=\pi/2\f$ [r_unit].
    double ray_eq2() const ; 
    
    /// Coordinate radius for fluid 2 at \f$\phi=\pi/2\f$, \f$\theta=\pi/2\f$ [r_unit].
    double ray_eq2_pis2() const ; 
    
    /// Coordinate radius for fluid 2 at \f$\phi=\pi\f$, \f$\theta=\pi/2\f$ [r_unit].
    double ray_eq2_pi() const ; 
    
    /// Coordinate radius for fluid 2 at \f$\theta=0\f$ [r_unit]. 
    double ray_pole2() const ; 
    
    /// Baryon mass of fluid 1
    double mass_b1() const ;
    
    /// Baryon mass of fluid 2
    double mass_b2() const ;
    
    virtual double mass_b() const ;	///< Total Baryon mass
    virtual double mass_g() const ;	///< Gravitational mass
    virtual double angu_mom() const ;	///< Angular momentum
    
    /** Error on the virial identity GRV2.
     * Given by the integral Eq. (4.6) in 
     * [Bonazzola, Gougoulhon, Salgado, Marck, A\&A \b 278 , 421 (1993)].
     */
    virtual double grv2() const ;		
    
    /** Error on the virial identity GRV3.
     *  The error is computed as the integral defined
     *  by Eq. (43) of [Gourgoulhon and Bonazzola, 
     *  \a Class. \a Quantum \a Grav. \b 11 , 443 (1994)] divided by
     *  the integral of the matter terms.
     * 
     *  @param ost output stream to give details of the computation;
     *		if set to 0x0 [default value], no details will be
     *		given.
     *   
     */
    virtual double grv3(ostream* ost = 0x0) const ;	
    
    virtual double r_circ2() const ;       ///< Circumferential radius for fluid 2
    virtual double area2() const ;         ///< Surface area for fluid 2
    virtual double mean_radius2() const ;  ///< Mean radius for fluid 2
    virtual double aplat2() const ;        ///< Flatening r_pole/r_eq for fluid 2

    /** Quadrupole moment.
     *  The quadrupole moment \e Q is defined according to Eq. (11) of
     *  [Pappas and Apostolatos, \a Physical \a Review \a Letters 
     *  \b 108, 231104 (2012)]. This is a corrected version of the quadrupole
     *  moment defined by [Salgado, Bonazzola, Gourgoulhon and Haensel,
     *  \a Astron. \a Astrophys. \b 291 , 155 (1994)]. 
     *  Following this definition, \f$Q = {\bar Q } - 4/3 (1/4 + b) M^3 \f$, 
     *  where \f$ {\bar Q } \f$ is defined as the negative of the (wrong) quadrupole 
     *  moment defined in Eq. (7) of [Salgado, Bonazzola, Gourgoulhon and Haensel,
     *  \a Astron. \a Astrophys. \b 291 , 155 (1994)], \e b is defined by
     *  Eq. (3.37) of [Friedman and Stergioulas, \a Rotating \a Relativistic 
     *  \a Stars, Cambridge Monograph on mathematical physics] and \e M is 
     *  the gravitational mass of the star.
     */
    virtual double mom_quad() const ;	
    
    /** Part of the quadrupole moment.
     *  This term \f$ {\bar Q } \f$ is defined by Laarakkers and Poisson,
     *  \a Astrophys. \a J. \b 512 , 282 (1999).
     *  Note that \f$ {\bar Q }\f$ is the negative of the (wrong) quadrupole moment 
     *  defined in Eq. (7) of [Salgado, Bonazzola, Gourgoulhon and Haensel,
     *  \a Astron. \a Astrophys. \b 291 , 155 (1994)]. 
     */
    virtual double mom_quad_old() const ;
    
    /** Part of the quadrupole moment.
     *  \e B_o is defined as \f$bM^2\f$, where \e b is given by Eq. (3.37) of 
     *  [Friedman and Stergioulas, \a Rotating \a Relativistic \a Stars, 
     *  Cambridge Monograph on mathematical physics] and \e M is the 
     *  the gravitational mass of the star. 
     */
    virtual double mom_quad_Bo() const ;
    
    virtual double angu_mom_1() const ;	///< Angular momentum of fluid 1
    virtual double angu_mom_2() const ;	///< Angular momentum of fluid 2
    
    /** Quantities used to study the different fluid couplings:
     *  \f$\tilde{I}_X\f$, \f$\tilde{\varepsilon}_X\tilde{I}_X\f$ 
     *  and \f$\tilde{\omega}_X\tilde{I}_X\f$ .
	  *  These terms are defined in Eqs. (C5) - (C7) of Sourie et al., 
     *  \a MNRAS \a \b 464 , 4641-4657 (2017). 
     *  Note that these quantities only make sense in the 
     *  slow-rotation approximation and to first order in the lag.
     */    
    virtual double coupling_mominert_1() const ;
    virtual double coupling_mominert_2() const ;
    virtual double coupling_entr() const ;
    virtual double coupling_LT_1() const ;
    virtual double coupling_LT_2() const ;

    // Computational routines
    // ----------------------
  public: 
    /** Computes the hydrodynamical quantities relative to the Eulerian
     *  observer from those in the fluid frame.
     *
     *  The calculation is performed starting from the quantities
     *  \c ent , \c ent2 , \c ener , \c press , \c K_nn , \c K_np , \c K_pp
     *  and \c a_car, which are supposed to be up to date. From these,  
     *  the following fields are updated:
     *  \c delta_car , \c gam_euler , \c gam_euler2 , \c ener_euler , 
     *  \c s_euler , \c sphph_euler  and \c j_euler .
     */
    virtual void hydro_euler() ; 
    
    /** Computes the proper baryon and energy densities, as well as
     *  pressure and the coefficients Knn, Knp and Kpp, from the enthalpies 
     *  and both velocities.
     */
    virtual void equation_of_state() ; 
    
    /** Computes an equilibrium configuration.
     *  
     *  @param ent_c  [input] Central enthalpy for fluid 1 
     *  @param ent_c2 [input] Central enthalpy for fluid 2
     *  @param omega0  [input] Requested angular velocity for fluid 1
     *  @param omega20  [input] Requested angular velocity for fluid 2
     *  @param ent_limit [input] 1-D \c Tbl  of dimension \c nzet  which
     *				defines the enthalpy for fluid 1 at the 
     *                            outer boundary of each domain
     *  @param ent2_limit [input] 1-D \c Tbl  of dimension \c nzet  which
     *				defines the enthalpy for fluid 2 at the 
     *                            outer boundary of each domain
     *  @param icontrol [input] Set of integer parameters (stored as a
     *			    1-D \c Itbl  of size 5) to control the 
     *			    iteration: 
     *	\li \c icontrol(0) = mer_max  : maximum number of steps 
     *	\li \c icontrol(1) = mer_rot  : step at which the rotation is 
     *				      switched on 
     *	\li \c icontrol(2) = mer_change_omega  : step at which the rotation
     *			  velocity is changed to reach the final one  
     *	\li \c icontrol(3) = mer_fix_omega  :  step at which the final
     *			    rotation velocity must have been reached  
     *	\li \c icontrol(4) = mermax_poisson  : maximum number of steps in 
     *				\c Map_et::poisson  
     *  @param control [input] Set of parameters (stored as a 
     *			    1-D \c Tbl  of size 5) to control the 
     *			    iteration: 
     *	\li \c control(0) = precis  : threshold on the enthalpy relative 
     *				change for ending the computation 
     *	\li \c control(1) = omega_ini  : initial angular velocity, 
     *			    switched on only if \c mer_rot < 0 , 
     *			    otherwise 0 is used  
     *	\li \c control(2) = omega2_ini  : initial angular velocity, 
     *			    switched on only if \c mer_rot < 0 , 
     *			    otherwise 0 is used  
     *	\li \c control(3) = relax  : relaxation factor in the main 
     *				   iteration  
     *	\li \c control(4) = relax_poisson  : relaxation factor in 
     *				   \c Map_et::poisson  
     *  @param diff [output]   1-D \c Tbl  of size 8 for the storage of 
     *			    some error indicators : 
     *	    \li \c diff(0)  : Relative change in the enthalpy field 1
     *			      between two successive steps 
     *	    \li \c diff(1)  : Relative change in the enthalpy field 2
     *			      between two successive steps 
     *	    \li \c diff(2)  : Relative error in the resolution of the
     *			    Poisson equation for \c nuf    
     *	    \li \c diff(3)  : Relative error in the resolution of the
     *			    Poisson equation for \c nuq    
     *	    \li \c diff(4)  : Relative error in the resolution of the
     *			    Poisson equation for \c dzeta    
     *	    \li \c diff(5)  : Relative error in the resolution of the
     *			    Poisson equation for \c tggg    
     *	    \li \c diff(6)  : Relative error in the resolution of the
     *			    equation for \c shift  (x comp.)   
     *	    \li \c diff(7)  : Relative error in the resolution of the
     *			    equation for \c shift  (y comp.)   
     */
    void equilibrium_bi(double ent_c, double ent_c2, double omega0, 
			double omega20, const Tbl& ent_limit, 
			const Tbl& ent2_limit, const Itbl& icontrol, 
			const Tbl& control, Tbl& diff,
			int mer_mass, double mbar1_wanted, double mbar2_wanted, 
			double aexp_mass);
  };
  
}
#endif
