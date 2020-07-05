/*
 *  Definition of Lorene class Excision_surf, friend class of Spheroid.
 *
 */

/*
 *   Copyright (c) 2008  Jose-Luis Jaramillo & Nicolas Vasset
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

#ifndef __EXCISIONSURF_H_ 
#define __EXCISIONSURF_H_ 

/*
 * $Header: /cvsroot/Lorene/C++/Include/excision_surf.h,v 1.7 2014/10/13 08:52:34 j_novak Exp $
 *
 */
#include "metric.h"
#include "spheroid.h"

namespace Lorene {
/**
 * Surface where boundary conditions for quantities in the bulk will be calculated
 * It relies on geometrical properties of the associated Spheroid\ingroup (star)
 * (*** WARNING! under development***)
 * 
 */
class Excision_surf {


    // Data : 
    // -----
 protected:
  /// The associated Spheroid object
  Spheroid sph ;
  
  ///The value of the conformal factor on the 3-slice
  Scalar conf_fact ; 
  
  /** The lapse defined on the 3 slice*/
  Scalar lapse ;
  
  /** The Shift 3-vector on the slice **/
  Vector shift ;
  
  /** The 3-d metric on the slice*/
  Metric gamij ; 
  
  /** The 3-d extrinsic curvature on the slice */
  Sym_tensor Kij ;
  
  /** The time step for evolution in parabolic drivers */
  double delta_t;

  /** The internal number of timesteps for one iteration. */
  double no_of_steps;

  /** The 2d expansion, directly evolved from the initial excision with Einstein Equations */
  Scalar expa;
  
  /** The time derivative of the expansion, derived from Einstein equations and arbitrary evolution */
  Scalar dt_expa;
  

    // Derived data : 
    // ------------
 protected:

	mutable Scalar* p_get_BC_conf_fact_1 ; ///< Source of Neumann boundary condition on \f$ \psi \f$, 
	mutable Scalar* p_get_BC_lapse_1 ;   ///< Source of Dirichlet boundary condition of \f$ N \f$
	mutable Vector* p_get_BC_shift_1 ; ///< Source of Dirichlet BC for the shift vector  \f$ \beta^{i} \f$
	mutable Scalar* p_get_BC_Npsi_1 ; ///<  Source of Neumann boundary condition on \f$ \psi \f$.
	mutable Scalar* p_get_BC_conf_fact_2 ; ///< Source of Neumann boundary condition on \f$ \psi \f$,  
	mutable Scalar* p_get_BC_conf_fact_3 ; ///< Source of Neumann boundary condition on \f$ \psi \f$, 
	mutable Scalar* p_get_BC_conf_fact_4 ; ///< Source of Birichlet boundary condition on \f$ \psi \f$, 
	mutable Scalar* p_get_BC_lapse_2 ;   ///< Source of Dirichlet boundary condition of \f$ N \f$
	mutable Scalar* p_get_BC_lapse_3 ;   ///< Source of Dirichlet condtion on  \f$ N \f$, based on einstein equations.
	mutable Scalar* p_get_BC_lapse_4 ;   ///< Source of Dirichlet condtion on  \f$ N \f$, based on einstein equations (conservation of isotropic gauge)
	mutable Scalar* p_derive_t_expa ;   ///< Computation of an updated expansion scalar.
	mutable Vector* p_get_BC_shift_2 ; ///< Source of Dirichlet BC for the shift vector  \f$ \beta^{i} \f$
	mutable Vector* p_get_BC_shift_3 ; ///< Source of Dirichlet BC for the shift vector  \f$ \beta^{i} \f$, partly derived from kinematical relation
	mutable Vector* p_get_BC_shift_4 ; ///< Source of Dirichlet BC for the shift vector  \f$ \beta^{i} \f$, partly from projection of Einstein Equations
	mutable Scalar* p_get_BC_Npsi_2 ; ///<  Source of Dirichlet boundary condition on \f$ N \psi \f$.
	mutable Scalar* p_get_BC_Npsi_3 ; ///<  Source of Dirichlet boundary condition on \f$ N \psi \f$.
	mutable Scalar* p_get_BC_Npsi_4 ; ///<  Source of Dirichlet boundary condition on \f$ N \psi \f$.
	mutable Scalar* p_get_BC_Npsi_5 ; ///<  Source of Neumann boundary condition on \f$ N \psi \f$.


    // Constructors - Destructor
    // -------------------------
    public:

	/** Constructor of an excision surface embedded in a 3-slice (\c Time_slice ) of 3+1 
	 * formalism. 
	 * This is done from the \c Time_slice data.
	 * @param h_in : the location of the surface \e r = h_in (WARNING:must be 
	                    defined on a mono-domain angular grid)
	 * @param gij : the 3-metric on the 3slice
	 * @param Kij : the extrinsic curvature of the 3-slice 
	 *                 (covariant representation)
	 * @param timestep : time interval associated with the parabolic-driven boundary conditions.
	 * @param int_nos : Number of iterations to be done during timestep.
	 */
	Excision_surf(const Scalar& h_in, const Metric& gij, const Sym_tensor& Kij2, const Scalar& ppsi, const Scalar& nn, const Vector& beta, double timestep, int int_nos) ;
	Excision_surf(const Excision_surf& ) ;		///< Copy constructor

	/// Constructor from a file (see \c sauve(FILE*) )
	Excision_surf(FILE* ) ;    		

	virtual ~Excision_surf() ;			///< Destructor
 

    // Memory management
    // -----------------
    protected:
	/// Deletes all the derived quantities
	virtual void del_deriv() const ; 
	
	/// Sets to \c 0x0 all the pointers on derived quantities
	void set_der_0x0() const ; 



    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another Excision_surf
	void operator=(const Excision_surf&) ;	


	// Functions/Accessors
	//---------------------

	/** Computes the parameters for the hyperbolic evolution in set_expa_hyperb(), so that the expansion has a C1 matching 
	 * with initial data. Sets also values for expa() and dt_expa() accordingly with initial conditions */
	void get_evol_params_from_ID(double alpha, double beta, double gamma, Scalar& Ee, Vector& Jj, Sym_tensor& Ss) ;

	///Sets a new value for expansion rescaled over lapse (and its derivative), obtained by parabolic evolution.
	void set_expa_parab(double c_theta_lap, double c_theta_fin, Scalar& expa_fin) ;
	/**Sets a new value for expansion rescaled over lapse (and its derivative), obtained by hyperbolic evolution.
	 * Parameters for the hyperbolic driver are determined by the function Excision_surf::get_evol_params_from_ID() 
	 * so that the expansion stays of regularity $C^{1}$ throughout.
	 * All manipulated quantities are 2-dimensional.
	 */
	void set_expa_hyperb(double alph0, double beta0, double gamma0) ;


    // Accessors
    // ---------
    public:
	/// Returns the spheroid
	const Spheroid& get_sph() const {return sph; }; 
	
	/// Returns the conformal factor associated with the surface
	const Scalar& get_conf_fact() const {return conf_fact; } ;

	/// Returns the lapse function
	const Scalar& get_lapse() const {return lapse ; } ;

	/// Returns the shift vector field.
	const Vector& get_shift() const {return shift ; } ;

	/// Returns the symmetric tensor \f$ gamij \f$
	const Metric& get_gamij() const {return gamij ; } ;

	/// returns the 3-d extrinsic curvature \f$ K_{ij}\f$
	const Sym_tensor& get_Kij() const {return Kij ; } ;

	/// Returns the timestep used for evolution.
	double get_delta_t() const {return delta_t ;};

	/// Returns the internal number of timesteps for one iteration.
	double get_no_of_steps() const {return no_of_steps ;};

	/// Returns the assumed expansion associated to the excised surface at t.
	const Scalar & get_expa() const {return expa; };

	/// Returns the assumed time derivative of the expansion at t, evolved by functions of this class; 
	const Scalar& get_dt_expa() const {return dt_expa;} ;

	/// Sets a new spheroid from data
	Spheroid& set_sph() {del_deriv() ; return sph ;};

	/// Sets the value of the conformal factor
	Scalar& set_conf_fact() {del_deriv() ; return conf_fact ; } ;

	/// Sets the lapse function
	Scalar& set_lapse() {del_deriv() ; return lapse ; } ;

	/// Sets the shift vector field
	Vector& set_shift() {del_deriv() ; return shift ; } ;

	/// Sets the 3d metric of the TimeSlice
	Metric& set_gamij() {del_deriv() ; return gamij ; } ;

	/// Sets the extrinsic curvature
	Sym_tensor& set_Kij() {del_deriv() ; return Kij ; } ;

	double set_delta_t() {del_deriv() ; return delta_t ; } ;

	double set_no_of_steps() {del_deriv() ; return no_of_steps ; } ;

	/// Sets the expansion function on the surface at time t (considering to protect this function)

	Scalar& set_expa() {del_deriv(); return expa;};

	/// Sets the time derivative of the expansion function on the surface at time t (considering to protect this function)

	Scalar& set_dt_expa() {del_deriv(); return dt_expa;};
	
    // Computational functions
    // -----------------------
    public:
	/// Source for a Neumann BC on the conformal factor. If boolean isMOTS is false, it is based on expansion value of the spheroid or the value of exppa; it is based on zero expansion if isMOTS is true. 
	const Scalar& get_BC_conf_fact_1(bool isMOTS = false) const ;
// Source for an arbitrary Dirichlet BC on the lapse	
	const Scalar& get_BC_lapse_1(double value) const ;
// Source for a global Dirichlet BC on the shift, imposing a conformal Killing symmetry on \f$ \varphi \f$.
	const Vector& get_BC_shift_1(double Omega) const ;
// Source for a Dirichlet arbitrary BC on (N*Psi1)
	const Scalar& get_BC_Npsi_1(double value) const ;
	/// Source for the Dirichlet BC on the conformal factor, based on a parabolic driver for the conformal factor
	const Scalar& get_BC_conf_fact_2(double c_psi_lap, double c_psi_fin, Scalar& expa_fin) const ;
/// Source for the Neumann BC on the conformal factor, based on a parabolic driver for the expansio
	const Scalar& get_BC_conf_fact_3(double c_theta_lap, double c_theta_fin, Scalar& expa_fin) const ;
/// Source for the Dirchlet BC on the conformal factor, based on the consistency condition derived from the trace
	const Scalar& get_BC_conf_fact_4() const ;
/// Source for Dirichlet BC on the lapse, based on a parabolic driver towards arbitrary constant value
	const Scalar& get_BC_lapse_2(double lapse_fin, double c_lapse_lap, double c_lapse_fi) const ;
/// Source for Dirichlet BC on the lapse, based on einstein equations
	const Scalar& get_BC_lapse_3(Scalar& dttheta, Scalar& Ee, Vector& Jj, Sym_tensor& Sij, bool sph_sym = true) const ;
	/// Source for Dirichlet BC on the lapse, based on einstein equations (conservation of isotropic gauge)
	const Scalar& get_BC_lapse_4 (Scalar& old_nn, Vector& beta_point, Sym_tensor& strain_tens) const ;
	/// Forms the prospective time derivative for the expansion using projected Einstein equations. Does NOT modify the member dt_expa: do it by hand!
	const Scalar& derive_t_expa(Scalar& Ee, Vector& Jj, Sym_tensor& Sij) const;
	/// Source for a Dirichlet BC on the shift, based on a Parabolic driver; no assumptions are made except a global conformal Killing symmetry.
	const Vector& get_BC_shift_2(double c_bb_lap, double c_bb_fin, double c_V_lap, double epsilon) const ;
	/// Source for a Dirichlet BC on the shift, based on a Parabolic driver; Radial part is dealt with using a kinematical relation.
	const Vector& get_BC_shift_3(Scalar& dtpsi, double c_V_lap, double epsilon) const ;
	/// Source for a Dirichlet BC on the shift, based on a Parabolic driver; Radial part is dealt with using projection of Einstein Equations.
	const Vector & get_BC_shift_4(Scalar& dttheta, Scalar& Ee, Vector& Jj, Sym_tensor& Sij, double c_V_lap, double epsilon, bool sph_sym= true) const ;
	/// Source for the Dirichlet BC on (N*Psi1), based on a parabolic driver.
	const Scalar& get_BC_Npsi_2(double value, double c_npsi_lap, double c_npsi_fin) const ;        
	/// Source for the Dirichlet BC on (N*Psi1), with Kerr_Schild-like form for the lapse boundary.
	const Scalar& get_BC_Npsi_3(double n_0, double beta) const ;      
	/// Source for a Dirichlet BC on (N*Psi1), fixing a constant surface gravity in space and time.
	const Scalar& get_BC_Npsi_4(double Kappa) const ;      
	/// Source for a Neumann BC on (N*Psi1), fixing a constant surface gravity in space and time.
	const Scalar& get_BC_Npsi_5(double Kappa) const ; 
    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a file
    
	/// Display
	friend ostream& operator<<(ostream& , const Spheroid& ) ;	

};

ostream& operator<<(ostream& , const Spheroid& ) ;


}
#endif
