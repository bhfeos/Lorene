 /*
 *  Definition of Lorene class Vector
 *
 */

/*
 *   Copyright (c) 2003  Eric Gourgoulhon & Jerome Novak
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

#ifndef __VECTOR_H_ 
#define __VECTOR_H_ 

/*
 * $Id: vector.h,v 1.42 2014/10/13 08:52:37 j_novak Exp $
 * $Log: vector.h,v $
 * Revision 1.42  2014/10/13 08:52:37  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.41  2008/12/03 10:18:56  j_novak
 * Method 6 is now the default for calls to vector Poisson solver.
 *
 * Revision 1.40  2008/10/29 08:19:08  jl_cornou
 * Typo in the doxygen documentation + spectral bases for pseudo vectors added
 * and curl
 *
 * Revision 1.39  2008/08/27 08:45:21  jl_cornou
 * Implemented routines to solve dirac systems for divergence-free vectors
 *
 * Revision 1.38  2007/12/21 16:06:16  j_novak
 * Methods to filter Tensor, Vector and Sym_tensor objects.
 *
 * Revision 1.37  2007/05/11 11:38:16  n_vasset
 * Added poisson_boundary2.C routine
 *
 * Revision 1.36  2007/01/16 15:05:59  n_vasset
 * New constructor (taking a Scalar in mono-domain angular grid for
 * boundary) for function sol_elliptic_boundary
 *
 * Revision 1.35  2005/06/09 07:56:25  f_limousin
 * Implement a new function sol_elliptic_boundary() and
 * Vector::poisson_boundary(...) which solve the vectorial poisson
 * equation (method 6) with an inner boundary condition.
 *
 * Revision 1.34  2005/02/16 15:00:05  m_forot
 * Add visu_streamlime function
 *
 * Revision 1.33  2005/02/14 13:01:48  j_novak
 * p_eta and p_mu are members of the class Vector. Most of associated functions
 * have been moved from the class Vector_divfree to the class Vector.
 *
 * Revision 1.32  2004/05/25 14:54:01  f_limousin
 * Change arguments of method poisson with parameters.
 *
 * Revision 1.31  2004/05/09 20:54:22  e_gourgoulhon
 * Added method flux (to compute the flux accross a sphere).
 *
 * Revision 1.30  2004/03/29 11:57:13  e_gourgoulhon
 * Added methods ope_killing and ope_killing_conf.
 *
 * Revision 1.29  2004/03/28 21:25:14  e_gourgoulhon
 * Minor modif comments (formula for V^\theta in Vector_divfree).
 *
 * Revision 1.28  2004/03/27 20:59:55  e_gourgoulhon
 * Slight modif comment (doxygen \ingroup).
 *
 * Revision 1.27  2004/03/22 13:12:44  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.26  2004/03/10 12:52:19  f_limousin
 * Add a new argument "method" in poisson method.
 *
 * Revision 1.25  2004/03/03 09:07:02  j_novak
 * In Vector::poisson(double, int), the flat metric is taken from the mapping.
 *
 * Revision 1.24  2004/02/26 22:45:44  e_gourgoulhon
 * Added method derive_lie.
 *
 * Revision 1.23  2004/02/22 15:47:45  j_novak
 * Added 2 more methods to solve the vector Poisson equation. Method 1 is not
 * tested yet.
 *
 * Revision 1.22  2004/02/21 16:27:53  j_novak
 * Modif. comments
 *
 * Revision 1.21  2004/02/16 17:40:14  j_novak
 * Added a version of poisson with the flat metric as argument (avoids
 * unnecessary calculations by decompose_div)
 *
 * Revision 1.20  2003/12/14 21:47:24  e_gourgoulhon
 * Added method visu_arrows for visualization through OpenDX.
 *
 * Revision 1.19  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.18  2003/11/03 22:31:02  e_gourgoulhon
 * Class Vector_divfree: parameters of methods set_vr_eta_mu and set_vr_mu
 * are now all references.
 *
 * Revision 1.17  2003/11/03 14:02:17  e_gourgoulhon
 * Class Vector_divfree: the members p_eta and p_mu are no longer declared
 * "const".
 *
 * Revision 1.16  2003/10/20 15:12:17  j_novak
 * New method Vector::poisson
 *
 * Revision 1.15  2003/10/20 14:44:49  e_gourgoulhon
 * Class Vector_divfree: added method poisson().
 *
 * Revision 1.14  2003/10/20 09:32:10  j_novak
 * Members p_potential and p_div_free of the Helmholtz decomposition
 * + the method decompose_div(Metric).
 *
 * Revision 1.13  2003/10/17 16:36:53  e_gourgoulhon
 * Class Vector_divfree: Added new methods set_vr_eta_mu and set_vr_mu.
 *
 * Revision 1.12  2003/10/16 21:36:01  e_gourgoulhon
 * Added method Vector_divfree::update_vtvp().
 *
 * Revision 1.11  2003/10/16 15:25:00  e_gourgoulhon
 * Changes in documentation.
 *
 * Revision 1.10  2003/10/16 14:21:33  j_novak
 * The calculation of the divergence of a Tensor is now possible.
 *
 * Revision 1.9  2003/10/13 20:49:26  e_gourgoulhon
 * Corrected typo in the comments.
 *
 * Revision 1.8  2003/10/13 13:52:39  j_novak
 * Better managment of derived quantities.
 *
 * Revision 1.7  2003/10/12 20:33:15  e_gourgoulhon
 * Added new derived class Vector_divfree (preliminary).
 *
 * Revision 1.6  2003/10/06 13:58:45  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.5  2003/10/05 21:07:27  e_gourgoulhon
 * Method std_spectral_base() is now virtual.
 *
 * Revision 1.4  2003/10/03 14:08:03  e_gourgoulhon
 * Added constructor from Tensor.
 *
 * Revision 1.3  2003/09/29 13:48:17  j_novak
 * New class Delta.
 *
 * Revision 1.2  2003/09/29 12:52:56  j_novak
 * Methods for changing the triad are implemented.
 *
 * Revision 1.1  2003/09/26 08:07:32  j_novak
 * New class vector
 *
 *
 * $Header: /cvsroot/Lorene/C++/Include/vector.h,v 1.42 2014/10/13 08:52:37 j_novak Exp $
 *
 */

namespace Lorene {
class Vector_divfree ;

			//-------------------------//
			//       class Vector      //
			//-------------------------//
			

/**
 * Tensor field of valence 1. \ingroup (tensor)
 * 
 */
class Vector: public Tensor {

    // Derived data : 
    // ------------
    protected:
        /** The potential \f$\phi\f$ giving the gradient part in the Helmholtz 
	 * decomposition of any 3D vector \f$\vec{V}: \quad \vec{V} = 
	 * \vec{\nabla} \phi + \vec{\nabla} \wedge \vec{\psi}\f$.
	 * Only in the case of contravariant vectors.
	 */
        mutable Scalar* p_potential[N_MET_MAX] ;

	/** The divergence-free vector \f$\vec{W} =  \vec{\nabla} \wedge 
	 * \vec{\psi}\f$ of the Helmholtz decomposition of any 3D vector 
	 *\f$\vec{V}: \quad \vec{V} = \vec{\nabla} \phi + \vec{\nabla} 
	 *\wedge \vec{\psi}\f$. Only in the case of contravariant vectors.
	 */
	mutable Vector_divfree* p_div_free[N_MET_MAX] ;

 	/** Field \f$\eta\f$ such that the angular components \f$(V^\theta, V^\varphi)\f$
	 * of the vector are written:
	 * \f[
	 *	V^\theta =   {\partial \eta \over \partial\theta} -
	 *		 {1\over\sin\theta} {\partial \mu \over \partial\varphi} 
	 * \f] 
	 * \f[
	 *	V^\varphi =  {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} 
	 * \f] 
	 */
	mutable Scalar* p_eta ;
	
	/** Field \f$\mu\f$ such that the angular components \f$(V^\theta, V^\varphi)\f$
	 * of the vector are written:
	 * \f[
	 *	V^\theta =  {\partial \eta \over \partial\theta} -
	 *	 {1\over\sin\theta} {\partial \mu \over \partial\varphi} 
	 * \f] 
	 * \f[
	 *	V^\varphi =  {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} 
	 * \f] 
	 */
	mutable Scalar* p_mu ;

	/** Field \f$A\f$ defined by 
	 * \f[ 
	 *     A = {\partial \eta \over \ partial r} + { \eta \over r} - {V^r \over r}
	 *  \f]
	 * Insensitive to the longitudinal part of the vector, related to the curl.
	 */
	mutable Scalar* p_A ;
 	
   // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param tipe  the type \c COV  for a covariant vector (1-form) 
	 *		and \c CON  for a contravariant one
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *		    the vector components are defined 
	 */
	Vector(const Map& map, int tipe, const Base_vect& triad_i) ;

	/** Standard constructor with the triad passed as a pointer.
	 * 
	 * @param map   the mapping 
	 * @param tipe  the type \c COV  for a covariant vector (1-form) 
	 *		and \c CON  for a contravariant one
	 * @param triad_i  pointer on the vectorial basis (triad) 
	 * with respect to which the vector components are defined 
	 */
	Vector(const Map& map, int tipe, const Base_vect* triad_i) ;

	Vector(const Vector& a) ;       ///< Copy constructor

	/** Constructor from a \c Tensor .
	 *  The \c Tensor  must be of valence one.
	 */
	Vector(const Tensor& a) ;	

	/** Constructor from a file (see \c Tensor::sauve(FILE*) ).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param fich  file which has been created by 
	 *			    the function \c sauve(FILE*) .
	 */
	Vector(const Map& map, const Base_vect& triad_i, FILE* fich) ;

	virtual ~Vector() ;			///< Destructor

 
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const;	///< Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 

	/**
	 * Logical destructor of the derivatives depending on the i-th
	 * element of \c met_depend  in the class \c Vector.
	 */	
	virtual void del_derive_met(int) const ;

	/**
	 * Sets all the i-th components of \c met_depend  in the 
	 * class \c Vector  (\c p_potential , etc...) to 0x0.
	 */
	void set_der_met_0x0(int) const ;


    // Mutators / assignment
    // ---------------------
    public:
	/** Sets a new vectorial basis (triad) of decomposition and modifies
	 *  the components accordingly. 
	 */
	virtual void change_triad(const Base_vect& ) ; 

	/// Assignment from a Vector
	virtual void operator=(const Vector& a) ;	
	
	/// Assignment from a Tensor
	virtual void operator=(const Tensor& a) ;	

	/** Defines the components through potentials \f$\eta\f$ and \f$\mu\f$ 
	 *  (see members \c p_eta  and \c p_mu ), 
	 *  as well as the \f$V^r\f$ component of the vector. 
	 *
	 *	@param vr_i [input] component \f$V^r\f$ of the vector
	 *	@param eta_i [input] angular potential \f$\eta\f$
	 *	@param mu_i [input] angular potential \f$\mu\f$
	 *
	 */
	void set_vr_eta_mu(const Scalar& vr_i, const Scalar& eta_i,
		const Scalar& mu_i) ; 

	/**Makes the Helmholtz decomposition (see documentation of
	 * \c p_potential ) of \c this  with respect to a given
	 * \c Metric , only in the case of contravariant vectors.
	 */
	void decompose_div(const Metric&) const ;

	/**Returns the potential in the Helmholtz decomposition.
	 *
	 * It first makes the Helmholtz decomposition (see documentation of
	 * \c p_potential ) of \c this  with respect to a given
	 * \c Metric  and then returns \f$\phi\f$. Only in the case
	 * of contravariant vectors.
	 */
	const Scalar& potential(const Metric& ) const ;
	
	/**Returns the div-free vector in the Helmholtz decomposition.
	 *
	 * It first makes the Helmholtz decomposition (see documentation of
	 * \c p_potential ) of \c this  with respect to a given
	 * \c Metric  and then returns \f$\vec{W}\f$. Only in the case
	 * of contravariant vectors.	
	 */
	const Vector_divfree& div_free(const Metric& ) const;

	/** Applies exponential filters to all components 
	 * (see \c Scalar::exponential_filter_r ). Does a loop for Cartesian 
	 * components, and works in terms of the r-component, \f$\eta\f$ and
	 * \f$\mu\f$ for spherical components.
	 */
	virtual void  exponential_filter_r(int lzmin, int lzmax, int p, 
			    double alpha= -16.) ;

	/** Applies exponential filters to all components 
	 * (see \c Scalar::exponential_filter_ylm ). Does a loop for Cartesian 
	 * components, and works in terms of the r-component, \f$\eta\f$ and
	 * \f$\mu\f$ for spherical components. 
	 */
	virtual void exponential_filter_ylm(int lzmin, int lzmax, int p, 
			    double alpha= -16.) ;

	
    // Accessors
    // ---------
    public:
	Scalar& set(int ) ; ///< Read/write access to a component

	const Scalar& operator()(int ) const; ///<Readonly access to a component

	/**
	 * Returns the position in the \c Scalar  array \c cmp  of a 
	 * component given by its index.  
	 *
	 * @return position in the \c Scalar  array \c cmp   
	 * corresponding to the index given in \c idx . \c idx 
	 * must be a 1-D \c Itbl  of size 1, the element of which 
	 * must be one of the spatial indices 1, 2 or 3. 
	 */
	virtual int position(const Itbl& idx) const {
	  assert (idx.get_ndim() == 1) ;
	  assert (idx.get_dim(0) == 1) ;
	  assert ((idx(0) >= 1) && (idx(0) <= 3)) ;

	  return (idx(0) - 1) ;
    
	} ;

	/**
	 * Returns the index of a component given by its position in the 
	 * \c Scalar  array \c cmp . 
	 *
	 * @return the index is stored in an 1-D array (\c Itbl ) of
	 *         size 1 giving its value for the component located at 
	 *         the position \c place  in the \c Scalar  array 
	 *         \c cmp . The element of this \c Itbl  
	 *	   corresponds to a spatial index 1, 2 or 3. 
	 */
	virtual Itbl indices(int place) const {
	  assert((place>=0) && (place<3)) ;
	  
	  Itbl res(1) ;
	  res = place + 1;
	  return res ;

	};
	
	/**
	 * Sets the standard spectal bases of decomposition for each component.
	 *
	 */
	virtual void std_spectral_base() ; 

	/**
	 * Sets the standard spectal bases of decomposition for each component for a pseudo_vector.
	 *
	 */
	virtual void pseudo_spectral_base() ; 

	/** Gives the field \f$\eta\f$ such that the angular components 
	 * \f$(V^\theta, V^\varphi)\f$ of the vector are written:
	 * \f[
	 *	V^\theta =  {\partial \eta \over \partial\theta} -
	 *	 {1\over\sin\theta} {\partial \mu \over \partial\varphi} 
	 * \f] 
	 * \f[
	 *	V^\varphi =  {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} 
	 * \f] 
	 */
	virtual const Scalar& eta() const ;
	
	/** Gives the field \f$\mu\f$ such that the angular components 
	 * \f$(V^\theta, V^\varphi)\f$ of the vector are written:
	 * \f[
	 *	V^\theta =  {\partial \eta \over \partial\theta} -
	 *	 {1\over\sin\theta} {\partial \mu \over \partial\varphi}
	 * \f] 
	 * \f[
	 *	V^\varphi =  {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} 
	 * \f] 
	 */
	virtual const Scalar& mu() const ;


	/** Gives the field \f$A\f$ defined by 
	 * \f[ 
	 *     A = {\partial \eta \over \partial r} + { \eta \over r} - {V^r \over r}
	 *  \f]
	 *  Related to the curl, A is insensitive to the longitudinal part of the vector.
	 */
	virtual const Scalar& A() const ;

	
	/** Computes the components \f$V^\theta\f$ and \f$V^\varphi\f$ from the
	 *  potential \f$\eta\f$ and  \f$\mu\f$, according to:
	 * \f[
	 *	V^\theta =  {\partial \eta \over \partial\theta} - 
	 *		{1\over\sin\theta} {\partial \mu \over \partial\varphi}
	 * \f] 
	 * \f[
	 *	V^\varphi =  {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} 
	 * \f] 
	 */
	void update_vtvp() ;

	
    // Differential operators/ PDE solvers
    // -----------------------------------
    public:
	/**The divergence of \c this  with respect to a \c Metric .
	 * The \c Vector  is assumed to be contravariant.
	 */
	const Scalar& divergence(const Metric&) const ; 

	/** The curl of \c this with respect to a (flat) \c Metric .
	 *  The \c Vector is assumed to be contravariant.
	 */
	const Vector_divfree curl() const ;

    /** Computes the Lie derivative of \c this  with respect to some
     *  vector field \c v 
     */
    Vector derive_lie(const Vector& v) const ; 
        
    /** Computes the Killing operator associated with a given metric.
     *  The Killing operator is defined by \f$ D^i V^j + D^j V^i \f$
     *  for a contravariant vector and by \f$ D_i V_j + D_j V_i \f$
     *  for a covariant vector.
     * @param gam metric with respect to which the covariant derivative
     *  \f$ D_i \f$ is defined. 
     */
     Sym_tensor ope_killing(const Metric& gam) const ; 
     
    /** Computes the conformal Killing operator associated with a given metric.
     *  The conformal Killing operator is defined by 
     *  \f$ D^i V^j + D^j V^i - \frac{2}{3} D_k V^k \, \gamma^{ij} \f$
     *  for a contravariant vector and by 
     *  \f$ D_i V_j + D_j V_i - \frac{2}{3} D^k V_k \, \gamma_{ij}\f$
     *  for a covariant vector.
     * @param gam metric \f$\gamma_{ij}\f$ with respect to which the covariant 
     *  derivative \f$ D_i \f$ is defined. 
     */
     Sym_tensor ope_killing_conf(const Metric& gam) const ; 

    /**Solves the vector Poisson equation with \c *this  as a source.
     * 
     * The equation solved is \f$\Delta N^i +\lambda \nabla^i 
     * \nabla_k N^k = S^i\f$.
     * \c *this  must be given with \c dzpuis  = 4.
     * It uses the Helmholtz decomposition (see documentation of
     * \c p_potential ), with a flat metric, deduced from the triad.
     *
     * @param lambda [input] \f$\lambda\f$.
     * @param method [input] method used to solve the equation 
     * (see Vector::poisson(double, Metric_flat, int) for details).
     *
     * @return the solution \f$N^i\f$.
     */
    Vector poisson(double lambda, int method = 6) const ;
     
    /**Solves the vector Poisson equation with \c *this  as a source.
     * 
     * The equation solved is \f$\Delta N^i +\lambda \nabla^i 
     * \nabla_k N^k = S^i\f$.
     * \c *this  must be given with \c dzpuis  = 4.
     * It uses the Helmholtz decomposition (see documentation of
     * \c p_potential ), with the flat metric \c met_f  given 
     * in argument.
     *
     * @param lambda [input] \f$\lambda\f$.
     * @param met_f [input] the flat metric for the Helmholtz decomposition.
     * @param method [input] method used to solve the equation:
     *        \li 0 : It uses the Helmholtz decomposition (see documentation of
     *            \c p_potential ), with the flat metric \c met_f  given 
     *            in argument (the default).
     *        \li 1 : It solves, first for the divergence (calculated using 
     *            \c met_f ), then the \e r -component, the \f$\eta\f$ 
     *            potential, and fianlly the \f$\mu\f$ potential (see documentation
     *            of \c Vector_div_free .
     *        \li 2 : The sources is transformed to cartesian components and the 
     *            equation is solved using Shibata method (see Granclement 
     *            \e et \e al. JCPH 2001.
     *        \li 6 : Solves for the \e r -component and \f$ \eta \f$ together in a
     *            system, and for the \f$ \mu \f$ potential (which decouples). The
     *            solution is then built from these fields through the method
     *            \c Vector::set_vr_eta_mu(). It is the default method.
     *
     * @return the solution \f$N^i\f$.
     */
    Vector poisson(double lambda, const Metric_flat& met_f, int method = 6) const ;
    
    /**Solves the vector Poisson equation with \c *this  as a source
     * and parameters controlling the solution.
     * 
     * The equatiopn solved is \f$\Delta N^i +\lambda \nabla^i 
     * \nabla_k N^k = S^i\f$.
     * \c *this  must be given with \c dzpuis  = 4.
     * It uses the Helmholtz decomposition (see documentation of
     * \c p_potential ), with a flat metric, deduced from the triad.
     *
     * @param lambda [input] \f$\lambda\f$.
     *   @param par [input/output] possible parameters
     *   @param uu [input/output] solution \e u  with the 
     *              boundary condition \e u =0 at spatial infinity. 
     */
    
    Vector poisson(const double lambda, Param& par,
		   int method = 6) const ;
    
    /**Solves the vector Poisson equation with \c *this  as a source 
     * with a boundary condition on the excised sphere.
     * 
     * The equation solved is \f$\Delta N^i +\lambda \nabla^i 
     * \nabla_k N^k = S^i\f$.
     * \c *this  must be given with \c dzpuis  = 4.
     * It uses the Helmholtz decomposition (see documentation of
     * \c p_potential )
     *
     * @param lambda [input] \f$\lambda\f$.
     * @param resu [output] the solution \f$N^i\f$.
     */
    void poisson_boundary(double lambda, const Mtbl_cf& limit_vr, 
			  const Mtbl_cf& limit_eta, const Mtbl_cf& limit_mu, 
			  int num_front, double fact_dir, double fact_neu,
			  Vector& resu) const ;
 

    /**Alternative to previous poisson_boundary method for vectors ; 
     * this uses method 6 for vectorial solving, updated version 
     * (as in the poisson_vector_block routine). 
     * Boundary arguments are here required as scalar fields.
     */
    void poisson_boundary2(double lam, Vector& resu, Scalar boundvr, 
			   Scalar boundeta, Scalar boundmu, double dir_vr, 
			   double neum_vr, double dir_eta, double neum_eta, 
			   double dir_mu, double neum_mu ) const ;

 
    Vector poisson_dirichlet(double lambda, const Valeur& limit_vr, 
			  const Valeur& limit_vt, const Valeur& limit_vp, 
			  int num_front) const ;

     /**Solves the vector Poisson equation with \c *this  as a source 
     * with a boundary condition on the excised sphere.
     * 
     * The equation solved is \f$\Delta N^i +\lambda \nabla^i 
     * \nabla_k N^k = S^i\f$.
     * \c *this  must be given with \c dzpuis  = 4.
     * It uses the Helmholtz decomposition (see documentation of
     * \c p_potential )
     *
     * @param lambda [input] \f$\lambda\f$.
     * @param resu [output] the solution \f$N^i\f$.
     */
    Vector poisson_neumann(double lambda, const Valeur& limit_vr, 
			  const Valeur& limit_vt, const Valeur& limit_vp, 
			  int num_front) const ;

    /**Solves the vector Poisson equation with \c *this  as a source 
     * with a boundary condition on the excised sphere.
     * 
     * The equation solved is \f$\Delta N^i +\lambda \nabla^i 
     * \nabla_k N^k = S^i\f$.
     * \c *this  must be given with \c dzpuis  = 4.
     * It uses the Helmholtz decomposition (see documentation of
     * \c p_potential )
     *
     * @param lambda [input] \f$\lambda\f$.
     * @param resu [output] the solution \f$N^i\f$.
     */
    Vector poisson_robin(double lambda, const Valeur& limit_vr, 
			 const Valeur& limit_vt, const Valeur& limit_vp, 
			 double fact_dir, double fact_neu, 
			 int num_front) const ;
 

    /** Computes the flux of the vector accross a sphere \e r = const.
     *
     *  @param radius radius of the sphere \e S on which the flux is
     *      to be taken; the center of \e S is assumed to be the 
     *      center of the mapping (member \c mp).
     *      \c radius can take the value \c __infinity (to get the flux
     *      at spatial infinity).
     *  @param met metric \f$ \gamma \f$ giving the area element 
     *      of the sphere
     *  @return \f$ \oint_S V^i ds_i \f$, where \f$ V^i \f$ is the vector
     *  represented by \c *this and \f$ ds_i \f$ is the area element 
     *  induced on \e S by \f$ \gamma \f$.
     */
    double flux(double radius, const Metric& met) const ; 

    void poisson_block(double lambda, Vector& resu) const ;

        // Graphics
        // --------
 public:
  /** 3D visualization via OpenDX.
   *
   * @param xmin [input] defines with \c xmax  the x range of the visualization box 
   * @param xmax [input] defines with \c xmin  the x range of the visualization box 
   * @param ymin [input] defines with \c ymax  the y range of the visualization box 
   * @param ymax [input] defines with \c ymin  the y range of the visualization box 
   * @param zmin [input] defines with \c zmax  the z range of the visualization box 
   * @param zmax [input] defines with \c zmin  the z range of the visualization box 
   * @param title [input] title for the graph (for OpenDX legend)
   * @param filename [input] name for the file which will be the input for 
   *    OpenDX; the default 0x0 is transformed into "vector_arrows"
   * @param start_dx [input] determines whether OpenDX must be launched (as a
   *     subprocess) to view the field; if set to \c false , only input files
   *     for future usage of OpenDX are created 
   * @param nx [input] number of points in the x direction (uniform sampling)   
   * @param ny [input] number of points in the y direction (uniform sampling)   
   * @param nz [input] number of points in the z direction (uniform sampling)   
   *
   */
    void visu_arrows(double xmin, double xmax, double ymin, double ymax,
    double zmin, double zmax, const char* title0 = 0x0, 
    const char* filename0 = 0x0, bool start_dx = true, int nx = 8, int ny = 8, 
    int nz = 8) const ;    

    void visu_streamline(double xmin, double xmax, double ymin, double ymax,
    double zmin, double zmax, const char* title0 = 0x0, 
    const char* filename0 = 0x0, bool start_dx = true, int nx = 8, int ny = 8, 
    int nz = 8) const ;   
        
     
 
};



			//---------------------------------//
			//       class Vector_divfree      //
			//---------------------------------//
			

/**
 * Divergence-free vectors. \ingroup (tensor)
 *
 * This class is designed to store divergence-free vectors,
 * with the component expressed in a orthonormal spherical basis
 * \f$(e_r,e_\theta,e_\varphi)\f$.
 *
 * 
 */
class Vector_divfree: public Vector {

    // Data : 
    // -----
    protected:
	/// Metric with respect to which the divergence is defined
	const Metric* const met_div ; 
	
    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *		    the vector components are defined 
	 * @param met the metric with respect to which the divergence is defined
	 */
	Vector_divfree(const Map& map, const Base_vect& triad_i, 
		const Metric& met) ;

	Vector_divfree(const Vector_divfree& ) ;       ///< Copy constructor

	/** Constructor from a file (see \c Tensor::sauve(FILE*) ).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param met the metric with respect to which the divergence is defined
	 * @param fich  file which has been created by 
	 *			    the function \c sauve(FILE*) .
	 */
	Vector_divfree(const Map& map, const Base_vect& triad_i, 
		const Metric& met, FILE* fich) ;

	virtual ~Vector_divfree() ;			///< Destructor

 
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const;	///< Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------

	public:
	/// Assignment from another \c Vector_divfree 
	void operator=(const Vector_divfree& a) ;	
	
	/// Assignment from a \c Vector 
	virtual void operator=(const Vector& a) ;	
	
	/// Assignment from a \c Tensor 
	virtual void operator=(const Tensor& a) ;	
	
	/** Sets the angular potentials \f$\mu\f$ (see member
	 *  \c p_mu ), and the \f$V^r\f$ component
	 *  of the vector. The potential \f$\eta\f$ is then deduced from
	 *  \f$V^r\f$ by the divergence-free condition. 
	 *  The components \f$V^\theta\f$ and \f$V^\varphi\f$ are updated consistently
	 *  by a call to the method \c update_vtvp() .
	 *
	 *	@param vr_i [input] component \f$V^r\f$ of the vector
	 *	@param mu_i [input] angular potential \f$\mu\f$
	 *
	 */
	void set_vr_mu(const Scalar& vr_i, const Scalar& mu_i) ; 
	
	/** Defines the components through \f$V_r\f$, \f$\eta\f$ and \f$\mu\f$.
	 * (see members \c p_eta and \c p_mu ),
	 *
	 *	@param vr_i [input] r-component of the vector 
	 *	@param eta_i [input] Angular potential \f$\eta\f$
	 * 	@param mu_i [input] Angular potential \f$\mu\f$
	 *
	 */
	void set_vr_eta_mu(const Scalar& vr_i, const Scalar& eta_i, const Scalar& mu_i) ;

	/** Defines the components through potentials \f$A\f$ and \f$\mu\f$.
	 * (see members \c p_A and \c p_mu ),
	 * 
	 *	@param A_i [input] Angular potential \f$A\f$
	 * 	@param mu_i [input] Angular potential \f$\mu\f$
	 *
	 */
	void set_A_mu(const Scalar& A_i, const Scalar& mu_i, const Param* par_bc) ;	
	
	// Computational methods
	// ---------------------
	public:
	/** Gives the field \f$\eta\f$ such that the angular components 
	 * \f$(V^\theta, V^\varphi)\f$ of the vector are written:
	 * \f[
	 *	V^\theta =  {\partial \eta \over \partial\theta} -
	 *	 {1\over\sin\theta} {\partial \mu \over \partial\varphi} 
	 * \f] 
	 * \f[
	 *	V^\varphi =  {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} 
	 * \f] 
	 */
	virtual const Scalar& eta() const ;
	
	/** Computes the solution of a vectorial Poisson equation
	 *  with \c *this  \f$= \vec{V}\f$ as a source:
	 * \f[
	 *    \Delta \vec{W} = \vec{V}
	 * \f] 
	 * 
	 * @return solution \f$\vec{W}\f$ of the above equation with the boundary
	 *	condition \f$\vec{W}=0\f$ at spatial infinity.
	 */
	Vector_divfree poisson() const ; 

	/** Computes the components \f$V^r\f$ and \f$\eta\f$ from the potential
	 *  A and the divergence-free condition, according to : 
	 * \f[ 
	 *      {\partial \eta \over \ partial r} + { \eta \over r} 
	 *						- {V^r \over r} = A
	 * \f]
	 * \f[
	 * 	{\partial V^r \over \partial r} + {2 V^r \over r} 
	 * 		+ {1 \over r} \Delta_{\vartheta \varphi} \eta = 0
	 * \f]
	 */
	void update_etavr() ;
	
	/** Computes the solution of a vectorial Poisson equation
	 *  with \c *this  \f$= \vec{V}\f$ as a source:
	 * \f[
	 *    \Delta \vec{W} = \vec{V}
	 * \f] 
	 * 
	 * @return solution \f$\vec{W}\f$ of the above equation with the boundary
	 *	condition \f$\vec{W}=0\f$ at spatial infinity.
	 */
	Vector_divfree poisson(Param& par) const ; 
	protected:
	
	/** Solves a system of two-coupled first-order PDEs obtained from the 
	 *  divergence-free condition and the requirement that the potential A 
	 *  has a given value. The system reads : 
	 * \f[ 
	 *     {\partial \eta \over \partial r} + {\eta \over r} - {V^r \over r} = A
	 * \f]
	 * \f[ 
	 *     {\partial V^r \over \partial r} + {2 V^r \over r} 
	 * 		+ {1 \over r}\Delta_{\vartheta \varphi}\eta = 0
	 * \f]
	 */
	void sol_Dirac_A(const Scalar& aaa, Scalar& eta, Scalar& vr,
              const Param* par_bc = 0x0) const ;	

	/** Solves via a tau method a system of two-coupled first-order PDEs 
	 *  obtained from the divergence-free condition and the requirement
	 *  that the potential A has a given value. The system reads : 
	 * \f[ 
	 *     {\partial \eta \over \partial r} + {\eta \over r} - {V^r \over r} = A
	 * \f]
	 * \f[ 
	 *     {\partial V^r \over \partial r} + {2 V^r \over r} 
	 * 		+ {1 \over r}\Delta_{\vartheta \varphi}\eta = 0
	 * \f]
	 */
	void sol_Dirac_A_tau(const Scalar& aaa, Scalar& eta, Scalar& vr,
              const Param* par_bc = 0x0) const ;

       /** Solves via a poisson method a system of two-coupled first-order PDEs 
	 *  obtained from the divergence-free condition and the requirement
	 *  that the potential A has a given value. The system reads : 
	 * \f[ 
	 *     {\partial \eta \over \partial r} + {\eta \over r} - {V^r \over r} = A
	 * \f]
	 * \f[ 
	 *     {\partial V^r \over \partial r} + {2 V^r \over r} 
	 * 		+ {1 \over r}\Delta_{\vartheta \varphi}\eta = 0
	 * \f]
	 */
	void sol_Dirac_A_poisson(const Scalar& aaa, Scalar& eta, Scalar& vr,
              const Param* par_bc = 0x0) const ;

	/** Solves a one-domain system of two-coupled first-order PDEs obtained from the 
	 *  divergence-free condition and the requirement that the potential A 
	 *  has a given value. The system reads : 
	 * \f[ 
	 *     {\partial \eta \over \partial r} + {\eta \over r} - {V^r \over r} = A
	 * \f]
	 * \f[ 
	 *     {\partial V^r \over \partial r} + {2 V^r \over r} 
	 * 		+ {1 \over r}\Delta_{\vartheta \varphi}\eta = 0
	 * \f]
	 */
	void sol_Dirac_A_1z(const Scalar& aaa, Scalar& eta, Scalar& vr,
              const Param* par_bc = 0x0) const ;	
	
};





}
#endif
