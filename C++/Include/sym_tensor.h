/*
 *  Definition of Lorene class Sym_tensor, 
 *  as well as derived classes Sym_tensor_trans and Sym_tensor_tt
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
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

#ifndef __SYM_TENSOR_H_ 
#define __SYM_TENSOR_H_ 

/*
 * $Id: sym_tensor.h,v 1.50 2019/08/16 08:47:36 j_novak Exp $
 * $Log: sym_tensor.h,v $
 * Revision 1.50  2019/08/16 08:47:36  j_novak
 * *** empty log message ***
 *
 * Revision 1.49  2014/10/13 08:52:36  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.48  2010/10/11 10:23:03  j_novak
 * Removed methods Sym_tensor_trans::solve_hrr() and Sym_tensor_trans::set_WX_det_one(), as they are no longer relevant.
 *
 * Revision 1.47  2008/12/05 08:46:19  j_novak
 * New method Sym_tensor_trans_aux::set_tt_part_det_one.
 *
 * Revision 1.46  2008/12/03 10:18:56  j_novak
 * Method 6 is now the default for calls to vector Poisson solver.
 *
 * Revision 1.45  2008/08/20 14:39:53  n_vasset
 * New Dirac solvers handling degenerate elliptic operators on excised spacetimes.
 *
 * Revision 1.44  2007/12/21 16:06:16  j_novak
 * Methods to filter Tensor, Vector and Sym_tensor objects.
 *
 * Revision 1.43  2007/11/27 15:48:52  n_vasset
 * New member p_tilde_c for class Sym_tensor
 *
 * Revision 1.42  2007/05/04 16:43:50  n_vasset
 * adding of functions sol_Dirac_BC2 and sol_Dirac_A2
 *
 * Revision 1.41  2006/10/24 13:03:17  j_novak
 * New methods for the solution of the tensor wave equation. Perhaps, first
 * operational version...
 *
 * Revision 1.40  2006/08/31 12:13:21  j_novak
 * Added an argument of type Param to Sym_tensor_trans::sol_  rac_A().
 *
 * Revision 1.39  2006/06/20 12:07:13  j_novak
 * Improved execution speed for sol_Dirac_tildeB...
 *
 * Revision 1.38  2006/06/14 10:04:19  j_novak
 * New methods sol_Dirac_l01, set_AtB_det_one and set_AtB_trace_zero.
 *
 * Revision 1.37  2006/06/13 13:30:12  j_novak
 * New members sol_Dirac_A and sol_Dirac_tildeB (see documentation).
 *
 * Revision 1.36  2006/06/12 13:37:23  j_novak
 * Added bounds in l (multipolar momentum) for Sym_tensor_trans::solve_hrr.
 *
 * Revision 1.35  2006/06/12 07:42:28  j_novak
 * Fields A and tilde{B} are defined only for l>1.
 *
 * Revision 1.34  2006/06/12 07:27:18  j_novak
 * New members concerning A and tilde{B}, dealing with the transverse part of the
 * Sym_tensor.
 *
 * Revision 1.33  2005/11/28 14:45:14  j_novak
 * Improved solution of the Poisson tensor equation in the case of a transverse
 * tensor.
 *
 * Revision 1.32  2005/09/16 13:58:10  j_novak
 * New Poisson solver for a Sym_tensor_trans.
 *
 * Revision 1.31  2005/09/07 16:47:42  j_novak
 * Removed method Sym_tensor_trans::T_from_det_one
 * Modified Sym_tensor::set_auxiliary, so that it takes eta/r and mu/r as
 * arguments.
 * Modified Sym_tensor_trans::set_hrr_mu.
 * Added new protected method Sym_tensor_trans::solve_hrr
 *
 * Revision 1.30  2005/04/08 08:22:04  j_novak
 * New methods set_hrr_mu_det_one() and set_WX_det_one(). Not tested yet...
 *
 * Revision 1.29  2005/04/06 15:43:58  j_novak
 * New method Sym_tensor_trans::T_from_det_one(...).
 *
 * Revision 1.28  2005/04/04 15:25:22  j_novak
 * Added new members www, xxx, ttt and the associated methods.
 *
 * Revision 1.27  2005/04/01 14:28:31  j_novak
 * Members p_eta and p_mu are now defined in class Sym_tensor.
 *
 * Revision 1.26  2005/01/03 08:34:58  f_limousin
 * Come back to the previous version.
 *
 * Revision 1.25  2005/01/03 08:15:39  f_limousin
 * The first argument of the function trace_from_det_one(...) is now
 * a Sym_tensor_trans instead of a Sym_tensor_tt (because of a
 * compilation error with some compilators).
 *
 * Revision 1.24  2004/12/28 14:21:46  j_novak
 * Added the method Sym_tensor_trans::trace_from_det_one
 *
 * Revision 1.23  2004/12/28 10:37:22  j_novak
 * Better way of enforcing zero divergence.
 *
 * Revision 1.22  2004/06/14 20:44:44  e_gourgoulhon
 * Added argument method_poisson to Sym_tensor::longit_pot and
 * Sym_tensor::transverse.
 *
 * Revision 1.21  2004/05/25 14:57:20  f_limousin
 * Add parameters in argument of functions transverse, longit_pot,
 * set_tt_trace, tt_part and set_khi_mu for the case of a Map_et.
 *
 * Revision 1.20  2004/05/24 13:44:54  e_gourgoulhon
 * Added parameter dzp to method Sym_tensor_tt::update.
 *
 * Revision 1.19  2004/04/08 16:37:54  e_gourgoulhon
 * Sym_tensor_tt::set_khi_mu: added argument dzp (dzpuis of resulting h^{ij}).
 *
 * Revision 1.18  2004/03/30 14:01:19  j_novak
 * Copy constructors and operator= now copy the "derived" members.
 *
 * Revision 1.17  2004/03/29 16:13:06  j_novak
 * New methods set_longit_trans and set_tt_trace .
 *
 * Revision 1.16  2004/03/22 13:12:43  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.15  2004/03/03 13:54:16  j_novak
 * Error in comments corrected.
 *
 * Revision 1.14  2004/03/03 13:16:20  j_novak
 * New potential khi (p_khi) and the functions manipulating it.
 *
 * Revision 1.13  2004/02/26 22:45:13  e_gourgoulhon
 * Added method derive_lie.
 *
 * Revision 1.12  2004/02/18 18:43:22  e_gourgoulhon
 * Method trace() renamed the_trace() in order to avoid
 * any confusion with new method Tensor::trace().
 *
 * Revision 1.11  2004/01/04 20:49:06  e_gourgoulhon
 * Sym_tensor is now a derived class of Tensor_sym.
 * Suppressed methods Sym_tensor::indices and Sym_tensor::position:
 *  they are now implemented at the Tensor_sym level.
 *
 * Revision 1.10  2003/11/27 16:05:11  e_gourgoulhon
 * Changed return value of methods transverse( ) and longit_pot( ).
 *
 * Revision 1.9  2003/11/26 21:56:21  e_gourgoulhon
 * Class Sym_tensor: added the members p_transverse and p_longit_pot,
 * and the associated methods transverse( ), longit_pot( ),
 * del_deriv_met( ) and set_der_met_0x0( ).
 *
 * Revision 1.8  2003/11/07 16:54:23  e_gourgoulhon
 * Added method Sym_tensor_tt::poisson().
 *
 * Revision 1.7  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.6  2003/11/05 15:26:31  e_gourgoulhon
 * Modif documentation.
 *
 * Revision 1.5  2003/11/04 22:57:26  e_gourgoulhon
 * Class Sym_tensor_tt: method set_eta_mu renamed set_rr_eta_mu
 *    method update_tp() renamed update()
 *    added method set_rr_mu.
 *
 * Revision 1.4  2003/11/03 22:29:54  e_gourgoulhon
 * Class Sym_tensor_tt: added functions set_eta_mu and update_tp.
 *
 * Revision 1.3  2003/11/03 17:09:30  e_gourgoulhon
 * Class Sym_tensor_tt: added the methods eta() and mu().
 *
 * Revision 1.2  2003/10/28 21:22:51  e_gourgoulhon
 * Class Sym_tensor_trans: added methods trace() and tt_part().
 *
 * Revision 1.1  2003/10/27 10:45:19  e_gourgoulhon
 * New derived classes Sym_tensor_trans and Sym_tensor_tt.
 *
 *
 * $Header: /cvsroot/Lorene/C++/Include/sym_tensor.h,v 1.50 2019/08/16 08:47:36 j_novak Exp $
 *
 */

namespace Lorene {
class Sym_tensor_trans ;
class Sym_tensor_tt ;




			//---------------------------------//
			//        class Sym_tensor         //
			//---------------------------------//
			
/**
 * Class intended to describe valence-2 symmetric tensors.
 * The storage and the calculations are different and quicker than with an 
 * usual \c Tensor .
 * 
 * The valence must be 2. \ingroup (tensor)
 *
 */
class Sym_tensor : public Tensor_sym {

    // Derived data : 
    // ------------
    protected:
	/** Array of the transverse part \f${}^t T^{ij}\f$ of the tensor with respect 
	 * to various metrics, transverse meaning divergence-free with respect
	 * to a metric. Denoting \c *this  by \f$T^{ij}\f$, we then have
	 * \f[
	 *		T^{ij} = {}^t T^{ij} + \nabla^i W^j + \nabla^j W^i  
	 *		\qquad\mbox{with}\quad \nabla_j {}^t T^{ij} = 0 
	 *\f]
	 * where \f$\nabla_i\f$ denotes the covariant derivative with respect
	 * to the given metric and \f$W^i\f$ is the vector potential of the
	 * longitudinal part of \f$T^{ij}\f$ (member \c p_longit_pot  below)
	 */
	mutable Sym_tensor_trans* p_transverse[N_MET_MAX] ;

	/** Array of the vector potential of the
	 * longitudinal part of the tensor with respect 
	 * to various metrics (see documentation of member 
	 * \c p_transverse 
	 */
	mutable Vector* p_longit_pot[N_MET_MAX] ;

	/** Field \f$\eta\f$ such that the components \f$(T^{r\theta}, T^{r\varphi})\f$
	 * of the tensor are written (has only meaning with spherical components!):
	 * \f[
	 *	T^{r\theta} =  {1\over r} \left( {\partial \eta \over \partial\theta} -
	 *	{1\over\sin\theta} {\partial \mu \over \partial\varphi} \right) 
	 *\f] 
	 * \f[
	 *	T^{r\varphi} =  {1\over r} \left( {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} \right)
	 *\f] 
	 */
	mutable Scalar* p_eta ;
	
	/** Field \f$\mu\f$ such that the components \f$(T^{r\theta}, T^{r\varphi})\f$
	 * of the tensor are written (has only meaning with spherical components!):
	 * \f[
	 *	T^{r\theta} =  {1\over r} \left( {\partial \eta \over \partial\theta} -
	 *	 {1\over\sin\theta} {\partial \mu \over \partial\varphi} \right) 
	 *\f] 
	 * \f[
	 *	T^{r\varphi} =  {1\over r} \left( {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} \right)
	 *\f] 
	 */
	mutable Scalar* p_mu ;

	/** Field \e W such that the components \f$T^{\theta\theta}, 
	 * T^{\varphi\varphi}\f$ and \f$T^{\theta\varphi}\f$
	 * of the tensor are written (has only meaning with spherical components!):
	 * \f[
	 * \frac{1}{2}\left(T^{\theta\theta} - T^{\varphi\varphi} \right) 
	 * = \frac{\partial^2 W}{\partial\theta^2} - \frac{1}{\tan
	 * \theta} \frac{\partial W}{\partial \theta} - \frac{1}{\sin^2 \theta} 
	 * \frac{\partial^2 W}{\partial \varphi^2} - 2\frac{\partial}{\partial \theta} 
	 * \left( \frac{1}{\sin \theta} \frac{\partial X}{\partial \varphi} \right) ,
	 *\f] 
	 * \f[
	 *  T^{\theta\varphi} = \frac{\partial^2 X}{\partial\theta^2} - \frac{1}{\tan
	 * \theta} \frac{\partial X}{\partial \theta} - \frac{1}{\sin^2 \theta} 
	 * \frac{\partial^2 X}{\partial \varphi^2} + 2\frac{\partial}{\partial \theta} 
	 * \left( \frac{1}{\sin \theta} \frac{\partial W}{\partial \varphi} \right) .
	 *\f] 
	 */
	mutable Scalar* p_www ;

	/** Field \e X such that the components \f$T^{\theta\theta}, 
	 * T^{\varphi\varphi}\f$ and \f$T^{\theta\varphi}\f$
	 * of the tensor are written (has only meaning with spherical components!):
	 * \f[
	 * \frac{1}{2}\left(T^{\theta\theta} - T^{\varphi\varphi} \right) 
	 * = \frac{\partial^2 W}{\partial\theta^2} - \frac{1}{\tan
	 * \theta} \frac{\partial W}{\partial \theta} - \frac{1}{\sin^2 \theta} 
	 * \frac{\partial^2 W}{\partial \varphi^2} - 2\frac{\partial}{\partial \theta} 
	 * \left( \frac{1}{\sin \theta} \frac{\partial X}{\partial \varphi} \right) ,
	 *\f] 
	 * \f[
	 *  T^{\theta\varphi} = \frac{\partial^2 X}{\partial\theta^2} - \frac{1}{\tan
	 * \theta} \frac{\partial X}{\partial \theta} - \frac{1}{\sin^2 \theta} 
	 * \frac{\partial^2 X}{\partial \varphi^2} + 2\frac{\partial}{\partial \theta} 
	 * \left( \frac{1}{\sin \theta} \frac{\partial W}{\partial \varphi} \right) .
	 *\f] 
	 */
	mutable Scalar* p_xxx ;

	/// Field \e T defined as \f$ T = T^{\theta\theta} + T^{\varphi\varphi} \f$.
	mutable Scalar* p_ttt ;

	/** Field \e A defined from \e X and \f$\mu\f$ insensitive to the 
	 * longitudinal part of the \c Sym_tensor (only for \f$\ell \geq 2\f$).
	 * Its definition reads \f[
	 * A = \frac{\partial X}{\partial r} - \frac{\mu}{r^2}.
	 * \f] */
	mutable Scalar* p_aaa ;

	/** Field \f$ \tilde{B}\f$ defined from \f$ h^{rr}, \eta, W\f$ and \e h
	 * insensitive to the longitudinal part of the \c Sym_tensor.
	 * It is defined for each multipolar momentum \f$\ell \geq 2\f$ by
	 * \f[ 
	 * \tilde{B} = (\ell + 2) \frac{\partial W}{\partial r} + \ell(\ell + 2)
	 * \frac{W}{r} - \frac{2\eta}{r^2} + \frac{(\ell +2)T}{2r(\ell + 1)}
	 * + \frac{1}{2(\ell + 1)} \frac{\partial T}{\partial r} - \frac{h^{rr}}
	 * {(\ell + 1)r}.
	 * \f]
	 */
	mutable Scalar* p_tilde_b ;

	/** Field \f$ \tilde{C}\f$ defined from \f$ h^{rr}, \eta, W\f$ and \e h
	 * insensitive to the longitudinal part of the \c Sym_tensor.
	 * It is defined for each multipolar momentum \f$\ell \geq 2\f$ by
	 * \f[ 
	 * \tilde{C} = - (\ell - 1) \frac{\partial W}{\partial r} + (\ell + 1)(\ell - 1)
	 * \frac{W}{r} - \frac{2\eta}{r^2} + \frac{(\ell - 1)T}{2r\ell}
	 * - \frac{1}{2 \ell } \frac{\partial T}{\partial r} - \frac{h^{rr}}
	 * {\ell r}.
	 * \f]
	 */
	mutable Scalar* p_tilde_c ;


     



    // Constructors - Destructor :
    // -------------------------
	
    public:
	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param tipe  1-D array of integers (class \c Itbl ) of size 2 
	 *		containing the type 
	 *		of each index, \c COV  for a covariant one 
	 *		and \c CON  for a contravariant one,  with the 
	 *		following storage convention: 
	 *			\li \c tipe(0)  : type of the first index 
	 *			\li \c tipe(1)  : type of the second index 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 */
	Sym_tensor(const Map& map, const Itbl& tipe, const Base_vect& triad_i) ;

	/** Standard constructor when both indices are of the same type.
	 * 
	 * @param map   the mapping 
	 * @param tipe  the type of the indices.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 * 
	 */
	Sym_tensor(const Map& map, int tipe, const Base_vect& triad_i) ;

	Sym_tensor(const Sym_tensor& a) ; ///< Copy constructor

	/** Constructor from a \c Tensor .
	 *  The symmetry of the input tensor is assumed but is not checked.
	 */
	Sym_tensor(const Tensor& a) ;
	
	/** Constructor from a file (see \c sauve(FILE*) ).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param fich  file which has been used by 
	 *			    the function \c sauve(FILE*) .
	 */
	Sym_tensor(const Map& map, const Base_vect& triad_i, FILE* fich) ;

	virtual ~Sym_tensor() ;    ///< Destructor

      

    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const;	///< Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 

	/** Logical destructor of the derivatives depending on the i-th
	 *  element of \c met_depend  specific to the
	 *  class \c Sym_tensor  (\c p_transverse , etc...).
	 */	
	virtual void del_derive_met(int i) const ;

	/** Sets all the i-th components of \c met_depend  specific to the
	 * class \c Sym_tensor  (\c p_transverse , etc...) to 0x0.
	 */
	void set_der_met_0x0(int i) const ;


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another \c Sym_tensor 
	virtual void operator=(const Sym_tensor& a) ;

	/// Assignment to a \c Tensor_sym 
	virtual void operator=(const Tensor_sym& a) ;

	/**
	 * Assignment to a \c Tensor .
	 * 
	 * The symmetry is assumed but not checked.
	 */
	virtual void operator=(const Tensor& a) ;

	/**
	 * Assigns the derived members \c p_longit_pot and \c p_transverse
	 *  and updates the components accordingly.
	 * (see the documentation of these derived members for details)
	 */
	void set_longit_trans( const Vector& v, const Sym_tensor_trans& a) ;

	/** 
	 * Assigns the component \f$ T^{rr} \f$ and the derived members 
	 * \c p_eta , \c p_mu , \c p_www, \c p_xxx and \c p_ttt ,
	 * fro, their values and \f$ \eta / r\f$, \f$\mu / r \f$.
	 * It updates the other components accordingly.
	 */
	void set_auxiliary( const Scalar& trr, const Scalar& eta_over_r, const
			    Scalar& mu_over_r, const Scalar& www, const Scalar&
			    xxx, const Scalar& ttt ) ;

	/** Applies exponential filters to all components 
	 * (see \c Scalar::exponential_filter_r ). Does a loop for Cartesian 
	 * components, and works in terms of the rr-component, \f$\eta\f$,
	 * \f$\mu\f$, \c W, \c X, \c T for spherical components.
	 */
	virtual void  exponential_filter_r(int lzmin, int lzmax, int p, 
			    double alpha= -16.) ;

	/** Applies exponential filters to all components 
	 * (see \c Scalar::exponential_filter_ylm ). Does a loop for Cartesian 
	 * components, and works in terms of the r-component, \f$\eta\f$,
	 * \f$\mu\f$, \c W, \c X, \c T for spherical components. 
	 */
	virtual void exponential_filter_ylm(int lzmin, int lzmax, int p, 
			    double alpha= -16.) ;

    // Computation of derived members
    // ------------------------------
    public:


	/**Returns the divergence of \c this  with respect to a \c Metric .
	 * The indices are assumed to be contravariant.
	 */
	const Vector& divergence(const Metric&) const ; 

        /** Computes the Lie derivative of \c this  with respect to some
         *  vector field \c v 
         */
        Sym_tensor derive_lie(const Vector& v) const ; 

	/** Computes the transverse part \f${}^t T^{ij}\f$ of the tensor with respect 
	 * to a given metric, transverse meaning divergence-free with respect
	 * to that metric. Denoting \c *this  by \f$T^{ij}\f$, we then have
	 * \f[
	 *		T^{ij} = {}^t T^{ij} + \nabla^i W^j + \nabla^j W^i  
	 *		\qquad\mbox{with}\quad \nabla_j {}^t T^{ij} = 0 
	 *\f]
	 * where \f$\nabla_i\f$ denotes the covariant derivative with respect
	 * to the given metric and \f$W^i\f$ is the vector potential of the
	 * longitudinal part of \f$T^{ij}\f$ (function \c longit_pot()  below)
         * @param gam metric with respect to the transverse decomposition 
         *      is performed
         * @param par parameters for the vector Poisson equation
         * @param method_poisson type of method for solving the vector
         *      Poisson equation to get the longitudinal part (see 
         *      method \c Vector::poisson)
	 */
	const Sym_tensor_trans& transverse(const Metric& gam, Param* par = 0x0,
                int method_poisson = 6) const ; 

	/** Computes the vector potential \f$W^i\f$ of
	 * longitudinal part of the tensor (see documentation of
	 * method \c transverse() above).
         * @param gam metric with respect to the transverse decomposition 
         *      is performed
         * @param par parameters for the vector Poisson equation
         * @param method_poisson type of method for solving the vector
         *      Poisson equation to get the longitudinal part (see 
         *      method \c Vector::poisson)
	 */
	const Vector& longit_pot(const Metric& gam, Param* par = 0x0,
                int method_poisson = 6) const ; 
	
	/// Gives the field \f$\eta\f$ (see member \c p_eta ).
	virtual const Scalar& eta(Param* par = 0x0) const ;

	/// Gives the field \f$\mu\f$ (see member \c p_mu ).
	const Scalar& mu(Param* par = 0x0) const ;

	/// Gives the field \e W (see member \c p_www ).
	const Scalar& www() const ;

	/// Gives the field \e X (see member \c p_xxx ).
	const Scalar& xxx() const ;

	/// Gives the field \e T (see member \c p_ttt ).
	const Scalar& ttt() const ;

	/** Gives the field \e A (see member \c p_aaa ).
	 * @param output_ylm a flag to control the spectral decomposition 
	 * base of the result: if true (default) the spherical harmonics base 
	 * is used.
	 */
	const Scalar& compute_A(bool output_ylm = true, Param* par = 0x0) const ;

	/** Gives the field \f$\tilde{B}\f$ (see member \c p_tilde_b ).
	 * @param output_ylm a flag to control the spectral decomposition 
	 * base of the result: if true (default) the spherical harmonics base 
	 * is used.
	 */
	const Scalar& compute_tilde_B(bool output_ylm = true, Param* par = 0x0) const ;

	/** Gives the field \f$\tilde{B}\f$ (see member \c p_tilde_b )
	 * associated with the TT-part of the \c Sym_tensor .
	 * @param output_ylm a flag to control the spectral decomposition 
	 * base of the result: if true (default) the spherical harmonics base 
	 * is used.
	 */
	Scalar compute_tilde_B_tt(bool output_ylm = true, Param* par = 0x0) const ;

	/** Gives the field \f$\tilde{C}\f$ (see member \c p_tilde_c ).
	 * @param output_ylm a flag to control the spectral decomposition 
	 * base of the result: if true (default) the spherical harmonics base 
	 * is used.
	 */
	const Scalar& compute_tilde_C(bool output_ylm = true, Param* par = 0x0) const ;





 protected:
	/** Computes \f$\tilde{B}\f$ (see \c Sym_tensor::p_tilde_b ) from its
	 * transverse-traceless part and the trace.
	 */
	Scalar get_tilde_B_from_TT_trace(const Scalar& tilde_B_tt_in, const Scalar&
	    trace) const ;
	
    // Mathematical operators
    // ----------------------
 protected:
	/**
	 * Returns a pointer on the inverse of the \c Sym_tensor  
	 * (seen as a matrix).
	 */
	Sym_tensor* inverse() const ;

    // Friend classes
    //-----------------
	friend class Metric ;
 
} ;


			//---------------------------------//
			//    class Sym_tensor_trans       //
			//---------------------------------//
			

/**
 * Transverse symmetric tensors of rank 2. \ingroup (tensor)
 *
 * This class is designed to store transverse (divergence-free) 
 * symmetric contravariant tensors of rank 2,
 * with the component expressed in an orthonormal spherical basis
 * \f$(e_r,e_\theta,e_\varphi)\f$.
 *
 * 
 */
class Sym_tensor_trans: public Sym_tensor {

    // Data : 
    // -----
    protected:
	/// Metric with respect to which the divergence and the trace are defined
	const Metric* const met_div ; 
	
	/// Trace with respect to the metric \c *met_div  
	mutable Scalar* p_trace ; 
	
	/// Traceless part with respect to the metric \c *met_div  
	mutable Sym_tensor_tt* p_tt ;
	
    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *		    the tensor components are defined 
	 * @param met the metric with respect to which the divergence is defined
	 */
	Sym_tensor_trans(const Map& map, const Base_vect& triad_i, 
		const Metric& met) ;

	Sym_tensor_trans(const Sym_tensor_trans& ) ;       ///< Copy constructor

	/** Constructor from a file (see \c Tensor::sauve(FILE*) ).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param met the metric with respect to which the divergence is defined
	 * @param fich  file which has been used by 
	 *			    the function \c sauve(FILE*) .
	 */
	Sym_tensor_trans(const Map& map, const Base_vect& triad_i, 
		const Metric& met, FILE* fich) ;

	virtual ~Sym_tensor_trans() ;			///< Destructor

 
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const;	///< Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 


    // Accessors
    // ---------
        public:
	/** Returns the metric with respect to which the divergence 
	 *  and the trace are defined.
	 */
	const Metric& get_met_div() const {return *met_div ; } ;

    // Mutators / assignment
    // ---------------------

	public:
	/// Assignment to another \c Sym_tensor_trans 
	virtual void operator=(const Sym_tensor_trans& a) ;	
	
	/// Assignment to a \c Sym_tensor 
	virtual void operator=(const Sym_tensor& a) ;	
	
	/// Assignment to a \c Tensor_sym 
	virtual void operator=(const Tensor_sym& a) ;

	/// Assignment to a \c Tensor 
	virtual void operator=(const Tensor& a) ;	
	
	/**
	 * Assigns the derived members \c p_tt and \c p_trace
	 *  and updates the components accordingly.
	 * (see the documentation of these derived members for details)
	 */
	void set_tt_trace(const Sym_tensor_tt& a, const Scalar& h, 
			  Param* par = 0x0) ;

	// Computational methods
	// ---------------------
	/// Returns the trace of the tensor with respect to metric \c *met_div 
	const Scalar& the_trace() const ; 
	
	/** Returns the transverse traceless part of the tensor, 
	 * the trace being defined
	 * with respect to metric \c *met_div 
	 */
	const Sym_tensor_tt& tt_part(Param* par = 0x0) const ; 

 protected:
	/** Solves a system of two coupled first-order PDEs obtained from 
	 * the divergence-free condition (Dirac gauge) and the requirement that
	 * the potential \e A (see \c Sym_tensor::p_aaa ) has a given value.
	 * The system reads: \f{eqnarray*}
	 * \frac{\partial \tilde{\mu}}{\partial r}  + \frac{3\tilde{\mu}}{r} + \left( 
	 * \Delta_{\theta\varphi } + 2\right) X &=& 0;\\
	 * \frac{\partial X}{\partial r} - \frac{\tilde{\mu}}{r} &=& A. \f}
	 * Note that this is solved only for \f$\ell \geq 2\f$ and that 
	 * \f$\tilde{\mu} = \mu / r\f$ (see \c Sym_tensor::p_mu ).
	 *
	 * @param aaa [input] the source \e A
	 * @param tilde_mu [output] the solution \f$\tilde{\mu}\f$
	 * @param xxx [output] the solution \e X
	 * @param par_bc [input] \c Param to control the boundary conditions
	 */
	void sol_Dirac_A(const Scalar& aaa, Scalar& tilde_mu, Scalar& xxx,
			 const Param* par_bc = 0x0) const ;


	/** Solves a system of three coupled first-order PDEs obtained from 
	 * divergence-free conditions (Dirac gauge) and the requirement that
	 * the potential \f$\tilde{B}\f$ (see \c Sym_tensor::p_tilde_b ) has 
	 * a given value. The system reads: \f{eqnarray*}
	 * \frac{\partial T^{rr}}{r} + \frac{3T^{rr}}{r} +\frac{1}{r}
	 *  \Delta_{\theta\varphi } \tilde{\eta} &=& \frac{h}{r};\\
	 * \frac{\partial \tilde{\eta}}{\partial r} + \frac{3\tilde{\eta}}{r} -
	 * \frac{T^{rr}}{2r} + \left( \Delta_{\theta\varphi } + 2\right) 
	 * \frac{W}{r} &=& -\frac{h}{2r};\\
	 * (\ell + 2) \frac{\partial W}{\partial r} + \ell(\ell + 2)
	 * \frac{W}{r} - \frac{2\tilde{\eta}}{r} + \frac{(\ell +2)T}{2r(\ell + 1)}
	 * + \frac{1}{2(\ell + 1)} \frac{\partial T}{\partial r} - \frac{T^{rr}}
	 * {(\ell + 1)r} &=& \tilde{B} - \frac{1}{2(\ell +1)} \frac{\partial h}
	 * {\partial r} - \frac{\ell +2}{\ell +1} \frac{h}{2r}.\f}
	 * Note that \f$\tilde{\eta} = \eta / r\f$ (for definitions, see derived
	 * members of \c Sym_tensor).
	 *
	 * @param tilde_b [input] the source \f$\tilde{B}\f$
	 * @param hh [input] the trace of the tensor
	 * @param hrr [output] the \e rr component of the result
	 * @param tilde_eta [output] the solution \f$\tilde{\eta}\f$
	 * @param www [output] the solution \e W
	 * @param par_bc [input] \c Param to control the boundary conditions
	 * @param par_mat [input/output] \c Param in which the operator matrix is
	 *                stored.
	 */
	void sol_Dirac_tilde_B(const Scalar& tilde_b, const Scalar& hh, Scalar& hrr,
			       Scalar& tilde_eta, Scalar& www, Param* par_bc=0x0,
			       Param* par_mat=0x0) const ;

	/** Solves the same system as \c Sym_tensor_trans::sol_Dirac_tilde_B
	 * but only for \f$\ell=0,1\f$. In these particular cases, \e W =0
	 * the system is simpler and homogeneous solutions are different.
	 */
	void sol_Dirac_l01(const Scalar& hh, Scalar& hrr, Scalar& tilde_eta,
			   Param* par_mat) const ;



 public:


        
	/** Same resolution as sol_Dirac_A, but with inner boundary conditions added. 
	 *For now, only Robyn-type boundary conditions on \f$\frac {\mu}  {r} \f$ can be imposed.
	 */


	void sol_Dirac_Abound(const Scalar& aaa, Scalar& tilde_mu, Scalar& x_new,
		 	  Scalar bound_mu, const Param* par_bc);
 
        
	/** Same resolution as sol_Dirac_Abound, but here the boundary conditions
	 * are the degenerate elliptic conditions encontered when solving the
	 * Kerr problem.
	 */


	void sol_Dirac_A2(const Scalar& aaa, Scalar& tilde_mu, Scalar& x_new,
		 	  Scalar bound_mu, const Param* par_bc);       

	/** Same resolution as sol_Dirac_tilde_B, but with inner boundary conditions added.
         *  The difference is here, one has to put B and C values in (and not only  \f$\tilde{B}\f$).
         * For now, only Robyn-type boundary conditions on \f$ h^{rr} \f$ can be imposed.
	 */

 	void sol_Dirac_BC2(const Scalar& bb, const Scalar& cc, const Scalar& hh, 
	 			Scalar& hrr, Scalar& tilde_eta, Scalar& ww, Scalar bound_eta,double dir, double neum, double rhor, Param* par_bc, Param* par_mat); 
 

	/** Same resolution as sol_Dirac_Abound, but here the boundary conditions
	 * are the degenerate elliptic conditions encontered when solving the
	 * Kerr problem.
	 */

 	void sol_Dirac_BC3(const Scalar& bb, const Scalar& hh, 
	 			Scalar& hrr, Scalar& tilde_eta, Scalar& ww, Scalar bound_hrr, Scalar bound_eta, Param* par_bc, Param* par_mat); 



	 // Solving the electric system for l=0 and l=1 only (simpler case), with boundary conditions imposed by the degenerate elliptic system.
   
 	void sol_Dirac_l01_bound(const Scalar& hh, Scalar& hrr, Scalar& tilde_eta, Scalar& bound_hrr, Scalar& bound_eta, Param* par_mat) ;

	// Provisory: just for compilation, to be removed
 	void sol_Dirac_l01_2(const Scalar& hh, Scalar& hrr, Scalar& tilde_eta, Param* par_mat) ;


	/** Assigns the derived member \c p_tt and computes the trace so that 
	 * \c *this + the flat metric has a determinant equal to 1. It then
	 * updates the components accordingly, with a \c dzpuis = 2. 
	 * This function makes an 
	 * iteration until the relative difference in the trace between 
	 * two steps is lower than \c precis . 
	 *
	 * @param htt the transverse traceless part; all components must have
	 *            dzpuis = 2.
	 * @param precis relative difference in the trace computation to end
	 *               the iteration.
	 * @param it_max maximal number of iterations.
	 */
	void trace_from_det_one(const Sym_tensor_tt& htt, 
				double precis = 1.e-14, int it_max = 100) ;

	/** Assigns the \e rr component and the derived member \f$\mu\f$.
	 * Other derived members are deduced from the divergence-free 
	 * condition. Finally, it computes \c T (\c Sym_tensor::p_ttt )  so that 
	 * \c *this + the flat metric has a determinant equal to 1. It then
	 * updates the components accordingly. This function makes an 
	 * iteration until the relative difference in \c T between 
	 * two steps is lower than \c precis . 
	 *
	 * @param hrr the \e rr component of the tensor,
	 * @param mu_in the \f$\mu\f$ potential,
	 * @param precis relative difference in the trace computation to end
	 *               the iteration.
	 * @param it_max maximal number of iterations.
	 */
	void set_hrr_mu_det_one(const Scalar& hrr, const Scalar& mu_in,
				double precis = 1.e-14, int it_max = 100) ;

	/** Assignes the TT-part of the tensor.
	 * The trace is deduced from the divergence-free condition, through the 
	 * Dirac system on \f$ \tilde{B} \f$, so that 
	 * \c *this + the flat metric has a determinant equal to 1. It then
	 * updates the components accordingly. This function makes an 
	 * iteration until the relative difference in the trace between 
	 * two steps is lower than \c precis . 
	 * @param hijtt  the TT part for \c this.
	 * @param h_prev a pointer on a guess for the trace of \c *this; if
	 *               null, then the iteration starts from 0.
	 * @param precis relative difference in the trace computation to end
	 *               the iteration.
	 * @param it_max maximal number of iterations.
	 */
	void set_tt_part_det_one(const Sym_tensor_tt& hijtt, const 
				 Scalar* h_prev = 0x0, Param* par_mat = 0x0, 
				 double precis = 1.e-14, int it_max = 100) ;

	/** Assigns the derived member \c A and computes \f$\tilde{B}\f$ 
	 * from its TT-part (see \c Sym_tensor::compute_tilde_B_tt() ).
	 * Other derived members are deduced from the divergence-free 
	 * condition. Finally, it computes the trace so that 
	 * \c *this + the flat metric has a determinant equal to 1. It then
	 * updates the components accordingly. This function makes an 
	 * iteration until the relative difference in the trace between 
	 * two steps is lower than \c precis . 
	 *
	 * @param a_in the \c A potential (see \c Sym_tensor::p_aaa )
	 * @param tbtt_in the TT-part of \f$\tilde{B}\f$ potential 
	 *   (see \c Sym_tensor::p_tilde_b and \c Sym_tensor::compute_tilde_B_tt() )
	 * @param h_prev a pointer on a guess for the trace of \c *this; if
	 *               null, then the iteration starts from 0.
	 * @param precis relative difference in the trace computation to end
	 *               the iteration.
	 * @param it_max maximal number of iterations.
	 */
	void set_AtBtt_det_one(const Scalar& a_in, const Scalar& tbtt_in, 
			       const Scalar* h_prev = 0x0, Param* par_bc = 0x0,
			       Param* par_mat = 0x0, double precis = 1.e-14, 
			       int it_max = 100) ;

	/** Assigns the derived members \c A , \f$\tilde{B}\f$ and the trace.
	 * Other derived members are deduced from the divergence-free condition.
	 *
	 * @param a_in the \c A potential (see \c Sym_tensor::p_aaa )
	 * @param tb_in the \f$\tilde{B}\f$ potential (see \c Sym_tensor::p_tilde_b )
	 * @param trace the trace of the \c Sym_tensor.
	 */
	void set_AtB_trace(const Scalar& a_in, const Scalar& tb_in, const 
			   Scalar& trace, Param* par_bc = 0x0, Param* par_mat = 0x0) ;

	/** Computes the solution of a tensorial transverse Poisson equation
	 *  with \c *this  \f$= S^{ij}\f$ as a source:
	 * \f[
	 *    \Delta h^{ij} = S^{ij}.
	 *\f] 
	 * In particular, it makes an iteration on the trace of the result, using
	 * \c Sym_tensor::set_WX_det_one.
	 * 
	 * @param h_guess a pointer on a guess for the trace of the result; it is
	 *                passed to \c Sym_tensor::set_WX_det_one.
	 * @return solution \f$h^{ij}\f$ of the above equation with the boundary
	 *	condition \f$h^{ij}=0\f$ at spatial infinity.
	 */
	Sym_tensor_trans poisson(const Scalar* h_guess = 0x0) const ; 
} ; 
	

			//------------------------------//
			//    class Sym_tensor_tt       //
			//------------------------------//
			

/**
 * Transverse and traceless symmetric tensors of rank 2.
 *
 * This class is designed to store transverse (divergence-free) 
 * and transverse symmetric contravariant tensors of rank 2,
 * with the component expressed in an orthonormal spherical basis
 * \f$(e_r,e_\theta,e_\varphi)\f$.\ingroup (tensor)
 *
 * 
 */
class Sym_tensor_tt: public Sym_tensor_trans {

    // Data : 
    // -----

    protected:
	/** Field \f$\chi\f$ such that the component \f$h^{rr} = \frac{\chi}{r^2}\f$.
	 */
	mutable Scalar* p_khi ;
	
	
    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *		    the tensor components are defined 
	 * @param met the metric with respect to which the divergence is defined
	 */
	Sym_tensor_tt(const Map& map, const Base_vect& triad_i, 
		const Metric& met) ;

	Sym_tensor_tt(const Sym_tensor_tt& ) ;       ///< Copy constructor

	/** Constructor from a file (see \c Tensor::sauve(FILE*) ).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param met the metric with respect to which the divergence is defined
	 * @param fich  file which has been used by 
	 *			    the function \c sauve(FILE*) .
	 */
	Sym_tensor_tt(const Map& map, const Base_vect& triad_i, 
		const Metric& met, FILE* fich) ;

	virtual ~Sym_tensor_tt() ;			///< Destructor

 
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const;	///< Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------

	public:
	/// Assignment to another \c Sym_tensor_tt 
	virtual void operator=(const Sym_tensor_tt& a) ;	
	
	/// Assignment to a \c Sym_tensor_trans 
	virtual void operator=(const Sym_tensor_trans& a) ;	
	
	/// Assignment to a \c Sym_tensor 
	virtual void operator=(const Sym_tensor& a) ;	
	
	/// Assignment to a \c Tensor_sym 
	virtual void operator=(const Tensor_sym& a) ;

	/// Assignment to a \c Tensor 
	virtual void operator=(const Tensor& a) ;	
	
	/** Sets the component \f$h^{rr}\f$, as well as the angular potentials 
	 * \f$\eta\f$ and \f$\mu\f$ (see members
	 *  \c p_eta  and \c p_mu ). 
	 *  The other components are updated consistently
	 *  by a call to the method \c update() .
	 *
	 *	@param hrr [input] value of \f$h^{rr}\f$
	 *	@param eta_i [input] angular potential \f$\eta\f$
	 *	@param mu_i [input] angular potential \f$\mu\f$
	 *
	 */
	void set_rr_eta_mu(const Scalar& hrr, const Scalar& eta_i, 
						const Scalar& mu_i) ; 

	/** Sets the component \f$h^{rr}\f$, as well as the angular potential
	 * \f$\mu\f$ (see member \c p_mu ). 
	 * The angular potential \f$\eta\f$ (member \c p_eta ) is deduced from
	 * the divergence free condition. 
	 * The other tensor components are updated consistently
	 * by a call to the method \c update() .
	 *
	 *	@param hrr [input] value of \f$h^{rr}\f$
	 *	@param mu_i [input] angular potential \f$\mu\f$
	 *
	 */
	void set_rr_mu(const Scalar& hrr, const Scalar& mu_i) ; 
	
	
	/** Sets the component \f$\chi\f$, as well as the angular potentials 
	 * \f$\eta\f$ and \f$\mu\f$ (see members \c p_khi ,
	 *  \c p_eta  and \c p_mu ). 
	 *  The components are updated consistently
	 *  by a call to the method \c update() .
	 *
	 *	@param khi_i [input] value of \f$\chi\f$
	 *	@param eta_i [input] angular potential \f$\eta\f$
	 *	@param mu_i [input] angular potential \f$\mu\f$
	 *
	 */
	void set_khi_eta_mu(const Scalar& khi_i, const Scalar& eta_i, 
						const Scalar& mu_i) ; 
		
	/** Sets the component \f$\chi\f$, as well as the angular potential
	 * \f$\mu\f$ (see member \c p_khi  and \c p_mu ). 
	 * The angular potential \f$\eta\f$ (member \c p_eta ) is deduced from
	 * the divergence free condition. 
	 * The tensor components are updated consistently
	 * by a call to the method \c update() .
	 *
	 *	@param khi_i [input] value of \f$\chi\f$
	 *	@param mu_i [input] angular potential \f$\mu\f$
         *      @param dzp [input] \c dzpuis parameter of the resulting
         *                      tensor components
	 *
	 */
	void set_khi_mu(const Scalar& khi_i, const Scalar& mu_i, int dzp = 0,
			Param* par1 = 0x0, Param* par2 = 0x0, 
			Param* par3 = 0x0) ; 

	/** Assigns the derived members \c A and \f$\tilde{B}\f$.
	 * Other derived members are deduced from the divergence-and trace-free 
	 * conditions.
	 *
	 * @param a_in the \c A potential (see \c Sym_tensor::p_aaa )
	 * @param tb_in the \f$\tilde{B}\f$ potential (see \c Sym_tensor::p_tilde_b )
	 */
	void set_A_tildeB(const Scalar& a_in, const Scalar& tb_in, Param* par_bc = 0x0,
			  Param* par_mat = 0x0) ;

	// Computational methods
	// ---------------------
	
	public:
	/** Gives the field \f$\chi\f$ such that the component 
	 * \f$h^{rr} = \frac{\chi}{r^2}\f$.
	 */
	const Scalar& khi() const ;
	
	/// Gives the field \f$\eta\f$ (see member \c p_eta ).
	virtual const Scalar& eta(Param* par = 0x0) const ;

	protected:
	/** Computes the components \f$h^{r\theta}\f$, \f$h^{r\varphi}\f$,
	 * \f$h^{\theta\theta}\f$, \f$h^{\theta\varphi}\f$ and \f$h^{\varphi\varphi}\f$,
	 *  from \f$h^{rr}\f$ and the potentials \f$\eta\f$ and \f$\mu\f$.
         *  @param dzp \c dzpuis parameter of the result, i.e. of the 
         *      components \f$ h^{ij} \f$.
	 */
	void update(int dzp, Param* par1 = 0x0, Param* par2 = 0x0) ;

	public:
	/** Computes the solution of a tensorial TT Poisson equation
	 *  with \c *this  \f$= S^{ij}\f$ as a source:
	 * \f[
	 *    \Delta h^{ij} = S^{ij}
	 *\f] 
	 * 
	 * @param dzfin [input] the \c dzpuis for all the components of the result 
	 *        (see the documentation for \c Scalar ).
	 * @return solution \f$h^{ij}\f$ of the above equation with the boundary
	 *	condition \f$h^{ij}=0\f$ at spatial infinity.
	 */
	Sym_tensor_tt poisson(int dzfin = 2) const ; 
	


} ; 
	



}
#endif
