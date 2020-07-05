/*
 *  Definition of Lorene classes Ope_elementary
 *
 */

/*
 *   Copyright (c) 2003 Philippe Grandclement
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

#ifndef __OPE_ELEMENTARY_H_ 
#define __OPE_ELEMENTARY_H_ 

/*
 * $Id: ope_elementary.h,v 1.12 2014/10/13 08:52:36 j_novak Exp $
 * $Log: ope_elementary.h,v $
 * Revision 1.12  2014/10/13 08:52:36  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.11  2007/05/06 10:48:08  p_grandclement
 * Modification of a few operators for the vorton project
 *
 * Revision 1.10  2007/04/24 09:04:11  p_grandclement
 * Addition of an operator for the vortons
 *
 * Revision 1.9  2004/12/23 16:30:14  j_novak
 * New files and class for the solution of the rr component of the tensor Poisson
 * equation.
 *
 * Revision 1.8  2004/08/24 09:14:40  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.7  2004/06/22 08:49:57  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * Revision 1.6  2004/06/14 15:07:10  j_novak
 * New methods for the construction of the elliptic operator appearing in
 * the vector Poisson equation (acting on eta).
 *
 * Revision 1.5  2004/05/10 15:28:21  j_novak
 * First version of functions for the solution of the r-component of the
 * vector Poisson equation.
 *
 * Revision 1.4  2004/03/23 14:54:45  j_novak
 * More documentation
 *
 * Revision 1.3  2004/03/09 09:15:11  p_grandclement
 * Correction pour le raccord...
 *
 * Revision 1.2  2004/03/05 09:18:48  p_grandclement
 * Addition of operator sec_order_r2
 *
 * Revision 1.1  2003/12/11 14:57:00  p_grandclement
 * I had forgotten the .h (sorry folks...)
 *
 * Revision 1.2  2001/12/11 06:44:41  e_gourgoulhon
 *
 * $Header: /cvsroot/Lorene/C++/Include/ope_elementary.h,v 1.12 2014/10/13 08:52:36 j_novak Exp $
 *
 */

#include "matrice.h"

namespace Lorene {
/**
 * Basic class for elementary elliptic operators.
 *
 * Such objects describe a type of elliptic operator, in a given domain and 
 * for a given spherical harmonics.
 * They are called by the general elliptic solver. \c Ope_elementary 
 * objects 
 * know how to compute the approriate particular solution for a given source 
 * and the associated homogeneous ones.
 * 
 * The class \c Ope_elementary is an abstract one: 
 * it cannot be instanciated. 
 * Specific implementation of coordinate mappings will be performed by derived
 * classes of \c Ope_elementary . \ingroup (ellip)
 * 
 **/

class Ope_elementary {

 protected:

  int nr ; ///< Number of radial points
  int base_r ; ///< Radial basis of decomposition
  double alpha ; ///< Parameter \f$\alpha\f$ of the associated mapping.
  double beta ; ///< Parameter \f$\beta\f$ of the associated mapping.
  
  /**
   * Pointer on the matrix representation of the operator.
   **/
  mutable Matrice* ope_mat ;
  /**
   * Pointer on the banded-matrix of the operator.
   **/
  mutable Matrice* ope_cl ;
  /**
   * Pointer on the non-degenerated matrix of the operator.
   **/
  mutable Matrice* non_dege ;

  /**
   * Value of the first homogeneous solution at the outer boundary.
   **/
  mutable double s_one_plus ;
  /**
   * Value of the first homogeneous solution at the inner boundary.
   **/
  mutable double s_one_minus ;
  /**
   * Value of the derivative  of the 
   * first homogeneous solution at the outer boundary.
   **/
  mutable double ds_one_plus ;
  /**
   * Value of the derivative  of the 
   * first homogeneous solution at the inner boundary.
   **/
  mutable double ds_one_minus ;
  
  /**
   * Value of the second homogeneous solution at the outer boundary.
   **/
  mutable double s_two_plus ;
  /**
   * Value of the second homogeneous solution at the inner boundary.
   **/
  mutable double s_two_minus ;
  /**
   * Value of the derivative  of the 
   * second homogeneous solution at the outer boundary.
   **/
  mutable double ds_two_plus ;
  /**
   * Value of the derivative  of the 
   * second homogeneous solution at the inner boundary.
   **/
  mutable double ds_two_minus ;
 
  /**
   * Value of the particular solution at the inner boundary.
   **/
  mutable double sp_minus ;
  /**
   * Value of the particular solution at the outer boundary.
   **/
  mutable double sp_plus ;
  /**
   * Value of the derivative of the particular solution at the inner boundary.
   **/
  mutable double dsp_minus ;
  /**
   * Value of the derivative of the particular solution at the outer boundary.
   **/
  mutable double dsp_plus ;

  // Constructors destructor
 protected:
  /**
   * Standard constructor, protected because the class is an abstract one.
   *
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   **/
  explicit Ope_elementary (int nbr , int baser , double alf, double eta) ;
  Ope_elementary (const Ope_elementary&) ; ///< Constructor by copy

 public:
  virtual ~Ope_elementary() ; ///< Destructor

 public: 
  /**
   * Returns the value of the first homogeneous solution at the inner
   * boundary.
   **/
  double val_sh_one_minus() const {return s_one_minus ;} ; 
  /**
   * Returns the value of the first homogeneous solution at the outer
   * boundary.
   **/
  double val_sh_one_plus() const {return s_one_plus ;} ; 
  /**
   * Returns the value of the derivative of the 
   * first homogeneous solution at the inner
   * boundary.
   **/
  double der_sh_one_minus() const {return ds_one_minus ;} ; 
  /**
   * Returns the value of the derivative of the 
   * first homogeneous solution at the outer
   * boundary.
   **/
  double der_sh_one_plus() const {return ds_one_plus ;} ;
 
  /**
   * Returns the value of the second homogeneous solution at the inner
   * boundary.
   **/
  double val_sh_two_minus() const {return s_two_minus ;} ;
  /**
   * Returns the value of the second homogeneous solution at the outer
   * boundary.
   **/
  double val_sh_two_plus() const {return s_two_plus ;} ;
  /**
   * Returns the value of the derivative of the 
   * second homogeneous solution at the inner
   * boundary.
   **/
  double der_sh_two_minus() const {return ds_two_minus ;} ;
  /**
   * Returns the value of the derivative of the 
   * second homogeneous solution at the outer
   * boundary.
   **/
  double der_sh_two_plus() const {return ds_two_plus ;} ;

  /**
   * Returns the value of the particular solution at the inner boundary.
   **/
  double val_sp_minus() const {return sp_minus ;} ;
  /**
   * Returns the value of the particular solution at the outer boundary.
   **/
  double val_sp_plus() const {return sp_plus ;} ; 
  /**
   * Returns the value of the derivative 
   * particular solution at the inner boundary.
   **/
  double der_sp_minus() const {return dsp_minus ;} ;
  /**
   * Returns the value of the derivative 
   * particular solution at the outer boundary.
   **/
  double der_sp_plus() const {return dsp_plus ;} ;

  /// Returns \c alpha .
  double get_alpha() const {return alpha ;} ; 

  /// Returns \c beta}.
  double get_beta() const {return beta ;} ; 

  /// Returns \c base_r}. 
  int get_base_r() const {return base_r ;} ;

  /// Returns the matrix representation.
   Matrice get_ope_mat() {
	if (ope_mat ==0x0)
		do_ope_mat() ;
	return *ope_mat ;
}

  /// Returns the banded matrix representation.
   Matrice get_ope_cl() {
	if (ope_cl ==0x0)
		do_ope_cl() ;
	return *ope_cl ;
} 
 
  /// Returns the non degenerate matrix representation.
   Matrice get_non_dege() {
	if (non_dege ==0x0)
		do_non_dege() ;
	return *non_dege ;
}
 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const = 0 ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const = 0 ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const = 0 ;  
 
 public:
  /**
   * Computes the particular solution, given the source \c so .
   **/
  virtual Tbl get_solp(const Tbl& so) const = 0 ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const = 0 ;
  /**
   * Increases the quatum number \e l by one unit.
   **/
  virtual void inc_l_quant() = 0 ;
} ;

/**
 * Class for the operator of the Poisson equation (i.e. the Laplacian !).
 *
 * It is implemented in every type of domains.
 **/

class Ope_poisson : public Ope_elementary {

 protected:
  int l_quant ; ///< quantum number
  int dzpuis ; ///< the associated dzpuis, if in the compactified domain. 
  
 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   * @param lq [input] quantum number \e l .
   * @param dz [input] dzpuis of the source.
   **/
  Ope_poisson (int nbr, int baser, double alf, double bet, int lq, int dz) ;
  Ope_poisson (const Ope_poisson&) ; ///< Constructor by copy
  virtual ~Ope_poisson() ; ///< Destructor

  /// Returns the associated dzpuis, if in the compactified domain.
  int get_dzpuis() {return dzpuis ;} ;

  /// Returns the quantum number \e l
  int get_lquant() {return l_quant;} ;

 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the particular solution, given the source \c so .
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number \e l by one unit.
   **/
  virtual void inc_l_quant() ;
  /**
   * Decreases the quatum number \e l by one unit.
   **/
  virtual void dec_l_quant() ;
} ;

/**
 * Class for the Helmholtz operator \f$\Delta - m^2\f$ (\f$m > 0\f$).
 * 
 * It is implemented only in the shells and in the compactified domain.
 **/
class Ope_helmholtz_minus : public Ope_elementary {

 protected: 
  int lq ; ///< The quantum number \e l
  double masse ; ///< The mass parameter \e m .
 
 public:
   /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param lq [input] the quatum number \e l.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   * @param mas [input] mass parameter \e m .
   **/
  Ope_helmholtz_minus (int nbr, int baser, int lq, double alf, double bet, 
		       double mas) ;
  Ope_helmholtz_minus (const Ope_helmholtz_minus&) ; ///< Constructor by copy
  virtual ~Ope_helmholtz_minus() ; ///< Destructor
  
 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the particular solution, given the source \c so .
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number \e l by one unit (CURRENTLY NOT IMPLEMENTED)
   **/
  virtual void inc_l_quant() ;
} ;

/**
 * Class for the Helmholtz operator \f$\Delta + m^2\f$ (\e m > 0).
 * 
 * It is implemented only in the shells.
 **/
class Ope_helmholtz_plus : public Ope_elementary {

 protected:
  int lq ; ///< The quantum number \e l
  double masse ; ///< The mass parameter \e m .

 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.  
   * @param lq [input] the quatum number \e l.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   * @param mas [input] mass parameter \e m .
   **/
  Ope_helmholtz_plus (int nbr, int baser, int lq, double alf, double bet, 
		      double mas) ;
  Ope_helmholtz_plus (const Ope_helmholtz_plus&) ; ///< Constructor by copy
  virtual ~Ope_helmholtz_plus() ; ///< Destructor
  
 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the particular solution, given the source \c so .
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number \e l  by one unit (CURRENTLY NOT IMPLEMENTED)
   **/
  virtual void inc_l_quant() ;
} ;

/**
 * Class for operator of the type 
 * \f$ a r^2 \partial^2 / \partial r^2 + b r \partial / \partial r + c\f$.
 * 
 * It is implemented only in the shells.
 **/
class Ope_sec_order_r2 : public Ope_elementary {

 protected:

  double a_param ; ///< The parameter \e a .
  double b_param ; ///< The parameter \e b .
  double c_param ; ///< The parameter \e c .

 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   * @param a [input] parameter \e a  .
   * @param b [input] parameter \e b .    
   * @param c [input] parameter \e c .
   **/

  Ope_sec_order_r2 (int nbr, int baser, double alf, double bet, 
		      double a, double b, double c) ;

  Ope_sec_order_r2 (const Ope_sec_order_r2&) ; ///< Constructor by copy
  virtual ~Ope_sec_order_r2() ; ///< Destructor
  
 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the particular solution, given the source \c so .
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number \e l  by one unit (CURRENTLY NOT IMPLEMENTED)
   **/
  virtual void inc_l_quant() ;
} ;

/**
 * Class for operator of the type 
 * \f$ a \partial^2 / \partial r^2 + b \partial / \partial r + c\f$.
 * 
 * It is implemented only in the shells.
 **/
class Ope_sec_order : public Ope_elementary {

 protected:

  double a_param ; ///< The parameter \f$a\f$.
  double b_param ; ///< The parameter \f$b\f$.
  double c_param ; ///< The parameter \f$c\f$.

 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   * @param a [input] parameter \f$a\f$ .
   * @param b [input] parameter \f$b\f$ .    
   * @param c [input] parameter \f$c\f$ .
   **/

  Ope_sec_order (int nbr, int baser, double alf, double bet, 
		      double a, double b, double c) ;

  Ope_sec_order (const Ope_sec_order&) ; ///< Constructor by copy
  virtual ~Ope_sec_order() ; ///< Destructor
  
 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:

  /**
   * Computes the particular solution, given the source \c so .
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number \f$l\f$ by one unit (CURRENTLY NOT IMPLEMENTED)
   **/
  virtual void inc_l_quant() ;
} ;

/**
 * Class for the operator of the \e r component of the vector
 * Poisson equation. The operator reads \f$\Delta + \frac{2}{r} 
 *\frac{\partial}{\partial r} 
 * + \frac{2}{r^2} \f$ in all domains, for \f$ l \not= 0 \f$; and to 
 * \f$\frac{\partial^2}{\partial r^2} + \frac{2}{r} \frac{\partial}
 * {\partial r} - \frac{2}{r^2} \f$ 
 * in all domains otherwise.
 *
 * It is implemented in every type of domain.
 **/
class Ope_pois_vect_r : public Ope_poisson {

 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   * @param lq [input] quantum number \e l .
   * @param dz [input] dzpuis of the source.
   **/
  Ope_pois_vect_r (int nbr, int baser, double alf, double bet, int lq, int dz) ;
  Ope_pois_vect_r (const Ope_pois_vect_r&) ; ///< Constructor by copy
  virtual ~Ope_pois_vect_r() ; ///< Destructor

 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
} ;

/**
 * Class for the operator of the \e rr component of the divergence-free
 * tensor Poisson equation. The operator reads \f$\Delta + \frac{4}{r} 
 *\frac{\partial}{\partial r} 
 * + \frac{6}{r^2} \f$ in all domains, for \f$l \geq 2\f$.
 *
 * It is implemented in every type of domain.
 **/
class Ope_pois_tens_rr : public Ope_poisson {

 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   * @param lq [input] quantum number \e l .
   * @param dz [input] dzpuis of the source.
   **/
  Ope_pois_tens_rr (int nbr, int baser, double alf, double bet, int lq, int dz) ;
  Ope_pois_tens_rr (const Ope_pois_tens_rr&) ; ///< Constructor by copy
  virtual ~Ope_pois_tens_rr() ; ///< Destructor

 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
} ;

/**
 * Class for the operator of the Poisson equation in 2D.
 *
 * It is implemented in every type of domains.
 **/

class Ope_poisson_2d : public Ope_elementary {

 protected:
  int l_quant ; ///< quantum number
  int dzpuis ; ///< the associated dzpuis, if in the compactified domain. 
  
 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   * @param lq [input] quantum number \e l .
   * @param dz [input] dzpuis of the source.
   **/
  Ope_poisson_2d (int nbr, int baser, double alf, double bet, int lq, int dz) ;
  Ope_poisson_2d (const Ope_poisson_2d&) ; ///< Constructor by copy
  virtual ~Ope_poisson_2d() ; ///< Destructor

  /// Returns the associated dzpuis, if in the compactified domain.
  int get_dzpuis() {return dzpuis ;} ;

  /// Returns the quantum number \e l
  int get_lquant() {return l_quant;} ;

 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the particular solution, given the source \c so .
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number \e l by one unit.
   **/
  virtual void inc_l_quant() ;
  /**
   * Decreases the quatum number \e l by one unit.
   **/
  virtual void dec_l_quant() ;

} ;

/**
 * Class for the operator of the Helmholtz equation in 2D.
 *
 **/

class Ope_helmholtz_minus_2d : public Ope_elementary {

 protected:
  int l_quant ; ///< quantum number
  double masse ; ///< The mass term.
  int dzpuis ; ///< the associated dzpuis, if in the compactified domain. 
  
 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   * @param lq [input] quantum number \e l .
   * @param masse [input] mass term.
   * @param dz [input] dzpuis of the source.
   **/
  Ope_helmholtz_minus_2d (int nbr, int baser, double alf, double bet, int lq, double masse, int dz) ;
  Ope_helmholtz_minus_2d (const Ope_helmholtz_minus_2d&) ; ///< Constructor by copy
  virtual ~Ope_helmholtz_minus_2d() ; ///< Destructor

  /// Returns the associated dzpuis, if in the compactified domain.
  int get_dzpuis() {return dzpuis ;} ;

  /// Returns the quantum number \e l
  int get_lquant() {return l_quant;} ;

  /// Returns the mass term
  double get_masse() {return masse;} ;

 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the particular solution, given the source \c so .
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number \e l by one unit.
   **/
  virtual void inc_l_quant() ;
  /**
   * Decreases the quatum number \e l by one unit.
   **/
  virtual void dec_l_quant() ;
} ;

/**
 * Class for the operator of the Poisson equation in pseudo 1d.
 *
 **/

class Ope_poisson_pseudo_1d : public Ope_elementary {

 protected:
  int l_quant ; ///< quantum number
    
 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   * @param lq [input] quantum number \e l .
   **/
  Ope_poisson_pseudo_1d (int nbr, int baser, double alf, double bet, int lq) ;
  Ope_poisson_pseudo_1d (const Ope_poisson_pseudo_1d&) ; ///< Constructor by copy
  virtual ~Ope_poisson_pseudo_1d() ; ///< Destructor

  /// Returns the quantum number \e l
  int get_lquant() {return l_quant;} ;

 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the particular solution, given the source \c so .
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number \e l by one unit.
   **/
  virtual void inc_l_quant() ;
  /**
   * Decreases the quatum number \e l by one unit.
   **/
  virtual void dec_l_quant() ;

} ;


/**
 * Class for the operator of the modified Helmholtz equation in pseudo-1d.
 *
 * It is implemented only in the external domain
 **/

class Ope_helmholtz_minus_pseudo_1d : public Ope_elementary {

 protected:
  int l_quant ; ///< quantum number
  double masse ; ///< The mass term.
  int dzpuis ; ///< the associated dzpuis, if in the compactified domain. 
  
 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   * @param lq [input] quantum number \e l .
   * @param masse [input] mass term.
   * @param dz [input] dzpuis of the source.
   **/
  Ope_helmholtz_minus_pseudo_1d (int nbr, int baser, double alf, double bet, 
				 int lq, double masse, int dz) ;
  Ope_helmholtz_minus_pseudo_1d (const Ope_helmholtz_minus_pseudo_1d&) ; ///< Constructor by copy
  virtual ~Ope_helmholtz_minus_pseudo_1d() ; ///< Destructor

  /// Returns the associated dzpuis, if in the compactified domain.
  int get_dzpuis() {return dzpuis ;} ;

  /// Returns the quantum number \e l
  int get_lquant() {return l_quant;} ;

  /// Returns the mass term
  double get_masse() {return masse;} ;

 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the particular solution, given the source \c so .
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number \e l by one unit.
   **/
  virtual void inc_l_quant() ;
  /**
   * Decreases the quatum number \e l by one unit.
   **/
  virtual void dec_l_quant() ;
} ;

/**
 * Class for the operator appearing for the vortons
 *
 * It is implemented in the shells and the compactified domain
 **/

class Ope_vorton : public Ope_elementary {

 protected:
  int l_quant ; ///< quantum number
  int dzpuis ; ///< the associated dzpuis, if in the compactified domain. 
  
 public:
  /**
   * Standard constructor.
   * 
   * @param nbr [input] number of radial points.
   * @param baser [input] radial basis of decomposition.
   * @param alf [input] parameter \f$\alpha\f$ of the mapping.
   * @param bet [input] parameter \f$\beta\f$ of the mapping.
   * @param lq [input] quantum number \e l .
   * @param dz [input] dzpuis of the source.
   **/
  Ope_vorton (int nbr, int baser, double alf, double bet, int lq, int dz) ;
  Ope_vorton (const Ope_vorton&) ; ///< Constructor by copy
  virtual ~Ope_vorton() ; ///< Destructor

  /// Returns the associated dzpuis, if in the compactified domain.
  int get_dzpuis() {return dzpuis ;} ;

  /// Returns the quantum number \e l
  int get_lquant() {return l_quant;} ;

 private:
  /**
   * Computes the matrix of the operator.
   **/
  virtual void do_ope_mat() const ;
  /**
   * Computes the banded-matrix of the operator.
   **/
  virtual void do_ope_cl() const ;
  /**
   * Computes the non-degenerated matrix of the operator.
   **/
  virtual void do_non_dege() const ;  
  
 public:
  /**
   * Computes the particular solution, given the source \c so .
   **/
  virtual Tbl get_solp(const Tbl& so) const ;
  /**
   * Computes the homogeneous solutions(s).
   **/
  virtual Tbl get_solh() const ;
  /**
   * Increases the quatum number \e l by one unit.
   **/
  virtual void inc_l_quant() ;
  /**
   * Decreases the quatum number \e l by one unit.
   **/
  virtual void dec_l_quant() ;
} ;

}
#endif

