#ifndef _LINEAR_H_
#define _LINEAR_H_

#include"cParticles.h"

#ifdef CHOLMOD
 #include <Eigen/CholmodSupport>
#else
// #include <Eigen/IterativeLinearSolvers>
 #include <Eigen/SparseCholesky>
//#include <Eigen/SparseQR>
#endif

#include <iostream>
#include <vector>

#include <unsupported/Eigen/SparseExtra>

using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
//using Eigen::ConjugateGradient;

//const Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "];");

class linear {
 public:  // constructor
 linear(Triangulation& TT) : T(TT) {}

  void fill_matrices();

  void solve_for_weights( );
  void ps_equation(const FT dt );
  void p_equation(const FT dt );
  void s_equation(const FT dt );
  void u_star( void );
  void reset_p( void );
  void reset_s( void );
  void u_add_grad_ps( const FT dt ) ;
  void u_add_grad_p( const FT dt );
  void u_add_grad_s( const FT dt ) ;
  
  void DD_scalar_vfield(const vfield_list::take from , const sfield_list::take to );  
  VectorXd DD_scalar_vfield(const vfield_list::take from );
  void DD_times_sfield(const sfield_list::take from ,
		       VectorXd& Dx,VectorXd& Dy);
  void MM_times_sfield(const sfield_list::take from ,
		       VectorXd& Dx,VectorXd& Dy);
  void copy(const sfield_list::take from, sfield_list::take to  );


  
private:

  
  Triangulation& T; // Its triangulation

  typedef   SparseMatrix<FT>  SpMat;
  typedef Eigen::Triplet<FT> triplet;

  SpMat Dx, Dy;
  SpMat Mx, My;

  SpMat DD, DM, MD, MM;
  SpMat DDMM;

  VectorXd field_to_vctr(const sfield_list::take sf );
  void vctr_to_field(const VectorXd& vv, const sfield_list::take sf  );
  void vfield_to_vctrs(const vfield_list::take vf , VectorXd& vx, VectorXd& vy ) ;
  void vctrs_to_vfield(const VectorXd& vx, const VectorXd& vy , const vfield_list::take vf ) ;

#define DIRECT_SOLVER

#ifdef DIRECT_SOLVER
#ifdef CHOLMOD
  Eigen::CholmodSupernodalLLT<SpMat> DDMM_solver;
#else
  Eigen::SimplicialLDLT<SpMat> DDMM_solver;
  Eigen::SimplicialLDLT<SpMat> DD_solver;
  Eigen::SimplicialLDLT<SpMat> MM_solver;
#endif
#else
  Eigen::BiCGSTAB<SpMat> solver_stiffp1;
#endif


};


#endif
