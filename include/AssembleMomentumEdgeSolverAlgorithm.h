/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleMomentumEdgeSolverAlgorithm_h
#define AssembleMomentumEdgeSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;
template <typename T> class PecletFunction;

class AssembleMomentumEdgeSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleMomentumEdgeSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleMomentumEdgeSolverAlgorithm();
  virtual void initialize_connectivity();
  virtual void execute();
  
  double van_leer(
    const double &dqm,
    const double &dqp,
    const double &small);

  const bool meshMotion_;
  const double includeDivU_;

  VectorFieldType *velocityRTM_;
  VectorFieldType *velocity_;
  VectorFieldType *coordinates_;
  GenericFieldType *dudx_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  ScalarFieldType *viscosity_lam_;
  VectorFieldType *edgeAreaVec_;
  ScalarFieldType *massFlowRate_;
  ScalarFieldType *dualNodalVolume_;

  // peclet function specifics
  PecletFunction<double>* pecletFunction_;

  // Coefficient of kinetic energy perturbation from user
  double coeff_kk_;
  
  // Target eigenvalue and relative distance perturbation from user
  double BinvXt_[3];
  double deltaB_;
  
  // Eigenvector basis perturbation from user
  int eigenvectorPermutation_;
  
  // momentum perturbation function members
  void diagonalize( const double (&A)[3][3], double (&Q)[3][3], double (&D)[3][3]);
  void sort(double (&Q)[3][3], double (&D)[3][3]);
  void perturb(double (&Q)[3][3], double (&D)[3][3]);
  void form_perturbed_stress( const double (&D)[3][3], const double (&Q)[3][3], double (&A)[3][3]);
  void matrix_matrix_multiply( const double (&A)[3][3], const double (&B)[3][3], double (&C)[3][3]);
  
  // momentum perturbation data members
  double S_[3][3];
  double A_[3][3];
  double D_[3][3];
  double Q_[3][3];
  double Q0_[3][3];
  double resolved_lambda_[3];
  double filter;
  double nu_sgs;
  double total_kk;
  double resolved_kk;
  double modeled_kk;
  const double Cw_;

};

} // namespace nalu
} // namespace Sierra

#endif
