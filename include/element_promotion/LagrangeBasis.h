/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef LagrangeBasis_h
#define LagrangeBasis_h

#include <vector>

namespace sierra{
namespace nalu{

class Lagrange1D
{
public:
  Lagrange1D(const double* nodeLocs, int order);

  Lagrange1D(std::vector<double> nodeLocs);

  Lagrange1D(int order);

  virtual ~Lagrange1D();

  double interpolation_weight(double x, unsigned nodeNumber) const;

  double derivative_weight(double x, unsigned nodeNumber) const;

private:
  void set_lagrange_weights();
  std::vector<double> lagrangeWeights_;
  std::vector<double> nodeLocs_;
};

class LagrangeBasis
{
public:
  LagrangeBasis(
    const std::vector<std::vector<int>>&  indicesMap,
    const std::vector<double>& nodeLocs
  );

  virtual ~LagrangeBasis();

  std::vector<double> eval_basis_weights(
    const std::vector<double>& intgLoc) const;

  std::vector<double> eval_deriv_weights(
    const std::vector<double>& intgLoc) const;

  double tensor_lagrange_derivative(
    unsigned dimension,
    const double* x,
    const int* node_ordinals,
    unsigned derivativeDirection
  ) const;

  double tensor_lagrange_interpolant(unsigned dimension, const double* x, const int* node_ordinals) const;

  std::vector<std::vector<int>> indicesMap_;
  const Lagrange1D basis1D_;
  unsigned numNodes1D_;
  unsigned polyOrder_;
  const unsigned dim_;
};


} // namespace nalu
} // namespace Sierra

#endif
