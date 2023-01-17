/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef ScalarGaussianAuxFunction_h
#define ScalarGaussianAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class ScalarGaussianAuxFunction : public AuxFunction
{
public:

  ScalarGaussianAuxFunction(
    const std::vector<double> &params);

  virtual ~ScalarGaussianAuxFunction() {}
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;
  
private:
  int coordDirection_;
  double a_;
  double b_;
  double c_;
};

} // namespace nalu
} // namespace Sierra

#endif
