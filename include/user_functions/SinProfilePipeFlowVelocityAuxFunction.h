/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SinProfilePipeFlowVelocityAuxFunction_h
#define SinProfilePipeFlowVelocityAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class SinProfilePipeFlowVelocityAuxFunction : public AuxFunction
{
public:

  SinProfilePipeFlowVelocityAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    const std::vector<double> &params);

  virtual ~SinProfilePipeFlowVelocityAuxFunction() {}
  
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
  double umag_;
  double a_;
  double pertx_;
  double perty_;
  double pertz_;
  const double pi_;
};

} // namespace nalu
} // namespace Sierra

#endif
