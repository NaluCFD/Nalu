/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef PowerlawPipeVelocityAuxFunction_h
#define PowerlawPipeVelocityAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class PowerlawPipeVelocityAuxFunction : public AuxFunction
{
public:

  PowerlawPipeVelocityAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    const std::vector<double> &params);

  virtual ~PowerlawPipeVelocityAuxFunction() {}
  
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
  double uzCenterline_;
  double xCentroid_;
  double yCentroid_;
  double R_;
  double powerExponent_;
};

} // namespace nalu
} // namespace Sierra

#endif
