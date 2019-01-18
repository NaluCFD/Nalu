/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef PowerlawVelocityAuxFunction_h
#define PowerlawVelocityAuxFunction_h

#include "AuxFunction.h"

#include <vector>

namespace sierra{
namespace nalu{

class PowerlawVelocityAuxFunction : public AuxFunction
{
public:

  PowerlawVelocityAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    const std::vector<double> &params);

  virtual ~PowerlawVelocityAuxFunction() {}
  
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
  double uKnown_;
  double yKnown_;
  double alpha_;
};

} // namespace nalu
} // namespace Sierra

#endif
