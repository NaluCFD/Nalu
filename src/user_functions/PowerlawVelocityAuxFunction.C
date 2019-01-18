/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "user_functions/PowerlawVelocityAuxFunction.h"

// basic c++
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace sierra{
namespace nalu{

PowerlawVelocityAuxFunction::PowerlawVelocityAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> &params) :
  AuxFunction(beginPos, endPos),
  uKnown_(10.0),
  yKnown_(5.0),
  alpha_(1.0/7.0)
{
  // extract required params
  if (params.size() != 3) {
    throw std::runtime_error("PowerlawVelocityAuxFunction requires three parameters");
  }
  else {
    uKnown_ = params[0];
    yKnown_ = params[1];
    alpha_ = 1.0/params[2];
  }
}

void
PowerlawVelocityAuxFunction::do_evaluate(
  const double *coords,
  const double /*time*/,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  for(unsigned p=0; p < numPoints; ++p) {

    // ux = f(y)    
    double cY = std::max(1.0e-16, coords[1]);
    
    // standard power law 
    const double uX = uKnown_*std::pow(cY/yKnown_, alpha_);

    fieldPtr[0] = uX;
    fieldPtr[1] = 0.0;
    fieldPtr[2] = 0.0;
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
