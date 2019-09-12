/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/PulseVelocityAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

PulseVelocityAuxFunction::PulseVelocityAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> &params) :
  AuxFunction(beginPos, endPos),
  ux_(-10.0),
  uy_(1.0),
  uz_(0.01),
  pulseTime_(1.0)
{
  if ( params.size() != 4 || params.empty() )
    throw std::runtime_error("PulseVelocityAuxFunction: requires u, v, w, time");
  
  // extract them
  ux_ = params[0];
  uy_ = params[1];
  uz_ = params[2];
  pulseTime_ = params[3];
}

void
PulseVelocityAuxFunction::do_evaluate(
  const double *coords,
  const double time,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  const double fac = time > pulseTime_ ? 0.0 : 1.0;
  for(unsigned p=0; p < numPoints; ++p) {
        
    fieldPtr[0] = ux_*fac;
    fieldPtr[1] = uy_*fac;
    fieldPtr[2] = uz_*fac;

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
