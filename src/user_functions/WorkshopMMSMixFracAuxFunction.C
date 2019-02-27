/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/WorkshopMMSMixFracAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

WorkshopMMSMixFracAuxFunction::WorkshopMMSMixFracAuxFunction() :
  AuxFunction(0,1),
  u_(1.0),
  v_(0.0),
  sigma_(0.5)
{
  // does nothing
}

void
WorkshopMMSMixFracAuxFunction::do_evaluate(
  const double *coords,
  const double time,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  for(unsigned p=0; p < numPoints; ++p) {
    const double x = coords[0];
    const double y = coords[1];
    const double xtu = x-time*u_;
    const double ytv = y-time*v_;
    fieldPtr[0] = std::exp(-(xtu*xtu + ytv*ytv)/sigma_/sigma_);
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
