/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/SinVofAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

SinVofAuxFunction::SinVofAuxFunction(
  const std::vector<double> &/*params*/) :
  AuxFunction(0, 1),
  coordDirection_(0),
  heightDirection_(1),
  lambda_(1.0e-4),
  aNot_(0.1*lambda_),
  pi_(std::acos(-1.0))
{
  // extract params
  /*
    if ( params.size() != 2 || params.empty() )
    throw std::runtime_error("SinVofAuxFunction: requires coordinate direction and lambda");
    
    // extract parameters
    coordDirection_ = params[0];
    coordDirection_ = params[0];
    a_ = params[1]; 
  */
}


void
SinVofAuxFunction::do_evaluate(
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
    const double height = aNot_*std::cos(2.0*pi_*coords[coordDirection_]/lambda_);
    fieldPtr[0] = 1.0;
    if ( coords[heightDirection_] >= height )
      fieldPtr[0] = 0.0;
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
