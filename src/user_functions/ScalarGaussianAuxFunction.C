/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/ScalarGaussianAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

ScalarGaussianAuxFunction::ScalarGaussianAuxFunction(
  const std::vector<double> &params) :
  AuxFunction(0, 1),
  coordDirection_(0),
  a_(1.0),
  b_(0.0),
  c_(0.2)
{
  // extract params
  if ( params.size() != 4 || params.empty() )
    throw std::runtime_error("ScalarGaussianAuxFunction: requires coordinate direction, scale, l, a, and sigma");
  
  // extract parameters
  coordDirection_ = params[0];
  a_ = params[1]; // height
  b_ = params[2]; // center of peak
  c_ = params[3]; // standard deviation
}


void
ScalarGaussianAuxFunction::do_evaluate(
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
    
    //const double theValue = std::sin(pi_*coords[0]);
    const double bracket = (coords[coordDirection_] - b_);

    fieldPtr[0] = a_*std::exp(-0.5*bracket*bracket/(c_*c_));
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
