/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/ConvectingTaylorVortexPressureAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

ConvectingTaylorVortexPressureAuxFunction::ConvectingTaylorVortexPressureAuxFunction() :
  AuxFunction(0,1),
  uNot_(1.0),
  vNot_(1.0),
  pNot_(1.0),
  visc_(0.001),
  pi_(std::acos(-1.0))
{
  // nothing to do
}


void
ConvectingTaylorVortexPressureAuxFunction::do_evaluate(
  const double *coords,
  const double t,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  for(unsigned p=0; p < numPoints; ++p) {

    double x = coords[0];
    double y = coords[1];
    const double omega = pi_*pi_*visc_;    

    fieldPtr[0] =  -(pNot_/4.)*(cos(2.*pi_*(x-uNot_*t)) + cos(2.*pi_*(y-vNot_*t)))*exp(-4.0*omega*t);

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
