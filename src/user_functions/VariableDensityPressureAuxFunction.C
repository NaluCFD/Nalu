/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/VariableDensityPressureAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

VariableDensityPressureAuxFunction::VariableDensityPressureAuxFunction() :
  AuxFunction(0,1),
  pnot_(1.0),
  a_(20.0),
  pi_(acos(-1.0))
{
  // does nothing
}

void
VariableDensityPressureAuxFunction::do_evaluate(
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

    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];
    
    fieldPtr[0] = -pnot_/4.0*( cos(2.0*a_*pi_*x) + 
                               cos(2.0*a_*pi_*y) + 
                               cos(2.0*a_*pi_*z));
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
