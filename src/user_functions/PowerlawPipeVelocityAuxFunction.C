/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/PowerlawPipeVelocityAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

// remove
#include <iostream>
#include <ostream>

namespace sierra{
namespace nalu{

PowerlawPipeVelocityAuxFunction::PowerlawPipeVelocityAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> &params) :
  AuxFunction(beginPos, endPos),
  uzCenterline_(1.0),
  xCentroid_(0.0),
  yCentroid_(0.0),
  R_(0.005),
  powerExponent_(1.0)
{
  if ( params.size() != 5 || params.empty() )
    throw std::runtime_error("PowerlawPipeVelocityAuxFunction: requires uzCenterline_, xCentroid_, yCentroid_, R_ and n");
  
  // extract them
  uzCenterline_ = params[0];
  xCentroid_ = params[1];
  yCentroid_ = params[2];
  R_ = params[3];
  powerExponent_ = 1.0/params[4];
}

void
PowerlawPipeVelocityAuxFunction::do_evaluate(
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

    const double x = coords[0] - xCentroid_;
    const double y = coords[1] - yCentroid_;

    const double r = std::sqrt(x*x + y*y);
    const double powArg = std::max(0.0, 1.0 - r/R_);
    
    fieldPtr[0] = 0.0;
    fieldPtr[1] = 0.0;
    fieldPtr[2] = uzCenterline_*std::pow(powArg, powerExponent_);
    
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
