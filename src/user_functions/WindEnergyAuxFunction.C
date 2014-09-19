/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/WindEnergyAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

WindEnergyAuxFunction::WindEnergyAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> theParams) :
  AuxFunction(beginPos, endPos),
  omega_(1.0)
{
  // nothing; note hard coded for 2D
  if (theParams.size() != 1)
    throw std::runtime_error("Wind_energy user function requires one parameter");
  omega_ = theParams[0];
}


void
WindEnergyAuxFunction::do_evaluate(
  const double *coords,
  const double /*time*/,
  const unsigned /*spatialDimension*/,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  for(unsigned p=0; p < numPoints; ++p) {

    double cX = coords[0];
    double cY = coords[1];

    fieldPtr[0] = -omega_*cY;
    fieldPtr[1] = +omega_*cX;
    
    fieldPtr += fieldSize;
    coords += fieldSize;
  }
}

} // namespace nalu
} // namespace Sierra
