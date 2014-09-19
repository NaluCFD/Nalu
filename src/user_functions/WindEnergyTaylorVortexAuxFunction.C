/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/WindEnergyTaylorVortexAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

WindEnergyTaylorVortexAuxFunction::WindEnergyTaylorVortexAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> &params) :
  AuxFunction(beginPos, endPos),
  centroidX_(0.0),
  centroidY_(0.0),
  rVortex_(0.0),
  beta_(0.0),
  uInf_(0.0)
{
  // check size and populate
  if ( params.size() != 5 )
    throw std::runtime_error("Realm::setup_initial_conditions: wind_energy_taylor_vortex requires 5 params: ");
  centroidX_ = params[0];
  centroidY_ = params[1];
  rVortex_   = params[2];
  beta_      = params[3];
  uInf_      = params[4];
}


void
WindEnergyTaylorVortexAuxFunction::do_evaluate(
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

    const double cX = coords[0];
    const double cY = coords[1];

    const double cDx = cX - centroidX_;
    const double cDy = cY - centroidY_;

    const double radius = std::sqrt(cDx*cDx+cDy*cDy)/rVortex_;
    const double factor = beta_/rVortex_*std::exp(0.5*(1.0-radius*radius));

    const double velX = uInf_-factor*cDy;
    const double velY = factor*cDx;

    fieldPtr[0] = velX;
    fieldPtr[1] = velY;
    
    fieldPtr += fieldSize;
    coords += fieldSize;
  }
}

} // namespace nalu
} // namespace Sierra
