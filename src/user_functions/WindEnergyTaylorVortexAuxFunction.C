/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/WindEnergyTaylorVortexAuxFunction.h>
#include <NaluEnv.h>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace sierra{
namespace nalu{

WindEnergyTaylorVortexAuxFunction::WindEnergyTaylorVortexAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> &params) :
  AuxFunction(beginPos, endPos),
  centroidX_(-2.5),
  centroidY_(0.0),
  rVortex_(0.25),
  beta_(15.0),
  uInf_(10.0),
  density_(1.0e-3),
  visc_(1.0e-4),
  Lx_(5.)
{
 //  check size and populate
  if ( params.size() != 5 && !params.empty() )
    throw std::runtime_error("Realm::setup_initial_conditions: wind_energy_taylor_vortex takes five parameters:"
        " centroidX, centroidY, initial vortex radius, utheta0, uInf.");

  if (!params.empty()) {

    centroidX_ = params[0];
    centroidY_ = params[1];
    rVortex_ = params[2];
    beta_ = params[3];
    uInf_ = params[4];
  }
  else {
    NaluEnv::self().naluOutputP0()
        << "wind_energy_taylor_vortex proceeding with default parameters."
        << "\n  centroidX: " << centroidX_
        << "\n  centroidY: " << centroidY_
        << "\n  rVortex: " << rVortex_
        << "\n  beta: " << beta_
        << "\n  uinf: " << uInf_
        << "\n density: " << density_
        << "\n visc: " << visc_
        << std::endl;
  }
}

void
WindEnergyTaylorVortexAuxFunction::do_evaluate(
  const double *coords,
  const double time,
  const unsigned /*spatialDimension*/,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  double xCentroid = centroidX_ + uInf_ * time;
  double travelDist = Lx_ - centroidX_;
  xCentroid = (xCentroid > travelDist) ? xCentroid - Lx_ : xCentroid;

  const double tHat = (visc_ > 1.0e-30) ? time * beta_ / rVortex_ : 0.0; // zero out viscous terms if mu ~ 0
  const double Re = (visc_ > 1.0e-30) ? density_ * beta_ * rVortex_ / visc_ : 0.0; // zero out viscous terms if mu ~ 0
  const double tFac = 1.0 + 2.0 * tHat / Re;
  for(unsigned p=0; p < numPoints; ++p) {
    const double cX = coords[0];
    const double cY = coords[1];
    const double cDx = cX - (centroidX_ + uInf_ * time);
    const double cDy = cY - centroidY_;

    const double radius = std::sqrt(cDx*cDx+cDy*cDy)/(rVortex_);
    const double factor = beta_ / (rVortex_ * tFac * tFac) * std::exp(0.5*(1.0 - radius * radius / tFac));

    const double velX = -factor * cDy + uInf_;
    const double velY =  factor * cDx;

    fieldPtr[0] = velX;
    fieldPtr[1] = velY;
    fieldPtr += fieldSize;
    coords += fieldSize;
  }
}

} // namespace nalu
} // namespace Sierra
