/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/WindEnergyTaylorVortexPressureAuxFunction.h>
#include <NaluEnv.h>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace sierra{
namespace nalu{

WindEnergyTaylorVortexPressureAuxFunction::WindEnergyTaylorVortexPressureAuxFunction(
  const std::vector<double> &params) :
  AuxFunction(0,1),
  centroidX_(-2.5),
  centroidY_(0.0),
  rVortex_(0.25),
  beta_(15.0),
  uInf_(10.0),
  density_(1.0e-3)
{
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
        << std::endl;
  }
}

void
WindEnergyTaylorVortexPressureAuxFunction::do_evaluate(
  const double *coords,
  const double /*time*/,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  // no time update.  Just used for IC
  const double coeff = -0.5 * std::exp(1.) * beta_ * beta_ * density_;
  for(unsigned p=0; p < numPoints; ++p) {
    const double cX = coords[0];
    const double cY = coords[1];
    const double cDx = cX - centroidX_;
    const double cDy = cY - centroidY_;
    const double radius = std::sqrt(cDx*cDx+cDy*cDy)/rVortex_;
    *fieldPtr = coeff * std::exp(-radius * radius);
    ++fieldPtr;
    coords += spatialDimension;
  }
}

WindEnergyTaylorVortexPressureGradAuxFunction::WindEnergyTaylorVortexPressureGradAuxFunction(
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
  visc_(1.0e-4)
{
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
WindEnergyTaylorVortexPressureGradAuxFunction::do_evaluate(
  const double *coords,
  const double time,
  const unsigned /*spatialDimension*/,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  const double tHat = time * beta_ / rVortex_;
  const double Re = density_ * beta_ * rVortex_ / visc_;
  const double tFac = 1.0 + 2.0 * tHat / Re;

  for(unsigned p=0; p < numPoints; ++p) {

    const double cX = coords[0];
    const double cY = coords[1];

    const double cDx = cX - (centroidX_ + uInf_ * time);
    const double cDy = cY - centroidY_;

    const double radius = std::sqrt(cDx*cDx+cDy*cDy)/rVortex_;
    const double utheta = beta_/(rVortex_ * tFac * tFac) * std::exp(0.5*(1.0-radius*radius/tFac));
    double gradFac = density_ * utheta * utheta;
    double theta = std::atan2(cDy, cDx);

    fieldPtr[0] = gradFac * radius * rVortex_ * std::cos(theta);
    fieldPtr[1] = gradFac * radius * rVortex_ * std::sin(theta);

    fieldPtr += fieldSize;
    coords += fieldSize;
  }
}

} // namespace nalu
} // namespace Sierra
