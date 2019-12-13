/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/ChannelFlowPerturbedPlugVelocityAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

ChannelFlowPerturbedPlugVelocityAuxFunction::ChannelFlowPerturbedPlugVelocityAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> &params) :
  AuxFunction(beginPos, endPos),
    u_m(1.0),		     // Max velocity
    Lx_(2.0*acos(-1.0)),     // Domain length in x-direction
    Ly_(2.0),                // Domain length in y-direction
    Lz_(2.0/3.0*acos(-1.0)), // Domain length in z-direction
    streamwiseD_(0),
    spanwiseD_(1),
    normalD_(2),
    pi_(acos(-1.0))
{
  if ( params.size() > 0 && params[0] > 0.0 )
    u_m = params[0];
  if ( params.size() > 1 && params[1] > 0.0 )
    Lx_ = params[1];
  if ( params.size() > 2 && params[2] > 0.0 )
    Ly_ = params[2];
  if ( params.size() > 3 && params[3] > 0.0 )
    Lz_ = params[3];
  if ( params.size() > 4 )
    streamwiseD_ = (int)params[4];
  if ( params.size() > 5 )
    spanwiseD_ = (int)params[5];
  if ( params.size() > 6 )
    normalD_ = (int)params[6];
}

void
ChannelFlowPerturbedPlugVelocityAuxFunction::do_evaluate(
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

    const double x = coords[streamwiseD_];
    const double y = coords[spanwiseD_];
    const double z = coords[normalD_];

    const double ubase = u_m*std::tanh(8.0*y)*std::tanh(8.0*(Ly_-y));
    fieldPtr[streamwiseD_] = ubase + 0.05*ubase*sin(8.0*pi_/Lx_*x)*cos(4.0*pi_/Lz_*z);
    fieldPtr[spanwiseD_] = 0.05*ubase*cos(8.0*pi_/Lx_*x)*sin(4.0*pi_/Lz_*z);
    fieldPtr[normalD_] = 0.05*ubase*sin(4.0*pi_/Lx_*x)*cos(2.0*pi_/Lz_*z);

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
