/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/SinProfilePipeFlowVelocityAuxFunction.h>
#include <algorithm>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

SinProfilePipeFlowVelocityAuxFunction::SinProfilePipeFlowVelocityAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> &params) :
  AuxFunction(beginPos, endPos),
  umag_(-10.0),
  a_(1.0),
  pertx_(0.01),
  perty_(-0.01),
  pertz_(1.0),
  pi_(acos(-1.0))
{
  if ( params.size() != 5 || params.empty() )
    throw std::runtime_error("SinProfilePipeFlowVelocityAuxFunction: requires umag, a, pertx, perty, pertz");
  
  // extract them
  umag_ = params[0];
  a_ = params[1];
  pertx_ = params[2];
  perty_ = params[3]; 
  pertz_ = params[4]; 
}

void
SinProfilePipeFlowVelocityAuxFunction::do_evaluate(
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
    
    const double z = coords[2];    
    const double uz = umag_*sin(pi_*a_*z);
    
    fieldPtr[0] = uz*pertx_;
    fieldPtr[1] = uz*perty_;
    fieldPtr[2] = uz*pertz_;

    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

} // namespace nalu
} // namespace Sierra
