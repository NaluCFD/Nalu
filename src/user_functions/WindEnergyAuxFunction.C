/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/WindEnergyAuxFunction.h>
#include <algorithm>
#include <PecletFunction.h>
#include <Realm.h>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

WindEnergyAuxFunction::WindEnergyAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> theParams,
  Realm &realm) :
  AuxFunction(beginPos, endPos),
  centroidX_(0.0),
  centroidY_(0.0),
  omega_(1.0),
  omegaBlend_(1.0),
  tanhFunction_(NULL)
{
  // nothing; note hard coded for 2D
  if (theParams.size() < 1 )
    throw std::runtime_error("Wind_energy user function requires at least one parameter");
  omega_ = theParams[0];
  
  // check for possible centroids
  if ( theParams.size() > 1)
    centroidX_ = theParams[1];
  if ( theParams.size() > 2)
    centroidY_ = theParams[2];
  
  // check for omega blending
  const std::string omegaName = "omega";
  if ( realm.get_tanh_functional_form(omegaName) == "tanh") {
    const double c1 = realm.get_tanh_trans(omegaName);
    const double c2 = realm.get_tanh_width(omegaName);
    tanhFunction_ = new TanhFunction(c1, c2);
  }
}

WindEnergyAuxFunction::~WindEnergyAuxFunction()
{
  if ( NULL != tanhFunction_ )
    delete tanhFunction_;
}

void
WindEnergyAuxFunction::setup(const double time)
{
  if ( NULL != tanhFunction_  ) {
    omegaBlend_ = tanhFunction_->execute(time);
  }
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

    double cX = coords[0] - centroidX_;
    double cY = coords[1] - centroidY_;

    fieldPtr[0] = -omega_*cY*omegaBlend_;
    fieldPtr[1] = +omega_*cX*omegaBlend_;
    
    fieldPtr += fieldSize;
    coords += fieldSize;
  }
}

} // namespace nalu
} // namespace Sierra
