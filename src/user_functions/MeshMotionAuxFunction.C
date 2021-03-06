/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/MeshMotionAuxFunction.h>
#include <MeshMotionInfo.h>
#include <PecletFunction.h>
#include <Realm.h>
#include <SolutionOptions.h>

// basic c++
#include <algorithm>
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

MeshMotionAuxFunction::MeshMotionAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<std::string> theStringParams,
  Realm &realm) :
  AuxFunction(beginPos, endPos),
  omegaMM_(NULL),
  centroidMM_(NULL),
  velMM_(NULL)
{
  // extract mesh info 
  if (theStringParams.size() < 1 )
    throw std::runtime_error("Mesh user function requires at least one string parameter");
  const std::string meshMotionBlockName = theStringParams[0];

  // find the mesh motion info object
  std::map<std::string, MeshMotionInfo *>::const_iterator iterFind;
  MeshMotionInfo *meshInfo = NULL;
  iterFind = realm.solutionOptions_->meshMotionInfoMap_.find(meshMotionBlockName);
  if ( iterFind != realm.solutionOptions_->meshMotionInfoMap_.end() )
    meshInfo = iterFind->second;
  else 
    throw std::runtime_error("MeshMotionAuxFunction::error() Can not find mesh motion block name " + meshMotionBlockName);       
  
  // set the data
  omegaMM_ = &meshInfo->bodyOmega_;
  centroidMM_ = &meshInfo->centroid_;
  dispMM_ = &meshInfo->bodyDispCC_;
  velMM_ = &meshInfo->bodyVel_;
}

MeshMotionAuxFunction::~MeshMotionAuxFunction()
{

  // Nothing required 

}

void
MeshMotionAuxFunction::setup(const double time)
{

  // Nothing required

}

void
MeshMotionAuxFunction::do_evaluate(
  const double *coords,
  const double /*time*/,
  const unsigned /*spatialDimension*/,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  double c[3] = {0.0,0.0,0.0};
  double u[3] = {0.0,0.0,0.0};
  for(unsigned p=0; p < numPoints; ++p) {

    // define distance from centroid and compute cross product
    for ( unsigned i = 0; i < fieldSize; ++i ) 
      c[i] = coords[i] - ((*centroidMM_)[i]+(*dispMM_)[i]);  
 
    cross_product(c, u);

    for ( unsigned i = 0; i < fieldSize; ++i ) 
      u[i] += (*velMM_)[i];
    
    for ( unsigned i = 0; i < fieldSize; ++i )
      fieldPtr[i] = u[i];
  
    fieldPtr += fieldSize;
    coords += fieldSize;
  }
}

void
MeshMotionAuxFunction::cross_product(double *c, double *u) const
{
  u[0] = (*omegaMM_)[1]*c[2] - (*omegaMM_)[2]*c[1];
  u[1] = (*omegaMM_)[2]*c[0] - (*omegaMM_)[0]*c[2];
  u[2] = (*omegaMM_)[0]*c[1] - (*omegaMM_)[1]*c[0];
}

} // namespace nalu
} // namespace Sierra
