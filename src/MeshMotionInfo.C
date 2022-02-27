/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <MeshMotionInfo.h>

// standard c++
#include <string>
#include <vector>
#include <array>

namespace sierra{
namespace nalu{


//==========================================================================
// Class Definition
//==========================================================================
// MeshMotionInfo - holder for user options for mesh motion
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
// Rotation constructor
MeshMotionInfo::MeshMotionInfo(
  std::vector<std::string> meshMotionBlock, 
  const double omega, 
  std::vector<double> centroid,
  std::vector<double> unitVec,
  const bool computeCentroid,
  const double theAngle)
  : meshMotionBlock_(meshMotionBlock), 
    omega_(omega), 
    centroid_(centroid),
    unitVec_(unitVec),
    computeCentroid_(computeCentroid),
    computeCentroidCompleted_(false),
    theAngle_(theAngle),
    sixDof_(false),
    bodyDispCC_(std::vector<double>(3,0.0)),
    bodyVel_(std::vector<double>(3,0.0)),
    bodyAccel_(std::vector<double>(3,0.0)),
    bodyAlpha_(std::vector<double>(3,0.0)),
    bodyMass_(0.0),
    bodyDen_(0.0)
{

}

// 6-DOF constructor
MeshMotionInfo::MeshMotionInfo(
    std::vector<std::string> meshMotionBlock,
    std::vector<std::string> forceSurface,
    std::vector<double> bodyDispCC,
    std::vector<double> bodyAngle,
    std::vector<double> bodyOmega,
    std::vector<double> bodyPrincInertia,
    std::vector<double> centroid,
    std::vector<double> bodyVel,
    const double bodyMass,
    const double bodyDen,
    std::vector<double> appliedForce, 
    const bool computeCentroid,
    const std::vector<std::array<double, 9>> tetherGeom)
    : meshMotionBlock_(meshMotionBlock),
      omega_(0.0),
      centroid_(centroid),
      computeCentroid_(computeCentroid),
      computeCentroidCompleted_(false),
      theAngle_(0.0),
      sixDof_(true),
      bodyDispCC_(bodyDispCC),
      bodyAngle_(bodyAngle),
      bodyOmega_(bodyOmega),
      bodyPrincInertia_(bodyPrincInertia),
      bodyVel_(bodyVel),
      bodyForce_(std::vector<double>(3,0.0)),
      bodyMom_(std::vector<double>(3,0.0)),
      appliedForce_(appliedForce),
      bodyAccel_(std::vector<double>(3,0.0)),
      bodyMass_(bodyMass),
      bodyDen_(bodyDen),
      forceSurface_(forceSurface)
{

  for (size_t i = 0; i < tetherGeom.size(); ++i) {
    this->tetherGeom_.emplace_back(std::array<double,9>{});
    for (size_t j = 0; j < 9; ++j)
      tetherGeom_.back()[j] = tetherGeom[i][j];

  }  

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
MeshMotionInfo::~MeshMotionInfo()
{
  // nothing to do
}
 
} // namespace nalu
} // namespace Sierra
