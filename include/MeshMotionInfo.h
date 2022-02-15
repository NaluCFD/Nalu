/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MeshMotionInfo_h
#define MeshMotionInfo_h

// standard c++
#include <string>
#include <vector>
#include <array>

namespace sierra{
namespace nalu{

class MeshMotionInfo
{
 public:
  
  // Rotation constructor
  MeshMotionInfo(
   std::vector<std::string> meshMotionBlock, 
   const double omega, 
   std::vector<double> centroid,
   std::vector<double> unitVec,
   const bool computeCentroid,
   const double theAngle_ = 0.0);

  // 6-DOF constructor
  MeshMotionInfo(
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
    std::vector<std::array<double,9>> tetherGeom);

  ~MeshMotionInfo();

  std::vector<std::string> meshMotionBlock_;
  const double omega_;
  std::vector<double> centroid_;
  std::vector<double> unitVec_;
  const double computeCentroid_;

  std::vector<std::array<double,9>> tetherGeom_;

  double computeCentroidCompleted_;
  const double theAngle_;


  // General 6-DOF motion
  const bool sixDof_;
  std::vector<double> bodyDispCC_;
  std::vector<double> bodyAngle_;
  std::vector<double> bodyOmega_;
  std::vector<double> bodyPrincInertia_;
  std::vector<double> bodyVel_;
  std::vector<double> bodyForce_;
  std::vector<double> bodyMom_;
  std::vector<double> appliedForce_;
  std::vector<double> bodyAccel_;
  std::vector<double> bodyAlpha_;
  const double bodyMass_; 
  const double bodyDen_;
  std::vector<std::string> forceSurface_;
  
  
};

} // namespace nalu
} // namespace Sierra

#endif
