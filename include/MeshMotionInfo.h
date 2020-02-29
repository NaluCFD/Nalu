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

namespace sierra{
namespace nalu{

class MeshMotionInfo
{
 public:
  MeshMotionInfo(
   std::vector<std::string> meshMotionBlock, 
   const double omega, 
   std::vector<double> centroid,
   std::vector<double> unitVec,
   const bool computeCentroid,
   const double theAngle_ = 0.0);

  ~MeshMotionInfo();

  std::vector<std::string> meshMotionBlock_;
  const double omega_;
  std::vector<double> centroid_;
  std::vector<double> unitVec_;
  const double computeCentroid_;
  double computeCentroidCompleted_;
  const double theAngle_;
};

} // namespace nalu
} // namespace Sierra

#endif
