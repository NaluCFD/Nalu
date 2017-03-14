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
   double omega, 
   std::vector<double> centroid,
   std::vector<double> unitVec);

  ~MeshMotionInfo();

  std::vector<std::string> meshMotionBlock_;
  double omega_;
  std::vector<double> centroid_;
  std::vector<double> unitVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
