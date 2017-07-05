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
MeshMotionInfo::MeshMotionInfo(
  std::vector<std::string> meshMotionBlock, 
  const double omega, 
  std::vector<double> centroid,
  std::vector<double> unitVec,
  const bool computeCentroid)
  : meshMotionBlock_(meshMotionBlock), 
    omega_(omega), 
    centroid_(centroid),
    unitVec_(unitVec),
    computeCentroid_(computeCentroid),
    computeCentroidCompleted_(false)
{
  // nothing to do
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
