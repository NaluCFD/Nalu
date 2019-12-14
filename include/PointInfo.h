/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef PointInfo_h
#define PointInfo_h

//==============================================================================
// Includes and forwards
//==============================================================================

#include <master_element/MasterElement.h>

// stk
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Ghosting.hpp>

// stk_search
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/SearchMethod.hpp>
#include <stk_search/CoarseSearch.hpp>

#include <vector>
#include <map>

namespace sierra {
namespace nalu {

typedef stk::search::IdentProc<uint64_t,int>  uint64IdentProc;
typedef stk::search::Point<double> Point;
typedef stk::search::Sphere<double> Sphere;
typedef stk::search::Box<double> Box;
typedef std::pair<Point,uint64IdentProc> boundingPoint;
typedef std::pair<Sphere,uint64IdentProc> boundingSphere;
typedef std::pair<Box,uint64IdentProc> boundingBox;

//=============================================================================
// Class Definition
//=============================================================================
// PointInfo
//=============================================================================
/**
 * * @par Description:
 * - class to manage pointinformation.
 *
 * @par Design Considerations:
 * -
 */
//=============================================================================
class PointInfo {

 public:

  // constructor and destructor
  PointInfo(
    boundingPoint bPoint,
    const uint64_t localPointId,
    Point &ipCoords,
    Point &pointCoords,
    const int nDim);
  ~PointInfo();

  boundingPoint bPoint_;
  // should be able to extract this below from bPoint, right?
  const uint64_t localPointId_;
  const Point ipCoordinates_;
  const Point pointCoordinates_;
  
  stk::mesh::Entity owningElement_;

  double bestX_;
  const double bestXRef_;

  int elemIsGhosted_;

  // master element for background mesh
  MasterElement *meSCS_;

  std::vector<double> isoParCoords_;
};

} // end sierra namespace
} // end Acon namespace

#endif
