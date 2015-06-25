/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef NonConformalInfo_h
#define NonConformalInfo_h

//==============================================================================
// Includes and forwards
//==============================================================================

#include <master_element/MasterElement.h>

// stk
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Ghosting.hpp>

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/SearchMethod.hpp>

#include <vector>
#include <map>

namespace stk {
namespace mesh {
typedef std::vector<Part *> PartVector;
}
}

namespace sierra {
namespace nalu {

class Realm;
class DgInfo;

typedef stk::search::IdentProc<uint64_t,int>  theKey;
typedef stk::search::Point<double> Point;
typedef stk::search::Box<double> Box;
typedef std::pair<Point,theKey> boundingPoint;
typedef std::pair<Box,theKey> boundingElementBox;

//=============================================================================
// Class Definition
//=============================================================================
// NonConformalInfo
//=============================================================================
/**
 * * @par Description:
 * - class to manage dg information.
 *
 * @par Design Considerations:
 * -
 */
//=============================================================================
class NonConformalInfo {

 public:

  // constructor and destructor
  NonConformalInfo(
    Realm & realm,
    const stk::mesh::Part *currentPart,
    const stk::mesh::Part *opposingPart,
    const double expandBoxPercentage,
    const std::string &searchMethodName,
    const bool clipIsoParametricCoords,
    const double searchTolerance);

  ~NonConformalInfo();

  void initialize();
  void construct_dgInfo_state();
  void find_possible_face_elements();
  void set_best_x();
  void determine_elems_to_ghost();
  void complete_search();
  void provide_diagnosis();

  Realm &realm_;
  const std::string name_;

  // master slave parts; slave part can be subsetted while master is not..
  const stk::mesh::Part *currentPart_;
  const stk::mesh::Part *opposingPart_;

  /* expand search box */
  double expandBoxPercentage_;

  stk::search::SearchMethod searchMethod_;

  /* clip isoparametric coordinates if they are out of bounds */
  const bool clipIsoParametricCoords_;

  /* allow for some finite search tolereance for bounding box */
  const double searchTolerance_;

  /* does the realm have mesh motion */
  const bool meshMotion_;

  /* bounding box data types for stk_search */
  std::vector<boundingPoint>      boundingPointVec_;
  std::vector<boundingElementBox> boundingFaceElementBoxVec_;

  /* vector of DgInfo */
  std::vector<std::vector<DgInfo *> > dgInfoVec_;

  /* save off product of search */
  std::vector<std::pair<theKey, theKey> > searchKeyPair_;

};

} // end sierra namespace
} // end Acon namespace

#endif
