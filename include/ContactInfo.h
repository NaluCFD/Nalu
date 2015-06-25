/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ContactInfo_h
#define ContactInfo_h

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
class HaloInfo;

typedef stk::search::IdentProc<uint64_t,int>  theKey;
typedef stk::search::Point<double> Point;
typedef stk::search::Box<double> Box;
typedef std::pair<Point,theKey> boundingPoint;
typedef std::pair<Box,theKey> boundingElementBox;

//=============================================================================
// Class Definition
//=============================================================================
// ContactInfo
//=============================================================================
/**
 * * @par Description:
 * - class to manage contact information.
 *
 * @par Design Considerations:
 * -
 */
//=============================================================================
class ContactInfo {

 public:

  // constructor and destructor
  ContactInfo(
    Realm & realm,
    const std::string &name,
    const double maxSearchRadius,
    const double minSearchRadius,
    const std::vector<std::string> &contactSearchBlockName,
    const double expandBoxPercentage,
    const stk::mesh::Part *contactSurfacePart,
    const std::string &searchMethodName,
    const bool clipIsoParametricCoords,
    const bool useHermiteInterpolation);

  ~ContactInfo();

  void initialize();
  void construct_halo_state();
  void populate_halo_mesh_velocity();
  void find_possible_elements();
  void set_best_x();
  void determine_elems_to_ghost();
  void complete_search();
  void dump_diagnosis();

  Realm &realm_;
  const std::string name_;

   /* user input variable for max R */
  const double maxSearchRadius_;

  /* user input variable for min R */
  const double minSearchRadius_;

  // opposing search and contact bc part
  stk::mesh::PartVector contactSearchBlock_;
  const stk::mesh::Part *contactSurfacePart_;

  /* master element for homogeneous block search type; one topo per contactinfo */
  MasterElement *meSCS_;

  /* expand search box */
  double expandBoxPercentage_;

  /* does the realm have mesh motion */
  const bool meshMotion_;

  stk::search::SearchMethod searchMethod_;

  /* clip isoparametric coordinates if they are out of bounds */
  const bool clipIsoParametricCoords_;

  /* use higher order Hermite interpolation */
  const bool useHermiteInterpolation_;

  /* bounding box data types for stk_search */
  std::vector<boundingPoint>      boundingPointVec_;
  std::vector<boundingElementBox> boundingElementBoxVec_;

  /* map of HaloInfo */
  std::map<uint64_t, HaloInfo *> haloInfoMap_;

  /* save off product of search */
  std::vector<std::pair<theKey, theKey> > searchKeyPair_;

};

} // end sierra namespace
} // end Acon namespace

#endif
