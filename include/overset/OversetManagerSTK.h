/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef OversetManagerSTK_h
#define OversetManagerSTK_h

#include "overset/OversetManager.h"

// stk_mesh
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

// stk_search
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/SearchMethod.hpp>

// STL
#include <vector>
#include <map>

// field types
typedef stk::mesh::Field<double>  ScalarFieldType;
typedef stk::mesh::Field<double, stk::mesh::Cartesian>  VectorFieldType;
typedef stk::mesh::Field<double, stk::mesh::SimpleArrayTag>  GenericFieldType;

// search types
typedef stk::search::IdentProc<uint64_t,int>  theKey;
typedef stk::search::Point<double> Point;
typedef stk::search::Box<double> Box;
typedef std::pair<Point,theKey> boundingPoint;
typedef std::pair<Box,theKey> boundingElementBox;

namespace stk {
  namespace io {
    class StkMeshIoBroker;
  }
  namespace mesh {
    class Part;
    class MetaData;
    class BulkData;
    class Ghosting;
    typedef std::vector<Part*> PartVector;
    struct Entity;
  }
}

namespace sierra {
namespace nalu{
class OversetInfo;
struct OversetUserData;

/** Overset connectivity with native STK hole cutting algorithm
 *
 *  Example usage:
 *
 *  ```
 *  - overset_boundary_condition: bc_overset
 *    overset_user_data:
 *      percent_overlap: 10.0
 *      background_block: block_1
 *      overset_block: block_2
 *      overset_surface: surface_6
 *      background_cut_block: block_3
 *      background_cut_surface: surface_101
 *  ```
 */
class OversetManagerSTK : public OversetManager
{
public:

  // constructor/destructor
  OversetManagerSTK(
    Realm & realm, 
    const OversetUserData &oversetUserData);

  virtual ~OversetManagerSTK();

  virtual void setup();

  // main method called for initialization
  virtual void initialize();
  
  // initialize ghosting data structures
  void initialize_ghosting();

  // set space for inactive part; intersection of overset with background mesh
  void declare_inactive_part();

  // set space for inactive part exposed surfaces
  void declare_background_surface_part();
  
  // set space for inner part
  void declare_inner_part();

  // define the bounding boxes for inactive and inner (based on % reduction of the overset mesh)
  void define_inactive_bounding_box();

  // define the high level overset bounding boxes
  void define_overset_bounding_boxes();

  // define the background mesh set of bounding boxes
  void define_background_bounding_boxes();

  // determine all of the inactive intersected elements (coarse search on inactive bounding box and backgroundBoxesVec)
  void determine_intersected_elements(
    std::vector<boundingElementBox> &boundingBoxVec,  
    std::vector<boundingElementBox> &boundingBoxesVec,  
    std::vector<stk::mesh::Entity > &elementVec);

  // remove all elements from internally managed parts
  void clear_parts();

  // add elements to the inactive part
  void populate_inactive_part();

  // add elements to the inner paer
  void populate_inner_part();

  // skin the inactive part to obtain a surface part
  void skin_exposed_surface_on_inactive_part();

  // push back on all surfaces that contain the orphan nodes
  void set_orphan_surface_part_vec();

  // create an OversetInfo object for each locally owned exposed node
  void create_overset_info_vec();

  // create an OversetInfo object for each locally owned node that is part of inactive, however, not inner part
  void create_fringe_info_vec();

  // orphan node within element search; product is a valid ghosting and oversetInfoVec completed
  void orphan_node_search();
  
  // set the element variable for intersected elements to unity
  void set_data_on_inactive_part();

  // set the nodal variable for fringe
  void set_data_on_fringe_part();

  // general coarse search method 
  void coarse_search( 
    std::vector<boundingPoint> &boundingPointVec,    
    std::vector<boundingElementBox> &boundingElementVec,
    std::vector<std::pair<theKey, theKey> > &searchKeyPair);

  // deal with ghosting objects/data for fine search purposes
  void manage_ghosting();
  
  // general complete search method (fine search)
  void complete_search( 
    std::vector<std::pair<theKey, theKey> > searchKeyPair,
    std::map<uint64_t, OversetInfo *> &oversetInfoMap);

  // data set at construction
  const OversetUserData &oversetUserData_;
  const stk::search::SearchMethod searchMethod_;
  int nDim_;
  // lots of detailed information on the search
  const bool oversetAlgDetailedOutput_;

  uint64_t needToGhostCount_; 

  // internal flag to ensure that we declare parts once
  bool firstInitialization_;
  
  // vector of elements to ghost
  stk::mesh::EntityProcVec elemsToGhost_;

  // search data structures
  std::vector<boundingElementBox> boundingElementInactiveBoxVec_;
  std::vector<boundingElementBox> boundingElementInactiveBoxVecInner_;
  std::vector<boundingElementBox> boundingElementOversetBoxesVec_;
  std::vector<boundingElementBox> boundingElementBackgroundBoxesVec_;
  std::vector<boundingPoint>      boundingPointVecBackground_;
  std::vector<boundingPoint>      boundingPointVecOverset_;
  std::vector<boundingPoint>      boundingPointVecInner_;

  /* save off product of search */
  std::vector<std::pair<theKey, theKey> > searchKeyPairBackground_;
  std::vector<std::pair<theKey, theKey> > searchKeyPairOverset_;
  std::vector<std::pair<theKey, theKey> > searchKeyPairInner_;
  
  // vector of elements intersected... will want to push to a part
  std::vector<stk::mesh::Entity > intersectedInactiveElementVec_;
  std::vector<stk::mesh::Entity > intersectedInactiveElementVecInner_;

  // hold a vector of parts that correspond to the exposed surface orphan points
  stk::mesh::PartVector orphanPointSurfaceVecOverset_;

  // map of node global id to overset info
  std::map<uint64_t, OversetInfo *> oversetInfoMapOverset_;
  std::map<uint64_t, OversetInfo *> oversetInfoMapBackground_;
  std::map<uint64_t, OversetInfo *> oversetInfoMapFringe_;

  // part and info that holds elements within the inner
  stk::mesh::Part* innerPart_{nullptr};  
};

} // namespace nalu
} // namespace Sierra

#endif
