/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level NaluUnit      */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef OversetManager_h
#define OversetManager_h

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

class OversetManager
{
public:

  // constructor/destructor
  OversetManager(
    Realm & realm, 
    const OversetUserData &oversetUserData);

  ~OversetManager();

  // main method called for initialization
  void initialize();
  
  // allow for manager to populate orhan nodal values
  void overset_orphan_node_field_update(
     stk::mesh::FieldBase *theField,
     const int sizeRow,
     const int sizeCol);

  // initialize ghosting data structures
  void initialize_ghosting();

  // set space for inactive part; intersection of overset with background mesh
  void declare_inactive_part();

  // set space for inactive part exposed surfaces
  void declare_background_surface_part();
  
  // define the high level overset bounding box (single in size for cutting)
  void define_overset_bounding_box();

  // define the high level overset bounding boxes
  void define_overset_bounding_boxes();

  // define the background mesh set of bounding boxes
  void define_background_bounding_boxes();

  // determine all of the intersected elements (coarse search on oversetBoxVec and backgroundBoxVec)
  void determine_intersected_elements();

  // add elements to the inactive part
  void populate_inactive_part();

  // skin the inactive part to obtain a surface part
  void skin_exposed_surface_on_inactive_part();

  // push back on all surfaces that contain the orphan nodes
  void set_orphan_surface_part_vec();

  // create an OversetInfo object for each locally owned exposed node
  void create_overset_info_vec();

  // orphan node within element search; product is a valid ghosting and oversetInfoVec completed
  void orphan_node_search();
  
  // set the element variable for intersected elements to unity
  void set_data_on_inactive_part();

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
  Realm &realm_;
  const OversetUserData &oversetUserData_;
  const stk::search::SearchMethod searchMethod_;
  int nDim_;
  // meta, bulk and io
  stk::mesh::MetaData *metaData_;
  stk::mesh::BulkData *bulkData_;
  // lots of detailed information on the search
  const bool oversetAlgDetailedOutput_;

  /* ghosting infrastructure */
  stk::mesh::Ghosting *oversetGhosting_;
  uint64_t needToGhostCount_; 

  // part associated with inactive elements
  stk::mesh::Part *inActivePart_;
  
  // part associated with exposed surfaces for inactive elements
  stk::mesh::Part *backgroundSurfacePart_;

  // vector of elements to ghost
  stk::mesh::EntityProcVec elemsToGhost_;

  // search data structures
  std::vector<boundingElementBox> boundingElementOversetBoxVec_;
  std::vector<boundingElementBox> boundingElementOversetBoxesVec_;
  std::vector<boundingElementBox> boundingElementBackgroundBoxesVec_;
  std::vector<boundingPoint>      boundingPointVecBackground_;
  std::vector<boundingPoint>      boundingPointVecOverset_;

  // map for background elements (used for intersection)
  std::map<uint64_t, stk::mesh::Entity> searchIntersectedElementMap_;

  /* save off product of search */
  std::vector<std::pair<theKey, theKey> > searchKeyPairBackground_;
  std::vector<std::pair<theKey, theKey> > searchKeyPairOverset_;
  
  // vector of elements intersected... will want to push to a part
  std::vector<stk::mesh::Entity > intersectedElementVec_;

  // hold a vector of parts that correspond to the exposed surface orphan points
  stk::mesh::PartVector orphanPointSurfaceVecOverset_;
  stk::mesh::PartVector orphanPointSurfaceVecBackground_;

  // vector of overset information: contains orphan node -> donor element 
  std::vector<OversetInfo *> oversetInfoVec_;

  // map of node global id to overset info
  std::map<uint64_t, OversetInfo *> oversetInfoMapOverset_;
  std::map<uint64_t, OversetInfo *> oversetInfoMapBackground_;

};

} // namespace naluUnit
} // namespace Sierra

#endif
