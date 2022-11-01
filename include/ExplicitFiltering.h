/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ExplicitFiltering_h
#define ExplicitFiltering_h

#include <NaluParsing.h>
#include <FieldTypeDef.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Ghosting.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/SearchMethod.hpp>

// basic c++
#include <string>
#include <vector>
#include <utility>

namespace sierra{
namespace nalu{

// common type defs
typedef stk::search::IdentProc<uint64_t,int> theKey;
typedef stk::search::Point<double> Point;
typedef stk::search::Box<double> Box;
typedef std::pair<Box,theKey> boundingElementBox;

struct ExplicitFilteringNames {
  std::string fieldName_;
  std::string expFieldName_;
  int fieldSize_;

ExplicitFilteringNames(std::string fieldName, std::string expFieldName, int fieldSize) 
: fieldName_(fieldName), expFieldName_(expFieldName), fieldSize_(fieldSize) {}
};

struct ExplicitFilteringFields {
  stk::mesh::Field<double> *theField_;
  stk::mesh::Field<double> *expField_;
  int fieldSize_;

ExplicitFilteringFields(stk::mesh::Field<double> *theField, stk::mesh::Field<double> *expField, int fieldSize) 
: theField_(theField), expField_(expField), fieldSize_(fieldSize) {}
};
 
class Realm;

class ExplicitFiltering
{
public:

  ExplicitFiltering(
    Realm &realm,
    const YAML::Node &node);
  ~ExplicitFiltering();

  // load all of the options
  void load(
    const YAML::Node & node);

  // setup part creation and nodal field registration (before populate_mesh())
  void setup();

  // setup part creation and nodal field registration (after populate_mesh())
  void initialize();

  // determine element bounding box in the mesh
  void populate_candidate_elements();

  // fill in the map that will hold point and ghosted elements
  void create_explicit_filter_point_info_map();

  // figure out the set of elements that belong in the custom ghosting data structure
  void determine_elems_to_ghost();

  // deal with custom ghosting
  void manage_ghosting();

  // populate vector of elements
  void complete_search();

  // populate nodal field and output norms (if appropriate)
  void execute();

  // general gather methods for scalar and vector (both double)
  void gather_field(
    const int &sizeOfField,
    double *elemAveragedField,
    const stk::mesh::FieldBase &stkField,
    stk::mesh::Entity const* elem_node_rels,
    const int &nodesPerElement);

  // element volume and scv volume populated
  double compute_volume(
    const int &nDim,
    stk::mesh::Entity elem,
    const stk::mesh::BulkData & bulkData);

  void increment_elem_mean(
    const int &sizeOfField,
    double *fieldToFill,
    const stk::mesh::FieldBase &stkField,
    stk::mesh::Entity const* elem_node_rels,
    const int &nodesPerElement,
    const double &elemVolume);

  // hold the realm
  Realm &realm_;

  // type of stk search
  const stk::search::SearchMethod searchMethod_;

  // the size of the filter
  Coordinates filterSize_;

  // custom ghosting
  stk::mesh::Ghosting *explicitFilteringGhosting_;

  // how many elements to ghost?
  uint64_t needToGhostCount_;

  // provide debug output
  bool debugOutput_;

  // vector of Names struct
  std::vector<ExplicitFilteringNames> explicitFilteringNamesVec_;

  // vector of Fields struct
  std::vector<ExplicitFilteringFields> explicitFilteringFieldsVec_;

  // vector of elements/poit processor to ghost
  stk::mesh::EntityProcVec elemsToGhost_;

  // save off product of search
  std::vector<std::pair<theKey, theKey> > searchKeyPair_;

  // bounding box data types for stk_search */
  std::vector<boundingElementBox> boundingFilterVec_;
  std::vector<boundingElementBox> boundingElementBoxVec_;

  // target names for set of bounding boxes
  std::vector<std::string> searchTargetNames_;

  // corresponding parts for targets
  stk::mesh::PartVector searchParts_;

  // map of point info objects to vector of element entities
  std::map<stk::mesh::Entity, std::vector<stk::mesh::Entity> > explicitFilteringMap_;

  // scratch space
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scv_volume_;
};


} // namespace nalu
} // namespace Sierra

#endif
