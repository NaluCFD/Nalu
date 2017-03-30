/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ActuatorLineFAST_h
#define ActuatorLineFAST_h

#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include "Actuator.h"

// FAST C++ API
#include "OpenFAST.H"

namespace sierra{
namespace nalu{

class Realm;

class ActuatorLineFASTInfo {
public:
  ActuatorLineFASTInfo();
  ~ActuatorLineFASTInfo();

  // for each type of probe, e.g., line of site, hold some stuff
  int processorId_;
  int numPoints_;
  std::string turbineName_;
  Coordinates epsilon_;
};

// class that holds all of the action... for each point, hold the current location and other useful info
class ActuatorLineFASTPointInfo {
 public:
  ActuatorLineFASTPointInfo(
			    size_t localId, Point centroidCoords, double searchRadius, Coordinates epsilon, double *velocity, fast::ActuatorNodeType nType, size_t globTurbId);
  ~ActuatorLineFASTPointInfo();
  size_t globTurbId_; // Global turbine number
  size_t localId_;
  Point centroidCoords_;
  double searchRadius_;
  Coordinates epsilon_;
  double bestX_;
  stk::mesh::Entity bestElem_;

  fast::ActuatorNodeType nodeType_;
  // mesh motion specifics
  double velocity_[3];

  std::vector<double> isoParCoords_;
  std::vector<stk::mesh::Entity> elementVec_;
};
 
 class ActuatorLineFAST: public Actuator
{
public:
  
  ActuatorLineFAST(
    Realm &realm,
    const YAML::Node &node);
  ~ActuatorLineFAST();
  
  // load all of the options
  void load(
    const YAML::Node & node);

  // setup part creation and nodal field registration (before populate_mesh())
  void setup();

  // allocate turbines to processors containing hub location
  void allocateTurbinesToProcs() ;
  
  // Allocate turbines to processors, initialize FAST and get location of actuator points
  void initialize();

  // setup part creation and nodal field registration (after populate_mesh())
  void update();

  // determine processor bounding box in the mesh
  void populate_candidate_procs();

  // determine element bounding box in the mesh
  void populate_candidate_elements();

  // fill in the map that will hold point and ghosted elements
  void create_actuator_line_point_info_map();

  // figure out the set of elements that belong in the custom ghosting data structure
  void determine_elems_to_ghost();

  // deal with custom ghosting
  void manage_ghosting();

  // manage rotation, now only in the y-z plane
  void set_current_coordinates(
    double *lineCentroid, double *centroidCoords, const double &omega, const double &currentTime);
  void set_current_velocity(
    double *lineCentroid, const double *centroidCoords, double *velocity, const double &omega);

  // populate vector of elements
  void complete_search();
    
  // populate nodal field and output norms (if appropriate)
  void execute();

  // support methods to gather data; scalar and vector
  void resize_std_vector( 
    const int &sizeOfField,
    std::vector<double> &theVector,   
    stk::mesh::Entity elem, 
    const stk::mesh::BulkData & bulkData);

  // general gather methods for scalar and vector (both double)
  void gather_field(
    const int &sizeOfField,
    double *fieldToFill, 
    const stk::mesh::FieldBase &stkField,
    stk::mesh::Entity const* elem_node_rels, 
    const int &nodesPerElement);

  void gather_field_for_interp(
    const int &sizeOfField,
    double *fieldToFill, 
    const stk::mesh::FieldBase &stkField,
    stk::mesh::Entity const* elem_node_rels, 
    const int &nodesPerElement);

  // element volume and scv volume populated
  double compute_volume( 
    const int &nDim,
    stk::mesh::Entity elem, 
    const stk::mesh::BulkData & bulkData);

  // interpolate field to point centroid
  void interpolate_field(
    const int &sizeOfField,
    stk::mesh::Entity elem, 
    const stk::mesh::BulkData & bulkData,
    double *isoParCoords,
    const double *fieldAtNodes,
    double *pointField);

  // centroid of the element
  void compute_elem_centroid( 
    const int &nDim,
    double *elemCentroid,
    const int &nodesPerElement);

  // distance from element centroid to point centroid
  double compute_distance( 
    const int &nDim,
    const double *elemCentroid,
    const double *pointCentroid);

  // compute the body force at an element given a
  // projection weighting.
  void compute_elem_force_given_weight(
    const int &nDim,
    const double &g,
    const double *pointForce,
    double *elemForce);

  // isotropic Gaussian projection function.
  double isotropic_Gaussian_projection(
    const int &nDim,
    const double &dis,
    const Coordinates &epsilon);

  // finally, perform the assembly
  void assemble_source_to_nodes(
    const int &nDim,
    stk::mesh::Entity elem,
    const stk::mesh::BulkData & bulkData,
    const double &elemVolume,
    const double *drag,
    const double &dragLHS,
    const double &gLocal,
    stk::mesh::FieldBase &actuator_source,
    stk::mesh::FieldBase &actuator_source_lhs,
    stk::mesh::FieldBase &g,
    const double &lhsFac);

  // hold the realm
  Realm &realm_;

  // type of stk search
  stk::search::SearchMethod searchMethod_;
  
  // custom ghosting
  stk::mesh::Ghosting *actuatorLineGhosting_;
  // how many elements to ghost?
  uint64_t needToGhostCount_;
  stk::mesh::EntityProcVec elemsToGhost_;

  // custom ghosting
  stk::mesh::Ghosting *actuatorLineForceGhosting_;
  // how many elements to ghost?
  uint64_t needToGhostForceCount_;
  stk::mesh::EntityProcVec elemsToGhostForce_;

  // does the actuator line move?
  bool actuatorLineMotion_;

  // everyone needs pi
  const double pi_;

  // save off product of search
  std::vector<std::pair<theKey, theKey> > searchKeyPair_;
  std::vector<std::pair<theKey, theKey> > searchKeyPairForce_;

  // bounding box data types for stk_search */
  std::vector<boundingSphere> boundingSphereVec_;
  std::vector<boundingElementBox> boundingElementBoxVec_;
  std::vector<boundingSphere> boundingHubSphereVec_;
  std::vector<boundingElementBox> boundingProcBoxVec_;
  std::vector<boundingSphere> boundingSphereForceVec_;

  // target names for set of bounding boxes
  std::vector<std::string> searchTargetNames_;
 
  // vector of averaging information
  std::vector<ActuatorLineFASTInfo *> actuatorLineInfo_; 

  // map of point info objects
  std::map<size_t, ActuatorLineFASTPointInfo *> actuatorLinePointInfoMap_;

  // map of point info objects
  std::map<size_t, ActuatorLineFASTPointInfo *> actuatorLineForcePointInfoMap_;

  // scratch space
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scv_volume_;
  std::vector<double> ws_velocity_;
  std::vector<double> ws_density_;
  std::vector<double> ws_viscosity_;

  // FAST C++ API handle
  fast::OpenFAST FAST;

};


} // namespace nalu
} // namespace Sierra

#endif
