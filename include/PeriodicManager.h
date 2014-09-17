/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef PeriodicManager_h
#define PeriodicManager_h

//==============================================================================
// Includes and forwards
//==============================================================================

#include <FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Ghosting.hpp>

#include <stk_search/BoundingBox.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>

#include <vector>
#include <list>
#include <map>

namespace sierra {
namespace nalu {

class Realm;

typedef stk::search::IdentProc<stk::mesh::EntityKey,unsigned> theEntityKey;
typedef stk::search::Point<double> Point;
typedef stk::search::Sphere<double> Sphere;
typedef std::pair<Sphere,theEntityKey> sphereBoundingBox;

class PeriodicManager {

 public:

  // constructor and destructor
  PeriodicManager(
    Realm & realm);

  ~PeriodicManager();

  void add_periodic_pair(
    stk::mesh::Part* meshPartsMaster,
    stk::mesh::Part* meshPartsSlave,
    const double &searchTolerance,
    const std::string &searchMethodName);

  void build_constraints();

  // holder for master += slave; slave = master
  void apply_constraints(
    stk::mesh::FieldBase *,
    const unsigned &sizeOfField,
    const bool &bypassFieldCheck,
    const bool &addSlaves = true,
    const bool &setSlaves = true);

  // find the max
  void apply_max_field(
    stk::mesh::FieldBase *,
    const unsigned &sizeOfField);

  void create_ghosting_object();

  stk::mesh::Ghosting * get_ghosting_object();

  const stk::mesh::PartVector &get_slave_part_vector();

  double get_search_time();

 private:

  void setup_gid_pairs(
    const stk::mesh::Part *partMaster,
    const stk::mesh::Part *partSlave,
    const stk::search::SearchMethod searchMethod);

  void master_slave_reduction();

  void create_slave_part_vector();

  void update_global_id_field();

  /* communicate periodicGhosting nodes */
  void
  periodic_parallel_communicate_field(
    stk::mesh::FieldBase *theField);

  /* communicate shared nodes and aura nodes */
  void
  parallel_communicate_field(
    stk::mesh::FieldBase *theField);

  Realm &realm_;
  double searchTolerance_;

  stk::mesh::Ghosting *periodicGhosting_;
  const std::string ghostingName_;
  double timerSearch_;

  // the data structures to hold master/slave information
  typedef std::pair<stk::mesh::Part*, stk::mesh::Part*> MeshPartPair;
  typedef std::pair<stk::mesh::Entity, stk::mesh::Entity> EntityPair;
  typedef std::vector<std::pair<theEntityKey,theEntityKey> > SearchKeyVector;

  // vector of master:slave part pairs
  std::vector<MeshPartPair> periodicPartPairs_;

  // vector of search types
  std::vector<stk::search::SearchMethod> searchMethodVec_;

  // vector of masterEntity:slaveEntity
  std::vector<EntityPair> masterSlaveCommunicator_;

  // vector of slave parts
  stk::mesh::PartVector slavePartVector_;

  // culmination of all searches
  SearchKeyVector searchKeyVector_;

 private:

  void add_slave_to_master(
    stk::mesh::FieldBase *theField,
    const unsigned &sizeOfField,
    const bool &bypassFieldCheck);

  void set_slave_to_master(
    stk::mesh::FieldBase *theField,
    const unsigned &sizeOfField,
    const bool &bypassFieldCheck);

};

} // namespace nalu
} // namespace sierra

#endif
