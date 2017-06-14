/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef UNITTESTREALM_H
#define UNITTESTREALM_H

#include "Simulation.h"
#include "Realm.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include "yaml-cpp/yaml.h"

namespace unit_test_utils {

YAML::Node get_default_inputs();

YAML::Node get_realm_default_node();

class NaluTest
{
public:
  NaluTest(YAML::Node&);

  sierra::nalu::Realm& create_realm(
    const YAML::Node&,
    const std::string realm_type="multi_physics");

  YAML::Node doc_;
  stk::ParallelMachine comm_;
  unsigned spatialDim_;
  sierra::nalu::Simulation sim_;

  stk::mesh::PartVector partVec_;

private:
  NaluTest() = delete;
  NaluTest(const NaluTest&) = delete;
};

}  // unit_test_utils

#endif /* UNITTESTREALM_H */
