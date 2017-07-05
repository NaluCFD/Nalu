/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef UNITTESTALGORITHMUTILS_H
#define UNITTESTALGORITHMUTILS_H

#include <gtest/gtest.h>
#include "UnitTestUtils.h"

#include "Algorithm.h"
#include "SupplementalAlgorithm.h"

namespace unit_test_algorithm_utils {

/** Driver class that mimics Assemble*SolverAlgorithm
 *
 * It is the caller's responsibility to populate `activeSuppAlgs_` and call the
 * `setup` method on the activated SuppAlg before calling the `execute` method
 * of this class.
 *
 */
class TestSupplementalAlgorithmDriver
{
public:
  TestSupplementalAlgorithmDriver(
    const stk::mesh::BulkData& bulk)
    : bulk_(bulk)
  {}

  void nodal_execute()
  {
    const stk::mesh::MetaData& meta = bulk_.mesh_meta_data();

    const stk::mesh::Selector selector = meta.locally_owned_part();

    const auto& buckets = bulk_.get_buckets(stk::topology::NODE_RANK,
                                            selector);

    double lhs_value = 0.0;
    double rhs_value = 0.0;

    kokkos_thread_team_bucket_loop(buckets, [&](stk::mesh::Entity node) {
      for (size_t i=0; i < activeSuppAlgs_.size(); ++i){
	lhs_value = 0.0;
	rhs_value = 0.0;

        activeSuppAlgs_[i]->node_execute(&lhs_value, &rhs_value, node);

        lhs_.push_back(lhs_value);
        rhs_.push_back(rhs_value);
      }
      });
  }

  inline double get_lhs_norm() {return unit_test_utils::vector_norm(lhs_, bulk_.parallel());}
  inline double get_rhs_norm() {return unit_test_utils::vector_norm(rhs_, bulk_.parallel());}

  std::vector<sierra::nalu::SupplementalAlgorithm*> activeSuppAlgs_;

private:
  const stk::mesh::BulkData& bulk_;
  std::vector<double> lhs_;
  std::vector<double> rhs_;
};

}

#endif /* UNITTESTALGORITHMUTILS_H */
