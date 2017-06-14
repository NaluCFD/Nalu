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

    kokkos_thread_team_bucket_loop(buckets, [&](stk::mesh::Entity node) {

      double lhs_value = 0.0;
      double rhs_value = 0.0;

      for (size_t i=0; i < activeSuppAlgs_.size(); ++i){
        activeSuppAlgs_[i]->node_execute(&lhs_value, &rhs_value, node);
        lhs_.push_back(lhs_value);
        rhs_.push_back(rhs_value);
      }
      });
  }

  /** Compute norm of the LHS/RHS
   *
   * Used to generate the gold values and testing
   */
  double norm(const std::vector<double> & vec)
  {
    stk::ParallelMachine comm = bulk_.parallel();

    size_t N = vec.size();
    size_t g_N = 0;
    double norm = 0.0;
    double g_norm = 0.0;

    for (int i = 0; i < N; ++i) {
      norm += vec[i]*vec[i];
    }
    stk::all_reduce_sum(comm, &N, &g_N, 1);
    stk::all_reduce_sum(comm, &norm, &g_norm, 1);
    g_norm = std::sqrt(g_norm/g_N);

    return g_norm;
  }

  inline double get_lhs_norm() {return norm(lhs_);}
  inline double get_rhs_norm() {return norm(rhs_);}

  std::vector<sierra::nalu::SupplementalAlgorithm*> activeSuppAlgs_;

private:
  const stk::mesh::BulkData& bulk_;
  std::vector<double> lhs_;
  std::vector<double> rhs_;
};

}

#endif /* UNITTESTALGORITHMUTILS_H */
