/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "UnitTestFieldUtils.h"

namespace unit_test_utils {

double field_max(const ScalarFieldType& field, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
{
  stk::ParallelMachine comm = bulk.parallel();
  const auto& buckets = bulk.get_buckets(stk::topology::NODE_RANK, selector);

  double maxVal = -1.0e16;
  double g_maxVal = -1.0e16;
  kokkos_thread_team_bucket_loop(buckets, [&](stk::mesh::Entity node) {
      double nodeVal = *stk::mesh::field_data(field, node);
      maxVal = std::max(maxVal, nodeVal);
    });
  stk::all_reduce_max(comm, &maxVal, &g_maxVal, 1);

  return g_maxVal;
}

double field_min(const ScalarFieldType& field, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
{
  stk::ParallelMachine comm = bulk.parallel();
  const auto& buckets = bulk.get_buckets(stk::topology::NODE_RANK, selector);

  double minVal = 1.0e16;
  double g_minVal = 1.0e16;
  kokkos_thread_team_bucket_loop(buckets, [&](stk::mesh::Entity node) {
      double nodeVal = *stk::mesh::field_data(field, node);
      minVal = std::min(minVal, nodeVal);
    });
  stk::all_reduce_min(comm, &minVal, &g_minVal, 1);

  return g_minVal;
}


double field_norm(const ScalarFieldType& field, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
{
  stk::ParallelMachine comm = bulk.parallel();
  const auto& buckets = bulk.get_buckets(stk::topology::NODE_RANK, selector);

  size_t N = 0;
  size_t g_N = 0;
  double norm = 0.0;
  double g_norm = 0.0;

  kokkos_thread_team_bucket_loop(buckets, [&](stk::mesh::Entity node) {
      double nodeVal = *stk::mesh::field_data(field, node);
      norm += nodeVal*nodeVal;
      N++;
    });
  stk::all_reduce_sum(comm, &N, &g_N, 1);
  stk::all_reduce_sum(comm, &norm, &g_norm, 1);
  g_norm = std::sqrt(g_norm/g_N);

  return g_norm;
}

}
