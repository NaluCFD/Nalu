/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "UnitTestFieldUtils.h"

namespace unit_test_utils {

double field_norm(const ScalarFieldType& field, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
{
  stk::ParallelMachine comm = bulk.parallel();
  const auto& buckets = bulk.get_buckets(stk::topology::NODE_RANK, selector);

  size_t N = 0;
  size_t g_N = 0;
  double norm = 0.0;
  double g_norm = 0.0;

  kokkos_thread_team_bucket_loop(buckets, [&](stk::mesh::Entity node) {
      double node_value = *stk::mesh::field_data(field, node);
      Kokkos::atomic_add(&N, (size_t)1);
      Kokkos::atomic_add(&norm, (node_value * node_value));
    });
  stk::all_reduce_sum(comm, &N, &g_N, 1);
  stk::all_reduce_sum(comm, &norm, &g_norm, 1);
  g_norm = std::sqrt(g_norm/g_N);

  return g_norm;
}

}
