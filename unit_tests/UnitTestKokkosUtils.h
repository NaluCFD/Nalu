#ifndef _UnitTestKokkosUtils_h_
#define _UnitTestKokkosUtils_h_

#include <master_element/MasterElement.h>
#include <Kokkos_Core.hpp>

typedef Kokkos::Schedule<Kokkos::Dynamic> DynamicScheduleType;
typedef typename Kokkos::TeamPolicy<typename Kokkos::DefaultExecutionSpace, DynamicScheduleType>::member_type TeamHandleType;

using DeviceShmem = Kokkos::DefaultExecutionSpace::scratch_memory_space;
template<typename T>
using SharedMemView = Kokkos::View<T, Kokkos::LayoutRight, DeviceShmem, Kokkos::MemoryUnmanaged>;
using DeviceTeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
using DeviceTeam = DeviceTeamPolicy::member_type;

inline DeviceTeamPolicy get_team_policy(const size_t sz, const size_t bytes_per_team,
    const size_t bytes_per_thread)
{
  DeviceTeamPolicy policy(sz, Kokkos::AUTO);
  return policy.set_scratch_size(0, Kokkos::PerTeam(bytes_per_team), Kokkos::PerThread(bytes_per_thread));
}

inline
SharedMemView<double*> get_shmem_view_1D(const TeamHandleType& team, size_t len)
{
  return Kokkos::subview(SharedMemView<double**>(team.team_shmem(), team.team_size(), len), team.team_rank(), Kokkos::ALL());
}

inline
SharedMemView<double**> get_shmem_view_2D(const TeamHandleType& team, size_t len1, size_t len2)
{
  return Kokkos::subview(SharedMemView<double***>(team.team_shmem(), team.team_size(), len1, len2), team.team_rank(), Kokkos::ALL(), Kokkos::ALL());
}

template<class OUTER_LOOP_BODY, class INNER_LOOP_BODY>
void bucket_loop_serial_only(const stk::mesh::BucketVector& buckets, const OUTER_LOOP_BODY& outer_loop_body, const INNER_LOOP_BODY& inner_loop_body)
{
    for(const stk::mesh::Bucket* bptr : buckets)
    {   
        const stk::mesh::Bucket& bkt = *bptr;
        stk::topology topo = bkt.topology();
        sierra::nalu::MasterElement* meSCS = unit_test_utils::get_surface_master_element(topo);

        outer_loop_body(topo,*meSCS);

        for(size_t j=0; j<bkt.size(); ++j)
        {
            inner_loop_body(bkt[j], topo, *meSCS);
        }
    }
}

template<class LOOP_BODY>
void kokkos_bucket_loop(const stk::mesh::BucketVector& buckets, LOOP_BODY inner_loop_body)
{
    Kokkos::parallel_for(buckets.size(), [&](const size_t& i)
    {
        const stk::mesh::Bucket& bkt = *buckets[i];
        for(size_t j=0; j<bkt.size(); ++j)
        {
            inner_loop_body(bkt[j]);
        }
    });
}

template<class LOOP_BODY>
void kokkos_thread_team_bucket_loop(const stk::mesh::BucketVector& buckets, LOOP_BODY inner_loop_body)
{
    Kokkos::parallel_for(Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>(buckets.size(), Kokkos::AUTO), KOKKOS_LAMBDA(const TeamHandleType& team)
    {
        const stk::mesh::Bucket& bkt = *buckets[team.league_rank()];
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, bkt.size()), [&](const size_t& j)
        {
            inner_loop_body(bkt[j]);
        });
    });
}

template<class LOOP_BODY>
void kokkos_thread_team_bucket_loop_with_topo(const stk::mesh::BucketVector& buckets,
                                    const LOOP_BODY& inner_loop_body)
{
    Kokkos::parallel_for(Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>(buckets.size(), Kokkos::AUTO), KOKKOS_LAMBDA(const TeamHandleType& team)
    {
        const stk::mesh::Bucket& bkt = *buckets[team.league_rank()];
        stk::topology topo = bkt.topology();
        sierra::nalu::MasterElement* meSCS = unit_test_utils::get_surface_master_element(topo);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, bkt.size()), [&](const size_t& j)
        {
            inner_loop_body(bkt[j], topo, *meSCS);
        });
    });
}

#endif

