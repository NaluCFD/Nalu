#include <gtest/gtest.h>
#include <limits>

#include <stk_util/parallel/Parallel.hpp>
#include <Kokkos_Core.hpp>

TEST(BasicKokkos, discover_execution_space)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    int proc = stk::parallel_machine_rank(comm);

    if (proc == 0) {
        std::cout << std::endl;

#ifdef KOKKOS_HAVE_SERIAL
        std::cout << "Kokkos::Serial is available." << std::endl;
#endif

#ifdef KOKKOS_HAVE_OPENMP
        std::cout << "Kokkos::OpenMP is available. (Control num-threads via env-var OMP_NUM_THREADS)" << std::endl;
#endif

#ifdef KOKKOS_HAVE_CUDA
        std::cout << "Kokkos::Cuda is available." << std::endl;
#endif
        std::cout << "Default execution space info: ";
        Kokkos::DefaultExecutionSpace::print_configuration(std::cout);

        std::cout << std::endl;
    }
}

TEST(BasicKokkos, simple_views_1D)
{
    const double tolerance = 0.0000001;
    const size_t N = 10;
    Kokkos::View<double*>::HostMirror host_view1D("host_view1D", N);
    for(size_t i=0; i<N; ++i) {
        host_view1D(i) = i+1;
    }

    Kokkos::View<double*> device_view1D = Kokkos::create_mirror_view(host_view1D);
    Kokkos::deep_copy(device_view1D, host_view1D);

    Kokkos::View<double*>::HostMirror host_view1D_2("host_view1D_2", N);
    Kokkos::deep_copy(host_view1D_2, device_view1D);

    for(size_t i=0; i<N; ++i) {
        EXPECT_NEAR(host_view1D(i), host_view1D_2(i), tolerance);
    }
}

TEST(BasicKokkos, simple_views_2D)
{
    const double tolerance = 0.0000001;
    const size_t N = 10;
    const size_t M = 20;
    Kokkos::View<double**>::HostMirror host_view2D("host_view2D", N, M);
    for(size_t i=0; i<N; ++i) {
        for(size_t j=0; j<M; ++j) {
            host_view2D(i,j) = i+j+1;
        }
    }

    Kokkos::View<double**> device_view2D = Kokkos::create_mirror_view(host_view2D);
    Kokkos::deep_copy(device_view2D, host_view2D);

    Kokkos::View<double**>::HostMirror host_view2D_2("host_view2D_2", N, M);
    Kokkos::deep_copy(host_view2D_2, device_view2D);

    for(size_t i=0; i<N; ++i) {
        for(size_t j=0; j<M; ++j) {
            EXPECT_NEAR(host_view2D(i,j), host_view2D_2(i,j), tolerance);
        }
    }
}

TEST(BasicKokkos, parallel_for)
{
    const double tolerance = 0.0000001;
    const size_t N = 10;
    const size_t M = 20;
    Kokkos::View<double**>::HostMirror host_view2D("host_view2D", N, M);

    for(size_t i=0; i<N; ++i) {
        for(size_t j=0; j<M; ++j) {
            host_view2D(i,j) = i+j+1;
        }
    }

    Kokkos::View<double**> device_view2D = Kokkos::create_mirror_view(host_view2D);
    Kokkos::deep_copy(device_view2D, host_view2D);

    Kokkos::parallel_for(N, KOKKOS_LAMBDA(const size_t& i) {
        for(size_t j=0; j<M; ++j) {
            device_view2D(i, j) *= 2;
        }
    });

    Kokkos::View<double**>::HostMirror host_result("host_result", N, M);
    Kokkos::deep_copy(host_result, device_view2D);

    for(size_t i=0; i<N; ++i) {
        for(size_t j=0; j<M; ++j) {
            EXPECT_NEAR(host_result(i,j), host_view2D(i,j), tolerance);
        }
    }
}

TEST(BasicKokkos, nested_parallel_for_thread_teams)
{
    const double tolerance = 0.0000001;
    const size_t N = 100;
    const size_t M = 200;
    Kokkos::View<double**>::HostMirror host_view2D("host_view2D", N, M);

    for(size_t i=0; i<N; ++i) {
        for(size_t j=0; j<M; ++j) {
            host_view2D(i,j) = i+j+1;
        }
    }

    Kokkos::View<double**> device_view2D = Kokkos::create_mirror_view(host_view2D);
    Kokkos::deep_copy(device_view2D, host_view2D);

    typedef Kokkos::Schedule<Kokkos::Dynamic> DynamicScheduleType;
    typedef typename Kokkos::TeamPolicy<typename Kokkos::DefaultExecutionSpace, DynamicScheduleType>::member_type TeamHandleType;

    Kokkos::parallel_for(Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>(N, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamHandleType& team) {
            size_t i = team.league_rank();
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, (size_t)0, M), KOKKOS_LAMBDA(const size_t& j) {
                device_view2D(i, j) *= 2;
            });
        });

    Kokkos::View<double**>::HostMirror host_result("host_result", N, M);
    Kokkos::deep_copy(host_result, device_view2D);

    for(size_t i=0; i<N; ++i) {
        for(size_t j=0; j<M; ++j) {
            EXPECT_NEAR(host_result(i,j), host_view2D(i,j), tolerance);
        }
    }
}

