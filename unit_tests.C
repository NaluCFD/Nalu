/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for InitGoogleTest, etc
#include <mpi.h>                        // for MPI_Comm_rank, MPI_Finalize, etc
#include <Kokkos_Core.hpp>
#include <stk_util/parallel/Parallel.hpp>

// can't use stk_unit_test_utils until Trilinos/stk is updated, configuration is changed...
// #include <stk_unit_test_utils/ParallelGtestOutput.hpp>

#include "include/NaluEnv.h"

int gl_argc = 0;
char** gl_argv = 0;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    //NaluEnv will call MPI_Finalize for us.
    sierra::nalu::NaluEnv::self();
    Kokkos::initialize(argc, argv);
    int returnVal = 0;

    // Create a dummy nested scope to ensure destructors are called before
    // Kokkos::finalize_all. The instances owning threaded Kokkos loops must be
    // cleared out before Kokkos::finalize is called.
    {
      testing::InitGoogleTest(&argc, argv);

      gl_argc = argc;
      gl_argv = argv;

// can't use stk_unit_test_utils until Trilinos/stk is updated, configuration is changed...
//    int procId = stk::parallel_machine_rank(MPI_COMM_WORLD);
//    stk::unit_test_util::create_parallel_output(procId);

      returnVal = RUN_ALL_TESTS();
    }

    Kokkos::finalize_all();

    //NaluEnv will call MPI_Finalize when the NaluEnv singleton is cleaned up,
    //which is after we return.
    return returnVal;
}

