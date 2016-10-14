/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for InitGoogleTest, etc
#include <mpi.h>                        // for MPI_Comm_rank, MPI_Finalize, etc

#include "include/NaluEnv.h"

int gl_argc = 0;
char** gl_argv = 0;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    //NaluEnv will call MPI_Finalize for us.
    sierra::nalu::NaluEnv::self();

    testing::InitGoogleTest(&argc, argv);

    gl_argc = argc;
    gl_argv = argv;

    int returnVal = RUN_ALL_TESTS();

    return returnVal;
}

