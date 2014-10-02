/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef NaluEnv_h
#define NaluEnv_h

#include <mpi.h>
#include <string>
#include <fstream>

namespace sierra{
namespace nalu{

class NaluEnv : std::ostream
{
public:

  NaluEnv();
  ~NaluEnv();

  static NaluEnv &self();

  MPI_Comm parallelCommunicator_;
  int pSize_;
  int pRank_;
  std::ostream *naluLogStream_;

  std::ostream & naluOutputP0();
  MPI_Comm parallel_comm();
  int parallel_size();
  int parallel_rank();
  void set_log_file_stream(std::ofstream *str);

  template<class T>
    NaluEnv& operator<<(T& thing) {
    if ( pRank_ )
      (*naluLogStream_) << thing;
    return *naluLogStream_;
  }

};

} // namespace nalu
} // namespace Sierra

#endif
