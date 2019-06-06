/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef NaluEnv_h
#define NaluEnv_h

#include <mpi.h>
#include <fstream>
#include <streambuf>

namespace sierra{
namespace nalu{
  
  class NaluEmptyStreamBuffer : public std::filebuf {
  public:
    int overflow(int c) {return c;}
  };

class NaluEnv
{
 public:

  NaluEnv();
  ~NaluEnv();

  static NaluEnv &self();

  MPI_Comm parallelCommunicator_;
  int pSize_;
  int pRank_;
  std::streambuf *stdoutStream_;
  std::ostream *naluLogStream_;
  std::ostream *naluParallelStream_;
  bool parallelLog_;
  std::string baseName_;
  
  NaluEmptyStreamBuffer naluEmptyStreamBuffer_;
  std::filebuf naluStreamBuffer_;
  std::filebuf naluParallelStreamBuffer_;

  std::ostream & naluOutputP0();
  std::ostream & naluOutput();

  MPI_Comm parallel_comm();
  int parallel_size();
  int parallel_rank();
  void set_log_file_stream(std::string naluLogName, bool pprint = false);
  void close_log_file_stream();
  double nalu_time();
  void set_base_name(std::string baseName);
  std::string get_base_name();
};

} // namespace nalu
} // namespace Sierra

#endif
